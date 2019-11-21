using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using viper.corelib.horizon;
using viper.corelib.interfaces;
using viper.corelib.patch;
using viper.corelib.spice;
using viper.corelib.utilities;

namespace viper.corelib.lighting
{
    public class SafeHavenGenerator
    {
        public ILunarHorizon MainWindow;
        public Rectangle Region;
        public bool WriteHorizons = true;

        public DateTime Start = new DateTime(2022, 12, 1);
        public DateTime Stop = new DateTime(2023, 4, 1);
        public TimeSpan Step = TimeSpan.FromHours(2);

        public int ObserverHeightInMeters = 0;
        public float EarthMultipathThreshold = 2f;

        public float SafeHavenSunThreshold = 0.5f;

        public void WriteSafeHavenGeotiffs(string path)
        {
            var outer_start = Start.AddDays(-28);
            var outer_stop = Stop.AddDays(28);
            GetLowEarthTimes(outer_start, outer_stop, out List<DateTime> times_start_to_stop, out List<int> indices_of_minima_earth_elevation);

            for (var i=0;i<indices_of_minima_earth_elevation.Count;i++)
                Console.WriteLine($"i={i} index={indices_of_minima_earth_elevation[i]} time={times_start_to_stop[indices_of_minima_earth_elevation[i]]}");
            var region_width = Region.Width;
            var region_height = Region.Height;

            var combined = new float[region_height, region_width];
            var month_images = Enumerable.Range(0, indices_of_minima_earth_elevation.Count).Select(i => new float[region_height, region_width]).ToList();

            var sun_vectors = times_start_to_stop.Select(time => CSpice.SunPosition(time)).ToList();
            var earth_vectors = times_start_to_stop.Select(time => CSpice.EarthPosition(time)).ToList();

            // Get the indexes in the time and vector arrays of Start and Stop
            var first_inside = times_start_to_stop.Select((time, index) => (time, index)).Where(pair => pair.time >= Start).Select(pair => pair.index).First();
            var last_inside = -1;
            for (var i = times_start_to_stop.Count - 1; i >= 0; i--)
                if (times_start_to_stop[i] <= Stop)
                {
                    last_inside = i;
                    break;
                }

            //var debugpt = new Point(Region.Left + 545, Region.Top + 490);
            //debugpt = new Point(Region.Left + 494, Region.Top + 438);
            var ids = TerrainPatch.EnumerateIds(TerrainPatch.CoveringIdRectangle(Region)).ToList();

            {
                var not_generated = ids.Where(id => !File.Exists(TerrainPatch.FromId(id, ObserverHeightInMeters).Path)).ToList();
                Console.WriteLine($"{not_generated.Count} patches haven't been generated yet.  Generating them.");
                var gpuProcessor = ViperEnvironment.Processor as CPUHorizons;
                var queue = not_generated.Select(id => TerrainPatch.FromId(id)).ToList();
                gpuProcessor.RunQueue(queue, writeHorizons: WriteHorizons, unloadHorizons: true);
            }

            Parallel.ForEach(ids,
                new ParallelOptions { MaxDegreeOfParallelism = ViperEnvironment.MaxDegreeOfParallelism },
                id =>
                {
                    var patch = GetPatchWithHorizons(id, ObserverHeightInMeters);
                    for (var row = 0; row < TerrainPatch.DefaultSize; row++)
                    {
                        var region_row = (patch.Line + row) - Region.Y;
                        if (region_row < 0 || region_row >= region_height) continue;
                        for (var col = 0; col < TerrainPatch.DefaultSize; col++)
                        {
                            //if ((patch.Line + row) != debugpt.Y || (patch.Sample + col) != debugpt.X)
                            //    continue;
                            var region_col = patch.Sample + col - Region.X;
                            if (region_col < 0 || region_col >= region_width) continue;
                            var has_sun = patch.GetLightCurve(row, col, sun_vectors).Select(sun => sun >= SafeHavenSunThreshold).ToList();
                            var has_comm = patch.GetEarthOverHorizonCurve(row, col, earth_vectors).Select(rise => rise >= EarthMultipathThreshold).ToList();
                            var month_steps = MaxShadowWithoutComm2(indices_of_minima_earth_elevation, has_sun, has_comm);
                            combined[region_row, region_col] = 2f * month_steps.Max();
                            for (var i = 0; i < month_steps.Count; i++)
                                month_images[i][region_row, region_col] = 2f * month_steps[i];
                        }
                    }
                });

            {
                var path_combined = LandingSiteDataset.AppendToFilename(path, "_combined");
                GeotiffHelper.WriteArrayAsGeotiff(combined, Region, path_combined);
                WriteSafeHavenOverlay(combined, path_combined);
            }

            for (var i = 0; i < month_images.Count; i++)
            {
                var time = times_start_to_stop[indices_of_minima_earth_elevation[i]];
                var path1 = LandingSiteDataset.AppendToFilename(path, "_" + time.ToString("yyyy-MM-dd"));
                GeotiffHelper.WriteArrayAsGeotiff(month_images[i], Region, path1);
                WriteSafeHavenOverlay(month_images[i], path1);
            }
        }

        // Return all times from start to stop.  Also return the indexes of the times at which the earth elevation reaches a minima
        void GetLowEarthTimes(DateTime outer_start, DateTime outer_stop, out List<DateTime> all_times, out List<int> indices_of_minima_earth_elevation)
        {
            var center = new Point(Region.Left + Region.Width / 2, Region.Top + Region.Height / 2);
            var centerPatch = ViperEnvironment.GetPatch(TerrainPatch.LineSampleToId(center));
            centerPatch.FillMatrices(ViperEnvironment.Terrain);

            // Generate the time / elevation pairs
            all_times = new List<DateTime>();
            for (var time = outer_start; time <= outer_stop; time += Step)
                all_times.Add(time);
            var earth_elevations = all_times.Select(t =>
            {
                centerPatch.GetAzEl(CSpice.EarthPosition(t), 0, 0, out float _, out float earth_elevation_rad);
                return earth_elevation_rad;
            }).ToList();

            // Find minima
            indices_of_minima_earth_elevation = new List<int>();
            for (var i = 1; i < earth_elevations.Count - 1; i++)
                if (Start <= all_times[i] && all_times[i] <= Stop && earth_elevations[i - 1] > earth_elevations[i] && earth_elevations[i] < earth_elevations[i + 1])
                    indices_of_minima_earth_elevation.Add(i);
        }

        // Return value in steps, not hours
        public int MaxShadowWithoutComm(List<bool> sun, List<bool> comm, int first_inside, int last_inside)
        {
            var no_sun = FalseIntervals(sun);  // will come back sorted
            var no_comm = FalseIntervals(comm);

            no_sun = CoalesceIntervals(no_sun);
            no_comm = CoalesceIntervals(no_comm);

            no_sun = no_sun.Where(i => !(i.Start > last_inside || i.Stop < first_inside)).ToList();
            no_comm = no_comm.Where(i => !(i.Start > last_inside || i.Stop < first_inside)).ToList();

            if (no_sun.Count == 0 || no_comm.Count == 0)
                return 1000;

            // For each no_comm, find the longest intersecting no_sun.  Call that hours;
            // find the minimum value of that, i.e., it's interesting if it has a short, longest darkness period during any cycle.
            var hours = no_comm.Min(nc =>
            {
                var h = LengthOfLongestNoSunInterval(nc, no_sun);
                return h == -1 ? 1000 : h;
            });
            return hours;
        }

        public static List<int> MaxShadowWithoutComm2(List<int> minimal_earth_elevation_indices, List<bool> has_sun, List<bool> has_comm)
        {
            var r = new List<int>();
            for (var i = 0; i < minimal_earth_elevation_indices.Count; i++)
            {
                var index = minimal_earth_elevation_indices[i];
                var first_no_comm = FirstTrue(has_comm, index, -1);
                if (first_no_comm < 0) first_no_comm = 0;
                var first = FirstTrue(has_sun, first_no_comm, -1);
                if (first < 0) first = 0;

                var last_no_comm = FirstTrue(has_comm, index, 1);
                if (last_no_comm < 0) last_no_comm = has_comm.Count - 1;
                var last = FirstTrue(has_sun, last_no_comm, 1);
                if (last < 0) last = has_sun.Count - 1;

                var shadow = CountShadow(has_sun, first, last);
                r.Add(shadow);
            }
            return r;

            int CountShadow(List<bool> has_sun1, int f, int l)
            {
                const int recharge_count = 12;  // 1 day
                var count_sun = 0;
                var count_shadow = 0;
                var max_shadow = 0;
                for (var i = f; i <= l; i++)
                {
                    if (has_sun1[i])
                    {
                        if (++count_sun >= recharge_count)
                        {
                            max_shadow = Math.Max(max_shadow, count_shadow);
                            count_shadow = 0;
                        }
                    }
                    else
                    {
                        count_shadow++;
                        count_sun = 0;
                    }
                }
                max_shadow = Math.Max(max_shadow, count_shadow);
                return max_shadow;
            }
        }

        public static List<Interval> FalseIntervals(List<bool> booleans)
        {
            var intervals = new List<Interval>();
            var start_of_false = FirstFalse(booleans, 0, 1);
            while (start_of_false > -1)
            {
                var end_of_false = FirstTrue(booleans, start_of_false, 1);
                if (end_of_false == -1)
                {
                    intervals.Add(new Interval { Start = start_of_false, Stop = booleans.Count - 1 });
                    break;
                }
                intervals.Add(new Interval { Start = start_of_false, Stop = end_of_false });
                start_of_false = FirstFalse(booleans, end_of_false, 1);
            }
            return intervals;
        }

        public static int FirstTrue(List<bool> booleans, int start, int step)
        {
            for (var index = start; index >= 0 && index < booleans.Count; index += step)
                if (booleans[index])
                    return index;
            return -1;
        }

        public static int FirstFalse(List<bool> l, int start, int step)
        {
            for (var index = start; index >= 0 && index < l.Count; index += step)
                if (!l[index])
                    return index;
            return -1;
        }

        List<Interval> CoalesceIntervals(List<Interval> false_intervals)
        {
            if (false_intervals.Count < 1) return false_intervals;
            const int gap_threshold = 12;
            var r = new List<Interval>();
            var interval = false_intervals[0];
            for (var i=1;i< false_intervals.Count; i++)
            {
                if (false_intervals[i].Start - interval.Stop > gap_threshold)
                {
                    r.Add(interval);
                    interval = false_intervals[i];
                }
                else
                {
                    interval.Stop = false_intervals[i].Stop;
                }
            }
            r.Add(interval);
            return r;
        }

        int LengthOfLongestNoSunInterval(Interval no_comm_interval, List<Interval> no_sun) => no_sun.Max(ns => ns.Intersects(no_comm_interval) ? ns.Length : -1);

        void WriteSafeHavenOverlay(float[,] ary, string path)
        {
            var map = new SafeHavenColorTable();
            Func<float, byte> func = v => (byte)Math.Min(255f,v);  // color entries map to hours
            using (var bmp = BitmapHelper.ToFormat8bppIndexed(ary, func))
            {
                map.LoadPalette(bmp);
                bmp.Save(LandingSiteDataset.AppendToFilename(path, "", ".png"), System.Drawing.Imaging.ImageFormat.Png);
                using (var bmp2 = map.AddLegend(bmp))
                    bmp2.Save(LandingSiteDataset.AppendToFilename(path, "_legend", ".png"), System.Drawing.Imaging.ImageFormat.Png);
            }
        }

        TerrainPatch GetPatchWithHorizons(Point id, int observer_height)
        {
            var patch = TerrainPatch.FromId(id, observer_height);
            if (File.Exists(patch.Path))
            {
                patch.Load();
                return patch;
            }
            lock (this)
            {
                var gpuProcessor = ViperEnvironment.Processor as GPUHorizons;
                var queue = new List<TerrainPatch>() { patch };
                gpuProcessor.RunQueue(queue, writeHorizons: WriteHorizons, unloadHorizons: false);
                return patch;
            }
        }

        public struct Interval
        {
            public int Start;
            public int Stop;

            public override string ToString() => $"[{Start},{Stop},{Hours}]";
            public int Length => Stop - Start;
            public float Hours => 2f * Length;
            public bool Intersects(Interval other) => !(Stop < other.Start || other.Stop < Start);
        }

    }
}
