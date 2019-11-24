using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using viper.corelib.horizon;
using viper.corelib.lighting;
using viper.corelib.math;
using viper.corelib.patch;
using viper.corelib.spice;
using viper.corelib.terrain;
using viper.corelib.utilities;

namespace calculator
{
    public class LightingCalculator
    {
        public const int Observer = 0;
        public const float EarthMultipathThreshold = 2f;

        public Rectangle IdRectangle { get; set; } = new Rectangle(147, 72, 1, 1);
        public bool WriteHorizons { get; set; } = false;
        public bool CalculatePSRPatch { get; set; } = false;
        public bool CalculateSafeHavenPatch { get; set; } = false;
        public bool CalculateAverageSunPatch { get; set; } = false;
        public bool CalculateAverageEarthPatch { get; set; } = false;
        public bool Verbose { get; set; } = false;
        public string MapRoot { get; set; } = "./sp/";

        public DateTime SafeHavenStart = new DateTime(2022, 12, 1);
        public DateTime SafeHavenStop = new DateTime(2023, 4, 1);
        public TimeSpan SafeHavenStep = TimeSpan.FromHours(2);

        protected List<Vector3d> sunvecs = null;
        protected List<Vector3d> earthvecs = null;

        protected void Setup()
        {
            ViperEnvironment.SetMapRoot(MapRoot);

            WriteHorizons = AppConfiguration.GetBoolean("WriteHorizons", WriteHorizons);
            CalculatePSRPatch = AppConfiguration.GetBoolean("CalculatePSRPatch", CalculatePSRPatch);
            CalculateSafeHavenPatch = AppConfiguration.GetBoolean("CalculateSafeHavenPatch", CalculateSafeHavenPatch);
            CalculateAverageSunPatch = AppConfiguration.GetBoolean("CalculateAverageSunPatch", CalculateAverageSunPatch);
            CalculateAverageEarthPatch = AppConfiguration.GetBoolean("CalculateAverageEarthPatch", CalculateAverageEarthPatch);
            Verbose = AppConfiguration.GetBoolean("Verbose", Verbose);

            var rectangle_string = AppConfiguration.Get("IdRectangle", null);
            if (rectangle_string != null)
            {
                try
                {
                    var strs = rectangle_string.Split(',');
                    if (strs.Length != 4) throw new Exception(@"Wrong number of ints in IdRectangle (should be 4).  Ignoring ...");
                    var ints = strs.Select(s => int.TryParse(s, out int r) ? null : (int?)r).ToList();
                    var index = ints.FindIndex(i => !i.HasValue);
                    if (index > -1) throw new Exception($"Failure to parse {ints[index].Value} in IdRectangle.  Ignoring ...");
                    IdRectangle = new Rectangle(ints[0].Value, ints[1].Value, ints[2].Value, ints[3].Value);
                } catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                }
            }

            // Initialization
            OSGeo.GDAL.Gdal.AllRegister();
            ViperEnvironment.Terrain = new InMemoryInt16Terrain();
            ViperEnvironment.Terrain.LoadSouth();
            ViperEnvironment.Processor = new GPUHorizons { Terrain = ViperEnvironment.Terrain };

            SpiceManager.GetSingleton();
        }

        protected void PrintSetup()
        {
            Console.WriteLine($"             WriteHorizons = {WriteHorizons}");
            Console.WriteLine($"         CalculatePSRPatch = {CalculatePSRPatch}");
            Console.WriteLine($"   CalculateSafeHavenPatch = {CalculateSafeHavenPatch}");
            Console.WriteLine($"  CalculateAverageSunPatch = {CalculateAverageSunPatch}");
            Console.WriteLine($"CalculateAverageEarthPatch = {CalculateAverageEarthPatch}");
            Console.WriteLine($"                   MapRoot = {ViperEnvironment.MapRoot}");
            Console.WriteLine($"               HorizonRoot = {ViperEnvironment.HorizonRoot}");
            Console.WriteLine($"              AvgPatchRoot = {ViperEnvironment.AvgPatchRoot}");
            Console.WriteLine($"              PSRPatchRoot = {ViperEnvironment.PSRPatchRoot}");
            Console.WriteLine($"            HavenPatchRoot = {ViperEnvironment.HavenPatchRoot}");
            Console.WriteLine($"               IdRectangle = {IdRectangle}");
        }

        public unsafe void Run()
        {
            Setup();
            if (Verbose)
                PrintSetup();

            if (WriteHorizons)
                Directory.CreateDirectory(ViperEnvironment.HorizonRoot);
            if (CalculatePSRPatch)
                Directory.CreateDirectory(ViperEnvironment.PSRPatchRoot);
            if (CalculateSafeHavenPatch)
                Directory.CreateDirectory(ViperEnvironment.HavenPatchRoot);
            if (CalculateAverageSunPatch)
                Directory.CreateDirectory(ViperEnvironment.AvgPatchRoot);
            if (CalculateAverageEarthPatch)
                Directory.CreateDirectory(ViperEnvironment.AvgPatchRoot);

            // Perform the actual work
            var queue = new BlockingCollection<TerrainPatch>(48);
            var producer = Task.Run(() =>
            {
                var patches = TerrainPatch.EnumerateIds(IdRectangle).Select(id => TerrainPatch.FromId(id, 0)).ToList();
                var tocalculate = new List<TerrainPatch>();
                var toload = new List<TerrainPatch>();
                foreach (var p in patches)
                    if (File.Exists(p.Path))
                        toload.Add(p);
                    else
                        tocalculate.Add(p);
                var stopwatch = Stopwatch.StartNew();
                if (toload.Count > 0)
                    Console.WriteLine("Starting to load {patches.Count} horizons");
                foreach (var p in toload)
                {
                    p.Load();
                    queue.Add(p);
                }
                if (tocalculate.Count > 0)
                    Console.WriteLine($"Starting to calculate horizons for {queue.Count} patches");
                var processor = ViperEnvironment.Processor as GPUHorizons;
                foreach (var patch in tocalculate)
                {
                    Console.WriteLine($"Generating or loading patch id={patch.Id}");
                    processor.ProcessNearField(patch);
                    patch.Matrices = null;      // Force the matrices to be regenerated so they're not left over from the near horizon
                    patch.FillPoints(ViperEnvironment.Terrain);  // Should be a no op
                    patch.FillMatricesRelativeToPoint(ViperEnvironment.Terrain, patch.Points[0][0]);
                    if (patch.Horizons == null)
                        patch.InitializeHorizons();
                    var patchDictionary = new Dictionary<int, TerrainPatch>();
                    var far_field = processor.GenerateFarField(patch.Id, patchDictionary);
                    processor.ProcessFarField(patch, far_field);// Here is the GPU Call
                    if (WriteHorizons)
                        patch.Write();
                    patch.Matrices = null;
                    queue.Add(patch);
                }
                Console.WriteLine($"Horizon generation or loading took {stopwatch.Elapsed}");
                queue.CompleteAdding();
            });

            // Consumer
            Parallel.ForEach(queue.GetConsumingEnumerable(), patch =>
            {
                Console.WriteLine($"Processing patch id={patch.Id}");
                patch.Matrices = null;
                patch.FillPointsAndMatrices(ViperEnvironment.Terrain);
                if (CalculatePSRPatch)
                    GeneratePSRPatch(patch);
                if (CalculateSafeHavenPatch)
                    GenerateSafeHavenPatch(patch);
                if (CalculateAverageSunPatch)
                    GenerateAverageSunPatch(patch);
                if (CalculateAverageEarthPatch)
                    GenerateAverageEarthPatch(patch);
            });

            producer.Wait();
        }

        unsafe void GeneratePSRPatch(TerrainPatch patch)
        {
            var path = TileLightingProductManager.PSRPatchPath(patch.Line, patch.Sample);
            if (File.Exists(path))
                return;

            const int azimuth_steps = 360 * 4;
            const float twopi = (float)(Math.PI * 2d);
            var sunvec = new Vector3d[azimuth_steps];
            var elevation = new float[azimuth_steps];
            Array.Fill(elevation, float.MinValue);
            var interval = StudyInterval.MazaricoCycle;
            var step = TimeSpan.FromHours(2);
            for (var time = interval.Start; time < interval.Stop; time += step)
            {
                var sun = ViperEnvironment.SunPosition(time);
                patch.GetAzEl(sun, 0, 0, out float sun_azimuth_rad, out float sun_elevation_rad);
                var azimuth = (int)(sun_azimuth_rad * azimuth_steps / twopi);
                if (sun_elevation_rad > elevation[azimuth])
                    sunvec[azimuth] = sun;
            }

            var bmp = new Bitmap(TerrainPatch.DefaultSize, TerrainPatch.DefaultSize, System.Drawing.Imaging.PixelFormat.Format8bppIndexed);
            var bmp_data = bmp.LockBits(new Rectangle(0, 0, TerrainPatch.DefaultSize, TerrainPatch.DefaultSize), ImageLockMode.WriteOnly, bmp.PixelFormat);
            for (var row = 0; row < TerrainPatch.DefaultSize; row++)
            {
                var rowptr = (byte*)(bmp_data.Scan0 + row * bmp_data.Stride);
                for (var col = 0; col < TerrainPatch.DefaultSize; col++)
                {
                    rowptr[col] = 1;
                    var horizon = patch.Horizons[row][col];
                    for (var azimuth = 0; azimuth < azimuth_steps; azimuth++)
                    {
                        var sun = sunvec[azimuth];
                        patch.GetAzEl(sun, col, row, out float sun_azimuth_rad, out float sun_elevation_rad);
                        var val = horizon.SunFraction2(sun_azimuth_rad * 180f / 3.141592653589f, sun_elevation_rad * 180f / 3.141592653589f);
                        if (val > 0f)
                            rowptr[col] = 0;
                    }
                }
            }
            bmp.UnlockBits(bmp_data);
            var palette = bmp.Palette;
            palette.Entries[0] = Color.FromArgb(0, 0, 0, 0);
            palette.Entries[1] = Color.Red;
            bmp.Palette = palette;
            bmp.Save(path, ImageFormat.Png);
        }

        void GenerateSafeHavenPatch(TerrainPatch patch)
        {
            const float SafeHavenSunThreshold = 0.5f;

            var outer_start = SafeHavenStart.AddDays(-28);
            var outer_stop = SafeHavenStop.AddDays(28);
            GetLowEarthTimes(outer_start, outer_stop, out List<DateTime> times_start_to_stop, out List<int> indices_of_minima_earth_elevation);

            for (var i = 0; i < indices_of_minima_earth_elevation.Count; i++)
                Console.WriteLine($"i={i} index={indices_of_minima_earth_elevation[i]} time={times_start_to_stop[indices_of_minima_earth_elevation[i]]}");
            var combined = new float[TerrainPatch.DefaultSize, TerrainPatch.DefaultSize];
            var month_images = Enumerable.Range(0, indices_of_minima_earth_elevation.Count).Select(i => new float[TerrainPatch.DefaultSize, TerrainPatch.DefaultSize]).ToList();

            var sun_vectors = times_start_to_stop.Select(time => CSpice.SunPosition(time)).ToList();
            var earth_vectors = times_start_to_stop.Select(time => CSpice.EarthPosition(time)).ToList();

            for (var row = 0; row < TerrainPatch.DefaultSize; row++)
            {
                for (var col = 0; col < TerrainPatch.DefaultSize; col++)
                {
                    var has_sun = patch.GetLightCurve(row, col, sun_vectors).Select(sun => sun >= SafeHavenSunThreshold).ToList();
                    var has_comm = patch.GetEarthOverHorizonCurve(row, col, earth_vectors).Select(rise => rise >= EarthMultipathThreshold).ToList();
                    var month_steps = MaxShadowWithoutComm2(indices_of_minima_earth_elevation, has_sun, has_comm);
                    combined[row, col] = 2f * month_steps.Max();
                    for (var i = 0; i < month_steps.Count; i++)
                        month_images[i][row, col] = 2f * month_steps[i];
                }
            }

            {
                var path = TileLightingProductManager.HavenPath(patch.Line, patch.Sample, SafeHavenStart,SafeHavenStop,SafeHavenStep, Observer);
                var path_combined = LandingSiteDataset.AppendToFilename(path, "_combined");
                GeotiffHelper.WriteArrayAsGeotiff(combined, patch.Bounds, path_combined);
            }

            for (var i = 0; i < month_images.Count; i++)
            {
                var time = times_start_to_stop[indices_of_minima_earth_elevation[i]];
                var path = TileLightingProductManager.HavenPath(patch.Line, patch.Sample, time, Observer);
                //var path1 = LandingSiteDataset.AppendToFilename(path, "_" + time.ToString("yyyy-MM-dd"));
                GeotiffHelper.WriteArrayAsGeotiff(month_images[i], patch.Bounds, path);
            }

            // Return all times from start to stop.  Also return the indexes of the times at which the earth elevation reaches a minima
            void GetLowEarthTimes(DateTime outer_start, DateTime outer_stop, out List<DateTime> all_times, out List<int> indices_of_minima_earth_elevation)
            {
                // Generate the time / elevation pairs
                all_times = new List<DateTime>();
                for (var time = outer_start; time <= outer_stop; time += SafeHavenStep)
                    all_times.Add(time);
                var earth_elevations = all_times.Select(t =>
                {
                    patch.GetAzEl(CSpice.EarthPosition(t), 0, 0, out float _, out float earth_elevation_rad);
                    return earth_elevation_rad;
                }).ToList();

                // Find minima
                indices_of_minima_earth_elevation = new List<int>();
                for (var i = 1; i < earth_elevations.Count - 1; i++)
                    if (SafeHavenStart <= all_times[i] && all_times[i] <= SafeHavenStop && earth_elevations[i - 1] > earth_elevations[i] && earth_elevations[i] < earth_elevations[i + 1])
                        indices_of_minima_earth_elevation.Add(i);
            }

            static List<int> MaxShadowWithoutComm2(List<int> minimal_earth_elevation_indices, List<bool> has_sun, List<bool> has_comm)
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

                static int CountShadow(List<bool> has_sun1, int f, int l)
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

            static int FirstTrue(List<bool> booleans, int start, int step)
            {
                for (var index = start; index >= 0 && index < booleans.Count; index += step)
                    if (booleans[index])
                        return index;
                return -1;
            }

            static int FirstFalse(List<bool> l, int start, int step)
            {
                for (var index = start; index >= 0 && index < l.Count; index += step)
                    if (!l[index])
                        return index;
                return -1;
            }
        }

        unsafe void GenerateAverageSunPatch(TerrainPatch patch)
        {
            var interval = StudyInterval.MazaricoCycle;
            var path = TileLightingProductManager.AverageSunPath(patch.Line, patch.Sample, interval.Start, interval.Stop, interval.Step, Observer);
            if (File.Exists(path))
                return;
            PrepareSunAndEarthVectors();
            using (var bmp = patch.GenerateAverageSun(sunvecs))
            {
                bmp.Save(path, ImageFormat.Png);
            }
        }

        void GenerateAverageEarthPatch(TerrainPatch patch)
        {
            var path = TileLightingProductManager.PSRPatchPath(patch.Line, patch.Sample);
            if (File.Exists(path))
                return;
            PrepareSunAndEarthVectors();
            using (var bmp = patch.GenerateAverageEarth(earthvecs, EarthMultipathThreshold))
            {
                bmp.Save(path, ImageFormat.Png);
            }
        }

        void PrepareSunAndEarthVectors()
        {
            if (sunvecs != null) return;
            sunvecs = new List<Vector3d>();
            earthvecs = new List<Vector3d>();
            var interval = StudyInterval.MazaricoCycle;
            for (var time = interval.Start; time < interval.Stop; time += interval.Step)
            {
                sunvecs.Add(ViperEnvironment.SunPosition(time));
                earthvecs.Add(ViperEnvironment.EarthPosition(time));
            }
        }
    }
}
