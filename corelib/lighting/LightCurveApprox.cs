using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using viper.corelib.math;
using viper.corelib.patch;
using viper.corelib.utilities;

namespace viper.corelib.lighting
{
    /// <summary>
    /// Calculate light curves, sun and earth elevations.
    /// </summary>
    public class LightCurveApprox
    {
        const float PI = 3.141592653589f;
        const float RadToDeg = 180f / PI;

        private Dictionary<(int,int,int), TerrainPatch> _patches = ViperEnvironment.PatchCache; // new Dictionary<Point, TerrainPatch>();

        /*
        public PointPairList GetLightCurve(int line, int sample, DateTime start, DateTime stop, TimeSpan step)
        {
            var result = new PointPairList();
            Action<DateTime, float> action = (time, sun) => result.Add(new PointPair(time.ToOADate(), sun));
            MapOverLightCurve(line, sample, start, stop, step, action);
            return result;
        }

        public PointPairList GetSunElevationCurve(int line, int sample, DateTime start, DateTime stop, TimeSpan step)
        {
            var result = new PointPairList();
            Action<DateTime, float> action = (time, sun) => result.Add(new PointPair(time.ToOADate(), sun));
            MapOverSunElevationCurve(line, sample, start, stop, step, action);
            return result;
        }

        public PointPairList GetSunHorizonElevationCurve(int line, int sample, DateTime start, DateTime stop, TimeSpan step)
        {
            var result = new PointPairList();
            Action<DateTime, float> action = (time, sun) => result.Add(new PointPair(time.ToOADate(), sun));
            MapOverSunHorizonElevationCurve(line, sample, start, stop, step, action);
            return result;
        }

        // earth elevation in degrees
        public PointPairList GetEarthElevationCurve(int line, int sample, DateTime start, DateTime stop, TimeSpan step)
        {
            var result = new PointPairList();
            Action<DateTime, float> action = (time, elevation) => result.Add(new PointPair(time.ToOADate(), elevation));
            MapOverEarthElevationCurve(line, sample, start, stop, step, action);
            return result;
        }

        public PointPairList GetEarthHorizonElevationCurve(int line, int sample, DateTime start, DateTime stop, TimeSpan step)
        {
            var result = new PointPairList();
            Action<DateTime, float> action = (time, elevation) => result.Add(new PointPair(time.ToOADate(), elevation));
            MapOverEarthHorizonElevationCurve(line, sample, start, stop, step, action);
            return result;
        }
        */

        public float ShortestNight(int line, int sample, DateTime start, DateTime stop, float minWorkingHoursOfSun = 48f, float sunThreshold = 0.8f)
        {
            var step = new TimeSpan(1, 0, 0);
            var intervals = new List<ShortNightInterval>();
            var validState = false;
            var inSun = false;
            DateTime lastTime = new DateTime();
            Action<DateTime, float> action = (time, sun) =>
            {
                var light = sun >= sunThreshold;
                if (!validState)
                {
                    validState = true;
                    inSun = light;
                    lastTime = time;
                }
                else
                {
                    if (inSun != light)
                    {
                        intervals.Add(new ShortNightInterval { InSun = inSun, Start = lastTime, Stop = time });
                        inSun = light;
                        lastTime = time;
                    }
                }
            };
            MapOverLightCurve(line, sample, start, stop, step, action);

            return 0;  // not finished

        }

        public TerrainPatch GetPatch(Point id, int observer = 0)
        {
            lock (this)
            {
                if (_patches.TryGetValue((id.X,id.Y,observer), out TerrainPatch patch))
                    return patch;
                patch = TerrainPatch.FromId(id, observer);
                var horizons_path = patch.Path;
                if (!File.Exists(horizons_path))
                    return null;
                patch.Read(horizons_path);
                patch.FillMatrices(ViperEnvironment.Terrain);
                _patches.Add((id.X, id.Y, observer), patch);
                return patch;
            }
        }

        public float CumulativeSun(int line, int sample, DateTime start, DateTime stop, TimeSpan step)
        {
            var sum = 0d;
            MapOverLightCurve(line, sample, start, stop, step, (time, sun) => sum += sun);
            var totalSunHours = sum * step.TotalHours;
            var maxSunHours = (stop - start).TotalHours;
            return (float)(totalSunHours / maxSunHours);
        }

        public void MapOverLightCurve(int line, int sample, DateTime start, DateTime stop, TimeSpan step, Action<DateTime, float> action)
        {
            var id = TerrainPatch.LineSampleToId(line, sample);
            var patch = GetPatch(id);
            if (patch == null)
                return;

            patch.FillPointsAndMatrices(terrain.InMemoryInt16Terrain.Singleton);

            var y_offset = line - patch.Line;
            var x_offset = sample - patch.Sample;
            var mat = patch.Matrices[y_offset][x_offset];
            var horizon = patch.Horizons[y_offset][x_offset];
            if (!horizon.IsDegrees)
                lock (horizon)
                    horizon.ConvertSlopeToDegrees();

            var cache = SunVectorCache.GetSingleton();
            var temp = new Vector3d();
            for (var time = start; time <= stop; time += step)
            {
                var sunvec = cache.SunPosition(time);
                TerrainPatch.Transform(ref sunvec, ref mat, ref temp);

                var sun_x = temp[0];
                var sun_y = temp[1];
                var sun_z = temp[2];
                var azimuth_rad = Math.Atan2(sun_y, sun_x) + Math.PI;  // [0,2PI]
                var azimuth_deg = (float)(azimuth_rad * 180d / Math.PI);

                var alen = Math.Sqrt(sun_x * sun_x + sun_y * sun_y);
                var slope = sun_z / alen;
                var elevation_deg = ((float)Math.Atan(slope)) * 180f / 3.141592653589f;
                var sun = horizon.SunFraction2(azimuth_deg, elevation_deg);

                //Console.WriteLine($"time={time} elevation_deg={elevation_deg} sun={sun}");

                action(time, sun);
            }
        }

        public void MapOverSunElevationCurve(int line, int sample, DateTime start, DateTime stop, TimeSpan step, Action<DateTime, float> action)
        {
            var id = TerrainPatch.LineSampleToId(line, sample);
            var patch = GetPatch(id);
            if (patch == null)
                return;

            patch.FillPointsAndMatrices(terrain.InMemoryInt16Terrain.Singleton);

            var y_offset = line - patch.Line;
            var x_offset = sample - patch.Sample;
            var mat = patch.Matrices[y_offset][x_offset];

            var cache = SunVectorCache.GetSingleton();
            var temp = new Vector3d();
            for (var time = start; time <= stop; time += step)
            {
                var sunvec = cache.SunPosition(time);
                TerrainPatch.Transform(ref sunvec, ref mat, ref temp);

                var sun_x = temp[0];
                var sun_y = temp[1];
                var sun_z = temp[2];
                var azimuth_rad = Math.Atan2(sun_y, sun_x) + Math.PI;  // [0,2PI]

                var alen = Math.Sqrt(sun_x * sun_x + sun_y * sun_y);
                var slope = sun_z / alen;
                var elevation_deg = ((float)Math.Atan(slope)) * 180f / 3.141592653589f;

                action(time, elevation_deg);
            }
        }

        public void MapOverSunHorizonElevationCurve(int line, int sample, DateTime start, DateTime stop, TimeSpan step, Action<DateTime, float> action)
        {
            var id = TerrainPatch.LineSampleToId(line, sample);
            var patch = GetPatch(id);
            if (patch == null)
                return;

            patch.FillPointsAndMatrices(terrain.InMemoryInt16Terrain.Singleton);

            var y_offset = line - patch.Line;
            var x_offset = sample - patch.Sample;
            var mat = patch.Matrices[y_offset][x_offset];
            var horizon = patch.Horizons[y_offset][x_offset];
            if (!horizon.IsDegrees)
                lock (horizon)
                    horizon.ConvertSlopeToDegrees();

            var cache = SunVectorCache.GetSingleton();
            var temp = new Vector3d();
            for (var time = start; time <= stop; time += step)
            {
                var sunvec = cache.SunPosition(time);
                TerrainPatch.Transform(ref sunvec, ref mat, ref temp);

                var sun_x = temp[0];
                var sun_y = temp[1];
                var sun_z = temp[2];
                var azimuth_rad = (float)Math.Atan2(sun_y, sun_x) + PI;  // [0,2PI]

                var alen = Math.Sqrt(sun_x * sun_x + sun_y * sun_y);
                var slope = sun_z / alen;
                var elevation_deg = (float)Math.Atan(slope) * RadToDeg;

                var over_horizon_deg = horizon.OverHorizon(azimuth_rad, elevation_deg);  // deg

                action(time, over_horizon_deg);
            }
        }

        public void MapOverEarthElevationCurve(int line, int sample, DateTime start, DateTime stop, TimeSpan step, Action<DateTime, float> action)
        {
            var id = TerrainPatch.LineSampleToId(line, sample);
            var patch = GetPatch(id);
            if (patch == null)
                return;

            patch.FillPointsAndMatrices(terrain.InMemoryInt16Terrain.Singleton);

            var y_offset = line - patch.Line;
            var x_offset = sample - patch.Sample;
            var mat = patch.Matrices[y_offset][x_offset];

            var cache = SunVectorCache.GetSingleton();
            var temp = new Vector3d();
            for (var time = start; time <= stop; time += step)
            {
                var earthvec = cache.EarthPosition(time);
                TerrainPatch.Transform(ref earthvec, ref mat, ref temp);

                var earth_x = temp[0];
                var earth_y = temp[1];
                var earth_z = temp[2];
                var azimuth_rad = Math.Atan2(earth_y, earth_x) + Math.PI;  // [0,2PI]

                var alen = Math.Sqrt(earth_x * earth_x + earth_y * earth_y);
                var slope = earth_z / alen;
                var elevation_deg = (float)Math.Atan(slope) * RadToDeg;

                action(time, elevation_deg);
            }
        }

        public void MapOverEarthHorizonElevationCurve(int line, int sample, DateTime start, DateTime stop, TimeSpan step, Action<DateTime, float> action)
        {
            var id = TerrainPatch.LineSampleToId(line, sample);
            var patch = GetPatch(id);
            if (patch == null)
                return;

            patch.FillPointsAndMatrices(terrain.InMemoryInt16Terrain.Singleton);

            var y_offset = line - patch.Line;
            var x_offset = sample - patch.Sample;
            var mat = patch.Matrices[y_offset][x_offset];
            var horizon = patch.Horizons[y_offset][x_offset];
            if (!horizon.IsDegrees)
                lock (horizon)
                    horizon.ConvertSlopeToDegrees();

            var cache = SunVectorCache.GetSingleton();
            var temp = new Vector3d();
            for (var time = start; time <= stop; time += step)
            {
                var earthvec = cache.EarthPosition(time);
                TerrainPatch.Transform(ref earthvec, ref mat, ref temp);

                var earth_x = temp[0];
                var earth_y = temp[1];
                var earth_z = temp[2];
                var azimuth_rad = Math.Atan2(earth_y, earth_x) + Math.PI;  // [0,2PI]
                var azimuth_deg = (float)(azimuth_rad * 180d / Math.PI);

                var alen = Math.Sqrt(earth_x * earth_x + earth_y * earth_y);
                var slope = earth_z / alen;
                var elevation_deg = (float)Math.Atan(slope) * RadToDeg;

                var over_horizon_deg = horizon.OverHorizon((float)azimuth_rad, elevation_deg);  // deg

                action(time, over_horizon_deg);
            }
        }


        public IEnumerable<(DateTime time, float sun)> LightCurveGenerator(int line, int sample, DateTime start, DateTime stop, TimeSpan step)
        {
            var id = TerrainPatch.LineSampleToId(line, sample);
            var patch = GetPatch(id);
            if (patch == null)
                yield break;

            var y_offset = line - patch.Line;
            var x_offset = sample - patch.Sample;
            var mat = patch.Matrices[y_offset][x_offset];
            var horizon = patch.Horizons[y_offset][x_offset];
            if (!horizon.IsDegrees)
                lock (horizon)
                    horizon.ConvertSlopeToDegrees();

            var cache = SunVectorCache.GetSingleton();
            var temp = new Vector3d();
            for (var time = start; time <= stop; time += step)
            {
                var sunvec = cache.SunPosition(time);
                TerrainPatch.Transform(ref sunvec, ref mat, ref temp);

                var sun_x = temp[0];
                var sun_y = temp[1];
                var sun_z = temp[2];
                var azimuth_rad = Math.Atan2(sun_y, sun_x) + Math.PI;  // [0,2PI]
                var azimuth_deg = (float)(azimuth_rad * 180d / Math.PI);

                var alen = Math.Sqrt(sun_x * sun_x + sun_y * sun_y);
                var slope = sun_z / alen;
                var elevation_deg = ((float)Math.Atan(slope)) * 180f / 3.141592653589f;
                var sun = horizon.SunFraction2(azimuth_deg, elevation_deg);

                yield return (time: time, sun: sun);
            }
        }

        public IEnumerable<(DateTime time, bool inSun)> InSunGenerator(int line, int sample, DateTime start, DateTime stop, TimeSpan step, float threshold = 0.8f) =>
            LightCurveGenerator(line, sample, start, stop, step).Select(t => (time: t.time, inSun: t.sun >= threshold));

        public long[,] GenerateLongestNightDataset(Rectangle rect, DateTime start, DateTime stop, TimeSpan step, float threshold)
        {
            var width = rect.Width;
            var height = rect.Height;
            var pixel_data = new long[height, width];
            var options = new ParallelOptions { MaxDegreeOfParallelism = ViperEnvironment.MaxDegreeOfParallelism };
            Parallel.For(0, height, options, y =>
            {
                for (var x = 0; x < width; x++)
                {
                    var during_day = false;
                    var dark_start = start;
                    var longest_night_so_far = new TimeSpan();
                    foreach (var (date, in_sun) in InSunGenerator(y + rect.Top, x + rect.Left, start, stop, step, threshold))
                    {
                        if (during_day)
                        {
                            if (!in_sun)  // going into darkness
                                dark_start = date;
                        }
                        else
                        {
                            if (in_sun)
                            {
                                var nightspan = date - dark_start;
                                if (nightspan > longest_night_so_far)
                                    longest_night_so_far = nightspan;
                            }
                        }
                        during_day = in_sun;
                    }
                    var ticks = longest_night_so_far.Ticks;
                    pixel_data[y, x] = ticks;
                }
            });
            return pixel_data;
        }
    }

    public class ShortNightInterval
    {
        public bool InSun;
        public DateTime Start;
        public DateTime Stop;
    }

    public class TimeImageDataset
    {
        public Rectangle Bounds;
        public long[,] Data;

        public TimeImageDataset() { }
        public TimeImageDataset(string filename) => Read(filename);

        public void Read(string filename)
        {
            using (var reader = new BinaryReader(File.Open(filename, FileMode.Open)))
            {
                Bounds = new Rectangle(reader.ReadInt32(), reader.ReadInt32(), reader.ReadInt32(), reader.ReadInt32());
                var height = Bounds.Height;
                var width = Bounds.Width;
                var data = new long[height, width];
                for (var y = 0; y < height; y++)
                    for (var x = 0; x < width; x++)
                        data[y, x] = reader.ReadInt64();
                Data = data;
            }
        }

        public void Write(string filename)
        {
            var data = Data;
            var height = data.GetLength(0);
            var width = data.GetLength(1);
            Debug.Assert(Bounds.Width == width && Bounds.Height == height);
            using (var writer = new BinaryWriter(File.Open(filename, FileMode.Create)))
            {
                writer.Write(Bounds.Left);
                writer.Write(Bounds.Top);
                writer.Write(Bounds.Width);
                writer.Write(Bounds.Height);
                for (var y = 0; y < height; y++)
                    for (var x = 0; x < width; x++)
                        writer.Write(data[y, x]);
            }
        }

        public IEnumerable<long> Enumerate() => Data.Enumerate();
    }
}