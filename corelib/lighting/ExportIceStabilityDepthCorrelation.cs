using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using CsvHelper;
using CsvHelper.Configuration;
using viper.corelib.math;
using viper.corelib.patch;
using viper.corelib.terrain;
using viper.corelib.utilities;

namespace viper.corelib.lighting
{
    public class ExportIceStabilityDepthCorrelation
    {
        public static string LightingCacheFilename = "sun-contrib-cache.bin"; //"lighting-cache.bin";
        public static string IceStabilityDepthCorrelationFilename = "ice-stability-depth-correlation.csv";

        private Dictionary<(int,int,int), TerrainPatch> _patches = ViperEnvironment.PatchCache; // new Dictionary<Point, TerrainPatch>();
        DateTime _start = new DateTime(2018, 1, 1);
        DateTime _stop = new DateTime(2018 + 19, 1, 1);
        TimeSpan _step = TimeSpan.FromHours(2);

        public static string IceDepthCsvPath = @"C:\RP\maps\inputs_20m\HermiteA\Hermite\Other\HermiteA6000_Siegler.txt";

        Dictionary<Point, bool> _precache;  // Id -> bool, Cache answer about whether the point has horizons already calculated
        Dictionary<Point, TerrainPatch> _cache;

        public List<TerrainPatch> ExportTemperatureProfiles()
        {
            _precache = new Dictionary<Point, bool>();
            _cache = new Dictionary<Point, TerrainPatch>();

            var sieglers = ReadSieglerProduct2().ToList();
            Console.WriteLine($"There are {sieglers.Count} lines in the siegler product.");

            // Filter for sieglers near 12420,17549
            var center = new Point( 12420, 17549);
            sieglers = sieglers.Where(s => Distance(s.Pixel, center) < 30d).ToList();
            Console.WriteLine($"There are {sieglers.Count} pixels close to the test center");

            //var sieglers2 = sieglers.Where(s => s.Distance < d).ToList();
            //Console.WriteLine($"filtered by distance {d} count={sieglers2.Count} fraction={sieglers2.Count / (double)sieglers.Count}");
            //sieglers = sieglers.Where(s=>s.Depth<=1f && s.Depth>=0f).ToList();

            var group1 = sieglers.GroupBy(s => DepthEquivalenceClass(s.Depth)).ToList();
            var sieglers2 = new List<Siegler>();
            foreach (var g in group1)
                sieglers2.AddRange(g.Take(10));
            Console.WriteLine($"bunched by percent [0-1] count={sieglers2.Count} fraction={sieglers2.Count / (double)sieglers.Count}");

            var groups = sieglers2.GroupBy(s => s.Pixel).ToList();
            Console.WriteLine($"There are {groups.Count} grouped pixels within the map");
            var consolidated = groups.Select(g => new ConsolidatedLine(g.ToList())).ToList();

            consolidated = consolidated.Take(100).ToList();  // smaller set for now

            var filtered = consolidated.Where(HorizonsHaveBeenCalculated).ToList();
            Console.WriteLine($"There are {filtered.Count} pixels with horizons already calculated");
            Console.WriteLine($"There are {_cache.Count} tiles involved.");

            Console.WriteLine($"Beginning to load patches");
            Parallel.ForEach(_cache.Values, patch => GetPatch(patch.Id));
            Console.WriteLine($"Finished loading patches");

            var expected_count = consolidated.Count;
            var count = 0;
            var stopwatch = new Stopwatch();
            var starttime = DateTime.Now;
            stopwatch.Start();

            var queue = new BlockingCollection<ConsolidatedLine>(10000);
            Task.Run(() =>
            {
                Parallel.ForEach(consolidated, new ParallelOptions { MaxDegreeOfParallelism = Math.Max(Environment.ProcessorCount - 1, 1) }, c =>
                {
                    var r = DoCalculation(c);
                    if (r)
                    {
                        queue.Add(c);
                        Console.WriteLine(@"[{0:D4},{1:D4}].Lighting.Max()={2:F8}  Depth={3:F8}", c.Pixel.X, c.Pixel.Y, c.Lighting.Max(), c.Depth);
                    }
                    else
                        Console.WriteLine($"Calculation failed for [{c.Pixel.X},{c.Pixel.Y}]");
                });
                queue.CompleteAdding();
            });

            using (var fs = File.Create(LightingCacheFilename, 4096, FileOptions.None))
            using (var bw = new BinaryWriter(fs))
                foreach (var c in queue.GetConsumingEnumerable())
                {
                    c.Write(bw);
                    Math.DivRem(count++, 100, out int rem);
                    if (rem == 0)
                    {
                        var ticks_per_solution = stopwatch.Elapsed.Ticks / count;
                        var timespan_per_solution = new TimeSpan(ticks_per_solution);
                        var estimated_elapsed_ticks = ticks_per_solution * expected_count;
                        var estimated_stop = starttime.AddTicks(estimated_elapsed_ticks);
                        Console.WriteLine($"Processed {count} pixels.  {timespan_per_solution} per solution.  Est stop={estimated_stop}");
                    }
                }

            stopwatch.Stop();
            Console.WriteLine($"The lighting calculations took {stopwatch.Elapsed}");

            return _cache.Values.ToList();
        }

        IEnumerable<Siegler> ReadSieglerProduct2()
        {
            using (var tr = new StreamReader("ice-stability-depth-export.csv"))  // export of ice stability depth image
            {
                string file_line;
                tr.ReadLine();
                while ((file_line = tr.ReadLine())!=null)
                {
                    var split = file_line.Split(',');
                    var lat = double.Parse(split[0]);
                    var lon = double.Parse(split[1]);
                    var depth = float.Parse(split[2]);
                    var lat_rad = lat * Math.PI / 180d;
                    var lon_rad = lon * Math.PI / 180d;
                    InMemoryInt16Terrain.GetLineSample(lat_rad, lon_rad, out int line, out int sample);
                    //InMemoryTerrainManager.GetLatLonDegrees(line, sample, out double lat1, out double lon1);
                    if (line < 0 || sample < 0 || line >= InMemoryInt16Terrain.Samples || sample >= InMemoryInt16Terrain.Samples)
                        continue;
                    var pixel = new Point(sample, line);
                    yield return new Siegler { Lat = lat, Lon = lon, Depth = depth, Pixel = pixel };
                }
            }
        }

        IEnumerable<Siegler> ReadSieglerProduct()
        {
            var config = new Configuration { AllowComments = true, HasHeaderRecord = false, IgnoreBlankLines = true, Delimiter = "\t" };
            using (var tr = new StreamReader(IceDepthCsvPath))
            using (var csv = new CsvReader(tr, config))
            {
                while (csv.Read())
                {
                    // Assume lat,lon are at the center of the triangle
                    var depth = csv.GetField<float>(3);
                    var lon = -csv.GetField<float>(7);  // Note the - sign
                    var lat = csv.GetField<float>(8);
                    var x = csv.GetField<double>(4);
                    var y = csv.GetField<double>(5);
                    var z = csv.GetField<double>(6);

                    var lat1 = lat * Math.PI / 180d;
                    var lon1 = lon * Math.PI / 180d;

                    //var lon2 = Math.Atan2(y, x);
                    //if (lon2 < -Math.PI) lon2 += 2d * Math.PI;
                    //var d2 = Math.Sqrt(x * x + y * y);
                    //var lat2 = Math.Atan2(z, d2);

                    //lat2 = lat * Math.PI/180d;
                    //lon2 = lon * Math.PI / 180d;

                    InMemoryInt16Terrain.GetLineSample(lat1, lon1, out int line, out int sample);
                    if (line < 0 || sample < 0 || line >= InMemoryInt16Terrain.Samples || sample >= InMemoryInt16Terrain.Samples)
                        continue;

                    // Distance from tri center to one corner
                    var delta = Vector3d.Subtract(new Vector3d(x, y, z), TerrainPatch.MoonRadius * LatLonToPoint(lat1, lon1)).Length;
                    var pixel = new Point(sample, line);
                    yield return new Siegler { Lat = lat1, Lon = lon1, Depth = depth, Pixel = pixel, Distance = delta };
                }
            }
        }

        int DepthEquivalenceClass(float d)
        {
            const float maxdepth = 2.5f;
            if (d < 0f) d = 0f;
            var frac = d / maxdepth;
            return (int)(frac * 10);  // 10 classes
        }

        double Distance(Point a, Point b) => Math.Sqrt((a.X - b.X) * (a.X - b.X) + (a.Y - b.Y) * (a.Y - b.Y));

        // Cache by id
        private bool HorizonsHaveBeenCalculated(ConsolidatedLine consolidated)
        {
            var pixel = consolidated.Pixel;
            var id = TerrainPatch.LineSampleToId(pixel.Y, pixel.X);
            lock (_precache)
            {
                if (_precache.TryGetValue(id, out bool hasHorizons))
                    return hasHorizons;
                var patch = TerrainPatch.FromId(id);
                var horizons_path = patch.Path;
                var result = File.Exists(horizons_path);
                if (result)
                    _cache.Add(id, patch);
                _precache.Add(id, result);
                return result;
            }
        }

        Vector3d LatLonToPoint(double lat_rad, double lon_rad)
        {
            var z = Math.Sin(lat_rad);
            var d = Math.Cos(lat_rad);
            return new Vector3d(Math.Cos(lon_rad) * d, Math.Sin(lon_rad) * d, z);
        }

        private bool DoCalculation(ConsolidatedLine c)
        {
            var normal = GetSurfaceNormal(c.Pixel.Y, c.Pixel.X);
            c.NormalX = normal.X;
            c.NormalY = normal.Y;
            c.NormalZ = normal.Z;
            var steps = (_stop.Ticks - _start.Ticks) / _step.Ticks;
            var result = new List<float>((int)steps + 10);
            if (true)
            {
                if (!MapOverSunContribution(c.Pixel.Y, c.Pixel.X, _start, _stop, _step, (time, sun_contrib) => result.Add(sun_contrib)))
                    return false;
            }
            else
            {
                if (!MapOverLightCurve(c.Pixel.Y, c.Pixel.X, _start, _stop, _step, (time, light) => result.Add(light)))
                    return false;
            }
            c.Lighting = result.ToArray();
            return true;
        }

        public TerrainPatch GetPatch(Point id, int observer = 0)
        {
            lock (this)
            {
                if (_patches.TryGetValue((id.X, id.Y, observer), out TerrainPatch patch))
                    return patch;
                patch = TerrainPatch.FromId(id);
                var horizons_path = patch.Path;
                if (!File.Exists(horizons_path))
                    return null;
                patch.Read(horizons_path);
                patch.FillMatrices(ViperEnvironment.Terrain);
                _patches.Add((id.X, id.Y, observer), patch);
                return patch;
            }
        }

        /// <summary>
        /// Map over the fraction of the sun's disk visible
        /// </summary>
        /// <param name="line"></param>
        /// <param name="sample"></param>
        /// <param name="start"></param>
        /// <param name="stop"></param>
        /// <param name="step"></param>
        /// <param name="action"></param>
        /// <returns></returns>
        public bool MapOverLightCurve(int line, int sample, DateTime start, DateTime stop, TimeSpan step, Action<DateTime, float> action)
        {
            var id = TerrainPatch.LineSampleToId(line, sample);
            var patch = GetPatch(id);
            if (patch == null)
                return false;

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
                var sunfrac = horizon.SunFraction2(azimuth_deg, elevation_deg);

                //if (line == 17522 && sample == 12388)
                //{
                //    Console.WriteLine($"{line},{sample}, time={time} elevation_deg={elevation_deg} sun={sunfrac}");
                //}

                action(time, sunfrac);
            }
            return true;
        }

        /// <summary>
        /// Map over the fraction of the sun's disk visible times Cos(theta) - theta is the angle between the sun and normal
        /// </summary>
        /// <param name="line"></param>
        /// <param name="sample"></param>
        /// <param name="start"></param>
        /// <param name="stop"></param>
        /// <param name="step"></param>
        /// <param name="action"></param>
        /// <returns></returns>
        public bool MapOverSunContribution(int line, int sample, DateTime start, DateTime stop, TimeSpan step, Action<DateTime, float> action)
        {
            var id = TerrainPatch.LineSampleToId(line, sample);
            var patch = GetPatch(id);
            if (patch == null)
                return false;

            var cross = GetSurfaceNormal(line, sample);

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
                var sunfrac = horizon.SunFraction2(azimuth_deg, elevation_deg);

                //if (line == 17522 && sample == 12388)
                //{
                //    Console.WriteLine($"{line},{sample}, time={time} elevation_deg={elevation_deg} sun={sunfrac}");
                //}

                sunvec.NormalizeFast();  // Normalize before calculating surface temp
                var dot = Vector3d.Dot(sunvec, cross);
                var sun_contrib = dot < 0d ? 0d : sunfrac * dot;
                //if (sun_contrib > 0)
                //    Console.WriteLine(sun_contrib);
                action(time, (float) sun_contrib);
            }
            return true;
        }

        /// <summary>
        /// Map over estimated surface temperature
        /// </summary>
        /// <param name="line"></param>
        /// <param name="sample"></param>
        /// <param name="start"></param>
        /// <param name="stop"></param>
        /// <param name="step"></param>
        /// <param name="action"></param>
        /// <returns></returns>
        public bool MapOverTemperatureCurve(int line, int sample, DateTime start, DateTime stop, TimeSpan step, Action<DateTime, float> action)
        {
            var id = TerrainPatch.LineSampleToId(line, sample);
            var patch = GetPatch(id);
            if (patch == null)
                return false;

            var cross = GetSurfaceNormal(line, sample);

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
                var sunfrac = horizon.SunFraction2(azimuth_deg, elevation_deg);

                //if (line == 17522 && sample == 12388)
                //{
                //    Console.WriteLine($"{line},{sample}, time={time} elevation_deg={elevation_deg} sun={sunfrac}");
                //}

                sunvec.NormalizeFast();  // Normalize before calculating surface temp
                var tmp = SurfaceTemperature(sunvec, cross, sunfrac);
                //Console.WriteLine(tmp);
                action(time, tmp);
            }
            return true;
        }

        Vector3d GetSurfaceNormal(int line, int sample)
        {
            // Surface Normal
            var terrain = InMemoryInt16Terrain.Singleton;

            var span = 1f;                // step in pixels
            var vectorStep = span * 20d;  // step in meters
            var p1p2 = new Vector3d(vectorStep, 0d, 0d);
            var p1p3 = new Vector3d(0d, vectorStep, 0d);
            var cross = new Vector3d();

            var h1 = terrain.InterpolatedElevation(line, sample);
            var h2 = terrain.InterpolatedElevation(line + span, sample);
            var h3 = terrain.InterpolatedElevation(line, sample - span);

            p1p2.Z = h2 - h1;   // Z in meters
            p1p3.Z = h3 - h1;

            Vector3d.Cross(ref p1p2, ref p1p3, out cross);
            cross.NormalizeFast();
            return cross;
        }

        float SurfaceTemperature(Vector3d sunvec, Vector3d normal, float sunfrac)
        {
            const float albedo = 0.15f;
            const float emissivity = 0.975f;
            const float insolation = 1336f; //  W/m^2
            const float stefan_boltzmann = 5.670367e-8f;  // W / m^2 / K^4

            var dot = Vector3d.Dot(sunvec, normal);

            //var angle = (float)(Math.Acos(dot) * 180d / Math.PI);
            //return angle;
            //if (sunfrac == 1f && Math.Abs(angle - 88f) < 0.1f)
            //    Console.WriteLine("here");

            //if (dot < 0d)
            //    if (sunfrac > 0f)
            //        Console.WriteLine("here");

            if (dot < 0d)
                return 40f;
            var temp = (float)Math.Pow((insolation * (1f - albedo) * dot * sunfrac) / (emissivity * stefan_boltzmann), 0.25d);
            if (temp < 40f)
                return 40f;

            return temp;
        }

        public class Siegler
        {
            public double Lat;
            public double Lon;
            public double Distance;
            public Point Pixel;
            public float Depth;
        }

        public class OutputLine : Siegler
        {
            public float Estimate;  // originally longest cumulative time
            public OutputLine(Siegler s)
            {
                Lat = s.Lat;
                Lon = s.Lon;
                Depth = s.Depth;
            }
        }

        public class ConsolidatedLine
        {
            public const UInt32 Marker = 0xAAAAAAAA;
            public Point Pixel;
            public float Depth;
            public float StDev;
            public double NormalX;
            public double NormalY;
            public double NormalZ;
            public float[] Lighting;

            public int ByteLength => 48 + 4 * (Lighting != null ? Lighting.Length : 0);

            public ConsolidatedLine() { }
            public ConsolidatedLine(List<Siegler> lst)
            {
                Pixel = lst[0].Pixel;
                Depth = lst.Average(s => s.Depth);
                StDev = (float)Math.Sqrt(lst.Sum(s => Math.Pow(s.Depth - Depth, 2)) / (lst.Count - 1));
            }

            public void Write(BinaryWriter bw)
            {
                Debug.Assert(Lighting != null);
                bw.Write(Marker);
                bw.Write(Pixel.X);
                bw.Write(Pixel.Y);
                bw.Write(Depth);
                bw.Write(StDev);
                bw.Write(NormalX);
                bw.Write(NormalY);
                bw.Write(NormalZ);
                bw.Write(Lighting.Length);
                for (var i = 0; i < Lighting.Length; i++)
                    bw.Write(Lighting[i]);
            }

            public ConsolidatedLine Read(BinaryReader r, int expected_len = 0)
            {
                if (r.ReadUInt32() != Marker)
                    return null;
                Pixel = new Point(r.ReadInt32(), r.ReadInt32());
                Depth = r.ReadSingle();
                StDev = r.ReadSingle();
                NormalX = r.ReadDouble();
                NormalY = r.ReadDouble();
                NormalZ = r.ReadDouble();
                var buflen = r.ReadInt32();
                if (expected_len > 0 && expected_len != buflen)
                    return null;
                if (Lighting == null || Lighting.Length != buflen)
                    Lighting = new float[buflen];
                var lighting = Lighting;
                for (var i = 0; i < buflen; i++)
                    lighting[i] = r.ReadSingle();
                return this;
            }

            public static ConsolidatedLine ReadFrom(BinaryReader r, int expected_len = 0) => (new ConsolidatedLine()).Read(r, expected_len);
        }

        internal void WriteIceDepthCorrelation()
        {
            using (var tw = new StreamWriter(IceStabilityDepthCorrelationFilename))
            using (var fs = File.OpenRead(LightingCacheFilename))
            using (var br = new BinaryReader(fs))
            {
                tw.WriteLine("x,y,siegler_depth,max_temp");
                try
                {
                    while (true)
                    {
                        var l = ConsolidatedLine.ReadFrom(br);
                        tw.WriteLine($"{l.Pixel.X},{l.Pixel.Y},{l.Depth},{l.Lighting.BoxFilter(3).Max()},{l.Lighting.BoxFilter(6).Max()},{l.Lighting.BoxFilter(9).Max()},{l.Lighting.BoxFilter(12).Max()}");                        
                    }
                } catch (Exception) { }
            }
        }
    }

    public static class Extensions
    {
        public static IEnumerable<float> BoxFilter(this IEnumerable<float> stream, int width)
        {
            var buf = new float[width];
            var ptr = width-1;
            var sum = 0f;
            foreach (var v in stream)
            {
                sum = sum + v - buf[ptr];
                buf[ptr++] = v;
                if (ptr >= width)
                    ptr = 0;
                yield return sum / width;  // doesn't clip the first few responses, so they're low
            }
        }
    }
}
