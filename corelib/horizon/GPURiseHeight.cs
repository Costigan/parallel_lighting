using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using ILGPU;
using ILGPU.Runtime;
using System.Drawing;
using viper.corelib.terrain;
using viper.corelib.math;
using viper.corelib.patch;
using viper.corelib.spice;
using viper.corelib.lighting;
using viper.corelib.utilities;

namespace viper.corelib.horizon
{
    public class GPURiseHeight : CPUHorizons
    {
        public const float MoonRadius = 1737.4f;
        //public const float MetersToMoonRadiusUnits = 1f / MoonRadiusUnitsToMeters;
        //public const float MoonRadiusUnitsToMeters = 1737400f;
        //public const float PixelsToMoonRadiusUnits = 20f / MoonRadiusUnitsToMeters;
        public const float PixelsToKM = 0.020f;

        public const float S0 = 15199.5f;             // PDS SAMPLE_PROJECTION_OFFSET
        public const float L0 = 15199.5f;

        public const float Scale = 20f / 1000f;
        public const float LonP = 0f;

        //public static float LonFactor = 1f;  
        //public static float LatP;
        public const float LonFactor = 1f;           // appropriate for south
        public const float LatP = -GPUMath.PI / 2f;

        public const int TerrainSize = InMemoryInt16Terrain.Samples;

        public AcceleratorType DefaultAcceleratorType = AcceleratorType.Cuda;
        public Rectangle? BoundingBox = null;

        List<Vector3d> _sunvecs = null;

        public const int cpu_slope_size = 360 * 4;
        float[] cpu_slope = new float[cpu_slope_size];
        float[] cpu_rise = new float[TerrainPatch.DefaultSize * TerrainPatch.DefaultSize];

        const int small_range_step_count = 1024;
        float[] cpu_range;
        int[] cpu_range_limit;

        private static readonly Index2 Dimension = new Index2(TerrainPatch.DefaultSize, TerrainPatch.DefaultSize);

        public void RunQueue(List<TerrainPatch> queue, Action<List<TerrainPatch>> queue_saver = null)
        {
            Debug.Assert(LonFactor == (float)InMemoryInt16Terrain.LonFactor);
            Debug.Assert(LatP == (float)InMemoryInt16Terrain.LatP);
            Terrain = ViperEnvironment.Terrain;

            var stopwatchOuter = new Stopwatch();
            stopwatchOuter.Start();
            if (StatusReceiver != null) StatusReceiver.IsRunning = true;
            Debug.Assert(ViperEnvironment.Terrain != null);
            if (queue == null)
                queue = ReadShadowCalculationQueue();

            var run_time = PatchLight.EstimateCalculationTime(queue);
            Console.WriteLine($"Queue length          = {queue.Count}");
            Console.WriteLine($"Estimated run time    = {run_time}");
            Console.WriteLine($"Estimated finish time = {DateTime.Now + run_time}");

            var patchDictionary = new Dictionary<int, TerrainPatch>();
            BoundingBox = null;

            var interval = StudyInterval.MetonicCycle;
            interval.Step = TimeSpan.FromMinutes(20);
            PrepareSlopeArray(interval);
            WriteSlopeArray();
            FillRangeArray();
            cpu_range_limit = new int[cpu_slope.Length];

            using (var context = new Context())
            {
                if (!GetAccelerator(AcceleratorType.Cuda, out AcceleratorId aid))  // DefaultAcceleratorType
                    return;
                using (var accelerator = Accelerator.Create(context, aid))
                {
                    using (var gpu_terrain = accelerator.Allocate<short>(ViperEnvironment.Terrain.Data.Length))
                    using (var gpu_range = accelerator.Allocate<float>(cpu_range.Length))
                    using (var gpu_range_limit = accelerator.Allocate<int>(cpu_range_limit.Length))
                    using (var gpu_slope = accelerator.Allocate<float>(cpu_slope.Length))
                    using (var gpu_rise = accelerator.Allocate<float>(cpu_rise.Length))
                    {
                        gpu_terrain.CopyFrom(ViperEnvironment.Terrain.Data, 0, 0, ViperEnvironment.Terrain.Data.Length);
                        gpu_range.CopyFrom(cpu_range, 0, 0, cpu_range.Length);
                        gpu_slope.CopyFrom(cpu_slope, 0, 0, cpu_slope.Length);

                        var kernel1 = accelerator.LoadAutoGroupedStreamKernel<Index2, int, int, int, int, ArrayView<short>, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<int>>(NearFieldKernel1);

                        while (queue.Count > 0)
                        {
                            var center = queue[0];
                            Console.WriteLine($"Starting patch id=[{center.Id.X},{center.Id.Y}] ...");
                            try
                            {
                                var stopwatch = Stopwatch.StartNew();

                                for (var i = 0; i < cpu_rise.Length; i++)
                                    cpu_rise[i] = float.MinValue;
                                FillRangeLimit(center);
                                gpu_range_limit.CopyFrom(cpu_range_limit, 0, 0, cpu_range_limit.Length);
                                gpu_rise.CopyFrom(cpu_rise, 0, 0, cpu_rise.Length);

                                var launchDimension = new Index2(TerrainPatch.DefaultSize, TerrainPatch.DefaultSize);
                                //var launchDimension = new Index2(1, 1);

                                var batch = 32;
                                for (var slope_idx_start = 0; slope_idx_start < cpu_slope.Length; slope_idx_start += batch)
                                {
                                    var slope_idx_stop = Math.Min(slope_idx_start + batch, cpu_slope.Length);
                                    kernel1(launchDimension, center.Line, center.Sample, slope_idx_start, slope_idx_stop, gpu_terrain, gpu_range, gpu_slope, gpu_rise, gpu_range_limit);
                                    accelerator.Synchronize();  // Needed?
                                }

                                gpu_rise.CopyTo(cpu_rise, 0, 0, cpu_rise.Length);

                                stopwatch.Stop();

                                DumpRise(center, cpu_rise);
                                WriteSlopeArray();

                                // Update queue on disk
                                queue.RemoveAt(0);
                                queue_saver?.Invoke(queue);
                                Console.WriteLine($"  {queue.Count} patches remaining.  Est. completion at {DateTime.Now.AddMinutes(queue.Count * stopwatch.Elapsed.TotalMinutes)}");
                            }
                            catch (Exception e1)
                            {
                                Console.WriteLine(e1);
                                Console.WriteLine(e1.StackTrace);
                            }
                            if (StatusReceiver != null) StatusReceiver.IsRunning = false;
                        }
                    }
                    Console.WriteLine($"Finished queue.  time={stopwatchOuter.Elapsed}.");
                }
            }
        }

        bool GetAccelerator(AcceleratorType t, out AcceleratorId aid)
        {
            aid = Accelerator.Accelerators.Where(id => id.AcceleratorType == t).FirstOrDefault();
            if (aid.AcceleratorType != t)
                Console.WriteLine(@"There is accelerator present of the desired type.  Doing nothing.");
            return aid.AcceleratorType == t;
        }

        void PrepareSlopeArray(StudyInterval i) => PrepareSlopeArray(i.Start, i.Stop, i.Step);
        void PrepareSlopeArray(DateTime start, DateTime stop, TimeSpan step)
        {
            const int bins = 4 * 360;
            var flip = -LonFactor;

            for (var i = 0; i < bins; i++)
                cpu_slope[i] = float.MaxValue;

            var elev = cpu_slope;
            for (var time = start; time <= stop; time += step)
            {
                var vec = CSpice.SunPosition(time);
                var x = vec[0];
                var y = vec[1];
                var z = vec[2];
                var alen = Math.Sqrt(x * x + y * y);
                var azimuth_rad = (float)(Math.Atan2(y, x) + Math.PI);  // [0,2PI]
                var azimuth_index = (int)((azimuth_rad / (2 * Math.PI)) * bins);
                var slope = flip * z / alen;
                if (slope < elev[azimuth_index])
                    elev[azimuth_index] = (float)slope;
            }

            if (false)
            {
                for (var i = 0; i < elev.Length; i++)
                    elev[i] = 0f;
            }
        }

        void FillRangeArray()
        {
            var step = 0.70710710678f;
            var lst = Enumerable.Range(0, small_range_step_count).Select(i => i * step).ToList();

            // Fill this list, then copy to the second range array

            // the largest range in _range1 is d.  Go from there to 2*d, taking a larger step.  Then bump the range and the step again.
            var bump_step = 4f;
            var bump_fence = 2f;
            var max = (float)TerrainPatch.DEM_size;

            step *= bump_step;
            var d = lst[lst.Count - 1] + step;
            var fence = d * bump_fence;
            while (d < max)
            {
                lst.Add(d);
                d += step;
                if (d > fence)
                {
                    step *= bump_step;
                    fence *= bump_fence;
                }
            }
            cpu_range = lst.ToArray();
        }

        void FillRangeLimit(TerrainPatch center)
        {
            var center_line = center.Line;
            var center_sample = center.Sample;
            Debug.Assert(cpu_range_limit.Length == cpu_slope_size);
            var min = (double)(2 * TerrainPatch.DefaultSize);
            var max = (double)(TerrainPatch.DEM_size - 2 * TerrainPatch.DefaultSize);
            for (var i = 0; i < cpu_range_limit.Length; i++)
            {
                var ray_rad = Math.PI * 2f * i / cpu_slope_size;  // 0 deg in ME frame points toward the earth
                var ray_cos = Math.Cos(ray_rad);  // varies the line
                var ray_sin = Math.Sin(ray_rad);  // varies the sample
                var farthest = 0;
                for (var j = 0; j < cpu_range.Length; j++)
                {
                    var range = cpu_range[j];
                    var caster_line = center_line + ray_sin * range;
                    if (caster_line < min || caster_line > max)
                        break;
                    var caster_sample = center_sample + ray_cos * range;
                    if (caster_sample < min || caster_sample > max)
                        break;
                    farthest = j;
                }
                cpu_range_limit[i] = farthest;
            }
        }

        static void NearFieldKernel1(
            Index2 index,      // Contains (sample, line)  First dimension changes fastest
            int target_line,   // row within the 128 x 128 target patch
            int target_sample, // column
            int slope_idx_start,
            int slope_idx_stop,

            ArrayView<short> gpu_terrain,
            ArrayView<float> gpu_range,
            ArrayView<float> gpu_slope,
            ArrayView<float> gpu_rise,
            ArrayView<int> gpu_range_limit)
        {
            var flip = -LonFactor;

            // From ILGPU source code: public int ComputeLinearIndex(Index2 dimension) => Y * dimension.X + X;
            var rise_idx = index.ComputeLinearIndex(Dimension);
            var rise = gpu_rise[rise_idx];

            var center_line = target_line + index.Y;
            var center_sample = target_sample + index.X;
            var center_idx = center_line * TerrainSize + center_sample;
            var center_height = 0.5f * gpu_terrain[center_idx];
            GetVector3d(center_line, center_sample, center_height, out float center_x, out float center_y, out float center_z);
            center_z *= flip;

            for (var slope_idx = slope_idx_start; slope_idx < slope_idx_stop; slope_idx++)
            {
                var slope = gpu_slope[slope_idx];

                // Work out the direction vector
                var ray_rad = GPUMath.PI * 2f * slope_idx / cpu_slope_size;  // 0 deg in ME frame points toward the earth
                var ray_cos = GPUMath.Cos(ray_rad);  // varies the line
                var ray_sin = GPUMath.Sin(ray_rad);  // varies the sample

                // iterate over the ranges
                var range_limit = gpu_range_limit[slope_idx];
                for (var range_idx = 0; range_idx < range_limit; range_idx++)
                {
                    var range = gpu_range[range_idx];
                    var caster_line = center_line + ray_sin * range;
                    var caster_sample = center_sample + ray_cos * range;

                    var x1 = (int)caster_sample;  // Calculate the caster point by interpolating between four points from the points array
                    var y1 = (int)caster_line;
                    int x2 = x1 + 1;
                    int y2 = y1 + 1;

                    var q11_idx = y1 * TerrainSize + x1;
                    var q11 = gpu_terrain[q11_idx];

                    var q12_idx = q11_idx + TerrainSize;
                    var q12 = gpu_terrain[q12_idx];

                    // First interpolation across rows (line)
                    var q1_line = q11 + (caster_line - y1) * (q12 - q11);

                    var q21_idx = q11_idx + 1;
                    var q21 = gpu_terrain[q21_idx];

                    var q22_idx = q11_idx + TerrainSize + 1;
                    var q22 = gpu_terrain[q22_idx];

                    // Second interpolation across rows
                    var q2_line = q21 + (caster_line - y1) * (q22 - q21);

                    // Interpolate across samples and convert to meters
                    var caster_height = q1_line + (caster_sample - x1) * (q2_line - q1_line);
                    caster_height *= 0.5f;

                    GetVector3d(caster_line, caster_sample, caster_height, out float caster_x, out float caster_y, out float caster_z);
                    caster_z *= flip;

                    var dx = caster_x - center_x;
                    var dy = caster_y - center_y;
                    var d = GPUMath.Sqrt(dx * dx + dy * dy); // horizontal distance in moon radius units

                    var light_ray_height = caster_z - slope * d;  // negative slope gets higher as light ray goes toward the center
                    var ray_rise_height = light_ray_height - center_z;  // moon radius units
                    var ray_rise_meters = ray_rise_height * 1000f;

                    // Alternative
                    //var dInMeters = d * 1000f;
                    //var deltaHeightInMeters = (caster_z - center_z) * 1000f;
                    //var rise2 = deltaHeightInMeters - dInMeters * slope;
                    
                    rise = GPUMath.Max(rise, ray_rise_meters);
                }
            }

            gpu_rise[rise_idx] = rise;
        }

        // Returns the position in km
        static void GetVector3d(float line, float sample, float height_meters, out float x, out float y, out float z)
        {
            var radius = MoonRadius + height_meters / 1000f;
            GetLatLon(line, sample, out float lat, out float lon);
            z = radius * GPUMath.Sin(lat);
            var c = radius * GPUMath.Cos(lat);
            x = c * GPUMath.Cos(lon);  // TODO: Not sure about these
            y = c * GPUMath.Sin(lon);
        }

        static void GetLatLon(float line, float sample, out float latitude, out float longitude)
        {
            var x = (sample - S0) * Scale;
            var y = (L0 - line) * Scale;
            var P = GPUMath.Sqrt(x * x + y * y);
            var C = 2f * GPUMath.Atan2(P, 2f * MoonRadius);
            latitude = GPUMath.Asin(GPUMath.Cos(C) * GPUMath.Sin(LatP) + (y == 0 ? 0 : y * GPUMath.Sin(C) * GPUMath.Cos(LatP) / P));
            longitude = LonP + GPUMath.Atan2(x, y * LonFactor);
        }

        #region testing

        void DumpRise(TerrainPatch patch, float[] cpu_rise)
        {
            Debug.Assert(cpu_rise.Length == TerrainPatch.DefaultSize * TerrainPatch.DefaultSize);
            using (var sw = new StreamWriter("rise.csv"))
            {
                for (var row=0;row< TerrainPatch.DefaultSize;row++)
                {
                    for (var col=0;col< TerrainPatch.DefaultSize;col++)
                    {
                        var idx = row * TerrainPatch.DefaultSize + col;
                        var rise = cpu_rise[idx];
                        sw.Write(col == 0 ? "{0}" : ",{0}", rise);
                    }
                    sw.WriteLine();
                }
            }

            using (var sw = new StreamWriter("elevation_plus_rise.csv"))
            {
                var line = patch.Line;
                var sample = patch.Sample;
                for (var row = 0; row < TerrainPatch.DefaultSize; row++)
                {
                    for (var col = 0; col < TerrainPatch.DefaultSize; col++)
                    {
                        var idx = row * TerrainPatch.DefaultSize + col;
                        var rise = cpu_rise[idx];
                        var elevation = 0.5f * Terrain.LineSampleToTerrainOffset(line + row, sample + col);
                        var v = elevation + rise;
                        sw.Write(col == 0 ? "{0}" : ",{0}", v);
                    }
                    sw.WriteLine();
                }
            }

            using (var sw = new StreamWriter("elevation_plus_rise2.csv"))
            {
                var line = patch.Line;
                var sample = patch.Sample;
                for (var row = 0; row < TerrainPatch.DefaultSize; row++)
                {
                    for (var col = 0; col < TerrainPatch.DefaultSize; col++)
                    {
                        var idx = col * TerrainPatch.DefaultSize + row;
                        var rise = cpu_rise[idx];
                        var elevation = 0.5f * Terrain.LineSampleToTerrainOffset(line + row, sample + col);
                        var v = elevation + rise;
                        sw.Write(col == 0 ? "{0}" : ",{0}", v);
                    }
                    sw.WriteLine();
                }
            }

            using (var sw = new StreamWriter("elevation.csv"))
            {
                var line = patch.Line;
                var sample = patch.Sample;
                for (var row = 0; row < TerrainPatch.DefaultSize; row++)
                {
                    for (var col = 0; col < TerrainPatch.DefaultSize; col++)
                    {
                        var idx = row * TerrainPatch.DefaultSize + col;
                        var rise = cpu_rise[idx];
                        var elevation = 0.5f * Terrain.LineSampleToTerrainOffset(line + row, sample + col);
                        var v = elevation;
                        sw.Write(col == 0 ? "{0}" : ",{0}", v);
                    }
                    sw.WriteLine();
                }
            }
        }

        void WriteSlopeArray()
        {
            using (var sw = new StreamWriter("sun_slope.csv"))
                for (var i = 0; i < cpu_slope.Length; i++)
                    sw.Write(i == 0 ? "{0}" : ",{0}", cpu_slope[i]);
        }

        #endregion

        #region Generate Rise Height Dataset

        public List<(TerrainPatch,float[])> GenerateRiseHeightPatches(List<TerrainPatch> queue)
        {
            var result = new List<(TerrainPatch, float[])>();

            Debug.Assert(LonFactor == (float)InMemoryInt16Terrain.LonFactor);
            Debug.Assert(LatP == (float)InMemoryInt16Terrain.LatP);
            Terrain = ViperEnvironment.Terrain;

            var stopwatchOuter = new Stopwatch();
            stopwatchOuter.Start();
            if (StatusReceiver != null) StatusReceiver.IsRunning = true;
            Debug.Assert(ViperEnvironment.Terrain != null);
            if (queue == null)
                queue = ReadShadowCalculationQueue();

            var patchDictionary = new Dictionary<int, TerrainPatch>();
            BoundingBox = null;

            var interval = StudyInterval.MetonicCycle;
            interval.Step = TimeSpan.FromMinutes(20);
            PrepareSlopeArray(interval);
            WriteSlopeArray();
            FillRangeArray();
            cpu_range_limit = new int[cpu_slope.Length];

            using (var context = new Context())
            {
                if (!GetAccelerator(AcceleratorType.Cuda, out AcceleratorId aid))  // DefaultAcceleratorType
                    return result;
                using (var accelerator = Accelerator.Create(context, aid))
                {
                    using (var gpu_terrain = accelerator.Allocate<short>(ViperEnvironment.Terrain.Data.Length))
                    using (var gpu_range = accelerator.Allocate<float>(cpu_range.Length))
                    using (var gpu_range_limit = accelerator.Allocate<int>(cpu_range_limit.Length))
                    using (var gpu_slope = accelerator.Allocate<float>(cpu_slope.Length))
                    using (var gpu_rise = accelerator.Allocate<float>(cpu_rise.Length))
                    {
                        gpu_terrain.CopyFrom(ViperEnvironment.Terrain.Data, 0, 0, ViperEnvironment.Terrain.Data.Length);
                        gpu_range.CopyFrom(cpu_range, 0, 0, cpu_range.Length);
                        gpu_slope.CopyFrom(cpu_slope, 0, 0, cpu_slope.Length);

                        var kernel1 = accelerator.LoadAutoGroupedStreamKernel<Index2, int, int, int, int, ArrayView<short>, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<int>>(NearFieldKernel1);

                        while (queue.Count > 0)
                        {
                            var center = queue[0];
                            Console.WriteLine($"Starting patch id=[{center.Id.X},{center.Id.Y}] ...");
                            try
                            {
                                var stopwatch = Stopwatch.StartNew();

                                for (var i = 0; i < cpu_rise.Length; i++)
                                    cpu_rise[i] = float.MinValue;
                                FillRangeLimit(center);
                                gpu_range_limit.CopyFrom(cpu_range_limit, 0, 0, cpu_range_limit.Length);
                                gpu_rise.CopyFrom(cpu_rise, 0, 0, cpu_rise.Length);

                                var launchDimension = new Index2(TerrainPatch.DefaultSize, TerrainPatch.DefaultSize);
                                //var launchDimension = new Index2(1, 1);

                                var batch = 32;
                                for (var slope_idx_start = 0; slope_idx_start < cpu_slope.Length; slope_idx_start += batch)
                                {
                                    var slope_idx_stop = Math.Min(slope_idx_start + batch, cpu_slope.Length);
                                    kernel1(launchDimension, center.Line, center.Sample, slope_idx_start, slope_idx_stop, gpu_terrain, gpu_range, gpu_slope, gpu_rise, gpu_range_limit);
                                    accelerator.Synchronize();  // Needed?  Note location.  Was right after kernel call
                                }
                                gpu_rise.CopyTo(cpu_rise, 0, 0, cpu_rise.Length);

                                stopwatch.Stop();

                                var copy = new float[cpu_rise.Length];
                                Array.Copy(cpu_rise, copy, cpu_rise.Length);
                                result.Add((center, copy));

                                // Update queue on disk
                                queue.RemoveAt(0);

                                Console.WriteLine($"  {queue.Count} patches remaining.  Est. completion at {DateTime.Now.AddMinutes(queue.Count * stopwatch.Elapsed.TotalMinutes)}");
                            }
                            catch (Exception e1)
                            {
                                Console.WriteLine(e1);
                                Console.WriteLine(e1.StackTrace);
                            }
                            if (StatusReceiver != null) StatusReceiver.IsRunning = false;
                        }
                    }
                    Console.WriteLine($"Finished queue.  time={stopwatchOuter.Elapsed}.");
                }
            }
            return result;
        }

        public void GenerateRiseHeight(string path, Rectangle region)
        {
            var queue = TerrainPatch.EnumerateIds(TerrainPatch.CoveringIdRectangle(region)).Select(id => TerrainPatch.FromId(id)).ToList();
            var pairs = GenerateRiseHeightPatches(queue);
            var ary = new float[region.Height, region.Width];

            var region_width = region.Width;
            var region_height = region.Height;
            foreach (var (patch, data) in pairs)
            {
                var xoffset = patch.Sample - region.Left;
                var yoffset = patch.Line - region.Top;
                for (var row_src = 0; row_src < TerrainPatch.DefaultSize; row_src++)
                {
                    var row_dst = row_src + yoffset;
                    if (row_dst >= region_height || row_dst < 0) continue;
                    for (var col_src = 0; col_src < TerrainPatch.DefaultSize; col_src++)
                    {                        
                        var col_dst = col_src + xoffset;
                        if (col_dst >= region_width || col_dst < 0) continue;  // This can be faster if I don't clip per pixel.  Later.

                        ary[row_dst, col_dst] = data[row_src * TerrainPatch.DefaultSize + col_src];
                    }
                }
            }

            GeotiffHelper.WriteArrayAsGeotiff(ary, region, path);
        }

        #endregion
    }
}