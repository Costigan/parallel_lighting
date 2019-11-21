using ILGPU;
using ILGPU.Runtime;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using viper.corelib.math;
using viper.corelib.patch;
using viper.corelib.spice;
using viper.corelib.utilities;

namespace viper.corelib.lighting
{
    public class RiseToFullSun
    {
        public static double MaximumLocalDistance = 230;  // 1d / Math.Tan(025d * Math.PI/180d)
        //public LunarHorizon MainWindow;

        float[] cpu_slope;  // toward sun
        float[] cpu_vector_x;
        float[] cpu_vector_y;

        const int small_range_step_count = 1024;
        float[] cpu_range;

        public void Generate()
        {
            FillAzimuthElevation();
            FillRangeArrays();
        }

        void FillAzimuthElevation()
        {
            var center_patch = GetCenterPatch();
            var center_pt = center_patch.PointInPatch(new Point(TerrainPatch.DEM_size / 2, TerrainPatch.DEM_size / 2));
            var mat = center_patch.Matrices[center_pt.Y][center_pt.X];
            const int per_degree = 4;
            const int count = 360 * per_degree;
            var slope_array = Enumerable.Range(0, count).Select(i => float.MinValue).ToArray();
            var vec_array = new PointF[count];
            var stop = TileLightingProductManager.MetonicCycleStop;
            var step = TileLightingProductManager.MetonicCycleStep;
            var temp = new Vector3d();
            for (var time = TileLightingProductManager.MetonicCycleStart; time <= stop; time += step)
            {
                var sun = CSpice.SunPosition(time);
                TerrainPatch.Transform(ref sun, ref mat, ref temp);
                var sun_x = temp[0];
                var sun_y = temp[1];
                var sun_z = temp[2];
                var azimuth_rad = Math.Atan2(sun_y, sun_x) + Math.PI;  // [0,2PI]
                var azimuth_deg = azimuth_rad * 180d / Math.PI;
                var index = (int)(azimuth_deg * per_degree);
                var alen = Math.Sqrt(sun_x * sun_x + sun_y * sun_y);
                var slope = (float)(sun_z / alen);
                if (slope > slope_array[index])
                {
                    slope_array[index] = slope;
                    vec_array[index] = new PointF((float)(sun_x / alen), (float)(sun_y / alen));
                }
            }
            var slope_list = new List<float>();
            var vec_list = new List<PointF>();
            for (var i = 0; i < slope_array.Length; i++)
            {
                var s = slope_array[i];
                if (s > float.MinValue)
                {
                    slope_list.Add(s);
                    vec_list.Add(vec_array[i]);
                }
            }
            cpu_slope = slope_array.ToArray();
            cpu_vector_x = vec_list.Select(p => p.X).ToArray();
            cpu_vector_y = vec_list.Select(p => p.Y).ToArray();
        }

        TerrainPatch GetCenterPatch()
        {
            var pt = new Point(TerrainPatch.DEM_size / 2, TerrainPatch.DEM_size / 2);
            var center = TerrainPatch.FromLineSample(pt.Y, pt.X);
            center.FillMatrices(ViperEnvironment.Terrain);
            return center;
        }

        void FillRangeArrays()
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

        public void GenerateHeightField(TerrainPatch patch)
        {
            using (var context = new Context())
            {
                AcceleratorId aid = Accelerator.Accelerators.Where(id => id.AcceleratorType == AcceleratorType.Cuda).FirstOrDefault();
                if (aid.AcceleratorType != AcceleratorType.Cuda)
                {
                    Console.WriteLine(@"There is no CUDA accelerator present.  Doing nothing.");
                    return;
                }
                using (var accelerator = Accelerator.Create(context, aid))
                using (var gpu_dem = accelerator.Allocate<short>(ViperEnvironment.Terrain.Data.Length))
                using (var gpu_range = accelerator.Allocate<float>(cpu_range.Length))
                using (var gpu_slope = accelerator.Allocate<float>(cpu_slope.Length))
                using (var gpu_rise = accelerator.Allocate<short>(TerrainPatch.DefaultSize * TerrainPatch.DefaultSize))
                {

                    gpu_dem.CopyFrom(ViperEnvironment.Terrain.Data, 0, 0, ViperEnvironment.Terrain.Data.Length);
                    gpu_range.CopyFrom(cpu_range, 0, 0, cpu_range.Length);

                    var launchDimension = new Index2(TerrainPatch.DefaultSize, TerrainPatch.DefaultSize);

                    var kernel1 = accelerator.LoadStreamKernel<Index2, ArrayView<short>, ArrayView<float>, ArrayView<short>, int, int>(RiseKernel1);

                }
            }
        }

        static void RiseKernel1(
            Index2 index,
            ArrayView<short> dem,
            ArrayView<float> range,
            ArrayView<short> rise,
            int patch_x,
            int patch_y
            )
        {

        }
    }
}
