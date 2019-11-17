using System;
using System.Linq;
using ILGPU;
using ILGPU.Runtime;
using viper.corelib.patch;

namespace viper.corelib.horizon
{
    public class GPUTest
    {
        public const double MoonRadius = 1737.4d;
        public const double S0 = 15199.5d;             // PDS SAMPLE_PROJECTION_OFFSET
        public const double L0 = 15199.5d;
        public const double Scale = 20d / 1000d;
        public const double LonP = 0d;

        // These are only correct for the north DEM
        public const double LonFactor = -1d;
        public const double LatP = Math.PI / 2;

        public void Test1()
        {
            using (var context = new Context())
            {
                var cudaid = Accelerator.Accelerators.Where(id => id.AcceleratorType == AcceleratorType.Cuda).FirstOrDefault();
                if (cudaid.AcceleratorType != AcceleratorType.Cuda)
                {
                    Console.WriteLine(@"There is no CUDA accelerator present.  Doing nothing.");
                    return;
                }
                using (var cuda = Accelerator.Create(context, cudaid))
                {
                    Console.WriteLine($"Performing operations on {cuda}");
                    PrintAcceleratorInfo(cuda);
                }
            }
        }

        static void PrintAcceleratorInfo(Accelerator accelerator)
        {
            Console.WriteLine($"Name: {accelerator.Name}");
            Console.WriteLine($"MemorySize: {accelerator.MemorySize}");
            Console.WriteLine($"MaxThreadsPerGroup: {accelerator.MaxNumThreadsPerGroup}");
            Console.WriteLine($"MaxSharedMemoryPerGroup: {accelerator.MaxSharedMemoryPerGroup}");
            Console.WriteLine($"MaxGridSize: {accelerator.MaxGridSize}");
            Console.WriteLine($"MaxConstantMemory: {accelerator.MaxConstantMemory}");
            Console.WriteLine($"WarpSize: {accelerator.WarpSize}");
            Console.WriteLine($"NumMultiprocessors: {accelerator.NumMultiprocessors}");
        }

        internal float[] Test2(TerrainPatch a, TerrainPatch o, short[] relativeHeights)
        {
            using (var context = new Context())
            {
                AcceleratorId aid;
                if (true)
                {
                    aid = Accelerator.Accelerators.Where(id => id.AcceleratorType == AcceleratorType.Cuda).FirstOrDefault();
                    if (aid.AcceleratorType != AcceleratorType.Cuda)
                    {
                        Console.WriteLine(@"There is no CUDA accelerator present.  Doing nothing.");
                        return null;
                    }
                }
                else
                {
                    aid = Accelerator.Accelerators.Where(id => id.AcceleratorType != AcceleratorType.Cuda).FirstOrDefault();
                }
                using (var accelerator = Accelerator.Create(context, aid))
                {
                    Console.WriteLine($"Performing operations on {accelerator}");
                    PrintAcceleratorInfo(accelerator);

                    System.Diagnostics.Debug.Assert(a.Height == 128 && a.Width == 128);

                    // Matrices
                    var cpu_matrixes_size = a.Height * a.Width * 12;
                    var cpu_matrixes = new double[cpu_matrixes_size];
                    {
                        var r = cpu_matrixes;
                        var height = a.Height;
                        var width = a.Width;
                        for (var l = 0; l < height; l++)
                            for (var s = 0; s < width; s++)
                            {
                                var mat = a.Matrices[l][s];
                                var pos = (l * width + s) * 12;

                                r[pos++] = mat.Row0.X;
                                r[pos++] = mat.Row1.X;
                                r[pos++] = mat.Row2.X;
                                r[pos++] = mat.Row3.X;

                                r[pos++] = mat.Row0.Y;
                                r[pos++] = mat.Row1.Y;
                                r[pos++] = mat.Row2.Y;
                                r[pos++] = mat.Row3.Y;

                                r[pos++] = mat.Row0.Z;
                                r[pos++] = mat.Row1.Z;
                                r[pos++] = mat.Row2.Z;
                                r[pos++] = mat.Row3.Z;
                            }
                    }

                    // Heights
                    var gpu_heights_size = o.Height * o.Width;

                    // Horizon
                    var cpu_horizon_size = a.Height * a.Width * 360 * 4;
                    var cpu_horizon = new float[cpu_horizon_size];
                    for (var i = 0; i < cpu_horizon_size; i++)
                        cpu_horizon[i] = -10f;

                    using (var gpu_matrices = accelerator.Allocate<double>(cpu_matrixes_size))
                    using (var gpu_heights = accelerator.Allocate<short>(gpu_heights_size))
                    using (var gpu_horizon = accelerator.Allocate<float>(cpu_horizon_size))
                    {
                        gpu_matrices.CopyFrom(cpu_matrixes, 0, 0, cpu_matrixes_size);
                        gpu_heights.CopyFrom(relativeHeights, 0, 0, gpu_heights_size);
                        gpu_horizon.CopyFrom(cpu_horizon, 0, 0, cpu_horizon_size);

                        var kernel = accelerator.LoadSharedMemoryStreamKernel1<
                            GroupedIndex, ArrayView<short>, ArrayView<double>, ArrayView<float>, int, int, int, ArrayView<float>>(Test2Kernel);

                        var groupSize = Math.Min(accelerator.MaxNumThreadsPerGroup, 128);
                        var dimension = new GroupedIndex(1, groupSize);

                        var stopwatch = new System.Diagnostics.Stopwatch();
                        stopwatch.Start();

                        var defaultStream = accelerator.DefaultStream;


                        for (var aline_init = 0; aline_init < 128; aline_init++)
                        {
                            kernel(dimension, gpu_heights.View, gpu_matrices.View, gpu_horizon.View, o.Line, o.Sample, aline_init);

                            accelerator.Synchronize();
                        }

                        stopwatch.Stop();
                        Console.WriteLine($"kernel time={stopwatch.Elapsed}");

                        // Copy out data
                        gpu_horizon.CopyTo(cpu_horizon, 0, 0, cpu_horizon_size);
                    }

                    return cpu_horizon;
                }
            }
        }

        static void Test2Kernel(
            GroupedIndex index,
            ArrayView<short> heights,
            ArrayView<double> matrices,
            ArrayView<float> horizon,
            int a_line,
            int a_sample,
            int aline_init,

            [SharedMemory(1440)]
            ArrayView<float> horizon_shared)
        {
            var idx = index.GroupIdx;
            const int patch_size = 128;

            // Do the calculation
            for (var aline = aline_init; aline == aline_init; aline++)  // NOTE: constant here should be passed in
                for (var asample = 0; asample < patch_size; asample++)
                {
                    // Copy horizon for a[line,sample] into shared memory
                    {
                        var dim = Group.Dimension.X;
                        var len = horizon_shared.Length;
                        var passes = (len + (dim - 1)) / dim;
                        var offset = (aline * patch_size + asample) * len;
                        for (var pass = 0; pass < passes; pass++)
                        {
                            var ptr = pass * dim + idx;
                            if (ptr < len)
                                horizon_shared[ptr] = horizon[ptr + offset];
                            else
                                idx = idx;
                        }
                    }

                    //Group.Barrier();

                    // Copy the matrix into registers
                    var pos = (aline * patch_size + asample) * 12;

                    var row0x = matrices[pos++];
                    var row1x = matrices[pos++];
                    var row2x = matrices[pos++];
                    var row3x = matrices[pos++];

                    var row0y = matrices[pos++];
                    var row1y = matrices[pos++];
                    var row2y = matrices[pos++];
                    var row3y = matrices[pos++];

                    var row0z = matrices[pos++];
                    var row1z = matrices[pos++];
                    var row2z = matrices[pos++];
                    var row3z = matrices[pos];

                    for (var oline = 0; oline < patch_size; oline++)
                    {
                        // osample = idx
                        var relz = 0.5d * heights[aline * patch_size + idx];
                        var radius = MoonRadius + relz / 1000d;
                        var line = a_line + aline;
                        var sample = a_sample + idx;

                        var map_x = (sample - S0) * Scale;
                        var map_y = (L0 - line) * Scale;
                        var P = Math.Sqrt(map_x * map_x + map_y * map_y);
                        var C = 2d * Math.Atan2(P, 2 * MoonRadius);
                        var latitude = Math.Asin(Math.Cos(C) * Math.Sin(LatP) + map_y * Math.Sin(C) * Math.Cos(LatP) / P);
                        var longitude = LonP + Math.Atan2(map_x, map_y * LonFactor);

                        var latdeg = latitude * 180d / Math.PI;
                        var londeg = longitude * 180d / Math.PI;

                        // Calculate the other point in ME frame
                        var z_me = radius * Math.Sin(latitude);
                        var c = radius * Math.Cos(latitude);
                        var x_me = c * Math.Cos(longitude);
                        var y_me = c * Math.Sin(longitude);

                        // Transform the point to the local frame
                        var x = x_me * row0x + y_me * row1x + z_me * row2x + row3x;
                        var y = x_me * row0y + y_me * row1y + z_me * row2y + row3y;
                        var z = x_me * row0z + y_me * row1z + z_me * row2z + row3z;

//                        if (idx == 0)
//                            relz = relz;

                        var azimuth = Math.Atan2(y, x) + Math.PI;  // [0,2 PI]
                        var alen = Math.Sqrt(x * x + y * y);
                        var slope = z / alen;
                        var slopef = (float)slope;

                        var horizon_index = (int)(0.5d + 1439 * (azimuth / (2d * Math.PI)));
                        Atomic.Max(horizon_shared.GetVariableView(horizon_index), slopef);
                        //horizon_shared[horizon_index] = 1f;
                    }

                    //Group.Barrier();

                    {
                        var dim = Group.Dimension.X;
                        var len = horizon_shared.Length;
                        var passes = (len + (dim - 1)) / dim;
                        var offset = (aline * patch_size + asample) * len;
                        for (var pass = 0; pass < passes; pass++)
                        {
                            var ptr = pass * dim + idx;
                            if (ptr < len)
                                horizon[ptr + offset] = horizon_shared[ptr];
                            else
                                idx = idx;
                        }
                    }
                }
        }

        internal float[] Test3(TerrainPatch a, TerrainPatch o, short[] relativeHeights)
        {
            using (var context = new Context())
            {
                AcceleratorId aid;
                if (true)
                {
                    aid = Accelerator.Accelerators.Where(id => id.AcceleratorType == AcceleratorType.Cuda).FirstOrDefault();
                    if (aid.AcceleratorType != AcceleratorType.Cuda)
                    {
                        Console.WriteLine(@"There is no CUDA accelerator present.  Doing nothing.");
                        return null;
                    }
                }
                else
                {
                    aid = Accelerator.Accelerators.Where(id => id.AcceleratorType != AcceleratorType.Cuda).FirstOrDefault();
                }
                using (var accelerator = Accelerator.Create(context, aid))
                {
                    Console.WriteLine($"Performing operations on {accelerator}");
                    PrintAcceleratorInfo(accelerator);

                    System.Diagnostics.Debug.Assert(a.Height == 128 && a.Width == 128);

                    // Matrices
                    var cpu_matrixes_size = a.Height * a.Width * 12;
                    var cpu_matrixes = new float[cpu_matrixes_size];
                    {
                        var r = cpu_matrixes;
                        var height = a.Height;
                        var width = a.Width;
                        for (var l = 0; l < height; l++)
                            for (var s = 0; s < width; s++)
                            {
                                var mat = a.Matrices[l][s];
                                var pos = (l * width + s) * 12;

                                r[pos++] = (float)mat.Row0.X;
                                r[pos++] = (float)mat.Row1.X;
                                r[pos++] = (float)mat.Row2.X;
                                r[pos++] = (float)mat.Row3.X;

                                r[pos++] = (float)mat.Row0.Y;
                                r[pos++] = (float)mat.Row1.Y;
                                r[pos++] = (float)mat.Row2.Y;
                                r[pos++] = (float)mat.Row3.Y;

                                r[pos++] = (float)mat.Row0.Z;
                                r[pos++] = (float)mat.Row1.Z;
                                r[pos++] = (float)mat.Row2.Z;
                                r[pos++] = (float)mat.Row3.Z;
                            }
                    }

                    // Heights
                    var gpu_heights_size = o.Height * o.Width;

                    // Horizon
                    var cpu_horizon_size = a.Height * a.Width * 360 * 4;
                    var cpu_horizon = new float[cpu_horizon_size];
                    for (var i = 0; i < cpu_horizon_size; i++)
                        cpu_horizon[i] = float.MinValue;

                    using (var gpu_matrices = accelerator.Allocate<float>(cpu_matrixes_size))
                    using (var gpu_heights = accelerator.Allocate<short>(gpu_heights_size))
                    using (var gpu_horizon = accelerator.Allocate<float>(cpu_horizon_size))
                    {
                        gpu_matrices.CopyFrom(cpu_matrixes, 0, 0, cpu_matrixes_size);
                        gpu_heights.CopyFrom(relativeHeights, 0, 0, gpu_heights_size);
                        gpu_horizon.CopyFrom(cpu_horizon, 0, 0, cpu_horizon_size);

                        var kernel = accelerator.LoadSharedMemoryStreamKernel1<GroupedIndex, ArrayView<short>, ArrayView<float>, ArrayView<float>, int, int, ArrayView<float>>(Test3Kernel);

                        var groupSize = Math.Min(accelerator.MaxNumThreadsPerGroup, 128);
                        var dimension = new GroupedIndex(128, groupSize);

                        var stopwatch = new System.Diagnostics.Stopwatch();
                        stopwatch.Start();

                        // Note that *no* value is passed for the shared-memory variable
                        // since shared memory is handled automatically inside the runtime
                        // and shared memory has to be initialized inside a kernel.
                        // The delegate type for this kernel would be:
                        // Action<GroupedIndex, ArrayView<int>, ArrayView<int>>.
                        kernel(dimension, gpu_heights.View, gpu_matrices.View, gpu_horizon.View, o.Line, o.Sample);

                        accelerator.Synchronize();

                        stopwatch.Stop();
                        Console.WriteLine($"kernel time={stopwatch.Elapsed}");

                        // Copy out data
                        gpu_horizon.CopyTo(cpu_horizon, 0, 0, cpu_horizon_size);
                    }

                    return cpu_horizon;
                }
            }
        }

        static void Test3Kernel(
            GroupedIndex index,
            ArrayView<short> heights,
            ArrayView<float> matrices,
            ArrayView<float> horizon,
            int a_line,
            int a_sample,

            [SharedMemory(1440)]
            ArrayView<float> horizon_shared)
        {
            var idx = index.GroupIdx;
            const int patch_size = 128;

            // Do the calculation
            var aline = index.GridIdx;
            for (var asample = 0; asample < patch_size; asample++)
            {
                // Copy horizon for a[line,sample] into shared memory
                {
                    var dim = Group.Dimension.X;
                    var len = horizon_shared.Length;
                    var passes = (len + (dim - 1)) / dim;
                    var offset = (aline * patch_size + asample) * len;
                    for (var pass = 0; pass < passes; pass++)
                    {
                        var ptr = pass * dim + idx;
                        if (ptr < len)
                            horizon_shared[ptr] = horizon[ptr + offset];
                        // Note warp divergence
                    }
                }

                Group.Barrier();

                // Copy the matrix into registers
                var pos = (aline * patch_size + asample) * 12;

                var row0x = matrices[pos++];
                var row1x = matrices[pos++];
                var row2x = matrices[pos++];
                var row3x = matrices[pos++];

                var row0y = matrices[pos++];
                var row1y = matrices[pos++];
                var row2y = matrices[pos++];
                var row3y = matrices[pos++];

                var row0z = matrices[pos++];
                var row1z = matrices[pos++];
                var row2z = matrices[pos++];
                var row3z = matrices[pos];

                for (var oline = 0; oline < patch_size; oline++)
                {
                    // osample = idx
                    var relz = 0.5d * heights[aline * patch_size + idx];
                    var radius = MoonRadius + relz / 1000d;
                    var line = a_line + aline;
                    var sample = a_sample + idx;

                    var map_x = (sample - S0) * Scale;
                    var map_y = (L0 - line) * Scale;
                    var P = Math.Sqrt(map_x * map_x + map_y * map_y);
                    var C = 2d * Math.Atan2(P, 2 * MoonRadius);
                    var latitude = Math.Asin(Math.Cos(C) * Math.Sin(LatP) + map_y * Math.Sin(C) * Math.Cos(LatP) / P);
                    var longitude = LonP + Math.Atan2(map_x, map_y * LonFactor);

                    var latdeg = latitude * 180d / Math.PI;
                    var londeg = longitude * 180d / Math.PI;

                    // Calculate the other point in ME frame
                    var z_me = radius * Math.Sin(latitude);
                    var c = radius * Math.Cos(latitude);
                    var x_me = c * Math.Cos(longitude);
                    var y_me = c * Math.Sin(longitude);

                    // Transform the point to the local frame
                    var x = x_me * row0x + y_me * row1x + z_me * row2x + row3x;
                    var y = x_me * row0y + y_me * row1y + z_me * row2y + row3y;
                    var z = x_me * row0z + y_me * row1z + z_me * row2z + row3z;

                    //                        if (idx == 0)
                    //                            relz = relz;

                    var azimuth = Math.Atan2(y, x) + Math.PI;  // [0,2 PI]
                    var alen = Math.Sqrt(x * x + y * y);
                    var slope = z / alen;
                    var slopef = (float)slope;

                    var horizon_index = (int)(0.5d + 1439 * (azimuth / (2d * Math.PI)));
                    Atomic.Max(horizon_shared.GetVariableView(horizon_index), slopef);
                    //horizon_shared[horizon_index] = 1f;
                }

                Group.Barrier();

                {
                    var dim = Group.Dimension.X;
                    var len = horizon_shared.Length;
                    var passes = (len + (dim - 1)) / dim;
                    var offset = (aline * patch_size + asample) * len;
                    for (var pass = 0; pass < passes; pass++)
                    {
                        var ptr = pass * dim + idx;
                        if (ptr < len)
                            horizon[ptr + offset] = horizon_shared[ptr];
                        // Note warp divergence
                    }
                }
            }
        }

        /*
            output.X = vec.X * mat.Row0.X + vec.Y * mat.Row1.X + vec.Z * mat.Row2.X + mat.Row3.X;
            output.Y = vec.X * mat.Row0.Y + vec.Y * mat.Row1.Y + vec.Z * mat.Row2.Y + mat.Row3.Y;
            output.Z = vec.X * mat.Row0.Z + vec.Y * mat.Row1.Z + vec.Z * mat.Row2.Z + mat.Row3.Z;
            */
    }
}
