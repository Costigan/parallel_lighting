using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using ILGPU;
using ILGPU.Runtime;
using System.Drawing;
using viper.corelib.patch;
using viper.corelib.terrain;
using viper.corelib.math;

namespace viper.corelib.horizon
{
    public class GPUHorizons : CPUHorizons
    {
        public AcceleratorType DefaultAcceleratorType = AcceleratorType.Cuda;
        public Rectangle? BoundingBox = null;

        // TODO: remove overHorizonCheck

        public override void RunQueue(List<TerrainPatch> queue, Action<List<TerrainPatch>> queue_saver = null, int spread = -1, bool overHorizonCheck = true, bool centerOnly = false, bool center = true, bool unloadHorizons = true, bool writeHorizons = true)
        {
            var stopwatchOuter = new Stopwatch();
            stopwatchOuter.Start();
            if (StatusReceiver != null) StatusReceiver.IsRunning = true;
            Debug.Assert(Terrain != null);
            centerOnly |= NearHorizonsOnly;
            if (queue == null)
                queue = ReadShadowCalculationQueue();
            if (spread < 0) spread = MaxSpread;

            var run_time = PatchLight.EstimateCalculationTime(queue);
            Console.WriteLine($"Queue length          = {queue.Count}");
            Console.WriteLine($"Estimated run time    = {run_time}");
            Console.WriteLine($"Estimated finish time = {DateTime.Now + run_time}");

            var patchDictionary = new Dictionary<int, TerrainPatch>();
            BoundingBox = null;

            while (queue.Count > 0)
            {
                var patch = queue[0];
                Console.WriteLine($"Starting patch id=[{patch.Id.X},{patch.Id.Y}] observer={patch.ObserverHeightInMeters} m ...");
                try
                {
                    var stopwatch = Stopwatch.StartNew();

                    // Override observer height if specified in generator
                    if (ObserverHeightInMeters.HasValue)
                        patch.ObserverHeightInMeters = ObserverHeightInMeters.Value;

                    var far_field = GenerateFarField(patch.Id, patchDictionary);
                    if (StatusReceiver != null)
                    {
                        if (StatusReceiver.ProcessingPatches == null)
                            StatusReceiver.ProcessingPatches = new List<TerrainPatch>();
                        StatusReceiver.ProcessingPatches.Add(patch);
                        //MapView.ProcessingPatches = new List<TerrainPatch> { patch };
                        //MapView.FarPatches = far_field;
                    }
                    Console.WriteLine($"  The far field contains {far_field.Count} patches");
                    //Console.WriteLine($"The patch cache contains {patchDictionary.Count} patches");

                    Parallel.Invoke(
                        () =>
                        {
                            var stopwatch2 = new Stopwatch();
                            stopwatch2.Start();
                            // Add Nearby horizon.  Horizons will be slope, not angles
                            if (center)
                                ProcessNearField(patch);
                            //Console.WriteLine($"  near horizon calculated or loaded in {stopwatch2.Elapsed}");
                        },
                        () =>
                        {
                            var stopwatch2 = new Stopwatch();
                            stopwatch2.Start();
                            foreach (var p in far_field)
                                p.FillPoints(Terrain);
                            //Console.WriteLine($"Finished loading far field patches in {stopwatch2.Elapsed}");
                        }
                        );

                    patch.Matrices = null;      // Force the matrices to be regenerated so they're not left over from the near horizon
                    patch.FillPoints(Terrain);  // Should be a no op
                    patch.FillMatricesRelativeToPoint(Terrain, patch.Points[0][0]);

                    if (patch.Horizons == null)
                        patch.InitializeHorizons();

                    // Here is the GPU Call
                    ProcessFarField(patch, far_field);

                    foreach (var p in far_field)
                    {
                        var r = new Rectangle(p.Id, new Size(1, 1));
                        if (BoundingBox.HasValue)
                            BoundingBox = Rectangle.Union(BoundingBox.Value, r);
                        else
                            BoundingBox = r;
                    }

                    /*
                    var towrite = MapView.FarPatches.Select(p => p.Id.Y * 512 + p.Id.X).OrderBy(p => p).ToList();
                    Console.WriteLine($"Far patches for {patch.Id}:");
                    {
                        var count = 0;
                        foreach (var far in towrite)
                        {
                            Console.Write($"{far} ");
                            if (count++ > 20)
                            {
                                Console.WriteLine();
                                count = 0;
                            }
                        }
                        Console.WriteLine();
                    }
                    */
                    patch.Path = patch.DefaultPath;
                    if (writeHorizons)
                        patch.Write();               // This converts the horizon format to angles
                    if (unloadHorizons)
                    {
                        patch.Horizons = null;
                        patch.Points = null;
                        patch.Matrices = null;
                        //patch.InitializeHorizons();  // Unload the horizon data (100MB)
                    }
                    stopwatch.Stop();
                    var seconds_per_patch = far_field.Count == 0 ? 0f : (stopwatch.ElapsedMilliseconds / 1000f) / far_field.Count;
                    Console.WriteLine($"  Finished [{patch.Line},{patch.Sample}] time={stopwatch.Elapsed}.  sec/patch={seconds_per_patch}");

                    // Update queue on disk
                    queue.RemoveAt(0);
                    queue_saver?.Invoke(queue);
                    Console.WriteLine($"  {queue.Count} patches remaining.  Est. completion at {DateTime.Now.AddMinutes(queue.Count)}");
                }
                catch (Exception e1)
                {
                    Console.WriteLine(e1);
                    Console.WriteLine(e1.StackTrace);
                }
                if (StatusReceiver != null) StatusReceiver.IsRunning = false;
            }
            Console.WriteLine($"Finished queue.  time={stopwatchOuter.Elapsed}.");
        }

        List<TerrainPatch> GenerateFarField(Point id, Dictionary<int, TerrainPatch> dict)
        {
            var center = new PatchLight { Id = id, Step = 1 };
            var other = new PatchLight { Id = new Point(0, 0), Step = 256 };
            var patches = new List<PatchLight>();
            center.GenerateFarField(other, patches);
            return patches.Select(p =>
            {
                var hash = p.IdHash;
                if (dict.TryGetValue(hash, out TerrainPatch tp))
                    return tp;
                else
                {
                    tp = p.ToTerrainPatch();
                    dict.Add(hash, tp);
                    return tp;
                }
            }).ToList();
        }

        /// <summary>
        /// Update the horizons of a patch based on a list of shadow casters.
        /// The horizons will be in slope, not angle, format
        /// </summary>
        /// <param name="target"></param>
        /// <param name="casters"></param>
        public void ProcessFarField(TerrainPatch target, List<TerrainPatch> casters)
        {
            Debug.Assert(Terrain != null);
            if (casters.Count < 1) return;
            using (var context = new Context())
            {
                if (!GetAccelerator(DefaultAcceleratorType, out AcceleratorId aid))
                    return;
                using (var accelerator = Accelerator.Create(context, aid))
                {
                    target.FillPoints(Terrain);
                    target.FillMatricesRelativeToPoint(Terrain, target.Points[0][0]);

                    // Matrices
                    var cpu_matrices_size = target.Height * target.Width * 12;
                    var basePoint = target.Points[0][0];
                    var cpu_matrices = MakeCPUMatrices(target);

                    // Horizon (load from target)
                    var cpu_horizon_size = target.Height * target.Width * Horizon.HorizonSamples;
                    var cpu_horizon = new int[cpu_horizon_size];
                    for (var line = 0; line < TerrainPatch.DefaultSize; line++)
                        for (var sample = 0; sample < TerrainPatch.DefaultSize; sample++)
                        {
                            var offset = (line * TerrainPatch.DefaultSize + sample) * Horizon.HorizonSamples;
                            var buffer = target.Horizons[line][sample].Buffer;
                            for (var i = 0; i < Horizon.HorizonSamples; i++)
                                cpu_horizon[i + offset] = SlopeToEncoding(buffer[i]);
                        }

                    // Caster points
                    var cpu_caster_points_size = casters[0].Width * casters[0].Height * 3;
                    var cpu_caster_points = new float[cpu_caster_points_size];

                    // test array
                    var cpu_test_array = new float[20];

                    using (var gpu_matrices = accelerator.Allocate<float>(cpu_matrices_size))
                    using (var gpu_horizon = accelerator.Allocate<int>(cpu_horizon_size))
                    using (var gpu_caster_points = accelerator.Allocate<float>(cpu_caster_points_size))
                    using (var gpu_test_array = accelerator.Allocate<float>(cpu_test_array.Length))
                    {
                        gpu_matrices.CopyFrom(cpu_matrices, 0, 0, cpu_matrices_size);
                        gpu_horizon.CopyFrom(cpu_horizon, 0, 0, cpu_horizon_size);

                        //var groupSize = accelerator.MaxNumThreadsPerGroup;
                        var launchDimension = new GroupedIndex2(
                            new Index2(128, 128), // (data.Length + groupSize - 1) / groupSize,  // Compute the number of groups (round up)
                            new Index2(1, 128));

                        var kernel1 = accelerator.LoadSharedMemoryStreamKernel1<GroupedIndex2, ArrayView<float>, ArrayView<float>, ArrayView<int>, ArrayView<float>, float, ArrayView<int>>(FarFieldKernel1);

                        //var stopwatch = new Stopwatch();
                        //stopwatch.Start();

                        var observer_height_in_km = target.ObserverHeightInMeters / 1000f;  // convert to km
                        foreach (var caster in casters)
                        {
                            caster.FillPoints(Terrain);
                            CopyPointsToCpuArray(caster, basePoint, cpu_caster_points);
                            gpu_caster_points.CopyFrom(cpu_caster_points, 0, 0, cpu_caster_points_size);

                            kernel1(launchDimension, gpu_caster_points, gpu_matrices, gpu_horizon, gpu_test_array, observer_height_in_km);  

                            accelerator.Synchronize();
                        }

                        // Copy out data
                        gpu_horizon.CopyTo(cpu_horizon, 0, 0, cpu_horizon_size);
                        gpu_test_array.CopyTo(cpu_test_array, 0, 0, cpu_test_array.Length);

                        //stopwatch.Stop();
                        //Console.WriteLine($"kernel time={stopwatch.Elapsed} cpu_horizon.Max()={cpu_horizon.Max()} cpu_horizon[0]={cpu_horizon[0]}");

                        // Update the horizons
                        for (var line = 0; line < TerrainPatch.DefaultSize; line++)
                            for (var sample = 0; sample < TerrainPatch.DefaultSize; sample++)
                            {
                                var offset = (line * TerrainPatch.DefaultSize + sample) * Horizon.HorizonSamples;
                                var buffer = target.Horizons[line][sample].Buffer;
                                for (var i = 0; i < Horizon.HorizonSamples; i++)
                                    buffer[i] = EncodingToSlope(cpu_horizon[i + offset]);
                            }

                        //Console.WriteLine($"  max slope={cpu_horizon.Select(EncodingToSlope).Max()}");
                    }
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

        internal static float[] MakeCPUMatrices(TerrainPatch target)
        {
            var cpu_matrixes_size = target.Height * target.Width * 12;
            var r = new float[cpu_matrixes_size];
            var height = target.Height;
            var width = target.Width;
            for (var l = 0; l < height; l++)
                for (var s = 0; s < width; s++)
                {
                    var mat = target.Matrices[l][s];
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
            return r;
        }

        private void CopyPointsToCpuArray(TerrainPatch caster, Vector3d basePoint, float[] cpu_caster_points)
        {
            var ptr = 0;
            for (var line = 0; line < TerrainPatch.DefaultSize; line++)
            {
                var row = caster.Points[line];
                for (var sample = 0; sample < TerrainPatch.DefaultSize; sample++)
                {
                    var vec = row[sample] - basePoint;
                    cpu_caster_points[ptr++] = (float)vec.X;
                    cpu_caster_points[ptr++] = (float)vec.Y;
                    cpu_caster_points[ptr++] = (float)vec.Z;
                }
            }
        }

        static void FarFieldKernel1(
            GroupedIndex2 index,
            ArrayView<float> points,
            ArrayView<float> matrices,
            ArrayView<int> horizon,
            ArrayView<float> test_array,
            float observer_height_in_km,

            [SharedMemory(1440)]
            ArrayView<int> horizon_shared)
        {
            var target_line = index.GridIdx.Y;
            var target_sample = index.GridIdx.X;
            var caster_line = index.GroupIdx.Y;
            Debug.Assert(index.GroupIdx.X == 1);

            // Copy horizon for a[target_line,target_sample] into shared memory
            {
                var dim = Group.Dimension.Y;
                var len = horizon_shared.Length;
                var passes = (len + (dim - 1)) / dim;
                var offset = (target_line * TerrainPatch.DefaultSize + target_sample) * len;
                for (var pass = 0; pass < passes; pass++)
                {
                    var ptr = pass * dim + caster_line;
                    if (ptr < len)  // divergence
                        horizon_shared[ptr] = horizon[ptr + offset];
                }
            }

            Group.Barrier();

            // Copy the matrix into registers
            var pos = (target_line * TerrainPatch.DefaultSize + target_sample) * 12;

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

            for (var caster_sample = 0; caster_sample < TerrainPatch.DefaultSize; caster_sample++)
            {
                // Fetch the other point in local frame
                var points_offset = (caster_line * TerrainPatch.DefaultSize + caster_sample) * 3;
                var x_patch = points[points_offset];
                var y_patch = points[points_offset + 1];
                var z_patch = points[points_offset + 2];

                // Transform the point to the local frame
                var x = x_patch * row0x + y_patch * row1x + z_patch * row2x + row3x;
                var y = x_patch * row0y + y_patch * row1y + z_patch * row2y + row3y;
                var z = x_patch * row0z + y_patch * row1z + z_patch * row2z + row3z;

                // Adjust for solar array height
                z -= observer_height_in_km;

                var azimuth = GPUMath.Atan2(y, x) + GPUMath.PI;  // [0,2 PI]
                var alen = GPUMath.Sqrt(x * x + y * y);
                var slope = z / alen;

                var slopem = slope > 2f ? 2f : slope;
                slopem = slopem < -2f ? -2f : slopem;
                slopem = slopem / 4f;

                var slopei = (int)(slopem * 1000000);

                var horizon_index = (int)(0.5f + 1439 * (azimuth / (2f * GPUMath.PI)));
                Atomic.Max(horizon_shared.GetVariableView(horizon_index), slopei);

                if (caster_sample == 0 && caster_line == 0 && target_line == 0 && target_sample == 0)
                {
                    test_array[0] = x_patch;
                    test_array[1] = y_patch;
                    test_array[2] = z_patch;
                    test_array[3] = x;
                    test_array[4] = y;
                    test_array[5] = z;
                    test_array[6] = slope;
                    test_array[7] = slopem;
                    test_array[8] = slopei;
                    test_array[9] = row3x;
                    test_array[10] = row3y;
                    test_array[11] = row3z;
                }
            }

            Group.Barrier();

            {
                var dim = Group.Dimension.Y;
                var len = horizon_shared.Length;
                var passes = (len + (dim - 1)) / dim;
                var offset = (target_line * TerrainPatch.DefaultSize + target_sample) * len;
                for (var pass = 0; pass < passes; pass++)
                {
                    var ptr = pass * dim + caster_line;
                    if (ptr < len)  // divergence
                        horizon[ptr + offset] = horizon_shared[ptr];
                }
            }
        }

        static int SlopeToEncoding(float slope) => (int)((slope / 4f) * 1000000);
        static float EncodingToSlope(int encoding) => 4f * encoding / 1000000f;  // Haven't looked at the numerics closely to see whether these are exact inverses.  Don't care enough yet.

        /// <summary>
        /// Generate the near horizon for the patch, overwriting what was there.
        /// The horizons will be in slope format, not angles
        /// </summary>
        /// <param name="target"></param>
        public void ProcessNearField(TerrainPatch target)
        {
            Debug.Assert(Terrain != null);
            Debug.Assert(target != null);
            if (target.Horizons == null)
                target.InitializeHorizons();
            using (var context = new Context())
            {
                if (!GetAccelerator(DefaultAcceleratorType, out AcceleratorId aid))
                    return;
                using (var accelerator = Accelerator.Create(context, aid))
                {
                    target.Matrices = null;  // Be sure!!
                    target.FillPoints(Terrain);
                    target.FillMatricesRelativeToPoint(Terrain, target.Points[0][0]);

                    // Matrices
                    var cpu_matrices_size = target.Height * target.Width * 12;
                    var basePoint = target.Points[0][0];
                    var cpu_matrices = MakeCPUMatrices(target);

                    // Horizon (load from target)
                    var cpu_horizon_size = target.Height * target.Width * Horizon.HorizonSamples;
                    var cpu_horizon = new float[cpu_horizon_size];

                    // Initialize to float.MinValue.  Maintain horizon as slope
                    for (var i = 0; i < cpu_horizon_size; i++)
                        cpu_horizon[i] = float.MinValue;

                    // Caster array
                    const int dem_size = TerrainPatch.DEM_size;
                    const int patch_size = TerrainPatch.DefaultSize;
                    var maxDistance = (float)TerrainPatch.MaximumLocalDistance;
                    var border = 2 + (int)Math.Ceiling(maxDistance); // the 2 is margin for the bilinear interpolation

                    var line_min = Math.Max(0, target.Line - border);
                    var line_max = Math.Min(dem_size - 1, target.Line + patch_size + border);  // 1 more than the highest index
                    var line_size = line_max - line_min;
                    var line_offset = target.Line - line_min;

                    var sample_min = Math.Max(0, target.Sample - border);
                    var sample_max = Math.Min(dem_size - 1, target.Sample + patch_size + border);  // 1 more than the highest index
                    var sample_size = sample_max - sample_min;
                    var sample_offset = target.Sample - sample_min;

                    var cpu_caster_points_size = line_size * sample_size * 3;
                    var cpu_caster_points = new float[cpu_caster_points_size];
                    FillNearCasterArray(cpu_caster_points, target.Points[0][0], line_min, line_max, sample_min, sample_max);

                    //DumpNearFieldTestCsv(cpu_caster_points, cpu_matrices);

                    var test_floats_size = Horizon.HorizonSamples * TerrainPatch.NearHorizonOversample;
                    var test_floats = new float[test_floats_size];

                    using (var gpu_matrices = accelerator.Allocate<float>(cpu_matrices_size))
                    using (var gpu_horizon = accelerator.Allocate<float>(cpu_horizon_size))
                    using (var gpu_caster_points = accelerator.Allocate<float>(cpu_caster_points_size))
                    using (var gpu_test_floats = accelerator.Allocate<float>(test_floats_size))
                    {
                        gpu_matrices.CopyFrom(cpu_matrices, 0, 0, cpu_matrices_size);
                        gpu_horizon.CopyFrom(cpu_horizon, 0, 0, cpu_horizon_size);
                        gpu_caster_points.CopyFrom(cpu_caster_points, 0, 0, cpu_caster_points_size);
                        gpu_test_floats.CopyFrom(test_floats, 0, 0, test_floats_size);

                        var groupSize = accelerator.MaxNumThreadsPerGroup;
                        Index launchDimension = Horizon.HorizonSamples * TerrainPatch.NearHorizonOversample;

                        const float d_min = 1f;
                        var d_max = (float)TerrainPatch.MaximumLocalDistance;
                        const float d_step = TerrainPatch.LocalStep;

                        var kernel1 = accelerator.LoadAutoGroupedStreamKernel<Index, int, int, int, int, int, int, float, float, float, float, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>>(NearFieldKernel1);

                        //Console.WriteLine(@"Launching near horizon kernels ... ");
                        if (true)
                        {
                            var observer_height_in_km = target.ObserverHeightInMeters / 1000f;
                            for (var target_line = 0; target_line < TerrainPatch.DefaultSize; target_line++)
                                for (var target_sample = 0; target_sample < TerrainPatch.DefaultSize; target_sample++)
                                    kernel1(launchDimension, target_line, target_sample, line_offset, sample_offset, line_size, sample_size, d_min, d_max, d_step, observer_height_in_km, gpu_caster_points, gpu_matrices, gpu_horizon, gpu_test_floats);
                        }
                        else
                        {
                            kernel1(launchDimension, 0, 0, line_offset, sample_offset, line_size, sample_size, d_min, d_max, d_step, target.ObserverHeightInMeters / 1000f, gpu_caster_points, gpu_matrices, gpu_horizon, gpu_test_floats);
                        }

                        gpu_horizon.CopyTo(cpu_horizon, 0, 0, cpu_horizon_size);
                        gpu_test_floats.CopyTo(test_floats, 0, 0, test_floats_size);

                        // Update the horizons
                        for (var line = 0; line < TerrainPatch.DefaultSize; line++)
                            for (var sample = 0; sample < TerrainPatch.DefaultSize; sample++)
                            {
                                var offset = (line * TerrainPatch.DefaultSize + sample) * Horizon.HorizonSamples;
                                var buffer = target.Horizons[line][sample].Buffer;
                                for (var i = 0; i < Horizon.HorizonSamples; i++)
                                    buffer[i] = cpu_horizon[i + offset];
                            }
                    }
                }
            }
        }

        private void DumpNearFieldTestCsv(float[] cpu_caster_points, float[] matrices)
        {
            var decimate = 379;
            // Dump the points first
            using (var sw = new StreamWriter("points.csv"))
            {
                var ptr = 0;
                var ctr = 0;
                while (ptr < cpu_caster_points.Length)
                    if (0 == (ctr++) % decimate)
                        sw.WriteLine($"{cpu_caster_points[ptr++]},{cpu_caster_points[ptr++]},{cpu_caster_points[ptr++]}");
                    else
                        ptr += 3;
            }

            //var pointList = new List<Point> { new Point(29, 93) };
            var pointList = new List<Point> { new Point(29, 93), new Point(0, 0) };
            for (var i = 0; i < pointList.Count; i++)
            {
                var point = pointList[i];
                // Dump transformed points
                using (var sw = new StreamWriter($"points-transformed-{i}.csv"))
                {
                    var line = point.Y;
                    var sample = point.X;
                    var offset = (line * TerrainPatch.DefaultSize + sample) * 12;
                    var mat = new Matrix4();
                    {
                        var r = matrices;
                        var o = offset;
                        mat.Row0.X = r[o++];
                        mat.Row1.X = r[o++];
                        mat.Row2.X = r[o++];
                        mat.Row3.X = r[o++];

                        mat.Row0.Y = r[o++];
                        mat.Row1.Y = r[o++];
                        mat.Row2.Y = r[o++];
                        mat.Row3.Y = r[o++];

                        mat.Row0.Z = r[o++];
                        mat.Row1.Z = r[o++];
                        mat.Row2.Z = r[o++];
                        mat.Row3.Z = r[o++];
                    }
                    var ptr = 0;
                    var ctr = 0;
                    while (ptr < cpu_caster_points.Length)
                    {
                        if (0 == (ctr++) % decimate)
                        {
                            var v = new Vector3(cpu_caster_points[ptr++], cpu_caster_points[ptr++], cpu_caster_points[ptr++]);
                            var v1 = Vector3.Transform(v, mat);
                            sw.WriteLine($"{v1.X},{v1.Y},{v1.Z}");
                        }
                        else
                            ptr += 3;
                    }
                }
            }
        }

        internal void FillNearCasterArray(float[] ary, Vector3d base_point, int line_min, int line_max, int sample_min, int sample_max)
        {
            var ptr = 0;
            var bx = base_point.X;
            var by = base_point.Y;
            var bz = base_point.Z;
            for (var line = line_min; line < line_max; line++)
                for (var sample = sample_min; sample < sample_max; sample++)
                {
                    var relz = Terrain.LineSampleToTerrainOffset(line, sample);
                    var radius = TerrainPatch.MoonRadius + relz / 1000d;
                    InMemoryInt16Terrain.GetLatLon(line, sample, out double lat, out double lon);
                    var z = radius * Math.Sin(lat);
                    var c = radius * Math.Cos(lat);
                    var x = c * Math.Cos(lon);
                    var y = c * Math.Sin(lon);

                    ary[ptr++] = (float)(x - bx);
                    ary[ptr++] = (float)(y - by);
                    ary[ptr++] = (float)(z - bz);
                }
        }

        // This kernel runs very quickly and is called 128 x 128 times.   128 x 128 runs takes 8 sec ... longer than the 5 sec timeout.
        // It would be better if more work were done in the kernel and the kernel were called fewer times, I think.
        static void NearFieldKernel1(
            Index index,       // Horizon angle as integer
            int target_line,   // row within the 128 x 128 target patch
            int target_sample, // column
            int line_offset,
            int sample_offset,
            int points_rows,
            int points_columns,
            float d_min,
            float d_max,
            float d_step,
            float observer_height_in_km,
            ArrayView<float> points,
            ArrayView<float> matrices,
            ArrayView<float> horizon,
            ArrayView<float> test_floats)
        {
            //if (index == 0)
            //{
            //    var delme = 1;
            //}

            // Copy the matrix into registers
            var pos = (target_line * TerrainPatch.DefaultSize + target_sample) * 12;

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

            // Work out the location of the central point in the point array
            float center_line = target_line + line_offset;
            float center_sample = target_sample + sample_offset;

            var caster_line_max = (float)points_rows - 1;
            var caster_sample_max = (float)points_columns - 1;

            // Work out the direction vector
            var ray_angle = 3.141592653589f * 2f * index / (TerrainPatch.NearHorizonOversample * Horizon.HorizonSamples);
            var ray_cos = GPUMath.Cos(ray_angle);
            var ray_sin = GPUMath.Sin(ray_angle);

            // Work out the direction in horizon space (the horizon offset).  This duplicates the calculation of the first point.
            // I'm duplicating code rather than refactoring, because there would be a lot of values to pass
            float azimuth_rad;
            {
                var d = 1f;
                var caster_line = center_line + ray_sin * d;
                var caster_sample = center_sample + ray_cos * d;

                var x1 = (int)caster_sample;  // Calculate the caster point by interpolating between four points from the points array
                var y1 = (int)caster_line;
                int x2 = x1 + 1;
                int y2 = y1 + 1;

                var q11_offset = 3 * (y1 * points_columns + x1);  // (y1, x1);
                var q11 = new Vector3(points[q11_offset], points[q11_offset + 1], points[q11_offset + 2]);

                var q12_offset = 3 * (y2 * points_columns + x1);  // (y2, x1);
                var q12 = new Vector3(points[q12_offset], points[q12_offset + 1], points[q12_offset + 2]);

                // First interpolation across rows (line)
                var q1_line = q11 + (caster_line - y1) * (q12 - q11);

                var q21_offset = 3 * (y1 * points_columns + x2);  // (y1, x2);
                var q21 = new Vector3(points[q21_offset], points[q21_offset + 1], points[q21_offset + 2]);

                var q22_offset = 3 * (y2 * points_columns + x2);  // (y2, x2);
                var q22 = new Vector3(points[q22_offset], points[q22_offset + 1], points[q22_offset + 2]);

                // Second interpolation across rows
                var q2_line = q21 + (caster_line - y1) * (q22 - q21);

                // Interpolate across samples
                var caster = q1_line + (caster_sample - x1) * (q2_line - q1_line);

                // Break out the coordinates
                var x_patch = caster.X;
                var y_patch = caster.Y;
                var z_patch = caster.Z;

                // Transform the point to the local frame
                var x = x_patch * row0x + y_patch * row1x + z_patch * row2x + row3x;
                var y = x_patch * row0y + y_patch * row1y + z_patch * row2y + row3y;
                //var z = x_patch * row0z + y_patch * row1z + z_patch * row2z + row3z;

                azimuth_rad = GPUMath.Atan2(y, x) + GPUMath.PI;  // [0,2PI]
            }

            var normalized_azimuth = (Horizon.HorizonSamples - 1) * azimuth_rad / (2d * Math.PI);  // range [0,1)
            var horizon_offset = (int)(0.5d + normalized_azimuth);
            Debug.Assert(horizon_offset >= 0 && horizon_offset < Horizon.HorizonSamples);

            // index is which direction in [target_line,target_sample]'s horizon we're working on, but could be oversamples
            var horizon_index = (target_line * TerrainPatch.DefaultSize + target_sample) * Horizon.HorizonSamples + horizon_offset;
            var highest_slope = horizon[horizon_index];

            //test_floats[index] = ray_sin;

            for (var d = d_min; d <= d_max; d += d_step)
            {
                // Generate the location of the caster point in the points array
                var caster_line = center_line + ray_sin * d;
                var caster_sample = center_sample + ray_cos * d;

                test_floats[index] = caster_sample;
                /*{
                    var idx = index == 2160 ? 7 : 11;
                    test_floats[idx] = d;
                    test_floats[idx + 1] = ray_angle;
                    test_floats[idx + 2] = caster_line - line_offset;
                    test_floats[idx + 3] = caster_sample - sample_offset;
                }*/

                if (caster_line < 1f || caster_line > caster_line_max || caster_sample < 1f || caster_sample > caster_sample_max)
                    break;

                // Calculate the caster point by interpolating between four points from the points array
                var x1 = (int)caster_sample;
                var y1 = (int)caster_line;
                int x2 = x1 + 1;
                int y2 = y1 + 1;

                var q11_offset = 3 * (y1 * points_columns + x1);  // (y1, x1);
                var q11 = new Vector3(points[q11_offset], points[q11_offset + 1], points[q11_offset + 2]);

                var q12_offset = 3 * (y2 * points_columns + x1);  // (y2, x1);
                var q12 = new Vector3(points[q12_offset], points[q12_offset + 1], points[q12_offset + 2]);

                // First interpolation across rows (line)
                var q1_line = q11 + (caster_line - y1) * (q12 - q11);

                var q21_offset = 3 * (y1 * points_columns + x2);  // (y1, x2);
                var q21 = new Vector3(points[q21_offset], points[q21_offset + 1], points[q21_offset + 2]);

                var q22_offset = 3 * (y2 * points_columns + x2);  // (y2, x2);
                var q22 = new Vector3(points[q22_offset], points[q22_offset + 1], points[q22_offset + 2]);

                // Second interpolation across rows
                var q2_line = q21 + (caster_line - y1) * (q22 - q21);

                // Interpolate across samples
                var caster = q1_line + (caster_sample - x1) * (q2_line - q1_line);

                // Break out the coordinates
                var x_patch = caster.X;
                var y_patch = caster.Y;
                var z_patch = caster.Z;

                // Transform the point to the local frame
                var x = x_patch * row0x + y_patch * row1x + z_patch * row2x + row3x;
                var y = x_patch * row0y + y_patch * row1y + z_patch * row2y + row3y;
                var z = x_patch * row0z + y_patch * row1z + z_patch * row2z + row3z;

                // Adjust for solar array height (this is temporary, and I'm not sure we want this in the final version)
                z -= observer_height_in_km;

                var alen = GPUMath.Sqrt(x * x + y * y);
                var slope = z / alen;

                if (slope > highest_slope)
                    highest_slope = slope;

                /*if (target_line==0 && target_sample==0 && index == 0)
                {
                    test_floats[0] = slope;
                    test_floats[1] = x_patch;
                    test_floats[2] = y_patch;
                    test_floats[3] = z_patch;
                    test_floats[4] = x;
                    test_floats[5] = y;
                    test_floats[6] = z;
                }*/
            }

            horizon[horizon_index] = highest_slope;
        }
    }
}
