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

        public string DEMPath { get; set; }

        protected List<Vector3d> sunvecs = null;
        protected List<Vector3d> earthvecs = null;

        protected void Setup()
        {
            WriteHorizons = AppConfiguration.GetBoolean("WriteHorizons", WriteHorizons);
            CalculatePSRPatch = AppConfiguration.GetBoolean("CalculatePSRPatch", CalculatePSRPatch);
            CalculateSafeHavenPatch = AppConfiguration.GetBoolean("CalculateSafeHavenPatch", CalculateSafeHavenPatch);
            CalculateAverageSunPatch = AppConfiguration.GetBoolean("CalculateAverageSunPatch", CalculateAverageSunPatch);
            CalculateAverageEarthPatch = AppConfiguration.GetBoolean("CalculateAverageEarthPatch", CalculateAverageEarthPatch);
            Verbose = AppConfiguration.GetBoolean("Verbose", Verbose);

            DEMPath = AppConfiguration.Get("DEMPath", DEMPath);

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
        }

        protected void PrintSetup()
        {
            Console.WriteLine($"             WriteHorizons = {WriteHorizons}");
            Console.WriteLine($"         CalculatePSRPatch = {CalculatePSRPatch}");
            Console.WriteLine($"   CalculateSafeHavenPatch = {CalculateSafeHavenPatch}");
            Console.WriteLine($"  CalculateAverageSunPatch = {CalculateAverageSunPatch}");
            Console.WriteLine($"CalculateAverageEarthPatch = {CalculateAverageEarthPatch}");
            Console.WriteLine($"                   DEMPath = {DEMPath}");
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
                (ViperEnvironment.Processor as GPUHorizons).RunQueue(tocalculate, unloadHorizons: false, writeHorizons: true);
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
            var sunvec = new Vector3d[azimuth_steps];
            var elevation = new float[azimuth_steps];
            Array.Fill(elevation, float.MinValue);
            var interval = StudyInterval.MazaricoCycle;
            var step = TimeSpan.FromHours(2);
            for (var time = interval.Start; time < interval.Stop; time += step)
            {
                var sun = ViperEnvironment.SunPosition(time);
                patch.GetAzEl(sun, 0, 0, out float sun_azimuth_rad, out float sun_elevation_rad);
                var azimuth = (int)((sun_azimuth_rad + Math.PI) / azimuth_steps);
                if (sun_elevation_rad > elevation[azimuth])
                    sunvec[azimuth] = sun;
            }

            var bmp = new Bitmap(TerrainPatch.DefaultSize, TerrainPatch.DefaultSize, System.Drawing.Imaging.PixelFormat.Format8bppIndexed);
            var bmp_data = bmp.LockBits(new Rectangle(0, 0, TerrainPatch.DefaultSize, TerrainPatch.DefaultSize), System.Drawing.Imaging.ImageLockMode.WriteOnly, bmp.PixelFormat);
            for (var row = 0; row < TerrainPatch.DefaultSize; row++)
            {
                var rowptr = (byte*)(bmp_data.Scan0 + row * bmp_data.Stride);
                for (var col = 0; col < TerrainPatch.DefaultSize; col++)
                {
                    var horizon = patch.Horizons[row][col];
                    for (var azimuth = 0; azimuth < azimuth_steps; azimuth++)
                    {
                        var sun = sunvec[azimuth];
                        patch.GetAzEl(sun, col, row, out float sun_azimuth_rad, out float sun_elevation_rad);
                        var val = horizon.SunFraction2(sun_azimuth_rad * 180f / 3.141592653589f, sun_elevation_rad * 180f / 3.141592653589f);
                        if (val > 0f)
                            rowptr[col] = 1;
                    }
                }
            }
            bmp.UnlockBits(bmp_data);
            var palette = bmp.Palette;
            palette.Entries[0] = Color.Black;
            palette.Entries[1] = Color.Red;
            bmp.Palette = palette;
            bmp.Save(path, ImageFormat.Png);
        }

        void GenerateSafeHavenPatch(TerrainPatch patch)
        {

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
