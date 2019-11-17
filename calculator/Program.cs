using Microsoft.Extensions.Configuration;
using System;
using System.Collections.Generic;
using System.IO;
using viper.corelib.horizon;
using viper.corelib.lighting;
using viper.corelib.patch;
using viper.corelib.spice;
using viper.corelib.terrain;
using viper.corelib.utilities;

namespace calculator
{
    class Program
    {
        static void Main(string[] args)
        {
            new Test2().Run();
        }
    }

    public class BaseTest
    {
        public static SpiceManager Spice;
        public static string OutputPath;
        public virtual void Run()
        {
            Init();
        }

        public void Init()
        {
            // Initialize paths and create directories if needed
            OutputPath = AppConfiguration.Get("OutputPath", ".");
            if (!Directory.Exists(OutputPath)) Directory.CreateDirectory(OutputPath);
            ViperEnvironment.HorizonsRoot = AppConfiguration.Get("HorizonsRoot", ".");
            if (!Directory.Exists(ViperEnvironment.HorizonsRoot)) Directory.CreateDirectory(ViperEnvironment.HorizonsRoot);

            Spice = SpiceManager.GetSingleton();

            // Load South;
            ViperEnvironment.Terrain = new InMemoryInt16Terrain();
            ViperEnvironment.Terrain.LoadSouth();

            ViperEnvironment.Processor= new GPUHorizons { Terrain = ViperEnvironment.Terrain };
        }
    }

    public class Test0 : BaseTest
    {
        public override void Run()
        {
            base.Run();
            var patch = TerrainPatch.FromId(147, 72);
            var queue = new List<TerrainPatch> { patch };
            var gpu = new GPUHorizons { Terrain = ViperEnvironment.Terrain };
            Console.WriteLine("Running the queue");
            gpu.RunQueue(queue, unloadHorizons: false, writeHorizons: false);
            Console.WriteLine("Finished running the queue");
        }
    }

    public class Test1 : BaseTest
    {
        public override void Run()
        {
            base.Run();
            var patch = TerrainPatch.FromId(147, 72);
            var safe_haven_gen = new SafeHavenGenerator { Region = patch.Bounds, WriteHorizons = true };
            safe_haven_gen.WriteSafeHavenGeotiffs(OutputPath);
        }
    }

    public class Test2 : BaseTest
    {
        public override void Run()
        {
            base.Run();

            var patch = TerrainPatch.FromId(147, 72);

            Console.WriteLine(@"Generating safe havens");
            var safe_haven_gen = new SafeHavenGenerator { Region = patch.Bounds, WriteHorizons = false };
            safe_haven_gen.WriteSafeHavenGeotiffs(OutputPath);

            Console.WriteLine(@"Generating average Earth and average Sun patches");
            var (start, stop, step) = StudyInterval.SiteStudyPeriod.GetInterval();
            var gen = new TileLightingProductManager
            {
                //MainWindow = MainWindow,
                Selection = patch.Bounds,
                IntervalStart = start,
                IntervalStop = stop,
                IntervalStep = step,
                ObserverHeightInMeters = 0,
                EarthMultipathThreshold = 2f,
                GenerateAverageSun = true,
                GenerateAverageEarth = true
            };
            gen.GenerateAverageSunEarthPatches();

            Console.WriteLine(@"Finished.");
        }
    }
}
