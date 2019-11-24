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
            Console.WriteLine("Hello from netcore3.0");

            // Initialization
            OSGeo.GDAL.Gdal.AllRegister();
            ViperEnvironment.Terrain = new InMemoryInt16Terrain();
            ViperEnvironment.Terrain.LoadSouth();
            ViperEnvironment.Processor = new GPUHorizons { Terrain = ViperEnvironment.Terrain };

            var calc = new LightingCalculator{ WriteHorizons = true, Verbose = true, CalculatePSRPatch = true, CalculateAverageSunPatch = true, CalculateAverageEarthPatch = true };
            calc.Run();
        }
    }
}
