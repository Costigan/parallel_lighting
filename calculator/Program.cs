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



            ViperEnvironment.SetMapRoot("./sp/");

            new LightingCalculator {
                Verbose = true,
                WriteHorizons = true,
                CalculatePSRPatch = true,
                CalculateSafeHavenPatch = true,
                CalculateAverageSunPatch = false,
                CalculateAverageEarthPatch = false
            }
            .Run();
        }
    }
}
