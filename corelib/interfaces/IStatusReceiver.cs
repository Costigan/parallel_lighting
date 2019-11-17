using System;
using System.Collections.Generic;
using System.Text;
using viper.corelib.patch;

namespace viper.corelib.interfaces
{
    public interface IStatusReceiver
    {
        bool IsRunning { get; set; }
        List<TerrainPatch> ProcessingPatches { get; set; }
        List<TerrainPatch> FarPatches { get; set; }
    }
}
