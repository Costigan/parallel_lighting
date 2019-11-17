using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using viper.corelib.horizon;
using viper.corelib.patch;

namespace viper.corelib.utilities
{
    public class ViperEnvironment
    {
        public static int MaxDegreeOfParallelism = Environment.ProcessorCount - 2;
        public static terrain.InMemoryInt16Terrain Terrain;

        public const int TerrainWidth20m = 30400;

        public static string TileRoot = AppConfiguration.Get("TileRoot", @"c:/RP/tiles");
        public static string SouthPsrImgPath = AppConfiguration.Get("SouthPsrImgPath", @"c:/RP/LOLA/LPSR_85S_060M_201608_DEM_byte.img");
        public static string NorthPsrImgPath = AppConfiguration.Get("NorthPsrImgPath", @"c:/RP/LOLA/LPSR_85N_060M_201608_DEM_byte.img");
        public static string SouthAvgSunImgPath = AppConfiguration.Get("SouthAvgSunImgPath", @"c:/RP/LOLA/AVGVISIB_85S_060M_201608_DEM.img");
        public static string NorthAvgSunImgPath = AppConfiguration.Get("NorthAvgSunImgPath", @"c:/RP/LOLA/AVGVISIB_85N_060M_201608_DEM.img");
        public static string SouthAvgEarthImgPath = AppConfiguration.Get("SouthAvgEarthImgPath", @"c:/RP/LOLA/AVGVISIB_85S_060M_201608_DEM_EARTH.img");
        public static string NorthAvgEarthImgPath = AppConfiguration.Get("NorthAvgEarthImgPath", @"c:/RP/LOLA/AVGVISIB_85N_060M_201608_DEM_EARTH.img"); // not right yet


        public static bool IsNorth => Terrain.IsNorth;
        public static string PoleDirectory => Terrain.IsNorth ? "np" : "sp";

        public static string MapRoot = AppConfiguration.Get("MapRoot", @"c:/RP/tiles/np");
        public static string HorizonsRoot = AppConfiguration.Get("HorizonsRoot", @"c:/RP/tiles/np/horizons");
        public static string AvgPatchRoot = AppConfiguration.Get("AvgPatchRoot", @"c:/RP/tiles/np/avg_patches");
        public static string HavenPatchRoot = AppConfiguration.Get("HavenPatchRoot", @"c:/RP/tiles/np/haven_patches");

        public const string HorizonDirectoryName = @"horizons";
        public const string AvgPatchDirectoryName = @"avg_patches";
        public const string HavenPatchDirectoryName = @"haven_patches";

        #region Patch Cache

        public static Dictionary<(int, int, int), TerrainPatch> PatchCache = new Dictionary<(int, int, int), TerrainPatch>();

        public static IEnumerable<TerrainPatch> LoadPatches(IEnumerable<Point> ids) => ids.Select(id => GetPatch(id)).Where(p => p != null);

        public static TerrainPatch GetPatch(Point id, int observer = 0)
        {
            if (PatchCache.TryGetValue((id.X, id.Y, observer), out TerrainPatch patch))
                return patch;
            var patch2 = TerrainPatch.FromId(id, observer);
            if (File.Exists(patch2.Path))
            {
                Console.WriteLine($"Reading patch {patch2.Path}");
                var patch3 = TerrainPatch.ReadFrom(patch2.Path);
                PatchCache[(id.X, id.Y, observer)] = patch3;
                return patch3;
            }
            return patch2;
        }

        public static TerrainPatch GetPatch(TerrainPatch patch, int observer = 0)
        {
            var id = patch.Id;
            if (PatchCache.TryGetValue((id.X, id.Y, observer), out TerrainPatch patch2))
                return patch2;
            if (patch.Horizons != null && patch.Horizons[0][0].PartialPoints().Count() > 0)
            {
                PatchCache[(id.X, id.Y, observer)] = patch;
                return patch;
            }
            if (File.Exists(patch.Path))
            {
                Console.WriteLine($"Reading patch {patch.Path}");
                var patch3 = TerrainPatch.ReadFrom(patch.Path);
                PatchCache[(id.X, id.Y, observer)] = patch3;
                return patch3;
            }
            return null;
        }

        #endregion

        #region GPU-related

        public static CPUHorizons Processor;

        #endregion

        #region Drawing and Fonts

        public static Font TimestampFont = new Font("Arial", 32, FontStyle.Regular);
        public static Font DefaultOverlayFont = new Font("Arial", 16, FontStyle.Regular);

        #endregion
    }
}
