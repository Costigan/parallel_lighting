using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace viper.corelib.spice
{
    /// <summary>
    /// Wraps spice and doesn't load it if running in 32-bit mode.
    /// </summary>
    public class SpiceManager
    {
        protected static SpiceManager _singleton;
        protected readonly SpiceMethods Methods;

        public SpiceManager()
        {
            Methods = Environment.Is64BitProcess && Environment.Is64BitOperatingSystem ? new FullSpiceMethods() : (SpiceMethods)new DummySpiceMethods();
        }

        public static SpiceManager GetSingleton() => _singleton != null ? _singleton : _singleton = new SpiceManager();

        public void SunEarthAzEl(DateTime time, out double sunAzimuth, out double sunElevation, out double earthAzimuth, out double earthElevation)
        {
            Methods.SunEarthAzEl(time, out sunAzimuth, out sunElevation, out earthAzimuth, out earthElevation);
        }

        public void DSNAzEl(DateTime time, double lat, double lon, double height, out double dsnaz, out double dsnel)
        {
            Methods.DSNAzEl(time, lat, lon, height, out dsnaz, out dsnel);
        }
    }

    public abstract class SpiceMethods
    {
        public static DateTime Epoch = new DateTime(2000, 1, 1, 11, 58, 55, 816);
        public const int EarthId = 399;
        public const int MoonId = 301;
        public const int SunId = 10;
        public const double MoonRadius = 1737.4d;

        public abstract void SunEarthAzEl(DateTime time, out double sunAzimuth, out double sunElevation, out double earthAzimuth, out double earthElevation);
        public abstract void DSNAzEl(DateTime time, double lat, double lon, double height, out double dsnaz, out double dsnel);

        public static double DateTimeToET(DateTime time)
        {
            return (time - Epoch).TotalSeconds + 3d;
            // The 3D accounts for leap seconds since 2000.  This is valid only for dates after Jul 1 2012.
        }
    }

    public class DummySpiceMethods : SpiceMethods
    {
        public DummySpiceMethods() { }

        public override void DSNAzEl(DateTime time, double lat, double lon, double height, out double dsnaz, out double dsnel)
        {
            dsnaz = dsnel = 0d;
        }

        public override void SunEarthAzEl(DateTime time, out double sunAzimuth, out double sunElevation, out double earthAzimuth, out double earthElevation)
        {
            sunAzimuth = sunElevation = earthAzimuth = earthElevation = 0d;
        }
    }

    public class FullSpiceMethods : SpiceMethods
    {
        public FullSpiceMethods()
        {
            FurnishSpiceKernels("StaticFiles", "kernels/metakernel.txt");
        }

        public void FurnishSpiceKernels(string spiceKernelRoot, string rootFile)
        {
            string metakernelPath = Path.Combine(spiceKernelRoot, rootFile);
            foreach (string line in File.ReadAllLines(metakernelPath).Where(line => !string.IsNullOrEmpty(line) && line[0] != ' ' && line[0] != '#'))
                Furnish(Path.Combine(spiceKernelRoot, CanonicalizeDirectorySeparators(line)));
        }

        public static void Furnish(string path)
        {
            Console.WriteLine(@"Spice furnishing {0}", path);
            CSpice.furnsh_c(path);
        }

        public static string CanonicalizeDirectorySeparators(string path)
        {
            path = path.Replace('/', Path.DirectorySeparatorChar);
            path = path.Replace('\\', Path.DirectorySeparatorChar);
            return path;
        }

        public override void SunEarthAzEl(DateTime time, out double sunAzimuth, out double sunElevation, out double earthAzimuth, out double earthElevation)
        {
            //sunAzimuth = sunElevation = earthAzimuth = earthElevation = 0d;
            //return;

            var et = DateTimeToET(time);
            var state = new double[6];
            double lt = 0d;
            double[] pos = { 0d, 0d, 0d };

            CSpice.spkgeo_c(SunId, et, "MOON_ME", MoonId, state, ref lt);

            CSpice.spkcpo_c("SUN", et, "SITE_TOPO", "CENTER", "NONE", pos, "MOON", "MOON_ME", state, ref lt);
            //var d = Math.Sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
            //Console.WriteLine(@"HERMITE: Sun= dist={3} vec=[{0},{1},{2}]", state[0], state[1], state[2], d);

            sunAzimuth = -Math.Atan2(state[1], state[0]);
            var flatd = Math.Sqrt(state[0] * state[0] + state[1] * state[1]);
            sunElevation = Math.Atan2(state[2], flatd);

            CSpice.spkcpo_c("EARTH", et, "SITE_TOPO", "CENTER", "NONE", pos, "MOON", "MOON_ME", state, ref lt);
            //d = Math.Sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
            //Console.WriteLine(@"HERMITE: Earth= dist={3} vec=[{0},{1},{2}]", state[0], state[1], state[2], d);

            earthAzimuth = -Math.Atan2(state[1], state[0]);
            flatd = Math.Sqrt(state[0] * state[0] + state[1] * state[1]);
            earthElevation = Math.Atan2(state[2], flatd);
        }

        public string[] DSNSiteNames = { "DSS-54", "DSS-24", "DSS-34" };
        public string[] DSNSiteFrames = { "DSS-54_TOPO", "DSS-24_TOPO", "DSS-34_TOPO" };
        public string[] DSNComplexNames = { "Madrid", "Goldstone", "Canberra" };
        public int[] DSNIds = { 399054, 399024, 399034 };

        public override void DSNAzEl(DateTime time, double lat, double lon, double height, out double dsnaz, out double dsnel)
        {
            const double stationHorizonMask = 10d * Math.PI / 190d;

            var siteLat = lat * Math.PI / 180d;
            var siteLon = lon * Math.PI / 180d;

            double[] traverseSite = { 0d, 0d, 0d };
            CSpice.pgrrec_c("MOON", siteLon, siteLat, 0d, MoonRadius, 0f, traverseSite);

            var pz = MoonRadius * Math.Sin(siteLat);
            var temp = MoonRadius * Math.Cos(siteLat);
            var px = temp * Math.Cos(siteLon);
            var py = temp * Math.Sin(siteLon);

            Debug.Assert(Math.Abs(px - traverseSite[0]) < 1d);
            Debug.Assert(Math.Abs(py - traverseSite[1]) < 1d);
            Debug.Assert(Math.Abs(pz - traverseSite[2]) < 1d);

            var et = DateTimeToET(time);
            var state = new double[6];
            var earthState = new double[6];
            var moonPositionFromStation = new double[3];
            double lt = 0d;

            CSpice.spkcpo_c("EARTH", et, "SITE_TOPO", "CENTER", "NONE", traverseSite, "MOON", "MOON_ME", earthState, ref lt);
            var earthAzimuth = -Math.Atan2(earthState[1], earthState[0]);
            var flatd1 = Math.Sqrt(earthState[0] * earthState[0] + earthState[1] * earthState[1]);
            var earthElevation = Math.Atan2(earthState[2], flatd1);

            var data = new List<DSNSiteData>();
            for (var i = 0; i < DSNSiteNames.Length; i++)
            {
                // Station to the traverse site
                CSpice.spkezp_c(MoonId, et, DSNSiteFrames[i], "NONE", DSNIds[i], moonPositionFromStation, ref lt);
                var moonElevation = Math.Atan2(moonPositionFromStation[2], Math.Sqrt(moonPositionFromStation[0] * moonPositionFromStation[0] + moonPositionFromStation[1] * moonPositionFromStation[1]));

                // Traverse site to the station
                CSpice.spkcpo_c(DSNSiteNames[i], et, "SITE_TOPO", "CENTER", "NONE", traverseSite, "MOON", "MOON_ME", state, ref lt);

                var flatd = Math.Sqrt(state[0] * state[0] + state[1] * state[1]);
                data.Add(new DSNSiteData
                {
                    SiteToStationAzimuth = -Math.Atan2(state[1], state[0]),
                    SiteToStationElevation = Math.Atan2(state[2], flatd),
                    StationToMoonElevation = moonElevation,
                    Site = DSNSiteNames[i]
                });
            }

            foreach (var d in data)
                Console.WriteLine("complex={0} elev@complex={1}  complexElev@traverse={2} complexAz@traverse={3}", d.Site, d.StationToMoonElevation, d.SiteToStationElevation, d.SiteToStationAzimuth);

            DSNSiteData bestData = null;
            const double highestElev = double.MinValue;
            foreach (var d in data)
            {
                if (d.StationToMoonElevation >= stationHorizonMask && d.SiteToStationElevation > highestElev)
                    bestData = d;
            }

            if (bestData == null)
            {
                Console.WriteLine(@"Couldn't find a suitable ground complex.  Using {0}", data[0].Site);
                bestData = data[0];
            }
            dsnaz = bestData.SiteToStationAzimuth;
            dsnel = bestData.SiteToStationElevation;

            Console.WriteLine(@"moon center");
            Console.WriteLine(@"  earth az={0}", earthAzimuth);
            Console.WriteLine(@"  earth el={0}", earthElevation);
            Console.WriteLine(@"To {0}", bestData.Site);
            Console.WriteLine(@"  earth az={0}", bestData.SiteToStationAzimuth);
            Console.WriteLine(@"  earth el={0}", bestData.SiteToStationElevation);
        }

        internal class DSNSiteData
        {
            internal double SiteToStationAzimuth;
            internal double SiteToStationElevation;
            internal double StationToMoonElevation;
            internal string Site;
        }
    }
}

