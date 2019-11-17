using System;
using System.Drawing;
using System.IO;
using System.IO.MemoryMappedFiles;
using System.Runtime.CompilerServices;
using viper.corelib.math;
using viper.corelib.utilities;

// I'm working on handling multiple terrain formats and memory management approaches in a uniform way.
// However, the changes needed for that aren't done yet.

namespace viper.corelib.terrain
{
    public abstract class TerrainManager : IDisposable
    {
        public static TerrainManager Factory(bool isNorth = false, bool inMemory = true, bool isFloat = false)
        {
            if (inMemory && isFloat)
                return isNorth ? new InMemorySingleTerrain().LoadNorth() : new InMemorySingleTerrain().LoadSouth();
            if (inMemory && !isFloat)
                return isNorth ? new InMemoryInt16Terrain().LoadNorth() : new InMemoryInt16Terrain().LoadSouth();
            if (!inMemory && isFloat)
                throw new NotImplementedException();
            if (!inMemory && !isFloat)
                throw new NotImplementedException();
            throw new ArgumentException();
        }

        public abstract void Dispose();
    }

    public class InMemoryInt16Terrain : TerrainManager
    {
        //private const string TerrainImageFile = @"C:\UVS\svn\src\TestData\Terrain\ldem_64.img";
        public static string NorthDEM = AppConfiguration.Get("NorthDEM20", @"C:\RP\maps\DEM\ldem_80n_20m.horizon.img");

        //public static string SouthDEM = @"C:\RP\maps\DEM\ldem_80s_20m.img";  // Not right yet
        public static string SouthDEM = AppConfiguration.Get("SouthDEM20", @"C:\RP\maps\DEM\ldem_80s_20m_fixed.horizon.img");  // Not right yet

        public static InMemoryInt16Terrain Singleton = null;

        public short[] Data { get; private set; }

        // Conversion to lat/lon
        public const double S0 = 15199.5d;             // PDS SAMPLE_PROJECTION_OFFSET
        public const double L0 = 15199.5d;             // PDS LINE_PROJECTION_OFFSET
        public static double LonFactor = 1d;           // appropriate for south
        public const double Scale = 20d / 1000d;
        public const double LonP = 0d;
        public static double LatP = -Math.PI / 2;      // appropriate for south
        public const double MoonRadius = 1737.4d;
        public const double RadiusInMeters = 1737400.0d;
        public const int Samples = 30400;
        public bool IsNorth = false;

        // Constants

        public virtual InMemoryInt16Terrain LoadNorth()
        {
            IsNorth = true;
            LonFactor = -1d;
            LatP = Math.PI / 2;
            Open(NorthDEM);
            Singleton = this;
            return this;
        }

        public virtual InMemoryInt16Terrain LoadSouth()
        {
            IsNorth = false;
            LonFactor = 1d;
            LatP = -Math.PI / 2;
            Open(SouthDEM);
            Singleton = this;
            return this;
        }

        public static bool IsInitialized => Singleton != null;

        protected short _min = short.MaxValue;
        protected short _max = short.MinValue;
        protected bool _minMaxCalculated;

        public unsafe void Open(string lola_dem_path)
        {
            Console.WriteLine($"Loading terrain {lola_dem_path}");
            var fi = new FileInfo(lola_dem_path);
            var byte_count = fi.Length;
            Data = new short[byte_count / 2];
            fixed (short* shortptr = &Data[0])
            {
                var ptr = (byte*)shortptr;
                using (var memory_stream = new UnmanagedMemoryStream(ptr, 0, byte_count, FileAccess.Write))
                using (var file_stream = new FileStream(lola_dem_path, FileMode.Open, FileAccess.Read))
                    file_stream.CopyTo(memory_stream);
            }
        }

        public void Close() => Dispose(true);

        public short Min
        {
            get
            {
                if (_minMaxCalculated) return _min;
                GetMinMax();
                return _min;
            }
        }

        public short Max
        {
            get
            {
                if (_minMaxCalculated) return _max;
                GetMinMax();
                return _max;
            }
        }

        /// <summary>
        /// Return latitude and longitude in radians
        /// </summary>
        /// <param name="line"></param>
        /// <param name="sample"></param>
        /// <param name="latitude_rad"></param>
        /// <param name="longitude_rad"></param>
        /*public static void GetLatLon(int line, int sample, out double latitude, out double longitude)
        {
            const double R = MoonRadius;
            var x = (sample - S0 - 1d) * Scale;
            var y = (1d - L0 - line) * Scale;
            var P = Math.Sqrt(x * x + y * y);
            var C = 2 * Math.Atan2(P, 2 * R);
            latitude = Math.Asin(Math.Cos(C) * Math.Sin(LatP) + y * Math.Sin(C) * Math.Cos(LatP) / P);
            longitude = LonP + Math.Atan2(x, y * LonFactor);
        }*/

        /// <summary>
        /// Return latitude and longitude in radians
        /// </summary>
        /// <param name="line"></param>
        /// <param name="sample"></param>
        /// <param name="latitude"></param>
        /// <param name="longitude"></param>
        public static void GetLatLon(double line, double sample, out double latitude_rad, out double longitude_rad)
        {
            const double R = MoonRadius;
            var x = (sample - S0) * Scale;
            var y = (L0 - line) * Scale;
            var P = Math.Sqrt(x * x + y * y);
            var C = 2d * Math.Atan2(P, 2 * R);
            latitude_rad = Math.Asin(Math.Cos(C) * Math.Sin(LatP) + (y == 0 ? 0 : y * Math.Sin(C) * Math.Cos(LatP) / P));
            longitude_rad = LonP + Math.Atan2(x, y * LonFactor);
            //Console.WriteLine($"Line={line} Sample={sample} x={x} y={y} lat={latitude} lon={longitude}");
        }

        public static void GetLatLonDegrees(double line, double sample, out double latitude, out double longitude)
        {
            GetLatLon(line, sample, out double lat, out double lon);
            latitude = lat * 180d / Math.PI;
            longitude = lon * 180d / Math.PI;
        }

        /// <summary>
        /// Get Lat/Lon for a point in moon_me frame
        /// </summary>
        /// <param name="p">Position in MOON_ME frame</param>
        /// <param name="lat">Latitude of p in radians</param>
        /// <param name="lon">Longitude of p in radians [-Pi,Pi]</param>
        public static void MoonMEToLatLon(Vector3d p, out double lat, out double lon)
        {
            var d = Math.Sqrt(p.X * p.X + p.Y * p.Y);
            lat = Math.Atan2(p.Z, d);
            lon = Math.Atan2(p.Y, p.X);
        }

        public static void MoonMEToLatLonDeg(Vector3d p, out double lat, out double lon)
        {
            MoonMEToLatLon(p, out lat, out lon);
            lat *= 180d / Math.PI;
            lon *= 180d / Math.PI;
        }

        public Vector3d LineSampleHeightToMoonME(int line, int sample, float height)
        {
            var radius = MoonRadius + height / 1000d;
            GetLatLon(line, sample, out double lat, out double lon);
            var z = radius * Math.Sin(lat);
            var c = radius * Math.Cos(lat);
            var x = c * Math.Cos(lon);  // TODO: Not sure about these
            var y = c * Math.Sin(lon);
            return new Vector3d(x, y, z);
        }

        public static Point GetLineSampleDegrees(double lat_deg, double lon_deg)
        {
            GetLineSampleDegrees(lat_deg, lon_deg, out int line, out int sample);
            return new Point(sample, line);
        }

        public static void GetLineSampleDegrees(double lat_deg, double lon_deg, out int line, out int sample) =>
            GetLineSample(lat_deg * Math.PI / 180d, lon_deg * Math.PI / 180d, out line, out sample);

        // Checked against GetLatLon
        public static void GetLineSample(double lat_rad, double lon_rad, out int line, out int sample)
        {
            double x, y;
            if (InMemoryInt16Terrain.Singleton.IsNorth)
            {
                var temp1 = MoonRadius * Math.Tan(Math.PI / 4d - lat_rad / 2d);
                var temp2 = lon_rad - LonP;
                x = 2d * temp1 * Math.Sin(temp2);
                y = -2d * temp1 * Math.Cos(temp2);
            }
            else
            {
                var temp1 = MoonRadius * Math.Tan(Math.PI / 4d + lat_rad / 2d);
                var temp2 = lon_rad - LonP;
                x = 2d * temp1 * Math.Sin(temp2);
                y = 2d * temp1 * Math.Cos(temp2);
            }

            //x = (Sample - S0 - 1) * Scale
            //y = (1 - L0 - Line) * Scale => Line = 1 - L0 - y/Scale
            sample = (int)(x / Scale + S0 + 1d + 0.5d);

            // line = (int)(1d - L0 - y / Scale);   // Should be this per the notes
            line = -(int)(1d - L0 - y / Scale);  //NOTE: I'm not sure why there's a - sign here

            // I don't understand this correction
            const int mid = 30400 / 2;
            line = mid - (line - mid);
        }

        public static void GetLineSampleD(double lat_rad, double lon_rad, out double line, out double sample)
        {
            double x, y;
            if (InMemoryInt16Terrain.Singleton.IsNorth)
            {
                var temp1 = MoonRadius * Math.Tan(Math.PI / 4d - lat_rad / 2d);
                var temp2 = lon_rad - LonP;
                x = 2d * temp1 * Math.Sin(temp2);
                y = -2d * temp1 * Math.Cos(temp2);
            }
            else
            {
                var temp1 = MoonRadius * Math.Tan(Math.PI / 4d + lat_rad / 2d);
                var temp2 = lon_rad - LonP;
                x = 2d * temp1 * Math.Sin(temp2);
                y = 2d * temp1 * Math.Cos(temp2);
            }

            //x = (Sample - S0 - 1) * Scale
            //y = (1 - L0 - Line) * Scale => Line = 1 - L0 - y/Scale
            sample = x / Scale + S0;

            // line = (int)(1d - L0 - y / Scale);   // Should be this per the notes
            //line = -(1d - L0 - y / Scale);  //NOTE: I'm not sure why there's a - sign here
            line = -y / Scale + L0;

            // I don't understand this correction
            //const int mid = 30400 / 2;
            //line = mid - (line - mid);
        }

        public float InterpolatedElevation(float line, float sample)
        {
            var x1 = (int)sample;
            var y1 = (int)line;
            int x2 = x1 + 1;
            int y2 = y1 + 1;
            float q11 = LineSampleToTerrainOffset(y1, x1);
            float q12 = LineSampleToTerrainOffset(y2, x1);
            float q21 = LineSampleToTerrainOffset(y1, x2);
            float q22 = LineSampleToTerrainOffset(y2, x2);
            // From https://en.wikipedia.org/wiki/Bilinear_interpolation, note denominator is 1 in this case, so is omitted
            float r = (q11 * (x2 - sample) * (y2 - line) + q21 * (sample - x1) * (y2 - line) + q12 * (x2 - sample) * (line - y1) + q22 * (sample - x1) * (line - y1));
            //Console.WriteLine(@"  x={0:G7} y={1:G7} x1={2} y1={3} x2={4} y2={5} q11={6:G7} q12={7:G7} q21={8:G7} q22={9:G7}", x, y, x1, y1, x2, y2, q11, q12, q21, q22);
            return r;
        }

        public Vector3d GetPointInME(int line, int sample)
        {
            var relz = LineSampleToTerrainOffset(line, sample);
            var radius = MoonRadius + relz / 1000d;
            GetLatLon(line, sample, out double lat, out double lon);
            var z = radius * Math.Sin(lat);
            var c = radius * Math.Cos(lat);
            var x = c * Math.Cos(lon);  // TODO: Not sure about these
            var y = c * Math.Sin(lon);
            return new Vector3d(x, y, z);
        }

        public Vector3d GetPointInME(float line, float sample)
        {
            var relz = InterpolatedElevation(line, sample);
            var radius = MoonRadius + relz / 1000d;
            GetLatLon(line, sample, out double lat, out double lon);
            var z = radius * Math.Sin(lat);
            var c = radius * Math.Cos(lat);
            var x = c * Math.Cos(lon);  // TODO: Not sure about these
            var y = c * Math.Sin(lon);
            return new Vector3d(x, y, z);
        }

        public Vector3d GetPointInME(PointF p) => GetPointInME(p.Y, p.X);

        public float LineSampleToTerrainOffset(int line, int sample)
        {
            const int height = 30400;
            const int width = 30400;
            const int stride = width;

            System.Diagnostics.Debug.Assert(line >= 0 && line < height && sample >= 0 && sample < width);
            //if (!(line >= 0 && line < height && sample >= 0 && sample < width))
            //    Console.WriteLine("here");

            int index = line * stride + sample;
            short h = Data[index];
            return 0.5f * h;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public short Fetch(int line, int sample) => Data[line * 30400 + sample];

        public void GetMinMax()
        {
            var min = short.MaxValue;
            var max = short.MinValue;
            var length = Data.Length;
            for (var i = 0; i < length; i++)
            {
                var v = Data[i];
                if (v < min) min = v;
                if (v > max) max = v;
            }
            _max = max;
            _min = min;
        }

        #region IDisposable Support
        private bool disposedValue = false; // To detect redundant calls

        protected virtual void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (disposing)
                {
                    // TODO: dispose managed state (managed objects).
                }

                // TODO: free unmanaged resources (unmanaged objects) and override a finalizer below.
                Data = null;

                disposedValue = true;
            }
        }

        // TODO: override a finalizer only if Dispose(bool disposing) above has code to free unmanaged resources.
        // ~InMemoryTerrainManager() {
        //   // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
        //   Dispose(false);
        // }

        // This code added to correctly implement the disposable pattern.
        public override void Dispose()
        {
            // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
            Dispose(true);
            // TODO: uncomment the following line if the finalizer is overridden above.
            // GC.SuppressFinalize(this);
        }
        #endregion
    }

    public class InMemorySingleTerrain : TerrainManager
    {
        //private const string TerrainImageFile = @"C:\UVS\svn\src\TestData\Terrain\ldem_64.img";
        public static string NorthDEM = utilities.AppConfiguration.Get("NorthDEM20", @"C:\RP\maps\DEM\ldem_80n_20m.horizon.img");

        //public static string SouthDEM = @"C:\RP\maps\DEM\ldem_80s_20m.img";  // Not right yet
        public static string SouthDEM = utilities.AppConfiguration.Get("SouthDEM20", @"C:\RP\maps\DEM\ldem_80s_20m_fixed.horizon.img");  // Not right yet

        public static InMemorySingleTerrain Singleton = null;

        public float[] Data { get; private set; }

        // Conversion to lat/lon
        public const double S0 = 15199.5d;             // PDS SAMPLE_PROJECTION_OFFSET
        public const double L0 = 15199.5d;             // PDS LINE_PROJECTION_OFFSET
        public static double LonFactor = 1d;           // appropriate for south
        public const double Scale = 20d / 1000d;
        public const double LonP = 0d;
        public static double LatP = -Math.PI / 2;      // appropriate for south
        public const double MoonRadius = 1737.4d;
        public const double RadiusInMeters = 1737400.0d;
        public const int Samples = 30400;
        public bool IsNorth = false;

        // Constants

        public virtual InMemorySingleTerrain LoadNorth()
        {
            IsNorth = true;
            LonFactor = -1d;
            LatP = Math.PI / 2;
            Open(NorthDEM);
            Singleton = this;
            return this;
        }

        public virtual InMemorySingleTerrain LoadSouth()
        {
            IsNorth = false;
            LonFactor = 1d;
            LatP = -Math.PI / 2;
            Open(SouthDEM);
            Singleton = this;
            return this;
        }

        public static bool IsInitialized => Singleton != null;

        protected float _min = short.MaxValue;
        protected float _max = short.MinValue;
        protected bool _minMaxCalculated;

        public unsafe void Open(string lola_dem_path)
        {
            Console.WriteLine($"Loading terrain {lola_dem_path}");
            var fi = new FileInfo(lola_dem_path);
            var byte_count = fi.Length;
            Data = new float[byte_count / 4];
            fixed (float* shortptr = &Data[0])
            {
                var ptr = (byte*)shortptr;
                using (var memory_stream = new UnmanagedMemoryStream(ptr, 0, byte_count, FileAccess.Write))
                using (var file_stream = new FileStream(lola_dem_path, FileMode.Open, FileAccess.Read))
                    file_stream.CopyTo(memory_stream);
            }
        }

        public void Close() => Dispose(true);

        public float Min
        {
            get
            {
                if (_minMaxCalculated) return _min;
                GetMinMax();
                return _min;
            }
        }

        public float Max
        {
            get
            {
                if (_minMaxCalculated) return _max;
                GetMinMax();
                return _max;
            }
        }

        /// <summary>
        /// Return latitude and longitude in radians
        /// </summary>
        /// <param name="line"></param>
        /// <param name="sample"></param>
        /// <param name="latitude_rad"></param>
        /// <param name="longitude_rad"></param>
        /*public static void GetLatLon(int line, int sample, out double latitude, out double longitude)
        {
            const double R = MoonRadius;
            var x = (sample - S0 - 1d) * Scale;
            var y = (1d - L0 - line) * Scale;
            var P = Math.Sqrt(x * x + y * y);
            var C = 2 * Math.Atan2(P, 2 * R);
            latitude = Math.Asin(Math.Cos(C) * Math.Sin(LatP) + y * Math.Sin(C) * Math.Cos(LatP) / P);
            longitude = LonP + Math.Atan2(x, y * LonFactor);
        }*/

        /// <summary>
        /// Return latitude and longitude in radians
        /// </summary>
        /// <param name="line"></param>
        /// <param name="sample"></param>
        /// <param name="latitude"></param>
        /// <param name="longitude"></param>
        public static void GetLatLon(double line, double sample, out double latitude_rad, out double longitude_rad)
        {
            const double R = MoonRadius;
            var x = (sample - S0) * Scale;
            var y = (L0 - line) * Scale;
            var P = Math.Sqrt(x * x + y * y);
            var C = 2d * Math.Atan2(P, 2 * R);
            latitude_rad = Math.Asin(Math.Cos(C) * Math.Sin(LatP) + (y == 0 ? 0 : y * Math.Sin(C) * Math.Cos(LatP) / P));
            longitude_rad = LonP + Math.Atan2(x, y * LonFactor);
            //Console.WriteLine($"Line={line} Sample={sample} x={x} y={y} lat={latitude} lon={longitude}");
        }

        public static void GetLatLonDegrees(double line, double sample, out double latitude, out double longitude)
        {
            GetLatLon(line, sample, out double lat, out double lon);
            latitude = lat * 180d / Math.PI;
            longitude = lon * 180d / Math.PI;
        }

        /// <summary>
        /// Get Lat/Lon for a point in moon_me frame
        /// </summary>
        /// <param name="p">Position in MOON_ME frame</param>
        /// <param name="lat">Latitude of p in radians</param>
        /// <param name="lon">Longitude of p in radians [-Pi,Pi]</param>
        public static void MoonMEToLatLon(Vector3d p, out double lat, out double lon)
        {
            var d = Math.Sqrt(p.X * p.X + p.Y * p.Y);
            lat = Math.Atan2(p.Z, d);
            lon = Math.Atan2(p.Y, p.X);
        }

        public static void MoonMEToLatLonDeg(Vector3d p, out double lat, out double lon)
        {
            MoonMEToLatLon(p, out lat, out lon);
            lat *= 180d / Math.PI;
            lon *= 180d / Math.PI;
        }

        public Vector3d LineSampleHeightToMoonME(int line, int sample, float height)
        {
            var radius = MoonRadius + height / 1000d;
            GetLatLon(line, sample, out double lat, out double lon);
            var z = radius * Math.Sin(lat);
            var c = radius * Math.Cos(lat);
            var x = c * Math.Cos(lon);  // TODO: Not sure about these
            var y = c * Math.Sin(lon);
            return new Vector3d(x, y, z);
        }

        public static Point GetLineSampleDegrees(double lat_deg, double lon_deg)
        {
            GetLineSampleDegrees(lat_deg, lon_deg, out int line, out int sample);
            return new Point(sample, line);
        }

        public static void GetLineSampleDegrees(double lat_deg, double lon_deg, out int line, out int sample) =>
            GetLineSample(lat_deg * Math.PI / 180d, lon_deg * Math.PI / 180d, out line, out sample);

        // Checked against GetLatLon
        public static void GetLineSample(double lat_rad, double lon_rad, out int line, out int sample)
        {
            double x, y;
            if (InMemoryInt16Terrain.Singleton.IsNorth)
            {
                var temp1 = MoonRadius * Math.Tan(Math.PI / 4d - lat_rad / 2d);
                var temp2 = lon_rad - LonP;
                x = 2d * temp1 * Math.Sin(temp2);
                y = -2d * temp1 * Math.Cos(temp2);
            }
            else
            {
                var temp1 = MoonRadius * Math.Tan(Math.PI / 4d + lat_rad / 2d);
                var temp2 = lon_rad - LonP;
                x = 2d * temp1 * Math.Sin(temp2);
                y = 2d * temp1 * Math.Cos(temp2);
            }

            //x = (Sample - S0 - 1) * Scale
            //y = (1 - L0 - Line) * Scale => Line = 1 - L0 - y/Scale
            sample = (int)(x / Scale + S0 + 1d + 0.5d);

            // line = (int)(1d - L0 - y / Scale);   // Should be this per the notes
            line = -(int)(1d - L0 - y / Scale);  //NOTE: I'm not sure why there's a - sign here

            // I don't understand this correction
            const int mid = 30400 / 2;
            line = mid - (line - mid);
        }

        public static void GetLineSampleD(double lat_rad, double lon_rad, out double line, out double sample)
        {
            double x, y;
            if (InMemoryInt16Terrain.Singleton.IsNorth)
            {
                var temp1 = MoonRadius * Math.Tan(Math.PI / 4d - lat_rad / 2d);
                var temp2 = lon_rad - LonP;
                x = 2d * temp1 * Math.Sin(temp2);
                y = -2d * temp1 * Math.Cos(temp2);
            }
            else
            {
                var temp1 = MoonRadius * Math.Tan(Math.PI / 4d + lat_rad / 2d);
                var temp2 = lon_rad - LonP;
                x = 2d * temp1 * Math.Sin(temp2);
                y = 2d * temp1 * Math.Cos(temp2);
            }

            //x = (Sample - S0 - 1) * Scale
            //y = (1 - L0 - Line) * Scale => Line = 1 - L0 - y/Scale
            sample = x / Scale + S0;

            // line = (int)(1d - L0 - y / Scale);   // Should be this per the notes
            //line = -(1d - L0 - y / Scale);  //NOTE: I'm not sure why there's a - sign here
            line = -y / Scale + L0;

            // I don't understand this correction
            //const int mid = 30400 / 2;
            //line = mid - (line - mid);
        }

        public float InterpolatedElevation(float line, float sample)
        {
            var x1 = (int)sample;
            var y1 = (int)line;
            int x2 = x1 + 1;
            int y2 = y1 + 1;
            float q11 = LineSampleToTerrainOffset(y1, x1);
            float q12 = LineSampleToTerrainOffset(y2, x1);
            float q21 = LineSampleToTerrainOffset(y1, x2);
            float q22 = LineSampleToTerrainOffset(y2, x2);
            // From https://en.wikipedia.org/wiki/Bilinear_interpolation, note denominator is 1 in this case, so is omitted
            float r = (q11 * (x2 - sample) * (y2 - line) + q21 * (sample - x1) * (y2 - line) + q12 * (x2 - sample) * (line - y1) + q22 * (sample - x1) * (line - y1));
            //Console.WriteLine(@"  x={0:G7} y={1:G7} x1={2} y1={3} x2={4} y2={5} q11={6:G7} q12={7:G7} q21={8:G7} q22={9:G7}", x, y, x1, y1, x2, y2, q11, q12, q21, q22);
            return r;
        }

        public Vector3d GetPointInME(int line, int sample)
        {
            var relz = LineSampleToTerrainOffset(line, sample);
            var radius = MoonRadius + relz / 1000d;
            GetLatLon(line, sample, out double lat, out double lon);
            var z = radius * Math.Sin(lat);
            var c = radius * Math.Cos(lat);
            var x = c * Math.Cos(lon);  // TODO: Not sure about these
            var y = c * Math.Sin(lon);
            return new Vector3d(x, y, z);
        }

        public Vector3d GetPointInME(float line, float sample)
        {
            var relz = InterpolatedElevation(line, sample);
            var radius = MoonRadius + relz / 1000d;
            GetLatLon(line, sample, out double lat, out double lon);
            var z = radius * Math.Sin(lat);
            var c = radius * Math.Cos(lat);
            var x = c * Math.Cos(lon);  // TODO: Not sure about these
            var y = c * Math.Sin(lon);
            return new Vector3d(x, y, z);
        }

        public Vector3d GetPointInME(PointF p) => GetPointInME(p.Y, p.X);

        public float LineSampleToTerrainOffset(int line, int sample)
        {
            const int height = 30400;
            const int width = 30400;
            const int stride = width;
            System.Diagnostics.Debug.Assert(line >= 0 && line < height && sample >= 0 && sample < width);

            int index = line * stride + sample;
            return Data[index] / 1000f;  // product is in kilometers, return value is in meters
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public float Fetch(int line, int sample) => Data[line * 30400 + sample];

        public void GetMinMax()
        {
            var min = float.MaxValue;
            var max = float.MinValue;
            var length = Data.Length;
            for (var i = 0; i < length; i++)
            {
                var v = Data[i];
                if (v < min) min = v;
                if (v > max) max = v;
            }
            _max = max;
            _min = min;
        }

        #region IDisposable Support
        private bool disposedValue = false; // To detect redundant calls

        protected virtual void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (disposing)
                {
                    // TODO: dispose managed state (managed objects).
                }

                // TODO: free unmanaged resources (unmanaged objects) and override a finalizer below.
                Data = null;

                disposedValue = true;
            }
        }

        // TODO: override a finalizer only if Dispose(bool disposing) above has code to free unmanaged resources.
        // ~InMemoryTerrainManager() {
        //   // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
        //   Dispose(false);
        // }

        // This code added to correctly implement the disposable pattern.
        public override void Dispose()
        {
            // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
            Dispose(true);
            // TODO: uncomment the following line if the finalizer is overridden above.
            // GC.SuppressFinalize(this);
        }
        #endregion
    }

    /// <summary>
    ///  This is a memory mapped version that came from the traverse planner.  This will be folded into the new implementation.
    /// </summary>
    public class MemoryMappedInt16Terrain
    {
        static string PolarStereographic = @"PROJCS[""Moon_North_Pole_Stereographic"",
    GEOGCS[""Moon 2000"",
        DATUM[""D_Moon_2000"",
            SPHEROID[""Moon_2000_IAU_IAG"",1737400.0,0.0]],
        PRIMEM[""Greenwich"",0],
        UNIT[""Decimal_Degree"",0.0174532925199433]],
    PROJECTION[""Stereographic""],
    PARAMETER[""False_Easting"",0],
    PARAMETER[""False_Northing"",0],
    PARAMETER[""Central_Meridian"",0],
    PARAMETER[""Scale_Factor"",1],
    PARAMETER[""Latitude_Of_Origin"",90],
    UNIT[""Meter"",1]]";

        //private const string TerrainImageFile = @"C:\UVS\svn\src\TestData\Terrain\ldem_64.img";
        const string TerrainImageFile = @"C:\RP\maps\DEM\ldem_80n_20m.img";
        MemoryMappedViewAccessor _accessor;
        MemoryMappedFile _mmf;

        double[] _affineTransform;
        OSGeo.OSR.CoordinateTransformation _pixelToLineSample;

        public static MemoryMappedInt16Terrain CachedTerrainManager;

        public MemoryMappedInt16Terrain(string referenceImage)
        {
            Open(TerrainImageFile);

            using (OSGeo.GDAL.Dataset dataset = OSGeo.GDAL.Gdal.Open(referenceImage, OSGeo.GDAL.Access.GA_ReadOnly))
            {
                _affineTransform = new double[6];
                dataset.GetGeoTransform(_affineTransform);

                string projection = dataset.GetProjectionRef();

                var src = new OSGeo.OSR.SpatialReference(projection);
                var dst = new OSGeo.OSR.SpatialReference(PolarStereographic);
                _pixelToLineSample = new OSGeo.OSR.CoordinateTransformation(src, dst);
            }
        }

        public void Open(string filename)
        {
            _mmf = MemoryMappedFile.CreateFromFile(filename, FileMode.Open, "MissionManager", 0L, MemoryMappedFileAccess.Read);

            // This maps the whole file.  There is an accessor creation method
            // that maps only a section of the file, which could be used for larger
            // files, but it would need to know what the lat/lon are first.
            //_accessor = _mmf.CreateViewAccessor();
            _accessor = _mmf.CreateViewAccessor(0, 0L, MemoryMappedFileAccess.Read);
        }

        public void Close()
        {
            Dispose(true);
        }

        public void Dispose(bool ignore)
        {
            if (_mmf == null) return;
            _mmf.Dispose();
            _mmf = null;
            _accessor = null;
        }

        /// <summary>
        /// Read Terrain Height
        /// </summary>
        /// <param name="line"></param>
        /// <param name="sample"></param>
        /// <returns>terrain height in m relative to the reference radius (1737.4 km)</returns>
        public float LineSampleToTerrainOffset(int line, int sample)
        {
            const int height = 30400;
            const int width = 30400;
            const int stride = width * 2;

            if (line < 0 || line >= height || sample < 0 || sample >= width)
                throw new Exception(@"line,sample out-of-bounds");

            int index = line * stride + sample * 2;
            short h = _accessor.ReadInt16(index);
            return 0.5f * h;
        }

        public short LineSampleToRawTerrainOffset(int line, int sample)
        {
            const int height = 30400;
            const int width = 30400;
            const int stride = width * 2;

            if (line < 0 || line >= height || sample < 0 || sample >= width)
                throw new Exception(@"line,sample out-of-bounds");

            int index = line * stride + sample * 2;
            short h = _accessor.ReadInt16(index);
            return h;
        }

        public float PixelToTerrainOffset(int pixelX, int pixelY)
        {
            double x = pixelX * _affineTransform[1] + _affineTransform[0];
            double y = pixelY * _affineTransform[5] + _affineTransform[3];

            var o = new double[3];
            _pixelToLineSample.TransformPoint(o, x, y, 0d);

            double x1 = o[0];
            double y1 = o[1];

            int line = 30400 / 2 - (int)(y1 / 20d + 0.5d);
            int sample = 30400 / 2 + (int)(x1 / 20d + 0.5d);

            return LineSampleToTerrainOffset(line, sample);
        }

        public short PixelToRawTerrainOffset(int pixelX, int pixelY)
        {
            double x = pixelX * _affineTransform[1] + _affineTransform[0];
            double y = pixelY * _affineTransform[5] + _affineTransform[3];

            var o = new double[3];
            _pixelToLineSample.TransformPoint(o, x, y, 0d);

            double x1 = o[0];
            double y1 = o[1];

            int line = 30400 / 2 - (int)(y1 / 20d + 0.5d);
            int sample = 30400 / 2 + (int)(x1 / 20d + 0.5d);

            return LineSampleToRawTerrainOffset(line, sample);
        }
    }
}
