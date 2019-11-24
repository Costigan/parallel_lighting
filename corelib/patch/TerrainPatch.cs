using viper.corelib.math;
using viper.corelib.terrain;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using viper.corelib.utilities;

namespace viper.corelib.patch
{
    /// <summary>
    /// Terrain Patches know about the DEM and coordinate systems
    /// </summary>
    public class TerrainPatch : IDisposable
    {
        public const double DEM_Step = 20d;
        public const int DEM_size = 30400;
        public const int DefaultSize = 128;
        public const double MoonRadius = 1737.4d;
        public const int NearHorizonOversample = 2;

        public const double RayAngleStart = 0d;
        public const double RayAngleStop = Math.PI * 2d; // 10d * Math.PI / 180d;  // Math.PI * 2d;
        public const double RayAngleStep = Math.PI * 2d / (NearHorizonOversample * Horizon.HorizonSamples);

        /// <summary>
        /// These constants control the extent of the local horizon.  This distance (230 map units) is
        /// the distance at which the map pixels are less than 0.25 deg apart, which is the spacing of the current horizon mask
        /// </summary>
        public static double MaximumLocalDistance = 230;  // 1d / Math.Tan(025d * Math.PI/180d)
        public const float LocalStep = .70710698f;

        public int Line { get; set; }
        public int Sample { get; set; }
        public int Width { get; set; } = DefaultSize;
        public int Height { get; set; } = DefaultSize;
        public int Step { get; set; } = 1;
        public int ObserverHeightInMeters { get; set; }  // Meters.  Int to restrict options

        public Rectangle Bounds => new Rectangle(Sample, Line, Width, Height);

        protected string _Path = null;
        public string Path
        {
            get
            {
                return _Path ?? DefaultPath;
            }
            set
            {
                _Path = value;
            }
        }

        public bool UseCenterCache = true;

        // These are the ids of the patches that horizon filtering says can cast shadows on the current patch.  OLD CODE
        public List<Point> ShadowCasters = new List<Point>();

        public Point Id => new Point(Sample / DefaultSize, Line / DefaultSize);

        /// <summary>
        /// Surface points in MOON_ME frame
        /// Access like so Points[h][w][0] for X value;
        /// h corresponds to a y or a line value
        /// w corresponds to an x or a sample value
        /// </summary>
        public Vector3d[][] Points;
        public Matrix4d[][] Matrices;
        public Horizon[][] Horizons;
        public Vector3d[] TestVectors;
        protected Point TallestPoint;

        protected bool _centerMatrixIsValid;
        protected Matrix4d _centerMatrix;
        public Bitmap SunShadowBitmap;
        public Bitmap EarthShadowBitmap;
        protected float _lastShadowAzimuth = float.MinValue;
        protected float _lastShadowElevation = float.MinValue;

        public override string ToString() => $"<TerrainPatch l={Line} s={Sample} step={Step} h={ObserverHeightInMeters}";

        public override bool Equals(object obj)
        {
            if (!(obj is TerrainPatch o))
                return false;
            return Id.Equals(o.Id) && Step.Equals(o.Step) && Width == o.Width && Height == o.Height;
        }

        public override int GetHashCode()
        {
            unchecked
            {
                int hash = 17;
                hash = hash * 23 + Sample.GetHashCode();
                hash = hash * 23 + Line.GetHashCode();
                hash = hash * 23 + Step.GetHashCode();
                hash = hash * 23 + Width.GetHashCode();
                hash = hash * 23 + Height.GetHashCode();
                return hash;
            }
        }

        public static TerrainPatch FromId(int x, int y, int observer = 0)
            => new TerrainPatch { Line = y * DefaultSize, Sample = x * DefaultSize, Width = DefaultSize, Height = DefaultSize, ObserverHeightInMeters = observer };
        public static TerrainPatch FromId(Point id, int observer = 0) => FromId(id.X, id.Y, observer);
        public static TerrainPatch FromLineSample(int line, int sample, int observer = 0) => new TerrainPatch { Line = line, Sample = sample, Width = DefaultSize, Height = DefaultSize, ObserverHeightInMeters = observer };
        public static Point LineSampleToId(int line, int sample) => new Point(sample / DefaultSize, line / DefaultSize);

        public static Point LineSampleToId(Point p) => new Point(p.X / DefaultSize, p.Y / DefaultSize);
        public static Point IdToLineSample(Point p) => new Point(p.X * DefaultSize, p.Y * DefaultSize);

        public static Rectangle CoveringIdRectangle(Rectangle pixels)
        {
            var ul = LineSampleToId(pixels.Location);
            var lr = LineSampleToId(pixels.Bottom, pixels.Right);
            //return new Rectangle(ul.X, ul.Y, lr.X - ul.X + 1, lr.Y - ul.Y + 1);
            return new Rectangle(ul.X, ul.Y, lr.X - ul.X, lr.Y - ul.Y);
        }

        public static IEnumerable<Point> EnumerateIds(Rectangle ids)
        {
            for (var line = ids.Top; line < ids.Bottom; line++)
                for (var sample = ids.Left; sample < ids.Right; sample++)
                    yield return new Point(sample, line);
        }

        public Point PointInPatch(Point p) => new Point(p.X - Sample, p.Y - Line);
        public Point PointInPatchPointInDEM(Point p) => new Point(p.X + Sample, p.Y + Line);

        #region Accessors

        public Horizon GetHorizon(int line, int sample) => Horizons[line][sample];

        #endregion

        /// <summary>
        /// Return tile whose position is x,y relative to this tile in id coordinates
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public TerrainPatch RelativeTo(int x, int y)
        {
            var newline = Math.Max(0, Line + y * Height);
            var newheight = Math.Min(y < 0 ? Line : DEM_size, newline + Height) - newline;
            var newsample = Math.Max(0, Sample + x * Width);
            var newwidth = Math.Min(x < 0 ? Sample : DEM_size, newsample + Width) - newsample;
            // NOTE: This is where a cache fetch would occur
            return new TerrainPatch { Line = newline, Sample = newsample, Width = newwidth, Height = newheight };
        }

        public TerrainPatch RelativeTo(Point p) => RelativeTo(p.X, p.Y);

        /// <summary>
        /// Return a stream of Ids that are exactly range patches away from the current patch in one axis.
        /// </summary>
        /// <param name="range"></param>
        /// <returns></returns>
        public IEnumerable<Point> SurroundingIds(int range)
        {
            const int MaxId = DEM_size / 128;
            var id = Id;
            // Fix the X distance
            if (id.X - range >= 0)
            {
                var x1 = id.X - range;
                var y1 = Math.Max(id.Y - range, 0);
                var y2 = Math.Min(id.Y + range, MaxId - 1);
                for (var y = y1; y <= y2; y++)
                    yield return new Point(x1, y);
            }
            if (id.X + range < MaxId)
            {
                var x1 = id.X + range;
                var y1 = Math.Max(id.Y - range, 0);
                var y2 = Math.Min(id.Y + range, MaxId - 1);
                for (var y = y1; y <= y2; y++)
                    yield return new Point(x1, y);
            }
            // Now, fix the Y.  Make sure not to repeat the corners.
            if (id.Y - range >= 0)
            {
                var y1 = id.Y - range;
                var x1 = Math.Max(id.X - range + 1, 0);
                var x2 = Math.Min(id.X + range - 1, MaxId - 1);
                for (var x = x1; x <= x2; x++)
                    yield return new Point(x, y1);
            }
            if (id.Y + range < MaxId)
            {
                var y1 = id.Y + range;
                var x1 = Math.Max(id.X - range + 1, 0);
                var x2 = Math.Min(id.X + range - 1, MaxId - 1);
                for (var x = x1; x <= x2; x++)
                    yield return new Point(x, y1);
            }
        }

        public IEnumerable<TerrainPatch> SurroundingPatches(int range) => SurroundingIds(range).Select(p => TerrainPatch.FromId(p));

        public void FillPoints(InMemoryInt16Terrain m)
        {
            if (Points != null)
                return;
            Points = new Vector3d[Height][];
            var maxRadius = double.MinValue;
            for (var h = 0; h < Height; h++)
            {
                int line = Line + h * Step;
                var row = new Vector3d[Width];
                Points[h] = row;
                for (var w = 0; w < Width; w++)
                {
                    var sample = Sample + w * Step;
                    var relz = line < DEM_size && sample < DEM_size ? m.LineSampleToTerrainOffset(line, sample) : -12000f;  // guaranteed to be less than any real terrain
                    var radius = MoonRadius + relz / 1000d;  // km
                    if (radius > maxRadius)
                    {
                        maxRadius = radius;
                        TallestPoint = new Point(w, h);
                    }
                    InMemoryInt16Terrain.GetLatLon(line, sample, out double lat, out double lon);
                    var z = radius * Math.Sin(lat);
                    var c = radius * Math.Cos(lat);
                    var x = c * Math.Cos(lon);  // TODO: Not sure about these
                    var y = c * Math.Sin(lon);
                    row[w] = new Vector3d(x, y, z);  // km
                }
            }

            // Artificially increase maxRadius by 200 meters
            maxRadius += 0.2d;

            // Set up the test vectors
            TestVectors = PoleLocations().Select(loc =>
            {
                int line = Line + loc.Item1;
                var sample = Sample + loc.Item2;
                var radius = maxRadius;
                InMemoryInt16Terrain.GetLatLon(line, sample, out double lat, out double lon);
                var z = radius * Math.Sin(lat);
                var c = radius * Math.Cos(lat);
                var x = c * Math.Cos(lon);
                var y = c * Math.Sin(lon);
                return new Vector3d(x, y, z);
            }).ToArray();

            //Console.WriteLine($"FillPoints for [{Line},{Sample}] {GetHashCode()}");
        }

        protected IEnumerable<Tuple<int, int>> PoleLocations()
        {
            yield return new Tuple<int, int>(0, 0);
            yield return new Tuple<int, int>(0, Width - 1);
            yield return new Tuple<int, int>(Height - 1, 0);
            yield return new Tuple<int, int>(Height - 1, Width - 1);
            yield return new Tuple<int, int>(Height / 2, Width / 2);
        }

        public void FillPointsAndMatrices(InMemoryInt16Terrain m)
        {
            if (Matrices != null)
                return;
            Points = new Vector3d[Height][];
            Matrices = new Matrix4d[Height][];
            var maxRadius = double.MinValue;
            for (var h = 0; h < Height; h++)
            {
                int line = Line + h;
                var row = new Vector3d[Width];
                Points[h] = row;
                var mrow = new Matrix4d[Width];
                Matrices[h] = mrow;
                for (var w = 0; w < Width; w++)
                {
                    var sample = Sample + w;
                    var relz = m.LineSampleToTerrainOffset(line, sample);
                    var radius = MoonRadius + relz / 1000d;

                    if (radius > maxRadius)
                    {
                        maxRadius = radius;
                        TallestPoint = new Point(w, h);
                    }

                    InMemoryInt16Terrain.GetLatLon(line, sample, out double lat, out double lon);
                    var z = radius * Math.Sin(lat);
                    var c = radius * Math.Cos(lat);
                    var x = c * Math.Cos(lon);  // TODO: Not sure about these
                    var y = c * Math.Sin(lon);
                    var vec = new Vector3d(x, y, z);
                    row[w] = vec;

                    mrow[w] = MatrixFromLatLon(lat, lon, vec);
                }
            }
            if (Horizons == null)
            {
                Horizons = new Horizon[Height][];
                for (var h = 0; h < Height; h++)
                {
                    var hrow = new Horizon[Width];
                    Horizons[h] = hrow;
                    for (var w = 0; w < Width; w++)
                        hrow[w] = new Horizon();
                }
            }
        }

        public void FillMatrices(InMemoryInt16Terrain m)
        {
            if (Matrices != null)
                return;
            Matrices = new Matrix4d[Height][];
            for (var h = 0; h < Height; h++)
            {
                int line = Line + h;
                var mrow = new Matrix4d[Width];
                Matrices[h] = mrow;
                for (var w = 0; w < Width; w++)
                {
                    var sample = Sample + w;
                    var relz = m.LineSampleToTerrainOffset(line, sample);
                    var radius = MoonRadius + relz / 1000d;
                    InMemoryInt16Terrain.GetLatLon(line, sample, out double lat, out double lon);
                    var z = radius * Math.Sin(lat);
                    var c = radius * Math.Cos(lat);
                    var x = c * Math.Cos(lon);  // TODO: Not sure about these
                    var y = c * Math.Sin(lon);
                    var vec = new Vector3d(x, y, z);

                    mrow[w] = MatrixFromLatLon(lat, lon, vec);
                }
            }
        }

        public void FillMatricesRelativeToPoint(InMemoryInt16Terrain m, Vector3d base_point)
        {
            if (Matrices != null)
                return;
            Matrices = new Matrix4d[Height][];
            for (var h = 0; h < Height; h++)
            {
                int line = Line + h;
                var mrow = new Matrix4d[Width];
                Matrices[h] = mrow;
                for (var w = 0; w < Width; w++)
                {
                    var sample = Sample + w;
                    var relz = m.LineSampleToTerrainOffset(line, sample);
                    var radius = MoonRadius + relz / 1000d;
                    InMemoryInt16Terrain.GetLatLon(line, sample, out double lat, out double lon);
                    var z = radius * Math.Sin(lat);
                    var c = radius * Math.Cos(lat);
                    var x = c * Math.Cos(lon);  // TODO: Not sure about these
                    var y = c * Math.Sin(lon);
                    var vec = new Vector3d(x, y, z);

                    var mat = MatrixFromLatLon(lat, lon, vec - base_point);

                    //var vz = new Vector3d(0d, 0d, 1d);
                    //var vz1 = Vector3d.Transform(vz, mat);
                    //var vz1_length = vz1.Length;

                    mrow[w] = mat;
                }
            }
        }

        public Matrix4d MatrixFromLatLon(double lat, double lon, Vector3d vec)
        {
            var zaxis = new Vector3d(0d, 0d, 1d);
            var yaxis = new Vector3d(0d, 1d, 0d);
            var a = Matrix4d.CreateFromAxisAngle(zaxis, -lon);
            var b = Matrix4d.CreateFromAxisAngle(yaxis, -(Math.PI / 2 - lat));
            var c = Matrix4d.CreateTranslation(-vec);
            return c * a * b;  // working
        }

        public void MergeHorizonFrom(TerrainPatch other)
        {
            if (other.Horizons == null)
                return;
            if (Horizons == null)
                InitializeHorizons();
            if (Horizons == null || Horizons.Length < 1 || Horizons.Length != other.Horizons.Length || Horizons[0].Length != other.Horizons[0].Length)
                throw new Exception("Asked to merge incompatible horizon arrays");
            for (var line = 0; line < Horizons.Length; line++)
            {
                var rowa = Horizons[line];
                var rowb = other.Horizons[line];
                for (var sample = 0; sample < rowa.Length; sample++)
                {
                    var arya = rowa[sample].Buffer;
                    var aryb = rowb[sample].Buffer;
                    for (var i = 0; i < arya.Length; i++)
                        if (aryb[sample] > arya[sample])
                            arya[sample] = aryb[sample];
                }
            }
        }

        protected void MergeHorizonInto(short[] primary, short[] secondary)
        {
            Debug.Assert(primary.Length == secondary.Length);
            lock (primary)
            {
                for (var i = 0; i < primary.Length; i++)
                {
                    if (secondary[i] > primary[i])
                        primary[i] = secondary[i];
                }
            }
        }

        public void UpdateHorizon(TerrainPatch o, Horizon horizon_buffer = null)
        {
            var deltaSample = Math.Abs(Sample - o.Sample);
            var deltaLine = Math.Abs(Line - o.Line);

            if (horizon_buffer == null)
                horizon_buffer = new Horizon();
            var temp = new Vector3d();
            var oPoints = o.Points;
            var oHeight = o.Height;
            var oWidth = o.Width;
            for (var line = 0; line < Height; line++)
            {
                var row = Matrices[line];
                for (var sample = 0; sample < Width; sample++)
                {
                    //DEBUG
                    //if (line != 0 || sample != 0)
                    //    continue;

                    var m = row[sample];
                    var pixel_horizon = Horizons[line][sample];

                    // Initialize local, temporary horizon buffer
                    horizon_buffer.Init();

                    // Iterate through the other patch
                    for (var oline = 0; oline < oHeight; oline++)
                    {
                        var orow = oPoints[oline];
                        for (var osample = 0; osample < oWidth; osample++)
                        {
                            var opoint = orow[osample];
                            Transform(ref opoint, ref m, ref temp);

                            // temp contains the remote point in the horizon frame
                            var x = temp[0];
                            var y = temp[1];
                            var z = temp[2];
                            var azimuth = Math.Atan2(y, x) + Math.PI;  // [0,2PI]
                            var alen = Math.Sqrt(x * x + y * y);
                            var slope = z / alen;
                            horizon_buffer.Merge(azimuth, slope);
                        }
                    }

                    // Merge into global horizon buffer for the pixel.  Trying to hold this lock minimally
                    pixel_horizon.Merge(horizon_buffer);
                }
            }
            AddShadowCaster(o.Id);
        }

        private const int InhibitHorizonClipDistance = 6;

        public bool IsOverHorizon(TerrainPatch o)
        {
            // This is a hack to ignore this check for patches that are closer than some limit
            {
                var dx = Id.X - o.Id.X;
                var dy = Id.Y - o.Id.Y;
                if (dx * dx + dy * dy <= InhibitHorizonClipDistance * InhibitHorizonClipDistance)
                    return true;// true means that we want to check this shadow caster
            }
            var temp = new Vector3d();
            var oPoints = o.Points;
            var oHeight = o.Height;
            var oWidth = o.Width;
            var m = Matrices[TallestPoint.Y][TallestPoint.X];
            var pixel_horizon = Horizons[TallestPoint.Y][TallestPoint.X];

            var isDegrees = pixel_horizon.IsDegrees;
            //Debug.Assert(!isDegrees);

            var testVectors = o.TestVectors;
            for (var i = 0; i < testVectors.Length; i++)
            {
                var opoint = testVectors[i];
                Transform(ref opoint, ref m, ref temp);

                // temp contains the remote point in the horizon frame
                var x = temp[0];
                var y = temp[1];
                var z = temp[2];
                var azimuth_rad = Math.Atan2(y, x) + Math.PI;  // [0,2PI]
                var alen = Math.Sqrt(x * x + y * y);

                var elevation_deg_or_slope = isDegrees ? Math.Atan2(z, alen) * 180d / Math.PI : z / alen;

                if (pixel_horizon.IsOverHorizon(azimuth_rad, elevation_deg_or_slope))
                    return true;
            }
            return false;
        }

        //NOTE: I'm not sure what the basePoint is here.  The other version doesn't need it.
        public bool IsOverHorizon(TerrainPatch o, Vector3d basePoint)
        {
            if (Matrices == null)
                Console.WriteLine(@"Shouldn't get here");
            // This is a hack to ignore this check for patches that are closer than some limit
            {
                var dx = Id.X - o.Id.X;
                var dy = Id.Y - o.Id.Y;
                if (dx * dx + dy * dy <= InhibitHorizonClipDistance * InhibitHorizonClipDistance)
                    return true;// true means that we want to check this shadow caster
            }
            var temp = new Vector3d();
            var oPoints = o.Points;
            var oHeight = o.Height;
            var oWidth = o.Width;
            var testVectors = o.TestVectors;

            foreach (var terrain_pt in HorizonTestPoints())
            {
                var m = Matrices[terrain_pt.Y][terrain_pt.X];
                var pixel_horizon = Horizons[terrain_pt.Y][terrain_pt.X];
                Debug.Assert(!pixel_horizon.IsDegrees);

                for (var i = 0; i < testVectors.Length; i++)
                {
                    var opoint = testVectors[i] - basePoint;
                    Transform(ref opoint, ref m, ref temp);

                    // temp contains the remote point in the horizon frame
                    var x = temp[0];
                    var y = temp[1];
                    var z = temp[2];
                    var azimuth_rad = Math.Atan2(y, x) + Math.PI;  // [0,2PI]
                    var alen = Math.Sqrt(x * x + y * y);
                    var slope = z / alen;

                    if (pixel_horizon.IsOverHorizon(azimuth_rad, slope))
                        return true;
                }
            }
            return false;
        }

        public IEnumerable<Point> HorizonTestPoints()
        {
            yield return TallestPoint;
            yield return new Point(0, 0);
            yield return new Point(Height - 1, 0);
            yield return new Point(Height - 1, Width - 1);
            yield return new Point(0, Width - 1);
        }

        public static void Transform(ref Vector3d vec, ref Matrix4d mat, ref Vector3d output)
        {
            output.X = vec.X * mat.Row0.X + vec.Y * mat.Row1.X + vec.Z * mat.Row2.X + mat.Row3.X;
            output.Y = vec.X * mat.Row0.Y + vec.Y * mat.Row1.Y + vec.Z * mat.Row2.Y + mat.Row3.Y;
            output.Z = vec.X * mat.Row0.Z + vec.Y * mat.Row1.Z + vec.Z * mat.Row2.Z + mat.Row3.Z;
        }

        public static void Transform(ref Vector4d vec, ref Matrix4d mat, ref Vector4d output)
        {
            output.X = vec.X * mat.Row0.X + vec.Y * mat.Row1.X + vec.Z * mat.Row2.X + mat.Row3.X;
            output.Y = vec.X * mat.Row0.Y + vec.Y * mat.Row1.Y + vec.Z * mat.Row2.Y + mat.Row3.Y;
            output.Z = vec.X * mat.Row0.Z + vec.Y * mat.Row1.Z + vec.Z * mat.Row2.Z + mat.Row3.Z;
            output.W = vec.X * mat.Row0.Z + vec.Y * mat.Row1.Z + vec.Z * mat.Row2.Z + mat.Row3.Z;
        }

        public void ParallelAddNearHorizon(InMemoryInt16Terrain m)
        {
            if (Horizons == null)
                InitializeHorizons();
            if (UseCenterCache && File.Exists(DefaultCenterPath))
            {
                var p = ReadFrom(DefaultCenterPath);
                MergeHorizonFrom(p);
                return;
            }
            FillPointsAndMatrices(m);
            Parallel.For(0, Height, new ParallelOptions { MaxDegreeOfParallelism = ViperEnvironment.MaxDegreeOfParallelism }, l =>
            {
                for (var s = 0; s < Width; s++)
                    AddNearHorizon(m, l, s, Horizons);
            });
            if (UseCenterCache)
                Write(DefaultCenterPath);
        }

        public void AddNearHorizon(InMemoryInt16Terrain tm, int l, int s, Horizon[][] horizons)
        {
            //Console.WriteLine($"[{l},{s}]");
            var temp = new Vector3d();
            var mat = Matrices[l][s];
            var horizon = horizons[l][s];
            Debug.Assert(!horizon.IsDegrees);
            var ld = (float)(Line + l);
            var sd = (float)(Sample + s);
            const float dem_bound = DEM_size - 2f;
            if (RayAngleStart > 0d || RayAngleStop < Math.PI * 2d)
            {
                const double elev = 80d * Math.PI / 180d;
                for (var ang = 0d; ang < RayAngleStart; ang += RayAngleStep)
                    horizon.Merge(ang, elev);
                for (var ang = RayAngleStop; ang < Math.PI * 2d; ang += RayAngleStep)
                    horizon.Merge(ang, elev);
            }
            for (var ray = RayAngleStart; ray < RayAngleStop; ray += RayAngleStep)
            {
                var raycos = (float)Math.Cos(ray);
                var raysin = (float)Math.Sin(ray);
                for (var distance = 1f; distance < MaximumLocalDistance; distance += LocalStep)
                {
                    var l1 = ld + raysin * distance;
                    var s1 = sd + raycos * distance;
                    if (l1 < 1f || s1 < 1f || l1 > dem_bound || s1 > dem_bound)
                        break;  // If we step outside, we'll be outside until the end of this ray, so quit now
                    var e = tm.InterpolatedElevation(l1, s1);
                    var radius = MoonRadius + e / 1000d;
                    InMemoryInt16Terrain.GetLatLon(l1, s1, out double lat, out double lon);
                    var z = radius * Math.Sin(lat);
                    var c = radius * Math.Cos(lat);
                    var x = c * Math.Cos(lon);  // TODO: Not sure about these
                    var y = c * Math.Sin(lon);
                    var vec = new Vector3d(x, y, z);
                    Transform(ref vec, ref mat, ref temp);

                    // temp contains the remote point in the horizon frame
                    var rx = temp[0];
                    var ry = temp[1];
                    var rz = temp[2];
                    var azimuth = Math.Atan2(ry, rx) + Math.PI;  // [0,2PI]
                    var alen = Math.Sqrt(rx * rx + ry * ry);
                    var slope = Math.Atan2(rz, alen);

                    horizon.Merge(azimuth, slope);
                }
            }
        }

        public void InitializeHorizons()
        {
            // This does not change the existing arrays but replaces them
            Horizons = new Horizon[Height][];
            for (var h = 0; h < Height; h++)
            {
                Horizons[h] = new Horizon[Width];
                for (var w = 0; w < Width; w++)
                    Horizons[h][w] = new Horizon();
            }
        }

        public void LoadHorizons(float[] a)
        {
            if (a.Length != Horizon.HorizonSamples * Width * Height)
                return;
            Horizons = new Horizon[Height][];
            for (var line = 0; line < Height; line++)
            {
                var row = new Horizon[Width];
                Horizons[line] = row;
                for (var sample = 0; sample < Width; sample++)
                    row[sample] = Horizon.From(a, (line * Height + sample) * Horizon.HorizonSamples);
            }
        }

        public void ConvertSlopeToDegrees()
        {
            if (Horizons == null)
                return;
            for (var line = 0; line < TerrainPatch.DefaultSize; line++)
                for (var sample = 0; sample < TerrainPatch.DefaultSize; sample++)
                    Horizons[line][sample].ConvertSlopeToDegrees();
        }

        public void AddShadowCaster(Point p)
        {
            lock (ShadowCasters)
            {
                if (!ShadowCasters.Contains(p))
                    ShadowCasters.Add(p);
            }
        }

        #region Persist to Disk

        // Ignoring this
        public static string GeneratorPostfix = "gpu";  // must be 3 chars

        public string DefaultCenterFilename => string.Format(@"{1:D5}_{2:D5}_{0}_center.bin", GeneratorPostfix, Line, Sample);
        public string DefaultCenterPath => System.IO.Path.Combine(ViperEnvironment.HorizonRoot, DefaultCenterFilename);

        public string DefaultFilename => string.Format(@"horizon_{0:D5}_{1:D5}_{2:D3}.bin", Line, Sample, ObserverHeightInMeters);
        public string DefaultPath => System.IO.Path.Combine(ViperEnvironment.HorizonRoot, DefaultFilename);

        public static bool TryParseFilename(string filename, out int line, out int sample, out int observer_height)
        {
            line = sample = observer_height = 0;
            return filename.Length == 27 &&
                ".bin".Equals(filename.Substring(23, 4)) &&
                int.TryParse(filename.Substring(8, 5), out line) &&
                int.TryParse(filename.Substring(14, 5), out sample) &&
                int.TryParse(filename.Substring(20, 3), out observer_height);
        }

        public void Read(string path)
        {
            path = path ?? Path;
            var filename = System.IO.Path.GetFileName(path);
            if (!TryParseFilename(filename, out int line, out int sample, out int observer_height))
                throw new Exception(@"Couldn't parse filename when reading horizons");

            using (var fs = File.OpenRead(path))
                Read(fs, line, sample, observer_height, filename);
        }

        public void Read(Stream stream, int line = -1, int sample = -1, int observer_height = -1, string filename = null)
        {
            //using (var gz = new GZipStream(gzippedStream, CompressionMode.Decompress))
            using (var br = new BinaryReader(stream))
            {
                Line = br.ReadInt32();
                if (line > -1 && Line != line) throw new Exception($"Line in this horizon file doesn't match the line in the filename: {filename}");
                Sample = br.ReadInt32();
                if (sample > -1 && Sample != sample) throw new Exception($"Sample in this horizon file doesn't match the Sample in the filename: {filename}");
                Width = br.ReadInt32();
                Height = br.ReadInt32();
                Step = br.ReadInt32();
                ObserverHeightInMeters = br.ReadInt32();
                //if (observer_height > -1 && ObserverHeight != observer_height) throw new Exception($"ObserverHeight in this horizon file doesn't match the ObserverHeight in the filename: {filename}");
                var flag = br.ReadBoolean();
                if (flag)
                {
                    Horizons = new Horizon[Height][];
                    for (var h = 0; h < Height; h++)
                    {
                        Horizons[h] = new Horizon[Width];
                        for (var w = 0; w < Width; w++)
                        {
                            var horizon = new Horizon { IsDegrees = true };  // Always stored in degrees when on disk
                            var buf = horizon.Buffer;
                            for (var z = 0; z < Horizon.HorizonSamples; z++)
                                buf[z] = br.ReadSingle();
                            Horizons[h][w] = horizon;
                        }
                    }
                }
            }
            /*
            var spath = System.IO.Path.ChangeExtension(path, "json");
            if (File.Exists(spath))
            {
                using (var sr = new StreamReader(spath))
                using (var jr = new Newtonsoft.Json.JsonTextReader(sr))
                {
                    var serializer = new Newtonsoft.Json.JsonSerializer();
                    ShadowCasters = serializer.Deserialize<List<Point>>(jr);
                }
            }*/
        }

        public static TerrainPatch ReadFrom(string path)
        {
            var p = new TerrainPatch { Path = path };
            p.Read(path);
            return p;
        }

        /// <summary>
        /// A shell is a terrain patch with no data
        /// </summary>
        /// <param name="path"></param>
        /// <returns></returns>
        public static TerrainPatch ReadShell(string path)
        {
            var filename = System.IO.Path.GetFileName(path);
            if (!TryParseFilename(filename, out int line, out int sample, out int observer_height))
                return null;
            return new TerrainPatch { Path = path, Line = line, Sample = sample, Height = DefaultSize, Width = DefaultSize, ObserverHeightInMeters = observer_height };
        }

        /// <summary>
        /// Assume the current terrain patch is a shell with a valid path.  Load the data.
        /// </summary>
        public void Load(string path = null)
        {
            if (path == null)
                path = Path;
            if (!File.Exists(path))
                return;
            var p = ReadFrom(path);
            Line = p.Line;
            Sample = p.Sample;
            Width = p.Width;
            Height = p.Height;
            Step = p.Step;
            ObserverHeightInMeters = p.ObserverHeightInMeters;
            //Points = p.Points;
            //Matrices = p.Matrices;
            Horizons = p.Horizons;
            ShadowCasters = p.ShadowCasters;
            //TestVectors = p.TestVectors;
            p.Points = null;
            p.Matrices = null;
            p.Horizons = null;
        }

        public void Write(string path = null)
        {
            path = path ?? Path;
            using (var fs = File.Create(path))
            //using (var gz = new GZipStream(fs, CompressionMode.Compress))
            using (var bw = new BinaryWriter(fs))
            {
                bw.Write(Line);
                bw.Write(Sample);
                bw.Write(Width);
                bw.Write(Height);
                bw.Write(Step);
                bw.Write(ObserverHeightInMeters);
                bw.Write(Horizons != null);
                if (Horizons != null)
                {
                    for (var y = 0; y < Height; y++)
                    {
                        for (var x = 0; x < Width; x++)
                        {
                            var horizons = Horizons[y][x];
                            var buf = horizons.GetSlopeDegreeBuffer();
                            Debug.Assert(buf.Length == Horizon.HorizonSamples);
                            for (var z = 0; z < Horizon.HorizonSamples; z++)
                                bw.Write(buf[z]);
                        }
                    }
                }
            }
            if (ShadowCasters.Count > 0)
            {
                var spath = System.IO.Path.ChangeExtension(path, "json");
                using (var sr = new StreamWriter(spath))
                using (var jr = new Newtonsoft.Json.JsonTextWriter(sr))
                {
                    var serializer = new Newtonsoft.Json.JsonSerializer();
                    serializer.Serialize(jr, ShadowCasters);
                }
            }
        }

        #endregion

        #region Drawing

        public void Draw(Graphics g, DisplayTransform t, Pen pen)
        {
            var p1 = t[new PointF(Sample - 0.5f, Line - 0.5f)];
            var p2 = t[new PointF(Sample + Step * Width + 1f, Line + Step * Height + 1f)];
            g.DrawRectangle(pen, p1.X, p1.Y, p2.X - p1.X, p2.Y - p1.Y);
        }

        public void Draw(Graphics g, DisplayTransform t, Brush brush)
        {
            var p1 = t[new PointF(Sample - 0.5f, Line - 0.5f)];
            var p2 = t[new PointF(Sample + Step * Width + 1f, Line + Step * Height + 1f)];
            g.FillRectangle(brush, p1.X, p1.Y, p2.X - p1.X, p2.Y - p1.Y);
        }

        public void Draw(Graphics g, DisplayTransform t, Bitmap b)
        {
            var p1 = t[new PointF(Sample - 0.5f, Line - 0.5f)];
            var p2 = t[new PointF(Sample + Step * Width + 1f, Line + Step * Height + 1f)];
            g.DrawImage(b, p1.X, p1.Y, p2.X - p1.X, p2.Y - p1.Y);
        }

        public bool IsVisible(DisplayTransform t, Rectangle window)
        {
            var p1 = t[new PointF(Sample - 0.5f, Line - 0.5f)];
            var p2 = t[new PointF(Sample + Width + 1f, Line + Height + 1f)];
            var patchRectangle = new Rectangle((int)p1.X, (int)p1.Y, (int)(p2.X - p1.X), (int)(p2.Y - p1.Y));
            patchRectangle.Inflate(2, 2);
            return window.IntersectsWith(patchRectangle);
        }

        public Rectangle PixelBox(DisplayTransform t)
        {
            var p1 = t[new PointF(Sample, Line)];
            var p2 = t[new PointF(Sample + Width, Line + Height)];
            return new Rectangle((int)(p1.X + 0.5f), (int)(p1.Y + 0.5f), (int)(p2.X - p1.X + 0.5f), (int)(p2.Y - p1.Y + 0.5f));
        }

        #endregion

        #region Rendering Hillshade

        public unsafe Bitmap GetHillshade()
        {
            if (_lastShadowAzimuth == 1000f)  // This constant means hillshade was the last render
                return SunShadowBitmap;

            const double elevation = 30d * Math.PI / 180d;
            const double azimuth = -Math.PI / 2;
            var Zenith_rad = Math.PI / 2d - elevation;
            var Azimuth_rad = Math.PI - azimuth;  // Convert from map azimuth to hillshade azimuth

            var terrain = ViperEnvironment.Terrain;
            var dem = new short[Height, Width];
            for (var iy = 0; iy < Height; iy++)
            {
                for (var ix = 0; ix < Width; ix++)
                    dem[iy, ix] = terrain.Fetch(Line + iy, Sample + ix);
            }

            EnsureSunShadowBitmap();
            var bmp = SunShadowBitmap;
            var width = bmp.Width;
            var height = bmp.Height;
            var bmpData = bmp.LockBits(new Rectangle(0, 0, width, height), System.Drawing.Imaging.ImageLockMode.WriteOnly, bmp.PixelFormat);
            for (var y = 0; y < height; y++)
            {
                var rowptr = (byte*)(bmpData.Scan0 + y * bmpData.Stride);
                for (var x = 0; x < width; x++)
                {
                    SlopesForHillshade(x, y, dem, Width, Height, 20f, out float dzdx, out float dzdy);

                    double Aspect_rad;
                    if (dzdx != 0f)
                    {
                        Aspect_rad = Math.Atan2(dzdy, -dzdx);
                        if (Aspect_rad < 0f)
                            Aspect_rad = Math.PI * 2d + Aspect_rad;
                    }
                    else
                    {
                        if (dzdy > 0f)
                            Aspect_rad = Math.PI / 2d;
                        else if (dzdy < 0)
                            Aspect_rad = 1.5d * Math.PI;
                        else
                            Aspect_rad = 0d;
                    }

                    // Should be this.  There's some inversion in the slopes (I think)
                    //const double z_factor = 1d;
                    //var Slope_rad = Math.Atan(z_factor * Math.Sqrt(dzdx * dzdx + dzdy * dzdy));

                    var Slope_rad = -Math.Atan(Math.Sqrt(dzdx * dzdx + dzdy * dzdy));
                    var Hillshade = (Math.Cos(Zenith_rad) * Math.Cos(Slope_rad)) + (Math.Sin(Zenith_rad) * Math.Sin(Slope_rad) * Math.Cos(Azimuth_rad - Aspect_rad));

                    Hillshade = Math.Max(0d, Hillshade);
                    *(rowptr++) = (byte)(255 * Hillshade);
                }
            }
            bmp.UnlockBits(bmpData);

            _lastShadowAzimuth = 1000f;
            return SunShadowBitmap;
        }

        internal static void SlopesForHillshade(int x, int y, short[,] dem, int width, int height, float cellsize, out float dzdx, out float dzdy)
        {
            var xp = Math.Min(x + 1, width - 1);
            var xm = Math.Max(x - 1, 0);
            var yp = Math.Min(y + 1, height - 1);
            var ym = Math.Max(y - 1, 0);

            var a = dem[ym, xm];
            var b = dem[ym, x];
            var c = dem[ym, xp];
            var d = dem[y, xm];
            var e = dem[y, x];
            var f = dem[y, xp];
            var g = dem[yp, xm];
            var h = dem[yp, x];
            var i = dem[yp, xp];

            dzdx = ((c + f + f + i) - (a + d + d + g)) / (8f * cellsize);
            dzdy = ((g + h + h + i) - (a + b + b + c)) / (8f * cellsize);
        }
        #endregion

        #region Rendering Shadows

        public Bitmap GetShadows(Vector3d shadow_caster)
        {
            GetAzEl(shadow_caster, out float azimuth_rad, out float elevation_rad);
            return GetShadows(azimuth_rad * 180f / 3.141592653589f, elevation_rad * 180f / 3.141592653589f);
        }

        protected void EnsureSunShadowBitmap()
        {
            if (SunShadowBitmap != null)
                return;
            SunShadowBitmap = new Bitmap(Width, Height, PixelFormat.Format8bppIndexed);
            var palette = SunShadowBitmap.Palette;
            for (var i = 0; i < 256; i++)
                palette.Entries[i] = Color.FromArgb(i, i, i);
            SunShadowBitmap.Palette = palette;
        }

        public Bitmap GetShadows(float azimuth_deg, float elevation_deg)
        {
            EnsureSunShadowBitmap();
            if (Math.Abs(_lastShadowAzimuth - azimuth_deg) > 0.001f || Math.Abs(_lastShadowElevation - elevation_deg) > 0.001f)
            {
                RenderShadows(SunShadowBitmap, azimuth_deg, elevation_deg);
                _lastShadowAzimuth = azimuth_deg;
                _lastShadowElevation = elevation_deg;
                Console.WriteLine($"render shadows for {Id} az={_lastShadowAzimuth} el={_lastShadowElevation}");
            }
            return SunShadowBitmap;
        }

        internal unsafe void RenderShadows(Bitmap bmp, float azimuth_deg, float elevation_deg)
        {
            if (bmp == null || Horizons == null)
                return;
            var width = bmp.Width;
            var height = bmp.Height;
            var bmpData = bmp.LockBits(new Rectangle(0, 0, width, height), System.Drawing.Imaging.ImageLockMode.WriteOnly, bmp.PixelFormat);
            for (var row = 0; row < height; row++)
            {
                var rowptr = (byte*)(bmpData.Scan0 + row * bmpData.Stride);
                for (var col = 0; col < width; col++)
                {
                    var horizon = Horizons[row][col];
                    rowptr[col] = (byte)(255f * horizon.SunFraction(azimuth_deg, elevation_deg));
                }
            }
            bmp.UnlockBits(bmpData);
        }

        /// <summary>
        /// Generate the azimuth and elevation of a point expressed in rectangular coordinates in Mean Earth frame
        /// </summary>
        /// <param name="point_me"></param>
        /// <param name="azimuth_rad"></param>
        /// <param name="elevation_rad"></param>
        public void GetAzEl(Vector3d point_me, out float azimuth_rad, out float elevation_rad)
        {
            if (!_centerMatrixIsValid)
            {
                var h = Height / 2;
                var w = Width / 2;
                int line = Line + h;
                var sample = Sample + w;
                var relz = ViperEnvironment.Terrain.LineSampleToTerrainOffset(line, sample);   //Note: Ugh!  Bad modularity
                var radius = MoonRadius + relz / 1000d;
                InMemoryInt16Terrain.GetLatLon(line, sample, out double lat, out double lon);
                var z = radius * Math.Sin(lat);
                var c = radius * Math.Cos(lat);
                var x = c * Math.Cos(lon);  // TODO: Not sure about these
                var y = c * Math.Sin(lon);
                var vec = new Vector3d(x, y, z);
                _centerMatrix = MatrixFromLatLon(lat, lon, vec);
            }
            {
                var temp = new Vector3d();
                Transform(ref point_me, ref _centerMatrix, ref temp);
                var x = temp[0];
                var y = temp[1];
                var z = temp[2];
                azimuth_rad = (float)(Math.Atan2(y, x) + Math.PI);  // [0,2PI]
                var alen = Math.Sqrt(x * x + y * y);
                elevation_rad = (float)Math.Atan2(z, alen);
            }
        }

        Vector3d _last_sun_shadow_caster;

        public unsafe Bitmap GetSunShadowsCarefully(Vector3d shadow_caster)
        {
            if (SunShadowBitmap == null)
            {
                SunShadowBitmap = new Bitmap(DefaultSize, DefaultSize, System.Drawing.Imaging.PixelFormat.Format8bppIndexed);
                var palette = SunShadowBitmap.Palette;
                for (var i = 0; i < 256; i++)
                    palette.Entries[i] = Color.FromArgb(i, i, i);
                SunShadowBitmap.Palette = palette;
            }
            if (Horizons == null)
                return SunShadowBitmap;
            if (_last_sun_shadow_caster.Equals(shadow_caster))
                return SunShadowBitmap;
            var bmp = SunShadowBitmap;
            var width = bmp.Width;
            var height = bmp.Height;
            var bmpData = bmp.LockBits(new Rectangle(0, 0, width, height), System.Drawing.Imaging.ImageLockMode.WriteOnly, bmp.PixelFormat);
            for (var row = 0; row < height; row++)
            {
                var rowptr = (byte*)(bmpData.Scan0 + row * bmpData.Stride);
                for (var col = 0; col < width; col++)
                {
                    var horizon = Horizons[row][col];
                    GetAzEl(shadow_caster, col, row, out float azimuth_rad, out float elevation_rad);
                    rowptr[col] = (byte)(255f * horizon.SunFraction2(azimuth_rad * 180f / 3.141592653589f, elevation_rad * 180f / 3.141592653589f));
                }
            }
            bmp.UnlockBits(bmpData);

            _last_sun_shadow_caster = shadow_caster;
            return SunShadowBitmap;
        }

        Vector3d _last_earth_shadow_caster;

        public unsafe Bitmap GetEarthShadows(Vector3d shadow_caster)
        {
            if (EarthShadowBitmap == null)
            {
                EarthShadowBitmap = new Bitmap(DefaultSize, DefaultSize, System.Drawing.Imaging.PixelFormat.Format8bppIndexed);
                var palette = EarthShadowBitmap.Palette;
                palette.Entries[0] = Color.Black;
                palette.Entries[1] = Color.DarkRed;
                palette.Entries[2] = Color.IndianRed;
                palette.Entries[3] = Color.Gold;
                palette.Entries[4] = Color.Gold;
                palette.Entries[5] = Color.Yellow;
                palette.Entries[6] = Color.Yellow;
                palette.Entries[7] = Color.PaleGoldenrod;
                palette.Entries[8] = Color.PaleGoldenrod;
                for (var i = 9; i < 256; i++)
                    palette.Entries[i] = Color.White;
                EarthShadowBitmap.Palette = palette;
            }
            if (Horizons == null)
                return SunShadowBitmap;
            if (_last_earth_shadow_caster.Equals(shadow_caster))
                return EarthShadowBitmap;
            var bmp = EarthShadowBitmap;
            var width = bmp.Width;
            var height = bmp.Height;
            var bmpData = bmp.LockBits(new Rectangle(0, 0, width, height), System.Drawing.Imaging.ImageLockMode.WriteOnly, bmp.PixelFormat);
            for (var row = 0; row < height; row++)
            {
                var rowptr = (byte*)(bmpData.Scan0 + row * bmpData.Stride);
                for (var col = 0; col < width; col++)
                {
                    var horizon = Horizons[row][col];

                    GetAzEl(shadow_caster, col, row, out float azimuth_rad, out float elevation_rad);
                    Debug.Assert(azimuth_rad >= 0f && azimuth_rad <= 2f * 3.141592653589f);

                    var elevation_deg = elevation_rad * 180f / 3.141592653589f;
                    var over_horizon = horizon.OverHorizon(azimuth_rad, elevation_rad * 180f / 3.141592653589f);
                    var color_index = over_horizon <= 0f ? 0 : (int)Math.Ceiling(over_horizon);
                    rowptr[col] = (byte)color_index;
                }
            }
            bmp.UnlockBits(bmpData);

            _last_earth_shadow_caster = shadow_caster;
            return EarthShadowBitmap;
        }

        public void GetAzEl(Vector3d point_me, int ix, int iy, out float azimuth_rad, out float elevation_rad)
        {
            var temp = new Vector3d();
            Transform(ref point_me, ref Matrices[iy][ix], ref temp);
            var x = temp[0];
            var y = temp[1];
            var z = temp[2];
            azimuth_rad = (float)(Math.Atan2(y, x) + Math.PI);  // [0,2PI]
            var alen = Math.Sqrt(x * x + y * y);
            elevation_rad = (float)Math.Atan2(z, alen);
        }

        #endregion

        #region Average Sun and Earth Patches

        public unsafe Bitmap GenerateAverageSun(List<Vector3d> vectors = null)
        {
            FillPointsAndMatrices(ViperEnvironment.Terrain);
            var side = DefaultSize;
            var maxsun_scale = 255f / (vectors.Count * 255);
            var bitmap = new Bitmap(Width, Height, PixelFormat.Format8bppIndexed);
            var palette = bitmap.Palette;
            for (var i = 0; i < 256; i++)
                palette.Entries[i] = Color.FromArgb(i, i, i);
            bitmap.Palette = palette;
            var bitmap_rect = new Rectangle(0, 0, Width, Height);
            var bmpData = bitmap.LockBits(bitmap_rect, ImageLockMode.WriteOnly, bitmap.PixelFormat);
            for (var row = 0; row < side; row++)
            {
                //Console.WriteLine(row);
                var rowptr = (byte*)(bmpData.Scan0 + row * bmpData.Stride);
                for (var col = 0; col < side; col++)
                {
                    var horizon = Horizons[row][col];
                    if (!horizon.IsDegrees)
                        horizon.ConvertSlopeToDegrees();
                    var sum = 0;
                    //var time = 0;
                    for (var svi = 0; svi < vectors.Count; svi++)
                    {
                        GetAzEl(vectors[svi], col, row, out float sun_azimuth_rad, out float sun_elevation_rad);
                        var val = (int)(255f * horizon.SunFraction2(sun_azimuth_rad * 180f / 3.141592653589f, sun_elevation_rad * 180f / 3.141592653589f));
                        sum += val;
                        //time += Math.Sign(val);
                    }
                    rowptr[col] = (byte)(maxsun_scale * sum);
                }
            }
            bitmap.UnlockBits(bmpData);
            return bitmap;
        }

        public unsafe Bitmap GenerateAverageEarth(List<Vector3d> vectors, float threshold = 2)
        {
            FillPointsAndMatrices(ViperEnvironment.Terrain);
            var side = DefaultSize;
            var maxearth_scale = 255f / vectors.Count;
            var bitmap = new Bitmap(Width, Height, PixelFormat.Format8bppIndexed);
            var palette = bitmap.Palette;
            for (var i = 0; i < 256; i++)
                palette.Entries[i] = Color.FromArgb(i, i, i);
            bitmap.Palette = palette;
            var bitmap_rect = new Rectangle(0, 0, Width, Height);
            var bmpData = bitmap.LockBits(bitmap_rect, ImageLockMode.WriteOnly, bitmap.PixelFormat);
            for (var row = 0; row < side; row++)
            {
                var rowptr = (byte*)(bmpData.Scan0 + row * bmpData.Stride);
                for (var col = 0; col < side; col++)
                {
                    var horizon = Horizons[row][col];
                    if (!horizon.IsDegrees)
                        horizon.ConvertSlopeToDegrees();
                    var sum = 0;
                    for (var svi = 0; svi < vectors.Count; svi++)
                    {
                        GetAzEl(vectors[svi], col, row, out float azimuth_rad, out float elevation_rad);
                        var over_horizon = horizon.OverHorizon(azimuth_rad, elevation_rad * 180f / 3.141592653589f);
                        if (over_horizon >= threshold)
                            sum++;
                    }
                    rowptr[col] = (byte)(maxearth_scale * sum);
                }
            }
            bitmap.UnlockBits(bmpData);
            return bitmap;
        }

        #endregion

        public void MapOverLightCurve(int row, int col, List<Vector3d> sun_vectors, Action<int, float> action)
        {
            FillPointsAndMatrices(InMemoryInt16Terrain.Singleton);

            var mat = Matrices[row][col];
            var horizon = Horizons[row][col];
            if (!horizon.IsDegrees)
                lock (horizon)
                    horizon.ConvertSlopeToDegrees();

            var temp = new Vector3d();
            for (var i = 0; i < sun_vectors.Count; i++)
            {
                var sunvec = sun_vectors[i];
                Transform(ref sunvec, ref mat, ref temp);

                var sun_x = temp[0];
                var sun_y = temp[1];
                var sun_z = temp[2];
                var azimuth_rad = Math.Atan2(sun_y, sun_x) + Math.PI;  // [0,2PI]
                var azimuth_deg = (float)(azimuth_rad * 180d / Math.PI);

                var alen = Math.Sqrt(sun_x * sun_x + sun_y * sun_y);
                var slope = sun_z / alen;
                var elevation_deg = ((float)Math.Atan(slope)) * 180f / 3.141592653589f;
                var sun = horizon.SunFraction2(azimuth_deg, elevation_deg);

                action(i, sun);
            }
        }

        public float[] GetLightCurve(int row, int col, List<Vector3d> sun_vectors, float[] result = null)
        {
            if (result == null)
                result = new float[sun_vectors.Count];
            MapOverLightCurve(row, col, sun_vectors, (i, sun) => result[i] = sun);
            return result;
        }

        public void MapOverEarthElevation(int row, int col, List<Vector3d> vectors, Action<int, float> action)
        {
            FillPointsAndMatrices(InMemoryInt16Terrain.Singleton);

            var mat = Matrices[row][col];
            var horizon = Horizons[row][col];
            if (!horizon.IsDegrees)
                lock (horizon)
                    horizon.ConvertSlopeToDegrees();

            for (var i = 0; i < vectors.Count; i++)
            {
                GetAzEl(vectors[i], col, row, out float azimuth_rad, out float elevation_rad);
                var over_horizon = horizon.OverHorizon(azimuth_rad, elevation_rad * 180f / 3.141592653589f);
                action(i, over_horizon);
            }
        }

        public float[] GetEarthOverHorizonCurve(int row, int col, List<Vector3d> earth_vectors, float[] result = null)
        {
            if (result == null)
                result = new float[earth_vectors.Count];
            MapOverEarthElevation(row, col, earth_vectors, (i, over_horizon) => result[i] = over_horizon);
            return result;
        }

        #region Utilities

        internal void Describe(string point_name, int line, int sample)
        {
            double latd, lond;
            InMemoryInt16Terrain.GetLatLon(line, sample, out double latr, out double lonr);
            latd = latr * 180d / Math.PI;
            lond = lonr * 180d / Math.PI;
            Console.WriteLine($"{point_name} [l,s]=[{line},{sample}] [latd,lond]=[{latd},{lond}]");
        }

        #endregion

        #region IDisposable Support
        private bool disposedValue = false; // To detect redundant calls

        protected virtual void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (disposing)
                {
                    SunShadowBitmap?.Dispose();
                    SunShadowBitmap = null;
                    EarthShadowBitmap?.Dispose();
                    EarthShadowBitmap = null;
                }

                // TODO: free unmanaged resources (unmanaged objects) and override a finalizer below.
                // TODO: set large fields to null.

                disposedValue = true;
            }
        }

        // TODO: override a finalizer only if Dispose(bool disposing) above has code to free unmanaged resources.
        // ~TerrainPatch() {
        //   // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
        //   Dispose(false);
        // }

        // This code added to correctly implement the disposable pattern.
        public void Dispose()
        {
            // Do not change this code. Put cleanup code in Dispose(bool disposing) above.
            Dispose(true);
            // TODO: uncomment the following line if the finalizer is overridden above.
            // GC.SuppressFinalize(this);
        }
        #endregion

        #region Format Conversion

        public static void ConvertFile(string path)
        {
            var filename = System.IO.Path.GetFileName(path);
            using (var ms = new MemoryStream())
            {
                using (FileStream file = new FileStream(path, FileMode.Open, FileAccess.Read))
                    file.CopyTo(ms);
                ms.Position = 0;
                var patch = new TerrainPatch();
                patch.OldRead(ms);
                patch.Write();
            }
        }

        public void OldRead(string path)
        {
            path = path ?? Path;
            using (var fs = File.OpenRead(path))
                OldRead(fs);
        }

        public void OldRead(Stream stream)
        {
            using (var br = new BinaryReader(stream))
            {
                Line = br.ReadInt32();
                Sample = br.ReadInt32();
                Width = br.ReadInt32();
                Height = br.ReadInt32();
                Step = br.ReadInt32();
                var flag = br.ReadBoolean();
                if (flag)
                {
                    Horizons = new Horizon[Height][];
                    for (var h = 0; h < Height; h++)
                    {
                        Horizons[h] = new Horizon[Width];
                        for (var w = 0; w < Width; w++)
                        {
                            var horizon = new Horizon { IsDegrees = true };  // Always stored in degrees when on disk
                            var buf = horizon.Buffer;
                            for (var z = 0; z < Horizon.HorizonSamples; z++)
                                buf[z] = br.ReadSingle();
                            Horizons[h][w] = horizon;
                        }
                    }
                }
            }
            // There was code here to read ShadowCasters.  I'm ignoring that.
        }

        #endregion
    }
}
