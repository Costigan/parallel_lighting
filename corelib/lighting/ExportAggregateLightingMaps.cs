using OSGeo.GDAL;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Threading.Tasks;
using viper.corelib.math;
using viper.corelib.patch;
using viper.corelib.spice;
using viper.corelib.utilities;

namespace viper.corelib.lighting
{
    public class ExportAggregateLightingMaps
    {
        public static string Directory = ".";

        public Rectangle Bounds { get; set; }
        public Rectangle IdBounds { get; set; }

        List<TerrainPatch> _patches;

        public void DoExport()
        {
            GdalConfiguration.ConfigureGdal();
            var b = Bounds;
            Console.WriteLine($"Bounds: {b.Left} {b.Top} {b.Right} {b.Bottom}");
            var c = GetGeoExtent(b);
            Console.WriteLine($"Bounds: {c.Left} {c.Top} {c.Right} {c.Bottom}");
            var s = GetGdalTranslateCommand();
            Console.WriteLine(s);

            if (false)
            {
                Hack();
                return;
            }

            ReprojectionTest();
            return;
            LoadTiles();
            GeneratePSRImage();
        }

        void LoadTiles()
        {
            var s = TerrainPatch.DefaultSize;
            IdBounds = new Rectangle(Bounds.Left / s, Bounds.Top / s, Bounds.Width / s, Bounds.Height / s);
            _patches = new List<TerrainPatch>();
            for (var y = IdBounds.Top; y < IdBounds.Bottom; y++)
                for (var x = IdBounds.Left; x < IdBounds.Right; x++)
                {
                    var p = ViperEnvironment.GetPatch(new Point(x, y));
                    if (p != null)
                        _patches.Add(p);
                }
            Console.WriteLine($"{_patches.Count} patches were loaded");
            foreach (var p in _patches)
                Console.WriteLine(p.Id);
        }

        unsafe void GeneratePSRImage()
        {
            var side = TerrainPatch.DefaultSize;

            // Prepare bitmap
            var width = side * IdBounds.Width;
            var height = side * IdBounds.Height;
            var sum_image = new Int32[height, width];
            var time_image = new Int32[height, width];
            var bitmap_rect = new Rectangle(0, 0, width, height);

            // Generate sun vectors
            var start_time = new DateTime(2000, 1, 1);
            var stop_time = new DateTime(2019, 1, 1);
            var time_step = TimeSpan.FromHours(2);
            var sun_vectors = new List<Vector3d>();
            for (var time = start_time; time <= stop_time; time += time_step)
                sun_vectors.Add(CSpice.SunPosition(time));

            Console.WriteLine(sun_vectors.Count);

            System.Diagnostics.Debug.Assert(sun_vectors.Count <= 256 * 256 * 256);

            // Calculate psr image; parallelize over patches
            var options = new ParallelOptions { MaxDegreeOfParallelism = ViperEnvironment.MaxDegreeOfParallelism };  //TODO: Environment.ProcessorCount - 2 isn't working
            Parallel.ForEach(_patches, options,
                patch =>
                {
                    patch.FillPointsAndMatrices(ViperEnvironment.Terrain);
                    var base_row = (patch.Id.Y - IdBounds.Top) * side;
                    var base_col = (patch.Id.X - IdBounds.Left) * side;
                    for (var row = 0; row < side; row++)
                    {
                        var row1 = row + base_row;
                        for (var col = 0; col < side; col++)
                        {
                            var horizon = patch.Horizons[row][col];
                            var sum = 0;
                            var time = 0;
                            for (var svi = 0;svi<sun_vectors.Count;svi++)
                            {
                                var sun = sun_vectors[svi];
                                patch.GetAzEl(sun, col, row, out float sun_azimuth_rad, out float sun_elevation_rad);
                                var val = (int)(255f * horizon.SunFraction2(sun_azimuth_rad * 180f / 3.141592653589f, sun_elevation_rad * 180f / 3.141592653589f));
                                sum += val;
                                time += Math.Sign(val);
                            }
                            sum_image[row1,col + base_col] = sum;
                            time_image[row1, col + base_col] = time;
                        }
                    }
                });
        }

        string GetGdalTranslateCommand()
        {
            var b = Bounds;
            var c = GetGeoExtent(b);
            var xmin = Math.Min(c.Left, c.Right);
            var xmax = Math.Max(c.Left, c.Right);
            var ymin = Math.Min(c.Top, c.Bottom);
            var ymax = Math.Max(c.Top, c.Bottom);
            return $"gdalwarp -overwrite -t_srs \"+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=1737400 +b=1737400 +units=m +no_defs\" -ts {b.Width} {b.Height}  -of GTiff -te {xmin} {ymin} {xmax} {ymax} -r bilinear  PermShadows.85.1940S35.4360E.tex.tif temp.tif ";
        }

        public PointD PixelToCoordinate(int line, int sample)
        {
            const double S0 = 15199.5d;             // PDS SAMPLE_PROJECTION_OFFSET
            const double L0 = 15199.5d;             // PDS LINE_PROJECTION_OFFSET
            const double Scale = 20d;

            var x = (sample - S0) * Scale;
            var y = (L0 - line) * Scale;
            return new PointD(x, y);
        }

        RectangleD GetGeoExtent(Rectangle r)
        {
            var ul = PixelToCoordinate(r.Top, r.Left);  // line,sample
            var lr = PixelToCoordinate(r.Bottom, r.Right);
            return new RectangleD(ul, lr);
        }

        // [minx, scalex, 0, maxy, 0, -scaley]
        double[] GetAffineTransformOfBounds()
        {
            var r = GetGeoExtent(Bounds);
            return new double[] { r.Left, 20d, 0d, r.Top, 0d, -20d };
        }

        void ReprojectionTest()
        {
            GdalConfiguration.ConfigureGdal();
            const string gdalFormat = "GTiff";
            var driver = Gdal.GetDriverByName(gdalFormat);
            if (driver == null)
            {
                Console.WriteLine($"Can't load {gdalFormat} driver.  Not writing ice stability depth float tiff");
                return;
            }
            var driverMetadata = driver.GetMetadata(null);
            var canCreate = driverMetadata.Contains("DCAP_CREATE=YES");
            var canCreateCopy = driverMetadata.Contains("DCAP_CREATECOPY=YES");


            var src_path = @"C:\RP\maps\inputs_20m\Nobile_1\PermShadows.85.1940S35.4360E.tex.tif";
            var dst_path = "/vsimem/temp.tif";  //  @"C:\RP\maps\inputs_20m\Nobile_1\tempgen.tif";
            //var srs_proj4 = @"+proj=sterea +lat_0=-85.19257213508449 +lon_0=35.43540608703493 +k=1 +x_0=0 +y_0=0 +a=1737400 +b=1737400 +units=m +no_defs";
            //var dst_proj4 = @"+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=1737400 +b=1737400 +units=m +no_defs";

            var src_wkt = @"PROJCS['Moon Oblique Stereographic',GEOGCS['GCS_MOON',DATUM['D_MOON',SPHEROID['MOON',1737400,0]],PRIMEM['Reference_Meridian',0],UNIT['degree',0.0174532925199433]],PROJECTION['Oblique_Stereographic'],PARAMETER['latitude_of_origin',-85.19257213508449],PARAMETER['central_meridian',35.43540608703493],PARAMETER['scale_factor',1],PARAMETER['false_easting',0],PARAMETER['false_northing',0],UNIT['metre',1,AUTHORITY['EPSG','9001']]]";
            var dst_wkt = @"PROJCS['unnamed',GEOGCS['unnamed ellipse',DATUM['unknown',SPHEROID['unnamed',1737400,0]],PRIMEM['Greenwich',0],UNIT['degree',0.0174532925199433]],PROJECTION['Polar_Stereographic'],PARAMETER['latitude_of_origin',-90],PARAMETER['central_meridian',0],PARAMETER['scale_factor',1],PARAMETER['false_easting',0],PARAMETER['false_northing',0],UNIT['metre',1,AUTHORITY['EPSG','9001']]]";

            using (Dataset src_ds = Gdal.Open(src_path, Access.GA_ReadOnly))
            {
                if (src_ds == null)
                    throw new Exception($"GDAL couldn't open {src_path}");

                //File.Delete(dst_path);
                using (var dst_ds = driver.Create(dst_path, Bounds.Width, Bounds.Height, 1, DataType.GDT_Byte, null))
                {
                    var transform = GetAffineTransformOfBounds();
                    dst_ds.SetGeoTransform(transform);
                    Gdal.ReprojectImage(src_ds, dst_ds, ToWkt(src_wkt), ToWkt(dst_wkt), ResampleAlg.GRA_Bilinear, 0d, 0d, null, null, null);
                }
            }
        }

        void Hack()
        {
            var path = @"C:\RP\maps\inputs_20m\Nobile_1\temp1.tif";
            using (Dataset src_ds = Gdal.Open(path, Access.GA_ReadOnly))
            {
                if (src_ds == null)
                    throw new Exception($"GDAL couldn't open {path}");

                string Projection;
                double[] AffineTransform = new double[6];

                src_ds.GetGeoTransform(AffineTransform);
                Projection = src_ds.GetProjectionRef();
            }
        }

        string ToWkt(string s) => s.Replace('\'', '"');

    }
}
