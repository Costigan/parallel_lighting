using OSGeo.GDAL;
using System;
using System.Diagnostics;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;

namespace viper.corelib.utilities
{
    public class GeotiffHelper
    {
        public static void WriteArrayAsGeotiff(float[,] ary, Rectangle r, string path)
        {
            Debug.Assert(r.Height == ary.GetLength(0) && r.Width == ary.GetLength(1));
            if (!GetGdalDrivers())
                return;
            if (File.Exists(path))
                File.Delete(path);
            var height = r.Height;
            var width = r.Width;
            using (var dst_ds = _geotif_driver.Create(path, width, height, 1, DataType.GDT_Float32, null))
            {
                var band = dst_ds.GetRasterBand(1);
                var buf = new float[width];
                for (var row = 0; row < height; row++)
                {
                    //Buffer.BlockCopy(dem, row * width * 4, buf, 0, width * 4);
                    for (var col = 0; col < width; col++)
                        buf[col] = ary[row, col];
                    band.WriteRaster(0, row, width, 1, buf, width, height, 0, 0);
                }
                var transform = GetAffineTransformOfBounds(r);
                dst_ds.SetGeoTransform(transform);
                dst_ds.SetProjection(dst_wkt);
            }
        }

        public static string dst_wkt = ToWkt(@"PROJCS['unnamed',GEOGCS['unnamed ellipse',DATUM['unknown',SPHEROID['unnamed',1737400,0]],PRIMEM['Greenwich',0],UNIT['degree',0.0174532925199433]],PROJECTION['Polar_Stereographic'],PARAMETER['latitude_of_origin',-90],PARAMETER['central_meridian',0],PARAMETER['scale_factor',1],PARAMETER['false_easting',0],PARAMETER['false_northing',0],UNIT['metre',1,AUTHORITY['EPSG','9001']]]");

        public static PointD PixelToCoordinate(int line, int sample)
        {
            const double S0 = 15199.5d;             // PDS SAMPLE_PROJECTION_OFFSET
            const double L0 = 15199.5d;             // PDS LINE_PROJECTION_OFFSET
            const double Scale = 20d;

            var x = (sample - S0) * Scale;
            var y = (L0 - line) * Scale;
            return new PointD(x, y);
        }

        public static RectangleD GetGeoExtent(Rectangle r)
        {
            var ul = PixelToCoordinate(r.Top, r.Left);  // line,sample
            var lr = PixelToCoordinate(r.Bottom, r.Right);
            return new RectangleD(ul, lr);
        }

        // [minx, scalex, 0, maxy, 0, -scaley]
        public static double[] GetAffineTransformOfBounds(Rectangle bounds)
        {
            var r = GetGeoExtent(bounds);
            return new double[] { r.Left, 20d, 0d, r.Top, 0d, -20d };
        }

        public static string ToWkt(string s) => s.Replace('\'', '"');

        public unsafe static float[,] ReadGeotiffAsFloatArray(string path, float? no_data = null)
        {
            if (!GetGdalDrivers())
                return null;
            float[,] ary = null;
            using (var ds = Gdal.Open(path, Access.GA_ReadOnly))
            {
                var height = ds.RasterYSize;
                var width = ds.RasterXSize;
                ary = new float[height, width];
                fixed (float* ary_ptr = &ary[0,0])
                {
                    var iptr = new IntPtr(ary_ptr);
                    var band = ds.GetRasterBand(1);
                    band.GetNoDataValue(out double band_no_data, out int hasval);
                    band.ReadRaster(0, 0, width, height, iptr, width, height, DataType.GDT_Float32, 0, 0);

                    if (no_data.HasValue && hasval != 0)
                    {
                        var band_no_dataf = (float)band_no_data;
                        var no_dataf = no_data.Value;
                        for (var row = 0; row < ary.GetLength(0); row++)
                            for (var col = 0; col < ary.GetLength(1); col++)
                                if (ary[row, col] == band_no_dataf)
                                    ary[row, col] = no_dataf;
                    }
                }
            }
            return ary;
        }

        #region GDAL Driver Support

        static Driver _memory_driver = null;
        static Driver _geotif_driver = null;

        public static bool GetGdalDrivers()
        {
            if (_memory_driver != null)
                return true;
            GdalConfiguration.ConfigureGdal();
            const string gdalFormat = "GTiff";
            _memory_driver = Gdal.GetDriverByName("MEM");
            _geotif_driver = Gdal.GetDriverByName(gdalFormat);
            if (!(_memory_driver != null && _geotif_driver != null))
            {
                Console.WriteLine($"Can't load {(_geotif_driver == null ? "GTiff" : "MEM")} driver.");
                if (_memory_driver != null) _memory_driver.Dispose();
                if (_geotif_driver != null) _geotif_driver.Dispose();
                _memory_driver = _geotif_driver = null;
                return false;
            }
            return true;
        }

        public static Driver GetGeoTIFFDriver()
        {
            if (GetGdalDrivers())
                return _geotif_driver;
            throw new Exception(@"Couldn't open both GDAL Drivers");
        }

        public static Driver GetMemoryDriver()
        {
            if (GetGdalDrivers())
                return _memory_driver;
            throw new Exception(@"Couldn't open both GDAL Drivers");
        }


        #endregion

        #region Bitmaps

        public static void WriteBitmapAsGeotiff(Bitmap bmp, Rectangle bounds, string target_path)
        {
            if (!GetGdalDrivers())
                return;

            // set up MEM driver to read bitmap data
            int bandCount = 1;
            int pixelOffset = 1;
            DataType dataType = DataType.GDT_Byte;
            int[] band2ptroffset = null;
            switch (bmp.PixelFormat)
            {
                case PixelFormat.Format16bppGrayScale:
                    dataType = DataType.GDT_Int16;
                    bandCount = 1;
                    pixelOffset = 2;
                    break;
                case PixelFormat.Format24bppRgb:
                    dataType = DataType.GDT_Byte;
                    bandCount = 3;
                    pixelOffset = 3;
                    break;
                case PixelFormat.Format32bppArgb:
                    dataType = DataType.GDT_Byte;
                    //bandCount = 4;  // This was the band count, but the order is wrong
                    bandCount = 3;
                    band2ptroffset = new int[] { 3, 2, 1, 0 };
                    pixelOffset = 4;
                    break;
                case PixelFormat.Format48bppRgb:
                    dataType = DataType.GDT_UInt16;
                    bandCount = 3;
                    pixelOffset = 6;
                    break;
                case PixelFormat.Format64bppArgb:
                    dataType = DataType.GDT_UInt16;
                    bandCount = 4;
                    pixelOffset = 8;
                    break;
                case PixelFormat.Format8bppIndexed:
                    dataType = DataType.GDT_Byte;
                    bandCount = 1;
                    pixelOffset = 1;
                    break;
                default:
                    Console.WriteLine("Invalid pixel format " + bmp.PixelFormat.ToString());
                    break;
            }

            // Use GDAL raster reading methods to read the image data directly into the Bitmap
            BitmapData bitmapData = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);

            int stride = bitmapData.Stride;
            IntPtr buf = bitmapData.Scan0;

            try
            {
                Dataset ds = _memory_driver.Create("", bmp.Width, bmp.Height, 0, dataType, null);
                var transform = GetAffineTransformOfBounds(bounds);
                ds.SetGeoTransform(transform);
                ds.SetProjection(ToWkt(dst_wkt));
                //var color_interp = new ColorInterp[] { ColorInterp.GCI_AlphaBand, ColorInterp.GCI_RedBand, ColorInterp.GCI_GreenBand, ColorInterp.GCI_BlueBand };
                // add bands in a reverse order
                for (int i = 1; i <= bandCount; i++)
                {
                    // original
                    CPLErr err;
                    if (band2ptroffset == null)
                        err = ds.AddBand(dataType, new string[] { "DATAPOINTER=" + Convert.ToString(buf.ToInt64() + bandCount - i), "PIXELOFFSET=" + pixelOffset, "LINEOFFSET=" + stride });
                    else
                        err = ds.AddBand(dataType, new string[] { "DATAPOINTER=" + Convert.ToString(buf.ToInt64() + band2ptroffset[i]), "PIXELOFFSET=" + pixelOffset, "LINEOFFSET=" + stride });
                    if (err != CPLErr.CE_None)
                        Console.WriteLine($"Error adding band to {target_path}");
                }

                if (bmp.PixelFormat == PixelFormat.Format8bppIndexed)
                {
                    var ct1 = new ColorTable(PaletteInterp.GPI_RGB);
                    var band = ds.GetRasterBand(1);
                    var ct = band.GetRasterColorTable();
                    Console.WriteLine($"color table count={ct.GetCount()}");
                }

                if (File.Exists(target_path))
                    File.Delete(target_path);
                var ds_copy = _geotif_driver.CreateCopy(target_path, ds, 0, null, null, null);
                ds_copy.FlushCache();
                ds_copy.Dispose();
                ds.Dispose();
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.Message);
            }
            finally
            {
                bmp.UnlockBits(bitmapData);
            }
        }

        #endregion
    }
}
