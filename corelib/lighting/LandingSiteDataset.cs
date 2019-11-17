using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Linq;
using OSGeo.GDAL;
using viper.corelib.horizon;
using viper.corelib.interfaces;
using viper.corelib.math;
using viper.corelib.patch;
using viper.corelib.terrain;
using viper.corelib.utilities;

namespace viper.corelib.lighting
{
    /// <summary>
    /// Landing site datasets for Artemis
    /// </summary>
    public class LandingSiteDataset
    {
        public ILunarHorizon MainWindow;

        public const string SiteDirectory = @"landing_site_datasets";
        public static string SiteRoot;

        public static Dictionary<int, SiteDescriptor> SiteDictionary = new Dictionary<int, SiteDescriptor>() {
            // Southern
            { 0001, new SiteDescriptor{ Id=0001,Lon=222.69, Lat=-89.45, Name= "Mazarico 1"} },
            { 0004, new SiteDescriptor{ Id=0004,Lon=204.27, Lat=-89.78 , Name= "Mazarico 4"} },
            { 0007, new SiteDescriptor{ Id=0007, Lon=123.64, Lat=-88.81, Name= "Mazarico 7" } },
            { 0011, new SiteDescriptor{ Id=0011, Lon=291.90, Lat=-88.67, Name= "Mazarico 11" } },
            { 0101, new SiteDescriptor{ Id=0101, Lon=-35.8187, Lat=-83.6933, Name= "North Cabeus" } },
            { 0102, new SiteDescriptor{ Id=0102, Lon=32.1105, Lat=-85.4446, Name= "Nobile 1" } },
            { 0103, new SiteDescriptor{ Id=0103, Lon=2.0827, Lat=-85.9750, Name= "Malapert 1" } },
            { 0104, new SiteDescriptor{ Id=0104, Lon=-3.1406, Lat=-85.9698, Name= "Malapert 2" } },
            { 0105, new SiteDescriptor{ Id=0105, Lon = 57.1368, Lat = -87.1180, Name= "Shoemaker 1" } },
            { 0106, new SiteDescriptor{ Id=0106, Lon = -92.9704, Lat = -89.2588, Name= "de Gerlache South" } },
            { 0107, new SiteDescriptor{ Id=0107, Lon = -7.3766, Lat = -89.2804, Name= "Shirley1" } },
            { 0108, new SiteDescriptor{ Id=0108, Lon=-11.0677, Lat=-89.3781, Name="Shirley2" } },
            { 0201, new SiteDescriptor{ Id=0201, Lon=223.13621, Lat=-89.44075, Name="Spudis Point" } },
            { 0202, new SiteDescriptor{ Id=0202, Lon=78.2717, Lat=-88.7990, Name= "Shoemaker 2" } },
            { 0109, new SiteDescriptor{ Id=0109, IdRectangle = new Rectangle(149,66,158-149,74-66), Name="Nobile 2" } },
            { 0110, new SiteDescriptor{ Id=0110, IdRectangle = new Rectangle(100,68,108-100,76-68), Name="Haworth 1" } },
            { 0111, new SiteDescriptor{ Id=0111, IdRectangle = new Rectangle(100,80,107-100,87-80), Name="Haworth 2" } },
            { 0203, new SiteDescriptor{ Id=0109, IdRectangle = new Rectangle(149,69,157-149,77-69), Name="Nobile 3" } },

            // Northern
            { 1001, new SiteDescriptor{ Southern = false, Id=1001, IdRectangle = new Rectangle(122,150,128-122,154-150), Name="Peary_South" } },
            { 1002, new SiteDescriptor{ Southern = false, Id=1002, IdRectangle = new Rectangle(103,135,107-103,139-135), Name="Hermite A East" } },

        };

        public static StudyInterval NorthStudyInterval = new StudyInterval {
            Start = new DateTime(2023, 5, 1),
            Stop = new DateTime(2023, 9, 6),
            Step = TimeSpan.FromHours(2) };

        public static StudyInterval SouthStudyInterval = new StudyInterval
        {
            Start = new DateTime(2022, 12, 1),
            Stop = new DateTime(2023, 4, 1),
            Step = TimeSpan.FromHours(2)
        };

        public static IEnumerable<int> SiteKeys => SiteDictionary.Keys.Where(k => SiteDictionary[k].Southern == !InMemoryInt16Terrain.Singleton.IsNorth);
        public static IEnumerable<SiteDescriptor> Sites => SiteDictionary.Values.Where(s => s.Southern == !InMemoryInt16Terrain.Singleton.IsNorth);
        public static IEnumerable<int> SouthernSiteKeys => SiteDictionary.Keys.Where(k => SiteDictionary[k].Southern);
        public static IEnumerable<int> NorthernSiteKeys => SiteDictionary.Keys.Where(k => !SiteDictionary[k].Southern);
        public static IEnumerable<SiteDescriptor> SouthernSites => SiteDictionary.Values.Where(s => s.Southern);
        public static IEnumerable<SiteDescriptor> NorthernSites => SiteDictionary.Values.Where(s => s.Southern);

        public static float SiteRadius = 9f;  // km
        public static float SiteBoxRadius = 10f;  // km  400 pixels
        public const float MetersPerPixel = 20f;

        public static SiteDescriptor GetSiteDescriptor(int site) => SiteDictionary.TryGetValue(site, out SiteDescriptor s) ? s : null;

        public static Rectangle GetSiteBounds(int site)
        {
            if (!SiteDictionary.TryGetValue(site, out SiteDescriptor desc))
                return new Rectangle();
            return desc.Bounds;
        }

        public static SiteRectangle GetSiteRectangle(int site) => new SiteRectangle { Id = site.ToString(), Bounds = GetSiteBounds(site).ToRectangleF() };

        public static Point GetSiteCenter(int site)
        {
            var desc = SiteDictionary[site];  // will error if no site
            return InMemoryInt16Terrain.GetLineSampleDegrees(desc.Lat,desc.Lon);
        }

        public static Rectangle GetSiteInnerBounds(int site) => GetSiteBounds(site).Center().SpanFrom(new Point(80, 80));

        #region File Path Manipulation

        public static string SiteDatasetRoot(int site) => Path.Combine(SiteRoot, string.Format(@"site_{0:D3}", site));
        public static string SiteBaseFilename(int site) => string.Format(@"site_{0:D3}", site);

        #endregion

        #region Top Level

        public static void Init()
        {
            SiteRoot = Path.Combine(ViperEnvironment.MapRoot, SiteDirectory);
        }

        public void GenerateDatasets()
        {
            if (false)
            {
                foreach (var key in SiteKeys)
                    GenerateDatasetForStudySite(key);
            }
            else
            {
                GenerateDatasetForStudySite(1002);
            }
        }

        public void GenerateDatasetForStudySite(int site, bool regen = false)
        {
            Console.WriteLine($"** Writing dataset for site {site}");
            Init();
            var dataset_root = SiteDatasetRoot(site);
            if (!Directory.Exists(dataset_root))
                Directory.CreateDirectory(dataset_root);

            // Generate dem and slope-related maps
            if (true)
            {
                var r = GetSiteBounds(site);
                var base_path_9km = Path.Combine(SiteDatasetRoot(site), SiteBaseFilename(site)) + ".tif";
                GenerateDEMandSlope(site, base_path_9km, r, regen);
            }

            // Needed for PSR
            if (true)
            {
                GenerateAverageImages(site, height: 0, regen: regen, interval: StudyInterval.MazaricoCycle, sun: true);
            }

            // Generating other lighting-related maps
            if (false)  // normally false
            {
                //GenerateAverageImages(site, height: 0, regen: regen, interval: LightingInterval.MazaricoCycle, sun: true);
                //GenerateAverageImages(site, height: 2, regen: regen, interval: LightingInterval.MazaricoCycle, sun: true);
                GenerateAverageImages(site, height: 0, regen: regen, interval: StudyInterval.SiteStudyPeriod, sun: true);
                GenerateAverageImages(site, height: 2, regen: regen, interval: StudyInterval.SiteStudyPeriod, sun: true, earth: true, multipath_threshold: 2f);
                //GenerateAverageImages(site, height: 2, regen: regen, interval: LightingInterval.SiteStudyPeriod, sun: false, earth: true, multipath_threshold: 0f);

                //GenerateAverageImages(site, height: 0, regen: regen, interval: LightingInterval.MetonicCycle, sun: true);
                //GenerateAverageImages(site, height: 0, regen: regen, interval: LightingInterval.FourMetonicCycles, sun: true);
            }

            if (true)
            {
                GeneratePSRImage(site);
                CombineSunSlopePSR(site);
            }

            if (false)  // normally false
            {
                var gen = new TileLightingProductManager { };
                gen.WriteLongestNightsGeotiffs(site, interval: StudyInterval.SiteStudyPeriod, observer: 2);   
            }
            
            if (false)  // normally false
            {
                var gen = new TileLightingProductManager { };
                gen.WriteLongestNightsGeotiffs(site, interval: StudyInterval.ViperStudyPeriod1, observer: 1);
            }

            //if (false)  // normally false
            //    PlotHorizonAndLightCurveForBestPixel(site, interval: StudyInterval.SiteStudyPeriod, observer: 2);

            if (false)  // normally false
            {
                var gen = new TileLightingProductManager { };
                gen.GenerateHeightVsAverageSun(site);
            }

            if (false)  // normally false
            {
                var gen = new TileLightingProductManager { };
                gen.WritePixelCounts(site);
            }

            if (false)  // normally false
            {
                var gen = new TileLightingProductManager { };
                var interval = StudyInterval.SiteStudyPeriod;
                var observer = 10;
                gen.GenerateFullSunDurationMaps(site, interval.Start, interval.Stop, interval.Step, observer, 10f, 5);
            }

            if (true)
            {
                var base_path = Path.Combine(SiteDatasetRoot(site), SiteBaseFilename(site)) + ".tif";
                var path_20m = AppendToFilename(base_path, "_mosaic_20m", ".tif");
                GenerateHighResImage(GetSiteBounds(site), 20, path_20m);
                var path_1m = AppendToFilename(base_path, "_mosaic_1m", ".tif");
                GenerateHighResImage(GetSiteBounds(site), 1, path_1m);

                //using (var bmp = new Bitmap(path_1m))
                //    bmp.Save(AppendToFilename(path_1m, "", ".png"), ImageFormat.Png);
            }

            if (true)
            {
                GenerateIceStabilityImages(site);
            }

            if (true)
            {
                GenerateSafeHavenImages(site);
            }
        }

        #endregion

        #region DEM and slope images using GDAL

        // ignore regen
        private void GenerateDEMandSlope(int site, string base_path, Rectangle r, bool regen = false)
        {
            Debug.Assert(r.Width > 0);
            var dem = GetDEMArray(r);

            var dem_path = AppendToFilename(base_path, "_dem", ".tif");
            GeotiffHelper.WriteArrayAsGeotiff(dem, r, dem_path);

            //File.WriteAllText("color-slope.txt", "0 255 255 255\r\n90 0 0 0");

            var hillshade_path = AppendToFilename(base_path, "_hillshade", ".tif");
            RunGDALCommand($"gdaldem hillshade -of GTiff {dem_path} {hillshade_path}");
            WriteImageMetadata(hillshade_path, "Hillshade", "grayscale.map", "Default hillshade generated by GDAL");
            //FinishSiteImages(hillshade_png, "Hillshade");

            var slope_filename = AppendToFilename(base_path, "_slope", ".tif");
            RunGDALCommand($"gdaldem slope -of GTiff {dem_path} {slope_filename}");
            WriteImageMetadata(slope_filename, "Slope (deg)", "grayscale.map", "Slope in degrees generated by GDAL");
            var trafficability_path = AppendToFilename(base_path, "_trafficability", ".tif");
            RunGDALCommand($"gdaldem color-relief -of GTiff {slope_filename} gdallib/slope-color3.txt {trafficability_path}");
            WriteImageMetadata(trafficability_path, "Trafficability", "", "0-6 deg=green, 7-15 deg=yellow, >15 deg=red");
            //FinishSiteImages(slope_safety_filename, "Slope (safe landing/roving)");

            //var slope_colormap_filename = AppendToFilename(base_path, "_slope", ".png");
            //RunGDALCommand($"gdal_translate -ot Byte -scale 0 60 255 0 {slope_filename} {slope_colormap_filename}");
            //(int, int) minmax;
            //using (var bmp = Image.FromFile(slope_colormap_filename) as Bitmap)
            //    minmax = bmp.GetMinMax();
            //FinishSiteImages(slope_colormap_filename, "Slope (deg)", new ViridisColorTable { Low = 0, High = minmax.Item2 });

            var contour_shapefile = AppendToFilename(base_path, "_contour", ".shp");
            RunGDALCommand($"gdal_contour -a height {dem_path} {contour_shapefile} -i 100");
            var contour_overlay = AppendToFilename(base_path, "_contour_overlay", ".tif");
            File.Delete(contour_overlay);
            //RunGDALCommand($"gdal_translate -of GTiff -ot Byte {hillshade_filename} {contour_filename}");
            RunGDALCommand($"gdal_rasterize -of GTiff -ts {r.Width} {r.Height} -co alpha=yes -a height -add -init 0 0 0 0 -burn 255 0 0 255 -ot Byte {contour_shapefile} {contour_overlay}");

            var contour_path = AppendToFilename(base_path, "_contour", ".tif");
            var tmpfile = Path.ChangeExtension(contour_path, ".png");
            using (var overlay_bmp = Image.FromFile(contour_overlay))
            using (var hillshade_bmp = Image.FromFile(hillshade_path))
            using (var bmp = new Bitmap(r.Width, r.Height, PixelFormat.Format32bppArgb))
            using (var g = Graphics.FromImage(bmp))
            {
                g.DrawImageUnscaled(hillshade_bmp, 0, 0);
                g.DrawImageUnscaled(overlay_bmp, 0, 0);
                bmp.Save(tmpfile, ImageFormat.Png);
            }
            using (var bmp = Image.FromFile(tmpfile) as Bitmap)
                GeotiffHelper.WriteBitmapAsGeotiff(bmp, GetSiteBounds(site), contour_path);
            WriteImageMetadata(contour_path, "Topographic Map (100m contour interval)");

            File.Delete(tmpfile);
            File.Delete(contour_shapefile);
            File.Delete(contour_overlay);
            File.Delete(Path.ChangeExtension(contour_shapefile, ".shx"));
            File.Delete(Path.ChangeExtension(contour_shapefile, ".prj"));
            File.Delete(Path.ChangeExtension(contour_shapefile, ".dbf"));
            //FinishSiteImages(contour_filename, "Contour (100m steps)");
        }

        #endregion

        #region Lighting images

        public void GenerateAverageImages(int site, StudyInterval interval = null, int height = 0, bool regen = false, bool sun = false, bool earth = false, float multipath_threshold = 2f)
        {
            if (interval == null) interval = StudyInterval.SiteStudyPeriod;
            if (!SiteDictionary.TryGetValue(site, out SiteDescriptor desc))
                return;
            var rect = desc.Bounds;
            var root = SiteDatasetRoot(site);

            var (start, stop, step) = interval.GetInterval();

            var average_sun_path = Path.Combine(root, TileLightingProductManager.AverageSunFilename(site, start, stop, step, height));
            var average_earth_path = Path.Combine(root, TileLightingProductManager.AverageEarthFilename(site, start, stop, step, height, multipath_threshold));
            var average_sun_exists = !sun || File.Exists(average_sun_path);  // disable generation by pretending the file already exists
            var average_earth_exists = !earth || File.Exists(average_earth_path);

            var gen = new TileLightingProductManager
            {
                //MainWindow = MainWindow,
                Selection = rect,
                IntervalStart = start,
                IntervalStop = stop,
                IntervalStep = step,
                ObserverHeightInMeters = height,
                EarthMultipathThreshold = multipath_threshold,
                GenerateAverageSun = !average_sun_exists,
                GenerateAverageEarth = !average_earth_exists,
                SunExposurePath = average_sun_path,
                EarthExposurePath = average_earth_path
            };

            if (!average_sun_exists || !average_earth_exists)
                gen.GenerateAverageSunAndEarthImages();

            if (File.Exists(average_sun_path))
                WriteLightingGeotiffForRick(
                    bounds: rect,
                    path: average_sun_path,
                    title: "Average Illumination",
                    details: $"Observer={height}m start={start:d} stop={stop:d} step={step.TotalHours} hr",
                    colormap: "lighting_colormap_float.map");
            if (File.Exists(average_earth_path))
                WriteLightingGeotiffForRick(
                    bounds: rect,
                    path: average_earth_path,
                    title: "Average Earth Visibility (earth center)",
                    details: $"Observer={height}m multipath_mask={multipath_threshold:F0.00} start={start:d} stop={stop:d} step={step.TotalHours} hr",
                    colormap: "lighting_colormap_float.map");

            if (File.Exists(average_sun_path))
            {
                using (var writer = new StreamWriter(AppendToFilename(average_sun_path, "_histogram", ".txt")))
                using (var src = Image.FromFile(average_sun_path) as Bitmap)
                {
                    var buf = src.GetByteHistogram();
                    for (var i = 0; i < buf.Length; i++)
                        writer.WriteLine($"{100 * i / 255f}\t{buf[i]}");
                }
            }
        }

        void WriteLightingGeotiffForRick(Rectangle bounds, string path = null, string title = "empty title", string colormap = "lighting_colormap_float.map", string details = null)
        {
            float[,] ary;
            using (var bmp = Image.FromFile(path) as Bitmap)
                ary = bmp.ToByteArray2().Mult(1 / 255f);
            GeotiffHelper.WriteArrayAsGeotiff(ary, bounds, path.FilenameAppend("_geotiff", ".tif"));
            WriteImageMetadata(path, title, colormap, details);
        }

        public static void WriteImageMetadata(string image_path, string title = "empty title", string colormap = "lighting_colormap_float.map", string details = null)
        {
            using (var sw = new StreamWriter(image_path.FilenameAppend("_metadata", ".txt")))
            {
                sw.WriteLine(title);
                sw.WriteLine(colormap);
                if (details != null)
                    sw.WriteLine(details);
            }
        }

        void GeneratePSRImage(int site)
        {
            var src_path = Path.Combine(SiteDatasetRoot(site), SiteBaseFilename(site) + "_avg_sun_start19700101_stop20440101_step006_observer000.png");
            var dst_path = Path.Combine(SiteDatasetRoot(site), SiteBaseFilename(site) + "_psr.png");
            using (var src = Image.FromFile(src_path) as Bitmap)
            using (var dst = BitmapHelper.ToFormat8bppIndexed(src, v => (byte)(v == 0 ? 255 : 0)))
            {
                var p = dst.Palette;
                for (var i = 0; i < 255; i++) p.Entries[i] = Color.FromArgb(0, 0, 0, 0);
                p.Entries[255] = Color.Red;
                dst.Palette = p;
                dst.Save(dst_path, ImageFormat.Png);
            }
        }

        void CombineSunSlopePSR(int site)
        {
            var psr_path = Path.Combine(SiteDatasetRoot(site), SiteBaseFilename(site) + "_psr.png");
            var sun_path = Path.Combine(SiteDatasetRoot(site), SiteBaseFilename(site) + @"_avg_sun_start19700101_stop20440101_step006_observer000.png");
            var slp_path = @"slope_tmp.png";
            var target_path = Path.Combine(SiteDatasetRoot(site), SiteBaseFilename(site) + "_sunSlopePSR.png");

            var slope_path = Path.Combine(SiteDatasetRoot(site), SiteBaseFilename(site) + "_slope.tif");
            RunGDALCommand($"gdaldem color-relief -of PNG {slope_path} gdallib/slope-color4.txt {slp_path}");

            using (var sun = Image.FromFile(sun_path) as Bitmap)
            using (var slp = Image.FromFile(slp_path) as Bitmap)
            using (var psr = Image.FromFile(psr_path) as Bitmap)
            using (var target = new Bitmap(sun.Width, sun.Height, PixelFormat.Format32bppArgb))
            using (var g = Graphics.FromImage(target))
            {
                Debug.Assert(sun.Width == slp.Width && slp.Width == psr.Width && sun.Height == slp.Height && slp.Height == psr.Height);

                var p = sun.Palette;
                for (var i = 0; i < 256; i++) p.Entries[i] = Color.FromArgb(255 , i, i, i);
                sun.Palette = p;

                var dstRect = new Rectangle(0, 0, sun.Width, sun.Height);
                var _imageAttributes = new ImageAttributes();  // Holds the image attributes so it can be reused
                var _colorMatrix = new ColorMatrix();
                g.CompositingMode = System.Drawing.Drawing2D.CompositingMode.SourceOver;

                g.DrawImageUnscaled(slp, 0, 0);

                _colorMatrix.Matrix33 = 0.8f;
                _imageAttributes.SetColorMatrix(_colorMatrix);
                g.DrawImage(sun, dstRect, 0, 0, sun.Width, sun.Height, GraphicsUnit.Pixel, _imageAttributes);

                _colorMatrix.Matrix33 = 0.4f;
                _imageAttributes.SetColorMatrix(_colorMatrix);
                g.DrawImage(psr, dstRect, 0, 0, psr.Width, psr.Height, GraphicsUnit.Pixel, _imageAttributes);

                target.Save(target_path, ImageFormat.Png);
            }
        }

        #endregion

        #region GeoTiff Support

        #endregion

        #region Image decorations

        /// <summary>
        /// Draw scale circles
        /// </summary>
        /// <param name="path"></param>
        /// <param name="r1">radius in km</param>
        /// <param name="r2">radius in km</param>
        /// <param name="penWidth"></param>
        /// <param name="cmap"></param>
        public static void FinishDatasetPng(string path, float r1 = -1f, float r2 = -1f, Color? circleColor = null, float penWidth = 1, ListColorTable cmap = null)
        {
            Debug.Assert(".png".Equals(Path.GetExtension(path)));
            if (!circleColor.HasValue) circleColor = Color.White;
            var src_bmp = Image.FromFile(path) as Bitmap;
            if (cmap != null)
            {
                var p = src_bmp.Palette;
                cmap.LoadPalette(p);
                src_bmp.Palette = p;
            }
            var bmp32 = src_bmp.ToFormat32bppArgb();
            if (bmp32 != src_bmp)
            {
                using (var g = Graphics.FromImage(bmp32))
                    g.DrawImageUnscaled(src_bmp, 0, 0);
                src_bmp.Dispose();
            }

            var center = new Point(bmp32.Width / 2, bmp32.Height / 2);
            r1 *= 1000f / 20f;
            r2 *= 1000f / 20f;

            using (var g = Graphics.FromImage(bmp32))
            using (var brush = new SolidBrush(circleColor.Value))
            using (var pen = new Pen(brush, penWidth))
            {
                g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;
                if (r1 > 0f)
                    g.DrawEllipse(pen, center.X - r1, center.Y - r1, 2 * r1, 2 * r1);
                if (r2 > 0f)
                    g.DrawEllipse(pen, center.X - r2, center.Y - r2, 2 * r2, 2 * r2);
            }
            bmp32.Save("temp.png", ImageFormat.Png);
            bmp32.Dispose();
            File.Copy("temp.png", path, true);
        }

        void SaveColormapLegend(string path, Colormap map, int width = 60, int color_width = 30, int height = 512, int steps = 11, int y_border = 8, Func<Colormap, int, List<(float, string)>> func = null)
        {
            using (var bmp = map.MakeLegend(width, color_width, height, steps, y_border, func))
                bmp.Save(path, ImageFormat.Png);
        }

        public static void FinishSiteImages(string source_path, string title, Colormap map = null)
        {
            var radius1 = 9 * 1000 / 20;
            var radius2 = 1 * 1000 / 20;
            using (var src = Image.FromFile(source_path) as Bitmap)
            {
                if (map != null)
                {
                    var p = src.Palette;
                    map.LoadPalette(p);
                    src.Palette = p;
                }
                if (false)
                using (var src32bit = src.ToFormat32bppArgb())
                {
                    var center = new Point(src32bit.Width / 2, src32bit.Height / 2);

                    using (var g = Graphics.FromImage(src32bit))
                    //using (var brush = new SolidBrush(Color.FromArgb(200, 255, 100, 100)))
                    using (var pen = new Pen(Brushes.White, 2))
                    {
                        g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;
                        g.DrawEllipse(pen, center.X - radius1, center.Y - radius1, 2 * radius1, 2 * radius1);
                        g.DrawEllipse(pen, center.X - radius2, center.Y - radius2, 2 * radius2, 2 * radius2);
                    }
                    if (map == null)
                        src32bit.Save(AppendToFilename(source_path, "_9km", ".png"), ImageFormat.Png);
                    else
                        using (var src_w_legend = map.AddLegend(src32bit, title: title, steps: 11))
                            src_w_legend.Save(AppendToFilename(source_path, "_9km_legend", ".png"), ImageFormat.Png);
                }
            }

            using (var src = Image.FromFile(source_path) as Bitmap)
            {
                var smaller_span = (int)(1.4f * radius2);
                using (var smaller = new Bitmap(smaller_span * 2, smaller_span * 2, PixelFormat.Format32bppArgb))
                using (var g = Graphics.FromImage(smaller))
                using (var pen = new Pen(Brushes.White, 1))
                {
                    g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;
                    var center = new Point(smaller.Width / 2, smaller.Height / 2);
                    g.DrawImageUnscaled(src, smaller_span - src.Width / 2, smaller_span - src.Height / 2);
                    g.DrawEllipse(pen, center.X - radius2, center.Y - radius2, 2 * radius2, 2 * radius2);

                    //if (map == null)
                        smaller.Save(AppendToFilename(source_path, "_1km", ".png"), ImageFormat.Png);
                    //else
                    //    using (var src_w_legend = map.AddLegend(smaller, title: title, steps: 11))
                    //        src_w_legend.Save(AppendToFilename(source_path, "_1km_legend", ".png"), ImageFormat.Png);
                }
            }
        }

        #endregion

        #region Shell Commands

        public static void RunGDALCommand(string cmd)
        {
            var process = new Process();
            process.StartInfo.FileName = "cmd.exe";
            process.StartInfo.RedirectStandardInput = true;
            process.StartInfo.RedirectStandardOutput = true;
            process.StartInfo.RedirectStandardError = true;
            process.StartInfo.CreateNoWindow = true;
            process.StartInfo.UseShellExecute = false;
            process.Start();

            process.StandardInput.WriteLine("c:/Apps/gdal/SDKShell.bat");
            process.StandardInput.WriteLine(cmd);
            process.StandardInput.Flush();
            process.StandardInput.Close();
            var std_msgs = process.StandardOutput.ReadToEnd();

            process.WaitForExit();
            Console.WriteLine(std_msgs);
            var err_msgs = process.StandardError.ReadToEnd();
            Console.WriteLine($"stderr: {err_msgs}");
        }

        public static string AppendToFilename(string path, string postfix, string extension = null) => Path.Combine(Path.GetDirectoryName(path) ?? ".", Path.GetFileNameWithoutExtension(path) + postfix + (extension ?? Path.GetExtension(path)));

        #endregion

        #region Utilities

        float[,] NewFloatArray(Size s) => new float[s.Height, s.Width];

        float[,] GetDEMArray(Rectangle r)
        {
            var dem = NewFloatArray(r.Size);
            var terrain = ViperEnvironment.Terrain;
            for (var line = 0; line < r.Height; line++)
                for (var sample = 0; sample < r.Width; sample++)
                    dem[line, sample] = terrain.LineSampleToTerrainOffset(line + r.Top, sample + r.Left);
            return dem;
        }

        #endregion

        #region Horizon Plots
        /*

        void PlotHorizonAndLightCurveForBestPixel(int site, StudyInterval interval = null, int observer = 2)
        {
            //var center = GetSiteCenter(site);
            var center = GetPointWithMostAverageSun(site);

            var horizon_path = HorizonPlotPath(site);
            PlotHorizon(horizon_path, center.Y, center.X, interval, observer);

            var light_curve_path = LightCurvePlotPath(site);
            PlotLightCurve(light_curve_path, center.Y, center.X, interval, observer);
        }

        string HorizonPlotPath(int site) => Path.Combine(SiteDatasetRoot(site), SiteBaseFilename(site) + "_horizon.png");
        string LightCurvePlotPath(int site) => Path.Combine(SiteDatasetRoot(site), SiteBaseFilename(site) + "_light_curve.png");

        public static Point GetPointWithMostAverageSun(int site)
        {
            var rect = GetSiteBounds(site);
            var path = Path.Combine(SiteDatasetRoot(site), SiteBaseFilename(site) + "_avg_sun_start19700101_stop20440101_step006_observer000.png");
            byte[,] ary = null;
            using (var bmp = Image.FromFile(path) as Bitmap)
                ary = bmp.ToByteArray2();
            var point = new Point();
            var max = -1;
            for (var row = 0; row < ary.GetLength(0); row++)
                for (var col = 0; col < ary.GetLength(1); col++)
                {
                    var v = ary[row, col];
                    if (v > max)
                    {
                        max = v;
                        point = new Point(col, row);
                    }
                }
            return new Point(point.X + rect.X, point.Y + rect.Y);
        }

        void PlotHorizon(string path, int line, int sample, StudyInterval interval, int observer) => PlotHorizon(path,line, sample, interval.Start, interval.Stop, interval.Step, observer);
        void PlotHorizon(string path, int line, int sample, DateTime? start1 = null, DateTime? stop1 = null, TimeSpan? step1 = null, int observer = 0)
        {
            var start = start1 ?? StudyInterval.SiteStudyPeriod.Start;
            var stop = stop1 ?? StudyInterval.SiteStudyPeriod.Stop;
            var step = step1 ?? StudyInterval.SiteStudyPeriod.Step;
            var id = TerrainPatch.LineSampleToId(line, sample);
            var patch = ViperEnvironment.GetPatch(id, observer);
            var point = patch.PointInPatch(new Point(sample, line));
            var horizon = patch.Horizons[point.Y][point.X];

            var sunvecs = TileLightingProductManager.GetSunVectors(start, stop, step);
            var earthvecs = TileLightingProductManager.GetEarthVectors(start, stop, step);

            InMemoryInt16Terrain.GetLatLonDegrees(line, sample, out double lat, out double lon);
            Func<double, double> map_azimuth = azimuth_deg =>
             {
                 azimuth_deg -= lon;
                 if (azimuth_deg < 0d) azimuth_deg += 360d;                   // Put in [0,360)
                 if (azimuth_deg >= 360d) azimuth_deg -= 360d;
                 azimuth_deg -= 180d;
                 if (azimuth_deg < 0d) azimuth_deg += 360d;                   // Put in [0,360)
                 if (azimuth_deg >= 360d) azimuth_deg -= 360d;                 
                 return azimuth_deg - 180d;                                   // move to [-180,180)
             };

            var horizon_points = GetHorizonCurve(horizon, Color.Black, map_azimuth);

            patch.FillMatrices(ViperEnvironment.Terrain);
            var sun_points = TransformVectors(sunvecs, patch.Matrices[point.Y][point.X], Color.DarkGoldenrod, map_azimuth);
            var earth_points = TransformVectors(earthvecs, patch.Matrices[point.Y][point.X], Color.Blue, map_azimuth);

            using (var zg = new ZedGraphControl())
            using (var g = zg.CreateGraphics())
            {
                zg.MasterPane.ReSize(g, new RectangleF(0, 0, 1000, 600));
                AddCurves(zg.GraphPane, "Horizon", horizon_points, Color.Black, SymbolType.None);
                AddCurves(zg.GraphPane, "Sun", sun_points, Color.Gold, SymbolType.None);
                AddCurves(zg.GraphPane, "Earth", earth_points, Color.Blue, SymbolType.None);

                zg.GraphPane.XAxis.Scale.Min = -180;
                zg.GraphPane.XAxis.Scale.Max = 180;
                zg.GraphPane.YAxis.Scale.Min = -10;
                zg.GraphPane.YAxis.Scale.Max = 10;

                zg.GraphPane.Title.Text = $"line={line} sample={sample}";
                zg.GraphPane.Title.FontSpec = new FontSpec { Size = 8f };

                zg.GraphPane.XAxis.Title.Text = "Azimuth (deg)";
                zg.GraphPane.YAxis.Title.Text = "Elevation (deg)";

                zg.GraphPane.AxisChange();
                if (File.Exists(path)) File.Delete(path);
                using (var image = zg.MasterPane.GetImage())
                    image.Save(path, ImageFormat.Png);
                zg.GraphPane.CurveList.Clear();
            }
        }

        ZedGraphControl GetGraphPane() => _zedGraphPane ?? (_zedGraphPane = new ZedGraphControl());

        List<PointPairList> GetHorizonCurve(Horizon h, Color color, Func<double, double> map_azimuth)
        {
            var r = new List<PointPairList>();
            var points = new PointPairList();
            foreach (var pt in h.PartialPoints().Where(pt => pt.HasValue))
            {
                var azimuth_deg = map_azimuth(pt.Value.X);
                var elevation_deg = pt.Value.Y;
                if (points.Count > 0 && Math.Abs(azimuth_deg - points[points.Count - 1].X) > 20d)
                {
                    r.Add(points);
                    points = new PointPairList();
                }
                points.Add(azimuth_deg, elevation_deg);
            }
            if (points.Count > 0)
                r.Add(points);
            return r;
        }

        List<PointPairList> TransformVectors(List<Vector3d> vecs, Matrix4d mat, Color color, Func<double,double> map_azimuth)
        {
            var r = new List<PointPairList>();
            var points = new PointPairList();
            Vector3d temp = new Vector3d();
            for (var i = 0; i < vecs.Count; i++)
            {
                var vec = vecs[i];
                TerrainPatch.Transform(ref vec, ref mat, ref temp);
                var rx = temp[0];
                var ry = temp[1];
                var rz = temp[2];
                var azimuth_rad = Math.Atan2(ry, rx) + Math.PI;  // [0,2PI]
                var alen = Math.Sqrt(rx * rx + ry * ry);
                var elevation_rad = Math.Atan2(rz, alen);

                var azimuth_deg = azimuth_rad * 180 / Math.PI;
                azimuth_deg = map_azimuth(azimuth_deg);

                var elevation_deg = elevation_rad * 180 / Math.PI;

                if (points.Count>0 && Math.Abs(azimuth_deg - points[points.Count-1].X) > 20d)
                {
                    r.Add(points);
                    points = new PointPairList();
                }
                points.Add(azimuth_deg, elevation_deg);
            }
            if (points.Count > 0)
                r.Add(points);
            return r;
        }

        void AddCurves(GraphPane pane, string title, List<PointPairList> pts, Color color, SymbolType symbol)
        {
            for (var i = 0; i < pts.Count; i++)
                pane.AddCurve(i == 0 ? title : null, pts[i], color, symbol);
        }

        ZedGraphControl _zedGraphPane = null;

        void PlotLightCurve(string path, int line, int sample, StudyInterval interval, int observer) => PlotLightCurve(path, line, sample, interval.Start, interval.Stop, interval.Step, observer);
        void PlotLightCurve(string path, int line, int sample, DateTime? start1 = null, DateTime? stop1 = null, TimeSpan? step1 = null, int observer = 0)
        {
            var start = start1 ?? StudyInterval.SiteStudyPeriod.Start;
            var stop = stop1 ?? StudyInterval.SiteStudyPeriod.Stop;
            var step = step1 ?? StudyInterval.SiteStudyPeriod.Step;

            var gen = new LightCurveApprox();
            var points = gen.GetLightCurve(line, sample, start, stop, step);

            using (var zg = new ZedGraphControl())
            using (var g = zg.CreateGraphics())
            {
                zg.GraphPane.CurveList.Clear();
                zg.MasterPane.ReSize(g, new RectangleF(0, 0, 1000, 400));
                AddCurves(zg.GraphPane, "Sun", new List<PointPairList>() { points }, Color.DarkGoldenrod, SymbolType.None);

                zg.GraphPane.XAxis.Title.Text = "Time";
                zg.GraphPane.XAxis.Type = AxisType.DateAsOrdinal;
                zg.GraphPane.YAxis.Title.Text = "Sun Fraction";
                zg.GraphPane.YAxis.Type = AxisType.Linear;
                zg.GraphPane.Title.Text = "Sun fraction over time";

                zg.GraphPane.AxisChange();
                if (File.Exists(path)) File.Delete(path);
                using (var image = zg.MasterPane.GetImage())
                    image.Save(path, ImageFormat.Png);
                zg.GraphPane.CurveList.Clear();
            }
        }
        */
        #endregion

        #region HiRes Image

        //void GenerateHighResImage(int site) => GenerateHighResImage(GetSiteBounds(site));
        public static void GenerateHighResImage(Rectangle r, int metersPerPixel, string dst_path, bool compress = true)
        {
            //gdalwarp -overwrite -te -21290.000 -22310.000 -1290.000 -2310.000 -tr 1 1  LRO_NAC_AvgMosaic_SPole855_1mp.tif site001_1m.tif

            var src_path = AppConfiguration.Get("PolarMosaic", @"C:\RP\maps\LRO_NAC_AvgMosaic_SPole855_1mp.tif");

            //var xmin = -21290.000;
            //var xmax = -1290.000;
            //var ymin = -22310.000;
            //var ymax = -2310.000;

            var width = r.Width * (20d / metersPerPixel);
            var height = r.Height * (20d / metersPerPixel);
            var t = GetAffineTransformOfBounds(r);
            var xmin = t[0];
            var xmax = xmin + width * metersPerPixel;
            var ymax = t[3];
            var ymin = ymax - height * metersPerPixel;

            var of = GetOutputFormat(dst_path);
            var compression = compress ? "-co COMPRESS=JPEG -co TILED=YES" : "";

            var cmd = $"gdalwarp -overwrite -te {xmin} {ymin} {xmax} {ymax} -tr {metersPerPixel} {metersPerPixel} -of {of} {compression} {src_path} {dst_path}";
            Console.WriteLine(cmd);
            RunGDALCommand(cmd);

            string GetOutputFormat(string path)
            {
                switch (Path.GetExtension(path).ToLowerInvariant())
                {
                    case ".png": return "PNG";
                    default: return "GTiff";

                }
            }
        }

        private void GenerateHighResImage2(int site)
        {
            GdalConfiguration.ConfigureGdal();

            var bounds = GetSiteBounds(site);
            const string src_path = @"c:\Users\mshirley\Downloads\LRO_NAC_AvgMosaic_SPole855_1mp.tif";

            //var src_wkt = @"PROJCRS['PolarStereographic MOON',BASEGEOGCRS['GCS_MOON',DATUM['D_MOON',ELLIPSOID['MOON',1737400,0,LENGTHUNIT['metre',1,ID['EPSG',9001]]]],PRIMEM['Reference_Meridian',0,ANGLEUNIT['degree',0.0174532925199433,ID['EPSG',9122]]]],CONVERSION['unnamed',METHOD['Polar Stereographic (variant A)',ID['EPSG',9810]],PARAMETER['Latitude of natural origin',-90,ANGLEUNIT['degree',0.0174532925199433],ID['EPSG',8801]],PARAMETER['Longitude of natural origin',0,ANGLEUNIT['degree',0.0174532925199433],ID['EPSG',8802]],PARAMETER['Scale factor at natural origin',0,SCALEUNIT['unity',1],ID['EPSG',8805]],PARAMETER['False easting',0,LENGTHUNIT['metre',1],ID['EPSG',8806]],PARAMETER['False northing',0,LENGTHUNIT['metre',1],ID['EPSG',8807]]],CS[Cartesian,2],AXIS['easting',north,ORDER[1],LENGTHUNIT['metre',1]],AXIS['northing',north,ORDER[2],LENGTHUNIT['metre',1]]]";
            var src_wkt = @"PROJCS['unnamed',GEOGCS['GCS_MOON',DATUM['D_MOON',ELLIPSOID['MOON',1737400,0,LENGTHUNIT['metre',1,ID['EPSG',9001]]]],PRIMEM['Reference_Meridian',0,ANGLEUNIT['degree',0.0174532925199433,ID['EPSG',9122]]]],CONVERSION['unnamed',METHOD['Polar Stereographic (variant A)',ID['EPSG',9810]],PARAMETER['Latitude of natural origin',-90,ANGLEUNIT['degree',0.0174532925199433],ID['EPSG',8801]],PARAMETER['Longitude of natural origin',0,ANGLEUNIT['degree',0.0174532925199433],ID['EPSG',8802]],PARAMETER['Scale factor at natural origin',0,SCALEUNIT['unity',1],ID['EPSG',8805]],PARAMETER['False easting',0,LENGTHUNIT['metre',1],ID['EPSG',8806]],PARAMETER['False northing',0,LENGTHUNIT['metre',1],ID['EPSG',8807]]],CS[Cartesian,2],AXIS['easting',north,ORDER[1],LENGTHUNIT['metre',1]],AXIS['northing',north,ORDER[2],LENGTHUNIT['metre',1]]]";
            var dst_wkt = @"PROJCS['unnamed',GEOGCS['unnamed ellipse',DATUM['unknown',SPHEROID['unnamed',1737400,0]],PRIMEM['Greenwich',0],UNIT['degree',0.0174532925199433]],PROJECTION['Polar_Stereographic'],PARAMETER['latitude_of_origin',-90],PARAMETER['central_meridian',0],PARAMETER['scale_factor',1],PARAMETER['false_easting',0],PARAMETER['false_northing',0],UNIT['metre',1,AUTHORITY['EPSG','9001']]]";

            src_wkt = dst_wkt;  // for now

            var src_sr = new OSGeo.OSR.SpatialReference(ToWkt(src_wkt));
            var dst_sr = new OSGeo.OSR.SpatialReference(ToWkt(dst_wkt));
            var tx = new OSGeo.OSR.CoordinateTransformation(src_sr, dst_sr);

            using (Dataset src_ds = Gdal.Open(src_path, Access.GA_ReadOnly))
            {
                if (src_ds == null)
                    throw new Exception($"GDAL couldn't open {src_path}");

                var geo_t = new double[6];
                src_ds.GetGeoTransform(geo_t);
                var x_size = src_ds.RasterXSize;
                var y_size = src_ds.RasterYSize;
                double[] ul = new double[3], lr = new double[3];
                tx.TransformPoint(ul, geo_t[0], geo_t[3], 0d);
                tx.TransformPoint(lr, geo_t[0] + geo_t[1] * x_size, geo_t[3] + geo_t[5] * y_size, 0d);

                var driver = GeotiffHelper.GetGeoTIFFDriver();
                var dst_path=@"hack.tif";

                var pixel_spacing = 20;

                var width = (int)((lr[0]-ul[0])/ pixel_spacing);
                var height = (int)((ul[1] - lr[1]) / pixel_spacing);
                if (File.Exists(dst_path)) File.Delete(dst_path);
                using (var dst_ds = driver.Create(dst_path, bounds.Width * 20, bounds.Height * 20, 1, DataType.GDT_Byte, null))
                {
                    var transform = GetAffineTransformOfBounds(bounds);
                    var transform1 = GetGeoTransformFromFile( @"C:\Users\mshirley\Downloads\site001.tif");

                    var new_geo = new double[] { ul[0], pixel_spacing, geo_t[2], ul[1], geo_t[4], -pixel_spacing };
                    dst_ds.SetGeoTransform(transform1);  // compare with transform

                    dst_ds.SetProjection(ToWkt(dst_wkt));

                    Gdal.ReprojectImage(src_ds, dst_ds, ToWkt(src_wkt), ToWkt(dst_wkt), ResampleAlg.GRA_Bilinear, 0d, 0d, null, null, null);
                }
            }

            double[] GetGeoTransformFromFile(string path1)
            {
                var geo_t = new double[6];
                using (Dataset src_ds = Gdal.Open(path1, Access.GA_ReadOnly))
                {
                    if (src_ds == null)
                        throw new Exception($"GDAL couldn't open {path1}");

                    src_ds.GetGeoTransform(geo_t);
                    var x_size = src_ds.RasterXSize;
                    var y_size = src_ds.RasterYSize;
                    var proj1 = src_ds.GetProjection();
                }
                return geo_t;
            }
        }

        public static double[] GetAffineTransformOfBounds(Rectangle bounds)
        {
            var r = GetGeoExtent(bounds);
            return new double[] { r.Left, 20d, 0d, r.Top, 0d, -20d };
        }

        public static RectangleD GetGeoExtent(Rectangle r)
        {
            var ul = PixelToCoordinate(r.Top, r.Left);  // line,sample
            var lr = PixelToCoordinate(r.Bottom, r.Right);
            return new RectangleD(ul, lr);
        }

        public static utilities.PointD PixelToCoordinate(int line, int sample)
        {
            const double S0 = 15199.5d;             // PDS SAMPLE_PROJECTION_OFFSET
            const double L0 = 15199.5d;             // PDS LINE_PROJECTION_OFFSET
            const double Scale = 20d;

            var x = (sample - S0) * Scale;
            var y = (L0 - line) * Scale;
            return new utilities.PointD(x, y);
        }

        public static string ToWkt(string s) => s.Replace('\'', '"');

        #endregion

        public static List<TerrainPatch> GetPatchShells(Rectangle region, int observer_height, bool generate_horizons_if_needed)
        {
            if (observer_height < 0) observer_height = 0;
            var result = new List<TerrainPatch>();
            var queue = new List<TerrainPatch>();
            foreach (var id in TerrainPatch.EnumerateIds(TerrainPatch.CoveringIdRectangle(region)))
            {
                var patch = TerrainPatch.FromId(id, observer_height);
                (File.Exists(patch.Path) ? result : queue).Add(patch);
            }
            if (queue.Count > 0 && generate_horizons_if_needed)
            {
                Console.WriteLine($"Generating horizons for observer {observer_height}");
                var gpuProcessor = ViperEnvironment.Processor as CPUHorizons;
                gpuProcessor.RunQueue(queue, unloadHorizons: true);
                result.AddRange(queue);
            }
            return result;
        }

        IceStabilityTileGenerator ice_stability_tile_generator = null;
        void GenerateIceStabilityImages(int site)
        {
            var bounds = GetSiteBounds(site);
            if (ice_stability_tile_generator==null)
            {
                ice_stability_tile_generator = new IceStabilityTileGenerator { ChunkSize = 64 };
                ice_stability_tile_generator.Load();
            }

            var ice_depth_array = NewFloatArray(bounds.Size);
            ice_stability_tile_generator.FillArray(ice_depth_array, bounds);
            var base_path = Path.Combine(SiteDatasetRoot(site), SiteBaseFilename(site)) + ".tif";
            GeotiffHelper.WriteArrayAsGeotiff(ice_depth_array, bounds, base_path.FilenameAppend("_ice_stab_depth_geotiff", ".tif"));
            using (var bmp = BitmapHelper.ToFormat8bppIndexed(ice_depth_array, f => f <= 0f ? (byte)0 : (byte)(f * 100)))
            {
                var jet = new InvertedJetColorTable();
                var p = bmp.Palette;
                jet.LoadPalette(p, 255, 0, 200);
                p.Entries[0] = Color.White;
                for (var i = 200; i < 256; i++)
                    p.Entries[i] = Color.FromArgb(216, 216, 219);
                bmp.Palette = p;
                bmp.Save(base_path.FilenameAppend("_ice_stab_depth", ".png"), ImageFormat.Png);
                for (var i = 1; i < 50; i++) p.Entries[i] = Color.Green;
                for (var i = 50; i < 100; i++) p.Entries[i] = Color.Yellow;
                for (var i = 100; i < 256; i++) p.Entries[i] = Color.DarkGray;
                bmp.Palette = p;
                bmp.Save(base_path.FilenameAppend("_ice_stab_zone", ".png"), ImageFormat.Png);
            }
        }

        SafeHavenGenerator _safeHavenGenerator;
        void GenerateSafeHavenImages(int site)
        {
            var desc = GetSiteDescriptor(site);
            var interval = desc.StudyInterval;
            var path = Path.Combine(SiteDatasetRoot(site), SiteBaseFilename(site) + "_safe_havens") + ".tif";
            if (_safeHavenGenerator == null)
                _safeHavenGenerator = new SafeHavenGenerator
                {
                    //MainWindow = LunarHorizon.Singleton,
                    Region = GetSiteBounds(site),
                    EarthMultipathThreshold = 2f,
                    SafeHavenSunThreshold = 0.5f,
                    ObserverHeightInMeters = 0,
                    Start = interval.Start,
                    Stop = interval.Stop,
                    Step = interval.Step
                };
            _safeHavenGenerator.Region = GetSiteBounds(site);
            _safeHavenGenerator.WriteSafeHavenGeotiffs(path);
        }
    }

    public class StudyInterval
    {
        public DateTime Start { get; set; } = new DateTime(1970, 1, 1);
        public DateTime Stop { get; set; } = new DateTime(2044, 1, 1);
        public TimeSpan Step { get; set; } = TimeSpan.FromHours(6);

        public StudyInterval() { }
        public StudyInterval(DateTime start,DateTime stop,TimeSpan step) { Start = start;Stop = stop;Step = step; }

        public (DateTime, DateTime, TimeSpan) GetInterval() => (Start, Stop, Step);

        public static StudyInterval MetonicCycle = new StudyInterval { Start = new DateTime(2020, 1, 1), Stop = new DateTime(2020, 1, 1).Add(TimeSpan.FromDays(365.25 * 18.6)), Step = TimeSpan.FromHours(2) };
        public static StudyInterval FourMetonicCycles = new StudyInterval { Start = new DateTime(1970, 1, 1), Stop = new DateTime(2044, 1, 1), Step = TimeSpan.FromHours(2) };
        public static StudyInterval MazaricoCycle = new StudyInterval { Start = new DateTime(1970, 1, 1), Stop = new DateTime(2044, 1, 1), Step = TimeSpan.FromHours(6) };

        public static StudyInterval SiteStudyPeriod = new StudyInterval { Start = new DateTime(2024, 1, 1), Stop = new DateTime(2026, 1, 1), Step = TimeSpan.FromHours(6) };

        public static StudyInterval ViperStudyPeriod1 = new StudyInterval { Start = new DateTime(2022, 12, 1), Stop = new DateTime(2023, 4, 1), Step = TimeSpan.FromHours(2) };

    }

    public class SiteDescriptor
    {
        public bool Southern { get; set; } = true;
        public int Id { get; set; }
        public double Lat { get; set; }
        public double Lon { get; set; }
        public string Name { get; set; } = "Unnamed";

        // If non-zero, then show horizons from this point
        public int Line { get; set; }
        public int Sample { get; set; }

        protected Rectangle _Bounds;
        public Rectangle Bounds
        {
            get
            {
                if (_Bounds.Width > 0 || _Bounds.Height > 0)
                    return _Bounds;
                InMemoryInt16Terrain.GetLineSampleDegrees(Lat, Lon, out int line, out int sample);
                var span = (int)Math.Round(LandingSiteDataset.SiteBoxRadius * 1000f / LandingSiteDataset.MetersPerPixel);
                return new Rectangle(sample - span, line - span, span * 2, span * 2);
            }
            set
            {
                _Bounds = value;
            }
        }

        public Rectangle IdRectangle
        {
            set
            {
                var s = TerrainPatch.DefaultSize;
                var r = value;
                var left = r.Left * s;
                var top = r.Top * s;
                var right = r.Right * s;
                var bottom = r.Bottom * s;
                Bounds = new Rectangle(left, top, right - left, bottom - top);
            }
        }

        protected StudyInterval _studyInterval = null;
        public StudyInterval StudyInterval => _studyInterval != null ? _studyInterval : (Southern ? LandingSiteDataset.SouthStudyInterval : LandingSiteDataset.NorthStudyInterval);
    }

    public class SiteRectangle
    {
        public string Id { get; set; }
        public RectangleF Bounds { get; set; }
    }
}
