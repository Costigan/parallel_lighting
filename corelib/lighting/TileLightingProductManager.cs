using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using viper.corelib.horizon;
using viper.corelib.math;
using viper.corelib.patch;
using viper.corelib.spice;
using viper.corelib.terrain;
using viper.corelib.utilities;

namespace viper.corelib.lighting
{
    /// <summary>
    /// Manage a database of images of average sun and earth exposure.
    /// Create aggregate images.
    /// </summary>
    public class TileLightingProductManager
    {
        public Rectangle Selection;

        protected static Dictionary<(DateTime, DateTime, TimeSpan), List<Vector3d>> _sunVectorCache = new Dictionary<(DateTime, DateTime, TimeSpan), List<Vector3d>>();
        protected static Dictionary<(DateTime, DateTime, TimeSpan), List<Vector3d>> _earthVectorCache = new Dictionary<(DateTime, DateTime, TimeSpan), List<Vector3d>>();

        public static DateTime MetonicCycleStart = new DateTime(2000, 1, 1);
        public static DateTime MetonicCycleStop = new DateTime(2019, 1, 1);
        public static TimeSpan MetonicCycleStep = TimeSpan.FromHours(2);

        public DateTime IntervalStart = new DateTime(2000, 1, 1);
        public DateTime IntervalStop = new DateTime(2019, 1, 1);
        public TimeSpan IntervalStep = TimeSpan.FromHours(2);

        public int ObserverHeightInMeters = 0;
        public float EarthMultipathThreshold = 2f;

        public string SunExposurePath = @"sun_exposure.png";
        public string EarthExposurePath = @"earth_exposure.png";

        public bool GenerateAverageSun = true;
        public bool GenerateAverageEarth = true;

        public float[] LongestNightThresholds = new float[] { 1f, 0.5f, 0.1f };
        public int[] SolarArrayHeights = new int[] { 0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40, 44, 50, 60, 70, 80, 90, 100 };

        public List<TerrainPatch> GetPatchShells(Rectangle? region = null, int observer_height = -1, bool generate_horizons_if_needed = false)
        {
            if (!region.HasValue) region = Selection;
            var r = region.Value;
            if (observer_height < 0) observer_height = ObserverHeightInMeters;

            var result = new List<TerrainPatch>();
            var queue = new List<TerrainPatch>();
            var idr = TerrainPatch.CoveringIdRectangle(r);
            foreach (var id in TerrainPatch.EnumerateIds(idr))
            {
                var patch = TerrainPatch.FromId(id, observer_height);
                (File.Exists(patch.Path) ? result : queue).Add(patch);
            }
            if (queue.Count>0 && generate_horizons_if_needed)
            {
                Console.WriteLine($"Generating horizons for observer {observer_height}");
                var gpuProcessor = ViperEnvironment.Processor as GPUHorizons;
                gpuProcessor.RunQueue(queue, unloadHorizons: true);
                result.AddRange(queue);
            }
            return result;
        }

        IEnumerable<TerrainPatch> LoadShells(List<TerrainPatch> shells) => shells.Select(shell => ViperEnvironment.GetPatch(shell));

        public void GenerateAverageSunAndEarthImages()
        {
            // Generate the needed patches and write them to the filesystem.  This doesn't repeat work
            GenerateAverageSunEarthPatches();

            var shells = GetPatchShells(Selection, ObserverHeightInMeters);  // Generating horizons should have already happened
            if (shells.Count < 1)
                return;

            var map = new LightingColorTable100();  // new ViridisColorTable();

            if (GenerateAverageEarth)
            {
                var dict = new Dictionary<TerrainPatch, Bitmap>();
                foreach (var patch in shells)
                {
                    var path = AverageEarthPath(patch.Line, patch.Sample, IntervalStart, IntervalStop, IntervalStep, ObserverHeightInMeters, EarthMultipathThreshold);
                    if (File.Exists(path))
                    {
                        var bmp = Image.FromFile(path) as Bitmap;
                        Debug.Assert(bmp.PixelFormat == PixelFormat.Format8bppIndexed);
                        dict.Add(patch, bmp);
                    }                    
                }
                using (var earth_image = CombinePatchImages(Selection, dict))
                {
                    map.LoadPalette(earth_image);
                    earth_image.Save(EarthExposurePath, ImageFormat.Png);

                }
                //LandingSiteDataset.FinishSiteImages(EarthExposurePath, $"Average DTE Comm", map);
            }

            if (GenerateAverageSun)
            {
                var dict = new Dictionary<TerrainPatch, Bitmap>();
                foreach (var patch in shells)
                {
                    var path = AverageSunPath(patch.Line, patch.Sample, IntervalStart, IntervalStop, IntervalStep, ObserverHeightInMeters);
                    if (File.Exists(path))
                    {
                        var bmp = Image.FromFile(path) as Bitmap;
                        Debug.Assert(bmp.PixelFormat == PixelFormat.Format8bppIndexed);
                        dict.Add(patch, bmp);
                    }
                }
                using (var sun_image = CombinePatchImages(Selection, dict))
                {
                    map.LoadPalette(sun_image);
                    sun_image.Save(SunExposurePath, ImageFormat.Png);
                    using (var bmp = map.AddLegend(sun_image, steps: 11))
                        bmp.Save(LandingSiteDataset.AppendToFilename(SunExposurePath, "_legend", ".png"), ImageFormat.Png);
                }
                //LandingSiteDataset.FinishSiteImages(SunExposurePath, null, map);  //$"Average Illumination"
            }
        }

        public void GenerateAverageSunEarthPatches()
        {
            List<Vector3d> sunvecs = null, earthvecs = null;
            if (GenerateAverageSun)
                sunvecs = GetSunVectors(IntervalStart, IntervalStop, IntervalStep);
            if (GenerateAverageEarth)
                earthvecs = GetEarthVectors(IntervalStart, IntervalStop, IntervalStep);

            var shells = GetPatchShells(Selection, ObserverHeightInMeters, generate_horizons_if_needed: true);
            shells = shells
                .Where(shell =>
                    (GenerateAverageSun && !File.Exists(AverageSunPath(shell.Line, shell.Sample, IntervalStart, IntervalStop, IntervalStep, ObserverHeightInMeters)))
                 || (GenerateAverageEarth && !File.Exists(AverageEarthPath(shell.Line, shell.Sample, IntervalStart, IntervalStop, IntervalStep, ObserverHeightInMeters, EarthMultipathThreshold))))
                 .ToList();
            var patches = new List<TerrainPatch>();
            for (var i = 0; i < shells.Count; i++)
            {
                Console.Write($"Loading {i}/{shells.Count}: ");
                var patch = ViperEnvironment.GetPatch(shells[i].Id, ObserverHeightInMeters);
                Debug.Assert(patch != null);
                patches.Add(patch);
            }

            // Now, run in parallel
            var stopwatch = new Stopwatch();
            stopwatch.Start();
            var counter = 0;
            var count = (GenerateAverageSun ? patches.Count : 0) + (GenerateAverageEarth ? patches.Count : 0);
            var options = new ParallelOptions { MaxDegreeOfParallelism = ViperEnvironment.MaxDegreeOfParallelism };
            if (GenerateAverageSun)
                Parallel.ForEach(patches, options, patch =>
                {
                    var path = AverageSunPath(patch.Line, patch.Sample, IntervalStart, IntervalStop, IntervalStep, ObserverHeightInMeters);
                    if (!File.Exists(path))
                    {
                        using (var bmp = patch.GenerateAverageSun(sunvecs))
                            bmp.Save(path, ImageFormat.Png);
                        Console.WriteLine($"Wrote {Interlocked.Increment(ref counter)}/{count}: {path}");
                    }
                });
            if (GenerateAverageEarth)
                Parallel.ForEach(patches, options, patch =>
                {
                    var path = AverageEarthPath(patch.Line, patch.Sample, IntervalStart, IntervalStop, IntervalStep, ObserverHeightInMeters, EarthMultipathThreshold);
                    if (!File.Exists(path))
                    {
                        using (var bmp = patch.GenerateAverageEarth(earthvecs, EarthMultipathThreshold))
                            bmp.Save(path, ImageFormat.Png);
                        Console.WriteLine($"Wrote {Interlocked.Increment(ref counter)}/{count}: {path}");
                    }
                });
            stopwatch.Stop();
            Console.WriteLine($"patches calculated in {stopwatch.Elapsed}");
        }

        // For now, 8 bit only
        unsafe Bitmap CombinePatchImages(Rectangle selection, Dictionary<TerrainPatch, Bitmap> dict)
        {
            var selection_width = selection.Width;
            var selection_height = selection.Height;
            var bmp = new Bitmap(selection_width, selection_height, PixelFormat.Format8bppIndexed);
            var bmp_data = bmp.LockBits(new Rectangle(new Point(0, 0), selection.Size), ImageLockMode.WriteOnly, bmp.PixelFormat);
            foreach (var patch in dict.Keys)
            {
                var src_bmp = dict[patch];
                Debug.Assert(src_bmp.PixelFormat == PixelFormat.Format8bppIndexed);
                var xoffset = patch.Sample - selection.Left;
                var yoffset = patch.Line - selection.Top;
                var width = src_bmp.Width;
                var height = src_bmp.Height;
                var src_data = src_bmp.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, src_bmp.PixelFormat);
                for (var row = 0; row < height; row++)
                {
                    var bmp_row = row + yoffset;
                    if (bmp_row >= selection_height || bmp_row < 0) continue;
                    var bmp_rowptr = (byte*)(bmp_data.Scan0 + bmp_data.Stride * bmp_row);
                    var rowptr = (byte*)(src_data.Scan0 + src_data.Stride * row);
                    for (var col = 0; col < width; col++)
                    {
                        // This can be faster if I don't clip per pixel.  Later.
                        var bmp_col = col + xoffset;
                        if (bmp_col > selection_width || bmp_col < 0) continue;
                        bmp_rowptr[bmp_col] = rowptr[col];
                    }
                }
                src_bmp.UnlockBits(src_data);
            }

            bmp.UnlockBits(bmp_data);
            return bmp;
        }

        public void GenerateLongestNightPatches(Rectangle rect, DateTime start, DateTime stop, TimeSpan step, int observer)
        {
            var use_cache = true;
            var shells = GetPatchShells(rect, observer, generate_horizons_if_needed: true);
            if (use_cache)  // true lets cache be used
                shells = shells
                    .Where(shell =>
                        !File.Exists(LongestNightPath(shell.Line, shell.Sample, start, stop, step, observer)))
                     .ToList();
            var patches = new List<TerrainPatch>();
            for (var i = 0; i < shells.Count; i++)
            {
                Console.Write($"Loading {i}/{shells.Count}: ");
                var patch = ViperEnvironment.GetPatch(shells[i].Id, observer);
                Debug.Assert(patch != null);
                patches.Add(patch);
            }

            if (patches.Count == 0)
                return;

            List<Vector3d> sunvecs = GetSunVectors(start, stop, step);

            // Now, run in parallel
            var stopwatch = new Stopwatch();
            stopwatch.Start();
            var counter = 0;
            var count = patches.Count;
            var options = new ParallelOptions { MaxDegreeOfParallelism = ViperEnvironment.MaxDegreeOfParallelism };
            Parallel.ForEach(patches, options, patch =>
            {
                var path = LongestNightPath(patch.Line, patch.Sample, start, stop, step, observer);
                if (!use_cache || !File.Exists(path))
                {
                    var sunfrac = new float[sunvecs.Count];
                    WriteLongestNightsDataset(path, patch, sunvecs, sunfrac, LongestNightThresholds, step);
                    Console.WriteLine($"Wrote {Interlocked.Increment(ref counter)}/{count}: {path}");
                }
            });
            stopwatch.Stop();
            Console.WriteLine($"longest night patches calculated in {stopwatch.Elapsed}");
        }

        void WriteLongestNightsDataset(string path, TerrainPatch patch, List<Vector3d> sunvecs, float[] sunfrac, float[] thresholds, TimeSpan step)
        {
            using (var bw = new BinaryWriter(File.Open(path, FileMode.Create)))
            {
                var width = patch.Width;
                var height = patch.Height;
                bw.Write(patch.Line);
                bw.Write(patch.Sample);
                bw.Write(height);
                bw.Write(width);
                var count = sunvecs.Count;
                bw.Write(thresholds.Length);
                for (var i = 0; i < thresholds.Length; i++)
                    bw.Write(thresholds[i]);
                var buf = new float[count];
                var step_in_hours = (float)step.TotalHours;
                for (var band=0;band<thresholds.Length;band++)
                {
                    var threshold = thresholds[band];
                    for (var row = 0; row < height; row++)
                        for (var col = 0; col<width;col++)
                        {
                            patch.GetLightCurve(row, col, sunvecs, buf);
                            var isDay1 = !(buf[0] >= threshold);
                            var start_of_night = -1;
                            var longest_night = -1;
                            for (var i=0;i<count;i++)
                            {
                                var isDay2 = buf[i] >= threshold;
                                if (isDay1 != isDay2)
                                {
                                    if (isDay1) // Night starts
                                        start_of_night = i;
                                    else if (start_of_night > -1)      // Night ends.  Note this ignores the first partial night if things start in darkness
                                        longest_night = Math.Max(longest_night, i - start_of_night);
                                }
                                isDay1 = isDay2;
                            }
                            var longest_night_in_hours = longest_night < 0 ? count * step_in_hours : longest_night * step_in_hours;
                            Debug.Assert(longest_night_in_hours >= 0f);
                            bw.Write(longest_night_in_hours);

                            if (patch.Id.X == 111 && patch.Id.Y == 123 && longest_night_in_hours < 5 * 24)
                                Console.WriteLine($"here: longest_night_in_hours{longest_night_in_hours}");
                        }
                }
            }
            Debug.Assert((new FileInfo(path)).Length == 196640);
        }

        public void WriteLongestNightsGeotiffs(int site, StudyInterval interval = null, int observer = 0) => WriteLongestNightsGeotiffs(site, interval.Start, interval.Stop, interval.Step, observer);

        public void WriteLongestNightsGeotiffs(int site, DateTime start, DateTime stop, TimeSpan step, int observer)
        {
            var region = LandingSiteDataset.GetSiteBounds(site);
            GenerateLongestNightPatches(region, start, stop, step, observer);  // Ensure the patch files have been created
            var shells = GetPatchShells(region, observer, generate_horizons_if_needed: false);

            var bands = LongestNightThresholds.Select(t => new float[region.Height, region.Width]).ToArray();

            var region_width = region.Width;
            var region_height = region.Height;
            foreach (var patch in shells)
            {
                var path = LongestNightPath(patch.Line, patch.Sample, start, stop, step, observer);
                var src = LongestNightsForPatch.Read(path);
                Debug.Assert(src.BandCount == bands.Length);
                var xoffset = patch.Sample - region.Left;
                var yoffset = patch.Line - region.Top;
                for (var i=0;i<bands.Length;i++)
                {
                    var band_dst = bands[i];
                    var band_src = src.Bands[i];
                    for (var row_src = 0; row_src < TerrainPatch.DefaultSize; row_src++)
                    {
                        var row_dst = row_src + yoffset;
                        if (row_dst >= region_height || row_dst < 0) continue;
                        for (var col_src = 0; col_src < TerrainPatch.DefaultSize; col_src++)
                        {
                            // This can be faster if I don't clip per pixel.  Later.
                            var col_dst = col_src + xoffset;
                            if (col_dst >= region_width || col_dst < 0) continue;
                            band_dst[row_dst,col_dst] = band_src[row_src,col_src];
                        }
                    }
                }
            }

            for (var i=0;i<bands.Length;i++)
            {
                var path = LongestNightPath(site, start, stop, step, observer, LongestNightThresholds[i]);
                GeotiffHelper.WriteArrayAsGeotiff(bands[i], region, path);
                WriteLongestNightTestImage(bands[i], path);
            }
        }

        void WriteLongestNightTestImage(float[,] ary, string path)
        {
            /*
            var height = ary.GetLength(0);
            var width = ary.GetLength(1);
            var min = float.MaxValue;
            var max = float.MinValue;
            for (var row=0;row<height;row++)
                for (var col=0;col<width;col++)
                {
                    var v = ary[row, col];
                    min = Math.Min(min, v);
                    max = Math.Max(max, v);
                }
            var b = min;
            var a = 255f / (max - min);
            var func = v => (byte)((max - v + min) * a));
            */
            var map = new LongestNightColorTableV1();
            Func<float, byte> func = v => (byte)(v / 4f);  // four hour chunks
            using (var bmp = BitmapHelper.ToFormat8bppIndexed(ary, func))
            {
                map.LoadPalette(bmp);
                using (var bmp2 = map.AddLegend(bmp))
                    bmp2.Save(LandingSiteDataset.AppendToFilename(path, "_test", ".png"), ImageFormat.Png);
            }
        }

        internal void WritePixelCounts(int site)
        {
            var bins = BleacherBins().ToList();
            var p = StudyInterval.MazaricoCycle;
            var observers = new int[] { 0, 2 };
            var area = observers.Select(o =>
            {
                var ary = BitmapHelper.Image8bppToByteArray(AverageSunPath(site, p.Start, p.Stop, p.Step, o));
                return bins.Select(bin => CountBetweenPercent(ary, bin) * 20f * 20f).ToList();
            }).ToList();

            // Write table
            var table_path = Path.Combine(LandingSiteDataset.SiteDatasetRoot(site), string.Format(@"site_{0:D3}_sun_vs_surface_area.csv", site));
            using (var sw = new StreamWriter(table_path))
            {
                sw.WriteLine(@"low,high,surface_area_0m,surface_area_2m");
                for (var i = 0; i < bins.Count; i++)
                    sw.WriteLine($"{bins[i].Item1},{bins[i].Item2},{area[0][i]},{area[1][i]}");
            }

            var jet = new JetColorTable();
            List<Color> colors = null;
            foreach (var o in observers)
            {
                using (var bmp = Image.FromFile(AverageSunPath(site, p.Start, p.Stop, p.Step, o)) as Bitmap)
                {
                    jet.LoadPalette(bmp);
                    var pal = bmp.Palette;
                    if (colors==null) colors = BleacherBins256().Select(bin => pal.Entries[bin.Item2 - 1]).ToList();
                    for (var i = 0; i < 128; i++) pal.Entries[i] = Color.Black;
                    foreach (var bin in BleacherBins256())
                    {
                        var c = pal.Entries[bin.Item2 - 1];
                        for (var i = bin.Item1; i < bin.Item2; i++)
                            pal.Entries[i] = c;
                    }
                    bmp.Palette = pal;
                    var table_img_path = Path.Combine(LandingSiteDataset.SiteDatasetRoot(site), string.Format(@"site_{0:D3}_sun@{1}m_vs_surface_area_raw.png", site, o));
                    bmp.Save(table_img_path, ImageFormat.Png);

                    using (var sw = new StreamWriter("surface_area_colors.txt"))
                        foreach (var c in colors)
                            sw.WriteLine($"{c.R},{c.G},{c.B}");

                    var table_img_path_combined = Path.Combine(LandingSiteDataset.SiteDatasetRoot(site), string.Format(@"site_{0:D3}_sun@{1}m_vs_surface_area.png", site, o));
                    using (var bmp1 = CombineImages(bmp))
                        bmp1.Save(table_img_path_combined, ImageFormat.Png);
                }
            }

            int CountBetweenPercent(byte[,] a, (int, int) bin) => CountBetween(a, (int)(bin.Item1 * 255f / 100f), (int)(bin.Item2 * 255f / 100f));
            int CountBetween(byte[,] a, int low, int high) => ValuesOf(a).Count(v => v >= low && v < high);
            IEnumerable<int> ValuesOf(byte[,] a)
            {
                for (var i = 0; i < a.GetLength(0); i++)
                    for (var j = 0; j < a.GetLength(1); j++)
                        yield return a[i, j];
            }
            // Bins from 50%-100% in 5% steps
            IEnumerable<(int, int)> BleacherBins() => Enumerable.Range(0, 10).Select(i => (50 + 5 * i, 50 + 5 * (i + 1)));
            IEnumerable<(int, int)> BleacherBins256()
            {
                var last = 0;
                for (var i=0;i<11;i++)
                {
                    var v = (int)((50 + 5 * i) * 256f / 100f);
                    if (i > 0)
                        yield return (last, v);
                    last = v;
                }
            }

            Bitmap CombineImages(Bitmap bmp)
            {
                Debug.Assert(bmp.Width == 1000 && bmp.Height == 1000);
                var w = bmp.Width;
                var sep = 100;
                var row_height = 30;
                var drop = 20;
                var bot = row_height * 2 + drop;
                var center1 = new Point(w / 2, w / 2);
                var radius1 = 9 * 1000 / 20;

                var center2 = new Point(center1.X + w + sep, center1.Y);
                var radius2 = 1 * 1000 / 20;

                var w2 = (int)(1.4f * radius2);
                var ratio = w / (float)w2;

                var r = new Bitmap(bmp.Width * 2 + sep, bmp.Height + bot, PixelFormat.Format32bppArgb);
                using (var g = Graphics.FromImage(r))
                using (var circle_pen = new Pen(Brushes.White, 2f))
                {
                    g.FillRectangle(Brushes.White, 0, 0, w * 2 + sep, w);
                    g.DrawImageUnscaled(bmp, 0, 0);

                    g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.None;
                    g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.NearestNeighbor;
                    g.SetClip(new Rectangle(w + sep, 0, w, w), System.Drawing.Drawing2D.CombineMode.Replace);
                    g.DrawImage(bmp, new Rectangle(w + sep, 0, w, w), new Rectangle(center1.X - w2, center1.Y - w2, w2 * 2 + 1, w2 * 2 + 1), GraphicsUnit.Pixel);
                    g.ResetClip();

                    g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;
                    g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.High;
                    g.DrawEllipse(circle_pen, center1.X - radius1, center1.Y - radius1, 2 * radius1, 2 * radius1);
                    g.DrawEllipse(circle_pen, center1.X - radius2, center1.Y - radius2, 2 * radius2, 2 * radius2);

                    var radius3 = radius2 * (w / (float)w2) / 2f;
                    g.DrawEllipse(circle_pen, center2.X - radius3, center2.Y - radius3, 2 * radius3, 2 * radius3);

                    // Calculate the areas again
                    var ary = bmp.ToByteArray2();
                    var counts = bins.Select(bin => CountBetweenPercent(ary, bin) * 20 * 20).ToList();
                    counts.Insert(0, 1000 * 1000 * 20 * 20 - counts.Sum());
                    var copied_bins = bins.Select(a => a).ToList();  // Don't modify bins because it generates counts
                    if (counts.Count > copied_bins.Count)
                        copied_bins.Insert(0, (0, 50));
                    if (counts.Count > colors.Count)
                        colors.Insert(0, Color.Black);
                    // Draw the table
                    g.FillRectangle(Brushes.LightGray, 0, w+drop, r.Width, r.Height);  // background for counts
                    g.DrawLine(Pens.Black, 0, w + drop + row_height, r.Width, w + drop + row_height);
                    var cell_width = r.Width / (float)counts.Count;

                    for (var i=0;i< copied_bins.Count;i++)
                    {
                        var bin = copied_bins[i];
                        var count = counts[i];
                        var x1 = cell_width * i;
                        var x2 = cell_width * (i + 1);
                        var y1 = w + drop;
                        using (var brush = new SolidBrush(colors[i]))
                            g.FillRectangle(brush, x1, y1, x2 - x1 + 10, row_height);
                        g.DrawLine(Pens.Black, x1, y1, x1, y1 + w);
                        var percent_brush = i > 0 && i < 6 ? Brushes.Black : Brushes.White;
                        var font = ViperEnvironment.DefaultOverlayFont;
                        g.DrawString($"{bin.Item1}-{bin.Item2}%", font, percent_brush, x1 + 3, y1 + 3);
                        var str = FormatAreaMeasurement(counts[i]);
                        g.DrawString(str, font, Brushes.Black, x1 + 3, y1 + 3 + row_height);
                    }
                }
                return r;
            }
            string FormatAreaMeasurement(int square_meters) => square_meters >= 100000 ? (square_meters / 1000000f).ToString("0.000") + " km2" : square_meters.ToString() + " m2";
        }

        // Generate mazarizo cycle lighting patches for all existing horizon patches

        public void GenerateNewLightingPatches() =>
            GenerateLightingPatches(Directory.EnumerateFiles(ViperEnvironment.HorizonsRoot).Select(filename => TerrainPatch.ReadShell(filename)));

        public void GenerateLightingPatches(IEnumerable<TerrainPatch> patches, bool always = false)
        {
            const float SafeHavenSunThreshold = 0.5f;
            const int observer = 0;
            var sunvecs_mazarico = GetSunVectors(StudyInterval.MazaricoCycle);

            var isSouth = !InMemoryInt16Terrain.Singleton.IsNorth;
            var haven_start = isSouth ? new DateTime(2022, 12, 1) : new DateTime(2023, 5, 1);
            var haven_stop = isSouth ? new DateTime(2023, 4, 1) : new DateTime(2023, 9, 6);
            var haven_step = TimeSpan.FromHours(2);

            var sun_vectors = GetSunVectors(haven_start, haven_stop, haven_step);
            var earth_vectors = GetEarthVectors(haven_start, haven_stop, haven_step);
            var haven_times = new List<DateTime>();
            for (var time = haven_start; time <= haven_stop; time += haven_step)
                haven_times.Add(time);

            var map = new SafeHavenColorTable();
            Func<float, byte> func = v => (byte)Math.Min(255f, v);  // color entries map to hours

            Parallel.ForEach(patches,
                    new ParallelOptions { MaxDegreeOfParallelism = ViperEnvironment.MaxDegreeOfParallelism },
                    shell =>
                    {
                        bool any_processing = false;
                        // Get patch, filter by observer
                        var patch = TerrainPatch.ReadShell(shell.Path);
                        
                        if (always || patch.ObserverHeightInMeters != observer)
                            return;

                        var combined_haven_path = HavenPath(patch.Line, patch.Sample, haven_start, haven_stop, haven_step, observer);
                        if (!File.Exists(combined_haven_path))
                        {
                            any_processing = true;
                            patch.FillPointsAndMatrices(ViperEnvironment.Terrain);  // Fill both because both will be used in the avg sun calc
                            patch.Load();
                            var earth_elevations = earth_vectors.Select(v =>
                            {
                                patch.GetAzEl(v, 0, 0, out float _, out float earth_elevation_rad);
                                return earth_elevation_rad;
                            }).ToList();

                            // Find minima
                            var indices_of_minima_earth_elevation = new List<int>();
                            for (var i = 1; i < earth_vectors.Count - 1; i++)
                                if (haven_start <= haven_times[i] && haven_times[i] <= haven_stop && earth_elevations[i - 1] > earth_elevations[i] && earth_elevations[i] < earth_elevations[i + 1])
                                    indices_of_minima_earth_elevation.Add(i);

                            var combined = new float[TerrainPatch.DefaultSize, TerrainPatch.DefaultSize];
                            var month_images = Enumerable.Range(0, indices_of_minima_earth_elevation.Count).Select(i => new float[TerrainPatch.DefaultSize, TerrainPatch.DefaultSize]).ToList();

                            for (var row = 0; row < TerrainPatch.DefaultSize; row++)
                            {
                                for (var col = 0; col < TerrainPatch.DefaultSize; col++)
                                {
                                    var has_sun = patch.GetLightCurve(row, col, sun_vectors).Select(sun => sun >= SafeHavenSunThreshold).ToList();
                                    var has_comm = patch.GetEarthOverHorizonCurve(row, col, earth_vectors).Select(rise => rise >= EarthMultipathThreshold).ToList();
                                    var month_steps = SafeHavenGenerator.MaxShadowWithoutComm2(indices_of_minima_earth_elevation, has_sun, has_comm);
                                    combined[row, col] = 2f * month_steps.Max();
                                    for (var i = 0; i < month_steps.Count; i++)
                                        month_images[i][row, col] = 2f * month_steps[i];  // hours (steps are in units of 2 hours)
                                }
                            }

                            using (var bmp = BitmapHelper.ToFormat8bppIndexed(combined, func))
                            {
                                map.LoadPalette(bmp);
                                bmp.Save(LandingSiteDataset.AppendToFilename(combined_haven_path, "", ".png"), ImageFormat.Png);
                            }
                            for (var i = 0; i < month_images.Count; i++)
                                using (var bmp = BitmapHelper.ToFormat8bppIndexed(month_images[i], func))
                                {
                                    var path = HavenPath(patch.Line, patch.Sample, haven_times[indices_of_minima_earth_elevation[i]], observer);
                                    map.LoadPalette(bmp);
                                    bmp.Save(LandingSiteDataset.AppendToFilename(path, "", ".png"), ImageFormat.Png);
                                }
                        }

                        if (true)  // disable average sun here
                        {
                            Debug.Assert(patch.ObserverHeightInMeters == observer);
                            var m = StudyInterval.MazaricoCycle;
                            var mazarico_path = AverageSunPath(patch.Line, patch.Sample, m.Start, m.Stop, m.Step, observer);

                            if (always || !File.Exists(mazarico_path))
                            {
                                any_processing = true;
                                patch.FillPointsAndMatrices(ViperEnvironment.Terrain);  // Fill both because both will be used in the avg sun calc
                                if (patch.Horizons == null) patch.Load();
                                using (var bmp = patch.GenerateAverageSun(sunvecs_mazarico))
                                    bmp.Save(mazarico_path, ImageFormat.Png);
                            }
                        }

                        // throw the patch away.  Needed when we're enumerating the selected patches
                        patch.Horizons = null;
                        patch.Points = null;
                        patch.Matrices = null;

                        if (any_processing)
                            Console.WriteLine($"Processed {patch.DefaultPath}");
                        else
                            Console.WriteLine($"Skipped {patch.DefaultPath}");
                    });
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="site"></param>
        /// <param name="start"></param>
        /// <param name="stop"></param>
        /// <param name="step"></param>
        /// <param name="observer"></param>
        /// <param name="landing_slope_threshold"></param>
        /// <param name="error_span">Error circle radius in pixels</param>
        public unsafe void GenerateFullSunDurationMaps(int site, DateTime start, DateTime stop, TimeSpan step, int observer, float? landing_slope_threshold = null, int? error_span = null)
        {
            var rect = LandingSiteDataset.GetSiteBounds(site);
            byte[,,] in_sun = null;
            byte[,] max_days = null;

            var path = FullSunForTimeMapPath(site, start, stop, step, observer);
            path = path.FilenameAppend("_no_dte");
            if (File.Exists(path))
                return;
            if (max_days == null)
                GenerateMaxDays();
            var img = FillBitmap();

            var days_array = new int[] { 7, 8, 10, 14, 21, 28 };
            var colors = new Color[] { Color.Blue, Color.Plum, Color.Green, Color.LightBlue, Color.Pink, Color.Yellow };

            var p = img.Palette;
            for (var i = 0; i < 256; i++)
            {
                p.Entries[i] = Color.Black;
                for (var j = days_array.Length - 1; j >= 0; j--)
                {
                    if (i >= days_array[j])
                    {
                        p.Entries[i] = colors[j];
                        break;
                    }
                }
            }
            img.Palette = p;

            img.Save(path, ImageFormat.Png);

            if (landing_slope_threshold.HasValue && error_span.HasValue)
                WriteWithSlopeConstraint(img, landing_slope_threshold.Value, error_span.Value);

            img.Dispose();

            void GenerateMaxDays()
            {
                var width = rect.Width;
                var height = rect.Height;
                max_days = new byte[height, width];
                if (in_sun == null)
                {
                    in_sun = GenerateFullSunBitArray(site, start, stop, step, observer);
                    if (in_sun == null)
                        in_sun = ReadFullSunBitArray(FullSunBitArrayCachePath(site, start, start, step, observer));
                }
                var base_pt = TerrainPatch.IdToLineSample(TerrainPatch.LineSampleToId(rect.Location));
                var row_offset = rect.Y - base_pt.Y;
                var col_offset = rect.X - base_pt.X;
                var layers = in_sun.GetLength(0);
                var step_in_hours = (int)step.TotalHours;
                var options = new ParallelOptions { MaxDegreeOfParallelism = ViperEnvironment.MaxDegreeOfParallelism };
                Parallel.For(0, height, row =>
                {
                    var row1 = row + row_offset;
                    for (var col = 0; col < width; col++)
                    {
                        var col1 = col + col_offset;
                        var max_hours = -1;
                        var current_hours = 0;
                        var last_s = 0;
                        for (var layer = 0; layer < layers; layer++)
                        {
                            var s = in_sun[layer, row1, col1];
                            if (s == 1)
                            {
                                current_hours += step_in_hours;
                            }
                            else if (last_s == 1)
                            {
                                max_hours = Math.Max(max_hours, current_hours);
                                current_hours = 0;
                            }
                            last_s = s;
                        }
                        var longest_run_of_days = max_hours / 24;

                        //if (longest_run_of_days >= 80)
                        //{
                         //   Console.WriteLine($"days={longest_run_of_days}");
                        //}

                        max_days[row, col] = longest_run_of_days > 255 ? (byte)255 : (byte)longest_run_of_days;
                    }
                });
            }

            Bitmap FillBitmap()
            {
                var width = rect.Width;
                var height = rect.Height;
                var bmp = new Bitmap(width, height, PixelFormat.Format8bppIndexed);
                var bmpdata = bmp.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, bmp.PixelFormat);
                for (var row = 0; row < height; row++)
                {
                    var rowptr = (byte*)(bmpdata.Scan0 + row * bmpdata.Stride);
                    for (var col = 0; col < width; col++)
                        rowptr[col] = max_days[row, col];
                }
                bmp.UnlockBits(bmpdata);
                return bmp;
            }

            void WriteWithSlopeConstraint(Bitmap bmp, float slope_threshold, int span)
            {
                Debug.Assert(bmp.PixelFormat == PixelFormat.Format8bppIndexed);
                var slope_path = Path.Combine(LandingSiteDataset.SiteDatasetRoot(site), $"site_{site.ToString("000")}_slope.tif");
                var slope = GeotiffHelper.ReadGeotiffAsFloatArray(slope_path, 89f);  // no data -> steep
                Debug.Assert(bmp.Height == slope.GetLength(0) && bmp.Width == slope.GetLength(1));

                // Calculate offsets
                var d2 = span * span;
                var offsets = new List<Point>((2*span+1)*(2*span+1));
                for (var y = -span; y <= span; y++)
                    for (var x = -span; x <= span; x++)
                        if (d2 >= y * y + x * x)
                            offsets.Add(new Point(x, y));

                var count = offsets.Count;
                var height = slope.GetLength(0);
                var width = slope.GetLength(1);
                var legal = new byte[height, width];
                for (var row=0;row<height;row++)
                    for (var col=0;col<width;col++)
                    {
                        for (var i=0;i<count;i++)
                        {
                            var offset = offsets[i];
                            var row1 = row + offset.Y;
                            if (row1 < 0 || row1 >= height)
                                continue;
                            var col1 = col + offset.X;
                            if (col1 < 0 || col1 >= width)
                                continue;
                            if (slope[row1,col1] > slope_threshold)
                            {
                                legal[row, col] = 0;
                                goto continue_loop;
                            }
                        }
                        legal[row, col] = 1;
                        continue_loop: { }
                    }

                // Save the legal landing sites as an image
                using (var legal_bmp = BitmapHelper.ToFormat8bppIndexed(legal))
                {
                    var pal = legal_bmp.Palette;
                    pal.Entries[0] = Color.Black;
                    pal.Entries[1] = Color.Green;
                    legal_bmp.Palette = pal;
                    var legal_bmp_path = Path.Combine(LandingSiteDataset.SiteDatasetRoot(site), $"site_{site.ToString("000")}_legal_landing_sites_{(20*span).ToString("000")}m_{slope_threshold.ToString("00.0")}deg.png");
                    legal_bmp.Save(legal_bmp_path, ImageFormat.Png);
                }

                // Zero out all of the illegal landing sites
                var bmpdata = bmp.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, bmp.PixelFormat);
                for (var row = 0; row < height; row++)
                {
                    var rowptr = (byte*)(bmpdata.Scan0 + row * bmpdata.Stride);
                    for (var col=0;col<width;col++)
                        if (legal[row, col] == 0)
                            rowptr[col] = 0;
                }
                bmp.UnlockBits(bmpdata);

                // Write the new map
                var path_with_slope_constraint = FullSunForTimeMapPathWithSlope(site, start, stop, step, observer, span, slope_threshold);
                bmp.Save(path_with_slope_constraint, ImageFormat.Png);

                // Write overlaid on a hillshade
                using (var overlaid = new Bitmap(width, height, PixelFormat.Format32bppArgb))
                using (var hillshade = Image.FromFile(Path.Combine(LandingSiteDataset.SiteDatasetRoot(site), $"site_{site.ToString("000")}_hillshade.tif")) as Bitmap)
                using (var g = Graphics.FromImage(overlaid))
                {
                    Debug.Assert(hillshade.Width == width && hillshade.Height == height);
                    g.DrawImageUnscaled(hillshade, 0, 0);
                    g.CompositingMode = System.Drawing.Drawing2D.CompositingMode.SourceOver;
                    var imageatt = new ImageAttributes();
                    imageatt.SetColorMatrix(new ColorMatrix { Matrix33 = 0.5f });
                    g.DrawImage(bmp, new Rectangle(0, 0, width, height), 0, 0, width, height, GraphicsUnit.Pixel, imageatt);
                    var overlaid_path = LandingSiteDataset.AppendToFilename(path_with_slope_constraint, "_overlaid");
                    overlaid.Save(overlaid_path, ImageFormat.Png);
                }
            }
        }

        public byte[,,] GenerateFullSunBitArray(int site, DateTime start, DateTime stop, TimeSpan step, int observer)
        {
            const float threshold = 1f;
            var path = FullSunBitArrayCachePath(site, start, start, step, observer);
            if (File.Exists(path))
                return null;
            var rect = LandingSiteDataset.GetSiteBounds(site);
            GetPatchShells(rect, observer, generate_horizons_if_needed: true);  // generated horizons that haven't been generated before
            var idrect = TerrainPatch.CoveringIdRectangle(rect);
            var patches = TerrainPatch.EnumerateIds(idrect).Select(id => ViperEnvironment.GetPatch(id, observer)).ToList();
            Debug.Assert(patches.Count == idrect.Width * idrect.Height);

            var id_origin = patches[0].Id;
            var id_width = 1 + patches.Max(p => p.Id.X - id_origin.X);
            var id_height = 1 + patches.Max(p => p.Id.Y - id_origin.Y);
            var ary_width = TerrainPatch.DefaultSize * id_width;
            var ary_height = TerrainPatch.DefaultSize * id_height;

            List<Vector3d> sunvecs = GetSunVectors(start, stop, step);
            var earthvecs = GetEarthVectors(start, stop, step);
            Debug.Assert(sunvecs.Count == earthvecs.Count);

            var ary = new byte[sunvecs.Count, ary_height, ary_width];

            // Now, run in parallel
            var stopwatch = new Stopwatch();
            stopwatch.Start();
            var counter = 0;
            var count = patches.Count;
            var options = new ParallelOptions { MaxDegreeOfParallelism = ViperEnvironment.MaxDegreeOfParallelism };
            Parallel.ForEach(patches, options, patch =>
            {
                var sun = new float[sunvecs.Count];
                var earth_over_horizon = new float[sunvecs.Count];
                var width = patch.Width;
                var height = patch.Height;
                var row_offset = TerrainPatch.DefaultSize * (patch.Id.Y - id_origin.Y);
                var col_offset = TerrainPatch.DefaultSize * (patch.Id.X - id_origin.X);
                for (var row = 0; row < height; row++)
                {
                    var ary_row = row + row_offset;
                    for (var col = 0; col < width; col++)
                    {
                        var ary_col = col + col_offset;
                        patch.GetLightCurve(row, col, sunvecs, sun);

                        // Debugging
                        /*
                        var flag = true;
                        if (flag)
                        {
                            var threshold_count = (int)(80d / (step.TotalHours/24d));
                            var first = ContainsRunOverThreshold(sun, 1f, threshold_count);
                            if (first>=0)
                            {
                                Console.WriteLine($"Found one id={patch.Id} row={row} col={col} first={first} step_hours={step.TotalHours}");
                                WriteSunArray(sun, row, col, patch);
                                flag = false;
                            }
                        }*/

                        patch.GetEarthOverHorizonCurve(row, col, earthvecs, earth_over_horizon);
                        for (var i = 0; i < sun.Length; i++)
                            ary[i, ary_row, ary_col] = earth_over_horizon[i] < threshold ? (byte)0 : (sun[i] >= threshold ? (byte)1 : (byte)0);
                            //ary[i, ary_row, ary_col] = (sun[i] >= threshold ? (byte)1 : (byte)0);
                    }
                }
                Console.WriteLine($"Wrote {Interlocked.Increment(ref counter)}/{count}: {path}");
            });

            WriteFullSunBitArray(ary, path);
            stopwatch.Stop();
            Console.WriteLine($"Full Sun Bit Array generated in {stopwatch.Elapsed}");

            if (false)
            {
                var ary1 = ReadFullSunBitArray(path);
                for (var level = 0; level < ary.GetLength(0); level++)
                    for (var row = 0; row < ary.GetLength(1); row++)
                        for (var col = 0; col < ary.GetLength(2); col++)
                            if (ary[level, row, col] != ary1[level, row, col])
                                Console.WriteLine(@"Error in full sun bit array");
            }

            return ary;

            // Start of first run over c length
            int ContainsRunOverThreshold(float[] a, float t, int c)
            {
                int first = 0;
                for (var i=0;i<a.Length;i++)
                {
                    if (a[i] < t)
                        first = i;
                    if (i - first > c)
                        return first;
                }
                return -1;
            }

            void WriteSunArray(float[] a, int row, int col, TerrainPatch patch)
            {
                var pt = TerrainPatch.IdToLineSample(patch.Id);
                using (var sw = new StreamWriter(File.Create($"sun_{pt.Y + row:D3}_{pt.X + col:D3}.csv")))
                {
                    sw.WriteLine("sun");
                    for (var i = 0; i < a.Length; i++)
                        sw.WriteLine(a[i]);
                }
            }
        }

        void WriteFullSunBitArray(byte[,,] ary, string path)
        {
            Debug.Assert(ary.GetLength(2).Rem(8) == 0);
            using (var bw = new BinaryWriter(File.OpenWrite(path)))
            {
                var levels = ary.GetLength(0);
                var height = ary.GetLength(1);
                var width = ary.GetLength(2);
                bw.Write(levels);
                bw.Write(height);
                bw.Write(width);
                var width8 = width / 8;
                for (var level = 0; level < levels; level++)
                    for (var row = 0; row < height; row++)
                        for (var col8 = 0; col8 < width8; col8++)
                        {
                            var b = 0;
                            for (var bit = 0; bit < 8; bit++)
                            {
                                var v = ary[level, row, col8 * 8 + bit];
                                b |= (v & 1) << bit;
                            }
                            bw.Write((byte)b);
                        }
            }
        }

        byte[,,] ReadFullSunBitArray(string path)
        {
            using (var br = new BinaryReader(File.OpenRead(path)))
            {
                var levels = br.ReadInt32();
                var height = br.ReadInt32();
                var width = br.ReadInt32();
                Debug.Assert(width.Rem(8) == 0);
                var width8 = width / 8;
                var ary = new byte[levels, height, width];
                for (var level = 0; level < levels; level++)
                    for (var row = 0; row < height; row++)
                        for (var col8 = 0; col8 < width8; col8++)
                        {
                            var b = br.ReadByte();
                            for (var bit = 0;bit<8;bit++)
                                ary[level, row, col8 * 8 + bit] = (byte)((b >> bit) & 1);
                        }
                return ary;
            }
        }

        string FullSunBitArrayCachePath(int site, DateTime start, DateTime stop, TimeSpan step, int observer) 
            => Path.Combine(LandingSiteDataset.SiteDatasetRoot(site),$"site_{site.ToString("000")}_full_sun_bitmap_cache_{start.ToString("yyyyMMdd")}_{stop.ToString("yyyyMMdd")}_{step.ToString("hh")}_{observer.ToString("000")}.bin");

        string FullSunForTimeMapPath(int site, DateTime start, DateTime stop, TimeSpan step, int observer)
            => Path.Combine(LandingSiteDataset.SiteDatasetRoot(site), $"site_{site.ToString("000")}_full_sun_for_time_{start.ToString("yyyyMMdd")}_{stop.ToString("yyyyMMdd")}_{step.ToString("hh")}_{observer.ToString("000")}.png");

        string FullSunForTimeMapPathWithSlope(int site, DateTime start, DateTime stop, TimeSpan step, int observer, int span, float slope_threshold)
            => Path.Combine(LandingSiteDataset.SiteDatasetRoot(site), $"site_{site.ToString("000")}_full_sun_for_time_with_slope_{start.ToString("yyyyMMdd")}_{stop.ToString("yyyyMMdd")}_{step.ToString("hh")}_{observer.ToString("000")}_{(20 * span).ToString("000")}m_{slope_threshold.ToString("00.0")}deg.png");

        #region Path creators

        // These are used for the individual patch files

        public static string AverageSunFilename(int line, int sample, DateTime start, DateTime stop, TimeSpan step, int observer) =>
            string.Format("avg_sun_{0:D5}_{1:D5}_{2}_{3}_{4:D3}_{5:D3}.png",
                line, sample, start.ToString("yyyyMMdd"), stop.ToString("yyyyMMdd"), (int)Math.Round(step.TotalHours), observer);
        public static string AverageSunPath(int line, int sample, DateTime start, DateTime stop, TimeSpan step, int observer) =>
            Path.Combine(ViperEnvironment.AvgPatchRoot, AverageSunFilename(line, sample, start, stop, step, observer));

        public static string AverageEarthFilename(int line, int sample, DateTime start, DateTime stop, TimeSpan step, int observer, float threshold) =>
            string.Format("avg_earth_{0:D5}_{1:D5}_{2}_{3}_{4:D3}_{5:D3}_{6}.png",
                line, sample, start.ToString("yyyyMMdd"), stop.ToString("yyyyMMdd"), (int)Math.Round(step.TotalHours), observer, threshold.ToString("00.00"));
        public static string AverageEarthPath(int line, int sample, DateTime start, DateTime stop, TimeSpan step, int observer, float threshold) =>
            Path.Combine(ViperEnvironment.AvgPatchRoot, AverageEarthFilename(line, sample, start, stop, step, observer, threshold));

        public static string HavenFilename(int line, int sample, DateTime start, DateTime stop, TimeSpan step, int observer) =>
            string.Format("haven_interval_{0:D5}_{1:D5}_{2}_{3}_{4:D3}_{5:D3}.png",
                line, sample, start.ToString("yyyyMMdd"), stop.ToString("yyyyMMdd"), (int)Math.Round(step.TotalHours), observer);
        public static string HavenPath(int line, int sample, DateTime start, DateTime stop, TimeSpan step, int observer) =>
            Path.Combine(ViperEnvironment.HavenPatchRoot, HavenFilename(line, sample, start, stop, step, observer));

        public static string HavenFilename(int line, int sample, DateTime earth_minimum, int observer) =>
            string.Format("haven_overnight_{0:D5}_{1:D5}_{2}_{3:D3}.png",
                line, sample, earth_minimum.ToString("yyyyMMdd"), observer);
        public static string HavenPath(int line, int sample, DateTime earth_minimum, int observer) =>
            Path.Combine(ViperEnvironment.HavenPatchRoot, HavenFilename(line, sample, earth_minimum, observer));

        public static string LongestNightFilename(int line, int sample, DateTime start, DateTime stop, TimeSpan step, int observer) =>
            string.Format("long_night_{0:D5}_{1:D5}_{2}_{3}_{4:D3}_{5:D3}.bin",
                line, sample, start.ToString("yyyyMMdd"), stop.ToString("yyyyMMdd"), (int)Math.Round(step.TotalHours), observer);
        public static string LongestNightPath(int line, int sample, DateTime start, DateTime stop, TimeSpan step, int observer) =>
            Path.Combine(ViperEnvironment.AvgPatchRoot, LongestNightFilename(line, sample, start, stop, step, observer));

        // These are used for site files

        public static string AverageSunFilename(int site, DateTime start, DateTime stop, TimeSpan step, int observer) =>
            string.Format("site_{0:D3}_avg_sun_start{1}_stop{2}_step{3:D3}_observer{4:D3}.png",
                site, start.ToString("yyyyMMdd"), stop.ToString("yyyyMMdd"), (int)Math.Round(step.TotalHours), observer);
        public static string AverageSunPath(int site, DateTime start, DateTime stop, TimeSpan step, int observer) =>
            Path.Combine(LandingSiteDataset.SiteDatasetRoot(site), AverageSunFilename(site, start, stop, step, observer));

        public static string AverageEarthFilename(int site, DateTime start, DateTime stop, TimeSpan step, int observer, float threshold) =>
            string.Format("site_{0:D3}_avg_earth_start{1}_stop{2}_step{3:D3}_observer{4:D3}_threshold{5}.png",
                site, start.ToString("yyyyMMdd"), stop.ToString("yyyyMMdd"), (int)Math.Round(step.TotalHours), observer, threshold.ToString("00.00"));
        public static string AverageEarthPath(int site, DateTime start, DateTime stop, TimeSpan step, int observer, float threshold) =>
            Path.Combine(LandingSiteDataset.SiteDatasetRoot(site), AverageEarthFilename(site, start, stop, step, observer, threshold));

        public static string LongestNightFilename(int site, DateTime start, DateTime stop, TimeSpan step, int observer, float threshold) =>
            string.Format("site_{0:D3}_longest_night_start{1}_stop{2}_step{3:D3}_observer{4:D3}_threshold{5}.tif",
                site, start.ToString("yyyyMMdd"), stop.ToString("yyyyMMdd"), (int)Math.Round(step.TotalHours), observer, threshold.ToString("00.00"));
        public static string LongestNightPath(int site, DateTime start, DateTime stop, TimeSpan step, int observer, float threshold) =>
                Path.Combine(LandingSiteDataset.SiteDatasetRoot(site), LongestNightFilename(site, start, stop, step, observer, threshold));

        #endregion

        #region Height to sun attempt #1

        public void GenerateSolarArrayHeightImages(Point id, DateTime? start = null, DateTime? stop = null, TimeSpan? step = null)
        {
            var start1 = start.HasValue ? start.Value : new DateTime(2024, 1, 1);
            var stop1 = stop.HasValue ? stop.Value : new DateTime(2026, 1, 1);
            var step1 = step.HasValue ? step.Value : TimeSpan.FromHours(6);
            var sunvecs = GetSunVectors(start1, stop1, step1);

            var patches = SolarArrayHeights.Select(h => TerrainPatch.FromId(id, h)).ToList();
            var queue = patches.Where(p => !File.Exists(p.Path)).ToList();

            if (queue.Count > 0)
            {
                Console.WriteLine($"Generating horizons observer heights: {string.Join(",", queue.Select(p => p.ObserverHeightInMeters))}");
                var gpuProcessor = ViperEnvironment.Processor as GPUHorizons;
                gpuProcessor.RunQueue(queue, unloadHorizons: false);
                patches.AddRange(queue);
            }
            var options = new ParallelOptions { MaxDegreeOfParallelism = ViperEnvironment.MaxDegreeOfParallelism };
            var counter = 0;
            var count = patches.Count;
            var stopwatch = Stopwatch.StartNew();
            Parallel.ForEach(patches, options, patch =>
            {
                var path = AverageSunPath(patch.Line, patch.Sample, start1, stop1, step1, patch.ObserverHeightInMeters);
                if (!File.Exists(path))
                {
                    using (var bmp = patch.GenerateAverageSun(sunvecs))
                        bmp.Save(path, ImageFormat.Png);
                    Console.WriteLine($"Wrote {Interlocked.Increment(ref counter)}/{count}: {path}");
                }
            });
            Console.WriteLine($"Generating average sun images took {stopwatch.Elapsed}");
        }

        public void GenerateHeightVsAverageSun(int site, DateTime? start = null, DateTime? stop = null, TimeSpan? step = null)
        {
            var start1 = start.HasValue ? start.Value : new DateTime(2024, 1, 1);
            var stop1 = stop.HasValue ? stop.Value : new DateTime(2026, 1, 1);
            var step1 = step.HasValue ? step.Value : TimeSpan.FromHours(6);
            var sunvecs = GetSunVectors(start1, stop1, step1);

            var id = TerrainPatch.LineSampleToId(LandingSiteDataset.GetSiteCenter(site));

            GenerateSolarArrayHeightImages(id, start1, stop1, step1);

            var patches = SolarArrayHeights.Select(observer => TerrainPatch.FromId(id, observer)).ToList();
            var filenames = patches.Select(patch => (patch.ObserverHeightInMeters, AverageSunPath(patch.Line, patch.Sample, start1, stop1, step1, patch.ObserverHeightInMeters))).ToList();
            var arrays = filenames.Select(pair =>
            {
                using (var bmp = Image.FromFile(pair.Item2) as Bitmap)
                    return (pair.Item1, bmp.ToByteArray2());
            }).OrderBy(pair => pair.Item1).ToList();

            var twometer = arrays.First(pair => pair.Item1 == 2f).Item2;
            var best_10 = ByteArrayEnumeration(twometer).OrderBy(cell => -cell.Item2).Take(100).ToList();
            var path = Path.Combine(LandingSiteDataset.SiteDatasetRoot(site), string.Format(@"site_{0:D3}_height_vs_average_sun.csv", site));
            using (var sw = new StreamWriter(path))
            {
                sw.Write(@"height");
                for (var i = 0; i < best_10.Count; i++)
                {
                    var pt = best_10[i].Item1;
                    sw.Write($",cell_[{pt.X}_{pt.Y}]");
                }
                sw.WriteLine();
                for (var j = 0; j < arrays.Count; j++)
                {
                    sw.Write(arrays[j].Item1);
                    for (var i = 0; i < best_10.Count; i++)
                    {
                        var pt = best_10[i].Item1;
                        sw.Write($",{arrays[j].Item2[pt.Y, pt.X] / 255f}");
                    }
                    sw.WriteLine();
                }
            }
        }

        IEnumerable<(Point,byte)> ByteArrayEnumeration(byte[,] a)
        {
            for (var row = 0; row < a.GetLength(0); row++)
                for (var col = 0; col < a.GetLength(1); col++)
                    yield return (new Point(col, row), a[row, col]);
        }

        // 255 = no data
        // 254 = higher than any probe
        Bitmap GenerateSolarArrayHeightForFractionalSun(Rectangle r, float sunfrac, DateTime? start = null, DateTime? stop = null, TimeSpan? step = null, bool skip_empties = true)
        {
            var start1 = start.HasValue ? start.Value : new DateTime(2024, 1, 1);
            var stop1 = stop.HasValue ? stop.Value : new DateTime(2026, 1, 1);
            var step1 = step.HasValue ? step.Value : TimeSpan.FromHours(6);
            var terrain_size = TerrainPatch.DefaultSize;

            var reversed_heights = SolarArrayHeights.Reverse().ToList();
            var threshold = (byte)(255 * sunfrac);

            var idrect = TerrainPatch.CoveringIdRectangle(r);
            var target = new Bitmap(r.Width, r.Height, PixelFormat.Format8bppIndexed);
            target.SetPixels8Bit((row, col) => 255);
            foreach (var id in TerrainPatch.EnumerateIds(idrect))
            {
                foreach (var height in reversed_heights)
                {
                    var patch = TerrainPatch.FromId(id, height);
                    var path = AverageSunPath(patch.Line, patch.Sample, start1, stop1, step1, patch.ObserverHeightInMeters);
                    var path_exists = File.Exists(path);
                    if (!path_exists && skip_empties)
                        continue;
                    if (!path_exists)
                    {
                        GenerateSolarArrayHeightImages(id, start1, stop1, step1);
                        if (!File.Exists(path))
                        {
                            Console.WriteLine($"Tried to generate {path} but failed.  Skipping.");
                            goto continue_ids;
                        }
                    }
                    using (var source = Image.FromFile(path) as Bitmap)
                        target.CombineInto(source, new Point(id.X * terrain_size - r.Left, id.Y * terrain_size - r.Top),
                            (row, col, srcval, tarval)
                            =>
                            {
                                var result = srcval >= threshold ? (byte)height : tarval;
                                //if (srcval >= threshold)
                                //    Console.WriteLine("here");
                                return result;
                            });
                }
                continue_ids: { }
            }

            LoadSolarArrayHeightColorTable(target);

            return target;
        }

        void LoadSolarArrayHeightColorTable(Bitmap target)
        {
            var p = target.Palette;
            for (var i=0;i<10;i++)
                p.Entries[i] = Color.FromArgb(255-i*20, 255, 255 - i * 20);
            for (var i = 10; i < 20; i++)
                p.Entries[i] = Color.FromArgb(55, 255 - (i-10) * 20, 55);
            for (var i = 21; i < 255; i++)
                p.Entries[i] = Color.DarkSlateGray;
            p.Entries[255] = Color.FromArgb(255, 0, 0, 0);
            target.Palette = p;
        }

        public void WriteSolarArrayHeightForFractionalSun(string path, Rectangle r, float sunfrac, DateTime? start = null, DateTime? stop = null, TimeSpan? step = null, bool skip_empties = true)
        {
            using (var bmp = GenerateSolarArrayHeightForFractionalSun(r, sunfrac, start, stop, step, skip_empties))
            {
                bmp.Save(path, ImageFormat.Png);

                var histogram_path = LandingSiteDataset.AppendToFilename(path, "_histogram", ".csv");
                using (var sw = new StreamWriter(histogram_path))
                    foreach (var b in bmp.GetByteHistogram())
                        sw.WriteLine(b);
            }
        }

        #endregion

        #region Random

        public void CompareFincannon()
        {
            var line = 15503;
            var sample = 15065;
            var id = TerrainPatch.LineSampleToId(line, sample);
            var shell = TerrainPatch.FromId(id, 2);
            shell.Read(shell.DefaultPath);
            using (var sw = new StreamWriter(@"fincannon_comparison.csv"))
            {
                sw.WriteLine(@"Date, sun");
                var start = new DateTime(2023, 1, 1);
                var stop = new DateTime(2026, 1, 1);
                var step = TimeSpan.FromHours(1);
                var times = new List<DateTime>();
                for (var time = start; time <= stop; time += step)
                    times.Add(time);
                var sunvecs = times.Select(time => CSpice.SunPosition(time)).ToList();
                var light_curve = shell.GetLightCurve(line - shell.Line, sample - shell.Sample, sunvecs);
                for (var i = 0; i < times.Count; i++)
                {
                    sw.Write(times[i].ToString("MM/dd/yyyy HH:00:00"));
                    sw.WriteLine($", {light_curve[i]}");
                }
            }
        }

        public void CompareFincannon2()
        {
            var sites = new[] {
                (name: "Shackleton Rim Site", lat: -89.781, lon: 203.844 ),
                //(name: "Spudis Ridge Site",lat: -89.439,lon: 222.855 ),
                //(name: "De Gerlache Rim Site",lat: -88.688,lon: 291.679 ),
                //(name: "Plateau Site",lat:  -88.809, lon: 123.723 )
            };
            var site_data = LoadPatches(sites);
            foreach (var s in site_data)
                s.patch.FillPointsAndMatrices(terrain.InMemoryInt16Terrain.Singleton);

            var site_data1 = site_data.ToList();
            for (var l=-1;l<=1;l++)
                for(var s=-1;s<=1;s++)
                    if (l != 0 || s != 0)
                    {
                        var z = site_data[0];
                        site_data1.Add(($"{z.name}_{l}_{s}", z.lat, z.lon, z.line + l, z.sample + s, z.patch));
                    }
            site_data = site_data1.ToArray();     

            using (var sw = new StreamWriter(@"fincannon_comparison.csv"))
            {
                // Header
                sw.Write(@"Date");
                foreach (var s in site_data) sw.Write($",{s.name}");
                sw.WriteLine();

                var start = new DateTime(2023, 1, 1);
                var stop = new DateTime(2026, 1, 1);
                var step = TimeSpan.FromHours(1);
                var times = new List<DateTime>();
                for (var time = start; time <= stop; time += step)
                    times.Add(time);
                var sunvecs = times.Select(time => CSpice.SunPosition(time)).ToList();
                var curves = site_data.AsParallel().Select(s => s.patch.GetLightCurve(s.line - s.patch.Line, s.sample - s.patch.Sample, sunvecs));
                for (var i=0;i<times.Count;i++)
                {
                    sw.Write(times[i].ToString("MM/dd/yyyy HH:00:00"));
                    foreach (var c in curves)
                        sw.Write($",{c[i]}");
                    sw.WriteLine();
                }
            }

            (string name, double lat, double lon, int line, int sample, TerrainPatch patch)[] LoadPatches((string name, double lat, double lon)[] sites1)
            {
                var id2patch = new Dictionary<Point, TerrainPatch>();
                var result = new List<(string, double, double, int, int, TerrainPatch)>();
                foreach (var s in sites1)
                {
                    terrain.InMemoryInt16Terrain.GetLineSample(s.Item2 * Math.PI / 180d, s.Item3 * Math.PI / 180d, out int line, out int sample);
                    var id = TerrainPatch.LineSampleToId(line, sample);
                    if (id2patch.ContainsKey(id))
                        result.Add((s.Item1, s.Item2, s.Item3, line, sample, id2patch[id]));
                    else
                    {
                        var shell = TerrainPatch.FromId(id, 2);
                        if (!File.Exists(shell.DefaultPath))
                            continue;
                        shell.Read(shell.DefaultPath);
                        id2patch.Add(id, shell);
                        result.Add((s.name, s.lat, s.lon, line, sample, shell));
                    }
                }
                return result.ToArray();
            }
        }

        #endregion

        #region Utilities

        public static List<Vector3d> GetSunVectors(StudyInterval interval) => GetSunVectors(interval.Start, interval.Stop, interval.Step);
        public static List<Vector3d> GetSunVectors(DateTime start, DateTime stop, TimeSpan step)
        {
            lock (_sunVectorCache)
            {
                if (_sunVectorCache.TryGetValue((start, stop, step), out List<Vector3d> vectors))
                    return vectors;
                vectors = new List<Vector3d>();
                for (var time = start; time <= stop; time += step)
                    vectors.Add(CSpice.SunPosition(time));
                _sunVectorCache.Add((start, stop, step), vectors);
                return vectors;
            }
        }

        public static List<Vector3d> GetEarthVectors(StudyInterval interval) => GetEarthVectors(interval.Start, interval.Stop, interval.Step);
        public static List<Vector3d> GetEarthVectors(DateTime start, DateTime stop, TimeSpan step)
        {
            lock (_earthVectorCache)
            {
                if (_earthVectorCache.TryGetValue((start, stop, step), out List<Vector3d> vectors))
                    return vectors;
                vectors = new List<Vector3d>();
                for (var time = start; time <= stop; time += step)
                    vectors.Add(CSpice.EarthPosition(time));
                _earthVectorCache.Add((start, stop, step), vectors);
                return vectors;
            }
        }

        #endregion
    }

    public class LongestNightsForPatch
    {
        public int Line { get; set; }
        public int Sample { get; set; }
        public int BandCount { get; set; }
        public float[] BandThresholds { get; set; }
        public float[][,] Bands { get; set; }

        public LongestNightsForPatch() { }
        public static LongestNightsForPatch Read(string path)
        {
            using (var br = new BinaryReader(File.Open(path, FileMode.Open)))
            {
                var r = new LongestNightsForPatch();
                r.Line = br.ReadInt32();
                r.Sample = br.ReadInt32();
                var height = br.ReadInt32();
                var width = br.ReadInt32();
                Debug.Assert(height == TerrainPatch.DefaultSize && width == TerrainPatch.DefaultSize);
                var bandcount = br.ReadInt32();
                var band_thresholds = new float[bandcount];
                for (var i = 0; i < bandcount; i++)
                    band_thresholds[i] = br.ReadSingle();
                var bands = new float[bandcount][,];
                for (var i=0;i<bandcount;i++)
                {
                    var band = new float[height, width];
                    for (var row = 0; row < height; row++)
                        for (var col = 0; col < width; col++)
                            band[row, col] = br.ReadSingle();
                    bands[i] = band;
                }
                r.BandCount = bandcount;
                r.BandThresholds = band_thresholds;
                r.Bands = bands;
                return r;
            }
        }
    }
}
