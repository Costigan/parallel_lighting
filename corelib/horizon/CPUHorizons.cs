using viper.corelib.patch;
using viper.corelib.terrain;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using viper.corelib.utilities;
using viper.corelib.interfaces;

namespace viper.corelib.horizon
{
    public class CPUHorizons
    {
        public int? ObserverHeightInMeters;  // If this has a value, then override what's already in the patch

        public string ShadowCalculationQueueFilename = @"shadow_calculation_queue.json";
        public int MaxDegreeOfParallelism = GetParallelism();
        public int MaxSpread = 177; // 70;
        public bool FilterTentpoles = true;
        public bool NearHorizonsOnly = false;
        public IStatusReceiver StatusReceiver = null;
        //public views.MapView MapView;

        public InMemoryInt16Terrain Terrain;

        protected static int GetParallelism() => Math.Max(Environment.ProcessorCount - 2, 1);

        protected string _ShadowCalculationQueuePath;
        public string ShadowCalculationQueuePath
        {
            get
            {
                if (_ShadowCalculationQueuePath==null)
                    _ShadowCalculationQueuePath = Path.Combine(ViperEnvironment.MapRoot, ShadowCalculationQueueFilename);
                return _ShadowCalculationQueuePath;
            }
            set
            {
                _ShadowCalculationQueuePath = value;
            }
        }

        public virtual void RunQueue(List<TerrainPatch> queue, Action<List<TerrainPatch>> queue_saver = null, int spread = -1, bool overHorizonCheck = true, bool centerOnly = false, bool center = true, bool unloadHorizons = true, bool writeHorizons = true)
        {
            if (StatusReceiver != null) StatusReceiver.IsRunning = true;
            Debug.Assert(ViperEnvironment.Terrain != null);
            centerOnly |= NearHorizonsOnly;
            if (queue == null)
                queue = ReadShadowCalculationQueue();
            if (spread < 0) spread = MaxSpread;
            while (queue.Count > 0)
            {
                var patch = queue[0];
                Console.WriteLine($"Starting [{patch.Line},{patch.Sample}] ...");
                try
                {
                    var stopwatch = Stopwatch.StartNew();

                    // Override observer height if specified in generator
                    if (ObserverHeightInMeters.HasValue)
                        patch.ObserverHeightInMeters = ObserverHeightInMeters.Value;

                    var far_field = new List<TerrainPatch>();
                    if (StatusReceiver != null) StatusReceiver.ProcessingPatches = new List<TerrainPatch> { patch };
                    patch.FillPointsAndMatrices(Terrain);
                    if (center)
                        patch.ParallelAddNearHorizon(Terrain);
                    Console.WriteLine($"  near horizon calculated or loaded in {stopwatch.Elapsed}");

                    if (true)
                    {
                        var options = new ParallelOptions { MaxDegreeOfParallelism = MaxDegreeOfParallelism };
                        for (var i = 1; i < spread; i++)
                        {
                            var other1 = patch.SurroundingPatches(i).ToList();
                            Parallel.ForEach(other1, options, o => o.FillPoints(Terrain));
                            var other2 = overHorizonCheck ? other1.Where(patch.IsOverHorizon).ToList() : other1;
                            if (StatusReceiver != null)
                            {
                                far_field.AddRange(other2);
                                StatusReceiver.FarPatches = far_field;
                            }
                            Parallel.ForEach(other2, options, o => patch.UpdateHorizon(o));
                        }
                    }
                    if (writeHorizons)
                        patch.Write();
                    if (unloadHorizons)
                        patch.InitializeHorizons();  // Unload the horizon data (100MB)
                    stopwatch.Stop();
                    var seconds_per_patch = far_field.Count == 0 ? 0f : (stopwatch.ElapsedMilliseconds / 1000f) / far_field.Count;
                    Console.WriteLine($"  Finished [{patch.Line},{patch.Sample}] time={stopwatch.Elapsed}.  sec/patch={seconds_per_patch}");

                    // Update queue on disk
                    queue.RemoveAt(0);
                    queue_saver?.Invoke(queue);
                }
                catch (Exception e1)
                {
                    Console.WriteLine(e1);
                    Console.WriteLine(e1.StackTrace);
                }
                Console.WriteLine(@"Finished queue with CPU.");
                if (StatusReceiver != null) StatusReceiver.IsRunning = false;
            }
        }

        public void RunQueue()
        {
            RunQueue(ReadShadowCalculationQueue(), q => WriteShadowCalculationQueue(q));
        }

        #region Helper Functions

        // Right now, this queue is independent of oberver height.
        public List<TerrainPatch> ReadShadowCalculationQueue() => ReadShadowCalculationQueue(ShadowCalculationQueuePath);

        public static List<TerrainPatch> ReadShadowCalculationQueue(string path)
        {
            if (!File.Exists(path)) return new List<TerrainPatch>();
            using (var sr = new StreamReader(path))
            using (var jr = new JsonTextReader(sr))
            {
                var serializer = new JsonSerializer();
                return serializer.Deserialize<List<TerrainPatch>>(jr);
            }
        }

        public void WriteShadowCalculationQueue(List<TerrainPatch> q) => WriteShadowCalculationQueue(q, ShadowCalculationQueuePath);

        public static void WriteShadowCalculationQueue(List<TerrainPatch> q, string path)
        {
            using (var sr = new StreamWriter(path))
            using (var jr = new JsonTextWriter(sr))
            {
                var serializer = new JsonSerializer();
                serializer.Serialize(jr, q);
            }
        }

        #endregion

        #region Extend

        public void ExtendExistingPatches()
        {
            if (MaxSpread != 177)
                throw new Exception($"Unexpected value of MaxSpread={MaxSpread}");
            var filenames = (new DirectoryInfo(ViperEnvironment.HorizonsRoot)).EnumerateFiles("*.bin").Select(fi => fi.FullName).ToList();
            foreach (var filename in filenames)
            {
                var patch = TerrainPatch.ReadFrom(filename);
                var stopwatch = new Stopwatch();
                stopwatch.Start();
                Console.WriteLine($"Starting [{patch.Line},{patch.Sample}] ...");
                patch.FillPointsAndMatrices(ViperEnvironment.Terrain);
                try
                {
                    var far_field = new List<TerrainPatch>();
                    if (StatusReceiver != null) StatusReceiver.ProcessingPatches = new List<TerrainPatch> { patch };
                    var options = new ParallelOptions { MaxDegreeOfParallelism = MaxDegreeOfParallelism };
                    for (var i = 1; i < MaxSpread; i++)
                    {
                        var other1 = patch.SurroundingPatches(i).Where(p => !patch.ShadowCasters.Contains(p.Id)).ToList();
                        Parallel.ForEach(other1, options, o => o.FillPoints(ViperEnvironment.Terrain));
                        var other2 = other1.Where(patch.IsOverHorizon).ToList();
                        if (StatusReceiver != null)
                        {
                            far_field.AddRange(other2);
                            StatusReceiver.FarPatches = far_field;
                        }
                        Parallel.ForEach(other2, options, o => patch.UpdateHorizon(o));
                    }
                    patch.Write();
                    patch.InitializeHorizons();  // Unload the horizon data (100MB)
                    stopwatch.Stop();
                    var seconds_per_patch = far_field.Count == 0 ? 0f : (stopwatch.ElapsedMilliseconds / 1000f) / far_field.Count;
                    Console.WriteLine($"  Finished [{patch.Line},{patch.Sample}] time={stopwatch.Elapsed}.  sec/patch={seconds_per_patch}");
                }
                catch (Exception e1)
                {
                    Console.WriteLine(e1);
                    Console.WriteLine(e1.StackTrace);
                }
            }
        }

        #endregion
    }
}
