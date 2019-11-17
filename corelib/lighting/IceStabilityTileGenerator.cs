using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using viper.corelib.spice;
using viper.corelib.terrain;
using viper.corelib.utilities;

namespace viper.corelib.lighting
{
    public class IceStabilityTileGenerator
    {
        const int DefaultChunkSize = 1024;

        const string NP_Siegler_240m = @"StaticFiles/NP_Siegler_etal_2016_Ice.txt";
        const string SP_Siegler_240m = @"StaticFiles/SP_Siegler_etal_2016_Ice.txt";

        string NP_Siegler_HermiteA_20m = AppConfiguration.Get("NP_Siegler_HermiteA_20m", @"C:\RP\maps\inputs_20m\HermiteA\Hermite\Other\HermiteA6000_Siegler.txt");

        List<StabilityChunk> Chunks;
        readonly Brush[] brushes = { Brushes.Red, Brushes.Green, Brushes.Yellow, Brushes.DarkGray };
        //Brush[] _depthBrushes = Enumerable.Range(0, 256).Select(i => new SolidBrush(Color.FromArgb(i, i, i))).ToArray();

        public int ChunkSize = DefaultChunkSize;

        #region Load data

        public void Load()
        {
            Console.WriteLine(@"IceStabilityTileGenerator loading ...");
            var path = ViperEnvironment.IsNorth ? NP_Siegler_240m : SP_Siegler_240m;
            var triangles_240 = GetTriangles240m(path);
            //var triangles_20m = GetTriangles240m(NP_Siegler_HermiteA_20m);
            //triangles_240.AddRange(triangles_20m);
            Chunks = BuildChunks(triangles_240);
            Console.WriteLine(@"IceStabilityTileGenerator finished loading.");
        }

        private List<StabilityChunk> BuildChunks(List<IceStabilityTriangle> triangles)
        {
            const int DEMSize = ViperEnvironment.TerrainWidth20m;

            // Build 2d array of chunks
            var chunkwidth = (int)Math.Ceiling(DEMSize / (float)ChunkSize);
            var chunkArray = new StabilityChunk[chunkwidth, chunkwidth];
            for (var y = 0; y < chunkwidth; y++)
                for (var x = 0; x < chunkwidth; x++)
                    chunkArray[x, y] = new StabilityChunk { BoundingBox = new Rectangle(x * ChunkSize, y*ChunkSize, ChunkSize, ChunkSize), Triangles = new List<IceStabilityTriangle>() };

            // loop through the triangles, assigning them to chunks
            foreach (var t in triangles)
            {
                var box = t.BoundingBox();
                InsertIntoChunks(chunkArray, box.X / ChunkSize, box.Y / ChunkSize, t, box);
            }

            var chunks = new List<StabilityChunk>();
            for (var y = 0; y < chunkwidth; y++)
                for (var x = 0; x < chunkwidth; x++)
                    if (chunkArray[x, y].Triangles.Count > 0)
                        chunks.Add(chunkArray[x, y]);
            return chunks;
        }

        void InsertIntoChunks(StabilityChunk[,] ary, int x, int y, IceStabilityTriangle tri, Rectangle box)
        {
            InsertIntoChunkMaybe(ary, x, y, tri, box);
            InsertIntoChunkMaybe(ary, x - 1, y, tri, box);
            InsertIntoChunkMaybe(ary, x + 1, y, tri, box);

            InsertIntoChunkMaybe(ary, x, y-1, tri, box);
            InsertIntoChunkMaybe(ary, x - 1, y - 1, tri, box);
            InsertIntoChunkMaybe(ary, x + 1, y - 1, tri, box);

            InsertIntoChunkMaybe(ary, x, y + 1, tri, box);
            InsertIntoChunkMaybe(ary, x - 1, y + 1, tri, box);
            InsertIntoChunkMaybe(ary, x + 1, y + 1, tri, box);
        }

        void InsertIntoChunkMaybe(StabilityChunk[,] ary, int x, int y, IceStabilityTriangle tri, Rectangle box)
        {
            if (x>=0 && y>=0 && x<ary.GetLength(0) && y<ary.GetLength(1))
            {
                var chunk = ary[x, y];
                if (box.IntersectsWith(chunk.BoundingBox))
                    chunk.Triangles.Add(tri);
            }            
        }

        List<IceStabilityTriangle> GetTriangles240m(string path)
        {
            var triangles = new List<IceStabilityTriangle>();
            string line;
            using (var file = new System.IO.StreamReader(path))
                while ((line = file.ReadLine()) != null)
                {
                    var tokens = line.Split('\t');
                    var column = 0;
                    var x1 = double.Parse(tokens[column++]);
                    var y1 = double.Parse(tokens[column++]);
                    var z1 = double.Parse(tokens[column++]);
                    var x2 = double.Parse(tokens[column++]);
                    var y2 = double.Parse(tokens[column++]);
                    var z2 = double.Parse(tokens[column++]);
                    var x3 = double.Parse(tokens[column++]);
                    var y3 = double.Parse(tokens[column++]);
                    var z3 = double.Parse(tokens[column++]);
                    var lon = double.Parse(tokens[column++]);
                    var lat = double.Parse(tokens[column++]);
                    var elevation = double.Parse(tokens[column++]);
                    var num = (int)Math.Round(float.Parse(tokens[column++]));
                    var modern = float.Parse(tokens[column++]);
                    var paleo = float.Parse(tokens[column++]);

                    var p1 = new PointF();
                    var p2 = new PointF();
                    var p3 = new PointF();

                    ProjectPoint(new double[] { x1, y1, z1 }, ref p1);
                    ProjectPoint(new double[] { x2, y2, z2 }, ref p2);
                    ProjectPoint(new double[] { x3, y3, z3 }, ref p3);

                    var polygon = new PointF[] { p1, p2, p3 };

                    var m = modern;
                    var zone = m < 0 ? 0
                             : m < 0.5 ? 1
                             : m < 1 ? 2
                             : 3;
                    triangles.Add(new IceStabilityTriangle { Depth = modern, Zone = zone, Polygon = polygon });
                }
            return triangles;
        }

        private void ProjectPoint(double[] rectan, ref PointF p)
        {
            double r = 0, lat_rad = 0, lon_rad = 0;
            CSpice.reclat_c(rectan, ref r, ref lon_rad, ref lat_rad);
            InMemoryInt16Terrain.GetLineSample(lat_rad, lon_rad, out int line, out int sample);
            p.X = sample;
            p.Y = line;
        }

        #endregion

        #region structs

        public struct IceStabilityTriangle
        {
            public PointF[] Polygon;
            public float Depth;
            public int Zone;

            public Rectangle BoundingBox()
            {
                int minx = int.MaxValue, maxx = int.MinValue, miny = int.MaxValue, maxy = int.MinValue;
                if (Polygon!=null)
                for (var i=0;i< Polygon.Length;i++)
                    {
                        var p = Polygon[i];
                        minx = Math.Min(minx, (int)Math.Floor(p.X));
                        maxx = Math.Max(maxx, (int)Math.Ceiling(p.X));
                        miny = Math.Min(miny, (int)Math.Floor(p.Y));
                        maxy = Math.Max(maxy, (int)Math.Ceiling(p.Y));
                    }
                return new Rectangle(minx-1, miny-1, maxx - minx + 2, maxy - miny + 2);  // Some margin
            }

            public PointF AveragePoint()
            {
                System.Diagnostics.Debug.Assert(Polygon.Length > 0);
                var x = 0f;
                var y = 0f;
                for (var i=0;i<Polygon.Length;i++)
                {
                    x += Polygon[i].X;
                    y += Polygon[i].Y;
                }
                return new PointF(x / Polygon.Length, y / Polygon.Length);
            }
        }

        public struct StabilityChunk
        {
            public Rectangle BoundingBox;
            public List<IceStabilityTriangle> Triangles;
        }

        #endregion

        #region drawing

        // Note no graphics scaling here.  That is in the caller.
        public void Draw(Graphics g, Rectangle box)
        {
            //g.FillRectangle(Brushes.Green, 0, 0, 32000, 32000);
            int count = 0;
            foreach (var chunk in Chunks)
                if (box.IntersectsWith(chunk.BoundingBox))
                {
                    var tris = chunk.Triangles;
                    for (var i = 0; i < tris.Count; i++) // ***** Bug Hunt i++
                    {
                        var t = tris[i];
                        var brush = brushes[t.Zone];
                        //var d = t.Depth;
                        //var index = d < 0f ? 0 : (int)(255f * d / 2.5f);
                        //var brush = _depthBrushes[index];
                        g.FillPolygon(brush, t.Polygon);
                        count++;
                    }
                }
            Console.WriteLine($"Drawing {box} count={count}");
        }

        #endregion

        public void FillArray(float[,] ary, Rectangle box)
        {
            System.Diagnostics.Debug.Assert(ary.GetLength(0) == box.Height && ary.GetLength(1) == box.Width);
            var local_chunks = Chunks.Where(c => box.IntersectsWith(c.BoundingBox)).ToList();
            FillArray(ary, box, box, local_chunks);
        }

        public void FillArray(float[,] ary, Rectangle outer, Rectangle inner, List<StabilityChunk> chunks)
        {
            //Console.WriteLine(inner.ToString());
            if (inner.Height > 64)
            {
                //  i1, i3
                //  i2, i4
                var x2 = (inner.Left + inner.Right) / 2;
                var y2 = (inner.Top + inner.Bottom) / 2;
                var i1 = new Rectangle(inner.X, inner.Y, x2 - inner.X, y2 - inner.Y);
                var i2 = new Rectangle(inner.X, y2, x2 - inner.X, inner.Bottom - y2);
                var i3 = new Rectangle(x2, inner.Y, inner.Right - x2, y2 - inner.Y);
                var i4 = new Rectangle(x2, y2, inner.Right - x2, inner.Bottom - y2);
                FillArray(ary, outer, i1, chunks.Where(c => i1.IntersectsWith(c.BoundingBox)).ToList());
                FillArray(ary, outer, i2, chunks.Where(c => i2.IntersectsWith(c.BoundingBox)).ToList());
                FillArray(ary, outer, i3, chunks.Where(c => i3.IntersectsWith(c.BoundingBox)).ToList());
                FillArray(ary, outer, i4, chunks.Where(c => i4.IntersectsWith(c.BoundingBox)).ToList());
                return;
            }

            // Calculate the triangle centers for this set
            var centers = new List<(PointF, float)>();
            for (var i = 0; i < chunks.Count; i++)
            {
                var tris = chunks[i].Triangles;
                for (var j = 0; j < tris.Count; j++)
                {
                    var tri = tris[j];
                    centers.Add((tri.AveragePoint(), tri.Depth));
                }
            }

            for (var row=inner.Y;row<inner.Bottom; row++)
                for (var col=inner.X;col<inner.Right;col++)
                {
                    var best_distance = float.MaxValue;
                    var best_depth = 1.5f;
                    var pt = new PointF(col, row);
                    for (var i=0;i<centers.Count;i++)
                    {
                        var center = centers[i];
                        var d = pt.DistanceTo(center.Item1);
                        if (d < best_distance)
                        {
                            best_distance = d;
                            best_depth = center.Item2;
                        }
                    }
                    ary[row - outer.Y, col - outer.X] = best_depth;
                }
        }
    }
}
