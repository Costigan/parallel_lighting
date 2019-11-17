using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using viper.corelib.patch;
using viper.corelib.utilities;

namespace viper.corelib.terrain
{
    /// <summary>
    /// MeshGenerator can create and write a variable-resolution .ply file corresponding to a portion of a LOLA terrain
    /// To generate a mesh
    /// 1. Call Initialize()
    /// 2. Call AllocateMeshes()
    /// </summary>
    public class MeshGenerator
    {
        public const int Size = TerrainPatch.DefaultSize;
        public const int IdBoundMinusOne = TerrainPatch.DEM_size / TerrainPatch.DefaultSize - 1;

        public InMemoryInt16Terrain Terrain;

        public Dictionary<Point, MeshPatch> Allocated = new Dictionary<Point, MeshPatch>();
        public Dictionary<float, int> DistanceToZoom { get; set; } = new Dictionary<float, int> { { 3f, 1 }, { 6f, 2 }, { 12f, 8 }, { 16f, 32 }, { 1000f, 64 } };  // distance units are meshes
        public Dictionary<int, Color> ZoomToColor { get; set; } = new Dictionary<int, Color> { { 1, Color.White }, { 2, Color.LightGreen }, { 4, Color.LightSkyBlue }, { 8, Color.PaleVioletRed }, { 16, Color.LightGoldenrodYellow }, { 32, Color.CornflowerBlue }, { 64, Color.LightGray }, { 128, Color.Yellow } };

        public bool WriteColorVertices { get; set; } = true;  // default
        public bool WriteASCII { get; set; } = false;

        public List<MeshPatch> Patches = new List<MeshPatch>();

        public void GenerateMesh(List<Point> highRes, Rectangle boundingBox, Dictionary<int, Color> zoomToColor = null, Dictionary<float, int> distanceToZoom = null)
        {
            GeneratePatches(boundingBox);
            AllocateZooms(highRes, distanceToZoom, zoomToColor);
            GenerateVertices();
            GenerateFaces();
        }

        // Call Allocate the normal meshes and the edges
        public void GeneratePatches(Rectangle boundingBox)
        {
            Allocated = new Dictionary<Point, MeshPatch>();

            // Allocate normal meshes
            for (var row = boundingBox.Top; row < boundingBox.Bottom; row++)
                for (var col = boundingBox.Left; col < boundingBox.Right; col++)
                {
                    var mesh = MeshPatch.FromId(col, row);
                    mesh.IsEdge = false;
                    Allocated.Add(mesh.Id, mesh);
                }

            // Allocate edge meshes
            foreach (var m in Allocated.Values.ToList())
            {
                var below = m.RelativeId(0, 1);
                var toright = m.RelativeId(1, 0);
                var diagonal = m.RelativeId(1, 1);
                if (!Allocated.ContainsKey(below))
                    Allocated.Add(below, MeshPatch.FromId(below.X, below.Y, m.Zoom, true));
                if (!Allocated.ContainsKey(toright))
                    Allocated.Add(toright, MeshPatch.FromId(toright.X, toright.Y, m.Zoom, true));
                if (!Allocated.ContainsKey(diagonal))
                    Allocated.Add(diagonal, MeshPatch.FromId(diagonal.X, diagonal.Y, m.Zoom, true));
            }
        }

        private void AllocateZooms(List<Point> highRes, Dictionary<float, int> distanceToZoom, Dictionary<int, Color> zoomToColor)
        {
            if (distanceToZoom == null) distanceToZoom = DistanceToZoom;
            var dtoz = distanceToZoom.Keys.OrderBy(k => k).Select(k => Tuple.Create(k, DistanceToZoom[k])).ToArray();
            if (zoomToColor == null) zoomToColor = ZoomToColor;
            foreach (var mesh in Allocated.Values)
            {
                mesh.Zoom = mesh.IsEdge ? MeshPatch.Size : GetZoom(dtoz, highRes.Min(p => p.DistanceTo(mesh.Id)));
                mesh.Color = zoomToColor.ContainsKey(mesh.Zoom) ? zoomToColor[mesh.Zoom] : Color.White;
            }
        }

        protected int GetZoom(Tuple<float, int>[] dtoz, float distance)
        {
            for (var i = 0; i < dtoz.Length; i++)
                if (distance <= dtoz[i].Item1)
                    return dtoz[i].Item2;
            return 1;
        }

        public int GenerateVertices()
        {
            var count = 0;
            foreach (var m in Allocated.Values.Where(m => !m.IsEdge))
                count = m.AllocateCenterVertices(count);
            foreach (var m in Allocated.Values.Where(m => m.IsEdge))
                count = m.AllocateEdgeVertices(this, count);
            return count;
        }

        public void GenerateFaces()
        {
            foreach (var m in Allocated.Values)
                m.Faces = m.GenerateFacesInternal(this);
        }

        /// <summary>
        /// Return the vertices as an enumerable of Point(sample, line)
        /// </summary>
        /// <returns></returns>
        public IEnumerable<Point> EnumerateVertices()
        {
            if (Allocated != null)
                foreach (var m in Allocated.Values)
                    foreach (var vertex_offset in m.VertexOffsets)
                        yield return new Point(m.Sample + vertex_offset.X, m.Line + vertex_offset.Y);
        }

        public IEnumerable<int[]> EnumerateFaces()
        {
            if (Allocated != null)
                foreach (var m in Allocated.Values)
                    if (m.Faces != null)
                        foreach (var f in m.Faces)
                            yield return f;
        }

        #region Writing meshes to files

        public void WritePlyFile(InMemoryInt16Terrain manager, string path, bool? writeASCII = null, bool? writeColorVertices = null)
        {
            if (Allocated == null)
                return;
            var ascii = writeASCII.HasValue ? writeASCII.Value : WriteASCII;
            var color = writeColorVertices.HasValue ? writeColorVertices.Value : WriteColorVertices;

            var vertex_count = EnumerateVertices().Count();
            var face_count = EnumerateFaces().Count();

            var file = File.Open(path, FileMode.Create);
            var sw = new StreamWriter(file);
            using (var bw = new BinaryWriter(file))
            {
                using (var ms = new MemoryStream(500))
                using (var msw = new StreamWriter(ms))
                {
                    msw.WriteLine("ply");
                    msw.WriteLine(ascii ? "format ascii 1.0" : "format binary_little_endian 1.0");

                    msw.WriteLine($"element vertex {vertex_count}");
                    msw.WriteLine("property double x");
                    msw.WriteLine("property double y");
                    msw.WriteLine("property double z");

                    if (color)
                    {
                        msw.WriteLine("property uchar red");
                        msw.WriteLine("property uchar green");
                        msw.WriteLine("property uchar blue");
                    }

                    msw.WriteLine($"element face {face_count}");
                    msw.WriteLine("property list uchar int vertex_index");

                    msw.WriteLine("end_header");
                    msw.Flush();
                    var bytes = ms.GetBuffer();
                    if (ascii)
                    {
                        for (var i = 0; i < ms.Length; i++)
                            sw.Write((char)bytes[i]);
                    }
                    else
                    {
                        bw.Write(bytes, 0, (int)ms.Length);
                    }
                }

                sw.Flush();
                bw.Flush();

                var count_vertices = 0;

                // testing
                math.Vector3d offset;
                {
                    var vo = Allocated.Values.First();
                    var off = vo.VertexOffsets[0];
                    offset = MeshGenerator.LineSampleToMoonME(manager, vo.Line + off.Y, vo.Sample+off.X);
                }

                // Write the vertices
                foreach (var m in Allocated.Values)
                {
                    var colorFromTile = m.Color;
                    foreach (var vertex_offset in m.VertexOffsets)
                    {
                        var sample = m.Sample + vertex_offset.X;
                        var line = m.Line + vertex_offset.Y;
                        var pos = MeshGenerator.LineSampleToMoonME(manager, line, sample);

                        // Testing
                        pos = pos - offset;

                        count_vertices++;
                        if (ascii)
                        {
                            sw.Write($"{pos.X} {pos.Y} {pos.Z}");
                            if (color)
                                sw.WriteLine($" {colorFromTile.R} {colorFromTile.G} {colorFromTile.B}");
                        }
                        else
                        {
                            bw.Write(pos.X);
                            bw.Write(pos.Y);
                            bw.Write(pos.Z);
                            if (color)
                            {
                                bw.Write(colorFromTile.R);
                                bw.Write(colorFromTile.G);
                                bw.Write(colorFromTile.B);
                            }
                        }
                    }
                }

                System.Diagnostics.Debug.Assert(count_vertices == vertex_count);

                // write the faces
                if (ascii)
                {
                    foreach (var face in EnumerateFaces())
                    {
                        sw.Write(face.Length);
                        for (var i = 0; i < face.Length; i++)
                        {
                            var idx = face[i];
                            System.Diagnostics.Debug.Assert(idx >= 0 && idx < vertex_count);
                            sw.Write(' ');
                            sw.Write(idx);
                        }
                        sw.WriteLine();
                    }
                }
                else
                {
                    foreach (var face in EnumerateFaces())
                    {
                        bw.Write((byte)face.Length);
                        for (var i = 0; i < face.Length; i++)
                        {
                            var idx = face[i];
                            System.Diagnostics.Debug.Assert(idx >= 0 && idx < vertex_count);
                            bw.Write(idx);
                        }
                    }
                }

                bw.Flush();
                sw.Flush();
            }
        }

        #endregion

        #region Utility methods

        protected static int Clamp(int v, int min, int max) => Math.Min(Math.Max(v, min), max);
        protected static int Clamp(int v) => Clamp(v, 0, IdBoundMinusOne);

        public static math.Vector3d LineSampleToMoonME(InMemoryInt16Terrain m, int line, int sample)
        {
            var relz = m.LineSampleToTerrainOffset(line, sample);
            var radius = InMemoryInt16Terrain.MoonRadius + relz / 1000d;
            InMemoryInt16Terrain.GetLatLon(line, sample, out double lat, out double lon);
            var z = radius * Math.Sin(lat);
            var c = radius * Math.Cos(lat);
            var x = c * Math.Cos(lon);  // TODO: Not sure about these
            var y = c * Math.Sin(lon);
            return new math.Vector3d(x, y, z);
        }

        #endregion

        #region Testing Methods

        public bool EveryFacePointsUp()
        {
            var up = new math.Vector3d(0, 0, 1d);
            var vertices = EnumerateVertices().ToArray();
            foreach (var face in EnumerateFaces())
            {
                var tri = FaceToTriangle(face, vertices);
                System.Diagnostics.Debug.Assert(tri.Length == 3);
                var v1 = tri[1] - tri[0];
                var v2 = tri[2] - tri[0];
                var cross = math.Vector3d.Cross(v1, v2).Normalized();
                var dot = math.Vector3d.Dot(up, cross);
                if (dot <= 0d)
                    return false;
            }
            return true;
        }

        math.Vector3d[] FaceToTriangle(int[] face, Point[] vertices) => face.Select(v =>
        {
            var pt = vertices[v];
            return LineSampleToMoonME(Terrain, pt.Y, pt.X);
        }).ToArray();

        #endregion
    }
}
