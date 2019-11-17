using System;
using System.Collections.Generic;
using System.Drawing;
using viper.corelib.patch;

namespace viper.corelib.terrain
{
    /// <summary>
    /// MeshPatch represents the portion of a ply mesh that corresponds to a single terrain patch
    /// </summary>
    public class MeshPatch
    {
        public const int IdBoundMinusOne = TerrainPatch.DEM_size / TerrainPatch.DefaultSize - 1;
        public const int Size = TerrainPatch.DefaultSize;

        public Point Id;
        public int Zoom = 1;  // Zoom elt {1,2,4,8,16,32,64,128}
        public bool IsEdge;
        public Color Color = Color.White;

        public int FirstVertexIndex;
        public Point[] VertexOffsets;
        public int VertexCount => VertexOffsets.Length;  // Point offsets for each vertex within this mesh

        // These contain absolute vertex indices
        public int[,] OffsetToVertex;

        public int[][] Faces;
        public int FaceCount => Faces == null ? 0 : Faces.Length;

        public int Line => Id.Y * Size;
        public int Sample => Id.X * Size;

        public static int Clamp(int v, int min, int max) => Math.Min(Math.Max(v, min), max);
        public static int Clamp(int v) => Clamp(v, 0, IdBoundMinusOne);

        public Point RelativeId(int x, int y) => new Point(Id.X + x, Id.Y + y);
        public MeshPatch RelativeMesh(int x, int y) => new MeshPatch { Id = new Point(Clamp(Id.X + x), Clamp(Id.Y + y)) };

        public static MeshPatch FromId(int x, int y, int zoom = 1, bool isEdge = false) => new MeshPatch { Id = new Point(x, y), Zoom = zoom, IsEdge = isEdge };

        /// <summary>
        /// Allocate vertices for a full patch (with some zoom)
        /// </summary>
        /// <param name="firstVertexIndex"></param>
        /// <returns></returns>
        internal int AllocateCenterVertices(int firstVertexIndex)
        {
            System.Diagnostics.Debug.Assert(IsEdge == false);
            FirstVertexIndex = firstVertexIndex;
            var countPerSide = Size / Zoom;
            var vcount = countPerSide * countPerSide;
            VertexOffsets = new Point[vcount];
            var ptr = 0;
            for (var row = 0; row < Size; row += Zoom)
                for (var col = 0; col < Size; col += Zoom)
                    VertexOffsets[ptr++] = new Point(col, row);  // These are local offsets

            OffsetToVertex = new int[Size, Size];
            for (var row = 0; row < Size; row++)
                for (var col = 0; col < Size; col++)
                    OffsetToVertex[col, row] = -1;

            for (var i = 0; i < VertexOffsets.Length; i++)
            {
                var offset = VertexOffsets[i];
                OffsetToVertex[offset.X, offset.Y] = i + firstVertexIndex;
            }
            return FirstVertexIndex + VertexOffsets.Length;
        }

        /// <summary>
        /// Allocate vertices for an edge patch
        /// Call this after all normal patches have been allocated
        /// </summary>
        /// <param name="firstVertexIndex"></param>
        /// <returns></returns>
        internal int AllocateEdgeVertices(MeshGenerator g, int firstVertexIndex)
        {
            Zoom = Size;
            IsEdge = true;
            FirstVertexIndex = firstVertexIndex;
            var vertices = new List<Point>();
            // The upper left must be used if this is called in the first place
            vertices.Add(new Point(0, 0));
            // If there's a patch to the left, then fill the first column, stepping according to that patch's zoom
            if (g.Allocated.TryGetValue(RelativeId(-1, 0), out MeshPatch left_patch) && !left_patch.IsEdge)
            {
                var zoom = left_patch.Zoom;
                for (var i = zoom; i < Size; i += zoom)
                    vertices.Add(new Point(0, i));
                Zoom = Math.Min(Zoom, zoom);
            }
            // If there's a patch above, then fill the first column, stepping according to that patch's zoom
            if (g.Allocated.TryGetValue(RelativeId(0, -1), out MeshPatch above_patch) && !above_patch.IsEdge)
            {
                var zoom = above_patch.Zoom;
                for (var i = zoom; i < Size; i += zoom)
                    vertices.Add(new Point(i, 0));
                Zoom = Math.Min(Zoom, zoom);
            }

            if (Zoom == Size)
            {
                System.Diagnostics.Debug.Assert(vertices.Count == 1);  // Must be a diagonal.  Zoom shouldn't matter;
            }
            //System.Diagnostics.Debug.Assert(Zoom != Size);  // At least one of the two cases above must have been true

            VertexOffsets = vertices.ToArray();
            OffsetToVertex = new int[Size, Size];
            for (var row = 0; row < Size; row++)
                for (var col = 0; col < Size; col++)
                    OffsetToVertex[col, row] = -1;
            for (var i = 0; i < VertexOffsets.Length; i++)
            {
                var offset = VertexOffsets[i];
                OffsetToVertex[offset.X, offset.Y] = i + firstVertexIndex;
            }

            return FirstVertexIndex + VertexOffsets.Length;
        }

        int VertexIndex(int col, int row)
        {
            var index = OffsetToVertex[col, row];
            System.Diagnostics.Debug.Assert(index >= 0);
            return index;
        }

        bool IsVertexIndex(int col, int row) => OffsetToVertex[col, row] >= 0;

        public static bool DrawCenter = true, DrawDiagonal = true, DrawRight1 = true, DrawRight2 = true, DrawRight3 = true, DrawDown1 = true, DrawDown2 = true, DrawDown3 = true;
        internal int[][] GenerateFacesInternal(MeshGenerator g)
        {
            if (IsEdge)
            {
                return new int[0][];
            }
            // Normal case
            {
                var below = g.Allocated[RelativeId(0, 1)];
                var right = g.Allocated[RelativeId(1, 0)];
                var diagonal = g.Allocated[RelativeId(1, 1)];
                var edges = new List<int[]>();

                // Simple case
                for (var row = 0; row < Size; row += Zoom)
                    for (var col = 0; col < Size; col += Zoom)
                    {
                        var v1 = VertexIndex(col, row);
                        if (col + Zoom == Size)
                        {   // to right
                            if (row + Zoom == Size)
                            {   // diagonal
                                if (DrawDiagonal)
                                {
                                    if (false) // (Zoom == right.Zoom && Zoom == below.Zoom)
                                    {
                                        var right_row = right.Zoom * (row / right.Zoom);
                                        var v2 = right.VertexIndex(0, right_row);
                                        var below_col = below.Zoom * (col / below.Zoom);
                                        var v3 = below.VertexIndex(below_col, 0);
                                        var v4 = diagonal.VertexIndex(0, 0);
                                        edges.Add(new[] { v1, v3, v2 });
                                        edges.Add(new[] { v2, v3, v4 });
                                    }
                                    else
                                    {
                                        // Compare right
                                        if (Zoom == right.Zoom)
                                            edges.Add(new[] { v1, diagonal.VertexIndex(0, 0), right.VertexIndex(0, Size - right.Zoom) });
                                        else if (Zoom < right.Zoom && below.IsVertexIndex(col, 0))
                                        {
                                            var v2 = below.VertexIndex(col, 0);
                                            var v3 = diagonal.VertexIndex(0, 0);
                                            edges.Add(new[] { v1, v2, v3 });
                                        }
                                        else if (Zoom > right.Zoom)
                                        {
                                            for (var right_row = row; right_row < Size - right.Zoom; right_row += right.Zoom)
                                                edges.Add(new[] { v1, right.VertexIndex(0, right_row + right.Zoom), right.VertexIndex(0, right_row) });
                                            edges.Add(new[] { v1, diagonal.VertexIndex(0, 0), right.VertexIndex(0, Size - right.Zoom) });
                                        }

                                        // Compare below
                                        if (Zoom == below.Zoom)
                                            edges.Add(new[] { v1, below.VertexIndex(Size - below.Zoom, 0), diagonal.VertexIndex(0, 0) });
                                        else if (Zoom < below.Zoom && right.IsVertexIndex(0, row))
                                        {
                                            var v2 = right.VertexIndex(0, row);
                                            var v3 = diagonal.VertexIndex(0, 0);
                                            edges.Add(new[] { v1, v3, v2 });
                                        }
                                        else if (Zoom > below.Zoom)
                                        {
                                            for (var below_col = col; below_col < Size - below.Zoom; below_col += below.Zoom)
                                                edges.Add(new[] { v1, below.VertexIndex(below_col, 0), below.VertexIndex(below_col + below.Zoom, 0) });
                                            edges.Add(new[] { v1, below.VertexIndex(Size - below.Zoom, 0), diagonal.VertexIndex(0, 0) });
                                        }
                                    }
                                }
                            }
                            else
                            {   // to right only, same zoom
                                if (Zoom == right.Zoom)
                                {
                                    if (DrawRight1)
                                    {
                                        var v2 = right.VertexIndex(0, row);
                                        var v3 = VertexIndex(col, row + Zoom);
                                        var v4 = right.VertexIndex(0, row + Zoom);
                                        edges.Add(new[] { v1, v3, v2 });
                                        edges.Add(new[] { v2, v3, v4 });
                                    }
                                }
                                else if (Zoom > right.Zoom)
                                {   // right only, spacing getting tighter
                                    if (DrawRight2)
                                    {
                                        var v2 = VertexIndex(col, row + Zoom);
                                        var other_vertices = new List<int>(Zoom / right.Zoom + 1);
                                        for (var i = row; i < row + Zoom; i += right.Zoom)
                                            other_vertices.Add(right.VertexIndex(0, i));
                                        other_vertices.Add(row + Zoom == Size ? diagonal.VertexIndex(0, 0) : right.VertexIndex(0, row + Zoom));
                                        for (var i = 0; i < other_vertices.Count - 1; i++)
                                            edges.Add(new[] { v1, other_vertices[i + 1], other_vertices[i] });
                                        edges.Add(new[] { v1, v2, other_vertices[other_vertices.Count - 1] });
                                    }
                                }
                                else if (Zoom < right.Zoom)
                                {   // right only, spacing getting looser
                                    if (DrawRight3)
                                    {
                                        var v2 = VertexIndex(col, row + Zoom);
                                        var other_row = right.Zoom * (row / right.Zoom);
                                        if (row == other_row)
                                        {
                                            var v3 = right.VertexIndex(0, row);
                                            var v4 = other_row + right.Zoom == Size ? diagonal.VertexIndex(0, 0) : right.VertexIndex(0, other_row + right.Zoom);
                                            edges.Add(new[] { v1, v4, v3 });
                                            edges.Add(new[] { v1, v2, v4 });
                                        }
                                        else
                                        {
                                            var v3 = other_row + right.Zoom == Size ? diagonal.VertexIndex(0, 0) : right.VertexIndex(0, other_row + right.Zoom);
                                            edges.Add(new[] { v1, v2, v3 });
                                        }
                                    }
                                }
                                else
                                    System.Diagnostics.Debug.Assert(false, $"invalid pair of zooms: {Zoom} vs {right.Zoom}");
                            }
                        }
                        else
                        {   // Below
                            if (row + Zoom == Size)
                            {   // below only, same zoom
                                if (Zoom == below.Zoom)
                                {
                                    if (DrawDown1)
                                    {
                                        var v2 = VertexIndex(col + Zoom, row);
                                        var v3 = below.VertexIndex(col, 0);
                                        var v4 = below.VertexIndex(col + Zoom, 0);
                                        edges.Add(new[] { v1, v3, v2 });
                                        edges.Add(new[] { v2, v3, v4 });
                                    }
                                }
                                else if (Zoom > below.Zoom)
                                {   // below only, spacing getting tighter
                                    if (DrawDown2)
                                    {
                                        var v2 = VertexIndex(col + Zoom, row);
                                        var other_vertices = new List<int>(Zoom / below.Zoom + 1);
                                        for (var i = col; i < col + Zoom; i += below.Zoom)
                                            other_vertices.Add(below.VertexIndex(i, 0));
                                        other_vertices.Add(col + Zoom == Size ? diagonal.VertexIndex(0, 0) : below.VertexIndex(col + Zoom, 0));
                                        for (var i = 0; i < other_vertices.Count - 1; i++)
                                            edges.Add(new[] { v1, other_vertices[i], other_vertices[i + 1] });
                                        edges.Add(new[] { v1, other_vertices[other_vertices.Count - 1], v2 });
                                    }
                                }
                                else if (Zoom < below.Zoom)
                                {   // below only, spacing getting looser
                                    if (DrawDown3)
                                    {
                                        var v2 = VertexIndex(col + Zoom, row);
                                        var other_col = below.Zoom * (col / below.Zoom);
                                        if (col == other_col)
                                        {
                                            var v3 = below.VertexIndex(col, 0);
                                            var v4 = other_col + below.Zoom == Size ? diagonal.VertexIndex(0, 0) : below.VertexIndex(other_col + below.Zoom, 0);
                                            edges.Add(new[] { v1, v3, v4 });
                                            edges.Add(new[] { v1, v4, v2 });
                                        }
                                        else
                                        {
                                            var v3 = other_col + below.Zoom == Size ? diagonal.VertexIndex(0, 0) : below.VertexIndex(other_col + below.Zoom, 0);
                                            edges.Add(new[] { v1, v3, v2 });
                                        }
                                    }
                                }
                                else
                                    System.Diagnostics.Debug.Assert(false, $"invalid pair of zooms: {Zoom} vs {below.Zoom}");
                            }
                            else
                            {   // inside center patch
                                if (DrawCenter)
                                {
                                    var v2 = VertexIndex(col + Zoom, row);
                                    var v3 = VertexIndex(col, row + Zoom);
                                    var v4 = VertexIndex(col + Zoom, row + Zoom);
                                    edges.Add(new[] { v1, v3, v2 });
                                    edges.Add(new[] { v2, v3, v4 });
                                }
                            }
                        }
                    }
                return edges.ToArray();
            }
        }
    }
}
