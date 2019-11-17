using System;
using System.Collections.Generic;
using System.Drawing;

namespace viper.corelib.math
{
    public struct Triangle
    {
        public PointF A;
        public PointF B;
        public PointF C;

        public Triangle(PointF a, PointF b, PointF c)
        {
            A = a;
            B = b;
            C = c;
        }

        //public float Area => 0.5f * (A.X * B.Y - A.X * C.Y + B.X * C.Y - B.X * A.Y + C.X * A.Y - C.X * B.Y);
        public float Area => Math.Abs(0.5f * (A.X * (B.Y - C.Y) + B.X * (C.Y - A.Y) + C.X * (A.Y - B.Y)));

        public bool IsValid => !A.Equals(B) && !A.Equals(C) && !B.Equals(C);

        public bool Equivalent(Triangle o)
        {
            if (!IsValid || !o.IsValid) return false;
            if (A.Equals(o.A) && (B.Equals(o.B) && C.Equals(o.C) || B.Equals(o.C) && C.Equals(o.B))) return true;
            if (A.Equals(o.B) && (B.Equals(o.A) && C.Equals(o.C) || B.Equals(o.C) && C.Equals(o.A))) return true;
            if (A.Equals(o.C) && (B.Equals(o.A) && C.Equals(o.B) || B.Equals(o.B) && C.Equals(o.A))) return true;
            return false;
        }

        public override string ToString()
        {
            return string.Format(@"<triangle [{0},{1}], [{2},{3}], [{4},{5}]", A.X, A.Y, B.X, B.Y, C.X, C.Y);
        }

        public PointF[] Points => new PointF[] { A, B, C };
    }

    public class LabeledTriangle
    {
        public string Label;
        public Triangle Triangle;

        public override string ToString()
        {
            return string.Format(@"<LabeledTriangle {0} => {1}", Label, Triangle);
        }
    }

    public class TriangleClipper
    {
        public List<LabeledTriangle> Records = new List<LabeledTriangle>();

        public void AddLight(Triangle t, RectangleF r)
        {
            ClipLeft(t, r);
        }

        public void ClipLeft(Triangle t, RectangleF r)
        {
            //Record("ClipLeft", t);
            var outside = new PointF[3];
            var inside = new PointF[3];
            var outsidei = 0;
            var insidei = 0;

            // Split the triangle's points to the left and right sides of the square's left side
            var left = r.Left;
            if (t.A.X > left)
                inside[insidei++] = t.A;
            else
                outside[outsidei++] = t.A;
            if (t.B.X > left)
                inside[insidei++] = t.B;
            else
                outside[outsidei++] = t.B;
            if (t.C.X > left)
                inside[insidei++] = t.C;
            else
                outside[outsidei++] = t.C;

            if (insidei == 0)
                return;
            else if (outsidei == 0) // all outside
                ClipRight(t, r);
            else if (insidei == 1)  // inside has the single point.  We're going to pass a single triangle along
            {
                var p1 = inside[0];
                var p2 = ClipLineVertical(p1, outside[0], left);
                var p3 = ClipLineVertical(p1, outside[1], left);
                ClipRight(new Triangle(p1, p2, p3), r);
            }
            else
            {
                var p1 = ClipLineVertical(outside[0], inside[0], left);
                var p2 = ClipLineVertical(outside[0], inside[1], left);
                ClipRight(new Triangle(p1, inside[0], inside[1]), r);

                var thirdPoint = OppositeSidesOfLine(p2, inside[1], p1, inside[0]) ? 0 : 1;
                ClipRight(new Triangle(p1, p2, inside[thirdPoint]), r);
            }
        }

        public void ClipRight(Triangle t, RectangleF s)
        {
            //Record("ClipRight", t);
            var inside = new PointF[3];
            var outside = new PointF[3];
            var insidei = 0;
            var outsidei = 0;

            // Split the triangle's points to the left and right sides of the square's left side
            var right = s.Right;
            if (t.A.X > right)
                outside[outsidei++] = t.A;
            else
                inside[insidei++] = t.A;
            if (t.B.X > right)
                outside[outsidei++] = t.B;
            else
                inside[insidei++] = t.B;
            if (t.C.X > right)
                outside[outsidei++] = t.C;
            else
                inside[insidei++] = t.C;

            if (insidei == 0) // all outside
                return;
            else if (outsidei == 0)
                ClipTop(t, s);
            else if (insidei == 1)  // inside has the single point.  We're going to pass a single triangle along
            {
                var p1 = inside[0];
                var p2 = ClipLineVertical(p1, outside[0], right);
                var p3 = ClipLineVertical(p1, outside[1], right);
                ClipTop(new Triangle(p1, p2, p3), s);
            }
            else
            {
                var p1 = ClipLineVertical(inside[0], outside[0], right);
                var p2 = ClipLineVertical(inside[1], outside[0], right);
                ClipTop(new Triangle(p1, inside[0], inside[1]), s);

                var thirdPoint = OppositeSidesOfLine(p2, inside[0], p1, inside[1]) ? 1 : 0;
                ClipTop(new Triangle(p1, p2, inside[thirdPoint]), s);
            }
        }

        public void ClipTop(Triangle t, RectangleF s)
        {
            //Record("ClipTop", t);
            var outside = new PointF[3];
            var inside = new PointF[3];
            var outsidei = 0;
            var insidei = 0;

            // Split the triangle's points to the left and right sides of the square's left side
            var top = s.Top;
            if (t.A.Y > top) // y increases down
                inside[insidei++] = t.A;
            else
                outside[outsidei++] = t.A;
            if (t.B.Y > top)
                inside[insidei++] = t.B;
            else
                outside[outsidei++] = t.B;
            if (t.C.Y > top)
                inside[insidei++] = t.C;
            else
                outside[outsidei++] = t.C;

            if (insidei == 0) // all outside
                return;
            else if (outsidei == 0)
                ClipBottom(t, s);
            else if (insidei == 1)  // inside has the single point.  We're going to pass a single triangle along
            {
                var p1 = inside[0];
                var p2 = ClipLineHorizontal(p1, outside[0], top);
                var p3 = ClipLineHorizontal(p1, outside[1], top);
                ClipBottom(new Triangle(p1, p2, p3), s);
            }
            else
            {
                var p1 = ClipLineHorizontal(outside[0], inside[0], top);
                var p2 = ClipLineHorizontal(outside[0], inside[1], top);
                ClipBottom(new Triangle(p1, inside[0], inside[1]), s);

                var thirdPoint = OppositeSidesOfLine(p2, inside[0], p1, inside[1]) ? 1 : 0;
                ClipBottom(new Triangle(p1, p2, inside[thirdPoint]), s);
            }
        }

        public void ClipBottom(Triangle t, RectangleF s)
        {
            //Record("ClipBottom", t);
            var inside = new PointF[3];
            var outside = new PointF[3];
            var insidei = 0;
            var outsidei = 0;

            // Split the triangle's points to the left and right sides of the square's left side
            var bottom = s.Bottom;
            if (t.A.Y > bottom)
                outside[outsidei++] = t.A;
            else
                inside[insidei++] = t.A;
            if (t.B.Y > bottom)
                outside[outsidei++] = t.B;
            else
                inside[insidei++] = t.B;
            if (t.C.Y > bottom)
                outside[outsidei++] = t.C;
            else
                inside[insidei++] = t.C;

            if (insidei == 0) // all outside
                return;
            else if (outsidei == 0)
                Handler(t, s);
            else if (insidei == 1)  // inside has the single point.  We're going to pass a single triangle along
            {
                var p1 = inside[0];
                var p2 = ClipLineHorizontal(p1, outside[0], bottom);
                var p3 = ClipLineHorizontal(p1, outside[1], bottom);
                Handler(new Triangle(p1, p2, p3), s);
            }
            else
            {
                var p1 = ClipLineHorizontal(inside[0], outside[0], bottom);
                var p2 = ClipLineHorizontal(inside[1], outside[0], bottom);
                Handler(new Triangle(p1, inside[0], inside[1]), s);

                var thirdPoint = OppositeSidesOfLine(p2, inside[0], p1, inside[1]) ? 1 : 0;
                Handler(new Triangle(p1, p2, inside[thirdPoint]), s);
            }
        }

        public virtual void Handler(Triangle t, RectangleF r)
        {
            Records.Add(new LabeledTriangle { Label = "Handle", Triangle = t });
        }

        public PointF ClipLineVertical(PointF a, PointF b, float sideX)
        {
            var slope = (b.Y - a.Y) / (b.X - a.X);
            var sideY = a.Y + slope * (sideX - a.X);
            return new PointF(sideX, sideY);
        }

        public PointF ClipLineHorizontal(PointF a, PointF b, float sideY)
        {
            var slope = (b.X - a.X) / (b.Y - a.Y);
            var sideX = a.X + slope * (sideY - a.Y);
            return new PointF(sideX, sideY);
        }

        /// <summary>
        /// Determine whether two points are on opposite sides of a line
        /// </summary>
        /// <param name="a">One of the two points</param>
        /// <param name="b">The other point</param>
        /// <param name="p1">A point along the line</param>
        /// <param name="p2">Another point along the line</param>
        /// <returns></returns>
        public bool OppositeSidesOfLine(PointF a, PointF b, PointF p1, PointF p2)
        {
            //2 points (ax,ay)(ax,ay) and (bx,by)(bx,by) are on opposite sides of the line (x1,y1)→(x2,y2)(x1,y1)→(x2,y2)?
            //((y1−y2)(ax−x1)+(x2−x1)(ay−y1))((y1−y2)(bx−x1)+(x2−x1)(by−y1))<0.
            return ((p1.Y - p2.Y) * (a.X - p1.X) + (p2.X - p1.X) * (a.Y - p1.Y)) * ((p1.Y - p2.Y) * (b.X - p1.X) + (p2.X - p1.X) * (b.Y - p1.Y)) < 0f;
        }

        private void Record(string label, Triangle t)
        {
            Records.Add(new LabeledTriangle { Label = label, Triangle = t });
        }
    }

    public class LightingTriangleClipper : TriangleClipper
    {
        public PointF toP1;
        public PointF toP2;
        public PointF toP3;
        public PointF toP4;
        public float GridResolution;
        public int Width;
        public int Height;
        public float HalfSunDividedByGridCellArea;
        public float InvertedGridCellArea;
        public int GridCellX;
        public int GridCellY;
        public float[,] ShadowArray;

        public void AddLight(Triangle t, RectangleF r, int ix, int iy, float halfSun)
        {
            GridCellX = ix;
            GridCellY = iy;
            HalfSunDividedByGridCellArea = halfSun;
            //if (GridCellX == 100 && GridCellY == 100)
            //{
            //    Console.WriteLine(@"r={0}", r);
            //}
            ClipLeft(t, r);
        }

        public override void Handler(Triangle t, RectangleF r)
        {
            var increment = HalfSunDividedByGridCellArea * t.Area;
            ShadowArray[GridCellX, GridCellY] += increment;
        }
    }
}
