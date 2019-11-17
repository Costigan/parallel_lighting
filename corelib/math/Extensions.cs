using System;
using System.Collections.Generic;
using System.Drawing;

namespace viper.corelib.math
{
    public static class Extensions
    {
        public static float Distance(this PointF a, PointF b)
        {
            var dx = a.X - b.X;
            var dy = a.Y - b.Y;
            return (float)Math.Sqrt(dx * dx + dy * dy);
        }

        public static bool Equivalent(this PointF a, PointF b, float epsilon = 0.01f) => a.Distance(b) < epsilon;
        public static PointF Plus(this PointF a, PointF b) => new PointF(a.X + b.X, a.Y + b.Y);
        public static PointF Plus(this PointF a, SizeF b) => new PointF(a.X + b.Width, a.Y + b.Height);
        public static PointF Minus(this PointF a, PointF b) => new PointF(a.X - b.X, a.Y - b.Y);
        public static PointF Times(this PointF a, float b) => new PointF(a.X * b, a.Y * b);
        public static PointF Corner(RectangleF r) => new PointF(r.Right, r.Bottom);

        public static float Clamp(this float a, float low, float high) => a < low ? low : (a > high ? high : a);

        public static void ForEach<T>(this IEnumerable<T> enumeration, Action<T> action)
        {
            foreach (T item in enumeration)
                action(item);
        }

        public static double[] ToArray(this Vector3d v) => new double[] { v.X, v.Y, v.Z };

        public static void Populate<T>(this T[] arr, T value)
        {
            for (int i = 0; i < arr.Length; i++)
                arr[i] = value;
        }
    }
}
