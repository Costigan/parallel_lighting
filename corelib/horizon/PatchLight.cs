using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using viper.corelib.patch;

namespace viper.corelib.horizon
{
    public class PatchLight
    {
        public const int PatchSize = TerrainPatch.DefaultSize;
        public const int PatchDrawSize = 4;
        public const int DEMSize = 30400;
        public const int PatchesPerSize = DEMSize / PatchSize;

        public static bool DistanceFilter = true;
        public static int DistanceFilterThreshold = 90 * PatchSize;  // Pixel distance

        public static float SpacingSafetyFactor = 1f;
        public static bool Verbose = false;

        public Point Id;
        public int Step = 1;

        public static TimeSpan EstimateCalculationTime(IEnumerable<TerrainPatch> queue) => EstimateCalculationTime(queue.Select(p => p.Id));

        public static TimeSpan EstimateCalculationTime(IEnumerable<Point> ids)
        {
            var count = ids.Sum(id => (new PatchLight { Id = id }).GenerateFarField().Count);
            return TimeSpan.FromSeconds(count * .1625f);
        }

        public List<PatchLight> GenerateFarField()
        {
            var patches = new List<PatchLight>();
            var other = new PatchLight { Id = new Point(0, 0), Step = 256 };
            GenerateFarField(other, patches);
            return patches;
        }

        public void GenerateFarField(PatchLight other, List<PatchLight> patches)
        {
            if (other.Step < 1)
                return;
            if (DistanceFilter && MinimumPixelDistance(other) > DistanceFilterThreshold)
            {
                // Ignore this patch based on distance from the center
            }
            else if (Overlaps(other) || !AdequateResolution(other, 360 * 4))
            {
                if (Verbose)
                    Console.WriteLine($"Recur: {other}");
                var s = other.Step / 2;  // recur
                GenerateFarField( new PatchLight { Id = new Point(other.Id.X + 0, other.Id.Y + 0), Step = s }, patches);
                GenerateFarField( new PatchLight { Id = new Point(other.Id.X + s, other.Id.Y + 0), Step = s }, patches);
                GenerateFarField( new PatchLight { Id = new Point(other.Id.X + 0, other.Id.Y + s), Step = s }, patches);
                GenerateFarField( new PatchLight { Id = new Point(other.Id.X + s, other.Id.Y + s), Step = s }, patches);
            }
            else
            {
                if (Verbose)
                    Console.WriteLine($"Add:  {other}");
                patches.Add(other);
            }
        }

        public Rectangle Rectangle => new Rectangle(Id.X * PatchSize, Id.Y * PatchSize, (PatchSize - 1) * Step, (PatchSize - 1) * Step);

        // Does the other patch overlap this one?
        public bool Overlaps(PatchLight o)
        {
            if (o.Id.X + o.Step <= Id.X) return false;
            if (o.Id.Y + o.Step <= Id.Y) return false;
            if (o.Id.X >= Id.X + Step) return false;
            if (o.Id.Y >= Id.Y + Step) return false;
            return true;
        }

        // Is the spacing of the current patch tight enough when viewed from the center?
        public bool AdequateResolution(PatchLight other, int horizon_count)
        {
            var distance = MinimumPixelDistance(other);
            //var spacing_rad = 1f / distance;
            var spacing_rad = (float)other.Step / distance;
            var max_allowed_rad = 2f * 3.141592653f / horizon_count;
            var isAdequate = other.Step <= 1 || spacing_rad * SpacingSafetyFactor < max_allowed_rad;
            if (Verbose)
                Console.WriteLine($"distance: {this} -> {other} d={distance} spacing={spacing_rad} max={max_allowed_rad} adequate={isAdequate}");
            return isAdequate;
        }

        public int MinimumPixelDistance(PatchLight other)
        {
            var r1 = Rectangle;
            var r2 = other.Rectangle;

            if (r1.Top >= r2.Bottom)
            {
                if (r2.Right <= r1.Left)
                    return Distance(UpperLeft(r1), LowerRight(r2));
                if (r1.Right <= r2.Left)
                    return Distance(UpperRight(r1), LowerLeft(r2));
                return r1.Top - r2.Bottom;
            }
            else if (r1.Bottom <= r2.Top)
            {
                if (r2.Right <= r1.Left)
                    return Distance(LowerLeft(r1), UpperRight(r2));
                if (r1.Right <= r2.Left)
                    return Distance(LowerRight(r1), UpperLeft(r2));
                return r2.Top - r1.Bottom;
            }
            else
            {
                if (r2.Right <= r1.Left)
                    return r1.Left - r2.Right;
                if (r1.Left <= r2.Right)
                    return r2.Left - r1.Right;
                throw new Exception("Shouldn't get here due to overlap");
            }
        }

        Point UpperLeft(Rectangle r) => new Point(r.Left, r.Top);
        Point UpperRight(Rectangle r) => new Point(r.Right, r.Top);
        Point LowerLeft(Rectangle r) => new Point(r.Left, r.Bottom);
        Point LowerRight(Rectangle r) => new Point(r.Right, r.Bottom);
        int Distance(Point a, Point b) => (int)Math.Sqrt((a.X - b.X) * (a.X - b.X) + (a.Y - b.Y) * (a.Y - b.Y));

        public override string ToString() => $"<TerrainPatch x={Id.X} y={Id.Y} stride={Step}>";

        public int IdHash => Id.X | Id.Y << 10 | Step << 20;

        public TerrainPatch ToTerrainPatch()
        {
            var tp = TerrainPatch.FromId(Id);
            tp.Step = Step;
            return tp;
        }
    }
}
