using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Linq;

namespace viper.corelib.terrain
{
    public class Horizon
    {
        public const int HorizonSamples = 360 * 4;
        public const float HorizonSamplesF = 360 * 4f;
        public const double HorizonSamplesD = HorizonSamples - 1;
        public const float HorizonSamples_float = HorizonSamples - 1;
        private const double max_elevation_angle_deg = 20d * Math.PI / 180d;  // 20 degrees
        private const short max_elevation = short.MaxValue - 1;  // remember this -1
        private const double max_elevation_deg = (double)max_elevation;

        private float[] _buffer;

        public bool IsDegrees = false;

        public float[] Buffer => _buffer;

        public Horizon()
        {
            Init();
        }

        public Horizon Init()
        {
            if (_buffer==null)
                _buffer = new float[HorizonSamples];
            for (var i = 0; i < HorizonSamples; i++)
                _buffer[i] = float.MinValue;
            return this;
        }

        public static Horizon From(float[] a)
        {
            var h = new Horizon();
            h.Load(a);
            return h;
        }

        public static Horizon From(float[] a, int offset)
        {
            var h = new Horizon();
            h.Load(a, offset);
            return h;
        }

        public void Load(float[] a) => Array.Copy(a, _buffer, Math.Min(a.Length, _buffer.Length));
        public void Load(float[] a, int offset) => Array.Copy(a, offset, _buffer, 0, _buffer.Length);

        public void Merge(double azimuth, double slope)
        {
            // Simplify this
            var slopef = (float)slope;

            var normalized_azimuth = HorizonSamplesD * azimuth / (2d * Math.PI);  // range [0,1)
            var horizon_index = (int)(0.5d + normalized_azimuth);
            Debug.Assert(horizon_index >= 0 && horizon_index < HorizonSamples);

            if (slopef > _buffer[horizon_index])
            {
                _buffer[horizon_index] = slopef;
            }
        }

        public void Merge(Horizon other)
        {
            var other_horizon = other.Buffer;
            lock (_buffer)
            {
                for (var i = 0; i < HorizonSamples; i++)
                {
                    if (other_horizon[i] > _buffer[i])
                        _buffer[i] = other_horizon[i];
                }
            }
        }

        public void ConvertSlopeToDegrees()
        {
            if (IsDegrees)
                return;
            for (var i=0;i<_buffer.Length;i++)
                _buffer[i] = ((float)Math.Atan(_buffer[i])) * 180f / 3.141592653589f;
            IsDegrees = true;
        }

        public float[] GetSlopeDegreeBuffer() => IsDegrees ? _buffer : _buffer.Select(v => (float)(Math.Atan(v) * 180d / Math.PI)).ToArray();

        // First cut at this
        public float SunFraction(float azimuth_deg, float elevation_deg)
        {
            Debug.Assert(IsDegrees);
            var normalized_azimuth = azimuth_deg / 360f;  // range [0,1)
            var horizon_index = (int)(0.5d + HorizonSamplesD * normalized_azimuth);
            Debug.Assert(horizon_index >= 0 && horizon_index < HorizonSamples);
            var horizon_deg = _buffer[horizon_index];

            // It should be possible to optimize this more
            const float SunRadius = 0.25f;
            var delta = (elevation_deg - horizon_deg) / SunRadius;
            if (delta >= 1f)
                return 1f;
            if (delta <= -1f)
                return 0f;
            var frac = Math.Abs(delta);  // [0,1]
            var half_angle = Math.Acos(frac);
            var two_triangles = frac * Math.Sin(half_angle);
            //var wedge = (Math.PI) * (half_angle / Math.PI);  // radius=1, so simplified out
            var wedge = half_angle;
            var other = wedge - two_triangles;
            var ret = delta < 0f ? other : Math.PI * 2d - other;
            ret /= 2f * 3.141592653589f;
            return (float)ret;
        }

        const float SunHalfAngle_deg = 0.54f / 2f;
        const int half_circle_ticks = 8;
        static float[] half_circle = PremultiplyHalfCircle(MakeHalfCircle());
        static float max_photons = 2f * PremultiplyHalfCircle(MakeHalfCircle()).Sum();

        static float[] MakeHalfCircle() => Enumerable.Range(0, half_circle_ticks * 2).Select(i => (float)(Math.Sqrt(64d - (half_circle_ticks - 0.5d - i) * (half_circle_ticks - 0.5d - i)) / (half_circle_ticks))).ToArray();
        static float[] PremultiplyHalfCircle(float[] h) => h.Select(v => v * SunHalfAngle_deg).ToArray();

        public float SunFraction2(float az_deg, float el_deg)
        {
            const float BucketWidth_deg = 360f / HorizonSamplesF;
            const float BucketHalfWidth_deg = BucketWidth_deg / 2f;  // Used because we're interpolating between two buckets.

            const float frac_step_per_half_circle_index = SunHalfAngle_deg / BucketWidth_deg / half_circle_ticks;

            var sun_left_deg = az_deg - SunHalfAngle_deg - BucketHalfWidth_deg;   // 
            var sun_left_bucket_float = sun_left_deg * (HorizonSamplesF / 360f);
            var sun_left_bucket = (int)sun_left_bucket_float; // [0,1440)
            var frac = sun_left_bucket_float - sun_left_bucket;
            if (sun_left_bucket < 0) sun_left_bucket += HorizonSamples;  // [0,1440)

            // We're going to interpolate the horizon between two buckets.  sun_left_bucket is the index
            // of the left most of the two buckets for the first interpolation.  frac is the fraction of the
            // way that the sun's left edge between these two buckets

            var left_bucket_index = sun_left_bucket;
            var right_bucket_index = left_bucket_index + 1;
            if (right_bucket_index >= HorizonSamples) right_bucket_index -= HorizonSamples;

            var left_bucket_elevation_deg = _buffer[left_bucket_index];
            var right_bucket_elevation_deg = _buffer[right_bucket_index];
            var bucket_delta_deg = right_bucket_elevation_deg - left_bucket_elevation_deg;

            var photons = 0f;

            for (var i = 0; i < half_circle.Length; i++)
            {
                var horizon_elevation_deg = frac * bucket_delta_deg + left_bucket_elevation_deg;
                var sun_column_deg = half_circle[i];    // This is now pre-multiplied * SunHalfAngle_deg;
                var sun_top_deg = el_deg + sun_column_deg;

                if (horizon_elevation_deg >= sun_top_deg)
                    goto continue_loop;

                var angle_delta = sun_top_deg - horizon_elevation_deg;
                var sun_column_deg2 = sun_column_deg + sun_column_deg;  // Can't have more than this much light
                if (angle_delta > sun_column_deg2)
                    angle_delta = sun_column_deg2;

                photons += angle_delta;

                continue_loop:

                frac += frac_step_per_half_circle_index;
                if (frac < 1f)
                    continue;

                // Move the bucket
                left_bucket_index = right_bucket_index;
                right_bucket_index = left_bucket_index + 1;
                if (right_bucket_index >= HorizonSamples) right_bucket_index -= HorizonSamples;

                left_bucket_elevation_deg = right_bucket_elevation_deg;
                right_bucket_elevation_deg = _buffer[right_bucket_index];
                bucket_delta_deg = right_bucket_elevation_deg - left_bucket_elevation_deg;

                frac -= 1f;
            }

            var sun_fraction = photons / max_photons;

            return sun_fraction;
        }

        public bool IsOverHorizon(double azimuth_rad, double elevation_deg_or_slope)
        {
            var normalized_azimuth = azimuth_rad / (2d * Math.PI);  // range [0,1)
            var horizon_index = (int)(0.5d + HorizonSamplesD * normalized_azimuth);
            horizon_index = Modulo(horizon_index, 0, HorizonSamples);
            Debug.Assert(horizon_index >= 0 && horizon_index < HorizonSamples);
            var horizon_elevation = _buffer[horizon_index];
            return ((float)elevation_deg_or_slope) > horizon_elevation;
        }

        //public bool IsOverHorizon(double azimuth_rad, double elevation_deg) => OverHorizon((float)azimuth_rad, (float)elevation_deg) > 0f;

        // Elevation may be either in radians or degrees
        public float OverHorizon(float azimuth_rad, float target_elevation)
        {
            var normalized_azimuth = azimuth_rad / (2f * 3.141592653589f);  // range [0,1)
            var horizon_index_float = HorizonSamples_float * normalized_azimuth;
            var left = (int)horizon_index_float;
            var frac = horizon_index_float - left;
            var right = Modulo(left + 1, 0, HorizonSamples);
            var elevation_left = _buffer[left];
            var elevation_right = _buffer[right];
            var interpolated_horizon = elevation_left + frac * (elevation_right - elevation_left);
            return target_elevation - interpolated_horizon;
        }

        protected int Modulo(int i, int low, int high)
        {
            if (i < low)
                return i + (high - low);
            if (i >= high)
                return i - (high - low);
            return i;
        }

        /// <summary>
        /// Return stream of points or null.  Null means no data for some angles.
        /// </summary>
        /// <returns></returns>
        public IEnumerable<PointF?> PartialPoints()
        {
            var flag = false;  // true when the last azimuth had data
            for (var i = 0; i < HorizonSamples; i++)
            {
                var v = _buffer[i];
                if (v == float.MinValue)
                {
                    if (flag)
                    {
                        yield return null;
                        flag = false;
                    }
                }
                else
                {
                    yield return new PointF(360f * i / HorizonSamples, v);
                    flag = true;
                }
            }
        }
    }
}
