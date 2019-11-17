using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using viper.corelib.math;
using viper.corelib.spice;

namespace viper.corelib.lighting
{
    public class SunVectorCache
    {
        protected static SunVectorCache _singleton = new SunVectorCache();

        // Copied from LunarHorizon.cs
        // The 3D accounts for leap seconds since 2000.  This is valid only for dates after Jul 1 2012.
        public const int MaxCacheCount = 100000;
        public const int EarthId = 399;
        public const int MoonId = 301;
        public const int SunId = 10;
        protected static long _epoch = new DateTime(2000, 1, 1, 11, 58, 55, 816).Ticks;

        protected long _minimumStep = new TimeSpan(1, 0, 0).Ticks;
        protected Dictionary<long, Vector3d> _sunCache = new Dictionary<long, Vector3d>();
        protected Dictionary<long, Vector3d> _earthCache = new Dictionary<long, Vector3d>();

        public static SunVectorCache GetSingleton() => _singleton != null ? _singleton : _singleton = new SunVectorCache();

        public double[] _state = new double[6];

        public TimeSpan MinimumStep
        {
            get { return new TimeSpan(_minimumStep); }
            set { _minimumStep = value.Ticks; }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector3d SunPosition(DateTime time)
        {
            lock (this)
            {
                var ticks = time.Ticks;
                var rep = ticks - (ticks / _minimumStep);
                if (_sunCache.TryGetValue(rep, out Vector3d v))
                    return v;
                var et = (rep - _epoch) / 10000000L + 3d;
                double lt = 0d;
                CSpice.spkgeo_c(SunId, et, "MOON_ME", MoonId, _state, ref lt);
                v = new Vector3d(_state[0], _state[1], _state[2]);
                if (_sunCache.Count < MaxCacheCount)
                    _sunCache.Add(rep, v);
                return v;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Vector3d EarthPosition(DateTime time)
        {
            lock (this)
            {
                var ticks = time.Ticks;
                var rep = ticks - (ticks / _minimumStep);
                if (_earthCache.TryGetValue(rep, out Vector3d v))
                    return v;
                var et = (rep - _epoch) / 10000000L + 3d;
                double lt = 0d;
                CSpice.spkgeo_c(EarthId, et, "MOON_ME", MoonId, _state, ref lt);
                v = new Vector3d(_state[0], _state[1], _state[2]);
                if (_earthCache.Count < MaxCacheCount)
                    _earthCache.Add(rep, v);
                return v;
            }
        }
    }
}
