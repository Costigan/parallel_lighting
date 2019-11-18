using Microsoft.Extensions.Configuration;
using System.Configuration;
using System.Linq;

namespace viper.corelib.utilities
{
    public class AppConfiguration
    {
        public const string ConfigPath = @"config.json";
        public static string HostName = null;
        private static IConfigurationRoot _environment;

        /// <summary>
        /// Return a machine-specific value, if one exists.  Otherwise, look for a global version.
        /// </summary>
        /// <param name="key"></param>
        /// <param name="default_value"></param>
        /// <returns></returns>
        public static string Get(string key, string default_value = null)
        {
            InitializeVariables();
            var key2 = HostName + "_" + key;
            var r = _environment[key2];
            if (r != null) return r;
            return GetGlobal(key, default_value);
        }

        public static string GetGlobal(string key, string default_value = null)
        {
            InitializeVariables();
            var r = _environment[key];
            if (r != null) return r;
            if (default_value != null) return default_value;
            return null;
        }

        protected static string GetHostName()
        {
            if (HostName != null) return HostName;
            var n = System.Environment.MachineName.ToLower();
            if (n.Contains('.'))
                n = n.Substring(0, n.IndexOf('.'));
            HostName = n;
            return HostName;
        }

        protected static void InitializeVariables(string prefix = "Horizon_")
        {
            if (_environment != null)
                return;
            HostName = GetHostName();
            _environment = new ConfigurationBuilder()
                    .AddJsonFile(ConfigPath)
                    .AddEnvironmentVariables(prefix)
                    .Build();
        }
    }
}
