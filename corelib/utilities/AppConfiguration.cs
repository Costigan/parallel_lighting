using Microsoft.Extensions.Configuration;
using System.Configuration;
using System.Linq;

namespace viper.corelib.utilities
{
    public class AppConfiguration
    {
        protected static string Hostname = null;
        private static IConfigurationRoot _environment;

        public static string Get(string key, string default_value =null)
        {
            GetEnvVars();
            var r = _environment[key];
            if (null != r) return r;
            return default_value;

            //return GetForHost(key, default_value);
        }

        /// <summary>
        /// Return a machine-specific value, if one exists.  Otherwise, look for a global version.
        /// </summary>
        /// <param name="key"></param>
        /// <param name="default_value"></param>
        /// <returns></returns>
        public static string GetForHost(string key, string default_value = null)
        {
            var n = GetHostName();
            var key2 = n + "_" + key;
            var r = System.Configuration.ConfigurationManager.AppSettings[key2];
            if (r != null) return r;
            return GetGlobal(key, default_value);
        }

        public static string GetGlobal(string key, string default_value = null)
        {
            var r = ConfigurationManager.AppSettings[key];
            if (r != null) return r;
            if (default_value != null) return default_value;
            return null;
        }

        static string GetHostName()
        {
            if (Hostname != null) return Hostname;
            var n = System.Environment.MachineName.ToLower();
            if (n.Contains('.'))
                n = n.Substring(0, n.IndexOf('.'));
            Hostname = n;
            return Hostname;
        }

        public static void GetEnvVars(string prefix = "Horizon_")
        {
            if (_environment == null)
                _environment = new Microsoft.Extensions.Configuration.ConfigurationBuilder()
                    .AddEnvironmentVariables(prefix)
                    .Build();
        }
    }
}
