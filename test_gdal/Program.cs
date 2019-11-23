using System;

namespace test_gdal
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Hello World!");
            OSGeo.GDAL.Gdal.AllRegister();
            Console.WriteLine("After GDAL AllRegister");
        }
    }
}
