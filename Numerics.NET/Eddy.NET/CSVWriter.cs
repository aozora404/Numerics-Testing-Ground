using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CsvHelper;
using CsvHelper.Configuration;

namespace Eddy.NET
{
    public class PhysicsOut
    {
        public Vector Position { get; set; };
        public Vector Velocity { get; set; };
        public Vector Force { get; set; };

        public PhysicsOut()
        {
            Position = new Vector(0, 0, 0);
            Velocity = new Vector(0, 0, 0);
            Force = new Vector(0, 0, 0);
        }

        public set(Vector position, Vector velocity, Vector force){
            Position = position;
            Velocity = velocity;
            Force = force;
        }
    }
    
    private static string ConvertToCsv<T>(List<T> list)
    {
        using (var writer = new StringWriter())
        using (var csv = new CsvWriter(writer, new CsvConfiguration(CultureInfo.InvariantCulture)))
        {
            csv.WriteRecords(list);
            return writer.ToString();
        }
    }

    private static void SaveCsvToFile(string csv, string filePath)
    {
        File.WriteAllText(filePath, csv);
        Console.WriteLine("CSV data saved to file: " + filePath);
    }

    public static void SaveToFile<T>(List<T> list, string filePath)
    {
        SaveCsvToFile(ConvertToCsv<T>(list), filePath);
        Console.WriteLine("CSV data saved to file: " + filePath);
    }
}
