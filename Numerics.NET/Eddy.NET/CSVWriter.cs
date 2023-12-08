using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using CsvHelper;
using CsvHelper.Configuration;
using System.Globalization;

namespace Eddy.NET
{
    internal static class CSVWriter
    {
        private static string ConvertToCsv(List<PhysicsOut> list)
        {
            using (var writer = new StringWriter())
            using (var csv = new CsvWriter(writer, new CsvConfiguration(CultureInfo.InvariantCulture)))
            {
                csv.WriteRecords(list);
                return writer.ToString();
            }
        }

        private static string ArrayToCsv(double[,] data){
            string output = "";
        }

        private static void SaveCsvToFile(string csv, string filePath)
        {
            FileInfo file = new FileInfo(filePath);
            file.Directory.Create();
            File.WriteAllText(file.FullName, csv);
            Console.WriteLine("CSV data saved to file: " + filePath);
        }

        public static void SaveToFile(List<PhysicsOut> list, string filePath)
        {
            SaveCsvToFile(ConvertToCsv(list), filePath);
        }
    }    
}
