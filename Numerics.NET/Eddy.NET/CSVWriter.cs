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
        }
    }    
}
