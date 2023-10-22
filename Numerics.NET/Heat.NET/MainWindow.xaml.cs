using System;
using System.Collections.Generic;
using System.IO;
using Newtonsoft.Json;
using System.Windows;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Threading;

namespace Heat.NET
{
    public partial class MainWindow : Window
    {
        private List<double[,]> _snapshots;
        private int _currentSnapshotIndex = 0;
        private DispatcherTimer _timer = new DispatcherTimer();

        public MainWindow()
        {
            InitializeComponent();
            EnsureDataGenerated();
            _snapshots = JsonConvert.DeserializeObject<List<double[,]>>(File.ReadAllText("simulationData.json"));
            _timer.Interval = TimeSpan.FromMilliseconds(100);
            _timer.Tick += Timer_Tick;
        }

        private void EnsureDataGenerated()
        {
            if (!File.Exists("simulationData.json"))
            {
                var simulation = new Simulation();
                simulation.Run();
            }
        }

        private void PlayButton_Click(object sender, RoutedEventArgs e)
        {
            _timer.Start();
        }

        private void ReplayButton_Click(object sender, RoutedEventArgs e)
        {
            _currentSnapshotIndex = 0;
            _timer.Start();
        }

        private void Timer_Tick(object sender, EventArgs e)
        {
            if (_currentSnapshotIndex >= _snapshots.Count)
            {
                _timer.Stop();
                return;
            }
            double[,] currentData = _snapshots[_currentSnapshotIndex++];
            DisplayHeatMap(currentData);
        }

        private void DisplayHeatMap(double[,] data)
        {
            int width = data.GetLength(0);
            int height = data.GetLength(1);

            WriteableBitmap heatmap = new WriteableBitmap(width, height, 96, 96, PixelFormats.Bgra32, null);

            int[] pixelData = new int[width * height];

            double maxValue = 10; // You can adjust this according to your data range
            for (int i = 0; i < width; i++)
            {
                for (int j = 0; j < height; j++)
                {
                    double normalized = data[i, j] / maxValue;
                    Color heatmapColor = GetHeatMapColor(normalized);
                    int pixel = heatmapColor.B << 24 | heatmapColor.G << 16 | heatmapColor.R << 8 | 0xFF;
                    pixelData[i + j * width] = pixel;
                }
            }
            heatmap.WritePixels(new Int32Rect(0, 0, width, height), pixelData, width * sizeof(int), 0);
            HeatmapImage.Source = heatmap;
        }

        private Color GetHeatMapColor(double value)
        {
            // Simple blue to red gradient
            byte bluePart = 0;
            byte redPart = 0;

            if (value < 0.5)
            {
                bluePart = (byte)(255 - value * 2 * 255);
                redPart = (byte)(value * 2 * 255);
            }
            else
            {
                bluePart = 0;
                redPart = 255;
            }

            return Color.FromArgb(255, redPart, 0, bluePart);
        }
    }
}
