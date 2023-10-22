using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Eddy.NET
{
    internal class PDESolver
    {
        private const double Length = 1.0;
        private const int ResolutionSpace = 32;
        private const double Dx = Length / ResolutionSpace;

        private const double TimeStart = 0.0;
        private const double TimeEnd = 1.0;
        private const int ResolutionTime = 1000;
        private const double Dt = (TimeEnd - TimeStart) / ResolutionTime;

        private const double Omega = 1.5;
        private const double Tolerance = 10E-9;
        private const int MaxIterations = 5000;

        private const double Alpha = 0.1;
        private const double Beta = Alpha * Dt / (2 * Dx * Dx);

        private double[,] u;
        private double[,] uPrevious;

        public PDESolver()
        {
            InitializeArrays();
        }

        private void InitializeArrays()
        {
            u = new double[ResolutionSpace + 2, ResolutionSpace + 2];
            uPrevious = new double[ResolutionSpace + 2, ResolutionSpace + 2];

            for (int i = 0; i < ResolutionSpace + 2; i++)
            {
                for (int j = 0; j < ResolutionSpace + 2; j++)
                {
                    if (i > ResolutionSpace / 4 && i < 3 * ResolutionSpace / 4 &&
                        j > ResolutionSpace / 4 && j < 3 * ResolutionSpace / 4)
                    {
                        uPrevious[i, j] = 10;
                    }
                }
            }
        }

        public void Solve()
        {
            double currentTime = TimeStart;
            while (currentTime < TimeEnd)
            {
                SolveForOneTimeStep();
                currentTime += Dt;
                Console.WriteLine($"Time: {currentTime}");
            }
        }

        private void SolveForOneTimeStep()
        {
            Array.Copy(uPrevious, u, uPrevious.Length);

            bool isConverged = false;
            int iterationCount = 0;

            while (!isConverged && iterationCount < MaxIterations)
            {
                isConverged = UpdateUsingSOR();

                iterationCount++;
            }

            Console.WriteLine($"Converged in {iterationCount} iterations.");

            Array.Copy(u, uPrevious, uPrevious.Length);
        }

        private bool UpdateUsingSOR()
        {
            bool isConverged = true;
            double[,] delta = new double[ResolutionSpace, ResolutionSpace];

            Parallel.For(0, ResolutionSpace, i =>
            {
                for (int j = 0; j < ResolutionSpace; j++)
                {
                    double oldValue = u[i + 1, j + 1];
                    u[i + 1, j + 1] = (1 - Omega) * u[i + 1, j + 1]
                                    + Omega / (1 + 4 * Beta) * (uPrevious[i + 1, j + 1] + Beta *
                                      (u[i + 2, j + 1] + u[i, j + 1] + u[i + 1, j + 2] + u[i + 1, j]
                                      + uPrevious[i + 2, j + 1] + uPrevious[i, j + 1] + uPrevious[i + 1, j + 2] + uPrevious[i + 1, j] - 4 * uPrevious[i + 1, j + 1]));

                    delta[i, j] = Math.Abs(u[i + 1, j + 1] - oldValue);

                    if (delta[i, j] > Tolerance)
                    {
                        isConverged = false;
                    }
                }
            });

            return isConverged;
        }

        public void PrintResults()
        {
            for (int i = 0; i < ResolutionSpace + 2; i++)
            {
                for (int j = 0; j < ResolutionSpace + 2; j++)
                {
                    Console.Write($"{u[i, j]:f1} ");
                }
                Console.WriteLine();
            }
        }
    }
}
