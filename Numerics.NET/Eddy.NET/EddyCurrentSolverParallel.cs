using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Eddy.NET
{
    internal class EddyCurrentSolverParallel : IDifferentialEquationSolver
    {
        private double Dx;
        private double Dt;
        private double Alpha;
        private double Beta;

        private Vector2D[,] A;
        private Vector2D[,] APrevious;
        private Vector2D[,] A0;
        private Vector2D[,] A0Previous;

        private Vector2D F;
        private Vector2D s, v, a;

        private readonly SimulationSettings _settings;

        public EddyCurrentSolverParallel(SimulationSettings settings)
        {
            _settings = settings;
            Initialize();
        }

        public void Initialize()
        {
            A = new Vector2D[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            APrevious = new Vector2D[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            A0 = new Vector2D[_settings.ResolutionSpace, _settings.ResolutionSpace];
            A0Previous = new Vector2D[_settings.ResolutionSpace, _settings.ResolutionSpace];

            Dx = _settings.Length / _settings.ResolutionSpace;
            Dt = (_settings.TimeEnd - _settings.TimeStart) / _settings.ResolutionTime;
            Alpha = 1.0 / (_settings.MaterialConductivity * _settings.MaterialMagneticPermeability);
            Beta = Alpha * Dt / (2 * Dx * Dx);
        }

        public void Solve()
        {
            double currentTime = _settings.TimeStart;
            while (currentTime < _settings.TimeEnd)
            {
                SolveForOneTimeStep();
                currentTime += Dt;
                Console.WriteLine($"Time: {currentTime}");
            }
        }

        public void PrintResults()
        {
            // ... Logic to print the results of the simulation
        }


        private void EmbedA0()
        {
            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {
                    A0[i, j] = new Vector2D(1, 1);
                }
            });
        }

        private void SolveForOneTimeStep()
        {
            Array.Copy(APrevious, A, APrevious.Length);

            bool isConverged = false;
            int iterationCount = 0;

            while (!isConverged && iterationCount < _settings.MaxIterations)
            {
                isConverged = Update();

                iterationCount++;
            }

            Console.Write($"Converged in {iterationCount} iterations. ");

            Array.Copy(A, APrevious, APrevious.Length);
        }

        private bool Update()
        {
            bool isConverged = true;
            double[,] delta = new double[_settings.ResolutionSpace, _settings.ResolutionSpace];

            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {
                    Vector2D oldValue = A[i + 1, j + 1];
                    A[i + 1, j + 1] = (1 - _settings.Omega) * A[i + 1, j + 1]
                                    + _settings.Omega / (1 + 4 * Beta) * (APrevious[i + 1, j + 1] + Beta *
                                      (A[i + 2, j + 1] + A[i, j + 1] + A[i + 1, j + 2] + A[i + 1, j]
                                      + APrevious[i + 2, j + 1] + APrevious[i, j + 1] + APrevious[i + 1, j + 2] + APrevious[i + 1, j] - 4 * APrevious[i + 1, j + 1]));

                    delta[i, j] = (A[i + 1, j + 1] - oldValue).Magnitude();

                    if (delta[i, j] > _settings.Tolerance)
                    {
                        isConverged = false;
                    }
                }
            });

            return isConverged;
        }
    }

}
