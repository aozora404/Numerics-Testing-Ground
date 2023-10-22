using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Eddy.NET
{
    internal class DirectFieldSolverMultigrid : IDifferentialEquationSolver
    {
        private double Dx;
        private double Dt;
        private double Alpha;
        private double Beta;
        private double Mass;

        private bool[,] isMaterial;

        private Vector[,] B;
        private Vector[,] BPrevious;
        private Vector[,] B0;
        private Vector[,] B0Previous;
        private Vector[,] A;
        private Vector[,] APrevious;

        private Vector[,] J;
        private Vector[,] JPrevious;
        private Vector[][,] JGrids = new Vector[maxLevels][,];
        private const int numPreRelaxations = 3;
        private const int numPostRelaxations = 3;
        private const int maxLevels = 3;


        private Vector force;
        private Vector position, velocity, acceleration;

        private readonly SimulationSettings _settings;

        public DirectFieldSolverMultigrid(SimulationSettings settings)
        {
            _settings = settings;

            Dx = _settings.Length / _settings.ResolutionSpace;
            Dt = (_settings.TimeEnd - _settings.TimeStart) / _settings.ResolutionTime;
            Alpha = 1.0 / (_settings.MaterialConductivity * _settings.MaterialMagneticPermeability);
            Beta = Alpha * Dt / (2 * Dx * Dx);
            Mass = _settings.MaterialDensity * (4 / 3) * Math.PI * _settings.BallRadius * _settings.BallRadius * _settings.BallRadius;

            B = new Vector[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            BPrevious = new Vector[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            B0 = new Vector[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            B0Previous = new Vector[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            J = new Vector[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            JPrevious = new Vector[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            A = new Vector[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            APrevious = new Vector[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];


            isMaterial = new bool[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];

            Parallel.For(0, _settings.ResolutionSpace + 2, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace + 2; j++)
                {
                    B[i, j] = new Vector(0, 0, 0);
                    BPrevious[i, j] = new Vector(0, 0, 0);
                    B0[i, j] = new Vector(0, 0, 0);
                    B0Previous[i, j] = new Vector(0, 0, 0);
                    J[i, j] = new Vector(0, 0, 0);
                    JPrevious[i, j] = new Vector(0, 0, 0);
                    A[i, j] = new Vector(0, 0, 0);
                    APrevious[i, j] = new Vector(0, 0, 0);

                    double x = (i * 1.0 / _settings.ResolutionSpace) - 0.5;
                    double y = (j * 1.0 / _settings.ResolutionSpace) - 0.5;
                    double r = _settings.BallRadius / _settings.Length;
                    r = r * r;

                    if (x * x + y * y <= r)
                    {
                        isMaterial[i, j] = true;
                    }
                    else
                    {
                        isMaterial[i, j] = false;
                    }
                }
            });

            for (int level = 0; level < maxLevels; level++)
            {
                int size = _settings.ResolutionSpace / (int)Math.Pow(2, level) + 2;  // Halve the grid size for each coarser level
                JGrids[level] = new Vector[size, size];
            }

            force = new Vector(0, 0, 0);
            position = new Vector(0, 0, 0);
            velocity = new Vector(_settings.InitialVelocity, 0, 0);
            acceleration = new Vector(0, 0, 0);
        }

        public void Solve()
        {
            double currentTime = _settings.TimeStart;
            while (currentTime < _settings.TimeEnd)
            {
                EMStep();
                DynamicsStep();
                TimeStep();
                currentTime += Dt;
                Console.WriteLine($"Position: {position:f2} Velocity: {velocity:f2} Time: {currentTime:f2} ");
            }
        }

        public void PrintResults()
        {
            // ... Logic to print the results of the simulation
        }

        private void EMStep()
        {
            EmbedB0();

            bool isConverged = false;
            int iterationCount = 0;

            while (!isConverged && iterationCount < _settings.MaxIterations)
            {
                isConverged = UpdateB();

                iterationCount++;
            }

            Console.Write($"B Converged in {iterationCount} iterations. ");

            CalculateA();
            CalculateJ();

            Console.Write($"J Converged in {iterationCount} iterations. ");

            CalculateForce();
        }

        private void EmbedB0()
        {
            Vector relPosition = new Vector(0, 0, 0);
            Vector s = new Vector(0, 0, 0);
            Vector coilOrigin = new Vector(_settings.CoilDistance, 0, 0);

            Parallel.For(0, _settings.ResolutionSpace + 2, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace + 2; j++)
                {
                    relPosition.X = (i - _settings.ResolutionSpace / 2) * Dx;
                    relPosition.Y = (j - _settings.ResolutionSpace / 2) * Dx;

                    s = relPosition + position - coilOrigin;


                    if (s.Magnitude() > _settings.CoilRadius)
                    {
                        B0[i, j] = new Vector(0, 0, _settings.CoilMagneticFieldStrength);
                    }
                    else
                    {
                        B0[i, j] = new Vector(0, 0, 0);
                    }

                }
            });
        }

        private bool UpdateB()
        {
            bool isConverged = true;
            double delta;

            Parallel.For(1, _settings.ResolutionSpace + 1, i =>
            {
                for (int j = 1; j < _settings.ResolutionSpace + 1; j++)
                {
                    Vector oldValue = B[i, j];

                    if (isMaterial[i, j])
                    {
                        B[i, j] = (1 - _settings.Omega) * B[i, j]
                                     + _settings.Omega / (1 + 4 * Beta) * (BPrevious[i, j] - B0[i, j] + B0Previous[i, j] +
                                     Beta * (B[i + 1, j] + B[i - 1, j] + B[i, j + 1] + B[i, j - 1]
                                          + BPrevious[i + 1, j] + BPrevious[i - 1, j] + BPrevious[i, j + 1] + BPrevious[i, j - 1] - 4 * BPrevious[i, j]));
                    }
                    else
                    {
                        B[i, j] = new Vector(0, 0, 0);
                    }

                    delta = (B[i, j] - oldValue).Magnitude();

                    if (delta > _settings.Tolerance)
                    {
                        isConverged = false;
                    }
                }
            });

            return isConverged;
        }

        private void CalculateA()
        {
            Parallel.For(1, _settings.ResolutionSpace + 1, i =>
            {
                for (int j = 1; j < _settings.ResolutionSpace + 1; j++)
                {
                    if (isMaterial[i, j])
                    {
                        A[i, j] = 0.5 * _settings.MaterialConductivity / Dx * new Vector(B[i, j + 1].Z + B0[i, j + 1].Z - B[i, j - 1].Z - B0[i, j - 1].Z, B[i + 1, j].Z + B0[i + 1, j].Z - B[i - 1, j].Z - B0[i - 1, j].Z, 0);
                    }
                    else
                    {
                        A[i, j] = new Vector(0, 0, 0);
                    }
                }
            });
        }

        private void CalculateJ()
        {
            MultigridSolve(0);
        }

        private void MultigridSolve(int currentLevel)
        {
            bool isCon = false;
            
            if (currentLevel == maxLevels - 1)
            {
                while (!isCon)
                {
                    isCon = Relaxation(currentLevel);
                }
            }
            else
            {
                for (int i = 0; i < numPreRelaxations; i++)
                {
                    Relaxation(currentLevel);
                }

                Vector[,] residual = ComputeResidual(currentLevel);
                Vector[,] coarseResidual = Restrict(currentLevel);
                MultigridSolve(currentLevel + 1);
                Vector[,] fineCorrection = Prolongate(currentLevel + 1);
                ApplyCorrection(fineCorrection, currentLevel);

                for (int i = 0; i < numPostRelaxations; i++)
                {
                    Relaxation(currentLevel);
                }
            }
        }

        private bool Relaxation(int currentLevel)
        {
            int size = JGrids[currentLevel].GetLength(0);
            // ... Implement a few iterations of the Gauss-Seidel relaxation here for the specified grid level...
        }
        private Vector[,] ComputeResidual(int currentLevel)
        {
            // Compute the difference between the current solution and the desired one
            // This will involve looping over the grid and calculating the residual for each node
        }
        private Vector[,] Restrict(int currentLevel)
        {
            int coarseSize = JGrids[currentLevel + 1].GetLength(0) / 2;
            Vector[,] coarseResidual = new Vector[coarseSize, coarseSize];

            for (int i = 0; i < coarseSize; i++)
            {
                for (int j = 0; j < coarseSize; j++)
                {
                    int iFine = 2 * i;
                    int jFine = 2 * j;
                    coarseResidual[i, j] = 0.25 * (JGrids[currentLevel][iFine, jFine] + JGrids[currentLevel][iFine + 1, jFine] + JGrids[currentLevel][iFine, jFine + 1] + JGrids[currentLevel][iFine + 1, jFine + 1]);
                }
            }

            return coarseResidual;
        }

        private Vector[,] Prolongate(int currentLevel)
        {
            int fineSize = JGrids[currentLevel - 1].GetLength(0);
            Vector[,] fineSolution = new Vector[fineSize, fineSize];

            for (int i = 0; i < fineSize; i++)
            {
                for (int j = 0; j < fineSize; j++)
                {
                    int iCoarse = i / 2;
                    int jCoarse = j / 2;
                    fineSolution[i, j] = 0.25 * (JGrids[currentLevel][iCoarse, jCoarse] + JGrids[currentLevel][Math.Min(iCoarse + 1, JGrids[currentLevel].GetLength(0) - 1), jCoarse] + JGrids[currentLevel][iCoarse, Math.Min(jCoarse + 1, JGrids[currentLevel].GetLength(1) - 1)] + JGrids[currentLevel][Math.Min(iCoarse + 1, JGrids[currentLevel].GetLength(0) - 1), Math.Min(jCoarse + 1, JGrids[currentLevel].GetLength(1) - 1)]);
                }
            }

            return fineSolution;
        }

        private void ApplyCorrection(Vector[,] fineCorrection, int currentLevel)
        {
            for (int i = 0; i < JGrids[currentLevel].GetLength(0); i++)
            {
                for (int j = 0; j < JGrids[currentLevel].GetLength(0); j++)
                {
                    JGrids[currentLevel][i, j] += fineCorrection[i, j];
                }
            }
        }

        private void CalculateForce()
        {
            force = new Vector(0, 0, 0);
            for (int i = 1; i < _settings.ResolutionSpace + 1; i++)
            {
                for (int j = 1; j < _settings.ResolutionSpace + 1; j++)
                {
                    if (isMaterial[i, j])
                    {
                        force += Dx * Dx * Dx * J[i, j].Cross(B[i, j]);
                    }
                }
            }
            force *= 2;
        }

        private void DynamicsStep()
        {
            acceleration = (1.0 / Mass) * force;
            position += velocity * Dt;
            velocity += acceleration * Dt;
        }

        private void TimeStep()
        {
            Array.Copy(B, BPrevious, BPrevious.Length);
            Array.Copy(B0, B0Previous, B0Previous.Length);
            Array.Copy(J, JPrevious, JPrevious.Length);
            Array.Copy(A, APrevious, APrevious.Length);
        }
    }

}
