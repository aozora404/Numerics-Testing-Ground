using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Eddy.NET
{
    internal class DirectFieldSolverMk2 : IDifferentialEquationSolver
    {
        private double Dx;
        private double Dt;
        private double Mass;

        private bool[,] isMaterial;

        private Vector[,] B;
        private Vector[,] B0;
        private Vector[,] E;
        private Vector[,] E0;

        private Vector[,] currentDensity;
        private double[,] chargeDensity;
        private double[,] chargeDensityPrevious;

        private Vector force;
        private Vector position, velocity, acceleration;

        private readonly SimulationSettings _settings;

        public DirectFieldSolverMultigrid(SimulationSettings settings)
        {
            _settings = settings;

            Dx = _settings.Length / _settings.ResolutionSpace;
            Dt = (_settings.TimeEnd - _settings.TimeStart) / _settings.ResolutionTime;
            Mass = _settings.MaterialDensity * (4 / 3) * Math.PI * _settings.BallRadius * _settings.BallRadius * _settings.BallRadius;

            B = new Vector[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            B0 = new Vector[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            E = new Vector[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            E0 = new Vector[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            
            currentDensity = new Vector[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            chargeDensity = new double[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            chargeDensityPrevious = new double[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];

            isMaterial = new bool[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];

            Parallel.For(0, _settings.ResolutionSpace + 2, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace + 2; j++)
                {
                    B[i, j] = new Vector(0, 0, 0);
                    B0[i, j] = new Vector(0, 0, 0);
                    E[i, j] = new Vector(0, 0, 0);
                    E0[i, j] = new Vector(0, 0, 0);
                    
                    currentDensity[i, j] = new Vector(0, 0, 0);
                    chargeDensity[i, j] = 0;
                    chargeDensityPrevious[i, j] = 0;
                    
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
            CalculateField();
            CalculateCharge();
            CalculateForce();
        }

        private void CalculateField()
        {
            EmbedB0();
            EmbedE0();
            bool isConverged = false;
            int iterationCount = 0;

            while (!isConverged && iterationCount < _settings.MaxIterations)
            {
                isConverged = UpdateB();

                iterationCount++;
            }

            Console.Write($"B Converged in {iterationCount} iterations. ");

            isConverged = false;
            iterationCount = 0;

            while (!isConverged && iterationCount < _settings.MaxIterations)
            {
                isConverged = UpdateE();

                iterationCount++;
            }

            Console.Write($"E Converged in {iterationCount} iterations. ");
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

        private void EmbedE0()
        {
            Parallel.For(0, _settings.ResolutionSpace + 2, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace + 2; j++)
                {
                    E0[i,j] = velocity.Cross(B0[i,j])

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

                    B[i, j] = (1 - _settings.Omega) * B[i, j]
                                     + _settings.Omega / 4 * ((B[i + 1, j] + B[i - 1, j] + B[i, j + 1] + B[i, j - 1]
                                                             + B0[i + 1, j] + B0[i - 1, j] + B0[i, j + 1] + B0[i, j - 1] - 4 * B0[i, j])
                                                             + Dx/2 * _settings.MaterialMagneticPermeability * new Vector(0,0,(J[i+1,j].Y - J[i-1,j].Y) - (J[i,j+1].X-J[i,j-1].X)));
                    
                    delta = (B[i, j] - oldValue).Magnitude();

                    if (delta > _settings.Tolerance)
                    {
                        isConverged = false;
                    }
                }
            });

            return isConverged;
        }

        private bool UpdateE()
        {
            bool isConverged = true;
            double delta;

            Parallel.For(1, _settings.ResolutionSpace + 1, i =>
            {
                for (int j = 1; j < _settings.ResolutionSpace + 1; j++)
                {
                    Vector oldValue = E[i, j];

                    E[i, j] = (1 - _settings.Omega) * E[i, j]
                                     + _settings.Omega / 4 * ((E[i + 1, j] + E[i - 1, j] + E[i, j + 1] + E[i, j - 1]
                                                             + E0[i + 1, j] + E0[i - 1, j] + E0[i, j + 1] + E0[i, j - 1] - 4 * E0[i, j])
                                                             - Dx/2 * _settings.MaterialElectricPermittivity * new Vector(0,0,(J[i+1,j].Y - J[i-1,j].Y) - (J[i,j+1].X-J[i,j-1].X)));
                    
                    delta = (E[i, j] - oldValue).Magnitude();

                    if (delta > _settings.Tolerance)
                    {
                        isConverged = false;
                    }
                }
            });

            return isConverged;
        }

        private void CalculateJ()
        {
            MultigridSolve(0);
        }

        private void MultigridSolve(int currentLevel)
        {
            if (currentLevel == maxLevels - 1)
            {
                RelaxationToConvergence(currentLevel);
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
            Parallel.For(1, size + 1, i =>
            {
                for (int j = 1; j < size + 1; j++)
                {
                    if (isMaterial[i, j]) // TODO: Redo
                    {
                        JGrids[currentLevel][i, j] = (1 - omegaJ) * JGrids[currentLevel][i, j]
                                     + omegaJ / 4 * ((Dx * Dx / Dt) * (APrevious[i, j] - A[i, j]) + JGrids[currentLevel][i + 1, j] + JGrids[currentLevel][i - 1, j] + JGrids[currentLevel][i, j + 1] + JGrids[currentLevel][i, j - 1]);
                    }
                    else
                    {
                        JGrids[currentLevel][i, j] = new Vector(0, 0, 0);
                    }
                }
            });
        }

        private void RelaxationToConvergence(int currentLevel)
        {
            int size = JGrids[currentLevel].GetLength(0);
            bool isNotConverged = true;
            int iterCount = 0;
            while(isNotConverged || iterCount < 5000)
            {
                isNotConverged = false;
                Parallel.For(1, size + 1, i =>
                {
                    for (int j = 1; j < size + 1; j++)
                    {
                        Vector oldValue = JGrids[currentLevel][i, j];

                        if (isMaterial[i, j]) // TODO: Redo
                        {
                            JGrids[currentLevel][i, j] = (1 - omegaJ) * JGrids[currentLevel][i, j]
                                        + omegaJ / 4 * ((Dx * Dx / Dt) * (APrevious[i, j] - A[i, j]) + JGrids[currentLevel][i + 1, j] + JGrids[currentLevel][i - 1, j] + JGrids[currentLevel][i, j + 1] + JGrids[currentLevel][i, j - 1]);
                        }
                        else
                        {
                            JGrids[currentLevel][i, j] = new Vector(0, 0, 0);
                        }

                        delta = (JGrids[currentLevel][i, j] - oldValue).Magnitude();

                        if (delta > _settings.Tolerance)
                        {
                            isNotConverged = true;
                        }
                    }
                });
                iterCount++;
            }
        }

        private Vector[,] ComputeResidual(int currentLevel)
        {
            // Compute the difference between the current solution and the desired one
            // This will involve looping over the grid and calculating the residual for each node
            int size = JGrids[currentLevel].GetLength(0);
            Vector[,] res = new Vector[size,size];
            Parallel.For(1, size + 1, i =>
            {
                for (int j = 1; j < size + 1; j++)
                {
                    res[i,j] = new Vector(0,0,0);
                }
            });
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
