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
        private Vector[,] currentDensityPrevious;
        private double[,] chargeDensity;
        private double[,] chargeDensityPrevious;

        private Vector force;
        private Vector position, velocity, acceleration;

        private readonly SimulationSettings _settings;

        public DirectFieldSolverMk2(SimulationSettings settings)
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
            currentDensityPrevious = new Vector[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
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
            CalculateB();
            CalculateE();
        }

        private void CalculateCharge()
        {
            CalculateCurrentDensity(); 
            CalculateChargeDensity();
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
                    E0[i, j] = velocity.Cross(B0[i, j]);
                }
            });
        }

        private void CalculateB()
        {
            bool isConverged = false;
            int iterationCount = 0;

            while (!isConverged && iterationCount < _settings.MaxIterations)
            {
                isConverged = UpdateB();
                iterationCount++;
            }

            Console.Write($"B Converged in {iterationCount} iterations. ");
        }

        private void CalculateE()
        {
            bool isConverged = false;
            int iterationCount = 0;

            while (!isConverged && iterationCount < _settings.MaxIterations)
            {
                isConverged = UpdateE();
                iterationCount++;
            }

            Console.Write($"E Converged in {iterationCount} iterations. ");
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
                                                             + Dx/2 * _settings.MaterialMagneticPermeability * new Vector(0,0,(currentDensity[i+1,j].Y - currentDensity[i-1,j].Y) - (currentDensity[i,j+1].X-currentDensity[i, j - 1].X)));
                    
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
                                                             - Dx/2 * _settings.MaterialElectricPermittivity * new Vector(chargeDensity[i+1,j] - chargeDensity[i-1,j], chargeDensity[i,j+1] - chargeDensity[i,j-1], 0)
                                                             - (Dx * Dx / Dt) * _settings.MaterialMagneticPermeability * (currentDensity[i,j] - currentDensityPrevious[i,j]));
                    
                    delta = (E[i, j] - oldValue).Magnitude();

                    if (delta > _settings.Tolerance)
                    {
                        isConverged = false;
                    }
                }
            });

            return isConverged;
        }

        private void CalculateCurrentDensity()
        {
            Parallel.For(1, _settings.ResolutionSpace + 1, i =>
            {
                for (int j = 1; j < _settings.ResolutionSpace + 1; j++)
                {
                    currentDensity[i,j] = _settings.MaterialConductivity * (E[i,j] + E0[i,j])
                }
            });
        }

        private void CalculateChargeDensity()
        {
            Parallel.For(1, _settings.ResolutionSpace + 1, i =>
            {
                for (int j = 1; j < _settings.ResolutionSpace + 1; j++)
                {
                    chargeDensity[i,j] = chargeDensityPrevious[i,j] + (Dt/(2 * Dx)) * (currentDensity[i+1, j].X - currentDensity[i-1, j].X + currentDensity[i, j+1].Y - currentDensity[i, j-1].Y)
                }
            });
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
                        force += Dx * Dx * Dx * currentDensity[i, j].Cross(B[i, j]);
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
            Array.Copy(currentDensity, currentDensityPrevious, currentDensityPrevious.Length);
            Array.Copy(chargeDensity, chargeDensityPrevious, chargeDensityPrevious.Length);
        }
    }

}
