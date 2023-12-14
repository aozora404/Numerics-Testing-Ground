using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Eddy.NET
{
    internal class DirectFieldSolver : IDifferentialEquationSolver
    {
        private double Dx;
        private double Dt;
        private double Alpha;
        private double Beta;
        private double Mass;
        private double omegaJ;

        private bool[,] isMaterial;

        private Vector[,] B;
        private Vector[,] BPrevious;
        private Vector[,] B0;
        private Vector[,] B0Previous;
        private Vector[,] J;
        private Vector[,] JPrevious;
        private Vector[,] A;
        private Vector[,] APrevious;


        private Vector force;
        private Vector position, velocity, acceleration;

        private readonly SimulationSettings _settings;

        public DirectFieldSolver(SimulationSettings settings)
        {
            _settings = settings;

            Dx = _settings.Length / _settings.ResolutionSpace;
            Dt = (_settings.TimeEnd - _settings.TimeStart) / _settings.ResolutionTime;
            Alpha = 1.0 / (_settings.MaterialConductivity * _settings.MaterialMagneticPermeability);
            Beta = Alpha * Dt / (2 * Dx * Dx);
            Mass = _settings.MaterialDensity * (4 / 3) * Math.PI * _settings.BallRadius * _settings.BallRadius * _settings.BallRadius;
            omegaJ = 1.9;

            B = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];
            BPrevious = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];
            B0 = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];
            B0Previous = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];
            J = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];
            JPrevious = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];
            A = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];
            APrevious = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];


            isMaterial = new bool[_settings.ResolutionSpace, _settings.ResolutionSpace];

            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
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

            isConverged = false;
            iterationCount = 0;

            while (!isConverged && iterationCount < _settings.MaxIterations)
            {
                isConverged = CalculateJ();

                iterationCount++;
            }
            

            Console.Write($"J Converged in {iterationCount} iterations. ");

            CalculateForce();
        }

        private void EmbedB0()
        {
            Vector relPosition = new Vector(0, 0, 0);
            Vector s = new Vector(0, 0, 0);
            Vector coilOrigin = new Vector(_settings.CoilDistance, 0, 0);

            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {
                    relPosition.X = (i - _settings.ResolutionSpace / 2) * Dx;
                    relPosition.Y = (j - _settings.ResolutionSpace / 2) * Dx;

                    s = relPosition + position - coilOrigin;


                    if (s.Magnitude() > _settings.CoilRadius)
                    {
                        B0[i, j] =  new Vector(0, 0, _settings.CoilMagneticFieldStrength);
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

        private bool CalculateJ()
        {
            bool isConverged = true;
            double delta;

            Parallel.For(1, _settings.ResolutionSpace + 1, i =>
            {
                for (int j = 1; j < _settings.ResolutionSpace + 1; j++)
                {
                    Vector oldValue = J[i, j];

                    if (isMaterial[i, j])
                    {
                        J[i, j] = (1 - omegaJ) * J[i, j]
                                     + omegaJ / 4 * ((Dx * Dx / Dt) * (APrevious[i, j] - A[i, j]) + J[i + 1, j] + J[i - 1, j] + J[i, j + 1] + J[i, j - 1]);
                    }
                    else
                    {
                        J[i, j] = new Vector(0, 0, 0);
                    }

                    delta = (J[i, j] - oldValue).Magnitude();

                    if (delta > _settings.Tolerance)
                    {
                        isConverged = false;
                    }
                }
            });

            return isConverged;
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
