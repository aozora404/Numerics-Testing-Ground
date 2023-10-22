using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Eddy.NET
{
    internal class EddyCurrentSolver : IDifferentialEquationSolver
    {
        private double Dx;
        private double Dt;
        private double Alpha;
        private double Beta;
        private double Mass;

        private bool[,] isMaterial;

        private Vector2D[,] A;
        private Vector2D[,] APrevious;
        private Vector2D[,] A0;
        private Vector2D[,] A0Previous;

        private Vector2D force;
        private Vector2D position, velocity, acceleration;

        private readonly SimulationSettings _settings;

        public EddyCurrentSolver(SimulationSettings settings)
        {
            _settings = settings;

            Dx = _settings.Length / _settings.ResolutionSpace;
            Dt = (_settings.TimeEnd - _settings.TimeStart) / _settings.ResolutionTime;
            Alpha = 1.0 / (_settings.MaterialConductivity * _settings.MaterialMagneticPermeability);
            Beta = Alpha * Dt / (2 * Dx * Dx);
            Mass = _settings.MaterialDensity * (4 / 3) * Math.PI * _settings.BallRadius * _settings.BallRadius * _settings.BallRadius;

            A = new Vector2D[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            APrevious = new Vector2D[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            A0 = new Vector2D[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            A0Previous = new Vector2D[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];

            isMaterial = new bool[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];

            Parallel.For(0, _settings.ResolutionSpace + 2, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace + 2; j++)
                {
                    A[i, j] = new Vector2D(0, 0);
                    APrevious[i, j] = new Vector2D(0, 0);
                    A0[i, j] = new Vector2D(0, 0);
                    A0Previous[i, j] = new Vector2D(0, 0);

                    double x = (i * 1.0 / _settings.ResolutionSpace) - 0.5;
                    double y = (j * 1.0 / _settings.ResolutionSpace) - 0.5;
                    double r = _settings.BallRadius / _settings.Length;
                    r = r * r;

                    if (x*x + y*y <= r)
                    {
                        isMaterial[i,j] = true;
                    } else
                    {
                        isMaterial[i, j] = false;
                    }
                }
            });

            force = new Vector2D(0, 0);
            position = new Vector2D(0, 0);
            velocity = new Vector2D(_settings.InitialVelocity, 0);
            acceleration = new Vector2D(0, 0);

            
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
                Console.WriteLine($"Position: {position} Velocity: {velocity} Time: {currentTime} ");
            }
        }

        public void PrintResults()
        {
            // ... Logic to print the results of the simulation
        }

        private void EMStep()
        {
            EmbedA0();

            Array.Copy(APrevious, A, APrevious.Length);

            bool isConverged = false;
            int iterationCount = 0;

            while (!isConverged && iterationCount < _settings.MaxIterations)
            {
                isConverged = UpdateA();

                iterationCount++;
            }

            Console.Write($"Converged in {iterationCount} iterations. ");

            CalculateForce();
        }

        private void EmbedA0()
        {
            Vector2D relPosition = new Vector2D(0, 0);
            Vector2D s = new Vector2D(0, 0);
            Vector2D coilOrigin = new Vector2D(_settings.CoilDistance,0);

            Parallel.For(0, _settings.ResolutionSpace + 2, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace + 2; j++)
                {
                    relPosition.X = (i - _settings.ResolutionSpace/2) * Dx;
                    relPosition.Y = (j - _settings.ResolutionSpace/2) * Dx;

                    s = relPosition + position - coilOrigin;

                    
                    if(s.Magnitude() > _settings.CoilRadius)
                    {
                        A0[i, j] = 0.5 * _settings.CoilMagneticFieldStrength * (_settings.CoilRadius * _settings.CoilRadius / (s.Magnitude() * s.Magnitude())) * new Vector2D(-s.Y, s.X);
                    }
                    else
                    {
                        A0[i, j] = 0.5 * _settings.CoilMagneticFieldStrength * new Vector2D(-s.Y, s.X);
                    }

                }
            });
        }

        private bool UpdateA()
        {
            bool isConverged = true;
            double delta;

            Parallel.For(1, _settings.ResolutionSpace + 1, i =>
            {
                for (int j = 1; j < _settings.ResolutionSpace + 1; j++)
                {
                    Vector2D oldValue = A[i, j];

                    if (isMaterial[i, j])
                    {
                        A[i, j] = (1 - _settings.Omega) * A[i, j]
                                     + _settings.Omega / (1 + 4 * Beta) * (APrevious[i, j] - A0[i, j] + A0Previous[i, j] +
                                     Beta * (A[i + 1, j] + A[i - 1, j] + A[i, j + 1] + A[i, j - 1]
                                          + APrevious[i + 1, j] + APrevious[i - 1, j] + APrevious[i, j + 1] + APrevious[i, j - 1] - 4 * APrevious[i, j]));
                    }
                    else
                    {
                        /*
                        A[i, j] = (1 - _settings.Omega) * A[i, j]
                                     + (_settings.Omega / 8) * (A[i + 1, j] + A[i - 1, j] + A[i, j + 1] + A[i, j - 1]
                                                             + APrevious[i + 1, j] + APrevious[i - 1, j] + APrevious[i, j + 1] + APrevious[i, j - 1] - 4 * APrevious[i, j]);
                        */
                        A[i, j] = new Vector2D(0, 0);
                    }

                    delta = (A[i, j] - oldValue).Magnitude();

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
            Vector2D F = new Vector2D(0, 0);
            Vector2D J = new Vector2D(0, 0);
            double B;
            force = new Vector2D(0, 0);
            for (int i = 1; i < _settings.ResolutionSpace + 1; i++)
            {
                for (int j = 1; j < _settings.ResolutionSpace + 1; j++)
                {
                    if (isMaterial[i, j])
                    {
                        B = A[i + 1, j].Y - A[i - 1, j].Y + A0[i + 1, j].Y - A0[i - 1, j].Y - A[i, j + 1].X + A[i, j - 1].X - A0[i, j + 1].X + A0[i, j - 1].X;
                        J.X = -A[i, j].Y + APrevious[i, j].Y - A0[i, j].Y + A0Previous[i, j].Y;
                        J.Y = A[i, j].X - APrevious[i, j].X + A0[i, j].X - A0Previous[i, j].X;
                        F = _settings.MaterialConductivity * Dx * Dx / (2 * Dt) * J * B;
                        force += F;
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
            Array.Copy(A0, A0Previous, A0Previous.Length);
            Array.Copy(A, APrevious, APrevious.Length);
        }
    }

}
