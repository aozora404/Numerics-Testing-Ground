﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Input;

namespace Eddy.NET
{
    internal class Shit : IDifferentialEquationSolver
    {
        private double Dx;
        private double Dt;
        private double Mass;

        private bool[,] isMaterial;

        private Vector[,] B;
        private Vector[,] BPrevious;
        private Vector[,] B0;
        private Vector[,] B0Previous;

        private Vector[,] currentDensity;
        private Vector[,] currentDensityPrevious;

        private Vector force;
        private Vector position, velocity, acceleration;

        private readonly SimulationSettings _settings;

        private List<PhysicsOut> SimulationOutput;

        public Shit(SimulationSettings settings)
        {
            _settings = settings;

            Dx = _settings.Length / _settings.ResolutionSpace;
            Dt = (_settings.TimeEnd - _settings.TimeStart) / _settings.ResolutionTime;
            Mass = _settings.MaterialDensity * (4 / 3) * Math.PI * _settings.BallRadius * _settings.BallRadius * _settings.BallRadius;

            B = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];
            B0 = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];

            currentDensity = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];
            currentDensityPrevious = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];

            isMaterial = new bool[_settings.ResolutionSpace, _settings.ResolutionSpace];

            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {
                    B[i, j] = new Vector(0, 0, 0);
                    B0[i, j] = new Vector(0, 0, 0);

                    currentDensity[i, j] = new Vector(0, 0, 0);
                    currentDensityPrevious[i, j] = new Vector(0, 0, 0);

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

            SimulationOutput = new List<PhysicsOut>();
        }

        public void Solve()
        {
            double currentTime = _settings.TimeStart;
            while (currentTime < _settings.TimeEnd && !Console.KeyAvailable)
            {
                EMStep();
                DynamicsStep(currentTime);
                TimeStep();
                currentTime += Dt;
                Console.WriteLine($"Position: {position:f2} Velocity: {velocity:f2} Time: {currentTime:f2} Force: {force:f2}");
            }
        }

        public void PrintResults()
        {
            CSVWriter.SaveToFile(SimulationOutput, $@"C:\temp\DirectFieldSolverLegacy\out{DateTime.Now.ToString("yyyyMMddHHmmss")}.csv");
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
            CalculateB();
        }

        private void CalculateCharge()
        {
            CalculateCurrentDensity();
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


                    if (s.X > 0)
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

        private void CalculateCurrentDensity()
        {
            bool isConverged = false;
            int iterationCount = 0;

            while (!isConverged && iterationCount < _settings.MaxIterations)
            {
                isConverged = UpdateJ();
                iterationCount++;
            }

            Console.Write($"J Converged in {iterationCount} iterations. ");
        }

        private bool UpdateB()
        {
            bool isConverged = true;
            double delta;

            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {
                    Vector oldValue = B[i, j];

                    if (isMaterial[i, j])
                    {
                        B[i, j] = (1 - _settings.Omega) * B[i, j]
                                     + _settings.Omega / 4 * ((B[Math.Min(i + 1, _settings.ResolutionSpace - 1), j] + B[Math.Max(i - 1, 0), j] + B[i, Math.Min(j + 1, _settings.ResolutionSpace - 1)] + B[i, Math.Max(j - 1, 0)])
                                                             + Dx / 2 * _settings.VacuumMagneticPermeability * new Vector(0, 0, (currentDensity[Math.Min(i + 1, _settings.ResolutionSpace - 1), j].Y - currentDensity[Math.Max(i - 1, 0), j].Y) - (currentDensity[i, Math.Min(j + 1, _settings.ResolutionSpace - 1)].X - currentDensity[i, Math.Max(j - 1, 0)].X)));
                    }
                    else
                    {
                        B[i, j] = (1 - _settings.Omega) * B[i, j]
                                     + _settings.Omega / 4 * ((B[Math.Min(i + 1, _settings.ResolutionSpace - 1), j] + B[Math.Max(i - 1, 0), j] + B[i, Math.Min(j + 1, _settings.ResolutionSpace - 1)] + B[i, Math.Max(j - 1, 0)]));
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

        private bool UpdateJ()
        {
            bool isConverged = true;
            double delta;

            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {
                    Vector oldValue = currentDensity[i, j];

                    if (isMaterial[i, j])
                    {
                        currentDensity[i, j] = (1 - _settings.Omega) * currentDensity[i, j]
                                     + _settings.Omega / 4 * ((currentDensity[Math.Min(i + 1, _settings.ResolutionSpace - 1), j] + currentDensity[Math.Max(i - 1, 0), j] + currentDensity[i, Math.Min(j + 1, _settings.ResolutionSpace - 1)] + currentDensity[i, Math.Max(j - 1, 0)])
                                                             + Dx / 2 * _settings.VacuumMagneticPermeability * new Vector(0, 0, (currentDensity[Math.Min(i + 1, _settings.ResolutionSpace - 1), j].Y - currentDensity[Math.Max(i - 1, 0), j].Y) - (currentDensity[i, Math.Min(j + 1, _settings.ResolutionSpace - 1)].X - currentDensity[i, Math.Max(j - 1, 0)].X)));
                    }
                    else
                    {
                        currentDensity[i, j] = (1 - _settings.Omega) * currentDensity[i, j]
                                     + _settings.Omega / 4 * ((currentDensity[Math.Min(i + 1, _settings.ResolutionSpace - 1), j] + currentDensity[Math.Max(i - 1, 0), j] + currentDensity[i, Math.Min(j + 1, _settings.ResolutionSpace - 1)] + currentDensity[i, Math.Max(j - 1, 0)]));
                    }

                    delta = (currentDensity[i, j] - oldValue).Magnitude();

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
            for (int i = 0; i < _settings.ResolutionSpace; i++)
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {
                    if (isMaterial[i, j])
                    {
                        force += Dx * Dx * Dx * currentDensity[i, j].Cross(B[i, j] + B0[i, j]);
                    }
                }
            }
            force *= 2;
        }

        private void DynamicsStep(double time)
        {
            SimulationOutput.Add(new PhysicsOut(position, velocity, force, time));
            acceleration = (1.0 / Mass) * force;
            position += velocity * Dt;
            velocity += acceleration * Dt;
        }

        private void TimeStep()
        {
            Array.Copy(currentDensity, currentDensityPrevious, currentDensityPrevious.Length);
        }
    }

}