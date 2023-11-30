using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Input;

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

        private List<PhysicsOut> SimulationOutput;

        public DirectFieldSolverMk2(SimulationSettings settings)
        {
            _settings = settings;

            Dx = _settings.Length / _settings.ResolutionSpace;
            Dt = (_settings.TimeEnd - _settings.TimeStart) / _settings.ResolutionTime;
            Mass = _settings.MaterialDensity * (4.0 / 3.0) * Math.PI * _settings.BallRadius * _settings.BallRadius * _settings.BallRadius;

            B = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];
            B0 = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];
            E = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];
            E0 = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];
            
            currentDensity = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];
            currentDensityPrevious = new Vector[_settings.ResolutionSpace, _settings.ResolutionSpace];
            chargeDensity = new double[_settings.ResolutionSpace, _settings.ResolutionSpace];
            chargeDensityPrevious = new double[_settings.ResolutionSpace, _settings.ResolutionSpace];

            isMaterial = new bool[_settings.ResolutionSpace, _settings.ResolutionSpace];

            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {
                    B[i, j] = new Vector(0, 0, 0);
                    B0[i, j] = new Vector(0, 0, 0);
                    E[i, j] = new Vector(0, 0, 0);
                    E0[i, j] = new Vector(0, 0, 0);
                    
                    currentDensity[i, j] = new Vector(0, 0, 0);
                    currentDensityPrevious[i, j] = new Vector(0, 0, 0);
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

            SimulationOutput = new List<PhysicsOut>();
        }

        public void Solve()
        {
            int steps = 0;
            double currentTime = _settings.TimeStart;
            
            while (currentTime < _settings.TimeEnd && !Console.KeyAvailable)
            {
                EMStep();
                DynamicsStep(currentTime);
                TimeStep();
                currentTime += Dt;
                Console.WriteLine($"Position: {position:f2} mm Velocity: {velocity:f2} mm/s Force: {force:f2} mN Time: {currentTime * 1000:f2} ms");
                steps++;
                if ( steps % 5 == 0 )
                {
                    //UpdateVisuals();
                    steps = 0;
                }
            }
        }

        public void PrintResults()
        {
            CSVWriter.SaveToFile(SimulationOutput, $@"C:\temp\outDirectFieldSolverMk2{DateTime.Now.ToString("yyyyMMddHHmmss")}.csv");
        }

        private void EMStep()
        {
            LorentzShift(velocity);
            CalculateField();
            CalculateCharge();
            LorentzShift(-velocity);
            CalculateForce();
        }

        private void LorentzShift(Vector velocity)
        {
            Vector EParallel, E0Parallel, BParallel, B0Parallel;
            Vector EPerpendicular, E0Perpendicular, BPerpendicular, B0Perpendicular;
            Vector JParallel, JPerpendicular;
            Vector unitVelocity = velocity / velocity.Magnitude();
            double gamma = 1.0 / Math.Sqrt(1 - (velocity.Dot(velocity)) / (_settings.c * _settings.c));
            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {
                    EParallel = E[i, j].Dot(unitVelocity) * unitVelocity;
                    E0Parallel = E0[i, j].Dot(unitVelocity) * unitVelocity;
                    BParallel = B[i, j].Dot(unitVelocity) * unitVelocity;
                    B0Parallel = B0[i, j].Dot(unitVelocity) * unitVelocity;
                    JParallel = currentDensity[i,j].Dot(unitVelocity) * unitVelocity;

                    EPerpendicular = E[i, j] - EParallel;
                    E0Perpendicular = E0[i, j] - E0Parallel;
                    BPerpendicular = B[i, j] - BParallel;
                    B0Perpendicular = B0[i, j] - B0Parallel;
                    JPerpendicular = currentDensity[i,j] - JParallel;

                    E[i, j] = EParallel + gamma * (EPerpendicular + velocity.Cross(BPerpendicular));
                    E0[i, j] = E0Parallel + gamma * (E0Perpendicular + velocity.Cross(B0Perpendicular));
                    B[i, j] = BParallel + gamma * (BPerpendicular - 1.0/(_settings.c * _settings.c) * velocity.Cross(EPerpendicular));
                    B0[i, j] = B0Parallel + gamma * (B0Perpendicular - 1.0 / (_settings.c * _settings.c) * velocity.Cross(E0Perpendicular));
                    //currentDensity[i, j] = JPerpendicular + gamma * (JParallel - chargeDensity[i, j] * velocity);
                    //chargeDensity[i, j] = gamma * (chargeDensity[i, j] - velocity.Dot(JParallel) / (_settings.c * _settings.c));
                }
            });
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

            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {
                    relPosition.X = (i - _settings.ResolutionSpace / 2) * Dx;
                    relPosition.Y = (j - _settings.ResolutionSpace / 2) * Dx;

                    s = relPosition + position - coilOrigin;

                    if (s.Magnitude() < _settings.CoilRadius)
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
            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {
                    E0[i, j] = new Vector(0, 0, 0);
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

            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {
                    Vector oldValue = B[i, j];

                    B[i, j] = (1 - _settings.Omega) * B[i, j]
                                     + _settings.Omega / 4 * ((B[Math.Min(i + 1, _settings.ResolutionSpace - 1), j] + B[Math.Max(i - 1, 0), j] + B[i, Math.Min(j + 1, _settings.ResolutionSpace - 1)] + B[i, Math.Max(j - 1, 0)]
                                                             + B0[Math.Min(i + 1, _settings.ResolutionSpace - 1), j] + B0[Math.Max(i - 1, 0), j] + B0[i, Math.Min(j + 1, _settings.ResolutionSpace - 1)] + B0[i, Math.Max(j - 1, 0)] - 4 * B0[i, j])
                                                             + Dx / 2 * _settings.VacuumMagneticPermeability * new Vector(0, 0, (currentDensity[Math.Min(i + 1, _settings.ResolutionSpace - 1), j].Y - currentDensity[Math.Max(i - 1, 0), j].Y) - (currentDensity[i, Math.Min(j + 1, _settings.ResolutionSpace - 1)].X - currentDensity[i, Math.Max(j - 1, 0)].X)));
                    
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

            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {
                    Vector oldValue = E[i, j];

                    E[i, j] = (1 - _settings.Omega) * E[i, j]
                              + _settings.Omega / 4 * ((E[Math.Min(i + 1, _settings.ResolutionSpace - 1), j] + E[Math.Max(i - 1, 0), j] + E[i, Math.Min(j + 1, _settings.ResolutionSpace - 1)] + E[i, Math.Max(j - 1, 0)]
                                                      + E0[Math.Min(i + 1, _settings.ResolutionSpace - 1), j] + E0[Math.Max(i - 1, 0), j] + E0[i, Math.Min(j + 1, _settings.ResolutionSpace - 1)] + E0[i, Math.Max(j - 1, 0)] - 4 * E0[i, j])
                                                      //- Dx / (2 * _settings.VacuumElectricPermittivity) * new Vector(chargeDensity[Math.Min(i + 1, _settings.ResolutionSpace - 1), j] - chargeDensity[Math.Max(i - 1, 0), j], chargeDensity[i, Math.Min(j + 1, _settings.ResolutionSpace - 1)] - chargeDensity[i, Math.Max(j - 1, 0)], 0)
                                                      - (Dx * Dx / Dt) * _settings.VacuumMagneticPermeability * (currentDensity[i, j] - currentDensityPrevious[i, j])
                                                      );

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
            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {
                    if (isMaterial[i,j])
                    {
                        currentDensity[i, j] = _settings.MaterialConductivity * (E[i, j] + E0[i, j]);
                    }
                    else
                    {
                        currentDensity[i, j] = new Vector(0, 0, 0);
                    }
                    
                }
            });
        }

        private void CalculateChargeDensity()
        {
            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {

                    chargeDensity[i, j] = chargeDensityPrevious[i, j] + (Dt / (2 * Dx)) * (currentDensity[Math.Min(i + 1, _settings.ResolutionSpace - 1), j].X - currentDensity[Math.Max(i - 1, 0), j].X + currentDensity[i, Math.Min(j + 1, _settings.ResolutionSpace - 1)].Y - currentDensity[i, Math.Max(j - 1, 0)].Y);
                    /*
                    if (isMaterial[i, j])
                    {
                        
                    }
                    else
                    {
                        chargeDensity[i, j] = 0;
                    }
                    */
                    
                }
            });
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
                        force += Dx * Dx * Dx * (currentDensity[i, j].Cross(B[i, j] + B0[i, j]));
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
            Array.Copy(chargeDensity, chargeDensityPrevious, chargeDensityPrevious.Length);
        }
    }

}
