using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ILGPU;
using ILGPU.Runtime;
using ILGPU.Runtime.Cuda;

namespace Eddy.NET
{
    internal class DirectFieldSolverILGPU : IDifferentialEquationSolver
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

        private List<PhysicsOut> SimulationOutput;
        private PhysicsOut output;

        public DirectFieldSolverILGPU(SimulationSettings settings)
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

            using var context = Context.CreateDefault();
            using var accelerator = context.GetPreferredDevice(preferCPU: false).CreateAccelerator(context);
            

            while (currentTime < _settings.TimeEnd && !Console.KeyAvailable)
            {
                EMStep();
                DynamicsStep();
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
            CSVWriter.SaveToFile(SimulationOutput, $@"C:\temp\DirectFieldSolverILGPU\out{DateTime.Now.ToString("yyyyMMddHHmmss")}.csv");
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
                    currentDensity[i, j] = JPerpendicular + gamma * (JParallel - chargeDensity[i, j] * velocity);
                    chargeDensity[i, j] = gamma * (chargeDensity[i, j] - velocity.Dot(JParallel) / (_settings.c * _settings.c));
                }
            });
        }

        private void CalculateField()
        {
            EmbedB0();
            EmbedE0();
            CalculateJefimenko();
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
            Parallel.For(0, _settings.ResolutionSpace, i =>
            {
                for (int j = 0; j < _settings.ResolutionSpace; j++)
                {
                    E0[i, j] = new Vector(0, 0, 0);
                }
            });
        }

        private void CalculateJefimenko()
        {
            var kernel = accelerator.LoadAutoGroupedStreamKernel<Index2D, ArrayView<double>>(JefimenkoKernel);
            
            using var bufferE = accelerator.Allocate2D<Vector>(_settings.ResolutionSpace, _settings.ResolutionSpace);
            using var bufferB = accelerator.Allocate2D<Vector>(_settings.ResolutionSpace, _settings.ResolutionSpace);
            using var bufferChargeDensity = accelerator.Allocate2D<double>(chargeDensity);
            using var bufferCurrentDensity = accelerator.Allocate2D<Vector>(currentDensity);
            
            kernel();
            
            E = bufferE.GetAsArray2D();
            B = bufferB.GetAsArray2D();
        }

        static private void JefimenkoKernel(
            Index2D index,
            int length,
            ArrayView2D<Vector> E,
            ArrayView2D<Vector> B,
            ArrayView2D<double> chargeDensity,
            ArrayView2D<Vector> currentDensity
        )
        {
            var x = index.X;
            var y = index.Y;
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
                    if (isMaterial[i, j])
                    {
                        chargeDensity[i, j] = chargeDensityPrevious[i, j] + (Dt / (2 * Dx)) * (currentDensity[Math.Min(i + 1, _settings.ResolutionSpace - 1), j].X - currentDensity[Math.Max(i - 1, 0), j].X + currentDensity[i, Math.Min(j + 1, _settings.ResolutionSpace - 1)].Y - currentDensity[i, Math.Max(j - 1, 0)].Y);
                    }
                    else
                    {
                        chargeDensity[i, j] = 0;
                    }
                    
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

        private void DynamicsStep()
        {
            SimulationOutput.Add(new PhysicsOut(position, velocity, force, time));
            acceleration = (1.0 / Mass) * force;
            position += velocity * Dt;
            velocity += acceleration * Dt;
        }

        private void TimeStep()
        {
            Array.Copy(chargeDensity, chargeDensityPrevious, chargeDensityPrevious.Length);
        }
    }
}