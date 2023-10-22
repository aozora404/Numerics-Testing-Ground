using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Eddy.NET
{
    internal class SkinEffectSolver : IDifferentialEquationSolver
    {
        private double Dx;
        private double Dt;
        private double Mass;

        private double resistance;
        private double skinDepth;
        private double omega;
        private double power;
        

        private Vector2D force;
        private Vector2D position, velocity, acceleration;

        private readonly SimulationSettings _settings;

        public SkinEffectSolver(SimulationSettings settings)
        {
            _settings = settings;

            Dx = _settings.Length / _settings.ResolutionSpace;
            Dt = (_settings.TimeEnd - _settings.TimeStart) / _settings.ResolutionTime;

            Mass = _settings.MaterialDensity * (4 / 3) * Math.PI * _settings.BallRadius * _settings.BallRadius * _settings.BallRadius;
            
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
            // Check for near-zero velocity
            if (Math.Abs(velocity.X) < 1e-6)
            {
                power = 0;
                return;
            }

            omega = velocity.X / (_settings.BallRadius * 2);
            skinDepth = Math.Sqrt(2 / (omega * _settings.MaterialMagneticPermeability * _settings.MaterialConductivity));
            // Bound the skin depth
            skinDepth = Math.Min(skinDepth, _settings.BallRadius * 0.95);

            double denom = (_settings.BallRadius * (_settings.BallRadius - skinDepth));
            // Ensure the denominator isn't too small
            if (Math.Abs(denom) < 1e-6)
            {
                power = 0;
                return;
            }

            resistance = 1.0 / (4 * Math.PI * _settings.MaterialConductivity) * (skinDepth / denom);
            power = _settings.CoilMagneticFieldStrength * velocity.X * Math.PI * _settings.BallRadius * _settings.BallRadius;
            power *= power;
            // Ensure resistance isn't too small
            if (resistance < 1e-6)
            {
                power = 0;
            }
            else
            {
                power /= resistance;
            }

            CalculateForce();
        }



        private void CalculateForce()
        {
            force.X = -power / velocity.X;
        }

        private void DynamicsStep()
        {
            acceleration = (1.0 / Mass) * force;
            position += velocity * Dt;
            velocity += acceleration * Dt;
        }

    }

}
