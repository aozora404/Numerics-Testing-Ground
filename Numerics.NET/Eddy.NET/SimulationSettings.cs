using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Eddy.NET
{
    internal class SimulationSettings
    {
        // mmkgs

        public double Length = 10.0;
        public int ResolutionSpace = 256;

        public double TimeStart = 0.0;
        public double TimeEnd = 1.0;
        public int ResolutionTime = 50000;

        public double Omega = 1.0;
        public double Tolerance = 10E-9;
        public int MaxIterations = 5000;

        public double InitialVelocity = 10.0;
        public double BallRadius = 5.0;
        public string Material = "Aluminium";
        public double MaterialConductivity = 3.5e1;
        public double MaterialMagneticPermeability = 1.256665e-3;
        public double MaterialElectricPermittivity = 1;
        public double MaterialDensity = 2.7e-6;

        public double CoilRadius = 10.0;
        public double CoilMagneticFieldStrength = 0.01;
        public double CoilDistance = 15.1;
    }
}
