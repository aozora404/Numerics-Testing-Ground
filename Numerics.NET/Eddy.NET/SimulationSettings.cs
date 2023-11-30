﻿using System;
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
        public double TimeEnd = 0.1;
        public int ResolutionTime = 100000;

        public double Omega = 1.95;
        public double Tolerance = 10E-3;
        public int MaxIterations = 5000;

        public double VacuumMagneticPermeability = 1.256637e-3;
        public double VacuumElectricPermittivity = 8.854187e-21;
        public double c = 3e11;

        public double InitialVelocity = 250.0;
        public double BallRadius = 4.95;
        public string Material = "Aluminium";
        public double MaterialConductivity = 3.5e1;
        public double MaterialMagneticPermeability = 1.256665e-3;
        public double MaterialElectricPermittivity = 1000;
        public double MaterialDensity = 2.7e-6;

        public double CoilRadius = 10.0;
        public double CoilMagneticFieldStrength = 1.0;
        public double CoilDistance = 10.0;
    }
}
