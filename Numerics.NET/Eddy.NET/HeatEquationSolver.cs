using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Eddy.NET
{
    internal class HeatEquationSolver : IDifferentialEquationSolver
    {
        private double[,] u;
        private double[,] uPrevious;
        private readonly SimulationSettings _settings;

        public HeatEquationSolver(SimulationSettings settings)
        {
            _settings = settings;
            Initialize();
        }

        public void Initialize()
        {
            u = new double[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            uPrevious = new double[_settings.ResolutionSpace + 2, _settings.ResolutionSpace + 2];
            // ... Initialize the arrays based on the simulation settings
        }

        public void Solve()
        {
            // ... Logic to solve the heat equation
        }

        public void PrintResults()
        {
            // ... Logic to print the results of the simulation
        }
    }

}
