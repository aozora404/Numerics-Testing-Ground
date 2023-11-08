﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Eddy.NET
{
    internal interface IDifferentialEquationSolver
    {
        void Solve();
        void PrintResults();
    }
}