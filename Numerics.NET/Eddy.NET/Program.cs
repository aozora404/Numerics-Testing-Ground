using Eddy.NET;
using Sharprompt;
using System.Runtime;

public class Program
{
    public static void Main()
    {
        SimulationSettings settings = new SimulationSettings();

        Console.Write($@"Settings:
----- Simulation Settings -----
Grid Size = {settings.Length} mm
Resolution = {settings.ResolutionSpace} 
Time Start = {settings.TimeStart} s
Time End = {settings.TimeEnd} s
Time Step = {(settings.TimeEnd - settings.TimeStart) / settings.ResolutionTime} s
----- Solver Settings -----
Omega = {settings.Omega}
Tolerance = {settings.Tolerance}
Iteration Limit = {settings.MaxIterations}
----- Object Settings -----
Material: {settings.Material}
Mass = {settings.MaterialDensity * (4 / 3) * Math.PI * settings.BallRadius * settings.BallRadius * settings.BallRadius} kg
Radius = {settings.BallRadius} mm
Initial Velocity = {settings.InitialVelocity} mm/s
----- Coil Settings -----
Radius = {settings.CoilRadius} mm
Strength = {settings.CoilMagneticFieldStrength} T
Distance = {settings.CoilDistance} mm


        ");
        bool proceed = Prompt.Confirm("Proceed with these settings?", defaultValue: false);

        if (proceed)
        {
            IDifferentialEquationSolver solver;

            var simulationType = Prompt.Select("Choose a simulation method", new[] { "Potential", "Direct Field" });


            if (simulationType == "Potential")
            {
                solver = new EddyCurrentSolver(settings);
            }
            else
            {
                solver = new DirectFieldSolver(settings);
            }

            solver.Solve();
            solver.PrintResults();
        }
        else
        {
            Console.WriteLine("Simulation aborted.");
        }
    }
}
