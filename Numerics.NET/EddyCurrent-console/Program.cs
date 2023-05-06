﻿using NumericsEngine;
using NLib = NumericsEngine.NumericsLib;


// Settings

// Simulation Settings
double gridSize = 20.0;                             // mm
double cellSize = 0.1;                              // mm
double timeStep = 0.00004;                           // s
double timeEnd = 1.0;                               // s

// Object Settings
double initialVelocity = 10.0;                      // mm s-1
double mass = 0.01;                                 // kg
double objectRadius = 8.0;                          // mm
double materialConductivity = 3.5e1;                // mm-2 s3 kg-1 A2
double materialMagneticPermeability = 1.256665e-3;  // mm s-2 kg A-2

// Electromagnet Coil Settings
double coilRadius = 200.0;
double coilMagneticFieldStrength = -0.1;
double coilDistance = 208.0;

Console.WriteLine("Current Settings");
Console.WriteLine("----- Simulation Settings -----");
Console.WriteLine("Grid Size = {0} mm", gridSize);
Console.WriteLine("Cell Size = {0} mm", cellSize);
Console.WriteLine("Time Step = {0} s", timeStep);
Console.WriteLine("Time End = {0} s", timeEnd);
Console.WriteLine("----- Object Settings -----");
Console.WriteLine("Initial Velocity = {0} mm/s", initialVelocity);
Console.WriteLine("Mass = {0} kg", mass);
Console.WriteLine("Radius = {0} mm", objectRadius);
Console.WriteLine("Conductivity = {0:E} S/m", materialConductivity * 1e6);
Console.WriteLine("Permeability = {0:E} H/m", materialMagneticPermeability * 1e-3);
Console.WriteLine("----- Coil Settings -----");
Console.WriteLine("Radius = {0} mm", coilRadius);
Console.WriteLine("Strength = {0} T", coilMagneticFieldStrength);
Console.WriteLine("Distance = {0} mm", coilDistance);
Console.WriteLine("Press Enter to continue...");
Console.ReadLine();

/*************************** SIMULATION ********************************/

// Pre-processing
// Constants and Coefficients
Console.WriteLine("Calculating coefficients...");

double permeability = materialMagneticPermeability;
double permeabilityInverse = 1 / permeability;

double massInverse = 1 / mass;
double objectRadiusSquared = objectRadius * objectRadius;

double coilRadiusSquared = coilRadius * coilRadius;
Vector2D coilOrigin = new(coilDistance, 0);

double alpha = 0.5 * coilMagneticFieldStrength;

double cellArea = cellSize * cellSize;
double cellSizeInverse = 1 / cellSize;
int cellCount = (int)(gridSize / cellSize);

// Dummy Variable
Vector2D s;

// Initialize Structures
Console.WriteLine("Building structures...");

Vector2D[,] relativePosition = new Vector2D[cellCount + 2, cellCount + 2];

Vector2D[,] current = new Vector2D[cellCount + 2, cellCount + 2];
double[,] conductivity = new double[cellCount + 2, cellCount + 2];

double[,] magneticField = new double[cellCount + 2, cellCount + 2];


Vector2D[,] magneticPotential = new Vector2D[cellCount + 2, cellCount + 2];
Vector2D[,] magneticPotentialNext = new Vector2D[cellCount + 2, cellCount + 2];

Vector2D[,] magneticPotentialCurrentX = new Vector2D[cellCount + 2, cellCount + 2];
Vector2D[,] magneticPotentialCurrentXNext = new Vector2D[cellCount + 2, cellCount + 2];

Vector2D[,] magneticPotentialCurrentY = new Vector2D[cellCount + 2, cellCount + 2];
Vector2D[,] magneticPotentialCurrentYNext = new Vector2D[cellCount + 2, cellCount + 2];

Vector2D[,] magneticPotentialDelta = new Vector2D[cellCount + 2, cellCount + 2];
Vector2D[,] magneticPotentialDeltaNext = new Vector2D[cellCount + 2, cellCount + 2];


Vector2D[,] magneticPotentialEmbedded = new Vector2D[cellCount + 2, cellCount + 2];
Vector2D[,] magneticPotentialEmbeddedPrevious = new Vector2D[cellCount + 2, cellCount + 2];

Vector2D[,] magneticPotentialEmbeddedDelta = new Vector2D[cellCount + 2, cellCount + 2];
Vector2D[,] magneticPotentialEmbeddedLagrangian = new Vector2D[cellCount + 2, cellCount + 2];

Vector2D[,] magneticPotentialTotal = new Vector2D[cellCount + 2, cellCount + 2];

for (int y = 0; y < cellCount + 2; y++)
{
    for (int x = 0; x < cellCount + 2; x++)
    {
        relativePosition[y, x] = new Vector2D(0, 0);
        current[y, x] = new Vector2D(0, 0);
        magneticPotential[y, x] = new Vector2D(0, 0);
        magneticPotentialNext[y, x] = new Vector2D(0, 0);
        magneticPotentialCurrentX[y, x] = new Vector2D(0, 0);
        magneticPotentialCurrentXNext[y, x] = new Vector2D(0, 0);
        magneticPotentialCurrentY[y, x] = new Vector2D(0, 0);
        magneticPotentialCurrentYNext[y, x] = new Vector2D(0, 0);
        magneticPotentialDelta[y, x] = new Vector2D(0, 0);
        magneticPotentialDeltaNext[y, x] = new Vector2D(0, 0);
        magneticPotentialEmbedded[y, x] = new Vector2D(0, 0);
        magneticPotentialEmbeddedPrevious[y, x] = new Vector2D(0, 0);
        magneticPotentialEmbeddedDelta[y, x] = new Vector2D(0, 0);
        magneticPotentialEmbeddedLagrangian[y, x] = new Vector2D(0, 0);
        magneticPotentialTotal[y, x] = new Vector2D(0, 0);
    }
}

Vector2D force = new(0, 0);
Vector2D velocity = new(0, 0);
Vector2D velocityNext = new(0, 0);
Vector2D position = new(0, 0);
Vector2D positionNext = new(0, 0);

List<double> dataTime = new();
List<double> dataPosition = new();
List<double> dataVelocity = new();
List<double> dataForce = new();


// Set Initial Values
Console.WriteLine("Initializing field values...");

velocity.SetValue(initialVelocity, 0);

for (int y = 1; y <= cellCount; y++)
{
    for (int x = 1; x <= cellCount; x++)
    {
        relativePosition[y, x].SetValue(-gridSize / 2 + cellSize * x, -gridSize / 2 + cellSize * y);
    }
}

for (int y = 1; y <= cellCount; y++)
{
    for (int x = 1; x <= cellCount; x++)
    {
        if (relativePosition[y, x].MagnitudeSquared() < objectRadiusSquared)
        {
            conductivity[y, x] = materialConductivity;
        }
    }
}


for (int y = 1; y <= cellCount; y++)
{
    for (int x = 1; x <= cellCount; x++)
    {
        s = relativePosition[y, x] + position - coilOrigin;
        if (s.MagnitudeSquared() > coilRadiusSquared)
        {
            magneticPotentialEmbeddedPrevious[y, x] = alpha * (coilRadiusSquared / s.MagnitudeSquared()) * new Vector2D(-s.y, s.x);
        }
        else
        {
            magneticPotentialEmbeddedPrevious[y, x] = alpha * new Vector2D(-s.y, s.x);
        }
    }
}


// Solver
Console.WriteLine("Calculating...");

double t = 0;
double dt = timeStep;

while (t < timeEnd)
{
    // Calculate embedded magnetic potential
    for (int y = 1; y <= cellCount; y++)
    {
        for (int x = 1; x <= cellCount; x++)
        {
            s = relativePosition[y, x] + position - coilOrigin;
            double sMagnitude = s.MagnitudeSquared();
            if (s.MagnitudeSquared() > coilRadiusSquared)
            {
                magneticPotentialEmbedded[y, x] = alpha * (coilRadiusSquared / sMagnitude) * new Vector2D(-s.y, s.x);
                magneticPotentialEmbeddedDelta[y, x] = (1 / dt) * (magneticPotentialEmbedded[y, x] - magneticPotentialEmbeddedPrevious[y, x]);
                magneticPotentialEmbeddedLagrangian[y, x] = (1 / sMagnitude) * magneticPotentialEmbedded[y, x];
            }
            else
            {
                magneticPotentialEmbedded[y, x] = alpha * new Vector2D(-s.y, s.x);
                magneticPotentialEmbeddedDelta[y, x] = alpha * (new Vector2D(-velocity.y, velocity.x));
                magneticPotentialEmbeddedLagrangian[y, x] = (1 / sMagnitude) * magneticPotentialEmbedded[y, x];
            }
        }
    }

    // Induced magnetic potential
    for (int y = 1; y <= cellCount; y++)
    {
        for (int x = 1; x <= cellCount; x++)
        {
            if (conductivity[y, x] != 0)
            {
                magneticPotentialDeltaNext[y, x] = (permeabilityInverse / conductivity[y, x]) * (cellSizeInverse * new Vector2D(NLib.DivergenceFlux2D(magneticPotentialCurrentX, x, y), NLib.DivergenceFlux2D(magneticPotentialCurrentY, x, y))) + magneticPotentialEmbeddedLagrangian[y, x] - magneticPotentialEmbeddedDelta[y, x];
                magneticPotentialCurrentXNext[y, x] = cellSizeInverse * NLib.GradientFlux2D(NLib.GetXComponent(magneticPotential), x, y);
                magneticPotentialCurrentYNext[y, x] = cellSizeInverse * NLib.GradientFlux2D(NLib.GetYComponent(magneticPotential), x, y);
                magneticPotentialNext[y, x] = magneticPotential[y, x] + dt * magneticPotentialDelta[y, x];
            }
        }
    }

    // Total magnetic potential
    magneticPotentialTotal = NLib.AddVectorFields(magneticPotential, magneticPotentialEmbedded);

    // Current and magnetic field
    for (int y = 1; y <= cellCount; y++)
    {
        for (int x = 1; x <= cellCount; x++)
        {
            if (conductivity[y, x] != 0)
            {
                current[y, x] = -conductivity[y, x] * (magneticPotentialDelta[y, x] + magneticPotentialEmbeddedDelta[y, x]);
                magneticField[y, x] = cellSizeInverse * NLib.CurlFlux2D(magneticPotentialTotal, x, y);
            }
        }
    }

    // Dynamics
    force.SetValue(0, 0);
    for (int y = 1; y <= cellCount; y++)
    {
        for (int x = 1; x <= cellCount; x++)
        {
            if (conductivity[y, x] != 0)
            {
                force += magneticField[y, x] * new Vector2D(current[y, x].y, -current[y, x].x);
            }
        }
    }

    force = cellArea * force;
    velocityNext = velocity + dt * massInverse * force;
    positionNext = position + dt * velocity;

    // Time step
    magneticPotential = magneticPotentialNext;
    magneticPotentialCurrentX = magneticPotentialCurrentXNext;
    magneticPotentialCurrentY = magneticPotentialCurrentYNext;
    magneticPotentialDelta = magneticPotentialDeltaNext;
    magneticPotentialEmbeddedPrevious = magneticPotentialEmbedded;
    position = positionNext;
    velocity = velocityNext;

    t += dt;

    dataTime.Add(t);
    dataPosition.Add(position.x);
    dataVelocity.Add(velocity.x);
    dataForce.Add(force.x);

    Console.WriteLine("{0} s", t);
    Console.WriteLine("Position = ({0}, {1}) mm", position.x, position.y);
    Console.WriteLine("Velocity = ({0}, {1}) mm", velocity.x, velocity.y);
    Console.WriteLine("Force = ({0}, {1}) mm", force.x, force.y);
    Console.WriteLine();
}

// Post-processing

Console.WriteLine("Final Position = ({0}, {1}) mm", position.x, position.y);
Console.WriteLine("Final Velocity = ({0}, {1}) mm", velocity.x, velocity.y);