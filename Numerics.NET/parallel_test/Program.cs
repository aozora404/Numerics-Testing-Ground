using System;
using System.Threading.Tasks;
// Settings

double length = 1.0;
int resolutionSpace = 32;
double dx = length / resolutionSpace;

double timeStart = 0.0;
double timeEnd = 1.0;
int resolutionTime = 1000;
double t;
double dt = (timeEnd - timeStart) / resolutionTime;

double omega = 1.5;
double tolerance = 10E-9;

double alpha = 0.1;
double beta = alpha * dt / (2 * dx * dx);

double[,] u = new double[resolutionSpace + 2, resolutionSpace + 2];
double[,] uPrevious = new double[resolutionSpace + 2, resolutionSpace + 2];


// Dummy
double[,] delta = new double[resolutionSpace, resolutionSpace];
int iter;
bool finish;


// Init

for (int i = 0; i < resolutionSpace + 2; i++)
{
    for (int j = 0; j < resolutionSpace + 2; j++)
    {
        u[i, j] = 0;
        if (i > resolutionSpace / 4 && i < 3 * resolutionSpace / 4 && j > resolutionSpace / 4 && j < 3 * resolutionSpace / 4)
        {
            uPrevious[i, j] = 10;
        }
        else
        {
            uPrevious[i, j] = 0;
        }

    }
}

// Solver

Array.Copy(uPrevious, u, uPrevious.Length);

for (t = timeStart; t < timeEnd; t += dt)
{
    finish = false;
    iter = 0;

    // Poisson
    while (!finish)
    {
        finish = true;

        Parallel.For(0, resolutionSpace, i =>
        {
            for (int j = 0; j < resolutionSpace; j++)
            {
                delta[i, j] = u[i + 1, j + 1];

                u[i + 1, j + 1] = (1 - omega) * u[i + 1, j + 1]
                                + omega / (1 + 4 * beta) * (uPrevious[i + 1, j + 1] + beta * (u[i + 2, j + 1] + u[i, j + 1] + u[i + 1, j + 2] + u[i + 1, j]
                                                                                 + uPrevious[i + 2, j + 1] + uPrevious[i, j + 1] + uPrevious[i + 1, j + 2] + uPrevious[i + 1, j] - 4 * uPrevious[i + 1, j + 1]));

                delta[i, j] -= u[i + 1, j + 1];
                if (delta[i, j] < 0)
                {
                    delta[i, j] = -delta[i, j];
                }
                if (delta[i, j] > tolerance && finish)
                {
                    finish = false;
                }
            }
        });
        iter++;
    }

    Console.WriteLine("iter = {0}, t = {1}", iter, t);
    Array.Copy(u, uPrevious, uPrevious.Length);
}


for (int i = 0; i < resolutionSpace + 2; i++)
{
    for (int j = 0; j < resolutionSpace + 2; j++)
    {
        Console.Write("{0:f1} ", u[i, j]);
    }
    Console.WriteLine("");
}


