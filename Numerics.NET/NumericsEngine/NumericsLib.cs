using System;

namespace NumericsEngine
{
    public static class UniversalConstants
    {
        public const float pi = 3.141592653589793f;

        public const float e_0 = 8.8541878128f * 0.000000000001f; //10^-12
        public const float e_0Inverse = 1.1294090674f * 100000000000; // 10^11

        public const float m_0 = 1.2566370621f * 0.000001f; //10^-6
        public const float m_0Inverse = 7.9577471502f * 100000; //10^5

        public const float c = 299792458;
        public const float cSquared = c * c;

        public const float electronCharge = 1.602176634f * 0.0000000000000000001f; //10^-19
    }

    public static class NumericsLib
    {
        public static double[] GradientFlux(double[,,] ScalarFunction, int x, int y, int z)
        {
            return new double[3] { 0.5*(ScalarFunction[z, y, x + 1] - ScalarFunction[z, y, x - 1]),
                                   0.5*(ScalarFunction[z, y + 1, x] - ScalarFunction[z, y - 1, x]),
                                   0.5*(ScalarFunction[z + 1, y, x] - ScalarFunction[z - 1, y, x])};
        }

        public static double DivergenceFlux(double[,,,] VectorFunction, int x, int y, int z)
        {
            return 0.5 * (VectorFunction[z, y, x + 1, 0] - VectorFunction[z, y, x - 1, 0] +
                          VectorFunction[z, y + 1, x, 1] - VectorFunction[z, y - 1, x, 1] +
                          VectorFunction[z + 1, y, x, 2] - VectorFunction[z - 1, y, x, 2]);
        }

        public static double[] CurlFlux(double[,,,] VectorFunction, int x, int y, int z)
        {
            return new double[3] { 0.5*(VectorFunction[z, y + 1, x, 2] - VectorFunction[z, y - 1, x, 2]
                                      - VectorFunction[z + 1, y, x, 1] + VectorFunction[z - 1, y, x, 1]),
                                   0.5*(VectorFunction[z + 1, y, x, 0] - VectorFunction[z - 1, y, x, 0]
                                      - VectorFunction[z, y, x + 1, 2] + VectorFunction[z, y, x - 1, 2]),
                                   0.5*(VectorFunction[z, y, x + 1, 1] - VectorFunction[z, y, x - 1, 1]
                                      - VectorFunction[z, y + 1, x, 0] + VectorFunction[z, y - 1, x, 0])};
        }

        public static double[] GradientFlux2D(double[,] ScalarFunction, int x, int y)
        {
            return new double[2] { 0.5*(ScalarFunction[y, x + 1] - ScalarFunction[y, x - 1]),
                                   0.5*(ScalarFunction[y + 1, x] - ScalarFunction[y - 1, x])};
        }

        public static double DivergenceFlux2D(double[,,] VectorFunction, int x, int y)
        {
            return 0.5 * (VectorFunction[y, x + 1, 0] - VectorFunction[y, x - 1, 0] +
                          VectorFunction[y + 1, x, 1] - VectorFunction[y - 1, x, 1]);
        }

        public static double CurlFlux2D(double[,,] VectorFunction, int x, int y)
        {
            return  0.5*(VectorFunction[y, x + 1, 1] - VectorFunction[y, x - 1, 1]
                       - VectorFunction[y + 1, x, 0] + VectorFunction[y - 1, x, 0]);
        }


    }
}
