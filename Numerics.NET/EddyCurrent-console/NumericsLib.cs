namespace NumericsEngine
{
    public static class UniversalConstants
    {
        public const double pi = 3.141592653589793;

        public const double e_0 = 8.8541878128e-12;
        public const double e_0Inverse = 1.1294090674e11;

        public const double m_0 = 1.2566370621e-6;
        public const double m_0Inverse = 7.9577471502e5;

        public const double c = 2.99792458e8;
        public const double cSquared = c * c;

        public const double electronCharge = 1.602176634e-19;
    }

    public static class NumericsLib
    {
        public static Vector2D[,] AddVectorFields(Vector2D[,] fieldA, Vector2D[,] fieldB)
        {
            int length = fieldA.GetLength(0);
            Vector2D[,] result = new Vector2D[length, length];
            for (int y = 0; y < length; y++)
            {
                for (int x = 0; x < length; x++)
                {
                    result[y, x] = fieldA[y, x] + fieldB[y, x];
                }
            }
            return result;
        }
    }

    public class Vector2D
    {
        public double x;
        public double y;

        public Vector2D(double x, double y)
        {
            this.x = x;
            this.y = y;
        }

        public Vector2D Dot(Vector2D other)
        {
            return new Vector2D(this.x * other.x, this.y * other.y);
        }

        public double MagnitudeSquared()
        {
            return x * x + y * y;
        }

        public double Magnitude()
        {
            return Math.Sqrt(x * x + y * y);
        }

        public void SetValue(double x, double y)
        {
            this.x = x;
            this.y = y;
        }

        public static Vector2D operator +(Vector2D a, Vector2D b)
        {
            return new Vector2D(a.x + b.x, a.y + b.y);
        }

        public static Vector2D operator -(Vector2D a, Vector2D b)
        {
            return new Vector2D(a.x - b.x, a.y - b.y);
        }

        public static Vector2D operator -(Vector2D a)
        {
            return new Vector2D(-a.x, -a.y);
        }

        public static Vector2D operator *(double a, Vector2D b)
        {
            return new Vector2D(a * b.x, a * b.y);
        }
    }
}
