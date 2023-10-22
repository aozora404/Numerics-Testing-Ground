using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Eddy.NET
{
    public class Vector2D
    {
        public double X { get; set; } = 0;
        public double Y { get; set; } = 0;

        public Vector2D(double x, double y)
        {
            X = x;
            Y = y;
        }

        public static Vector2D operator +(Vector2D v1, Vector2D v2)
        {
            return new Vector2D(v1.X + v2.X, v1.Y + v2.Y);
        }

        public static Vector2D operator -(Vector2D v1, Vector2D v2)
        {
            return new Vector2D(v1.X - v2.X, v1.Y - v2.Y);
        }

        public static Vector2D operator *(double scalar, Vector2D v)
        {
            return new Vector2D(scalar * v.X, scalar * v.Y);
        }

        public static Vector2D operator *(Vector2D v, double scalar)
        {
            return scalar * v;
        }

        public double Dot(Vector2D other)
        {
            return X * other.X + Y * other.Y;
        }

        public double Magnitude()
        {
            return Math.Sqrt(X * X + Y * Y);
        }

        public override string ToString()
        {
            return $"({X:f2}, {Y:f2})";
        }
    }
}
