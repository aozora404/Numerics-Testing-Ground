using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Eddy.NET
{
    public class Vector
    {
        public double X { get; set; } = 0;
        public double Y { get; set; } = 0;
        public double Z { get; set; } = 0;

        public Vector(double x, double y, double z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        

        public static Vector operator +(Vector v1, Vector v2)
        {
            return new Vector(v1.X + v2.X, v1.Y + v2.Y, v1.Z + v2.Z);
        }

        public static Vector operator -(Vector v1, Vector v2)
        {
            return new Vector(v1.X - v2.X, v1.Y - v2.Y, v1.Z - v2.Z);
        }

        public static Vector operator *(double scalar, Vector v)
        {
            return new Vector(scalar * v.X, scalar * v.Y, scalar * v.Z);
        }

        public static Vector operator *(Vector v, double scalar)
        {
            return scalar * v;
        }

        public static Vector operator /(Vector v, double scalar)
        {
            return (1.0/scalar) * v;
        }

        public static Vector operator -(Vector v)
        {
            return -1 * v;
        }
        public double Dot(Vector other)
        {
            return X * other.X + Y * other.Y + Z * other.Z;
        }

        public Vector Cross(Vector other)
        {
            return new Vector(Y * other.Z - Z * other.Y, Z * other.X - X * other.Z, X * other.Y - Y * other.X);
        }

        public double Magnitude()
        {
            return Math.Sqrt(X * X + Y * Y + Z * Z);
        }

        public override string ToString()
        {
            return $"({X:f2}, {Y:f2}, {Z:f2})";
        }
    }
}
