using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Eddy.NET
{
    public class PhysicsOut
    {
        public double Time { get; set; }

        public double PositionX { get; set; }
        public double PositionY { get; set; }
        public double PositionZ { get; set; }

        public double VelocityX { get; set; }
        public double VelocityY { get; set; }
        public double VelocityZ { get; set; }

        public double ForceX { get; set; }
        public double ForceY { get; set; }
        public double ForceZ { get; set; }

        public PhysicsOut(Vector position, Vector velocity, Vector force, double time)
        {
            PositionX = position.X;
            PositionY = position.Y;
            PositionZ = position.Z;

            VelocityX = velocity.X;
            VelocityY = velocity.Y;
            VelocityZ = velocity.Z;

            ForceX = force.X;
            ForceY = force.Y;
            ForceZ = force.Z;

            Time = time;
        }
    }
}
