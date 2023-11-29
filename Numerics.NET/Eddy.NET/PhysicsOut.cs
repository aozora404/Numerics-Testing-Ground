using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Eddy.NET
{
    public class PhysicsOut
    {
        public Vector Position;
        public Vector Velocity;
        public Vector Force;

        public PhysicsOut()
        {
            Position = new Vector(0, 0, 0);
            Velocity = new Vector(0, 0, 0);
            Force = new Vector(0, 0, 0);
        }

        public void set(Vector position, Vector velocity, Vector force)
        {
            Position = position;
            Velocity = velocity;
            Force = force;
        }
    }
}
