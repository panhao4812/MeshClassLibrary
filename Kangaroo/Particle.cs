using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Kangaroo
{
    public class Particle
    {
        // Fields
        public Vector3d AngularVelocity;
        public double Mass;
        public Vector3d MoveSum;
        public Plane Orientation;
        public Point3d Position;
        public Plane StartOrientation;
        public Point3d StartPosition;
        public Vector3d TorqueSum;
        public double TorqueWeightSum;
        public Vector3d Velocity;
        public double WeightSum;

        // Methods
        public Particle(Plane Pl, double mass)
        {
            this.Position = Pl.Origin;
            this.StartPosition = Pl.Origin;
            this.MoveSum = Vector3d.Zero;
            this.WeightSum = 0.0;
            this.Velocity = Vector3d.Zero;
            this.Mass = mass;
            this.Orientation = Pl;
            this.StartOrientation = Pl;
            this.TorqueSum = Vector3d.Zero;
            this.AngularVelocity = Vector3d.Zero;
            this.TorqueWeightSum = 0.0;
        }

        public Particle(Point3d p, double m)
        {
            this.Position = p;
            this.StartPosition = p;
            this.MoveSum = Vector3d.Zero;
            this.WeightSum = 0.0;
            this.Velocity = Vector3d.Zero;
            this.Mass = m;
            this.TorqueSum = Vector3d.Zero;
            this.TorqueWeightSum = 0.0;
        }

        public void ClearForces()
        {
            this.MoveSum = Vector3d.Zero;
            this.WeightSum = 0.0;
            this.TorqueSum = Vector3d.Zero;
            this.TorqueWeightSum = 0.0;
        }
    }
  


}
