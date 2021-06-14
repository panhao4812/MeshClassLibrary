
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GeoTools;

namespace Kangaroo
{
    public interface IGoal
    {
        // Methods
        void Calculate(List<Particle> P);
        IGoal Clone();
        object Output(List<Particle> P);

        // Properties
        Plane[] InitialOrientation { get; set; }
        Vector3d[] Move { get; set; }
        string Name { get; set; }
        int[] PIndex { get; set; }
        Point3d[] PPos { get; set; }
        Vector3d[] Torque { get; set; }
        double[] TorqueWeighting { get; set; }
        double[] Weighting { get; set; }
    }
    public abstract class GoalObject : IGoal
    {
        // Methods
        protected GoalObject()
        {
        }

        public abstract void Calculate(List<Particle> p);
        public IGoal Clone() =>
            (base.MemberwiseClone() as IGoal);

        public Point3d[] GetCurrentPositions(List<Particle> p)
        {
            Point3d[] pointdArray = new Point3d[this.PIndex.Length];
            for (int i = 0; i < this.PIndex.Length; i++)
            {
                pointdArray[i] = p[this.PIndex[i]].Position;
            }
            return pointdArray;
        }

        public virtual object Output(List<Particle> p) =>
            null;

        // Properties
        public Plane[] InitialOrientation { get; set; }

        public Vector3d[] Move { get; set; }

        public string Name { get; set; }

        public int[] PIndex { get; set; }

        public Point3d[] PPos { get; set; }

        public Vector3d[] Torque { get; set; }

        public double[] TorqueWeighting { get; set; }

        public double[] Weighting { get; set; }
    }
}
