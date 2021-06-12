using System.Collections.Generic;
using GeoTools;
using System;

namespace Kangaroo
{
    public class DynamicWeight1d : GoalObject
    {
        public DynamicWeight1d(Point3d s, Point3d e, double WeightPerLength)
        {
            PPos = new Point3d[2] { s, e };
            Move = new Vector3d[2];
            Weighting = new double[2] { WeightPerLength, WeightPerLength };
        }

        public override void Calculate(List<Particle> p)
        {
            double CurrentLength = p[PIndex[1]].Position.DistanceTo(p[PIndex[0]].Position);
            var Weight = new Vector3d(0, 0, CurrentLength);
            Move[0] = Move[1] = 0.5 * Weight;
        }
    }
    public class LengthMultiple : GoalObject
    {
        public double Stiffness;
        public double Factor;

        public LengthMultiple()
        {
        }

        public LengthMultiple(int S, int E, double F, double k)
        {
            PIndex = new int[2] { S, E };
            Move = new Vector3d[2];
            Weighting = new double[2];
            Stiffness = k;
            Factor = F;
        }

        public LengthMultiple(Point3d S, Point3d E, double F, double k)
        {
            PPos = new Point3d[2] { S, E };
            Move = new Vector3d[2];
            Weighting = new double[2];
            Stiffness = k;
            Factor = F;
        }

        public override void Calculate(List<Particle> p)
        {
            Vector3d current = p[PIndex[1]].Position - p[PIndex[0]].Position;
            double LengthNow = current.Length;
            double RestLength = (Math.Round(LengthNow / Factor)) * Factor;
            double stretchfactor = 1.0 - RestLength / LengthNow;
            Vector3d SpringMove = 0.5 * current * stretchfactor;
            Move[0] = SpringMove;
            Move[1] = -SpringMove;
            Weighting[0] = 2 * Stiffness;
            Weighting[1] = 2 * Stiffness;
        }

    }
    public class Spring : GoalObject
    {
        public double RestLength, Stiffness;

        public Spring()
        {
        }

        public Spring(int s, int e, double l, double k)
        {
            PIndex = new int[2] { s, e };
            Move = new Vector3d[2];
            Weighting = new double[2];
            RestLength = l;
            Stiffness = k;
        }

        public Spring(Point3d s, Point3d e, double l, double k)
        {
            PPos = new Point3d[2] { s, e };
            Move = new Vector3d[2];
            Weighting = new double[2];
            RestLength = l;
            Stiffness = k;
        }
        public Spring(Line L, double length, double k)
        {
            base.PPos = new Point3d[] { L.From, L.To };
            base.Move = new Vector3d[2];
            base.Weighting = new double[2];
            this.RestLength = length;
            this.Stiffness = k;
        }
        public override void Calculate(List<Particle> p)
        {
            Vector3d current = p[PIndex[1]].Position - p[PIndex[0]].Position;
            double stretchfactor = 1.0 - RestLength / current.Length;
            Vector3d SpringMove = 0.5 * current * stretchfactor;
            Move[0] = SpringMove;
            Move[1] = -SpringMove;
            Weighting[0] = 2 * Stiffness;
            Weighting[1] = 2 * Stiffness;
        }

        public override object Output(List<Particle> p)
        {
            return new Line(p[PIndex[0]].Position, p[PIndex[1]].Position);
        }
    }
    /// <summary>
    /// Keep the length of a line between some upper and lower bounds
    /// When it is between these lengths no force is applied
    /// </summary>
    public class ClampLength : GoalObject
    {
        public double Upper;
        public double Lower;
        public double Stiffness;

        public ClampLength()
        {
        }

        public ClampLength(int S, int E, double U, double L, double k)
        {
            PIndex = new int[2] { S, E };
            Move = new Vector3d[2];
            Weighting = new double[2];
            Upper = U;
            Lower = L;
            Stiffness = k;
        }

        public ClampLength(Point3d S, Point3d E, double U, double L, double k)
        {
            PPos = new Point3d[2] { S, E };
            Move = new Vector3d[2];
            Weighting = new double[2]{k,k};
            Upper = U;
            Lower = L;
            Stiffness = k;
        }

        public override void Calculate(List<Particle> p)
        {
            Vector3d current = p[PIndex[1]].Position - p[PIndex[0]].Position;
            double LengthNow = current.Length;
            if (LengthNow > Upper)
            {
                double stretchfactor = 1.0 - Upper / LengthNow;
                Vector3d SpringMove = 0.5 * current * stretchfactor;
                Move[0] = SpringMove;
                Move[1] = -SpringMove;
                Weighting[0] = Stiffness;
                Weighting[1] = Stiffness;
            }
            else if (LengthNow < Lower)
            {
                double stretchfactor = 1.0 - Lower / LengthNow;
                Vector3d SpringMove = 0.5 * current * stretchfactor;
                Move[0] = SpringMove;
                Move[1] = -SpringMove;
                Weighting[0] = Stiffness;
                Weighting[1] = Stiffness;
            }
            else
            {
                Move[0] = Vector3d.Zero;
                Move[1] = Vector3d.Zero;
            }
        }   
    }
    public class ConstantTension : GoalObject
    {
        public double Strength;

        public ConstantTension(Point3d s, Point3d e, double k)
        {
            PPos = new Point3d[2] { s, e };
            Move = new Vector3d[2];
            Weighting = new double[2];
            Strength = 2 * k;
        }

        public override void Calculate(List<Particle> p)
        {
            Vector3d current = p[PIndex[1]].Position - p[PIndex[0]].Position;
            Move[0] = 0.5 * current;
            Move[1] = -0.5 * current;

            Weighting[0] = Weighting[1] = Strength / current.Length;
        }

        public override object Output(List<Particle> p)
        {
            return new Line(p[PIndex[0]].Position, p[PIndex[1]].Position);
        }
    }
    public class Direction : GoalObject
    {
        public Vector3d Dir;
        public double Strength;

        public Direction()
        {
        }

        public Direction(int Start, int End, Vector3d Direction, double K)
        {
            PIndex = new int[2] { Start, End };
            Move = new Vector3d[2];
            Weighting = new double[2] { K, K };
            Dir = Direction;
            Dir.Unitize();
            Strength = K;
        }

        public Direction(Point3d Start, Point3d End, Vector3d Direction, double K)
        {
            PPos = new Point3d[2] { Start, End };
            Move = new Vector3d[2];
            Weighting = new double[2] { K, K };
            Dir = Direction;
            Dir.Unitize();
            Strength = K;
        }

        public override void Calculate(List<Particle> p)
        {
            Point3d S = p[PIndex[0]].Position;
            Point3d E = p[PIndex[1]].Position;
            Vector3d V = E - S;
            Vector3d To = (V - (V * Dir) * Dir) * 0.5;

            Move[0] = To;
            Move[1] = -To;

            Weighting[0] = Strength;
            Weighting[1] = Strength;
        }

    }
    public class PlasticLength : GoalObject
    {
        public double RestLength;
        public double Limit;
        public double Stiffness;

        public PlasticLength(Point3d S, Point3d E, double Lim, double k)
        {
            PPos = new Point3d[2] { S, E };
            Move = new Vector3d[2];
            Weighting = new double[2];
            RestLength = S.DistanceTo(E);
            Limit = Lim;
            Stiffness = k;
        }

        public override void Calculate(List<Particle> p)
        {
            Vector3d current = p[PIndex[1]].Position - p[PIndex[0]].Position;
            double CurrentLength = current.Length;
            double Stretch = CurrentLength - RestLength;

            if (Stretch > Limit)
            {
                RestLength += Stretch - Limit;
            }
            if (-Stretch > Limit)
            {
                RestLength += Stretch + Limit;
            }

            double stretchfactor = 1.0 - RestLength / CurrentLength;
            Vector3d SpringMove = 0.5 * current * stretchfactor;
            Move[0] = SpringMove;
            Move[1] = -SpringMove;

            Weighting[0] = 2 * Stiffness;
            Weighting[1] = 2 * Stiffness;
        }
    }
}
