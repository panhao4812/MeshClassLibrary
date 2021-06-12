using GeoTools;
using System;
using System.Collections.Generic;


namespace Kangaroo
{
    public class TargetPlane : GoalObject
    {
        public double Strength;
        public Plane _plane;

        public TargetPlane()
        {
        }

        public TargetPlane(int P, double k, Plane Target)
        {
            PIndex = new int[1] { P };
            Move = new Vector3d[1];
            Weighting = new double[1];
            Strength = k;
            _plane = Target;
        }

        public override void Calculate(List<Particle> p)
        {
            Point3d ThisPt = p[PIndex[0]].Position;
            Move[0] = _plane.ClosestPoint(ThisPt) - ThisPt;
            Weighting[0] = Strength;
        }
    }
    public class MagnetSnap : GoalObject
    {
        public double Strength;
        public double Range;
        public double RangeSq;

        public MagnetSnap(List<Point3d> V, double R, double k)
        {
            int L = V.Count;
            PPos = V.ToArray();
            Move = new Vector3d[L];
            Weighting = new double[L];
            for (int i = 0; i < L; i++)
            {
                Weighting[i] = k;
            }
            Range = R;
            RangeSq = R * R;
            Strength = k;
        }

        public override void Calculate(List<Particle> p)
        {
            int L = PIndex.Length;
            double[] Xcoord = new double[L];
            for (int i = 0; i < L; i++)
            {
                Xcoord[i] = p[PIndex[i]].Position.X;
            }
            Array.Sort(Xcoord, PIndex);

            for (int i = 0; i < L; i++)
            {
                Move[i] = Vector3d.Zero;
                Weighting[i] = 0;
            }

            for (int i = 0; i < (PIndex.Length - 1); i++)
            {
                for (int j = 1; (i + j) < PIndex.Length; j++)
                {
                    int k = i + j;
                    Vector3d Separation = p[PIndex[k]].Position - p[PIndex[i]].Position;
                    if (Separation.X < Range)
                    {
                        if (Separation.SquareLength < RangeSq)
                        {
                            Move[i] += 0.5 * Separation;
                            Move[k] -= 0.5 * Separation;
                            Weighting[i] = Strength;
                            Weighting[k] = Strength;
                        }
                    }
                    else { break; }
                }
            }
        }

    }
    public class Unary : GoalObject
    {
        public Vector3d Force;        

        public Unary()
        {
        }

        public Unary(int u, Vector3d v)
        {
            PIndex = new int[1] { u };
            Move = new Vector3d[1];
            Weighting = new double[1];
            Force = v;
        }

        public Unary(Point3d P, Vector3d v)
        {
            PPos = new Point3d[1] { P };
            Move = new Vector3d[1];
            Weighting = new double[1];
            Force = v;
        }

        public override void Calculate(List<Particle> p)
        {            
            Move[0] = Force;                  
            Weighting[0] = 1.0;
        }
      
    }
    public class Anchor : GoalObject
    {
        public double Strength;
        public Point3d Pt;

        public Anchor()
        {
        }

        /// <summary>
        /// Construct a new Anchor object by particle index and target position
        /// </summary>
        /// <param name="Id">The integer index of the particle to anchor.</param>
        /// <param name="P">The target position to keep the particle at.</param>
        /// <param name="K">Strength of the Anchor. For an absolute anchor, you can use double.MaxValue here.</param>
        public Anchor(int Id, Point3d P, double k)
        {
            PIndex = new int[1] { Id };
            Move = new Vector3d[1];
            Weighting = new double[1] { k };
            Strength = k;
            Pt = P;
        }

        /// <summary>
        /// Construct a new Anchor object by position.
        /// </summary>        
        /// <param name="P">Particle starting position. Also used as the target position to keep the particle at.</param>
        /// <param name="K">Strength of the Anchor. For an absolute anchor, you can use double.MaxValue here.</param>
        public Anchor(Point3d P, double k)
        {
            PPos = new Point3d[1] { P };
            Move = new Vector3d[1];
            Weighting = new double[1] { k };
            Strength = k;
            Pt = P;
        }

        public override void Calculate(List<Particle> p)
        {
            Move[0] = Pt - p[PIndex[0]].Position;
            Weighting[0] = Strength;
        }
    }
    public class AnchorXYZ : GoalObject
    {
        public double Strength;
        public Point3d PtA;
        public bool xFix, yFix, zFix;

        public AnchorXYZ(Point3d P, bool X, bool Y, bool Z, double k)
        {
            PPos = new Point3d[1] { P };
            Move = new Vector3d[1];
            Weighting = new double[1] { k };
            Strength = k;
            PtA = P;
            xFix = X;
            yFix = Y;
            zFix = Z;
        }

        public override void Calculate(List<Particle> p)
        {
            Vector3d V = PtA - p[PIndex[0]].Position;
            if (!xFix) { V.X = 0; }
            if (!yFix) { V.Y = 0; }
            if (!zFix) { V.Z = 0; }
            Move[0] = V;
            Weighting[0] = Strength;
        }
    }
    public class Coincident : GoalObject
    {
        public double Strength;

        public Coincident(Point3d P0, Point3d P1, double Strength)
        {
            this.Strength = Strength;
            PPos = new Point3d[2] { P0, P1 };
            Move = new Vector3d[2];
            Weighting = new Double[2] { Strength, Strength };
        }

        public override void Calculate(List<Particle> p)
        {
            Move[0] = 0.5 * (p[PIndex[1]].Position - p[PIndex[0]].Position);
            Move[1] = -1 * Move[0];
        }
    }
    public class FloorPlane : GoalObject
    {
        public double Strength;

        public FloorPlane()
        {
        }

        public FloorPlane(double k)
        {
            PIndex = new int[1] { 0 };
            Move = new Vector3d[1];
            Weighting = new double[1];
            Strength = k;
        }

        public override void Calculate(List<Particle> p)
        {
            int L = p.Count;

            PIndex = new int[L];
            Move = new Vector3d[L];
            Weighting = new double[L];

            for (int i = 0; i < L; i++)
            {
                PIndex[i] = i;
                double Height = p[i].Position.Z;
                if (Height < 0)
                {
                    Move[i] = -Vector3d.ZAxis * Height;
                    Weighting[i] = Strength;
                }
                else
                {
                    Move[i] = Vector3d.Zero;
                    Weighting[i] = 0;
                }

            }
        }

    }
    public class PlasticAnchor : GoalObject
    {
        public Point3d Location;
        public double Limit;

        public PlasticAnchor(Point3d P, double R, double k)
        {
            PPos = new Point3d[1] { P };
            Move = new Vector3d[1];
            Weighting = new double[1] { k };
            Location = P;
            Limit = R;
        }

        public override void Calculate(List<Particle> p)
        {
            Point3d ThisPt = p[PIndex[0]].Position;
            Vector3d Between = Location - ThisPt;
            Move[0] = Between;
            double stretch = Between.Length - Limit;
            if (stretch > 0)
            {
                Between.Unitize();
                Between *= stretch;
                Location -= Between;
                Move[0] -= Between;
            }
        }
    }
    public class Transformer : GoalObject
    {
        public double Strength;
        public Transform XForm;
        public Transform Inverse;

        public Transformer()
        {
        }

        public Transformer(int P0, int P1, Transform T, double k)
        {
            PIndex = new int[2] { P0, P1 };
            Move = new Vector3d[2];
            Weighting = new double[2] { k, k };
            Strength = k;
            XForm = T;
            Inverse = T;
            T.Invert();
        }

        public Transformer(Point3d P0, Point3d P1, Transform T, double k)
        {
            PPos = new Point3d[2] { P0, P1 };
            Move = new Vector3d[2];
            Weighting = new double[2] { k, k };
            Strength = k;
            XForm = T;
            Inverse = T;
            T.Invert();
        }

        public override void Calculate(List<Particle> p)
        {

            Point3d PT0 = p[PIndex[0]].Position;
            PT0.Transform(XForm);
            Vector3d Match = PT0 - p[PIndex[1]].Position;
            Move[1] = 0.5 * Match;

            Point3d PT1 = p[PIndex[1]].Position;
            PT1.Transform(Inverse);
            Match = PT1 - p[PIndex[0]].Position;

            Move[0] = 0.5 * Match;
            Weighting[0] = Weighting[1] = Strength;
        }


    }
}
