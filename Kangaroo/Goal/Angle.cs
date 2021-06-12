using System;
using System.Collections.Generic;
using GeoTools;

namespace Kangaroo
{
    public class Angle : GoalObject
    {
        public double EI;
        public double RestAngle;

        public Angle()
        {
        }

        /// <summary>
        /// Construct a new Angle goal by particle index.
        /// </summary>
        /// <param name="Strength">Strength of this goal.</param>
        /// <param name="RA">Rest Angle.</param>
        /// <param name="P0">Start of the first line segment.</param>
        /// <param name="P1">End of the first line segment. This can be identical to P2 if the line segments are connected.</param>
        /// <param name="P2">Start of the second line segment. This can be identical to P1 if the line segments are connected.</param>
        /// <param name="P3">End of the second line segment.</param>
        public Angle(double Strength, double RA, int P0, int P1, int P2, int P3)
        {
          PIndex = new int[4]{P0,P1,P2,P3};
          Move = new Vector3d[4];
          Weighting = new double[4];
          EI = Strength;
          RestAngle = RA;
        }

        public Angle(Line L0, Line L1, double RA, double Strength)
        {
            PPos = new Point3d[4] { L0.From, L0.To, L1.From, L1.To };
            Move = new Vector3d[4];
            Weighting = new double[4];
            EI = Strength;
            RestAngle = RA;
        }

        public override void Calculate(List<Particle> p)
        {
          Point3d P0 = p[PIndex[0]].Position;
          Point3d P1 = p[PIndex[1]].Position;
          Point3d P2 = p[PIndex[2]].Position;
          Point3d P3 = p[PIndex[3]].Position;

          Vector3d V01 = P1 - P0;
          Vector3d V23 = P3 - P2;
          double top = 2* Math.Sin(Vector3d.VectorAngle(V01, V23) - RestAngle);
          double Lc = (V01 + V23).Length;
          double Sa = top / (V01.Length * Lc);
          double Sb = top / (V23.Length * Lc);
          
          Vector3d Perp = Vector3d.CrossProduct(V01, V23);
          Vector3d ShearA = Vector3d.CrossProduct(V01, Perp);
          Vector3d ShearB = Vector3d.CrossProduct(Perp, V23);

          ShearA.Unitize();
          ShearB.Unitize();

          ShearA *= Sa;
          ShearB *= Sb;

          Move[0] = ShearA;
          Move[1] = -ShearA;
          Move[2] = ShearB;
          Move[3] = -ShearB;

          Weighting[0] = EI;
          Weighting[1] = EI;
          Weighting[2] = EI;
          Weighting[3] = EI;
        }       
    }
    public class Angle2 : GoalObject
    {
        public double EI;
        public double RestAngle;

        public Angle2()
        {
        }

        /// <summary>
        /// Construct a new Angle goal by particle index.
        /// </summary>
        /// <param name="Strength">Strength of this goal.</param>
        /// <param name="RA">Rest Angle.</param>
        /// <param name="P0">Start of the first line segment.</param>
        /// <param name="P1">End of the first line segment. This can be identical to P2 if the line segments are connected.</param>
        /// <param name="P2">Start of the second line segment. This can be identical to P1 if the line segments are connected.</param>
        /// <param name="P3">End of the second line segment.</param>
        public Angle2(double Strength, double RA, int P0, int P1, int P2, int P3)
        {
            PIndex = new int[4] { P0, P1, P2, P3 };
            Move = new Vector3d[4];
            Weighting = new double[4];
            EI = Strength;
            RestAngle = RA;
        }

        public Angle2(Line L0, Line L1, double RA, double Strength)
        {
            PPos = new Point3d[4] { L0.From, L0.To, L1.From, L1.To };
            Move = new Vector3d[4];
            Weighting = new double[4];
            EI = Strength;
            RestAngle = RA;
        }

        public override void Calculate(List<Particle> p)
        {

            Point3d P0 = p[PIndex[0]].Position;
            Point3d P1 = p[PIndex[1]].Position;
            Point3d P2 = p[PIndex[2]].Position;
            Point3d P3 = p[PIndex[3]].Position;

            Vector3d V01 = P1 - P0;
            Vector3d V23 = P3 - P2;
            double top = 2 * Math.Sin(Vector3d.VectorAngle(V01, V23) - RestAngle);
            double Lc = (V01 + V23).Length;

            double Sa = V01.Length * top / Lc;
            double Sb = V23.Length * top / Lc;

            Vector3d Perp = Vector3d.CrossProduct(V01, V23);
            Perp.Unitize();
            Vector3d ShearA = Vector3d.CrossProduct(V01, Perp);
            Vector3d ShearB = Vector3d.CrossProduct(Perp, V23);

            ShearA *= 0.25 * Sa;
            ShearB *= 0.25 * Sb;

            Move[0] = ShearA;
            Move[1] = -ShearA;
            Move[2] = ShearB;
            Move[3] = -ShearB;

            Weighting[0] = Weighting[1] = 4 * EI / (V01.Length * V01.Length * V01.Length);
            Weighting[2] = Weighting[3] = 4 * EI / (V23.Length * V23.Length * V23.Length);

        }
    }
    public class AngleMultiple : GoalObject
    {
        public double EI;
        public double AngleFactor;

        public AngleMultiple()
        {
        }

        public AngleMultiple(double Strength, double AF, int P0, int P1, int P2, int P3)
        {
            PIndex = new int[4] { P0, P1, P2, P3 };
            Move = new Vector3d[4];
            Weighting = new double[4];
            EI = Strength;
            AngleFactor = AF;
        }

        public AngleMultiple(Line L0, Line L1, double AF, double Strength)
        {
            PPos = new Point3d[4] { L0.From, L0.To, L1.From, L1.To };
            Move = new Vector3d[4];
            Weighting = new double[4];
            EI = Strength;
            AngleFactor = AF;
        }

        public override void Calculate(List<Particle> p)
        {
            Point3d P0 = p[PIndex[0]].Position;
            Point3d P1 = p[PIndex[1]].Position;
            Point3d P2 = p[PIndex[2]].Position;
            Point3d P3 = p[PIndex[3]].Position;

            Vector3d V01 = P1 - P0;
            Vector3d V23 = P3 - P2;

            double AngleNow = Vector3d.VectorAngle(V01, V23);
            double RestAngle = (Math.Round(AngleNow / AngleFactor)) * AngleFactor;

            double top = Math.Sin(Vector3d.VectorAngle(V01, V23) - RestAngle);
            double Lc = (V01 + V23).Length;
            double Sa = top / (V01.Length * Lc);
            double Sb = top / (V23.Length * Lc);
            Vector3d Perp = Vector3d.CrossProduct(V01, V23);
            Vector3d ShearA = Vector3d.CrossProduct(V01, Perp);
            Vector3d ShearB = Vector3d.CrossProduct(Perp, V23);
            ShearA.Unitize();
            ShearB.Unitize();
            ShearA *= Sa * 0.5;
            ShearB *= Sb * 0.5;

            Move[0] = ShearA;
            Move[1] = -ShearA;
            Move[2] = ShearB;
            Move[3] = -ShearB;

            Weighting[0] = EI;
            Weighting[1] = EI;
            Weighting[2] = EI;
            Weighting[3] = EI;
        }
    }
    public class ClampAngle : GoalObject
    {
        public double EI;
        public double Upper;
        public double Lower;

        public ClampAngle(Line L0, Line L1, double Upp, double Low, double Strength)
        {
            PPos = new Point3d[4] { L0.From, L0.To, L1.From, L1.To };
            Move = new Vector3d[4];
            Weighting = new double[4];
            EI = Strength;
            Upper = Upp;
            Lower = Low;
        }

        public override void Calculate(List<Particle> p)
        {
            Point3d P0 = p[PIndex[0]].Position;
            Point3d P1 = p[PIndex[1]].Position;
            Point3d P2 = p[PIndex[2]].Position;
            Point3d P3 = p[PIndex[3]].Position;

            Vector3d V01 = P1 - P0;
            Vector3d V23 = P3 - P2;

            double CurrentAngle = Vector3d.VectorAngle(V01, V23);
            double RestAngle = 0;
            bool Active = false;
            if (CurrentAngle > Upper)
            {
                RestAngle = Upper;
                Active = true;
            }
            else if (CurrentAngle < Lower)
            {
                RestAngle = Lower;
                Active = true;
            }

            if (Active)
            {
                double top = Math.Sin(Vector3d.VectorAngle(V01, V23) - RestAngle);
                double Lc = (V01 + V23).Length;
                double Sa = top / (V01.Length * Lc);
                double Sb = top / (V23.Length * Lc);

                Vector3d Perp = Vector3d.CrossProduct(V01, V23);
                Vector3d ShearA = Vector3d.CrossProduct(V01, Perp);
                Vector3d ShearB = Vector3d.CrossProduct(Perp, V23);

                ShearA *= Sa * 0.5;
                ShearB *= Sb * 0.5;

                Move[0] = ShearA;
                Move[1] = -ShearA;
                Move[2] = ShearB;
                Move[3] = -ShearB;

                Weighting[0] = EI;
                Weighting[1] = EI;
                Weighting[2] = EI;
                Weighting[3] = EI;
            }
            else
            {
                Move[0] = Move[1] = Move[2] = Move[3] = Vector3d.Zero;
                Weighting[0] = Weighting[1] = Weighting[2] = Weighting[3] = 0;
            }
        }
    }
}
