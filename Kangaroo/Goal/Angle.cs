﻿using System;
using System.Collections.Generic;


namespace Kangaroo
{
    /// <summary>
    /// Make 2 directed lines form a given angle
    /// (This is the same as bending, except the 2 lines need not share a vertex)
    /// 
    /// The approach used here is described in the paper "Tensegrity spline beam and grid shell structures"
    /// by S.M.L. Adriaenssens, M.R.Barnes
    /// http://formfindinglab.princeton.edu/wp-content/uploads/2011/09/Spline-Beam.pdf
    /// </summary>
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
}