﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GeoTools;

namespace Kangaroo
{
    public class CoLinear : GoalObject
    {
        public double Strength;
        public double TargetArea;

        public CoLinear()
        {
        }

        public CoLinear(List<int> V, double k)
        {
            int L = V.Count;
            PIndex = V.ToArray();
            Move = new Vector3d[L];
            Weighting = new double[L];
            for (int i = 0; i < L; i++)
            {
                Weighting[i] = k;
            }
            Strength = k;
        }

        public CoLinear(List<Point3d> V, double k)
        {
            int L = V.Count;
            PPos = V.ToArray();
            Move = new Vector3d[L];
            Weighting = new double[L];
            for (int i = 0; i < L; i++)
            {
                Weighting[i] = k;
            }
            Strength = k;
        }

        public override void Calculate(List<Particle> p)
        {
            int L = PIndex.Length;
            Point3d[] Pts = new Point3d[L];

            for (int i = 0; i < L; i++)
            {
                Pts[i] = p[PIndex[i]].Position;
            }
            Line Ln;
            Line.TryFitLineToPoints(Pts, out Ln);

            for (int i = 0; i < L; i++)
            {
                Move[i] = Ln.ClosestPoint((Pts[i]), false) - Pts[i];
                Weighting[i] = Strength;
            }
        }

    }
    public class CoPlanar : GoalObject
    {
        public double Strength;
        public double TargetArea;

        public CoPlanar()
        {
        }

        public CoPlanar(List<int> V, double k)
        {
            int L = V.Count;
            PIndex = V.ToArray();
            Move = new Vector3d[L];
            Weighting = new double[L];
            for (int i = 0; i < L; i++)
            {
                Weighting[i] = k;
            }
            Strength = k;
        }

        public CoPlanar(List<Point3d> V, double k)
        {
            int L = V.Count;
            PPos = V.ToArray();
            Move = new Vector3d[L];
            Weighting = new double[L];
            for (int i = 0; i < L; i++)
            {
                Weighting[i] = k;
            }
            Strength = k;
        }

        public override void Calculate(List<Particle> p)
        {
            int L = PIndex.Length;
            Point3d[] Pts = new Point3d[L];

            for (int i = 0; i < L; i++)
            {
                Pts[i] = p[PIndex[i]].Position;
            }
            Plane Pl = new Plane();
            Plane.FitPlaneToPoints(Pts, out Pl);

            for (int i = 0; i < L; i++)
            {
                Move[i] = Pl.ClosestPoint(Pts[i]) - Pts[i];
                Weighting[i] = Strength;
            }
        }
     
    }
}
