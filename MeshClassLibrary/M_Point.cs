using Rhino;
using Rhino.Geometry;
using Rhino.DocObjects;
using Rhino.Collections;

using GH_IO;
using GH_IO.Serialization;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

using System;
using System.IO;
using System.Xml;
using System.Xml.Linq;
using System.Linq;
using System.Data;
using System.Drawing;
using System.Reflection;
using System.Collections;
using System.Windows.Forms;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using Rhino.Geometry.Collections;

namespace MeshClassLibrary
{
    public class M_Point
    {
        public Point3d pos;
        public List<M_Point> refpoints;
        private double tol = 0.01;
        public Vector3d N;
        public double order = 0;
        public List<M_Point> children = new List<M_Point>();
        /*for sort and lay functions
         in sort method  order means clockwise  angel
         in  lay order equals boolean */
        public M_Point(Point3d Pos)
        {
            refpoints = new List<M_Point>();
            pos = Pos;
        }
        public bool isDump(M_Point other)
        {
            return this.pos.DistanceTo(other.pos) < this.tol;
        }
        public void computeNormal()
        {
            N = new Vector3d(0, 0, 1);
        }
        public void computeNormal(Vector3d v)
        {
            N = v;
        }
        public void computeNormal(Surface s)
        {
            double u, v;
            s.ClosestPoint(this.pos, out u, out v);
            N = s.NormalAt(u, v);
        }
        public void computeNormal(Mesh s)
        {
            Point3d outpt; Vector3d outNormal;
            int output = s.ClosestPoint(this.pos, out outpt, out outNormal, double.MaxValue);
            if (output == -1) { N = new Vector3d(0, 0, 1); }
            else
            {
                N = outNormal;
            }
        }
        public void Sort()
        { //sort the refer points in clockwise order
            Plane p1 = new Plane(this.pos, this.N);
            Plane p2 = new Plane(new Point3d(0, 0, 0), new Vector3d(0, 0, 1));
            for (int i = 0; i < this.refpoints.Count; i++)
            {
                Point3d P = new Point3d(this.refpoints[i].pos);
                P.Transform(Transform.PlaneToPlane(p1, p2));
                Vector3d v = new Vector3d(P.X, P.Y, 0);
                double t = 0;
                if (P.Y >= 0) { t = Vector3d.VectorAngle(new Vector3d(1, 0, 0), v); }
                else { t = Math.PI * 2 - Vector3d.VectorAngle(new Vector3d(1, 0, 0), v); }
                this.refpoints[i].order = t;
            }
            this.refpoints.Sort(CompareDinosByLength);
        }
        private static int CompareDinosByLength(M_Point x, M_Point y)
        {
            if (x == null) { if (y == null) { return 0; } else { return -1; } }
            else
            {
                if (y == null) { return 1; }
                else
                {
                    if (x.order > y.order) return 1;
                    if (x.order == y.order) return 0;
                    if (x.order < y.order) return -1;
                    else return 0;
                }
            }
        }
        public void lay()
        {
            Sort();
            if (this.refpoints.Count < 2)
            {
                M_Point p = new M_Point(this.pos);
                p.refpoints.AddRange(this.refpoints);
                p.refpoints.Add(this);
                children.Add(p); children.Add(p);
            }
            for (int i = 0; i < this.refpoints.Count; i++)
            {
                M_Point p = new M_Point(this.pos);
                int after = i + 1; if (after > refpoints.Count - 1) after = 0;
                p.refpoints.Add(refpoints[i]);
                p.refpoints.Add(refpoints[after]);
                children.Add(p);
            }

        }
    }
}
