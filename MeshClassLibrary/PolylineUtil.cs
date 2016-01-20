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
    class PolyPlannar
    {
        public PolyPlannar() { }
        private MeshCreation mc = new MeshCreation();
        public Mesh Plannar(Polyline pl)
        {
            if (pl[0].DistanceTo(pl[pl.Count - 1]) < Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance)
            {
                pl.RemoveAt(pl.Count - 1);
            }
            Polyline pl2 = new Polyline(pl);
            Mesh mesh = new Mesh();
            if (pl.Count < 3) return mesh;
            if (pl.Count == 3) { mesh.Append(mc.MeshFromPoints(pl2[0], pl2[1], pl2[2])); return mesh; }
            while (pl2.Count >= 3)
            {
                int sign = -1; int before = -1; int after = -1;
                double maxAngle = double.MinValue;
                for (int i = 0; i < pl2.Count; i++)
                {
                    if (i == 0) before = pl2.Count - 1; else before = i - 1;
                    if (i == pl2.Count - 1) after = 0; else after = i + 1;
                    double t = 0;
                 
                    int type = isHullVertice(pl2[before], pl2[i], pl2[after], pl2);
                    t = MinimalAngle(pl2[before], pl2[i], pl2[after]) + Math.Abs(type) * Math.PI;
                 
                    if (t > maxAngle) { maxAngle = t; sign = i; }
                }           
                if (sign == 0) before = pl2.Count - 1; else before = sign - 1;
                if (sign == pl2.Count - 1) after = 0; else after = sign + 1;
                if (sign > -1)
                {
                    mesh.Append(mc.MeshFromPoints(pl2[before], pl2[sign], pl2[after]));
                    pl2.RemoveAt(sign);
                }
            }
            mesh.Weld(Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
            return  mesh;
        }
        public double MinimalAngle(Point3d p1, Point3d p2, Point3d p3)
        {
            double output = double.MaxValue;
            Vector3d v1 = p2 - p1; Vector3d v2 = p3 - p1;
            double t = Vector3d.VectorAngle(v1, v2);
            if (t < output) output = t;
            v1 = p1 - p2; v2 = p3 - p2;
            t = Vector3d.VectorAngle(v1, v2);
            if (t < output) output = t;
            v1 = p2 - p3; v2 = p1 - p3;
            t = Vector3d.VectorAngle(v1, v2);
            if (t < output) output = t;
            return output;
        }
        public int isHullVertice(int i, Polyline pl)
        {
            int before = -1; int after = -1;
            if (i == 0) before = pl.Count - 1; else before = i - 1;
            if (i == pl.Count - 1) after = 0; else after = i + 1;
            return isHullVertice(pl[before], pl[i], pl[after], pl);
        }
        public int isHullVertice(Point3d p1, Point3d p2, Point3d p3, Polyline pl)
        {
            double tol = 0.01;
            Vector3d v1 = p1 - p2, v2 = p2 - p3; Vector3d v3 = Vector3d.CrossProduct(v1, v2);
            Plane plane1 = new Plane(p2, v1, v3);
            Plane plane2 = new Plane(p2, v2, v3);
            int signP1 = 0, signP2 = 0, signP1_1 = 0, signP2_1 = 0;
            for (int i = 0; i < pl.Count; i++)
            {
                if (plane1.DistanceTo(pl[i]) > tol) signP1++;
                if (plane2.DistanceTo(pl[i]) > tol) signP2++;
                if (plane1.DistanceTo(pl[i]) <= tol && plane1.DistanceTo(pl[i]) >= -tol) signP1_1++;
                if (plane2.DistanceTo(pl[i]) <= tol && plane2.DistanceTo(pl[i]) >= -tol) signP2_1++;
            }
            int output = 0;
            //Print(signP1.ToString() + "///" + signP1_1.ToString() + "\n");
            // Print(signP2.ToString() + "///" + signP2_1.ToString());
            // Print(pl.Count.ToString());
            if (signP1 == 0 | (signP1 + signP1_1) == pl.Count) { output+=1; }
            if (signP2 == 0 | (signP2 + signP2_1) == pl.Count) { output-=1; }
            return output;
        }
    }
    class PolylineSmooth
    {
        public static Polyline cc_Subdivide(Polyline ptlist)
        {
            List<Point3d> ps2 = new List<Point3d>();
            if (ptlist.Count < 3)
            {
                return ptlist;
            }
            if (ptlist[0].DistanceTo(ptlist[ptlist.Count - 1]) > 0.001)
            {
                ps2.Add(ptlist[0]);
                Point3d pt = (ptlist[0] + ptlist[1]) / 2; ps2.Add(pt);
                for (int i = 1; i < ptlist.Count - 1; i++)
                {
                    Point3d p1 = new Point3d(ptlist[i - 1]);
                    Point3d p2 = new Point3d(ptlist[i]);
                    Point3d p3 = new Point3d(ptlist[i + 1]);
                    Point3d p4 = (p2 + p3) / 2;
                    ps2.Add(p1 * (1 / 6.0) + p2 * (2 / 3.0) + p3 * (1 / 6.0));
                    ps2.Add(p4);
                }
                ps2.Add(ptlist[ptlist.Count - 1]);
            }
            else
            {
                if (ptlist.Count < 4)
                {
                    return ptlist;
                }
                Point3d p1 = new Point3d(ptlist[ptlist.Count - 2]);
                Point3d p2 = new Point3d(ptlist[0]);
                Point3d p3 = new Point3d(ptlist[0 + 1]);
                Point3d p4 = (p2 + p3) / 2;
                ps2.Add(p1 * (1 / 6.0) + p2 * (2 / 3.0) + p3 * (1 / 6.0));
                ps2.Add(p4);
                for (int i = 1; i < ptlist.Count - 1; i++)
                {
                    p1 = new Point3d(ptlist[i - 1]);
                    p2 = new Point3d(ptlist[i]);
                    p3 = new Point3d(ptlist[i + 1]);
                    p4 = (p2 + p3) / 2;
                    ps2.Add(p1 * (1 / 6.0) + p2 * (2 / 3.0) + p3 * (1 / 6.0));
                    ps2.Add(p4);
                }
                ps2.Add(ps2[0]);
            }

            Polyline pl2 = new Polyline(ps2);
            return pl2;
        }
        public static Polyline cc_Subdivide(Polyline ptlist, int level)
        {
            if (level >= 1)
            {
                for (int i = 0; i < level; i++)
                {
                    ptlist = cc_Subdivide(ptlist);
                }
            }
            return ptlist;
        }
    }
}
