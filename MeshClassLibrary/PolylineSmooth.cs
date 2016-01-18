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
        public Mesh Plannar(Polyline pl)
        {
            if (pl[0].DistanceTo(pl[pl.Count - 1]) < Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance)
            {
                pl.RemoveAt(pl.Count - 1);
            }
            Mesh mesh = new Mesh();
            if (pl.Count < 3) return mesh;
            List<int> index = new List<int>();
            for (int i = 0; i < pl.Count; i++)
            {
                mesh.Vertices.Add(pl[i]);
                index.Add(i);
            }
            if (pl.Count == 3) { mesh.Faces.AddFace(0, 1, 2); return mesh; }
            while (index.Count >= 3)
            {
                int sign = -1; int before = -1; int after = -1;
                double maxAngle = double.MinValue;
                for (int i = 0; i < index.Count; i++)
                {
                    if (i == 0) before = index.Count - 1; else before = i - 1;
                    if (i == index.Count - 1) after = 0; else after = i + 1;
                    double t = MinimalAngle(pl[index[before]], pl[index[i]], pl[index[after]]);
                    if (t > maxAngle) { maxAngle = t; sign = i; }
                }
                if (sign == 0) before = index.Count - 1; else before = sign - 1;
                if (sign == index.Count - 1) after = 0; else after = sign + 1;
                if (sign > -1)
                {
                    mesh.Faces.AddFace(index[before], index[sign], index[after]);
                    index.RemoveAt(sign);
                }
            }
            return mesh;

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
