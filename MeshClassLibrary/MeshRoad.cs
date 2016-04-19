using Rhino.Geometry;

using System;
using System.Collections.Generic;

namespace MeshClassLibrary
{
   public class MeshRoad
    {
        public void FitPoly(ref Polyline pl, Polyline pl2)
        {
            if (pl.Count < 3 || pl.Count != pl2.Count) return;
            for (int i = 1; i < pl.Count - 1; i++)
            {
                if (pl[i - 1].Z > pl2[i - 1].Z && pl[i + 1].Z < pl2[i + 1].Z) pl[i] = pl2[i];
                if (pl[i - 1].Z < pl2[i - 1].Z && pl[i + 1].Z > pl2[i + 1].Z) pl[i] = pl2[i];
            }
        }
        public void FitEdge(ref Mesh mesh, Polyline pts, double tol)
        {
            MeshCreation mc = new MeshCreation();
            List<int> index = mc.MeshEdgeVerticeIndex(mesh);
            for (int i = 0; i < index.Count; i++)
            {
                Point3d pt = mesh.Vertices[index[i]];
                int sign = 0; double dist = double.MaxValue;
                for (int j = 0; j < pts.Count; j++)
                {
                    double t = pts[j].DistanceTo(pt);
                    if (t < dist) { dist = t; sign = j; }
                }
                mesh.Vertices.SetVertex(index[i], pts[sign]);
            }
            for (int i = 0; i < mesh.Vertices.Count; i++)
            {
                Point3d pt = mesh.Vertices[i];
                int sign = 0; double dist = double.MaxValue;
                for (int j = 0; j < pts.Count; j++)
                {
                    double t = pts[j].DistanceTo(pt);
                    if (t < dist) { dist = t; sign = j; }
                }
                if (dist < tol) mesh.Vertices.SetVertex(i, pts[sign]);
            }

        }
        public Polyline Project(Polyline pl, Brep sf)
        {
            Mesh[] meshes = Mesh.CreateFromBrep(sf);
            if (meshes.Length > 0) return Project(pl, meshes[0]);
            return null;
        }
        public Polyline Project(Polyline pl, Mesh mesh)
        {
            Polyline output = new Polyline();
            for (int i = 0; i < pl.Count; i++)
            {
                Point3d pt = pl[i];
                Line l = new Line(pt + new Vector3d(0, 0, 99999), pt + new Vector3d(0, 0, -99999));
                int[] ids;
                Point3d[] pts = Rhino.Geometry.Intersect.Intersection.MeshLine(mesh, l, out ids);
                if (pts.Length > 0) pt = pts[0];
                output.Add(pt);
            }
            return output;
        }
        public Polyline OffsetPolyline(Polyline pl, double dis)
        {
            Polyline output = new Polyline();
            NurbsCurve c = pl.ToNurbsCurve();
            for (int i = 1; i < pl.Count; i++)
            {
                Point3d pt = pl[i];
                double t = 0;
                c.ClosestPoint(pt, out t);
                Vector3d v = c.TangentAt(t);
                v = Vector3d.CrossProduct(v, Vector3d.ZAxis);
                v.Unitize();
                v *= dis;
                output.Add(pt + v);
            }
            output.Add(output[0]);
            return output;
        }
    }
}
