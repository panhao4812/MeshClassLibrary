using Rhino.Geometry;

using System;
using System.Collections.Generic;

namespace MeshClassLibrary
{
   public class MeshRoad
    {
        MeshCreation mc = new MeshCreation();
        public Mesh MeshUVRoad(Polyline pl1, Polyline pl2, double texturescale)
        {           
            Curve c1 = Curve.CreateControlPointCurve(pl1, 3);
            Curve c2 = Curve.CreateControlPointCurve(pl2, 3);
            Mesh mesh = new Mesh();
            // texturescale *= (c1.GetLength() + c2.GetLength()) / 2;
            for (int i = 0; i < pl1.Count; i++)
            {
                mesh.Vertices.Add(pl1[i]);
                mesh.Vertices.Add(pl2[i]);
                double t1 = -1, t2 = -1;
                c1.ClosestPoint(pl1[i], out t1);
                c2.ClosestPoint(pl2[i], out t2);
                if (i == 0)
                {
                    t1 = 0; t2 = 0;
                }
                else
                {
                    t1 = c1.GetLength(new Interval(c1.Domain.Min, t1)) / c1.GetLength();
                    t2 = t1;
                }        
                mesh.TextureCoordinates.Add(t1 * texturescale, 0);
                mesh.TextureCoordinates.Add(t2 * texturescale, 1);
                if (i > 0)
                {
                    mesh.Faces.AddFace((i - 1) * 2, i * 2, i * 2 + 1, (i - 1) * 2 + 1);

                }
            }
            mesh.UnifyNormals();
            return mesh;
        }
        public void FitPoly(ref Polyline pl, Polyline pl2)
        {
            if (pl.Count < 3 || pl.Count != pl2.Count) return;
            for (int i = 1; i < pl.Count - 1; i++)
            {
                if (pl[i - 1].Z > pl2[i - 1].Z && pl[i + 1].Z < pl2[i + 1].Z) pl[i] = pl2[i];
                if (pl[i - 1].Z < pl2[i - 1].Z && pl[i + 1].Z > pl2[i + 1].Z) pl[i] = pl2[i];
            }
        }
        public List<Mesh> MeshRoadWall(Polyline pl, Polyline pl2)
        {
            List<Mesh> meshes = new List<Mesh>();
            Polyline pl1a = new Polyline();
            Polyline pl2a = new Polyline();
            for (int i = 0; i < pl.Count; i++)
            {
                pl1a.Add(pl[i]);
                pl2a.Add(pl2[i]);
                if (pl[i].DistanceTo(pl2[i]) < 0.001 || i == pl.Count - 1)
                {
                    if (i != 0) meshes.Add(mc.MeshLoft(pl1a, pl2a, false, false));
                    pl1a.Clear();
                    pl2a.Clear();
                    pl1a.Add(pl[i]);
                    pl2a.Add(pl2[i]);
                }
            }
            return meshes;
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
