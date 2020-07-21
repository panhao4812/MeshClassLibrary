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
            Polyline pl3 = new Polyline();
            List<double> domain3 = new List<double>();
            double t = 0;
            for (int i = 0; i < pl1.Count; i++)
            {
                Point3d p3 = (pl1[i] + pl2[i]) / 2;
                if (i > 0) t += p3.DistanceTo(pl3[pl3.Count - 1]);
                domain3.Add(t);
                pl3.Add(p3);
            }
            Mesh mesh = new Mesh();
            for (int i = 0; i < pl3.Count; i++)
            {
                mesh.Vertices.Add(pl1[i]);
                mesh.Vertices.Add(pl2[i]);
                double t1 = domain3[i] / pl3.Length;
                // Print(t1.ToString());
                mesh.TextureCoordinates.Add(t1 * texturescale, 0);
                mesh.TextureCoordinates.Add(t1 * texturescale, 1);
                if (i > 0)
                {
                    mesh.Faces.AddFace((i - 1) * 2, i * 2, i * 2 + 1, (i - 1) * 2 + 1);
                }
            }
            mesh.Compact();
            mesh.UnifyNormals();
            return mesh;
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
                    if (i != 0)
                    {
                        meshes.Add(mc.MeshLoft(pl1a, pl2a, false, false));
                    }
                    pl1a.Clear();
                    pl2a.Clear();
                    pl1a.Add(pl[i]);
                    pl2a.Add(pl2[i]);
                }
            }
            return meshes;
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
        public void FitEdge(ref Mesh mesh, List<Polyline> pls, double tol)
        {
            if (pls.Count > 1)
            {
                for (int i = 1; i < pls.Count; i++)
                {
                    pls[0].AddRange(pls[i]);
                }
            }
            FitEdge(ref mesh, pls[0], 1);
            mesh.UnifyNormals();
        }
        public Polyline Project(Polyline pl, Brep sf)
        {
            Mesh[] meshes = Mesh.CreateFromBrep(sf, MeshingParameters.Default);
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
            Curve c = NurbsCurve.CreateControlPointCurve(pl, 3);

            for (int i = 1; i < pl.Count; i++)
            {
                Point3d pt = pl[i];
                double t = 0;
                c.ClosestPoint(pt, out t);
                Vector3d v = c.TangentAt(t);
                v = Vector3d.CrossProduct(v, Vector3d.ZAxis);
                v.Unitize();
                // Print(c.CurvatureAt(t).Length.ToString());
                v *= dis * (1 + c.CurvatureAt(t).Length * 60);
                output.Add(pt + v);
            }
            output.Add(output[0]);
            return output;
        }
        public void CreateMeshRoad(Polyline x, Mesh y, out Mesh road, out List<Mesh> bank,
            out List<Polyline> Ploutput)
        {
            double[] inputdata = { 1, 8, 8, 3 };
            CreateMeshRoad(x, y, inputdata, out road, out bank, out Ploutput);
        }
        public void CreateMeshRoad(Polyline x, Mesh y, double[] inputdata,
            out Mesh road, out List<Mesh> bank, out List<Polyline> Ploutput)
        {
            road = new Mesh();
            bank = new List<Mesh>();
            Ploutput = new List<Polyline>();
            Polyline pl1 = Project(x, y);
            //pl1 = MeshClassLibrary.PolylineSmooth.cc_Subdivide(pl1);
            pl1.Transform(Transform.Translation(0, 0, inputdata[0]));
            Polyline pl1a = OffsetPolyline(pl1, -inputdata[1]);
            Polyline pl1b = OffsetPolyline(pl1, inputdata[2]);
            pl1a = MeshClassLibrary.PolylineSmooth.cc_Subdivide(pl1a, (int)inputdata[3]);
            pl1b = MeshClassLibrary.PolylineSmooth.cc_Subdivide(pl1b, (int)inputdata[3]);
            Polyline pl1af = Project(pl1a, y);
            Polyline pl1bf = Project(pl1b, y);
            FitPoly(ref pl1af, pl1a);
            FitPoly(ref pl1bf, pl1b);
            bank.AddRange(MeshRoadWall(pl1a, pl1af));
            bank.AddRange(MeshRoadWall(pl1bf, pl1b));
            road = MeshUVRoad(pl1a, pl1b, 60);
            Ploutput.Add(pl1af);
            Ploutput.Add(pl1bf);
        }
        void filtCurve(ref Polyline pl, double tol, int step)
        {
            for (int i = 0; i < step; i++)
            {
                filtCurve(ref pl, tol);
            }
        }
        void filtCurve(ref Polyline pl, double tol)
        {
            if (pl[pl.Count - 1].DistanceTo(pl[0]) < 0.001) pl.RemoveAt(pl.Count - 1);
            //Print(pl.Count.ToString());
            double[] t1 = new double[pl.Count];
            for (int i = 0; i < pl.Count; i++)
            {
                int before = i - 1;
                int after = i + 1;
                if (before == -1) before = pl.Count - 1;
                if (after == pl.Count) after = 0;
                Vector3d v1 = pl[i] - pl[before]; v1.Unitize();
                Vector3d v2 = pl[after] - pl[i]; v2.Unitize();
                double t = Vector3d.VectorAngle(v1, v2);
                //Print(t.ToString());
                if (t > tol)
                {
                    double t3 = (t - tol) / Math.PI * 0.5;
                    if (t3 > 0.5) t3 = 0.5;
                    t1[before] += (pl[i].Z - pl[before].Z) * t3;
                    t1[after] += (pl[i].Z - pl[after].Z) * t3;
                }
            }
            for (int i = 0; i < pl.Count; i++)
            {
                pl[i] = pl[i] + new Vector3d(0, 0, t1[i]);
            }
        }
        public List<Mesh> MeshBay(Mesh mesh, int step)
        {

            List<Rface> faces = new List<Rface>();
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                faces.Add(new Rface());
            }

            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            for (int i = 0; i < el.Count; i++)
            {
                int[] index = el.GetConnectedFaces(i);
                if (index.Length == 1)
                {
                    faces[index[0]].age = 1;
                }
                if (index.Length == 2)
                {
                    faces[index[1]]._CollectedFaceIndex.Add(index[0]);
                    faces[index[0]]._CollectedFaceIndex.Add(index[1]);
                }
            }
            for (int k = 1; k < step + 1; k++)
            {
                for (int i = 0; i < faces.Count; i++)
                {
                    if (faces[i].age == k)
                    {
                        for (int j = 0; j < faces[i]._CollectedFaceIndex.Count; j++)
                        {
                            int index2 = faces[i]._CollectedFaceIndex[j];
                            if (faces[index2].age == 0) faces[index2].age = k + 1;
                        }
                    }
                }
            }
           
            Mesh mesh1 = new Mesh(); Mesh mesh2 = new Mesh();
            mesh1.Vertices.AddVertices(mesh.Vertices);
            mesh2.Vertices.AddVertices(mesh.Vertices);
            for (int i = 0; i < faces.Count; i++)
            {
                if (faces[i].age != 0)
                {
                    mesh2.Faces.AddFace(mesh.Faces[i]);

                }
                else
                {
                    mesh1.Faces.AddFace(mesh.Faces[i]);
                }
            }

            mesh1.Compact(); mesh1.UnifyNormals();
            mesh2.Compact(); mesh2.UnifyNormals();
            List<Mesh> output = new List<Mesh>();
            output.Add(mesh1);
            output.Add(mesh2);
            return output;
        }
        public class Rface
        {
            public Rface() { }
            public List<int> _CollectedFaceIndex = new List<int>();
            public int age = 0;
            public override string ToString()
            {
                string str = "age=>" + this.age.ToString() + "\n" + "collect=>";
                for (int i = 0; i < _CollectedFaceIndex.Count; i++)
                {
                    str += _CollectedFaceIndex[i].ToString() + "|";
                }
                return str;
            }
        }
    }
}
