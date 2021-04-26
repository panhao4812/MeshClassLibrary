using Rhino;
using Rhino.Geometry;

using System;
using System.Collections.Generic;

namespace MeshClassLibrary
{
    public class MeshSimplify
    {
        private class MeshSimplify_Point
        {
            public Point3d pos;
            public List<MeshSimplify_Point> refpoints;
            private double tol = 0.01;
            public Vector3d N;
            public double order = 0;
            public List<MeshSimplify_Point> children = new List<MeshSimplify_Point>();
            /*for sort and lay functions
             in sort method  order means clockwise  angel
             in  lay order equals boolean */
            public MeshSimplify_Point(Point3d Pos)
            {
                refpoints = new List<MeshSimplify_Point>();
                pos = Pos;
            }
            public bool isDump(MeshSimplify_Point other)
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
            private static int CompareDinosByLength(MeshSimplify_Point x, MeshSimplify_Point y)
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
                    MeshSimplify_Point p = new MeshSimplify_Point(this.pos);
                    p.refpoints.AddRange(this.refpoints);
                    p.refpoints.Add(this);
                    children.Add(p); children.Add(p);
                }
                for (int i = 0; i < this.refpoints.Count; i++)
                {
                    MeshSimplify_Point p = new MeshSimplify_Point(this.pos);
                    int after = i + 1; if (after > refpoints.Count - 1) after = 0;
                    p.refpoints.Add(refpoints[i]);
                    p.refpoints.Add(refpoints[after]);
                    children.Add(p);
                }

            }
        }
        private class MeshSimplify_Face
        {
            public Rhino.Geometry.MeshFace face = new MeshFace();
            public List<MeshSimplify_Face> reffaces = new List<MeshSimplify_Face>();
            public double order = 0;
            public MeshSimplify_Face() { }
            public MeshSimplify_Face(MeshFace f)
            {
                this.face = f;
            }
            public MeshSimplify_Face(int a, int b, int c, int d)
            {
                this.face = new MeshFace(a, b, c, d);
            }
            public MeshSimplify_Face(int a, int b, int c)
            {
                this.face = new MeshFace(a, b, c);
            }
            public void Add(MeshSimplify_Face f)
            {
                this.reffaces.Add(f);
            }
            public Point3d pos(List<Point3d> vs)
            {
                Point3d position;
                if (this.face.IsQuad)
                {
                    position = vs[this.face.A] + vs[this.face.B] + vs[this.face.C] + vs[this.face.D];
                    position /= 4;
                }
                else if (this.face.IsTriangle)
                {
                    position = vs[this.face.A] + vs[this.face.B] + vs[this.face.C];
                    position /= 3;
                }
                else position = new Point3d();
                return position;
            }
            public Point3d pos(List<Point3f> vs)
            {
                List<Point3d> ps = new List<Point3d>();
                for (int i = 0; i < vs.Count; i++)
                {
                    ps.Add((Point3d)vs[i]);

                }
                return pos(ps);
            }
            public Point3d pos(Rhino.Geometry.Collections.MeshVertexList vs)
            {
                List<Point3d> ps = new List<Point3d>();
                for (int i = 0; i < vs.Count; i++)
                {
                    ps.Add((Point3d)vs[i]);

                }
                return pos(ps);
            }
            public static List<List<MeshSimplify_Face>> Group(List<MeshSimplify_Face> faces)
            {
                for (int i = 0; i < faces.Count; i++)
                {
                    faces[i].order = 0;
                }
                int level = 1;
                for (int i = 0; i < faces.Count; i++)
                {
                    if (FindNext(faces[i], level)) { level++; }
                    ////////////////////////////////////////////////问题在于level
                }
                List<List<MeshSimplify_Face>> output = new List<List<MeshSimplify_Face>>();
                for (int i = 1; i < level; i++)
                {
                    List<MeshSimplify_Face> out1 = new List<MeshSimplify_Face>();
                    output.Add(out1);
                }

                for (int i = 0; i < faces.Count; i++)
                {
                    output[(int)faces[i].order - 1].Add(faces[i]);
                }
                return output;
            }
            public static bool FindNext(MeshSimplify_Face face, int level)
            {
                if (face.order == 0)
                {
                    face.order = level;
                    for (int i = 0; i < face.reffaces.Count; i++)
                    {
                        FindNext(face.reffaces[i], level);
                    }
                    return true;
                }
                return false;
            }
        }
        public MeshSimplify() { }
        private List<MeshSimplify_Point> preDivide(Mesh mesh)
        {
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            List<MeshSimplify_Point> PointList = new List<MeshSimplify_Point>();
            for (int i = 0; i < vs.Count; i++)
            {
                MeshSimplify_Point pt = new MeshSimplify_Point(vs[i]);
                if (vs.MeshVertexIndices(i).Length > 0)
                {
                    pt.N = mesh.Normals[vs.MeshVertexIndices(i)[0]];
                }
                else
                {
                    pt.computeNormal(mesh);
                }
                PointList.Add(pt);
            }
            for (int i = 0; i < vs.Count; i++)
            {
                int[] index = vs.ConnectedTopologyVertices(i);
                for (int j = 0; j < index.Length; j++)
                {
                    PointList[i].refpoints.Add(PointList[index[j]]);
                }
                PointList[i].Sort();
            }

            /////////////////////////////////////////////////////
            for (int i = 0; i < vs.Count; i++)
            {
                PointList[i].order = 0;
            }
            for (int i = 0; i < el.Count; i++)
            {
                if (el.GetConnectedFaces(i).Length == 1)
                {
                    PointList[el.GetTopologyVertices(i).I].order = 5;
                    PointList[el.GetTopologyVertices(i).J].order = 5;
                }
            }
            for (int i = 0; i < vs.Count; i++)
            {
                if (PointList[i].order == 5)
                {
                    if (PointList[i].refpoints.Count != 3)
                    {
                        PointList[i].order = 4;
                    }
                }
                else
                {
                    if (PointList[i].refpoints.Count != 4) PointList[i].order = 4;
                }
            }
            //////////////////////////////////////////////////////////
            for (int k = 0; k < PointList.Count; k++)
            {
                bool sign = true;
                for (int i = 0; i < PointList.Count; i++)
                {
                    if (PointList[i].order == 4)
                    {
                        sign = false;
                        PointList[i].order++;
                        if (PointList[i].refpoints.Count == 4)
                        {
                            if (PointList[i].refpoints[0].order == 5 && PointList[i].refpoints[2].order != 5) { PointList[i].refpoints[2].order = 1; }
                            if (PointList[i].refpoints[1].order == 5 && PointList[i].refpoints[3].order != 5) { PointList[i].refpoints[3].order = 1; }
                            if (PointList[i].refpoints[2].order == 5 && PointList[i].refpoints[0].order != 5) { PointList[i].refpoints[0].order = 1; }
                            if (PointList[i].refpoints[3].order == 5 && PointList[i].refpoints[1].order != 5) { PointList[i].refpoints[1].order = 1; }
                        }
                        else
                        {
                            for (int j = 0; j < PointList[i].refpoints.Count; j++)
                            {
                                PointList[i].refpoints[j].order++;
                            }
                        }
                    }
                }
                for (int i = 0; i < PointList.Count; i++)
                {
                    if (PointList[i].order > 0 && PointList[i].order < 4) PointList[i].order = 4;
                }
                if (sign) { break; }
            }
            return PointList;
        }
        public List<Line> MeshProfile(Mesh mesh)
        {
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            List<MeshSimplify_Point> PointList = preDivide(mesh);
            List<Line> output = new List<Line>();
            for (int i = 0; i < el.Count; i++)
            {
                int a = el.GetTopologyVertices(i).I;
                int b = el.GetTopologyVertices(i).J;
                if (PointList[a].order > 0 && PointList[b].order > 0)
                {
                    output.Add(el.EdgeLine(i));
                }
            }
            /*   List<Point3d> output1 = new List<Point3d>();
               List<string> output2 = new List<string>();
               for(int i = 0;i < PointList.Count;i++){
                 output1.Add(PointList[i].pos);
                 string str = PointList[i].order.ToString();
                 if(PointList[i].order == 0)str = "";
                 output2.Add(str);
               }*/
            return output;
        }
        public List<Mesh> MeshSeperate(Mesh mesh)
        {
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            List<MeshSimplify_Point> PointList = preDivide(mesh);
            List<MeshSimplify_Face> faces = new List<MeshSimplify_Face>();
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                faces.Add(new MeshSimplify_Face(mesh.Faces[i]));
            }
            for (int i = 0; i < el.Count; i++)
            {
                int a = el.GetTopologyVertices(i).I;
                int b = el.GetTopologyVertices(i).J;
                if (PointList[a].order == 0 || PointList[b].order == 0)
                {
                    int[] index = el.GetConnectedFaces(i);
                    if (index.Length == 2)
                    {
                        faces[index[0]].reffaces.Add(faces[index[1]]);
                        faces[index[1]].reffaces.Add(faces[index[0]]);
                    }
                }
            }
            List<Mesh> outputMesh = new List<Mesh>();
            List<List<MeshSimplify_Face>> outtemp = MeshSimplify_Face.Group(faces);
            for (int i = 0; i < outtemp.Count; i++)
            {
                Mesh meshout = new Mesh();
                meshout.Vertices.AddVertices(mesh.Vertices);
                for (int j = 0; j < outtemp[i].Count; j++)
                {
                    meshout.Faces.AddFace(outtemp[i][j].face);
                }
                meshout.Compact();
                outputMesh.Add(meshout);
            }
            return outputMesh;

            /*
    List<Point3d> outpos = new List<Point3d>();
    List<string> outsign = new  List<string>();
    for(int i = 0;i < faces.Count;i++){
      outpos.Add(faces[i].pos(mesh.Vertices));
      string str = i.ToString() + "/" + faces[i].reffaces.Count.ToString();
      outsign.Add(str);
    }
    B = outpos;C = outsign;
    */
        }
        public List<NurbsSurface> MeshSeperate2Nurbs(Mesh mesh)
        {
            MeshConvert conv = new MeshConvert();
            List<NurbsSurface> output = new List<NurbsSurface>();
            List<Mesh> input = MeshSeperate(mesh);

            input.ForEach(delegate (Mesh mesh1)
            {
                NurbsSurface Surf;
                conv.Mesh2Nurbs(mesh1, out Surf);
                output.Add(Surf);
            });
            return output;
        }
    }
    public class MeshFill
    {
        public MeshFill() { }
        public static Mesh MeshFromClosedPoly(List<Polyline> x)
        {
            Mesh mesh = new Mesh();
            for (int i = 0; i < x.Count; i++)
            {
                if (x[i].Count == 4)
                {
                    int n = mesh.Vertices.Count;
                    mesh.Vertices.Add(x[i][0]);
                    mesh.Vertices.Add(x[i][1]);
                    mesh.Vertices.Add(x[i][2]);
                    mesh.Faces.AddFace(new MeshFace(n, n + 1, n + 2));
                }
                else if (x[i].Count == 5)
                {
                    int n = mesh.Vertices.Count;
                    mesh.Vertices.Add(x[i][0]);
                    mesh.Vertices.Add(x[i][1]);
                    mesh.Vertices.Add(x[i][2]);
                    mesh.Vertices.Add(x[i][3]);
                    mesh.Faces.AddFace(new MeshFace(n, n + 1, n + 2, n + 3));
                }
            }
            mesh.Normals.ComputeNormals();
            return mesh;
        }
        public virtual void GetDirections(ref Vertice2 vertice)
        {
            //Vertice2.computeNormal
            return;
        }
        public List<Polyline> Remesh(List<Line> x)
        {
            List<Vertice2> vs; List<IndexPair> id;
            Vertice2.CreateCollection(x, out id, out vs);
            vs = Vertice2.CleanEdge(vs);
            for (int i = 0; i < vs.Count; i++)
            {
                vs[i].Sort(vs);
            }
            for (int i = 0; i < vs.Count; i++)
            {
                Vertice2 v = vs[i];
                GetDirections(ref v);
                vs[i] = v;
            }
            return Vertice2.Remesh(vs);
        }
    }
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
            pl1a = MeshClassLibrary.PolylineSmooth.Catmull_Clark(pl1a, (int)inputdata[3]);
            pl1b = MeshClassLibrary.PolylineSmooth.Catmull_Clark(pl1b, (int)inputdata[3]);
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