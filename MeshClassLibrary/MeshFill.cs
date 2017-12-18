using Rhino;
using Rhino.Geometry;
using System.Collections.Generic;

namespace MeshClassLibrary
{
    public class BasicFace
    {
        public List<int> refer = new List<int>();
        public int ID = -1;
        public bool dead = false;
        public double energy = 0;
        public BasicFace(int i)
        {
            this.ID = i;
        }
        public void Add(int i)
        {
            this.refer.Add(i);
        }
        public List<BasicFace> CreateFromMesh(Mesh x)
        {
            List<BasicFace> fs = new List<BasicFace>();
            for (int i = 0; i < x.Faces.Count; i++)
            {
                fs.Add(new BasicFace(i));
            }
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = x.TopologyEdges;
            for (int i = 0; i < el.Count; i++)
            {
                int[] parel = el.GetConnectedFaces(i);
                if (parel.Length == 2)
                {
                    fs[parel[0]].Add(parel[1]); fs[parel[1]].Add(parel[0]);
                }
            }
            return fs;
        }
    }
    public class Face1 : BasicFace
    {
        public Face1(int i) : base(i) { }
        List<double> Angles = new List<double>();
        public static List<Point3d> DisplayPos(List<Face1> fs, Mesh mesh)
        {
            List<Point3d> output = new List<Point3d>();
            fs.ForEach(delegate (Face1 f)
            {
                MeshFace f1 = mesh.Faces[f.ID];
                if (f1.IsQuad)
                {
                    Point3d p1 = mesh.Vertices[f1.A];
                    Point3d p2 = mesh.Vertices[f1.B];
                    Point3d p3 = mesh.Vertices[f1.C];
                    Point3d p4 = mesh.Vertices[f1.D];
                    Point3d cen = p1 + p2 + p3 + p4;
                    cen /= 4;
                    output.Add(cen);
                }
                else if (f1.IsTriangle)
                {
                    Point3d p1 = mesh.Vertices[f1.A];
                    Point3d p2 = mesh.Vertices[f1.B];
                    Point3d p3 = mesh.Vertices[f1.C];
                    Point3d cen = p1 + p2 + p3;
                    cen /= 3;
                    output.Add(cen);
                }
            });
            return output;
        }
        public static List<string> Displayenergy(List<Face1> fs)
        {
            List<string> output = new List<string>();
            fs.ForEach(delegate (Face1 f) { output.Add(f.energy.ToString()); });
            return output;
        }
        public static List<string> DisplayLife(List<Face1> fs)
        {
            List<string> output = new List<string>();
            fs.ForEach(delegate (Face1 f) { output.Add(f.dead.ToString()); });
            return output;
        }
        public static void CreateFromMesh(Mesh x, out List<Face1> fs)
        {
            x.UnifyNormals();
            x.FaceNormals.ComputeFaceNormals();
            fs = new List<Face1>();
            for (int i = 0; i < x.Faces.Count; i++)
            {
                fs.Add(new Face1(i));
            }
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = x.TopologyEdges;
            for (int i = 0; i < el.Count; i++)
            {
                int[] parel = el.GetConnectedFaces(i);
                if (parel.Length == 2)
                {
                    fs[parel[0]].Add(parel[1]); fs[parel[1]].Add(parel[0]);
                    Vector3d n1 = x.FaceNormals[parel[0]];
                    Vector3d n2 = x.FaceNormals[parel[1]];
                    double Angle1 = Vector3d.VectorAngle(n1, n2);
                    fs[parel[0]].Angles.Add(Angle1);
                    fs[parel[1]].Angles.Add(Angle1);
                }
            }
        }
        public static List<Face1> Fill(Mesh mesh, double angle)
        {
            List<Face1> fs;
            CreateFromMesh(mesh, out fs);
            for (int i = 0; i < fs.Count; i++)
            {
                FillNext(ref fs, i, angle, i + 1);
            }
            return fs;
        }
        public static void FillNext(ref List<Face1> fs, int index, double angle, int flag)
        {
            if (fs[index].energy != 0) return;
            fs[index].energy = flag;
            for (int i = 0; i < fs[index].refer.Count; i++)
            {
                int index2 = fs[index].refer[i];
                if (fs[index].Angles[i] <= angle)
                {
                    FillNext(ref fs, index2, angle, flag);
                }
            }

        }
        public static List<Mesh> MeshFill(Mesh mesh, double angle)
        {
            MeshCreation mc = new MeshCreation();
            List<Face1> fs = Fill(mesh, angle);
            List<Mesh> output = new List<Mesh>();
            for (int i = 0; i < fs.Count; i++)
            {
                Mesh mesh1 = new Mesh();
                for (int j=0;j< fs.Count; j++)
                {               
                    if (fs[j].energy == i+1)
                    {
                        MeshFace f1 = mesh.Faces[fs[j].ID];
                        if (f1.IsQuad)
                        {
                            Point3d p1 = mesh.Vertices[f1.A];
                            Point3d p2 = mesh.Vertices[f1.B];
                            Point3d p3 = mesh.Vertices[f1.C];
                            Point3d p4 = mesh.Vertices[f1.D];
                            mesh1.Append(mc.MeshFromPoints(p1, p2, p3, p4));
                        }
                        else if (f1.IsTriangle)
                        {
                            Point3d p1 = mesh.Vertices[f1.A];
                            Point3d p2 = mesh.Vertices[f1.B];
                            Point3d p3 = mesh.Vertices[f1.C];
                            mesh1.Append(mc.MeshFromPoints(p1, p2, p3));
                        }
                    }                   
                }
               if(mesh1.Vertices.Count>0) output.Add(mesh1);
            }
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
}
