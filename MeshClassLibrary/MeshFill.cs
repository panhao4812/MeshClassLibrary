using Rhino;
using Rhino.Geometry;
using System.Collections.Generic;

namespace MeshClassLibrary
{
   
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
