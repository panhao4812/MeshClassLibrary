using Rhino;
using Rhino.Geometry;

using System;
using System.Collections.Generic;

namespace MeshClassLibrary
{
    public class Face2
    {
        public static List<Line> MeshMaze2(Mesh x)
        {
            List<bool> sign;
            List<Face2> fs;
            Random rnd = new Random();
            fs = new List<Face2>();
            sign = new List<bool>();
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = x.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = x.TopologyVertices;
            List<Point3d> FaceC = new List<Point3d>();
            for (int i = 0; i < x.Faces.Count; i++)
            {
                Point3d f = new Point3d();
                if (x.Faces[i].IsQuad)
                {
                    f += x.Vertices[x.Faces[i].A];
                    f += x.Vertices[x.Faces[i].B];
                    f += x.Vertices[x.Faces[i].C];
                    f += x.Vertices[x.Faces[i].D];
                    f /= 4;
                }
                else if (x.Faces[i].IsTriangle)
                {
                    f += x.Vertices[x.Faces[i].A];
                    f += x.Vertices[x.Faces[i].B];
                    f += x.Vertices[x.Faces[i].C];
                    f /= 3;
                }
                FaceC.Add(f);
            }
            for (int i = 0; i < vs.Count; i++)
            {
                fs.Add(new Face2(i));
            }
            for (int i = 0; i < el.Count; i++)
            {
                sign.Add(true);
                IndexPair parel = el.GetTopologyVertices(i);
                fs[parel.I].EdgeIndex.Add(i); fs[parel.J].FaceIndex.Add(parel.I);
                fs[parel.J].EdgeIndex.Add(i); fs[parel.I].FaceIndex.Add(parel.J);
            }
            for (int i = 0; i < vs.Count; i++)
            {
                fs[i].WaveList(rnd);
            }
            int step = 0;
            for (int i = 0; i < fs.Count * 2; i++)
            {
                step = fs[step].FindNext(ref fs, ref sign);
                if (step == -1) break;
            }
            List<Line> output = new List<Line>();
            for (int i = 0; i < el.Count; i++)
            {
                if (sign[i])
                {
                    int[] index = el.GetConnectedFaces(i);
                    if (index.Length == 2)
                    {
                        Point3d p1 = FaceC[index[0]];
                        Point3d p2 = FaceC[index[1]];
                        output.Add(new Line(p1, p2));
                    }
                }
            }
            return output;
        }
        public static List<Line> MeshMaze1(Mesh x)
        {
            List<bool> sign;
            List<Face2> fs;
            Random rnd = new Random();
            fs = new List<Face2>();
            sign = new List<bool>();
            for (int i = 0; i < x.Faces.Count; i++)
            {
                fs.Add(new Face2(i));
            }
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = x.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = x.TopologyVertices;
            for (int i = 0; i < el.Count; i++)
            {
                sign.Add(true);
                int[] parel = el.GetConnectedFaces(i);
                if (parel.Length == 2)
                {
                    fs[parel[0]].EdgeIndex.Add(i); fs[parel[0]].FaceIndex.Add(parel[1]);
                    fs[parel[1]].EdgeIndex.Add(i); fs[parel[1]].FaceIndex.Add(parel[0]);
                }
            }
            for (int i = 0; i < x.Faces.Count; i++)
            {
                fs[i].WaveList(rnd);
            }
            int step = 0;
            for (int i = 0; i < fs.Count * 2; i++)
            {
                step = fs[step].FindNext(ref fs, ref sign);
                if (step == -1) break;
                //Print(step.ToString());
            }
            List<Line> output = new List<Line>();
            for (int i = 0; i < el.Count; i++)
            {
                if (sign[i]) output.Add(el.EdgeLine(i));
            }
            return output;
        }
        public List<int> EdgeIndex = new List<int>();
        public List<int> FaceIndex = new List<int>();
        public int ID = -1;
        public int parent = -1;
        public double energy = 0;
        public Face2(int i)
        {
            this.ID = i;
        }
        public bool WaveList(Random Rnd)
        {
            int dt = 0;
            if (EdgeIndex.Count != FaceIndex.Count) return false;
            if (this.EdgeIndex.Count <= 1) return false;
            if (this.FaceIndex.Count <= 1) return false;
            for (int i = 0; i <= this.EdgeIndex.Count; i++)
            {
                dt = Rnd.Next(FaceIndex.Count);
                int rep = EdgeIndex[0]; int rep2 = EdgeIndex[dt];
                EdgeIndex[0] = rep2;
                EdgeIndex[dt] = rep;
                rep = FaceIndex[0]; rep2 = FaceIndex[dt];
                FaceIndex[0] = rep2;
                FaceIndex[dt] = rep;
            }
            return true;
        }
        public int FindNext(ref List<Face2> fs, ref  List<bool> sign)
        {
            this.energy = 1;
            for (int i = 0; i < this.FaceIndex.Count; i++)
            {
                if (fs[this.FaceIndex[i]].energy == 0)
                {
                    fs[this.FaceIndex[i]].parent = this.ID;
                    sign[this.EdgeIndex[i]] = false;
                    return this.FaceIndex[i];
                }
            }
            return this.parent;
        }
    }
}