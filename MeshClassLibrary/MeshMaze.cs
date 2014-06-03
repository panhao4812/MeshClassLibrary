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
    public class BasicFace
    {
        public static List<Line> MeshMaze2(Mesh x)
        {
            List<bool> sign;
            List<BasicFace> fs;
            Random rnd = new Random();
            fs = new List<BasicFace>();
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
                fs.Add(new BasicFace(i));
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
            List<BasicFace> fs;
            Random rnd = new Random();
            fs = new List<BasicFace>();
            sign = new List<bool>();
            for (int i = 0; i < x.Faces.Count; i++)
            {
                fs.Add(new BasicFace(i));
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
        public BasicFace(int i)
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
        public int FindNext(ref List<BasicFace> fs, ref  List<bool> sign)
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
    public class MeshFollowLines
    {
        public MeshFollowLines() { }
        public List<Line> followlines1(Mesh mesh, List<double> t, double iso)
        {
            List<Line> ls = new List<Line>();
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                if (mesh.Faces[i].IsTriangle)
                {
                    Point3d p1 = mesh.Vertices[mesh.Faces[i].A];
                    Point3d p2 = mesh.Vertices[mesh.Faces[i].B];
                    Point3d p3 = mesh.Vertices[mesh.Faces[i].C];
                    double t1 = t[mesh.Faces[i].A];
                    double t2 = t[mesh.Faces[i].B];
                    double t3 = t[mesh.Faces[i].C];
                    Solve3Face(p1, p2, p3, t1, t2, t3, iso, ref ls);
                }
            }
            return ls;
        }
        public bool Solve3Face(Point3d p1, Point3d p2, Point3d p3, double t1, double t2, double t3, double iso, ref List<Line> lines)
        {
            if (t1 == iso && t2 == iso && t3 == iso) return false;

            if (t1 == iso)
            {
                if ((iso >= t2 && iso <= t3) || (iso >= t3 && iso <= t2))
                {
                    Point3d p;
                    if (t2 == t3) { p = (p2 + p3) / 2; }
                    else { p = p3 * (t2 - iso) / (t2 - t3) + p2 * (iso - t3) / (t2 - t3); }
                    lines.Add(new Line(p1, p));
                    return true;
                }
                return false;
            }
            if (t2 == iso)
            {
                if ((iso >= t1 && iso <= t3) || (iso >= t3 && iso <= t1))
                {
                    Point3d p;
                    if (t1 == t3) { p = (p1 + p3) / 2; }
                    else { p = p3 * (t1 - iso) / (t1 - t3) + p1 * (iso - t3) / (t1 - t3); }
                    lines.Add(new Line(p2, p));
                    return true;
                }
                return false;
            }
            if (t3 == iso)
            {
                if ((iso >= t2 && iso <= t1) || (iso >= t1 && iso <= t2))
                {
                    Point3d p;
                    if (t1 == t2) { p = (p1 + p2) / 2; }
                    p = p2 * (t1 - iso) / (t1 - t2) + p1 * (iso - t2) / (t1 - t2);
                    lines.Add(new Line(p3, p));
                    return true;
                }
                return false;
            }

            int square_idx = 0;
            if (t1 < iso) square_idx |= 1;
            if (t2 < iso) square_idx |= 2;
            if (t3 < iso) square_idx |= 4;
            int a = TriLine[square_idx, 0];
            int b = TriLine[square_idx, 1];
            if (a != -1 && b != -1)
            {
                List<Point3d> L = new List<Point3d>();
                L.Add(p2 * (t1 - iso) / (t1 - t2) + p1 * (iso - t2) / (t1 - t2));
                L.Add(p3 * (t2 - iso) / (t2 - t3) + p2 * (iso - t3) / (t2 - t3));
                L.Add(p1 * (t3 - iso) / (t3 - t1) + p3 * (iso - t1) / (t3 - t1));
                lines.Add(new Line(L[a], L[b]));
                return true;
            }
            return false;
        }
        public Mesh followlines2(Mesh mesh, List<double> t, double iso)
        {
            Mesh ls = new Mesh();
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                if (mesh.Faces[i].IsTriangle)
                {
                    Point3d p1 = mesh.Vertices[mesh.Faces[i].A];
                    Point3d p2 = mesh.Vertices[mesh.Faces[i].B];
                    Point3d p3 = mesh.Vertices[mesh.Faces[i].C];
                    double t1 = t[mesh.Faces[i].A];
                    double t2 = t[mesh.Faces[i].B];
                    double t3 = t[mesh.Faces[i].C];
                    Solve3Face(p1, p2, p3, t1, t2, t3, iso, ref ls);
                }
            }
            return ls;
        }
        public bool Solve3Face(Point3d p1, Point3d p2, Point3d p3, double t1, double t2, double t3, double iso, ref Mesh lines)
        {


            int square_idx = 0;
            if (t1 < iso) square_idx |= 1;
            if (t2 < iso) square_idx |= 2;
            if (t3 < iso) square_idx |= 4;
            int a = TriLine[square_idx, 0];
            int b = TriLine[square_idx, 1];
            int c = TriLine[square_idx, 2];
            int d = TriLine[square_idx, 3];
            if (a != -1)
            {
                Mesh mesh = new Mesh();
                if (a == -2)
                {
                    mesh.Vertices.Add(p1); mesh.Vertices.Add(p2); mesh.Vertices.Add(p3);
                    mesh.Faces.AddFace(0, 1, 2);
                    lines.Append(mesh);
                    return true;
                }
                List<Point3d> L = new List<Point3d>();
                L.Add(p2 * (t1 - iso) / (t1 - t2) + p1 * (iso - t2) / (t1 - t2));
                L.Add(p3 * (t2 - iso) / (t2 - t3) + p2 * (iso - t3) / (t2 - t3));
                L.Add(p1 * (t3 - iso) / (t3 - t1) + p3 * (iso - t1) / (t3 - t1));

                List<Point3d> L2 = new List<Point3d>();
                L2.Add(p1); L2.Add(p2); L2.Add(p3);
                mesh.Vertices.Add(L[a]);
                mesh.Vertices.Add(L[b]);
                mesh.Vertices.Add(L2[c]);
                if (d != -1)
                {
                    mesh.Vertices.Add(L2[d]);
                    mesh.Faces.AddFace(0, 1, 2, 3);
                }
                else { mesh.Faces.AddFace(0, 1, 2); }
                lines.Append(mesh);
                return true;
            }
            return false;
        }


        static int[,] TriLine = {
    {-1, -1} ,
    { 0, 2} ,
    { 0, 1} ,
    { 1, 2} ,
    { 1, 2} ,
    { 0, 1} ,
    { 0, 2} ,
    {-1, -1} ,
    };
    }


}