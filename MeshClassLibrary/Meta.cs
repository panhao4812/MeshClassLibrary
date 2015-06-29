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
    public static class MeshTopoVerticeConvert
    {
        static List<double> Data_Vertices2TopoVertice(Mesh mesh, List<double> t)
        {
            List<double> output = new List<double>();
            if (mesh.Vertices.Count != t.Count) return output;
             if (mesh.Vertices.Count != mesh.Vertices.Count) return output;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            for (int i = 0; i < vs.Count; i++)
            {
                int[] index = vs.MeshVertexIndices(i);
                if (index.Length < 1) return output;
                output.Add(t[index[0]]);
            }
            return output;
        }
        static List<double> Data_TopoVertices2Vertice(Mesh mesh, List<double> t)
        {
            List<double> output=new List<double>();
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            if (vs.Count != t.Count) return output;
            for (int i = 0; i < mesh.Vertices.Count; i++)
            {
                output.Add(0);
            }
            for (int i = 0; i < vs.Count; i++)
            {
                int[] index = vs.MeshVertexIndices(i);
                 if (index.Length < 1) return output;
                 for (int j = 0; j < index.Length; j++)
                 {
                     output[index[j]] = t[i];
                 }
            }
            return output;
        }
    }
    public class MeshFollowLines
    {
        public MeshFollowLines() { }
        public  List<Line> followlines(Mesh mesh, List<double> t, double iso)
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
                if (mesh.Faces[i].IsQuad)
                {
                    Point3d p1 = mesh.Vertices[mesh.Faces[i].A];
                    Point3d p2 = mesh.Vertices[mesh.Faces[i].B];
                    Point3d p3 = mesh.Vertices[mesh.Faces[i].C];
                    Point3d p4 = mesh.Vertices[mesh.Faces[i].D];
                    double t1 = t[mesh.Faces[i].A];
                    double t2 = t[mesh.Faces[i].B];
                    double t3 = t[mesh.Faces[i].C];
                    double t4 = t[mesh.Faces[i].D];
                    Solve4Face(p1, p2, p3, p4, t1, t2, t3, t4, iso, ref ls);
                }
            }
            return ls;
        }
        public bool Solve4Face(Point3d p1, Point3d p2, Point3d p3, Point3d p4, double t1, double t2, double t3, double t4,double iso, ref List<Line> lines)
        {
            int square_idx = 0;
            if (t1 < iso) square_idx |= 1;
            if (t2 < iso) square_idx |= 2;
            if (t3 < iso) square_idx |= 4;
            if (t4 < iso) square_idx |= 8;
            int a = QuadLine[square_idx, 0];
            int b = QuadLine[square_idx, 1];
            int c = QuadLine[square_idx, 2];
            int d = QuadLine[square_idx, 3];

            if (a != -1 && b != -1)
            {
                List<Point3d> L = new List<Point3d>();
                L.Add(p2 * (t1 - iso) / (t1 - t2) + p1 * (iso - t2) / (t1 - t2));
                L.Add(p2 * (t3 - iso) / (t3 - t2) + p3 * (iso - t2) / (t3 - t2));
                L.Add(p3 * (t4 - iso) / (t4 - t3) + p4 * (iso - t3) / (t4 - t3));
                L.Add(p1 * (t4 - iso) / (t4 - t1) + p4 * (iso - t1) / (t4 - t1));
                lines.Add(new Line(L[a], L[b]));
                if (c != -1 && d != -1)
                {
                    lines.Add(new Line(L[c], L[d]));
                }
                return true;
            }
            return false;
        }
        public bool Solve3Face(Point3d p1, Point3d p2, Point3d p3, double t1, double t2, double t3,double iso, ref List<Line> lines)
        {
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
        #region table
        static int[,] QuadLine = {
    {-1, -1, -1, -1 } ,
    {0,  3, -1, -1 },
    {0,  1, -1, -1 },
    {3,  1, -1, -1},
    {1,  2, -1, -1} ,
    {1,  2,  0,  3 },
    {0,  2, -1, -1 },
    {3,  2, -1, -1 },
    {3,  2, -1, -1 },
    {0,  2, -1, -1 },
    {3,  2,  0,  1 },
    {1,  2, -1, -1 },
    {3,  1, -1, -1 },
    {0,  1, -1, -1 },
    {0,  3, -1, -1 },
    {-1, -1, -1, -1}
    };
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
        #endregion
    }
    public class TriangleMeshFollowLines
    {
        public TriangleMeshFollowLines() { }
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
        public Mesh[] followlines2(Mesh mesh, List<double> t, double iso)
        {
            Mesh ls = new Mesh(); Mesh ls2 = new Mesh();
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
                    Solve3Face(p1, p2, p3, t1, t2, t3, iso, ref ls, ref ls2);
                }
            }
            Mesh[] output = { ls, ls2 };
            return output;
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
        public bool Solve3Face(Point3d p1, Point3d p2, Point3d p3, double t1, double t2, double t3, double iso, ref Mesh mesh, ref Mesh mesh2)
        {
            int n = mesh.Vertices.Count();
            int n2 = mesh2.Vertices.Count();
            int square_idx = 0;
            if (t1 < iso)
                square_idx |= 1;
            if (t2 < iso)
                square_idx |= 2;
            if (t3 < iso)
                square_idx |= 4;
            int a = TriLine2[square_idx, 0];
            int b = TriLine2[square_idx, 1];
            int c = TriLine2[square_idx, 2];
            int d = TriLine2[square_idx, 3];
            int e = TriLine2[square_idx, 4];
            int f = TriLine2[square_idx, 5];
            int g = TriLine2[square_idx, 6];
            int h = TriLine2[square_idx, 7];
            if (a == -1)
            {
                mesh2.Vertices.Add(p1);
                mesh2.Vertices.Add(p2);
                mesh2.Vertices.Add(p3);
                mesh2.Faces.AddFace(new MeshFace(0 + n2, 1 + n2, 2 + n2));
                return true;
            }
            else if (a == -2)
            {
                mesh.Vertices.Add(p1);
                mesh.Vertices.Add(p2);
                mesh.Vertices.Add(p3);
                mesh.Faces.AddFace(new MeshFace(0 + n, 1 + n, 2 + n));
                return true;
            }
            else
            {
                Point3d[] L = new Point3d[3];
                L[0] = p2 * (t1 - iso) / (t1 - t2) + p1 * (iso - t2) / (t1 - t2);
                L[1] = p3 * (t2 - iso) / (t2 - t3) + p2 * (iso - t3) / (t2 - t3);
                L[2] = p1 * (t3 - iso) / (t3 - t1) + p3 * (iso - t1) / (t3 - t1);
                Point3d[] L2 = new Point3d[3];
                L2[0] = p1;
                L2[1] = p2;
                L2[2] = p3;
                mesh.Vertices.Add(L[a]);
                mesh.Vertices.Add(L[b]);
                mesh.Vertices.Add(L2[c]);
                if (d != c) { mesh.Vertices.Add(L2[d]); mesh.Faces.AddFace(new MeshFace(0 + n, 1 + n, 2 + n, 3 + n)); }
                else { mesh.Faces.AddFace(new MeshFace(0 + n, 1 + n, 2 + n)); }
                mesh2.Vertices.Add(L[e]);
                mesh2.Vertices.Add(L[f]);
                mesh2.Vertices.Add(L2[g]);
                if (h != g) { mesh2.Vertices.Add(L2[h]); mesh2.Faces.AddFace(new MeshFace(0 + n2, 1 + n2, 2 + n2, 3 + n2)); }
                else { mesh2.Faces.AddFace(new MeshFace(0 + n2, 1 + n2, 2 + n2)); }
                return true;
            }
        }
        public Mesh[] followlines3(Mesh mesh, List<double> t, double iso)
        {
            Mesh[] mesh1 = new Mesh[mesh.Faces.Count];
            Mesh[] mesh2 = new Mesh[mesh.Faces.Count];
            System.Threading.Tasks.Parallel.For(0, mesh.Faces.Count, (i) =>
            {
                if (mesh.Faces[i].IsTriangle)
                {
                    Point3d p1 = mesh.Vertices[mesh.Faces[i].A];
                    Point3d p2 = mesh.Vertices[mesh.Faces[i].B];
                    Point3d p3 = mesh.Vertices[mesh.Faces[i].C];
                    double t1 = t[mesh.Faces[i].A];
                    double t2 = t[mesh.Faces[i].B];
                    double t3 = t[mesh.Faces[i].C];
                    Mesh[] ls = Solve3FaceMul(p1, p2, p3, t1, t2, t3, iso);
                    mesh1[i] = ls[0];
                    mesh2[i] = ls[1];
                }
            });
            Mesh out1 = new Mesh(), out2 = new Mesh();
            for (int i = 0; i < mesh1.Length; i++)
            {
                if (mesh1[i].Faces.Count != 0) out1.Append(mesh1[i]);
            }
            for (int i = 0; i < mesh2.Length; i++)
            {
                if (mesh2[i].Faces.Count != 0) out2.Append(mesh2[i]);
            }
            Mesh[] output = { out1, out2 };
            return output;
        }
        public Mesh[] Solve3FaceMul(Point3d p1, Point3d p2, Point3d p3, double t1, double t2, double t3, double iso)
        {
            Mesh mesh = new Mesh();
            Mesh mesh2 = new Mesh();
            int square_idx = 0;
            if (t1 < iso)
                square_idx |= 1;
            if (t2 < iso)
                square_idx |= 2;
            if (t3 < iso)
                square_idx |= 4;
            int a = TriLine2[square_idx, 0];
            int b = TriLine2[square_idx, 1];
            int c = TriLine2[square_idx, 2];
            int d = TriLine2[square_idx, 3];
            int e = TriLine2[square_idx, 4];
            int f = TriLine2[square_idx, 5];
            int g = TriLine2[square_idx, 6];
            int h = TriLine2[square_idx, 7];
            if (a == -1)
            {
                mesh2.Vertices.Add(p1);
                mesh2.Vertices.Add(p2);
                mesh2.Vertices.Add(p3);
                mesh2.Faces.AddFace(new MeshFace(0, 1, 2));
                Mesh[] output = { mesh, mesh2 };
                return output;
            }
            else if (a == -2)
            {
                mesh.Vertices.Add(p1);
                mesh.Vertices.Add(p2);
                mesh.Vertices.Add(p3);
                mesh.Faces.AddFace(new MeshFace(0, 1, 2));
                Mesh[] output = { mesh, mesh2 };
                return output;
            }
            else
            {
                List<Point3d> L = new List<Point3d>();
                L.Add(p2 * (t1 - iso) / (t1 - t2) + p1 * (iso - t2) / (t1 - t2));
                L.Add(p3 * (t2 - iso) / (t2 - t3) + p2 * (iso - t3) / (t2 - t3));
                L.Add(p1 * (t3 - iso) / (t3 - t1) + p3 * (iso - t1) / (t3 - t1));
                List<Point3d> L2 = new List<Point3d>();
                L2.Add(p1);
                L2.Add(p2);
                L2.Add(p3);
                mesh.Vertices.Add(L[a]);
                mesh.Vertices.Add(L[b]);
                mesh.Vertices.Add(L2[c]);
                if (d != c) { mesh.Vertices.Add(L2[d]); mesh.Faces.AddFace(new MeshFace(0, 1, 2, 3)); }
                else { mesh.Faces.AddFace(new MeshFace(0, 1, 2)); }
                mesh2.Vertices.Add(L[e]);
                mesh2.Vertices.Add(L[f]);
                mesh2.Vertices.Add(L2[g]);
                if (h != g) { mesh2.Vertices.Add(L2[h]); mesh2.Faces.AddFace(new MeshFace(0, 1, 2, 3)); }
                else { mesh2.Faces.AddFace(new MeshFace(0, 1, 2)); }
                Mesh[] output = { mesh, mesh2 };
                return output;
            }
        }
        #region table
        readonly  int[,] TriLine = {
    {-1, -1} ,
    { 0, 2} ,
    { 0, 1} ,
    { 1, 2} ,
    { 1, 2} ,
    { 0, 1} ,
    { 0, 2} ,
    {-1, -1} ,
    };
      readonly int[,] TriLine2 = { 
		{ -1, -1, -1, -1, -2, -2, -2, -2 },
		{ 0, 2, 0, 0, 2, 0, 1, 2 }, 
		{ 1, 0, 1, 1, 0, 1, 2, 0 },
		{ 1, 2, 0, 1, 2, 1, 2, 2 }, 
		{ 2, 1, 2, 2, 1, 2, 0, 1 },
		{ 0, 1, 2, 0, 1, 0, 1,1 }, 
		{ 2, 0, 1, 2, 0, 2, 0, 0 },
		{ -2, -2, -2, -2, -1, -1, -1, -1 }, };
# endregion
    }
    /// ////////////////////////////////
    public class MetaBall
    {
        //Old algorithm and there are some problems solving more than 100^3 box;
        int X, Y, Z;
        double[, ,] energy;
        Point3d[, , ,] EdgePoints;
        List<box> boxes = new List<box>();
        public static int initpt(int a, int b, int c, int maxa, int maxb, int maxc)
        {
            if (a < 0) return -1; if (b < 0) return -1; if (c < 0) return -1;
            if (a >= maxa) return -1;
            if (b >= maxb) return -1;
            if (c >= maxc) return -1;
            return c + b * maxc + a * maxc * maxb;
        }
        public static int initedge(int a, int b, int c, int d, int maxa, int maxb, int maxc)
        {   
            if (a < 0) return -1; if (b < 0) return -1; if (c < 0) return -1;
            if (a >= maxa) return -1;
            if (b >= maxb) return -1;
            if (c >= maxc) return -1;
            if (d < 0 || d >= 3) return -1;
            return d + c * 3 + b * maxc * 3 + a * maxc * maxb * 3;
        }
        public MetaBall(int x, int y, int z)
        {
            X = x; Y = y; Z = z;
            energy = new double[(x + 1), (y + 1), (z + 1)];
            EdgePoints = new Point3d[(x + 1), (y + 1), (z + 1), 3];
            for (int i = 0; i < x; i++)
            {
                for (int j = 0; j < y; j++)
                {
                    for (int k = 0; k < z; k++)
                    {
                        boxes.Add(new box(i, j, k));
                    }
                }
            }
        }
        public void AddEnergy(int x, int y, int z, double t)
        {
            energy[x, y, z] = t;
        }
        public void AddEnergy(double[, ,] t)
        {
            List<box> boxes = new List<box>();
            for (int i = 0; i <= X; i++)
            {
                for (int j = 0; j <= Y; j++)
                {
                    for (int k = 0; k <= Z; k++)
                    {
                        energy[i, j, k] = t[i, j, k];
                    }
                }
            }
        }
        public void AddEnergy(List<double> t)
        {
            int c = 0;
            for (int i = 0; i <= X; i++)
            {
                for (int j = 0; j <= Y; j++)
                {
                    for (int k = 0; k <= Z; k++)
                    {
                        if (c >= t.Count) { energy[i, j, k] = 0; }
                        else { energy[i, j, k] = t[c]; }
                        c++;
                    }
                }
            }
        }
        public Mesh Compute(double Iso)
        {
            Mesh meshoutput = new Mesh();
            double iso = Iso;
            for (int i = 0; i <= X; i++)
            {
                for (int j = 0; j <= Y; j++)
                {
                    for (int k = 0; k <= Z; k++)
                    {
                        if (i < X) { EdgePoints[i, j, k, 0] = caculateCp(i, j, k, i + 1, j, k, energy[i, j, k], energy[i + 1, j, k], iso); }
                        if (j < Y) { EdgePoints[i, j, k, 1] = caculateCp(i, j, k, i, j + 1, k, energy[i, j, k], energy[i, j + 1, k], iso); }
                        if (k < Z) { EdgePoints[i, j, k, 2] = caculateCp(i, j, k, i, j, k + 1, energy[i, j, k], energy[i, j, k + 1], iso); }
                    }
                }
            }
            for (int i = 0; i < boxes.Count; i++)
            {
                box b = boxes[i];
                Mesh mesh = new Mesh();
                int cubeindex = 0;
                if (energy[b.points[0, 0], b.points[0, 1], b.points[0, 2]] < iso) { cubeindex |= 1; }//;str += "1";}else{str += "0";}
                if (energy[b.points[1, 0], b.points[1, 1], b.points[1, 2]] < iso) { cubeindex |= 2; }//;str += "1";}else{str += "0";}
                if (energy[b.points[2, 0], b.points[2, 1], b.points[2, 2]] < iso) { cubeindex |= 4; }//;str += "1";}else{str += "0";}
                if (energy[b.points[3, 0], b.points[3, 1], b.points[3, 2]] < iso) { cubeindex |= 8; }//;str += "1";}else{str += "0";}
                if (energy[b.points[4, 0], b.points[4, 1], b.points[4, 2]] < iso) { cubeindex |= 16; }//;str += "1";}else{str += "0";}
                if (energy[b.points[5, 0], b.points[5, 1], b.points[5, 2]] < iso) { cubeindex |= 32; }//;str += "1";}else{str += "0";}
                if (energy[b.points[6, 0], b.points[6, 1], b.points[6, 2]] < iso) { cubeindex |= 64; }//;str += "1";}else{str += "0";}
                if (energy[b.points[7, 0], b.points[7, 1], b.points[7, 2]] < iso) { cubeindex |= 128; }//;str += "1";}else{str += "0";}
                int vs = 0;
                for (int l = 0; l <= 4; l++)
                {
                    if (triTable[cubeindex, l * 3] == -1) break;
                    int tr1 = triTable[cubeindex, l * 3];
                    int tr2 = triTable[cubeindex, l * 3 + 1];
                    int tr3 = triTable[cubeindex, l * 3 + 2];
                    Point3d p1 = EdgePoints[b.Cedges[tr1, 0], b.Cedges[tr1, 1], b.Cedges[tr1, 2], b.Cedges[tr1, 3]];
                    Point3d p2 = EdgePoints[b.Cedges[tr2, 0], b.Cedges[tr2, 1], b.Cedges[tr2, 2], b.Cedges[tr2, 3]];
                    Point3d p3 = EdgePoints[b.Cedges[tr3, 0], b.Cedges[tr3, 1], b.Cedges[tr3, 2], b.Cedges[tr3, 3]];
                    mesh.Vertices.Add(p1);
                    mesh.Vertices.Add(p2);
                    mesh.Vertices.Add(p3);
                    mesh.Faces.AddFace(vs * 3, vs * 3 + 1, vs * 3 + 2);
                    vs++;
                }
                meshoutput.Append(mesh);
            }
            return meshoutput;
        }
        private Point3d caculateCp(int x1, int x2, int x3,
          int y1, int y2, int y3,
          double e1, double e2, double level)
        {
            if ((e1 < level & e2 < level) | (e1 > level & e2 > level)) return default(Point3d);
            double mul = ((level - e1) / (e2 - e1));
            double p1x = (double)x1;
            double p1y = (double)x2;
            double p1z = (double)x3;
            double p2x = (double)y1;
            double p2y = (double)y2;
            double p2z = (double)y3;
            double x = p1x + (p2x - p1x) * mul;
            double y = p1y + (p2y - p1y) * mul;
            double z = p1z + (p2z - p1z) * mul;
            return new Point3d(x, y, z);
        }
        private Point3d caculateCp(Point3d p1, Point3d p2,
          double e1, double e2, double level)
        {
            if ((e1 < level & e2 < level) | (e1 > level & e2 > level)) return new Point3d();
            double mul = ((level - e1) / (e2 - e1));
            return p1 + (p2 - p1) * mul;
        }
        #region table
       
       readonly int[,] triTable =
      {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
      {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
      {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
      {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
      {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
      {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
      {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
      {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
      {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
      {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
      {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
      {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
      {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
      {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
      {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
      {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
      {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
      {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
      {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
      {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
      {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
      {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
      {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
      {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
      {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
      {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
      {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
      {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
      {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
      {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
      {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
      {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
      {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
      {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
      {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
      {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
      {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
      {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
      {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
      {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
      {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
      {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
      {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
      {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
      {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
      {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
      {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
      {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
      {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
      {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
      {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
      {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
      {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
      {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
      {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
      {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
      {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
      {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
      {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
      {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
      {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
      {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
      {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
      {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
      {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
      {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
      {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
      {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
      {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
      {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
      {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
      {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
      {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
      {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
      {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
      {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
      {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
      {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
      {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
      {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
      {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
      {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
      {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
      {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
      {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
      {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
      {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
      {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
      {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
      {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
      {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
      {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
      {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
      {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
      {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
      {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
      {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
      {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
      {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
      {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
      {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
      {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
      {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
      {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
      {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
      {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
      {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
      {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
      {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
      {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
      {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
      {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
      {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
      {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
      {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
      {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
      {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
      {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
      {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
      {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
      {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
      {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
      {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
      {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
      {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
      {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
      {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
      {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
      {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
      {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
      {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
      {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
      {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
      {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
      {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
      {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
      {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
      {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
      {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
      {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
      {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
      {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
      {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
      {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
      {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
      {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
      {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
      {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
      {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
      {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
      {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
      {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
      {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
      {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
      {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
      {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
      {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
      {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
      {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
      {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
      {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
      {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
      {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
      {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
      {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
      {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
      {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
      {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
      {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
      {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
      {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
      {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
      {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
      {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
      {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
      {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
      {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
      {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
      {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};
        #endregion
    }
    public struct box
    {
        public int[,] Cedges;// = new int[12,4];
        public int[,] points;// = new int[8,3];
        public box(int i, int j, int k)
        {
            int[,] Temp1 = {{i,j + 1,k},{i + 1,j + 1,k},{i + 1,j,k},{i,j,k},
          {i,j + 1,k + 1},{i + 1,j + 1,k + 1},{i + 1,j,k + 1},{i,j,k + 1}
          };
            int[,] Temp2 = {{i,j + 1,k,0}, {i + 1,j,k,1},{i,j,k,0},{i,j,k,1},
          {i,j + 1,k + 1,0},{i + 1,j,k + 1,1},{i,j,k + 1,0},{i,j,k + 1,1},
          {i,j + 1,k,2},{i + 1,j + 1,k,2},{i + 1,j,k,2},{i,j,k,2}
          };
            this.points = Temp1; this.Cedges = Temp2;
        }
    }
    ////////////////////////////////////
    struct BOX
    {
        public Point3f[] triangles;
        public int FaceCount;
    }
    class MetaBallMul
    {
        #region table
        readonly int[] edgeTable = {
      0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
      0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
      0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
      0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
      0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
      0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
      0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
      0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
      0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
      0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
      0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
      0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
      0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
      0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
      0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
      0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
      0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
      0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
      0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
      0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
      0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
      0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
      0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
      0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
      0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
      0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
      0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
      0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
      0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
      0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
      0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
      0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };
        readonly int[,] triTable =
      {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
      {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
      {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
      {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
      {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
      {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
      {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
      {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
      {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
      {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
      {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
      {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
      {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
      {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
      {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
      {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
      {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
      {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
      {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
      {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
      {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
      {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
      {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
      {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
      {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
      {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
      {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
      {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
      {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
      {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
      {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
      {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
      {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
      {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
      {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
      {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
      {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
      {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
      {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
      {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
      {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
      {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
      {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
      {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
      {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
      {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
      {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
      {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
      {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
      {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
      {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
      {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
      {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
      {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
      {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
      {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
      {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
      {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
      {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
      {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
      {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
      {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
      {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
      {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
      {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
      {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
      {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
      {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
      {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
      {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
      {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
      {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
      {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
      {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
      {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
      {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
      {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
      {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
      {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
      {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
      {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
      {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
      {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
      {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
      {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
      {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
      {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
      {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
      {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
      {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
      {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
      {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
      {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
      {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
      {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
      {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
      {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
      {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
      {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
      {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
      {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
      {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
      {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
      {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
      {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
      {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
      {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
      {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
      {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
      {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
      {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
      {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
      {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
      {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
      {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
      {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
      {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
      {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
      {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
      {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
      {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
      {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
      {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
      {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
      {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
      {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
      {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
      {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
      {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
      {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
      {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
      {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
      {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
      {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
      {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
      {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
      {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
      {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
      {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
      {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
      {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
      {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
      {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
      {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
      {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
      {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
      {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
      {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
      {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
      {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
      {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
      {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
      {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
      {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
      {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
      {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
      {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
      {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
      {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
      {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
      {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
      {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
      {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
      {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
      {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
      {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
      {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
      {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
      {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
      {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
      {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
      {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
      {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
      {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
      {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
      {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
      {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
      {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
      {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
      {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
      {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};
        #endregion
        public MetaBallMul(int x, int y, int z)
        {
            boxes = new BOX[x, y, z];
            for (int i = 0; i < x; i++)
            {
                for (int j = 0; j < x; j++)
                {
                    for (int k = 0; k < x; k++)
                    {
                        boxes[i, j, k] = new BOX();
                        boxes[i, j, k].triangles = new Point3f[15];
                    }
                }
            }
        }
        public BOX[, ,] boxes;
        public Mesh DrawMesh()
        {
            Mesh mesh = new Mesh();
            foreach (BOX b in this.boxes)
            {
                // Print(b.FaceCount.ToString());
                for (int ii = 0; ii < b.FaceCount; ii++)
                {
                    int n = mesh.Vertices.Count;
                    mesh.Vertices.Add(b.triangles[ii * 3]);
                    mesh.Vertices.Add(b.triangles[ii * 3 + 1]);
                    mesh.Vertices.Add(b.triangles[ii * 3 + 2]);
                    mesh.Faces.AddFace(n, n + 1, n + 2);
                }
            }
            return mesh;
        }
        public void Compute(int x, int y, int z, int maxx, int maxy, int maxz, double isolevel)
        {
            boxes[x, y, z].FaceCount = 0;
            Point3d[] p = new Point3d[8];
            double[] val = new double[8];
            double i = (double)x - (double)maxx / 2.5;
            double j = (double)y - (double)maxy / 2.5;
            double k = (double)z - (double)maxz / 2.5;
            p[0] = new Point3d(i, j + 1, k);
            p[1] = new Point3d(i + 1, j + 1, k);
            p[2] = new Point3d(i + 1, j, k);
            p[3] = new Point3d(i, j, k);
            p[4] = new Point3d(i, j + 1, k + 1);
            p[5] = new Point3d(i + 1, j + 1, k + 1);
            p[6] = new Point3d(i + 1, j, k + 1);
            p[7] = new Point3d(i, j, k + 1);
            for (int ii = 0; ii < 8; ii++)
            {
                double xa, ya, za;
                xa = p[ii].X; xa *= 2 / (double)maxx;
                ya = p[ii].Y; ya *= 2 / (double)maxy;
                za = p[ii].Z; za *= 2 / (double)maxz;
                val[ii] = Math.Sin(4 * xa) + Math.Sin(4 * ya) + Math.Sin(4 * za) + 4 * xa * ya * za;

            }

            int cubeindex = 0;
            if (val[0] < isolevel) cubeindex |= 1;
            if (val[1] < isolevel) cubeindex |= 2;
            if (val[2] < isolevel) cubeindex |= 4;
            if (val[3] < isolevel) cubeindex |= 8;
            if (val[4] < isolevel) cubeindex |= 16;
            if (val[5] < isolevel) cubeindex |= 32;
            if (val[6] < isolevel) cubeindex |= 64;
            if (val[7] < isolevel) cubeindex |= 128;
            Point3f[] vertlist = new Point3f[12];
            if (edgeTable[cubeindex] == 0) return;
            if ((edgeTable[cubeindex] & 1) != 0) vertlist[0] = VertexInterp(isolevel, p[0], p[1], val[0], val[1]);
            if ((edgeTable[cubeindex] & 2) != 0) vertlist[1] = VertexInterp(isolevel, p[1], p[2], val[1], val[2]);
            if ((edgeTable[cubeindex] & 4) != 0) vertlist[2] = VertexInterp(isolevel, p[2], p[3], val[2], val[3]);
            if ((edgeTable[cubeindex] & 8) != 0) vertlist[3] = VertexInterp(isolevel, p[3], p[0], val[3], val[0]);
            if ((edgeTable[cubeindex] & 16) != 0) vertlist[4] = VertexInterp(isolevel, p[4], p[5], val[4], val[5]);
            if ((edgeTable[cubeindex] & 32) != 0) vertlist[5] = VertexInterp(isolevel, p[5], p[6], val[5], val[6]);
            if ((edgeTable[cubeindex] & 64) != 0) vertlist[6] = VertexInterp(isolevel, p[6], p[7], val[6], val[7]);
            if ((edgeTable[cubeindex] & 128) != 0) vertlist[7] = VertexInterp(isolevel, p[7], p[4], val[7], val[4]);
            if ((edgeTable[cubeindex] & 256) != 0) vertlist[8] = VertexInterp(isolevel, p[0], p[4], val[0], val[4]);
            if ((edgeTable[cubeindex] & 512) != 0) vertlist[9] = VertexInterp(isolevel, p[1], p[5], val[1], val[5]);
            if ((edgeTable[cubeindex] & 1024) != 0) vertlist[10] = VertexInterp(isolevel, p[2], p[6], val[2], val[6]);
            if ((edgeTable[cubeindex] & 2048) != 0) vertlist[11] = VertexInterp(isolevel, p[3], p[7], val[3], val[7]);
            /* Create the triangle */
            for (int ii = 0; triTable[cubeindex, ii] != -1; ii += 3)
            {
                boxes[x, y, z].triangles[ii] = vertlist[triTable[cubeindex, ii]];
                boxes[x, y, z].triangles[ii + 1] = vertlist[triTable[cubeindex, ii + 1]];
                boxes[x, y, z].triangles[ii + 2] = vertlist[triTable[cubeindex, ii + 2]];
                boxes[x, y, z].FaceCount++;
            }

        }
        Point3f VertexInterp(double isolevel, Point3d p1, Point3d p2, double valp1, double valp2)
        {
            double mu;

            if (Math.Abs(isolevel - valp1) < 0.00001)
                return new Point3f((float)p1.X, (float)p1.Y, (float)p1.Z);
            if (Math.Abs(isolevel - valp2) < 0.00001)
                return new Point3f((float)p2.X, (float)p2.Y, (float)p2.Z);
            if (Math.Abs(valp1 - valp2) < 0.00001)
                return new Point3f((float)p1.X, (float)p1.Y, (float)p1.Z);
            Point3f p = new Point3f();
            mu = (isolevel - valp1) / (valp2 - valp1);
            p.X = (float)(p1.X + (p2.X - p1.X) * mu);
            p.Y = (float)(p1.Y + (p2.Y - p1.Y) * mu);
            p.Z = (float)(p1.Z + (p2.Z - p1.Z) * mu);
            return p;
        }
    }
    ////////////////////////////////////
}