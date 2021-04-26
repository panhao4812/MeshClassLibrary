using Rhino.Geometry;

using System;
using System.Linq;
using System.Collections.Generic;
using Grasshopper;
using Grasshopper.Kernel.Data;

namespace MeshClassLibrary
{
    public static class MeshTopoVerticeConvert
    {
        public static List<double> Data_Vertices2TopoVertice(Mesh mesh, List<double> t)
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
        public static List<double> Data_TopoVertices2Vertice(Mesh mesh, List<double> t)
        {
            List<double> output = new List<double>();
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
        public List<Line> followlines(Mesh mesh, List<double> t, double iso)
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
        public bool Solve4Face(Point3d p1, Point3d p2, Point3d p3, Point3d p4, double t1, double t2, double t3, double t4, double iso, ref List<Line> lines)
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
        public bool Solve3Face(Point3d p1, Point3d p2, Point3d p3, double t1, double t2, double t3, double iso, ref List<Line> lines)
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
        readonly int[,] TriLine = {
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
        /*Old algorithm and there are some problems solving more than 100^3 box;
         * use with ScaleField   ScaleField.Q=iso
         * ScaleField.ValueAt()=MetaBall.MathWoof()
         * |override MathWoof()|==>new()==>AddEnergy()==>Compute()
         */
        int X, Y, Z;
        public double[,,] energy;
        Point3d[,,,] EdgePoints;
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
        public virtual double MathWoof(int i, int j, int k)
        {
            return GyroidFunction(i / 5.0, j, k);
        }
        public double GyroidFunction(double i, double j, double k)
        {
            double xa = i, ya = j, za = k;
            //  energy = Math.Sin(4 * xa) + Math.Sin(4 * ya) + Math.Sin(4 * za) + 4 * xa * ya * za;
            return Math.Sin(xa) * Math.Cos(ya) + Math.Sin(ya) * Math.Cos(za) + Math.Sin(za) * Math.Cos(xa);
        }
        public double BlobsFunction(double i, double j, double k, List<Point3d> Pts)
        {
            Point3d P = new Point3d(i, j, k);
            double Sum = 0;
            foreach (Point3d Pt in Pts)
            {
                Sum += (1.0 / (P - Pt).SquareLength);
            }
            return Sum;
        }
        public double CurveFunction(double i, double j, double k, List<Curve> Crv)
        {
            Point3d P = new Point3d(i, j, k);
            double Sum = 0;
            double t;
            foreach (Curve C in Crv)
            {
                C.ClosestPoint(P, out t);
                Point3d Pt = C.PointAt(t);
                Sum += (1.0 / (P - Pt).Length);
            }
            return Sum;
        }
        public static double GetQ(double iso, Interval interval)
        {
            //iso is from 0 to 1;
            if (iso > 1) iso = 1;
            if (iso < 0) iso = 0;
            return iso * (interval.Length) + interval.Min;
        }
        public Interval AddEnergy()
        {
            double max = double.MinValue; double min = double.MaxValue;
            for (int i = 0; i <= energy.GetUpperBound(0); i++)
            {
                for (int j = 0; j <= energy.GetUpperBound(1); j++)
                {
                    for (int k = 0; k <= energy.GetUpperBound(2); k++)
                    {
                        double t1 = MathWoof(i, j, k);
                        if (t1 < min) min = t1;
                        if (t1 > max) max = t1;
                        energy[i, j, k] = t1;
                    }
                }
            }
            return new Interval(min, max);
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
                        if (i < X) { EdgePoints[i, j, k, 0] = VertexInterp(i, j, k, i + 1, j, k, energy[i, j, k], energy[i + 1, j, k], iso); }
                        if (j < Y) { EdgePoints[i, j, k, 1] = VertexInterp(i, j, k, i, j + 1, k, energy[i, j, k], energy[i, j + 1, k], iso); }
                        if (k < Z) { EdgePoints[i, j, k, 2] = VertexInterp(i, j, k, i, j, k + 1, energy[i, j, k], energy[i, j, k + 1], iso); }
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
        private Point3d VertexInterp(int x1, int x2, int x3,
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
        private Point3d VertexInterp(Point3d p1, Point3d p2,
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

        struct box
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
    }
    ////////////////////////////////////
    class MetaBallMul_BOX
    {
        public MetaBallMul_BOX(int x, int y, int z, int max1, int max2, int max3)
        {
            double i = (double)x - (double)max1 / 2;
            double j = (double)y - (double)max2 / 2;
            double k = (double)z - (double)max3 / 2;

            p[0] = new Point3d(i, j + 1, k);
            p[1] = new Point3d(i + 1, j + 1, k);
            p[2] = new Point3d(i + 1, j, k);
            p[3] = new Point3d(i, j, k);
            p[4] = new Point3d(i, j + 1, k + 1);
            p[5] = new Point3d(i + 1, j + 1, k + 1);
            p[6] = new Point3d(i + 1, j, k + 1);
            p[7] = new Point3d(i, j, k + 1);
        }
        public virtual void Compute(double scale)
        {
            for (int i = 0; i < 8; i++)
            {
                double xa, ya, za;
                xa = p[i].X; xa *= 3 / scale;
                ya = p[i].Y; ya *= 3 / scale;
                za = p[i].Z; za *= 3 / scale;
                val[i] = Math.Sin(4 * xa) + Math.Sin(4 * ya) + Math.Sin(4 * za) + 4 * xa * ya * za;
            }
        }
        public Point3d[] p = new Point3d[8];
        public double[] val = new double[8];
        public Point3d[] triangles = new Point3d[15];
        public int FaceCount = 0;
    }
    class MetaBallMul
    {
        public MetaBallMul(int x, int y, int z)
        {
            for (int i = 0; i < x; i++)
            {
                for (int j = 0; j < x; j++)
                {
                    for (int k = 0; k < x; k++)
                    {
                        boxes.Add(new MetaBallMul_BOX(i, j, k, x, y, z));
                    }
                }
            }
        }
        public List<MetaBallMul_BOX> boxes = new List<MetaBallMul_BOX>();
        public Mesh DrawMesh()
        {
            Mesh mesh = new Mesh();
            boxes.ForEach(delegate (MetaBallMul_BOX b)
            {
                for (int i = 0; i < b.FaceCount; i++)
                {
                    int n = mesh.Vertices.Count;
                    mesh.Vertices.Add(b.triangles[i * 3]);
                    mesh.Vertices.Add(b.triangles[i * 3 + 1]);
                    mesh.Vertices.Add(b.triangles[i * 3 + 2]);
                    mesh.Faces.AddFace(n, n + 1, n + 2);
                }
            });
            return mesh;
        }
        public void Compute(int ii, double isolevel)
        {
            MetaBallMul_BOX grid = boxes[ii];
            grid.FaceCount = 0;
            grid.Compute(100);
            int cubeindex = 0;
            if (grid.val[0] < isolevel) cubeindex |= 1;
            if (grid.val[1] < isolevel) cubeindex |= 2;
            if (grid.val[2] < isolevel) cubeindex |= 4;
            if (grid.val[3] < isolevel) cubeindex |= 8;
            if (grid.val[4] < isolevel) cubeindex |= 16;
            if (grid.val[5] < isolevel) cubeindex |= 32;
            if (grid.val[6] < isolevel) cubeindex |= 64;
            if (grid.val[7] < isolevel) cubeindex |= 128;
            Point3d[] vertlist = new Point3d[12];
            if (edgeTable[cubeindex] == 0) return;
            if ((edgeTable[cubeindex] & 1) != 0) vertlist[0] = VertexInterp(isolevel, grid.p[0], grid.p[1], grid.val[0], grid.val[1]);
            if ((edgeTable[cubeindex] & 2) != 0) vertlist[1] = VertexInterp(isolevel, grid.p[1], grid.p[2], grid.val[1], grid.val[2]);
            if ((edgeTable[cubeindex] & 4) != 0) vertlist[2] = VertexInterp(isolevel, grid.p[2], grid.p[3], grid.val[2], grid.val[3]);
            if ((edgeTable[cubeindex] & 8) != 0) vertlist[3] = VertexInterp(isolevel, grid.p[3], grid.p[0], grid.val[3], grid.val[0]);
            if ((edgeTable[cubeindex] & 16) != 0) vertlist[4] = VertexInterp(isolevel, grid.p[4], grid.p[5], grid.val[4], grid.val[5]);
            if ((edgeTable[cubeindex] & 32) != 0) vertlist[5] = VertexInterp(isolevel, grid.p[5], grid.p[6], grid.val[5], grid.val[6]);
            if ((edgeTable[cubeindex] & 64) != 0) vertlist[6] = VertexInterp(isolevel, grid.p[6], grid.p[7], grid.val[6], grid.val[7]);
            if ((edgeTable[cubeindex] & 128) != 0) vertlist[7] = VertexInterp(isolevel, grid.p[7], grid.p[4], grid.val[7], grid.val[4]);
            if ((edgeTable[cubeindex] & 256) != 0) vertlist[8] = VertexInterp(isolevel, grid.p[0], grid.p[4], grid.val[0], grid.val[4]);
            if ((edgeTable[cubeindex] & 512) != 0) vertlist[9] = VertexInterp(isolevel, grid.p[1], grid.p[5], grid.val[1], grid.val[5]);
            if ((edgeTable[cubeindex] & 1024) != 0) vertlist[10] = VertexInterp(isolevel, grid.p[2], grid.p[6], grid.val[2], grid.val[6]);
            if ((edgeTable[cubeindex] & 2048) != 0) vertlist[11] = VertexInterp(isolevel, grid.p[3], grid.p[7], grid.val[3], grid.val[7]);
            /* Create the triangle */

            for (int i = 0; triTable[cubeindex, i] != -1; i += 3)
            {
                grid.triangles[i] = vertlist[triTable[cubeindex, i]];
                grid.triangles[i + 1] = vertlist[triTable[cubeindex, i + 1]];
                grid.triangles[i + 2] = vertlist[triTable[cubeindex, i + 2]];
                grid.FaceCount++;
            }
        }
        Point3d VertexInterp(double isolevel, Point3d p1, Point3d p2, double valp1, double valp2)
        {
            double mu;
            Point3d p;
            if (Math.Abs(isolevel - valp1) < 0.00001)
                return (p1);
            if (Math.Abs(isolevel - valp2) < 0.00001)
                return (p2);
            if (Math.Abs(valp1 - valp2) < 0.00001)
                return (p1);
            mu = (isolevel - valp1) / (valp2 - valp1);
            p = p1 + (p2 - p1) * mu;
            return p;
        }
        #region table
        int[] edgeTable = {
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
        int[,] triTable =
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
    ////////////////////
    /*
 
    var rangePartitioner = System.Collections.Concurrent.Partitioner.Create(0, meta.boxes.Count);

    System.Threading.Tasks.Parallel.ForEach(rangePartitioner, (range, loopState) =>
      {
      for (int i = range.Item1; i < range.Item2; i++)
      {
        meta.Compute(i, iso);
      }});
    return meta.DrawMesh();
  
    */
    /*//If you want to have a go at defining your own fields, 
      * you’ll also need to include the gradient function.
      * For this I find Wolfram Alpha comes in very handy – you can just plug in Grad
      * (some function of xyz) and get the result. 
      * For example: grad(sin x * cos y + sin y * cos z + sin z * cos x)
       //*/
    public abstract class ScalarField
    {
        public double Q = 1;
        public void MeshPush(ref Mesh mesh)
        {

            List<Point3d> Pt = new List<Point3d>();
            for (int i = 0; i < mesh.Vertices.Count; i++)
            {
                Pt.Add(new Point3d(mesh.Vertices[i]));
            }
            System.Threading.Tasks.Parallel.For(0, Pt.Count, j =>
            {
                if (Pt[j].IsValid)
                {
                    {
                        double FValue = 1.0;
                        int counter = 0;
                        while (Math.Abs(FValue) > Acc && counter < Iter) //accuracy and max iterations
                        {
                            FValue = this.ValueAt(Pt[j]);
                            Vector3d Gradient = this.GradientAt(Pt[j]);
                            Gradient *= 1.0 / Gradient.SquareLength;
                            Pt[j] = Pt[j] - FValue * Gradient;
                            counter++;
                        }
                    }
                }
            });

            for (int i = 0; i < mesh.Vertices.Count; i++)
            {
                mesh.Vertices[i] = new Point3f((float)Pt[i].X, (float)Pt[i].Y, (float)Pt[i].Z);
            }
        }
        private double Acc = 0.001;
        public int Iter = 15;
        public abstract double ValueAt(Point3d p);
        public abstract Vector3d GradientAt(Point3d p);
    }
    public class ScalarField_Gyroid : ScalarField
    {
        public ScalarField_Gyroid() { }
        public ScalarField_Gyroid(double _Q)
        {
            Q = _Q;
        }
        public override double ValueAt(Point3d P)
        {
            return
              (Math.Sin(P.X) * Math.Cos(P.Y) +
              Math.Sin(P.Y) * Math.Cos(P.Z) +
              Math.Sin(P.Z) * Math.Cos(P.X)) + Q;
        }
        public override Vector3d GradientAt(Point3d P)
        {
            double x = P.X;
            double y = P.Y;
            double z = P.Z;
            return
              new Vector3d(
              Math.Cos(x) * Math.Cos(y) - Math.Sin(x) * Math.Sin(z),
              Math.Cos(y) * Math.Cos(z) - Math.Sin(x) * Math.Sin(y),
              Math.Cos(x) * Math.Cos(z) - Math.Sin(y) * Math.Sin(z));
        }
    }
    public class ScalarField_Blobs : ScalarField
    {
        public ScalarField_Blobs() { }
        public ScalarField_Blobs(double _Q, List<Point3d> P)
        {
            Q = _Q;
            Pts = P;
        }
        public List<Point3d> Pts = new List<Point3d>();
        public override double ValueAt(Point3d P)
        {
            double Sum = 0;
            foreach (Point3d Pt in Pts)
            {
                Sum += (1.0 / (P - Pt).SquareLength);
            }
            return Sum - Q;
        }
        public override Vector3d GradientAt(Point3d P)
        {
            var VSum = new Vector3d();
            foreach (Point3d Pt in Pts)
            {
                Vector3d V = P - Pt;
                double Value = Math.Pow((V.SquareLength), 2);
                VSum -= 2 * V / Value;
            }
            return VSum;
        }
    }
    public class ScalarField_Curve : ScalarField
    {
        public ScalarField_Curve() { }
        public List<Curve> Crv;
        public ScalarField_Curve(double _Q, List<Curve> C)
        {
            Q = _Q;
            Crv = C;
        }
        public override double ValueAt(Point3d P)
        {
            double Sum = 0;
            double t;
            foreach (Curve C in Crv)
            {
                C.ClosestPoint(P, out t);
                Point3d Pt = C.PointAt(t);
                Sum += (1.0 / (P - Pt).Length);
            }
            return Sum - Q;
        }
        public override Vector3d GradientAt(Point3d P)
        {
            double x = P.X;
            double y = P.Y;
            double z = P.Z;
            double t;
            var VSum = new Vector3d();
            foreach (Curve C in Crv)
            {
                C.ClosestPoint(P, out t);
                Point3d Pt = C.PointAt(t);
                Vector3d V = P - Pt;
                double Value = Math.Pow((V.SquareLength), 3 / 2);
                VSum -= new Vector3d(V.X / Value, V.Y / Value, V.Z / Value);
            }
            return VSum;
        }
    }
    public class Meta_Triangle
    {

        public class Matrix_Octahedra
        {
            public List<Mesh> mesh1 = new List<Mesh>();
            public List<Mesh> mesh2 = new List<Mesh>();
            public Matrix_Octahedra(int XCells, int YCells, int ZCells, double CellSize)
            {
                double CS = CellSize;
                double CH = 0.5 * CS;
                double CQ = 0.25 * CS;
                var Grids = new Point3d[2][,,];
                Grids[0] = new Point3d[XCells, YCells, ZCells];
                Grids[1] = new Point3d[XCells, YCells, ZCells];
                for (int i = 0; i < XCells; i++)
                {
                    for (int j = 0; j < YCells; j++)
                    {
                        for (int k = 0; k < ZCells; k++)
                        {
                            for (int g = 0; g < 2; g++)
                            {
                                Point3d P = new Point3d(i * CS, j * CS, k * CS) + g * new Vector3d(CH, CH, CH);
                                Mesh M = new Mesh();
                                Grids[g][i, j, k] = P;
                                #region mesh construct
                                int VCount;
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CH, P.Y + CQ, P.Z + 00);
                                M.Vertices.Add(P.X + CH, P.Y + 00, P.Z + CQ);
                                M.Vertices.Add(P.X + CH, P.Y - CQ, P.Z + 00);
                                M.Vertices.Add(P.X + CH, P.Y + 00, P.Z - CQ);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 3);
                                M.Faces.AddFace(VCount + 1, VCount + 2, VCount + 3);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y + CH, P.Z + 00);
                                M.Vertices.Add(P.X + 00, P.Y + CH, P.Z + CQ);
                                M.Vertices.Add(P.X - CQ, P.Y + CH, P.Z + 00);
                                M.Vertices.Add(P.X + 00, P.Y + CH, P.Z - CQ);
                                M.Faces.AddFace(VCount + 3, VCount + 1, VCount);
                                M.Faces.AddFace(VCount + 3, VCount + 2, VCount + 1);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y + 00, P.Z + CH);
                                M.Vertices.Add(P.X + 00, P.Y + CQ, P.Z + CH);
                                M.Vertices.Add(P.X - CQ, P.Y + 00, P.Z + CH);
                                M.Vertices.Add(P.X + 00, P.Y - CQ, P.Z + CH);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 3);
                                M.Faces.AddFace(VCount + 1, VCount + 2, VCount + 3);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X - CH, P.Y + CQ, P.Z + 00);
                                M.Vertices.Add(P.X - CH, P.Y + 00, P.Z + CQ);
                                M.Vertices.Add(P.X - CH, P.Y - CQ, P.Z + 00);
                                M.Vertices.Add(P.X - CH, P.Y + 00, P.Z - CQ);
                                M.Faces.AddFace(VCount + 3, VCount + 1, VCount);
                                M.Faces.AddFace(VCount + 3, VCount + 2, VCount + 1);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y - CH, P.Z + 00);
                                M.Vertices.Add(P.X + 00, P.Y - CH, P.Z + CQ);
                                M.Vertices.Add(P.X - CQ, P.Y - CH, P.Z + 00);
                                M.Vertices.Add(P.X + 00, P.Y - CH, P.Z - CQ);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 3);
                                M.Faces.AddFace(VCount + 1, VCount + 2, VCount + 3);

                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y + 00, P.Z - CH);
                                M.Vertices.Add(P.X + 00, P.Y + CQ, P.Z - CH);
                                M.Vertices.Add(P.X - CQ, P.Y + 00, P.Z - CH);
                                M.Vertices.Add(P.X + 00, P.Y - CQ, P.Z - CH);
                                M.Faces.AddFace(VCount + 3, VCount + 1, VCount);
                                M.Faces.AddFace(VCount + 3, VCount + 2, VCount + 1);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y + CQ, P.Z + CQ);
                                M.Vertices.Add(P.X + CH, P.Y + CQ, P.Z + 00);
                                M.Vertices.Add(P.X + CQ, P.Y + CH, P.Z + 00);
                                M.Vertices.Add(P.X + 00, P.Y + CH, P.Z + CQ);
                                M.Vertices.Add(P.X + 00, P.Y + CQ, P.Z + CH);
                                M.Vertices.Add(P.X + CQ, P.Y + 00, P.Z + CH);
                                M.Vertices.Add(P.X + CH, P.Y + 00, P.Z + CQ);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 6);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 1);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X - CQ, P.Y + CQ, P.Z + CQ);
                                M.Vertices.Add(P.X - CQ, P.Y + CH, P.Z + 00);
                                M.Vertices.Add(P.X - CH, P.Y + CQ, P.Z + 00);
                                M.Vertices.Add(P.X - CH, P.Y + 00, P.Z + CQ);
                                M.Vertices.Add(P.X - CQ, P.Y + 00, P.Z + CH);
                                M.Vertices.Add(P.X + 00, P.Y + CQ, P.Z + CH);
                                M.Vertices.Add(P.X + 00, P.Y + CH, P.Z + CQ);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 6);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 1);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X - CQ, P.Y - CQ, P.Z + CQ);
                                M.Vertices.Add(P.X - CH, P.Y - CQ, P.Z + 00);
                                M.Vertices.Add(P.X - CQ, P.Y - CH, P.Z + 00);
                                M.Vertices.Add(P.X + 00, P.Y - CH, P.Z + CQ);
                                M.Vertices.Add(P.X + 00, P.Y - CQ, P.Z + CH);
                                M.Vertices.Add(P.X - CQ, P.Y + 00, P.Z + CH);
                                M.Vertices.Add(P.X - CH, P.Y + 00, P.Z + CQ);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 6);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 1);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y - CQ, P.Z + CQ);
                                M.Vertices.Add(P.X + CQ, P.Y - CH, P.Z + 00);
                                M.Vertices.Add(P.X + CH, P.Y - CQ, P.Z + 00);
                                M.Vertices.Add(P.X + CH, P.Y + 00, P.Z + CQ);
                                M.Vertices.Add(P.X + CQ, P.Y + 00, P.Z + CH);
                                M.Vertices.Add(P.X + 00, P.Y - CQ, P.Z + CH);
                                M.Vertices.Add(P.X + 00, P.Y - CH, P.Z + CQ);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 6);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 1);
                                //******************
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y + CQ, P.Z - CQ);
                                M.Vertices.Add(P.X + CH, P.Y + CQ, P.Z - 00);
                                M.Vertices.Add(P.X + CQ, P.Y + CH, P.Z - 00);
                                M.Vertices.Add(P.X + 00, P.Y + CH, P.Z - CQ);
                                M.Vertices.Add(P.X + 00, P.Y + CQ, P.Z - CH);
                                M.Vertices.Add(P.X + CQ, P.Y + 00, P.Z - CH);
                                M.Vertices.Add(P.X + CH, P.Y + 00, P.Z - CQ);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 1);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 6);

                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X - CQ, P.Y + CQ, P.Z - CQ);
                                M.Vertices.Add(P.X - CQ, P.Y + CH, P.Z - 00);
                                M.Vertices.Add(P.X - CH, P.Y + CQ, P.Z - 00);
                                M.Vertices.Add(P.X - CH, P.Y + 00, P.Z - CQ);
                                M.Vertices.Add(P.X - CQ, P.Y + 00, P.Z - CH);
                                M.Vertices.Add(P.X + 00, P.Y + CQ, P.Z - CH);
                                M.Vertices.Add(P.X + 00, P.Y + CH, P.Z - CQ);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 1);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 6);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X - CQ, P.Y - CQ, P.Z - CQ);
                                M.Vertices.Add(P.X - CH, P.Y - CQ, P.Z - 00);
                                M.Vertices.Add(P.X - CQ, P.Y - CH, P.Z - 00);
                                M.Vertices.Add(P.X + 00, P.Y - CH, P.Z - CQ);
                                M.Vertices.Add(P.X + 00, P.Y - CQ, P.Z - CH);
                                M.Vertices.Add(P.X - CQ, P.Y + 00, P.Z - CH);
                                M.Vertices.Add(P.X - CH, P.Y + 00, P.Z - CQ);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 1);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 6);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y - CQ, P.Z - CQ);
                                M.Vertices.Add(P.X + CQ, P.Y - CH, P.Z - 00);
                                M.Vertices.Add(P.X + CH, P.Y - CQ, P.Z - 00);
                                M.Vertices.Add(P.X + CH, P.Y + 00, P.Z - CQ);
                                M.Vertices.Add(P.X + CQ, P.Y + 00, P.Z - CH);
                                M.Vertices.Add(P.X + 00, P.Y - CQ, P.Z - CH);
                                M.Vertices.Add(P.X + 00, P.Y - CH, P.Z - CQ);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 1);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 6);

                                #endregion

                                if (g == 0) mesh1.Add(M);
                                if (g == 1) mesh2.Add(M);
                            }
                        }
                    }
                }
            }
        }
        public class tribox
        {
            #region table
            readonly int[,] triTable = {
      {-1,-1,-1,-1},
      {0,2,1,-1},
      {3,5,0,-1},
      {1,3,5,2},
      {1,4,3,-1},
      {0,2,4,3},
      {4,5,0,1},
      {2,4,5,-1},
      //***********************
      {5,4,2,-1},
      {1,0,5,4},
      {3,4,2,0},
      {3,4,1,-1},
      {2,5,3,1},
      {0,5,3,-1},
      {1,2,0,-1},
      {-1,-1,-1,-1} };
            #endregion
            public edge[] Cedges;//= new edge[6];
            public point[] Points;//= new Point3d[4];
            public tribox(point[] pts, edge[] els)
            {
                Cedges = new edge[6];
                Points = new point[4];
                Cedges[0] = els[0]; Cedges[1] = els[1]; Cedges[2] = els[2];
                Cedges[3] = els[3]; Cedges[4] = els[4]; Cedges[5] = els[5];
                Points[0] = pts[0]; Points[1] = pts[1]; Points[2] = pts[2]; Points[3] = pts[3];
            }
            public tribox(point[] pts)
            {
                Cedges = new edge[6];
                Points = new point[4];
                Points[0] = pts[0]; Points[1] = pts[1]; Points[2] = pts[2]; Points[3] = pts[3];
                Cedges[0] = new edge(pts[0], pts[1]);
                Cedges[1] = new edge(pts[0], pts[2]);
                Cedges[2] = new edge(pts[0], pts[3]);
                Cedges[3] = new edge(pts[1], pts[2]);
                Cedges[4] = new edge(pts[2], pts[3]);
                Cedges[5] = new edge(pts[3], pts[1]);
            }
            public tribox(List<Point3d> pts)
            {
                Cedges = new edge[6];
                Points = new point[4];
                Points[0] = new point(pts[0]); Points[1] = new point(pts[1]);
                Points[2] = new point(pts[2]); Points[3] = new point(pts[3]);
                Cedges[0] = new edge(Points[0], Points[1]);
                Cedges[1] = new edge(Points[0], Points[2]);
                Cedges[2] = new edge(Points[0], Points[3]);
                Cedges[3] = new edge(Points[1], Points[2]);
                Cedges[4] = new edge(Points[2], Points[3]);
                Cedges[5] = new edge(Points[3], Points[1]);
            }
            public tribox(Point3f[] pts)
            {
                Cedges = new edge[6];
                Points = new point[4];
                Points[0] = new point(pts[0]); Points[1] = new point(pts[1]);
                Points[2] = new point(pts[2]); Points[3] = new point(pts[3]);
                Cedges[0] = new edge(Points[0], Points[1]);
                Cedges[1] = new edge(Points[0], Points[2]);
                Cedges[2] = new edge(Points[0], Points[3]);
                Cedges[3] = new edge(Points[1], Points[2]);
                Cedges[4] = new edge(Points[2], Points[3]);
                Cedges[5] = new edge(Points[3], Points[1]);
            }
            public void getValue(double val1, double val2, double val3, double val4)
            {
                this.Points[0].Value = val1;
                this.Points[1].Value = val2;
                this.Points[2].Value = val3;
                this.Points[3].Value = val4;

            }
            public bool Compute(out Mesh mesh, double isolevel)
            {
                for (int i = 0; i < this.Cedges.Length; i++)
                {
                    this.Cedges[i].VertexInterp(isolevel);
                }
                mesh = new Mesh();
                int sign = 0;
                if (Points[0].Value > isolevel) sign |= 1;
                if (Points[1].Value > isolevel) sign |= 2;
                if (Points[2].Value > isolevel) sign |= 4;
                if (Points[3].Value > isolevel) sign |= 8;
                if (sign == 0 || sign == 15) { return false; }
                else
                {
                    int a = triTable[sign, 0];
                    int b = triTable[sign, 1];
                    int c = triTable[sign, 2];
                    int d = triTable[sign, 3];
                    if (d == -1)
                    {
                        mesh.Vertices.Add(this.Cedges[a].Cpoint);
                        mesh.Vertices.Add(this.Cedges[b].Cpoint);
                        mesh.Vertices.Add(this.Cedges[c].Cpoint);
                        mesh.Faces.AddFace(0, 1, 2);
                    }
                    else
                    {
                        mesh.Vertices.Add(this.Cedges[a].Cpoint);
                        mesh.Vertices.Add(this.Cedges[b].Cpoint);
                        mesh.Vertices.Add(this.Cedges[c].Cpoint);
                        mesh.Vertices.Add(this.Cedges[d].Cpoint);
                        mesh.Faces.AddFace(0, 1, 2, 3);
                    }
                    return true;
                }
            }
            public static List<tribox> FromMesh(Mesh mesh)
            {
                List<tribox> output = new List<tribox>();
                if (mesh.Vertices.Count == 0) return output;
                Point3d cen = new Point3d();
                for (int i = 0; i < mesh.Vertices.Count; i++)
                {
                    cen += mesh.Vertices[i];
                }
                cen /= mesh.Vertices.Count;
                for (int i = 0; i < mesh.Faces.Count; i++)
                {
                    Point3f[] pts = new Point3f[4];
                    pts[0] = new Point3f((float)cen.X, (float)cen.Y, (float)cen.Z);
                    pts[1] = mesh.Vertices[mesh.Faces[i].A];
                    pts[2] = mesh.Vertices[mesh.Faces[i].B];
                    pts[3] = mesh.Vertices[mesh.Faces[i].C];
                    output.Add(new tribox(pts));
                }
                return output;
            }
        }
        public class point
        {
            public Point3f loc = new Point3f();
            public point(Point3d pt)
            {
                this.loc = new Point3f((float)pt.X, (float)pt.Y, (float)pt.Z);
            }
            public point(Point3f pt)
            {
                this.loc = pt;
            }
            public Point3d DrawPoint()
            {
                return new Point3d(this.loc);
            }
            public double Value = 0;
            public bool EqualTo(point pt)
            {
                return this.loc.Equals(pt.loc);
            }
        }
        public class edge
        {
            public Point3f Cpoint = new Point3f();
            public point From;
            public point To;
            public edge(int a, int b, List<point> pts)
            {
                From = pts[a]; To = pts[b];
            }
            public edge(point p1, point p2)
            {
                From = p1; To = p2;
            }
            public void VertexInterp(double isolevel)
            {
                if ((From.Value > isolevel && To.Value < isolevel) || (From.Value < isolevel && To.Value > isolevel))
                    VertexInterp(isolevel, From.loc, To.loc, From.Value, To.Value);
            }
            private void VertexInterp(double isolevel, Point3f p1, Point3f p2, double valp1, double valp2)
            {
                double mu;
                if (Math.Abs(isolevel - valp1) < 0.00001) Cpoint = p1;
                if (Math.Abs(isolevel - valp2) < 0.00001) Cpoint = p2;
                if (Math.Abs(valp1 - valp2) < 0.00001) Cpoint = p1;
                Point3f p = new Point3f();
                mu = (isolevel - valp1) / (valp2 - valp1);
                p.X = (float)(p1.X + (p2.X - p1.X) * mu);
                p.Y = (float)(p1.Y + (p2.Y - p1.Y) * mu);
                p.Z = (float)(p1.Z + (p2.Z - p1.Z) * mu);
                Cpoint = p;
            }
            private void VertexInterp(double isolevel, Point3d p1, Point3d p2, double valp1, double valp2)
            {
                double mu;
                if (Math.Abs(isolevel - valp1) < 0.00001)
                    Cpoint = new Point3f((float)p1.X, (float)p1.Y, (float)p1.Z);
                if (Math.Abs(isolevel - valp2) < 0.00001)
                    Cpoint = new Point3f((float)p2.X, (float)p2.Y, (float)p2.Z);
                if (Math.Abs(valp1 - valp2) < 0.00001)
                    Cpoint = new Point3f((float)p1.X, (float)p1.Y, (float)p1.Z);
                Point3f p = new Point3f();
                mu = (isolevel - valp1) / (valp2 - valp1);
                p.X = (float)(p1.X + (p2.X - p1.X) * mu);
                p.Y = (float)(p1.Y + (p2.Y - p1.Y) * mu);
                p.Z = (float)(p1.Z + (p2.Z - p1.Z) * mu);
                Cpoint = p;
            }
            public Line DrawLine(List<Point3f> Points)
            {
                Point3d p1 = new Point3d(From.loc);
                Point3d p2 = new Point3d(To.loc);
                return new Line(p1, p2);
            }
            public Line DrawLine(List<point> Points)
            {
                Point3d p1 = new Point3d(From.loc);
                Point3d p2 = new Point3d(To.loc);
                return new Line(p1, p2);
            }
            public bool EqualTo(edge el)
            {
                if (this.From.Equals(el.From) && this.To.Equals(el.To)) return true;
                if (this.From.Equals(el.To) && this.To.Equals(el.From)) return true;
                return false;
            }
            public bool IsValid()
            {
                return this.From.EqualTo(this.To);
            }
        }
        public Meta_Triangle() { }
        public Mesh Compute(double isolevel, int length, int width, int height, double cellsize)
        {
            Mesh mesh = new Mesh();
            Matrix_Octahedra mo = new Matrix_Octahedra(length, width, height, cellsize);
            List<Mesh> x = new List<Mesh>();
            x.AddRange(mo.mesh1); x.AddRange(mo.mesh2);
            List<tribox> boxes = new List<tribox>();
            for (int i = 0; i < x.Count; i++)
            {
                boxes.AddRange(tribox.FromMesh(x[i]));
            }
            //**************************************************
            Mesh[] meshes = new Mesh[boxes.Count];
            bool[] signs = new bool[boxes.Count];
            //for(int i = 0;i < boxes.Count;i++){
            System.Threading.Tasks.Parallel.For(0, boxes.Count, i =>
            {

                for (int j = 0; j < boxes[i].Points.Length; j++)
                {
                    Point3f P = boxes[i].Points[j].loc;
                    boxes[i].Points[j].Value =
                      (Math.Sin(P.X) * Math.Cos(P.Y) +
                      Math.Sin(P.Y) * Math.Cos(P.Z) +
                      Math.Sin(P.Z) * Math.Cos(P.X)) + 0.8;
                }
                Mesh mesh1;
                signs[i] = boxes[i].Compute(out mesh1, isolevel);
                meshes[i] = mesh1;

            });

            //**************************************************

            for (int i = 0; i < boxes.Count; i++)
            {
                if (signs[i]) mesh.Append(meshes[i]);
            }

            return mesh;
        }
    }
    public class MetaMCT
    {
        public MetaMCT(List<Point3d> c, List<double> rad, double Res)
        {
            this.use_charges = c;
            this.use_radii = rad;
            this.res = Res;
        }
        public MetaMCT(List<Point3d> c)
        {
            this.use_charges = c;
            use_radii = new List<double>();
            for (int i = 0; i < c.Count; i++)
            {
                use_radii.Add(100);
            }
        }
        private List<Point3d> use_charges;
        private List<double> use_radii;
        private double cA = -0.444444;
        private double cB = 1.888889;
        private double cC = -2.444444;
        public double res = 50;
        public double use_iso = 0.1;
        public Mesh Metaball(double iso)
        {
            this.use_iso = iso;
            Plane plane = new Plane();
            Plane.FitPlaneToPoints(use_charges, out plane);
            plane.Origin = use_charges[0];
            Interval xSize = new Interval(-use_radii[0], use_radii[0]);
            Box box = new Box(plane, xSize, xSize, xSize);
            int num27 = use_charges.Count - 1;
            for (int i = 1; i <= num27; i++)
            {
                plane.Origin = use_charges[i];
                box.Union(plane.PointAt(-use_radii[i], -use_radii[i], -use_radii[i]));
                box.Union(plane.PointAt(use_radii[i], use_radii[i], use_radii[i]));
            }
            box.Inflate(res);
            int xSet = (int)Math.Round((double)(box.X.Length / res));
            int ySet = (int)Math.Round((double)(box.Y.Length / res));
            int zSet = (int)Math.Round((double)(box.Z.Length / res));
            double xLength = box.X.Length / ((double)xSet);
            double yLength = box.Y.Length / ((double)ySet);
            double zLength = box.Z.Length / ((double)zSet);
            double xBase = xLength / 2.0;
            double yBase = yLength / 2.0;
            double zBase = zLength / 2.0;
            plane.Origin = box.GetCorners()[0];
            List<Point3d> list = new List<Point3d>();
            Mesh mesh = new Mesh();

            for (int j = 0; j <= xSet - 1; j++)
            {
                for (int m = 0; m <= ySet - 1; m++)
                {
                    for (int n = 0; n <= zSet - 1; n++)
                    {
                        Point3d item = plane.PointAt(xBase + (xLength * j), yBase + (yLength * m), zBase + (zLength * n));
                        int xSet1 = use_charges.Count - 1;
                        for (int k = 0; k <= xSet1; k++)
                        {
                            if (item.DistanceTo(use_charges[k]) < (use_radii[k] + res))
                            {
                                Plane pl = plane;
                                pl.Origin = item;
                                Mesh other = this.local_tet(pl, xBase, yBase, zBase);
                                mesh.Append(other);
                                list.Add(item);
                                break;
                            }
                        }
                    }
                }
            }
            mesh.Vertices.CombineIdentical(true, true);
            mesh.UnifyNormals();
            return mesh;
        }
        public Mesh[] Smooth(Mesh mesh, int smooth)
        {
            ////////////////////////////
            Mesh[] msh_smooth = mesh.SplitDisjointPieces();

            for (int k = 0; k < msh_smooth.Length; k++)
            {
                msh_smooth[k].Vertices.CombineIdentical(true, true);
                msh_smooth[k].FaceNormals.ComputeFaceNormals();
                msh_smooth[k].Normals.ComputeNormals();
                msh_smooth[k].UnifyNormals();
                for (int j = 0; j < smooth; j++)
                {
                    #region LaplaceFunction
                    Point3f[] new_v = msh_smooth[k].Vertices.ToArray<Point3f>();

                    for (int i = 0; i < msh_smooth[k].TopologyVertices.Count; i++)
                    {
                        int[] mv = msh_smooth[k].TopologyVertices.MeshVertexIndices(i);
                        int[] ctv = msh_smooth[k].TopologyVertices.ConnectedTopologyVertices(i);
                        Single new_x = 0;
                        Single new_y = 0;
                        Single new_z = 0;
                        foreach (int ctv_i in ctv)
                        {
                            int[] vtc = msh_smooth[k].TopologyVertices.MeshVertexIndices(ctv_i);

                            new_x += msh_smooth[k].Vertices[vtc[0]].X / (Single)ctv.Length;
                            new_y += msh_smooth[k].Vertices[vtc[0]].Y / (Single)ctv.Length;
                            new_z += msh_smooth[k].Vertices[vtc[0]].Z / (Single)ctv.Length;
                        }
                        new_v[mv[0]] = new Point3f(new_x, new_y, new_z);
                    }
                    for (int i = 0; i < msh_smooth[k].Vertices.Count; i++)
                    {
                        msh_smooth[k].Vertices[i] = new_v[i];
                    }

                    #endregion
                    msh_smooth[k].FaceNormals.ComputeFaceNormals();
                    msh_smooth[k].Normals.ComputeNormals();

                    for (int i = 0; i < new_v.Length; i++)
                    {
                        Point3d pointd3 = new_v[i];
                        Point3d pointd4 = pointd3 + (Vector3d)msh_smooth[k].Normals[i] * (res / 10);
                        double v1 = this.calculate_field(pointd3);
                        double v2 = this.calculate_field(pointd4);
                        Point3d pointd2 = new Point3d();
                        if (v2 != v1)
                        {
                            pointd2 = this.interp_vertex(pointd3, pointd4, v1, v2);
                        }
                        else
                        {
                            pointd2 = pointd3;
                        }
                        if (pointd3.DistanceTo(pointd2) > this.res)
                        {
                            new_v[i] = new Point3f((float)pointd3.X, (float)pointd3.Y, (float)pointd3.Z);
                        }
                        else
                        {
                            new_v[i] = new Point3f((float)pointd2.X, (float)pointd2.Y, (float)pointd2.Z);
                        }
                    }

                    for (int i = 0; i < msh_smooth[k].Vertices.Count; i++)
                    {
                        msh_smooth[k].Vertices[i] = new_v[i];
                    }
                }
            }
            //*/
            return msh_smooth;
        }
        private Point3d interp_vertex(Point3d p1, Point3d p2, double v1, double v2)
        {
            return new Point3d(p1 + ((Point3d)(((this.use_iso - v1) / (v2 - v1)) * (p2 - p1))));
        }
        private Mesh local_tet(Plane pl, double xBase, double yBase, double zBase)
        {
            List<Point3d> list = new List<Point3d> {
          pl.PointAt(-xBase, yBase, -zBase),
          pl.PointAt(xBase, yBase, -zBase),
          pl.PointAt(xBase, -yBase, -zBase),
          pl.PointAt(-xBase, -yBase, -zBase),
          pl.PointAt(-xBase, yBase, zBase),
          pl.PointAt(xBase, yBase, zBase),
          pl.PointAt(xBase, -yBase, zBase),
          pl.PointAt(-xBase, -yBase, zBase)
          };

            List<double> list2 = new List<double>();
            foreach (Point3d pointd in list)
            {
                double item = this.calculate_field(pointd);
                list2.Add(item);
            }
            DataTree<Point3d> tree = new DataTree<Point3d>();
            Point3d[] data = new Point3d[] { list[0], list[2], list[3], list[7] };
            tree.AddRange(data, new GH_Path(0));
            data = new Point3d[] { list[0], list[2], list[6], list[7] };
            tree.AddRange(data, new GH_Path(1));
            data = new Point3d[] { list[0], list[4], list[6], list[7] };
            tree.AddRange(data, new GH_Path(2));
            data = new Point3d[] { list[0], list[6], list[1], list[2] };
            tree.AddRange(data, new GH_Path(3));
            data = new Point3d[] { list[0], list[6], list[1], list[4] };
            tree.AddRange(data, new GH_Path(4));
            data = new Point3d[] { list[5], list[6], list[1], list[4] };
            tree.AddRange(data, new GH_Path(5));
            DataTree<double> tree2 = new DataTree<double>();
            tree2.AddRange(new double[] { list2[0], list2[2], list2[3], list2[7] }, new GH_Path(0));
            tree2.AddRange(new double[] { list2[0], list2[2], list2[6], list2[7] }, new GH_Path(1));
            tree2.AddRange(new double[] { list2[0], list2[4], list2[6], list2[7] }, new GH_Path(2));
            tree2.AddRange(new double[] { list2[0], list2[6], list2[1], list2[2] }, new GH_Path(3));
            tree2.AddRange(new double[] { list2[0], list2[6], list2[1], list2[4] }, new GH_Path(4));
            tree2.AddRange(new double[] { list2[5], list2[6], list2[1], list2[4] }, new GH_Path(5));
            Mesh mesh = new Mesh();
            foreach (GH_Path path in tree2.Paths)
            {
                List<Point3d> list3 = tree.Branch(path);
                List<double> list4 = tree2.Branch(path);
                int num2 = 0;
                if (list4[0] < this.use_iso)
                {
                    num2 |= 1;
                }
                if (list4[1] < this.use_iso)
                {
                    num2 |= 2;
                }
                if (list4[2] < this.use_iso)
                {
                    num2 |= 4;
                }
                if (list4[3] < this.use_iso)
                {
                    num2 |= 8;
                }
                Mesh other = new Mesh();
                switch (num2)
                {
                    case 1:
                    case 14:
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[1], list4[0], list4[1]));
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[2], list4[0], list4[2]));
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[3], list4[0], list4[3]));
                        other.Faces.AddFace(0, 1, 2);
                        break;

                    case 2:
                    case 13:
                        other.Vertices.Add(this.interp_vertex(list3[1], list3[0], list4[1], list4[0]));
                        other.Vertices.Add(this.interp_vertex(list3[1], list3[3], list4[1], list4[3]));
                        other.Vertices.Add(this.interp_vertex(list3[1], list3[2], list4[1], list4[2]));
                        other.Faces.AddFace(0, 1, 2);
                        break;

                    case 3:
                    case 12:
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[3], list4[0], list4[3]));
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[2], list4[0], list4[2]));
                        other.Vertices.Add(this.interp_vertex(list3[1], list3[3], list4[1], list4[3]));
                        other.Faces.AddFace(0, 1, 2);
                        other.Vertices.Add(this.interp_vertex(list3[1], list3[2], list4[1], list4[2]));
                        other.Faces.AddFace(2, 3, 1);
                        break;

                    case 4:
                    case 11:
                        other.Vertices.Add(this.interp_vertex(list3[2], list3[0], list4[2], list4[0]));
                        other.Vertices.Add(this.interp_vertex(list3[2], list3[1], list4[2], list4[1]));
                        other.Vertices.Add(this.interp_vertex(list3[2], list3[3], list4[2], list4[3]));
                        other.Faces.AddFace(0, 1, 2);
                        break;

                    case 5:
                    case 10:
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[1], list4[0], list4[1]));
                        other.Vertices.Add(this.interp_vertex(list3[2], list3[3], list4[2], list4[3]));
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[3], list4[0], list4[3]));
                        other.Faces.AddFace(0, 1, 2);
                        other.Vertices.Add(this.interp_vertex(list3[1], list3[2], list4[1], list4[2]));
                        other.Faces.AddFace(0, 3, 1);
                        break;

                    case 6:
                    case 9:
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[1], list4[0], list4[1]));
                        other.Vertices.Add(this.interp_vertex(list3[1], list3[3], list4[1], list4[3]));
                        other.Vertices.Add(this.interp_vertex(list3[2], list3[3], list4[2], list4[3]));
                        other.Faces.AddFace(0, 1, 2);
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[2], list4[0], list4[2]));
                        other.Faces.AddFace(0, 3, 2);
                        break;

                    case 7:
                    case 8:
                        other.Vertices.Add(this.interp_vertex(list3[3], list3[0], list4[3], list4[0]));
                        other.Vertices.Add(this.interp_vertex(list3[3], list3[2], list4[3], list4[2]));
                        other.Vertices.Add(this.interp_vertex(list3[3], list3[1], list4[3], list4[1]));
                        other.Faces.AddFace(0, 1, 2);
                        break;
                }
                mesh.Append(other);
            }
            mesh.Vertices.CombineIdentical(true, true);
            return mesh;
        }
        private double calculate_field(Point3d test_pt)
        {
            double num2 = 0.0;
            int num5 = this.use_charges.Count - 1;
            for (int i = 0; i <= num5; i++)
            {
                double x = test_pt.DistanceTo(this.use_charges[i]);
                if (x < this.use_radii[i])
                {
                    num2 += ((((this.cA * Math.Pow(x, 3.0)) / Math.Pow(this.use_radii[i], 3.0)) + ((this.cB * Math.Pow(x, 2.0)) / Math.Pow(this.use_radii[i], 2.0))) + ((this.cC * x) / this.use_radii[i])) + 1.0;
                }
            }
            if (num2 == 0.0)
            {
                num2 = -1.0;
            }
            return num2;
        }
    }
    public class ReactionDiffusion
    {
        //use with Meta
        List<Vertice3> vs;
        public ReactionDiffusion(Mesh x, Curve y)
        {
            Vertice3.CreateCollection(x, out vs);
            Random rnd = new Random();
            for (int i = 0; i < vs.Count; i++)
            {
                Point3d P1 = new Point3d(vs[i].pos.X, vs[i].pos.Y, 0);
                Plane p;
                if (y.TryGetPlane(out p))
                {
                    if (y.Contains(P1, p, Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance) == Rhino.Geometry.PointContainment.Inside)
                    {
                        vs[i].U = 0.5 * (rnd.NextDouble() * 2);
                        vs[i].V = 0.25 * (rnd.NextDouble() * 2);
                    }
                }
            }
        }
        public List<double> RunReactionDiffusion(Mesh x, Curve y, bool z, double iso)
        {
            for (int i = 0; i < vs.Count; i++)
            {
                vs[i].ComputeLaplation3(vs);
            }
            for (int i = 0; i < vs.Count; i++)
            {
                vs[i].ComputeUV1();
            }
            List<double> U2 = new List<double>();
            for (int i = 0; i < vs.Count; i++)
            {
                U2.Add(1 - vs[i].U);
            }
            return U2;
            /*
          A = mc.MeshTopoVerticeDisplay(x, U2);
          if(z) B = tmf.followlines3(x,
              MeshClassLibrary.MeshTopoVerticeConvert.Data_TopoVertices2Vertice(x, U2), iso);
        }
             * */
        }
    }
}