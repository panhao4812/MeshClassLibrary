﻿using Rhino;
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
    public class M_Point
    {
        public Point3d pos;
        public List<M_Point> refpoints;
        private double tol = 0.01;
        public Vector3d N;
        public double order = 0;
        public List<M_Point> children = new List<M_Point>();
        /*for sort and lay functions
         in sort method  order means clockwise  angel
         in  lay order equals boolean */
        public M_Point(Point3d Pos)
        {
            refpoints = new List<M_Point>();
            pos = Pos;
        }
        public bool isDump(M_Point other)
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
        private static int CompareDinosByLength(M_Point x, M_Point y)
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
                M_Point p = new M_Point(this.pos);
                p.refpoints.AddRange(this.refpoints);
                p.refpoints.Add(this);
                children.Add(p); children.Add(p);
            }
            for (int i = 0; i < this.refpoints.Count; i++)
            {
                M_Point p = new M_Point(this.pos);
                int after = i + 1; if (after > refpoints.Count - 1) after = 0;
                p.refpoints.Add(refpoints[i]);
                p.refpoints.Add(refpoints[after]);
                children.Add(p);
            }

        }
    }
    public class M_Face
    {
        public Rhino.Geometry.MeshFace face = new MeshFace();
        public List<M_Face> reffaces = new List<M_Face>();
        public double order = 0;
        public M_Face() { }
        public M_Face(MeshFace f)
        {
            this.face = f;
        }
        public M_Face(int a, int b, int c, int d)
        {
            this.face = new MeshFace(a, b, c, d);
        }
        public M_Face(int a, int b, int c)
        {
            this.face = new MeshFace(a, b, c);
        }
        public void Add(M_Face f)
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
        public static List<List<M_Face>> Group(List<M_Face> faces)
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
            List<List<M_Face>> output = new List<List<M_Face>>();
            for (int i = 1; i < level; i++)
            {
                List<M_Face> out1 = new List<M_Face>();
                output.Add(out1);
            }

            for (int i = 0; i < faces.Count; i++)
            {
                output[(int)faces[i].order - 1].Add(faces[i]);
            }
            return output;
        }
        public static bool FindNext(M_Face face, int level)
        {
            if (face.order == 0)
            {
                face.order = level;
                for (int i = 0; i < face.reffaces.Count; i++)
                {
                    FindNext(face.reffaces[i], level);
                }
                return true;
            } return false;
        }



    }

    public class Vertice : BasicVertice
    {
        /// static
        public Vertice(Point3d p):base(p){}
        public static List<Point3d> DisplayPos(List<Vertice> vs)
        {
            List<Point3d> output = new List<Point3d>();
            vs.ForEach(delegate(Vertice v) { output.Add(v.pos); });
            return output;
        }
        public static List<string> Displayenergy(List<Vertice> vs)
        {
            List<string> output = new List<string>();
            vs.ForEach(delegate(Vertice v) { output.Add(v.energy.ToString()); });
            return output;
        }
        public static List<string> DisplayLife(List<Vertice> vs)
        {
            List<string> output = new List<string>();
            vs.ForEach(delegate(Vertice v) { output.Add(v.dead.ToString()); });
            return output;
        }
        public static void CreateCollection(List<Line> x, out List<IndexPair> id, out  List<Vertice> vs)
        {
            id = new List<IndexPair>(); vs = new List<Vertice>();
            id.Add(new IndexPair(0, 1));
            vs.Add(new Vertice(x[0].From, 1));
            vs.Add(new Vertice(x[0].To, 0));
            for (int i = 1; i < x.Count; i++)
            {
                bool sign1 = true;
                bool sign2 = true;
                int a = 0, b = 0;
                for (int j = 0; j < vs.Count; j++)
                {
                    if (vs[j].equalTo(x[i].From)) { sign1 = false; a = j; }
                    if (vs[j].equalTo(x[i].To)) { sign2 = false; b = j; }
                    if (!sign1 && !sign2) { break; }
                }
                if (sign1) { vs.Add(new Vertice(x[i].From)); a = vs.Count - 1; }
                if (sign2) { vs.Add(new Vertice(x[i].To)); b = vs.Count - 1; }
                vs[a].Add(b); vs[b].Add(a);
                id.Add(new IndexPair(a, b));
            }
        }   
        /// ////////////////
        public Vertice(Point3d p, int index):base(p,index) { }
        public List<Polyline> edges = new List<Polyline>();
        public bool transferenergy(double percentage, ref List<Vertice> vs)
        {
            bool sign = false;
            if (!this.dead && this.energy != 0)
            {
                this.dead = true;
                for (int i = 0; i < this.refer.Count; i++)
                {
                    if (vs[this.refer[i]].energy == 0)
                    {
                        vs[this.refer[i]].energy = this.energy * percentage;
                        sign = true;
                    }
                }
            }
            return sign;
        }
        public void CrateEdges(List<Vertice> vs)
        {
            if (this.refer.Count == 3)
            {
                Point3d p1 = vs[this.refer[0]].pos; Vector3d v1 = p1 - this.pos; v1.Unitize();
                Point3d p2 = vs[this.refer[1]].pos; Vector3d v2 = p2 - this.pos; v2.Unitize();
                Point3d p3 = vs[this.refer[2]].pos; Vector3d v3 = p3 - this.pos; v3.Unitize();
                Plane p = new Plane(p1, p2, p3);
                Vector3d n = p.Normal; n.Unitize();
                Point3d N1 = this.pos + n * energy;
                Point3d N2 = this.pos - n * energy;
                Vector3d v0 = v1 + v2; v0.Unitize(); v0 *= this.energy;
                Point3d p12 = this.pos + v0;
                v0 = v2 + v3; v0.Unitize(); v0 *= this.energy;
                Point3d p23 = this.pos + v0;
                v0 = v3 + v1; v0.Unitize(); v0 *= this.energy;
                Point3d p31 = this.pos + v0;
                Polyline pl1 = new Polyline(); pl1.Add(N1); pl1.Add(p12); pl1.Add(N2); pl1.Add(p31);
                Polyline pl2 = new Polyline(); pl2.Add(N1); pl2.Add(p23); pl2.Add(N2); pl2.Add(p12);
                Polyline pl3 = new Polyline(); pl3.Add(N1); pl3.Add(p23); pl3.Add(N2); pl3.Add(p31);
                edges.Add(pl1); edges.Add(pl2); edges.Add(pl3);
            }
        } 
    }
    public class BasicVertice
    {
            /////////////////////basic
    public Point3d pos;
    public bool dead = false;
    public List<int> refer = new List<int>();
    public double energy = 0;
    public BasicVertice(Point3d p)
    {
      pos = new Point3d(p);
    }
    public BasicVertice(Point3d p, int index)
    {
      pos = new Point3d(p);
      this.refer.Add(index);
    }
    public void Add(int i)
    {
      this.refer.Add(i);
    }
    public bool equalTo(Point3d pt)
    {
      if (pos.DistanceTo(pt) < 0.01) { return true; }
      return false;
    }
    /// //////////////////static
    public static void CreateCollection(List<Line> x, out List<IndexPair> id, out  List<BasicVertice> vs)
    {
        id = new List<IndexPair>(); vs = new List<BasicVertice>();
      id.Add(new IndexPair(0, 1));
      vs.Add(new BasicVertice(x[0].From, 1));
      vs.Add(new BasicVertice(x[0].To, 0));
      for (int i = 1; i < x.Count; i++)
      {
        bool sign1 = true;
        bool sign2 = true;
        int a = 0, b = 0;
        for (int j = 0; j < vs.Count; j++)
        {
          if (vs[j].equalTo(x[i].From)) { sign1 = false; a = j; }
          if (vs[j].equalTo(x[i].To)) { sign2 = false; b = j; }
          if (!sign1 && !sign2) { break; }
        }
        if (sign1) { vs.Add(new BasicVertice(x[i].From)); a = vs.Count - 1; }
        if (sign2) { vs.Add(new BasicVertice(x[i].To)); b = vs.Count - 1; }
        vs[a].Add(b); vs[b].Add(a);
        id.Add(new IndexPair(a, b));
      }
    }
    public static List<Point3d> DisplayPos(List<BasicVertice> vs)
    {
      List<Point3d> output = new List<Point3d>();
      vs.ForEach(delegate(BasicVertice v) { output.Add(v.pos); });
      return output;
    }
    public static List<string> Displayenergy(List<BasicVertice> vs)
    {
      List<string> output = new List<string>();
      vs.ForEach(delegate(BasicVertice v) { output.Add(v.energy.ToString()); });
      return output;
    }
    public static List<string> DisplayLife(List<BasicVertice> vs)
    {
      List<string> output = new List<string>();
      vs.ForEach(delegate(BasicVertice v) { output.Add(v.dead.ToString()); });
      return output;
    }
    }
    /// <summary>
    /// Remesh class has bugs which are not fixed.Use Vertice2 instead.
    /// </summary>
    public class Vertice2 : BasicVertice
    {
        /// static
        public Vertice2(Point3d p) : base(p) { }
        public static List<Point3d> DisplayPos(List<Vertice2> vs)
        {
            List<Point3d> output = new List<Point3d>();
            vs.ForEach(delegate(Vertice2 v) { output.Add(v.pos); });
            return output;
        }
        public static List<string> Displayenergy(List<Vertice2> vs)
        {
            List<string> output = new List<string>();
            vs.ForEach(delegate(Vertice2 v) { output.Add(v.energy.ToString()); });
            return output;
        }
        public static List<string> DisplayLife(List<Vertice2> vs)
        {
            List<string> output = new List<string>();
            vs.ForEach(delegate(Vertice2 v) { output.Add(v.dead.ToString()); });
            return output;
        }

        public static List<string> DisplayRef(List<Vertice2> vs)
        {
            List<string> output = new List<string>();
            vs.ForEach(delegate(Vertice2 v)
            {
                string str = "";
                for (int i = 0; i < v.refer.Count; i++)
                {
                    str += v.refer[i].ToString() + "/";
                }
                output.Add(str);
            });
            return output;
        }
        public static void CreateCollection(List<Line> x, out List<IndexPair> id, out  List<Vertice2> vs)
        {
            id = new List<IndexPair>(); vs = new List<Vertice2>();
            id.Add(new IndexPair(0, 1));
            vs.Add(new Vertice2(x[0].From, 1));
            vs.Add(new Vertice2(x[0].To, 0));
            for (int i = 1; i < x.Count; i++)
            {
                bool sign1 = true;
                bool sign2 = true;
                int a = 0, b = 0;
                for (int j = 0; j < vs.Count; j++)
                {
                    if (vs[j].equalTo(x[i].From)) { sign1 = false; a = j; }
                    if (vs[j].equalTo(x[i].To)) { sign2 = false; b = j; }
                    if (!sign1 && !sign2) { break; }
                }
                if (sign1) { vs.Add(new Vertice2(x[i].From)); a = vs.Count - 1; }
                if (sign2) { vs.Add(new Vertice2(x[i].To)); b = vs.Count - 1; }
                vs[a].Add(b); vs[b].Add(a);
                id.Add(new IndexPair(a, b));
            }
        }
        /// ////////////////
        

        public Vertice2(Point3d p, int index) : base(p, index) { }
        Vector3d N = Vector3d.ZAxis;


     public static  List<Polyline> Remesh(List<Vertice2> vs)
        {
            List<Polyline> output = new List<Polyline>();
            List<List<IndexPair>> index = new List<List<IndexPair>>();
            List<List<bool>> sign = new List<List<bool>>();
            int TCount = 0;
            for (int i = 0; i < vs.Count; i++)
            {
                List<IndexPair> children = new List<IndexPair>();
                List<bool> sign2 = new List<bool>();
                for (int j = 0; j < vs[i].refer.Count; j++)
                {
                    TCount++;
                    int after = j + 1;
                    if (after == vs[i].refer.Count) after = 0;
                    children.Add(new IndexPair(vs[i].refer[j], vs[i].refer[after]));

                    //   Print(vs[i].refer[j].ToString() + "~" + i.ToString() + "~" + vs[i].refer[after].ToString());
                    sign2.Add(true);
                }
                index.Add(children);
                sign.Add(sign2);
            }
            ////////////////////////////////////////////////////////
            for (int i = 0; i < vs.Count; i++)
            {
                for (int j = 0; j < vs[i].refer.Count; j++)
                {
                    if (sign[i][j])
                    {
                        sign[i][j] = false;
                        Polyline pl = new Polyline();
                        pl.Add(vs[i].pos);
                        ///to find a start vertice to construct polyline
                        bool signST = true;
                        int before = i;
                        int next = index[before][j].J;
                        //string error = "";
                        for (int loop = 0; loop < TCount; loop++)
                        {
                            signST = false;
                           // error += next.ToString() + "-";
                            for (int k = 0; k < index[next].Count; k++)
                            {
                                if (index[next][k].I == before)
                                {
                                    if (sign[next][k])
                                    {
                                        sign[next][k] = false;
                                        pl.Add(vs[next].pos);
                                        before = next;
                                        next = index[before][k].J;
                                        signST = true;
                                        break;
                                    }
                                }
                            }
                            if (!signST) break;
                        }
                    //    Print(error);
                        /////////////////
                        if (pl.Count > 2) pl.Add(pl[0]);
                        output.Add(pl);
                    }//if
                }//j
            }//i
            return output;
        }

        public void Sort(List<Vertice2> vs)
        { //sort the refer points in clockwise order
            List<IndexPair> refpoints = new List<IndexPair>();
            if (this.refer.Count < 3) return;
            Plane p1 = new Plane(this.pos, this.N);
            Plane p2 = new Plane(new Point3d(0, 0, 0), new Vector3d(0, 0, 1));
            for (int i = 0; i < this.refer.Count; i++)
            {
                Point3d P = new Point3d(vs[this.refer[i]].pos);
                P.Transform(Transform.PlaneToPlane(p1, p2));
                Vector3d v = new Vector3d(P.X, P.Y, 0);
                double t = 0;
                if (P.Y >= 0) { t = Vector3d.VectorAngle(new Vector3d(1, 0, 0), v); }
                else { t = Math.PI * 2 - Vector3d.VectorAngle(new Vector3d(1, 0, 0), v); }
                refpoints.Add(new IndexPair(this.refer[i], (int)(t * 1000000)));
            }
            refpoints.Sort(CompareDinosByLength);
            this.refer.Clear();
            for (int i = 0; i < refpoints.Count; i++)
            {
                this.refer.Add(refpoints[i].I);
            }
        }
        private static int CompareDinosByLength(IndexPair x, IndexPair y)
        {

            if (x.J > y.J) return 1;
            if (x.J == y.J) return 0;
            if (x.J < y.J) return -1;
            else return 0;
        }
    }
    class M_Random
    {
        public int seed = -1;
        public M_Random() { }
        public double NextDouble()
        {
            seed++;
            if (seed == NL.Length) seed = 0;
            return NL[seed];
        }
        double[] NL = {0.84,0.13,0.14,0.23,0.76,0.54,0.67,0.59,0.91,0.74,
0.19,0.71,0.88,0.46,0.78,0.15,0.22,0.53,0.21,0.14,0.56,0.14,0.09,
0.75,0.33,0.67,0.26,0.24,0.22,0.04,0.03,0.96,0.71,0.13,0.43,0.94,
0.11,0.79,0.16,0.42,0.67,0.17,0.11,0.09,0.02,0.32,0.9,0.31,0.57,
0.05,0.94,0.6,0.6,0.82,0.04,0.7,0.05,0.39,0.9,0.09,0.28,0.43,0.37,
0.87,0.71,0.23,0,0.75,0.03,0.85,0.04,0.43,0.37,0.79,0.47,0.39,0.03,
1,0.73,0.02,0.77,0.95,0.67,0.17,0.11,0.43,0.36,0.89,0.09,0.73,0.89,
0.72,0.89,0.07,0.14,0.23,0.81,0.24,0.38,0.79,0.31,0.15,0.28,0.73,
0.01,0.51,0.23,0.81,0.35,0.65,0.67,0.05,0.67,0.88,0.32,0.33,0.77,
0.2,0.77,0.28,0.86,0.12,0.66,0.3,0.96,0.32,0.53,0.3,0.66,0.23,
0.58,0.79,0.62,0.93,0.7,0.62,0.67,0.94,0.16,0.6,0.2,0.56,0.53,
0.43,0.06,0.84,0.06,0.02,0.74,0.8,0.47,0.61,0.47,0.1,0.93,0.2,
0.49,0.99,0.77,0.7,0.97,0.92,0.15,0.12,0.07,0.88,0.43,0.73,0.17,
0.71,0.66,0.83,0.03,0.17,0.08,0.31,0.58,0.23,0.24,0.12,0.26,0.52,
0.56,0.86,0.77,0.97,0.32,0.52,1,0.51,0.12,0.68,0.17,0.46,0.62,
0.28,0.41,0.42,0.36,0.18,0.42,0.32,0.84,0.04,0.14,0.64,0.58,0.3,
0.02,0.63,0.62,0.26,0.74,0.65,0.44,0.46,0.37,0.29,0.35,0.1,0.56,
0.9,0.74,0.67,0.58,0.98,0.66,0.57,0.55,0.8,0.89,0.17,0.87,0.06,
0.7,0.94,0.68,0.52,0.72,0.13,0.4,0.02,0.5,0.37,0.89,0.86,0.94,
0.52,0.03,0.17,0.91,0.12,0.07,0.27,0.62,0.51,0.59,0.17,0.46,
0.16,0.98,0.01,0.75,0.23,0.73,0.45,0.4,0.68,0.95,0.5,0.78,0.85,
0.57,0.22,0.7,0.53,0.4,0.36,0.78,0.72,0.04,0.14,0.54,0.38,0.88,
0.77,0.1,0.6,0.44,0.19,0.35,0.5,0.06,0.56,0.15,0.39,0.27,0.28,
0.64,0.44,0.46,0.26,0.57,0.53,0.39,0.07,0.55,0.85,0.57,0.09,0.11,
0.23,0.39,0.74,0.12,0.85,0.46,0.37,0.35,0.96,0.35,0.8,0.24,0.77,
0.15,0.28,0.79,0.01,0.07,0.3,0.26,0.13,0.72,0.34,0.26,0.79,0.57,
0.02,0.99,0.81,0.22,0.25,0.03,0.35,0.07,0.12,0.11,0.32,0.45,0.3,
0.93,0.9,0.93,0.68,0.09,0.66,0.02,0.8,0.37,0.11,0.28,0.54,0.79,
0.26,0.83,0.99,0.51,0.05,0.48,0.33,0.28,0.45,0.38,0.53,0.74,0.1,
0.77,0.89,0.69,0.03,0.16,0.46,0.57,0.76,0.38,0.36,0.2,0.04,0.25,
0.6,0.77,0.77,0.64,0.88,0.53,0.68,0.46,0.76,0.52,0.08,0.62,0.06,
0.84,0.12,0.02,0.48,0.52,0.4,0.94,0.98,0.89,0.12,0.11,0.35,0.95,
0.82,0.97,0.02,0.89,0.47,0.79,0.47,0.8,0.88,0.56,0.51,0.8,0.5,0,
0.06,0.64,0.01,0.38,0.61,0.41,0.11,0.62,0.45,0.74,0.89,0.84,0.8,
0.1,0.27,0.71,0.65,0.66,0.29,0.93,0.72,0.7,0.44,0.87,0.04,0.29,
0.15,0.25,0.97,0.56,0.51,0.68,0.02,0.4,0.88,0.34,0.88,0.75,0.5,
0.94,0.84,0.2,0.52,0.29,1,0.64,0.99,0.37,0.54,0.16,0.91,0.85,
0.51,0.56,0.28,0.36,0.21,0.13,0.33,0.32,0.26,0.85,0.65,0.89,
0.23,0.21,0.81,0.4,0.22,0.93,0.83,0.9,0.16,0.35,0.09,0.52,0.18,
0.15,0.88,0.41,0.3,0.78,0.72,0.8,0.65,0.66,0.17,0.46,0.11,0.52,
0.13,0.75,0.41,0.18,0.68,0.99,0.55,0.63,0.06,0.78,0.82,0.59,0.15,
0.61,0.33,0.01,0.69,0.16,0.47,0.77,0.18,0.06,0.25,0.93,0.02,0.48,
0.13,0.85,0.25,0.57,0.05,0.35,0.29,0.71,0.79,0.09,0.49,0.98,0.67,
0.1,0.97,0.55,0.09,0.1,0.59,0.71,0.63,0.11,0.48,0.63,0.97,0.01,
0.99,0.35,0.34,0.08,0.49,0.49,0.16,0.2,0.86,0.7,0.38,0.5,0.73,
0.47,0.3,0.44,0.82,0.24,0.52,0.71,0.49,0.37,0.79,0.63,0.97,0.16,
0.34,0.31,0.85,0.03,0.37,0.61,0.6,0.04,0.36,0.94,0.37,0.72,0.6,0,
0.83,0.47,0.24,0.27,0.17,0.59,0.36,0.12,0.41,0.19,0.29,0.23,0.11,
0.26,0.52,0.62,0.55,0.71,0.11,0.33,0.14,0.84,0.35,0.83,0.33,0.77,
0.9,0.69,0.11,0.36,0.07,0.1,0.65,0.52,0.88,0.01,0.12,0.52,0.46,
0.38,0.79,0.22,0.91,0.65,0.74,0.14,0.5,0.34,0.52,0.74,0.39,0.66,
0.61,0.26,0.86,0.98,0.12,0.42,0.95,0.4,0.7,0.67,0.01,0.05,0.12,
0.18,0.59,0.59,0.38,0.51,0.5,0.03,0.25,0.73,0.54,0.92,0.94,0.7,
0.09,0.19,0.27,0.56,0.17,0.37,0.96,0.41,0.49,0.38,0.67,0.9,0.89,
0.71,0.58,0.07,0.68,0.12,0.21,0.86,0.53,0.56,0.55,0.92,0.96,0.01,
0.23,0.36,0.41,0.88,0.72,0.94,0.05,0.42,0.33,0.76,0.13,0.14,0.5,
0.64,0.09,0.71,0.69,0.21,0.92,0.48,0.61,0.79,0.45,0.18,0.05,0.41,
0.71,0.08,0.17,0.53,0.64,0.35,0.6,0.16,0.14,0.6,0.01,0.61,0.66,0.73,
0.86,0.47,0.38,0.82,0.94,0.54,0.63,0.57,0.77,0.82,0.87,0.34,1,0.49,
0.4,0.44,0.91,0.22,0.83,0.31,0.23,0.97,0.25,0.8,0.11,0.78,0.54,0.34,
0.5,0.49,0.7,0.09,0.55,0.19,0.62,0.14,0.42,0.63,0.25,0.51,0.79,0.13,
0.31,0.35,0.66,0.3,0.35,0.11,0.76,0.7,0.69,0.79,0.77,0.35,0.5,0.89,
0.22,0.57,0.71,0.16,1,0.29,0.07,0.28,0.12,0.78,0.8,0.8,0.87,0.26,
0.02,0.28,0.97,0.33,0.52,0.09,0.66,0.91,0.14,0.81,0.42,0.43,0.58,0.8,
0.8,0.91,0.31,0.2,0.69,0.73,0.92,0.84,0.93,0.09,0.51,0.5,0.06,0.03,
0.22,0.88,0.5,0.55,0.25,0.5,0.67,0.42,0.81,0.45,0.83,0.41,0.23,0.31,
0.43,0.89,0.74,0.57,0.71,0.27,0.49,0.21,0.46,0.6,0.12,0.14,0.34,0.18,
0.35,0.89,0.81,0.02,0.03,0.63,0.69,0.26,0.31,0.87,0.19,0.07,0.13,0.38,
0.1,0.87,0.37,0.28,0.5,0.61,0.42,0.03,0.35,0.94,0.78,0.79,0.54,0.01,
0.42,0.9,0.43,0.11,0.17,0.49,0.07,0.93,0.63,0.81,0.38,0.6,0.62,0.17,
0.58,0.86,0.38,0.64,0.15,0.1,0.12,0.6,0.23,0.84,0.64,0.72,0.76,0.32,
0.54,0.87,0.24,0.24,0.09,0.68,0.84,0.41,0.44,0.08,0.91,0.63,0.31,0.17,
0.23,0.56,0.9,0.9,0.99,0.25,0.45,0.48,0.56,0.14,0.64,0.44,0.89,0.82,
0.67,0.59,0.47,0.44,0.73,0.75,0.39,0.76,0.57,0.14,0.51,0.95,0.32,0.17,
0.42,0.3,0.73,0.51,0.79,0.95,0.37,0.68,0.94,0.74,0.73,0.51,0.87,0.06,
0.31,0.1,0.6,0.65,0.78,0.02,0.74};
    }
}
