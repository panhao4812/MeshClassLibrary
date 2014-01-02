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
}
