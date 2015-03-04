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
        public class tribox
    {
        #region table
        readonly int[,] triTable = {
        {-1,-1,-1,-1},
        {0,2,1,-1},
        {0,5,3,-1},
        {1,3,5,2},
        {1,4,3,-1},
        {0,2,4,3},
        {1,0,5,4},
        {5,4,2,-1},
//***********************
        {2,4,5,-1},
        {4,5,0,1},
        {3,4,2,0},
        {3,4,1,-1},
        {2,5,3,1},
        {3,5,0,-1},
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
        public void getvalue(double val1, double val2, double val3, double val4)
        {
            this.Points[0].value = val1;
            this.Points[1].value = val2;
            this.Points[2].value = val3;
            this.Points[3].value = val4;

        }
        public bool Compute(out Mesh mesh, double isolevel)
        {
            for (int i = 0; i < this.Cedges.Length; i++)
            {
                this.Cedges[i].VertexInterp(isolevel);
            }
            mesh = new Mesh();
            int sign = 0;
            if (Points[0].value > isolevel) sign |= 1;
            if (Points[1].value > isolevel) sign |= 2;
            if (Points[2].value > isolevel) sign |= 4;
            if (Points[3].value > isolevel) sign |= 8;
            if (sign == 0 || sign == 15) { return false; }
            else
            {
                int a = triTable[sign, 0];
                int b = triTable[sign, 1];
                int c = triTable[sign, 1];
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
                    mesh.Vertices.Add(this.Cedges[a].Cpoint);
                    mesh.Faces.AddFace(0, 1, 2, 3);
                }
                return true;
            }
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
            public double value = 0;
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
                if ((From.value > isolevel && To.value < isolevel) || (From.value < isolevel && To.value > isolevel))
                    VertexInterp(isolevel, From.loc, To.loc, From.value, To.value);
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
       
    }

