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
    public class TreeMap2D
    {
        public TreeMap2D() { }
        public List<string> __out = new List<string>();
        public List<Polyline> ComputeTreeMap(List<double> x, List<double> y)
        {
            double total = 0;
            for (int i = 0; i < x.Count; i++)
            {
                total += x[i];
            }
            List<double> output = new List<double>();
            for (int i = 0; i < x.Count; i++)
            {
                output.Add(x[i] * (y[0] * y[1]) / total);
            }
            TreeMap treemap1 = new TreeMap(y[0], y[1]);
            List<box> boxes = new List<box>();
            for (int i = 0; i < output.Count; i++)
            {
                boxes.Add(new box(output[i]));
            }
            __out.Add(treemap1.__Error);
            return treemap1.Run(boxes);
        }
        public class TreeMap : box
        {
            public TreeMap(double xEdge, double yEdge)
            {
                XEdge = xEdge;
                YEdge = yEdge;
                area = XEdge * YEdge;
                this.corner = new Point3d(xEdge, yEdge, 0);
            }
            public String __Error = "";
            public List<Polyline> Run(List<box> boxes)
            {
                List<Polyline> PL = new List<Polyline>();

                if (boxes.Count < 2) { __Error = "box Count<2"; return null; }
                double TotalAreal = 0;
                for (int i = 0; i < boxes.Count; i++) { TotalAreal += boxes[i].area; }
                if (TotalAreal > this.area * 1.1) { __Error = "area out of range"; return null; }
                boxes.Sort(CompareDinosByLength); boxes.Reverse();
                //test weather there exists any errors
                bool edgeside = true;
                //true==x false==y
                if (this.XEdge <= this.YEdge) { boxes[0].SetXEdge(this.XEdge); }
                else { boxes[0].SetYEdge(this.YEdge); edgeside = false; }
                //----------------------------------------------------
                List<box> Cache = new List<box>(); Cache.Add(boxes[0]);
                double ISO = boxes[0].XEdge + boxes[0].YEdge;


                for (int i = 1; i < boxes.Count; i++)
                {
                    double total_area = 0; double total_length = 0;
                    //------------------------------------------
                    for (int j = 0; j < Cache.Count; j++) { total_area += Cache[j].area; }
                    total_area += boxes[i].area;
                    double edge1;
                    if (edgeside) { edge1 = total_area / this.XEdge; } else { edge1 = total_area / this.YEdge; }
                    for (int j = 0; j < Cache.Count; j++)
                    {
                        if (edgeside) { Cache[j].SetYEdge(edge1); } else { Cache[j].SetXEdge(edge1); }
                        total_length += Cache[j].XEdge + Cache[j].YEdge;
                    }
                    if (edgeside) { boxes[i].SetYEdge(edge1); } else { boxes[i].SetXEdge(edge1); }
                    //-----------------------------
                    if (total_length <= ISO) { Cache.Add(boxes[i]); ISO = total_length + boxes[i].XEdge + boxes[i].YEdge; }
                    else
                    {
                        if (Cache.Count > 0)
                        {
                            for (int j = 0; j < Cache.Count; j++)
                            {
                                Cache[j].Goback();
                            }
                        }


                        for (int j = 0; j < Cache.Count; j++)
                        {
                            __Error += Cache[j].XEdge + "/" + Cache[j].YEdge + " ";
                            Polyline pl = Cache[j].getPolyline(this.corner);
                            PL.Add(pl);
                            if (edgeside) { this.corner = new Point3d(this.corner.X - Cache[j].XEdge, this.corner.Y, this.corner.Z); }
                            else { this.corner = new Point3d(this.corner.X, this.corner.Y - Cache[j].YEdge, this.corner.Z); }
                        }
                        __Error += "\r\n";


                        this.area -= total_area - boxes[i].area;
                        if (edgeside) this.SetXEdge(this.XEdge);
                        else this.SetYEdge(this.YEdge);
                        if (this.XEdge <= this.YEdge) { edgeside = true; boxes[i].SetXEdge(this.XEdge); }
                        else { edgeside = false; boxes[i].SetYEdge(this.YEdge); }
                        Cache.Clear(); Cache.Add(boxes[i]);
                        this.corner = new Point3d(XEdge, YEdge, 0);
                        ISO = boxes[i].XEdge + boxes[i].YEdge;
                    }
                }
                if (Cache.Count > 0)
                {
                    for (int j = 0; j < Cache.Count; j++)
                    {
                        __Error += Cache[j].XEdge + "/" + Cache[j].YEdge + " ";
                        Polyline pl = Cache[j].getPolyline(this.corner);
                        PL.Add(pl);
                        if (edgeside) { this.corner = new Point3d(this.corner.X - Cache[j].XEdge, this.corner.Y, this.corner.Z); }
                        else { this.corner = new Point3d(this.corner.X, this.corner.Y - Cache[j].YEdge, this.corner.Z); }
                    }
                    __Error += "\r\n";

                }
                return PL;
            }
            private static int CompareDinosByLength(box x, box y)
            {
                if (x == null) { if (y == null) { return 0; } else { return -1; } }
                else
                {
                    if (y == null) { return 1; }
                    else
                    {
                        if (x.area > y.area) return 1;
                        if (x.area == y.area) return 0;
                        if (x.area < y.area) return -1;
                        else return 0;
                    }
                }
            }
        }
        public class box
        {
            #region field
            public double XEdgeCache = 0;
            public double YEdgeCache = 0;
            public double XEdge;
            public double YEdge;
            public double area;
            public Point3d corner = new Point3d();
            #endregion
            #region init
            public box()
            {
                XEdge = 0;
                YEdge = 0;
                area = 0;
            }
            public box(double Area)
            {
                XEdge = 0;
                YEdge = 0;
                area = Area;
            }
            public box(double xEdge, double yEdge)
            {
                XEdge = xEdge;
                YEdge = yEdge;
                area = XEdge * YEdge;
            }
            #endregion
            #region method
            public void SetXEdge(double xEdge)
            {
                XEdgeCache = XEdge;
                YEdgeCache = YEdge;
                XEdge = xEdge;
                YEdge = area / XEdge;
            }
            public void SetYEdge(double yEdge)
            {
                XEdgeCache = XEdge;
                YEdgeCache = YEdge;
                YEdge = yEdge;
                XEdge = area / YEdge;
            }
            public void Goback()
            {
                XEdge = XEdgeCache;
                YEdge = YEdgeCache;
            }
            #endregion
            public Polyline getPolyline(Point3d p)
            {
                Polyline L = new Polyline();

                L.Add(p);
                L.Add(new Point3d(p.X - XEdge, p.Y, p.Z));
                L.Add(new Point3d(p.X - XEdge, p.Y - YEdge, p.Z));
                L.Add(new Point3d(p.X, p.Y - YEdge, p.Z));
                L.Add(p);
                return L;
            }
            public Polyline getPolyline()
            {
                Point3d p = this.corner;
                Polyline L = new Polyline();
                L.Add(p);
                L.Add(new Point3d(p.X - XEdge, p.Y, p.Z));
                L.Add(new Point3d(p.X - XEdge, p.Y - YEdge, p.Z));
                L.Add(new Point3d(p.X, p.Y - YEdge, p.Z));
                L.Add(p);
                return L;
            }
        }
    }
}

