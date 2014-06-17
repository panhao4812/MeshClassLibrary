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
   public class HullFrame
    {
       public HullFrame(){}
       public List<Line> ComputeVoronoi3d(List<Line> x, List<Point3d> y)
       {
           box hu = new box(x);
           hull[] hulls =new  hull[y.Count] ;
           /*
                 for (int ii = 0;ii < y.Count;ii++){
                   hull h = new hull(hu, y[ii]);
                   for(int i = 0;i < y.Count;i++){
                     if( i != ii && y[i].DistanceTo(y[ii]) < h.R * 2){
                       Point3d cen = new Point3d(y[ii]);cen += y[i];cen /= 2;
                       Vector3d v = y[ii] - y[i];
                       Plane plane = new Plane(cen, v);
                       h.intersect(plane);}
                   }
                   hulls.Add(h);
                 }
           */
           ///*
         //  System.Threading.Tasks.Parallel.ForEach(y, pt =>
          // {
       System.Threading.Tasks.Parallel.For(0, y.Count, (iii) => 
            {  
                Point3d pt=y[iii];
               hull h = new hull(hu, pt);
               for (int i = 0; i < y.Count; i++)
               {
                   double t = y[i].DistanceTo(pt);
                   if (t > 0.001 && t < h.R * 2)
                   {
                       Point3d cen = new Point3d(pt); cen += y[i]; cen /= 2;
                       Vector3d v = pt - y[i];
                       Plane plane = new Plane(cen, v);
                       h.intersect(plane);
                   }
               }
               hulls[iii]=h;
           });
           //  */
           List<Line> tree = new List<Line>();
           for (int k = 0; k < hulls.Length; k++)
           {
               hull h = hulls[k];
               for (int i = 0; i < h.edges.Count; i++)
               {
                   tree.Add(new Line(h.edges[i].p1.pos, h.edges[i].p2.pos));
               }
           }
           return tree;
       }
    public class vertex 
    {
        public Point3d pos;
        public int condition = -1;
        public double R = 0;
        public vertex() { }
        public vertex(Point3d pt, double radius)
        {
            this.pos=pt; R = radius;
        }
    }
    public class edge
    {
        public int condition = -1;
        public vertex p1;
        public vertex p2;
        public edge() { }
        public edge(vertex P1, vertex P2)
        {
            p1 = P1; p2 = P2;
        }
    }
    public class box
    {
        public List<int> lp1 = new List<int>();
        public List<int> lp2 = new List<int>();
        public List<Point3d> pts = new List<Point3d>();
        public box(List<Line> l)
        {
            pts.Add(new Point3d(l[0].From));
            pts.Add(new Point3d(l[0].To));
            lp1.Add(0); lp2.Add(1);
            for (int i = 1; i < l.Count; i++)
            {
                bool sign1 = false; bool sign2 = false;
                int a = -1; int b = -1;
                for (int j = 0; j < pts.Count; j++)
                {
                    if (pts[j].DistanceTo(l[i].From) < 0.0000001) { sign1 = true; a = j; }
                    if (pts[j].DistanceTo(l[i].To) < 0.000001) { sign2 = true; b = j; }
                    if (sign1 && sign2) break;
                }
                if (sign1 == false) { pts.Add(new Point3d(l[i].From)); a = pts.Count - 1; }
                if (sign2 == false) { pts.Add(new Point3d(l[i].To)); b = pts.Count - 1; }
                lp1.Add(a);
                lp2.Add(b);
            }
        }
    }
    public class hull
    {
        public double R = double.MaxValue;
        public Point3d center;
        public List<vertex> pts;
        public List<edge> edges;
        public hull(box hu, Point3d cen)
        {
            this.center = new Point3d(cen);
            this.pts = new List<vertex>();
            this.edges = new List<edge>();
            for (int i = 0; i < hu.pts.Count; i++)
            {
                this.pts.Add(new vertex(hu.pts[i], this.center.DistanceTo(hu.pts[i])));
            }
            for (int i = 0; i < hu.lp1.Count; i++)
            {
                this.edges.Add(new edge(pts[hu.lp1[i]], pts[hu.lp2[i]]));
            }
        }
        public void intersect(Plane p)
        {
            for (int i = 0; i < pts.Count; i++)
            {
                double db = p.DistanceTo(pts[i].pos);
                if (Math.Abs(db) < RhinoDoc.ActiveDoc.ModelAbsoluteTolerance) { pts[i].condition = 1; }
                else if (db > 0) { pts[i].condition = 2; }
                else if (db < 0) { pts[i].condition = 0; }
            }
            ///////////////////////
            int ii = 0;
            while (ii < edges.Count)
            {
                if (edges[ii].p1.condition == 0 && edges[ii].p2.condition == 0)
                {
                    edges.RemoveAt(ii);
                }
                else if (edges[ii].p1.condition == 1 && edges[ii].p2.condition == 0)
                {
                    edges.RemoveAt(ii);
                }
                else if (edges[ii].p1.condition == 1 && edges[ii].p2.condition == 1)
                {
                    edges.RemoveAt(ii);
                }
                else if (edges[ii].p1.condition == 0 && edges[ii].p2.condition == 1)
                {
                    edges.RemoveAt(ii);
                }
                else if (edges[ii].p1.condition == 0 && edges[ii].p2.condition == 2)
                {
                    double u ; Line line = new Line(edges[ii].p1.pos, edges[ii].p2.pos);
                    Rhino.Geometry.Intersect.Intersection.LinePlane(line, p, out u);
                    pts.Add(new vertex(line.PointAt(u), this.center.DistanceTo(line.PointAt(u))));
                    edges[ii].p1 = pts[pts.Count - 1];
                    ii++;
                }
                else if (edges[ii].p1.condition == 2 && edges[ii].p2.condition == 0)
                {
                    double u; Line line = new Line(edges[ii].p1.pos, edges[ii].p2.pos);
                    Rhino.Geometry.Intersect.Intersection.LinePlane(line, p, out u);
                    pts.Add(new vertex(line.PointAt(u), this.center.DistanceTo(line.PointAt(u))));
                    edges[ii].p2 = pts[pts.Count - 1];
                    ii++;
                }
                else { ii++; }
            }
            clearnull();
            //////////////////////////////////
            Transform w2p = Transform.PlaneToPlane(Plane.WorldXY, p);
            Transform p2w = Transform.PlaneToPlane(p,Plane.WorldXY );    
            Grasshopper.Kernel.Geometry.Node2List ls = new Grasshopper.Kernel.Geometry.Node2List();
            List<int> count = new List<int>();
            for (int i = 0; i < pts.Count; i++)
            {
                if (pts[i].condition == 1 || pts[i].condition == -1)
                {
                    pts[i].pos.Transform(w2p);
                    ls.Append(new Grasshopper.Kernel.Geometry.Node2(pts[i].pos.X, pts[i].pos.Y));
                    pts[i].pos.Transform(p2w);
                    count.Add(i);
                }
            }
            if (count.Count == 2) edges.Add(new edge(pts[count[0]], pts[count[1]]));
            else if (count.Count > 2)
            {
                List<int> count2 = new List<int>();
                Grasshopper.Kernel.Geometry.ConvexHull.Solver.Compute(ls, count2);
                for (int i = 0; i < count2.Count; i++)
                {
                    int c = i + 1; if (c == count2.Count) c = 0;
                    edges.Add(new edge(pts[count[count2[i]]], pts[count[count2[c]]]));
                }
            }
        }
        public void clearnull()
        {
            int i = 0;
            double max = 0;
            while (i < this.pts.Count)
            {
                if (this.pts[i].condition == 0) { this.pts.RemoveAt(i); }
                else
                {
                    if (max < this.pts[i].R) { max = this.pts[i].R; }
                    i++;
                }
            }
            this.R = max;
        }
    }
   public List<Line>  Offset3D( List<Polyline> x, double y){
           List<Line> output = new List<Line>();
           if (x.Count < 4) return output;
      List<Line> lines = breakPoly(x[0]);

      for(int i = 1;i < x.Count;i++){
        List<Line> ls = breakPoly(x[i]);
        //Print(ls.Count.ToString());
        for(int ii = 0;ii < ls.Count;ii++){
          bool sign = true;
          for(int j = 0;j < lines.Count;j++){
            if(isDumpLines(lines[j], ls[ii])){sign = false;break;}
          }
          //Print(sign.ToString());
          if(sign) lines.Add(ls[ii]);
        }
      }
      Point3d cen = new Point3d();
      for(int i = 0;i < lines.Count;i++){
        cen += lines[i].From;cen += lines[i].To;
      }
      // B = lines;
      cen /= 2 * lines.Count;
      HullFrame.box box = new HullFrame.box(lines);
      HullFrame.hull hull = new HullFrame.hull(box, cen);
      for(int i = 0;i < x.Count;i++){
        if(x[i].Count < 3){//Print("00001");
            return output;
        }
        Plane p = new Plane(x[i][0], x[i][1], x[i][2]);
        Vector3d v = cen - p.ClosestPoint(cen);
        v.Unitize(); p = new Plane(x[i][0], v);
        p.Transform(Transform.Translation(v * y));
        hull.intersect(p);
        hull.clearnull();
      }
    
      for(int i = 0;i < hull.edges.Count;i++){
        output.Add(new Line(hull.edges[i].p1.pos, hull.edges[i].p2.pos));
      }
      List<Point3d> pt = new List<Point3d>();
      for(int i = 0;i < hull.pts.Count;i++){
        pt.Add(hull.pts[i].pos);
      }
           return output;
       }
   public List<Line> breakPoly(Polyline pl)
   {
       List<Line> ls = new List<Line>();
       if (pl.Count < 1) return ls;
       for (int i = 1; i < pl.Count; i++)
       {
           ls.Add(new Line(pl[i], pl[i - 1]));
       }
       return ls;
   }
   public bool isDumpLines(Line l1, Line l2)
   {
       if ((l1.From.DistanceTo(l2.From) < RhinoDoc.ActiveDoc.ModelAbsoluteTolerance) && (l1.To.DistanceTo(l2.To) < RhinoDoc.ActiveDoc.ModelAbsoluteTolerance)) return true;
       if ((l1.From.DistanceTo(l2.To) < RhinoDoc.ActiveDoc.ModelAbsoluteTolerance) && (l1.To.DistanceTo(l2.From) < RhinoDoc.ActiveDoc.ModelAbsoluteTolerance)) return true;
       return false;
   }
    }
}
