using Rhino;
using Rhino.Geometry;

using System;
using System.Collections.Generic;
using System.Drawing;

namespace MeshClassLibrary
{//the Algorithm has bugs when more than 4 Points on a same plane.
 //Use sharpQhull Assemblies instead
    public class MeshConvexHull
    {
        public MeshConvexHull() { }
        private double dotProduct(Vector3d v1, Vector3d v2)
        {
            return (((v2.X * v1.X) + (v2.Y * v1.Y)) + (v2.Z * v1.Z));
        }
        private bool IsPtInMesh(Point3d point, Mesh mesh)
        {
            mesh.FaceNormals.ComputeFaceNormals();
            mesh.UnifyNormals();
            mesh.Normals.ComputeNormals();
            BoundingBox boundingBox = mesh.GetBoundingBox(true);
            if (!boundingBox.Contains(point))
            {
                return false;
            }
            Polyline points = new Polyline();
            points.Add(point);
            Vector3d vectord = new Vector3d(100.0, 100.0, 100.0);
            points.Add(boundingBox.Max + vectord);
            int[] faceIds = null;
            Point3d[] pointdArray = Rhino.Geometry.Intersect.Intersection.MeshPolyline(mesh, new PolylineCurve(points), out faceIds);
            if (pointdArray == null)
            {
                return false;
            }
            if ((pointdArray.Length % 2) == 0)
            {
                return false;
            }
            return true;
        }
        private List<int> seenFaces(Point3d point, Mesh mesh)
        {
            List<int> list2 = new List<int>();
            int num2 = mesh.Faces.Count - 1;
            for (int i = 0; i <= num2; i++)
            {
                Vector3d vectord = mesh.FaceNormals[i];
                Point3d faceCenter = mesh.Faces.GetFaceCenter(i);
                Vector3d vectord2 = this.vec2Pts(point, faceCenter);
                if (this.dotProduct(vectord2, vectord) < 0.001)
                {
                    list2.Add(i);
                }
            }
            return list2;
        }
        private Vector3d vec2Pts(Point3d pt1, Point3d pt2)
        {
            return new Vector3d(pt2.X - pt1.X, pt2.Y - pt1.Y, pt2.Z - pt1.Z);
        }
        private Mesh createNewFaces(Point3d pt, Mesh mesh)
        {
            List<Polyline> list2 = new List<Polyline>();
            list2.AddRange(mesh.GetNakedEdges());
            List<Line> list = new List<Line>();
            int num3 = list2.Count - 1;
            for (int i = 0; i <= num3; i++)
            {
                Polyline polyline = new Polyline();
                list.AddRange(list2[i].GetSegments());
            }
            Mesh other = new Mesh();
            int num4 = list.Count - 1;
            for (int j = 0; j <= num4; j++)
            {
                Mesh mesh4 = new Mesh();
                Line line = list[j];
                mesh4.Vertices.Add(line.From);
                mesh4.Vertices.Add(line.To);
                mesh4.Vertices.Add(pt);
                mesh4.Faces.AddFace(0, 1, 2);
                other.Append(mesh4);
            }
            mesh.Append(other);
            mesh.Vertices.CombineIdentical(true, true);
            mesh.UnifyNormals();
            return mesh;
        }
        private bool checkIfOnLine(Point3d pt1, Point3d pt2, Point3d pttested)
        {
            double num = pttested.DistanceTo(pt1);
            double num2 = pttested.DistanceTo(pt2);
            double num3 = pt1.DistanceTo(pt2);
            return ((num + num2) > num3);
        }
        private bool checkIfOnPlane(Point3d pt1, Point3d pt2, Point3d pt3, Point3d pttested)
        {
            Plane p = new Plane(pt1, pt2, pt3);
            // Print(p.DistanceTo(pttested).ToString());
            return Math.Abs(p.DistanceTo(pttested)) > RhinoDoc.ActiveDoc.ModelAbsoluteTolerance;
        }
        public Mesh ConvexHull(List<Point3d> list)
        {
            Mesh mesh = new Mesh();
            if (list.Count > 3)
            {
                List<Point3d> list2 = new List<Point3d>();
                list2.AddRange(list);
                for (int i = 0; i < list.Count; i++)
                {
                    bool flag = true;
                    int index = i;
                    if (mesh.Vertices.Count == 2)
                    {
                        flag = this.checkIfOnLine(mesh.Vertices[0], mesh.Vertices[1], list2[index]);
                    }
                    if (mesh.Vertices.Count == 3)
                    {
                        flag = checkIfOnPlane(mesh.Vertices[0], mesh.Vertices[1], mesh.Vertices[2], list2[index]);
                    }
                    if (flag)
                    {
                        mesh.Vertices.Add(list2[index]);
                    }

                    if (mesh.Vertices.Count == 4)
                    {
                        mesh.Faces.AddFace(0, 1, 2);
                        mesh.Faces.AddFace(1, 2, 3);
                        mesh.Faces.AddFace(2, 3, 0);
                        mesh.Faces.AddFace(3, 0, 1);
                        mesh.Vertices.CombineIdentical(true, true);
                        mesh.FaceNormals.ComputeFaceNormals();
                        mesh.UnifyNormals();
                        mesh.Normals.ComputeNormals();
                        break;
                    }
                }
                //Print(mesh.Vertices.Count.ToString());
                if (mesh.Vertices.Count < 4) { return default(Mesh); }
                if (list2.Count - 4 <= 0) { return mesh; }
                //
                for (int ii = 0; ii < list2.Count; ii++)
                {
                    int num2 = ii;
                    Point3d point = new Point3d(list2[num2]);
                    if (!this.IsPtInMesh(point, mesh))
                    {
                        List<int> list3 = new List<int>();
                        list3 = this.seenFaces(point, mesh);
                        list3.Sort();
                        for (int i = list3.Count - 1; i >= 0; i += -1)
                        {
                            mesh.Faces.RemoveAt(list3[i]);
                        }
                        mesh = this.createNewFaces(point, mesh);
                    }
                }

            }
            return mesh;
        }
    }
    public class PolylineHull2D
    {
        public PolylineHull2D() { }
        public Polyline Compute(List<Point3d> nodes)
        {
            Polyline pl = new Polyline();
            List<int> hull = new List<int>();
            if (Compute(nodes, ref hull))
            {
                for (int i = 0; i < hull.Count; i++)
                {
                    //Print(hull[i].ToString());
                    pl.Add(nodes[hull[i]]);
                }
            }
            return pl;
        }
        public bool Compute(List<Point3d> nodes, ref List<int> hull)
        {
            if (nodes == null)
            {
                //Print("points");
                nodes = new List<Point3d>();
            }
            if (hull == null)
            {
                // Print("hull");
                hull = new List<int>();
            }
            List<bool> list = new List<bool>();
            hull.Clear();
            list.Clear();
            hull.Capacity = nodes.Count;
            list.Capacity = nodes.Count;
            if (nodes.Count == 0)
            {
                return false;
            }
            if (nodes.Count == 1)
            {
                return false;
            }
            if (nodes.Count == 2)
            {
                hull.Add(0);
                hull.Add(1);
                return true;
            }
            int num9 = nodes.Count - 1;
            for (int i = 0; i <= num9; i++)
            {
                list.Add(false);
            }
            int num = -1;
            int item = -1;
            int num10 = nodes.Count - 1;
            for (int j = 0; j <= num10; j++)
            {
                if (nodes[j] != null)
                {
                    num = j;
                    item = j;
                    break;
                }
            }
            if (num < 0)
            {
                return false;
            }
            //Print(num.ToString() + "*0001");
            int num11 = nodes.Count - 1;
            for (int k = 1; k <= num11; k++)
            {
                if (nodes[k] != null)
                {
                    if (nodes[k].X < nodes[num].X)
                    {
                        num = k;
                    }
                    else if ((nodes[k].X == nodes[num].X) && (nodes[k].Y < nodes[num].Y))
                    {
                        num = k;
                    }
                }
            }
            item = num;
            //  Print(num.ToString() + "!!!!!");
            do
            {
                int num6 = -1;
                int num12 = nodes.Count - 1;
                for (int m = 0; m <= num12; m++)
                {
                    if ((nodes[m] != null) && (!list[m] && (m != item)))
                    {
                        if (num6 == -1)
                        {
                            num6 = m;
                        }
                        else
                        {
                            double num8 = CrossProduct(nodes[m], nodes[item], nodes[num6]);
                            if (num8 == 0.0)
                            {
                                if (DotProduct(nodes[item], nodes[m], nodes[m]) > DotProduct(nodes[item], nodes[num6], nodes[num6]))
                                {
                                    num6 = m;
                                }
                            }
                            else if (num8 < 0.0)
                            {
                                num6 = m;
                            }
                        }
                    }
                }
                item = num6;
                list[item] = true;
                hull.Add(item);
            } while (item != num);
            return true;
        }
        double CrossProduct(Point3d A, Point3d B, Point3d C)
        {
            return (((B.X - A.X) * (C.Y - A.Y)) - ((C.X - A.X) * (B.Y - A.Y)));
        }
        double DotProduct(Point3d A, Point3d B, Point3d C)
        {
            return (((B.X - A.X) * (C.X - A.X)) + ((B.Y - A.Y) * (C.Y - A.Y)));
        }
        List<Point3d> RemoveDupPts(List<Point3d> mypoints, double tolerance)
        {
            List<Point3d> list = new List<Point3d>();
            if (mypoints.Count > 0)
            {
                list.Add(mypoints[0]);
                for (int i = 1; i < mypoints.Count; i++)
                {
                    bool flag = true;
                    Point3d pointd = mypoints[i];
                    for (int j = 0; j < list.Count; j++)
                    {
                        if (OrthoClose(pointd, list[j], tolerance))
                        {
                            flag = false;
                            break;
                        }
                    }
                    if (flag)
                    {
                        list.Add(pointd);
                    }
                }
            }
            return list;
        }
        bool OrthoClose(Point3d Point1, Point3d Point2, double t)
        {
            return (((Math.Abs((double)(Point1.X - Point2.X)) < t) && (Math.Abs((double)(Point1.Y - Point2.Y)) < t)) && (Math.Abs((double)(Point1.Z - Point2.Z)) < t));
        }
    }
    public class HullFrame
    {
        public HullFrame() { }
        public List<Line> ComputeVoronoi3d(List<Line> x, List<Point3d> y)
        {
            box hu = new box(x);
            hull[] hulls = new hull[y.Count];
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
                Point3d pt = y[iii];
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
                hulls[iii] = h;
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
                this.pos = pt; R = radius;
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
                        double u; Line line = new Line(edges[ii].p1.pos, edges[ii].p2.pos);
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
                Transform p2w = Transform.PlaneToPlane(p, Plane.WorldXY);
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
        public List<Line> Offset3D(List<Polyline> x, double y)
        {
            List<Line> output = new List<Line>();
            if (x.Count < 4) return output;
            List<Line> lines = breakPoly(x[0]);

            for (int i = 1; i < x.Count; i++)
            {
                List<Line> ls = breakPoly(x[i]);
                //Print(ls.Count.ToString());
                for (int ii = 0; ii < ls.Count; ii++)
                {
                    bool sign = true;
                    for (int j = 0; j < lines.Count; j++)
                    {
                        if (isDumpLines(lines[j], ls[ii])) { sign = false; break; }
                    }
                    //Print(sign.ToString());
                    if (sign) lines.Add(ls[ii]);
                }
            }
            Point3d cen = new Point3d();
            for (int i = 0; i < lines.Count; i++)
            {
                cen += lines[i].From; cen += lines[i].To;
            }
            // B = lines;
            cen /= 2 * lines.Count;
            HullFrame.box box = new HullFrame.box(lines);
            HullFrame.hull hull = new HullFrame.hull(box, cen);
            for (int i = 0; i < x.Count; i++)
            {
                if (x[i].Count < 3)
                {//Print("00001");
                    return output;
                }
                Plane p = new Plane(x[i][0], x[i][1], x[i][2]);
                Vector3d v = cen - p.ClosestPoint(cen);
                v.Unitize(); p = new Plane(x[i][0], v);
                p.Transform(Transform.Translation(v * y));
                hull.intersect(p);
                hull.clearnull();
            }

            for (int i = 0; i < hull.edges.Count; i++)
            {
                output.Add(new Line(hull.edges[i].p1.pos, hull.edges[i].p2.pos));
            }
            List<Point3d> pt = new List<Point3d>();
            for (int i = 0; i < hull.pts.Count; i++)
            {
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
