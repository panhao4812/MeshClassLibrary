using Rhino.Geometry;

using System;
using System.Collections.Generic;
namespace MeshClassLibrary
{

    public class MeshLayOut
    {
        public MeshLayOut() { }
        public int X = 0;
        public int Y = 1;
        public double Dis = 3;
        private MeshConvexHull ch = new MeshConvexHull();
        public Mesh Mesh2DMinimalBox(Mesh mesh)
        {
            List<Point3d> x = new List<Point3d>();
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            for (int i = 0; i < vs.Count; i++)
            {
                x.Add(new Point3d(vs[i]));
            }
            Grasshopper.Kernel.Geometry.Node2List list = new Grasshopper.Kernel.Geometry.Node2List(x);
            Polyline pl = Grasshopper.Kernel.Geometry.ConvexHull.Solver.ComputeHull(list);
            double t = double.MaxValue;
            Transform xform = new Transform();
            for (int i = 0; i < pl.Count - 1; i++)
            {
                Vector3d Xaxis = pl[i + 1] - pl[i];
                Vector3d Yaxis = Vector3d.CrossProduct(Xaxis, Vector3d.ZAxis);
                Plane p = new Plane(pl[i], Xaxis, Yaxis);
                Polyline pl2 = new Polyline(pl);
                pl2.Transform(Transform.PlaneToPlane(p, Plane.WorldXY));
                Rhino.Geometry.BoundingBox box = pl2.BoundingBox;
                double area = (box.Max.X - box.Min.X) * (box.Max.Y - box.Min.Y);
                if (area < t) { t = area; xform = Transform.PlaneToPlane(p, Plane.WorldXY); }
            }
            mesh.Transform(xform);
            return mesh;
        }
        public Polyline MinimalBox(List<Point3d> x)
        {
            Grasshopper.Kernel.Geometry.Node2List list = new Grasshopper.Kernel.Geometry.Node2List(x);
            Polyline pl = Grasshopper.Kernel.Geometry.ConvexHull.Solver.ComputeHull(list);
            // List<Polyline> boxes = new List<Polyline>();
            Polyline output = new Polyline();
            double t = double.MaxValue;
            for (int i = 0; i < pl.Count - 1; i++)
            {
                Vector3d Xaxis = pl[i + 1] - pl[i];
                Vector3d Yaxis = Vector3d.CrossProduct(Xaxis, Vector3d.ZAxis);
                Plane p = new Plane(pl[i], Xaxis, Yaxis);
                Polyline pl2 = new Polyline(pl);
                pl2.Transform(Transform.PlaneToPlane(p, Plane.WorldXY));
                Rhino.Geometry.BoundingBox box = pl2.BoundingBox;
                Polyline pl3 = new Polyline();
                pl3.Add(box.Corner(false, false, false));
                pl3.Add(box.Corner(false, true, false));
                pl3.Add(box.Corner(true, true, false));
                pl3.Add(box.Corner(true, false, false));
                pl3.Add(box.Corner(false, false, false));
                double area = pl3[1].DistanceTo(pl3[0]) * pl3[1].DistanceTo(pl3[2]);
                if (area < t) { t = area; pl3.Transform(Transform.PlaneToPlane(Plane.WorldXY, p)); output = pl3; }
                // boxes.Add(pl3);
            }
            return output;
        }
        public void layout(ref List<Mesh> x, int DirAxis)
        {

            double t = 0;
            for (int i = 0; i < x.Count; i++)
            {

                Mesh mesh = Mesh2DMinimalBox(x[i]); BoundingBox box = mesh.GetBoundingBox(true);
                x[i].Transform(Transform.PlaneToPlane(new Plane(box.Center, Vector3d.ZAxis), Plane.WorldXY));
                double Legnth = box.Max.X - box.Min.X;
                double Width = box.Max.Y - box.Min.Y;

                if (DirAxis == X)
                {

                    if (Legnth > Width)
                    {
                        x[i].Transform(Transform.Rotation(Math.PI / 2, Point3d.Origin));
                        box = mesh.GetBoundingBox(true);
                    }

                    if (i > 0)
                    {

                        Mesh mesh2 = x[i - 1];
                        BoundingBox box2 = mesh2.GetBoundingBox(true);
                        t += box.Corner(false, true, true).DistanceTo(box.Corner(true, true, true)) / 2;
                        t += box2.Corner(false, true, true).DistanceTo(box2.Corner(true, true, true)) / 2;
                        t += Dis;
                    }
                    x[i].Transform(Transform.Translation(new Vector3d(t, 0, 0)));
                }
                else if (DirAxis == Y)
                {
                    if (Legnth < Width)
                    {
                        x[i].Transform(Transform.Rotation(Math.PI / 2, Point3d.Origin));
                        box = mesh.GetBoundingBox(true);
                    }

                    if (i > 0)
                    {

                        Mesh mesh2 = x[i - 1];
                        BoundingBox box2 = mesh2.GetBoundingBox(true);
                        t += box.Corner(true, false, true).DistanceTo(box.Corner(true, true, true)) / 2;
                        t += box2.Corner(true, false, true).DistanceTo(box2.Corner(true, true, true)) / 2;
                        t += Dis;
                    }
                    x[i].Transform(Transform.Translation(new Vector3d(0, t, 0)));
                }
            }
        }
        public Mesh Mesh3DMinimalBox(Mesh mesh)
        {
            List<Point3d> x = new List<Point3d>();
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            for (int i = 0; i < vs.Count; i++)
            {
                x.Add(new Point3d(vs[i]));
            }
            Mesh mesh2 = ch.ConvexHull(x);
            mesh2.Faces.ConvertQuadsToTriangles();
            double t = double.MaxValue;
            Transform xform = new Transform();
            for (int i = 0; i < mesh2.Faces.Count - 1; i++)
            {
                Point3d p1 = new Point3d(mesh2.Vertices[mesh2.Faces[i].A]);
                Point3d p2 = new Point3d(mesh2.Vertices[mesh2.Faces[i].B]);
                Point3d p3 = new Point3d(mesh2.Vertices[mesh2.Faces[i].C]);
                Plane p = new Plane(p1, p2, p3);
                Mesh mesh3 = new Mesh(); mesh3.Append(mesh);
                mesh3.Transform(Transform.PlaneToPlane(p, Plane.WorldXY));
                Rhino.Geometry.BoundingBox box = mesh3.GetBoundingBox(true);
                double area = (box.Max.X - box.Min.X) * (box.Max.Y - box.Min.Y) * (box.Max.Z - box.Min.Z);
                if (area < t) { t = area; xform = Transform.PlaneToPlane(p, Plane.WorldXY); }
            }
            mesh.Transform(xform);
            return mesh;
        }
    }
    public class MeshUnfold
    {
        public MeshUnfold() { }
        public Mesh Unfold(Mesh mesh)
        {
            List<face> Faces = new List<face>();
            List<edge> Edges = new List<edge>();

            mesh.Faces.ConvertQuadsToTriangles();
            mesh.UnifyNormals();
            mesh.Compact();
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            mesh.FaceNormals.ComputeFaceNormals();
            //Print(mesh.FaceNormals.Count.ToString());
            //  Print(mesh.Vertices.Count.ToString());
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                face f1 = new face(
                  new Point3d(mesh.Vertices[mesh.Faces[i].A]),
                  new Point3d(mesh.Vertices[mesh.Faces[i].B]),
                  new Point3d(mesh.Vertices[mesh.Faces[i].C]),
                  mesh.FaceNormals[i]
                  );

                Faces.Add(f1);
            }
            for (int i = 0; i < el.Count; i++)
            {
                int[] faceid = el.GetConnectedFaces(i);

                edge e1 = new edge(vs[el.GetTopologyVertices(i).I], vs[el.GetTopologyVertices(i).J]);
                if (faceid.Length == 1)
                {
                    e1.Faces = new face[1];
                    e1.Faces[0] = Faces[faceid[0]];
                }
                else if (faceid.Length > 1)
                {
                    e1.Faces = new face[2];
                    e1.Faces[0] = Faces[faceid[0]];
                    e1.Faces[1] = Faces[faceid[1]];
                }
                e1.ID = i;
                Edges.Add(e1);
            }
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                int[] edgeid = el.GetEdgesForFace(i);
                face f1 = Faces[i];
                f1.edges[0] = Edges[edgeid[0]];
                f1.edges[1] = Edges[edgeid[1]];
                f1.edges[2] = Edges[edgeid[2]];
                f1.ID = i;
            }


            /* List<  Mesh> output2 = new List<  Mesh>();
             for (int i = 0; i < Faces.Count; i++)
             {
               output2.Add(Faces[i].DrawFace());
             }
             B = output2;
         */
            face f = Faces[0];
            f.AddTransform(Transform.PlaneToPlane(new Plane(f.Center(), -f.Normal), Plane.WorldXY));

            FaceLoop(f);
            //Print(f.Pts[0].X.ToString() + "/" + f.Pts[0].Y.ToString() + "/" + f.Pts[0].Z.ToString());
            Mesh output = new Mesh();
            for (int i = 0; i < Faces.Count; i++)
            {
                output.Append(Faces[i].DrawFace());
            }
            return output;
        }
        public void FaceLoop(face f)
        {
            if (f.ID != -1)
            {
                ///////////
                // Print("00000002");
                f.ID = -1; f.Fold();
                for (int i = 0; i < f.edges.Length; i++)
                {
                    edge e = f.edges[i];
                    if (e.ID != -1)
                    {
                        e.ID = -1;
                        if (e.Faces.Length == 2)
                        {
                            face f2;
                            if (f.EqualTo(e.Faces[0])) { f2 = e.Faces[1]; } else { f2 = e.Faces[0]; }
                            if (f2.ID != -1)
                            {
                                //Print(f2.ID.ToString());
                                f2.AddTransform(f.Xform);
                                f2.AddTransform(e.fold(f2, e, f));
                                FaceLoop(f2);
                            }
                        }
                    }
                }
                //////
            }
        }
        public class face
        {
            public edge[] edges = new edge[3];
            public Point3d[] Pts = new Point3d[3];
            public int ID = -1;
            public Vector3d Normal;
            public Transform Xform = Rhino.Geometry.Transform.Identity;

            public face(List<Point3d> pts)
            {
                this.Pts = pts.ToArray();
                Normal = new Plane(this.Pts[0], this.Pts[1], this.Pts[2]).Normal;
            }
            public face(Point3d p1, Point3d p2, Point3d p3, Vector3d normal)
            {
                Pts[0] = p1; Pts[1] = p2; Pts[2] = p3;
                Normal = normal;
            }
            public face(Point3d p1, Point3d p2, Point3d p3)
            {
                Pts[0] = p1; Pts[1] = p2; Pts[2] = p3;
                Normal = new Plane(this.Pts[0], this.Pts[1], this.Pts[2]).Normal;
            }
            public face(Point3d[] pts)
            {
                this.Pts = pts;
                Normal = new Plane(this.Pts[0], this.Pts[1], this.Pts[2]).Normal;
            }
            public void Fold()
            {
                transform(this.Xform);
            }
            public void transform(Transform xform)
            {
                for (int i = 0; i < this.Pts.Length; i++)
                {
                    this.Pts[i].Transform(xform);
                }
            }
            public Point3d Center()
            {
                Point3d cen = new Point3d();
                if (this.Pts.Length == 3)
                {
                    cen += this.Pts[0];
                    cen += this.Pts[1];
                    cen += this.Pts[2];
                    cen /= 3;
                }
                return cen;
            }
            public Mesh DrawFace()
            {
                Mesh mesh = new Mesh();
                if (this.Pts.Length == 3)
                {
                    mesh.Vertices.Add(this.Pts[0]);
                    mesh.Vertices.Add(this.Pts[1]);
                    mesh.Vertices.Add(this.Pts[2]);
                    mesh.Faces.AddFace(0, 1, 2);
                    // mesh.Normals.ComputeNormals();
                }
                return mesh;
            }
            public void AddTransform(Transform xform)
            {
                this.Xform = Rhino.Geometry.Transform.Multiply(this.Xform, xform);
            }
            public bool EqualTo(face f1)
            {
                return this.ID == f1.ID;
            }
        }
        public class edge
        {
            public face[] Faces;
            public Point3d From;
            public Point3d To;
            public int ID = -1;

            public edge(Line L)
            {
                From = L.From; To = L.To;
            }
            public edge(Point3d p1, Point3d p2)
            {
                From = p1; To = p2;
            }
            public Line DrawLine()
            {
                return new Line(From, To);
            }
            public bool EqualTo(edge el)
            {
                /*
                 if (this.From.Equals(el.From) && this.To.Equals(el.To)) return true;
                     if (this.From.Equals(el.To) && this.To.Equals(el.From)) return true;
                     return false;
                 */
                return el.ID == this.ID;
            }
            public bool IsValid()
            {
                return this.From.Equals(this.To);
            }
            public Transform fold(face f1, edge e1, face f2)
            {
                Plane p1, p2;
                Point3d cen = this.From;
                Vector3d v = this.To - this.From;
                p1 = new Plane(cen, v, f1.Normal);
                p2 = new Plane(cen, v, f2.Normal);
                return Transform.PlaneToPlane(p1, p2);
            }
        }
        #region analysis
        public double FoldAngle(List<Point3d> pts)
        {
            if (pts.Count == 4) { return FoldAngle(pts[0], pts[1], pts[2], pts[3]); }
            if (pts.Count == 6) { return FoldAngle(pts[0], pts[1], pts[2], pts[3], pts[4], pts[5]); }
            return double.NaN;
        }
        public double FoldAngle(Point3d p1, Point3d p2, Point3d p3, Point3d p4)
        {
            //the triangles are (p1,p2,p3) and (p2,p4,p3)
            Plane plane1 = new Plane(p1, p2, p3);
            Plane plane2 = new Plane(p2, p4, p3);
            double t = plane1.DistanceTo(p4);
            double rad = Vector3d.VectorAngle(plane1.Normal, plane2.Normal);
            if (t < 0) rad *= -1;
            return rad;
        }
        public double FoldAngle(Point3d p1, Point3d p2, Point3d p3, Point3d p4, Point3d p5, Point3d p6)
        {
            Plane plane1 = new Plane(p1, p2, p3);
            Plane plane2 = new Plane(p4, p5, p6);
            Line l;
            if (!Rhino.Geometry.Intersect.Intersection.PlanePlane(plane1, plane2, out l)) return double.NaN;
            double max = l.DistanceTo(p4, false);
            double max2 = l.DistanceTo(p5, false);
            double max3 = l.DistanceTo(p6, false);
            double t = 0;
            if (max >= max2 && max >= max3) t = plane1.DistanceTo(p4);
            if (max2 >= max && max2 >= max3) t = plane1.DistanceTo(p5);
            if (max3 >= max2 && max3 >= max) t = plane1.DistanceTo(p6);
            double rad = Vector3d.VectorAngle(plane1.Normal, plane2.Normal);
            if (t < 0) rad *= -1;
            return rad;
        }
        #endregion    
    }
    public class MainCapsure
    {
        MeshUnfold mu = new MeshUnfold();
        MeshLayOut mlo = new MeshLayOut();
        MeshCreation mc = new MeshCreation();
        string MachineName = "";
        public MainCapsure()
        {
            MachineName = Global.GetComputerName();
        }
        public void Compute(List<Mesh> x, double y, double z,
          out List<Line> lines, out List<Mesh> meshes,
           out List<Rhino.Display.Text3d> t3d1,
            out List<Rhino.Display.Text3d> t3d2,
             out List<Rhino.Display.Text3d> t3d3)
        {

            lines = new List<Line>();
            meshes = new List<Mesh>();
            t3d1 = new List<Rhino.Display.Text3d>();
            t3d2 = new List<Rhino.Display.Text3d>();
            t3d3 = new List<Rhino.Display.Text3d>();
            if (!keys()) return;
            t3d1.AddRange(mc.MeshID(x));
            for (int i = 0; i < x.Count; i++)
            {
                Mesh mesh1 = mu.Unfold(x[i]);
                t3d2.AddRange(mc.FaceID(x[i]));
                meshes.Add(mesh1);
            }
            mlo.Dis = z;//modify the distance beteween the bounding
            mlo.layout(ref meshes, 0);
            for (int i = 0; i < x.Count; i++)
            {
                t3d2.AddRange(mc.FaceID(meshes[i]));

                lines.AddRange(createProfile(meshes[i], y));
            }       
            t3d1.AddRange(mc.MeshID(meshes));
            t3d3 = foldsign(x, meshes);
        }
        public List<Rhino.Display.Text3d> foldsign(List<Mesh> x, List<Mesh> y)
        {
            List<Rhino.Display.Text3d> t3d1 = new List<Rhino.Display.Text3d>();
            List<Point3d> pt = new List<Point3d>();
            List<double> rad = new List<double>();
            for (int i = 0; i < x.Count; i++)
            {
                Mesh mesh2 = y[i];
                mc.MeshClean(ref mesh2, 0.01);
                x[i].UnifyNormals();
                Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh2.TopologyEdges;

                x[i].FaceNormals.ComputeFaceNormals();
                for (int k = 0; k < el.Count; k++)
                {
                    if (el.GetConnectedFaces(k).Length != 2)
                    {
                    }
                    else
                    {
                        MeshFace f1 = x[i].Faces[el.GetConnectedFaces(k)[0]];
                        MeshFace f2 = x[i].Faces[el.GetConnectedFaces(k)[1]];
                        pt.Add(el.EdgeLine(k).PointAt(0.5));
                        double t = mu.FoldAngle(
                          x[i].Vertices[f1.A], x[i].Vertices[f1.B], x[i].Vertices[f1.C],
                          x[i].Vertices[f2.A], x[i].Vertices[f2.B], x[i].Vertices[f2.C]
                          );
                        t *= 180 / Math.PI;

                        rad.Add(t);
                    }
                }
            }
            for (int i = 0; i < pt.Count; i++)
            {
                double t1 = rad[i];
                t1 = Math.Round(t1);
                int t2 = (int)t1;

                if (t2 < -180) t2 = 0;
                if (t2 > 180) t2 = 0;

                t2 *= -1;
                t3d1.Add(new Rhino.Display.Text3d(t2.ToString(), new Plane(pt[i], Vector3d.ZAxis), 1));

            }
            return t3d1;
        }
        public List<Line> createProfile(Mesh x, double y)
        {
            mc.MeshClean(ref x, 0.01);
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = x.TopologyEdges;
            List<Line> ls = new List<Line>();
            for (int i = 0; i < el.Count; i++)
            {
                if (el.GetConnectedFaces(i).Length != 2)
                {
                    ls.Add(el.EdgeLine(i));
                }
                else
                {
                    Line l = el.EdgeLine(i);
                    Vector3d v = l.From - l.To;
                    v.Unitize(); v *= (double)y;
                    Point3d p1 = l.From - v; Point3d p2 = l.To + v;
                    ls.Add(new Line(p1, p2));
                }
            }
            return ls;
        }
        bool keys()
        {
            if (MachineName == "LAPTOP-P70FCQS1") return true;
            else if ( MachineName == "WIN-UVHSGJN4R8N") return true;
            else if (MachineName == "ZZDL-DESKTOP-MINI") return true;
            else if (MachineName == "Account-LWW") return true;
            else if ( MachineName == "USER-20181016WY") return true;
            else if (MachineName == "USER-201810169KD") return true;
            else if (MachineName == "USER-201810169YX") return true;
            else if (MachineName == "USER-20180727ZV") return true;
            else if (MachineName == "DESKTOP-C2JMJB0") return true;
            else return false;
        }
    }
}