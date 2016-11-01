using Rhino;
using Rhino.Geometry;

using System;
using System.Collections.Generic;

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
}
