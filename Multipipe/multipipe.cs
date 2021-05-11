using System;
using System.Collections.Generic;
using System.Linq;
using KPlankton;
using Rhino.Geometry;
using Rhino;
using Rhino.Collections;
using System.Threading.Tasks;

namespace Multipipe
{
    public class multipipe
    {
        public Mesh Default(List<Line> Lines, double radius)
        {
            List<Line> list5 = KangarooSolver.Util.RemoveDupLn2(Lines, 1E-05);
            List<Point3d> list4 = new List<Point3d>();
            List<double> list3 = new List<double>();
            list3.Add(radius);
            Mesh mesh = this.Fatten(list5, list4, list3, 1, 1, 1, false, 0, 1E-05, 0);
            return mesh;
        }
        public Mesh Default(List<Curve> Curves, double radius)
        {
            List<Line> list5 = CleanInput(Curves, 1E-05);
            List<Point3d> list4 = new List<Point3d>();
            List<double> list3 = new List<double>();
            list3.Add(radius);
            Mesh mesh = this.Fatten(list5, list4, list3, 1, 1, 1, false, 0, 1E-05, 0);
            return mesh;
        }
        public List<double> smoothValues(HalfEdgeGraph graph, List<Point3d> samplePts, List<double> sampleVals)
        {
            List<double> list = new List<double>();
            List<Point3d> positions = graph.GetPositions();
            Point3dList list3 = new Point3dList(samplePts);
            for (int i = 0; i < positions.Count; i++)
            {
                int num3 = list3.ClosestIndex(positions[i]);
                list.Add(sampleVals[num3]);
            }
            bool[] flagArray = new bool[list.Count];
            Point3dList list4 = new Point3dList(positions);
            for (int j = 0; j < samplePts.Count; j++)
            {
                int index = list4.ClosestIndex(samplePts[j]);
                flagArray[index] = true;
            }
            int num = graph.HalfEdges.Count * 2;
            double[] numArray = new double[graph.Nodes.Count];
            for (int k = 0; k < num; k++)
            {
                for (int m = 0; m < graph.Nodes.Count; m++)
                {
                    numArray[m] = 0.0;
                    foreach (Node node in graph.Nodes[m].GetNeighbours())
                    {
                        numArray[m] += list[node.Index];
                    }
                    numArray[m] *= 1.0 / ((double)graph.Nodes[m].GetValence());
                }
                for (int n = 0; n < graph.Nodes.Count; n++)
                {
                    if (!flagArray[n])
                    {
                        list[n] = numArray[n];
                    }
                }
            }
            return list;
        }
        public KPlanktonMesh[] FitCube(Vector3d[] v, Point3d center, double snap, double radius)
        {
            //特殊情况，绘制立方体节点
            Plane plane = new Plane();
            double num = 0.78539816339744828;
            Vector3d vectord = new Vector3d();
            Vector3d vectord2 = new Vector3d();
            int num2 = 0;
            bool[] flagArray = new bool[v.Length];
            int[] numArray1 = new int[v.Length];
            for (int i = 0; i < (v.Length - 1); i++)
            {
                for (int m = i + 1; m < v.Length; m++)
                {
                    double num5 = Vector3d.VectorAngle(v[i], v[m]);
                    if (Math.Abs(1.5707963267948966 - num5) < num)
                    {
                        if (num2 == 0)
                        {
                            vectord = v[i];
                            vectord2 = v[m];
                            num2++;
                        }
                        else
                        {
                            Vector3d vectord3 = vectord;
                            Vector3d vectord4 = vectord2;
                            vectord3.Unitize();
                            vectord4.Unitize();
                            Vector3d vectord5 = Vector3d.CrossProduct(v[i], v[m]);
                            vectord5.Unitize();
                            double num6 = (double)(vectord3 * v[i]);
                            double num7 = (double)(vectord3 * v[m]);
                            double num8 = (double)(vectord3 * vectord5);
                            double num9 = (double)(vectord4 * v[i]);
                            double num10 = (double)(vectord4 * v[m]);
                            double num11 = (double)(vectord4 * vectord5);
                            if ((Math.Abs(num6) > Math.Abs(num7)) && (Math.Abs(num6) > Math.Abs(num8)))
                            {
                                vectord = (num6 > 0.0) ? (vectord + v[i]) : (vectord - v[i]);
                                if (Math.Abs(num10) > Math.Abs(num11))
                                {
                                    vectord2 = (num10 > 0.0) ? (vectord2 + v[m]) : (vectord2 - v[m]);
                                }
                                else
                                {
                                    vectord2 = (num11 > 0.0) ? (vectord2 + vectord5) : (vectord2 - vectord5);
                                }
                            }
                            else if ((Math.Abs(num7) > Math.Abs(num6)) && (Math.Abs(num7) > Math.Abs(num8)))
                            {
                                vectord = (num7 > 0.0) ? (vectord + v[m]) : (vectord - v[m]);
                                if (Math.Abs(num9) > Math.Abs(num11))
                                {
                                    vectord2 = (num9 > 0.0) ? (vectord2 + v[i]) : (vectord2 - v[i]);
                                }
                                else
                                {
                                    vectord2 = (num11 > 0.0) ? (vectord2 + vectord5) : (vectord2 - vectord5);
                                }
                            }
                            else
                            {
                                vectord = (num8 > 0.0) ? (vectord + vectord5) : (vectord - vectord5);
                                if (Math.Abs(num9) > Math.Abs(num10))
                                {
                                    vectord2 = (num9 > 0.0) ? (vectord2 + v[i]) : (vectord2 - v[i]);
                                }
                                else
                                {
                                    vectord2 = (num10 > 0.0) ? (vectord2 + v[m]) : (vectord2 - v[m]);
                                }
                            }
                            num2++;
                            flagArray[i] = true;
                        }
                    }
                }
            }
            if (num2 <= 0)
            {
                return null;
            }
            vectord.Unitize();
            vectord2.Unitize();
            double u = 0.707 * radius;
            plane = new Plane(center, vectord - ((Vector3d)(0.5 * ((vectord * vectord2) * vectord2))), vectord2 - ((Vector3d)(0.5 * ((vectord2 * vectord) * vectord))));
            KPlanktonMesh source = new KPlanktonMesh();
            source.Vertices.Add(plane.PointAt(-u, -u, -u));
            source.Vertices.Add(plane.PointAt(u, -u, -u));
            source.Vertices.Add(plane.PointAt(-u, u, -u));
            source.Vertices.Add(plane.PointAt(u, u, -u));
            source.Vertices.Add(plane.PointAt(-u, -u, u));
            source.Vertices.Add(plane.PointAt(u, -u, u));
            source.Vertices.Add(plane.PointAt(-u, u, u));
            source.Vertices.Add(plane.PointAt(u, u, u));
            Vector3d[] vectordArray = new Vector3d[] { plane.XAxis, plane.YAxis, plane.ZAxis, -1.0 * plane.XAxis, -1.0 * plane.YAxis, -1.0 * plane.ZAxis };
            bool[] flagArray2 = new bool[6];
            int[][] numArray = new int[][] { new int[] { 5, 7, 3, 1 }, new int[] { 7, 6, 2, 3 }, new int[] { 6, 7, 5, 4 }, new int[] { 6, 4, 0, 2 }, new int[] { 4, 5, 1, 0 }, new int[] { 1, 3, 2, 0 } };
            KPlanktonMesh mesh2 = new KPlanktonMesh(source);
            for (int j = 0; j < v.Length; j++)
            {
                double num14 = Math.Cos(snap * 0.39);
                int index = -1;
                for (int n = 0; n < 6; n++)
                {
                    if (!flagArray2[n])
                    {
                        double num17 = (double)(v[j] * vectordArray[n]);
                        if (num17 > num14)
                        {
                            index = n;
                            num14 = num17;
                        }
                    }
                }
                if (index != -1)
                {
                    source.Faces.AddFace(numArray[index]);
                    flagArray2[index] = true;
                }
            }
            if (source.Faces.Count != v.Length)
            {
                return null;
            }
            for (int k = 0; k < 6; k++)
            {
                if (!flagArray2[k])
                {
                    mesh2.Faces.AddFace(numArray[k]);
                }
            }
            if (mesh2.Faces.Count == 0)
            {
                mesh2 = null;
            }
            return new KPlanktonMesh[] { source, mesh2 };
        }    
        public KPlanktonMesh[] TestVoroBall(List<Point3d> points, int sides, double minEdge, Point3d center, double radius, double endRadius)
        {
            //测试用的函数
            KPlanktonMesh[] output = new KPlanktonMesh[3];
            KPlanktonMesh mesh = new KPlanktonMesh();
            if (points.Count == 1)
            {
                Point3d[] pointdArray = points.ToArray();
                Polyline polyline = Polyline.CreateCircumscribedPolygon(new Circle(new Plane(center, (Vector3d)(center - pointdArray[0])), 0.66666666), sides);
                List<int> indices = new List<int>();
                for (int j = 0; j < sides; j++)
                {
                    mesh.Vertices.Add(polyline[j]);
                    indices.Add(j);
                }
                mesh.Faces.AddFace(indices);
            }
            else if (points.Count == 2)
            {
                Plane plane;
                Point3d[] pointdArray2 = points.ToArray();
                Vector3d a = pointdArray2[0] - pointdArray2[1];
                Vector3d b = ((pointdArray2[0] - center) + pointdArray2[1]) - center;
                Vector3d yDirection = Vector3d.CrossProduct(a, b);
                if (yDirection.IsValid && !yDirection.IsTiny(1E-06))
                {
                    plane = new Plane(center, Vector3d.CrossProduct(a, yDirection), yDirection);
                }
                else
                {
                    plane = new Plane(center, pointdArray2[1] - pointdArray2[0]);
                }
                double num3 = Vector3d.VectorAngle(pointdArray2[0] - center, pointdArray2[1] - center);
                mesh.Vertices.Add(plane.PointAt(1.0 / Math.Sin(num3 / 2.0), 0.0));
                mesh.Vertices.Add(plane.PointAt(0.0, 1.0));
                mesh.Vertices.Add(plane.PointAt(-1.0 / Math.Sin(num3 / 2.0), 0.0));
                mesh.Vertices.Add(plane.PointAt(0.0, -1.0));
                mesh.Faces.AddFace(0, 1, 2, 3);
                mesh.Faces.AddFace(0, 3, 2, 1);
            }
            else
            {
                if (Point3d.ArePointsCoplanar(points, 0.001))
                {
                    Point3d[] pointdArray3 = points.ToArray();
                    Plane plane2 = new Plane();
                    Circle circle = new Circle();
                    Circle.TryFitCircleToPoints(pointdArray3, out circle);
                    plane2 = circle.Plane;
                    plane2.Origin = center;
                    double angle = Vector3d.VectorAngle(plane2.XAxis, (Vector3d)(pointdArray3[0] - center), plane2);
                    plane2.Rotate(angle, plane2.ZAxis);
                    double[] source = new double[pointdArray3.Length];
                    int[] items = new int[pointdArray3.Length];
                    double num5 = 0.0;
                    for (int k = 0; k < pointdArray3.Length; k++)
                    {
                        source[k] = Vector3d.VectorAngle(plane2.XAxis, pointdArray3[k] - center, plane2);
                        items[k] = k;
                        num5 += Vector3d.VectorAngle(plane2.ZAxis, pointdArray3[k] - center);
                    }
                    source[0] = 0.0;
                    num5 *= 1.0 / pointdArray3.Length;
                    double w = 1.0 / Math.Sin(num5);
                    Array.Sort(source.ToArray(), pointdArray3);
                    Array.Sort(source, items);
                    mesh.Vertices.Add(plane2.PointAt(0.0, 0.0, w));
                    mesh.Vertices.Add(plane2.PointAt(0.0, 0.0, -w));
                    for (int i = 0; i < pointdArray3.Length; i++)
                    {
                        int index = (i + 1) % pointdArray3.Length;
                        Point3d vertex = pointdArray3[i];
                        double num10 = 0.5 * (source[index] - source[i]);
                        if (index == 0)
                        {
                            num10 = 0.5 * (6.2831853071795862 - source[i]);
                        }
                        vertex.Transform(Transform.Rotation((double)(1.0 * num10), plane2.ZAxis, center));
                        Vector3d vectord4 = vertex - center;
                        vectord4.Unitize();
                        vertex = center + vectord4;
                        double num11 = 0.5 * (Vector3d.VectorAngle(vectord4, (Vector3d)(pointdArray3[index] - center)) + Vector3d.VectorAngle(vectord4, (Vector3d)(pointdArray3[i] - center)));
                        vertex.Transform(Transform.Scale(center, 1.0 / Math.Sin(num11)));
                        mesh.Vertices.Add(vertex);
                    }
                    int[] numArray3 = new int[pointdArray3.Length];
                    for (int i = 0; i < items.Length; i++)
                    {
                        numArray3[i] = Array.IndexOf<int>(items, i);
                    }
                    for (int i = 0; i < items.Length; i++)
                    {
                        int[] numArray1 = new int[4];
                        numArray1[0] = 1;
                        numArray1[1] = (((numArray3[i] - 1) + pointdArray3.Length) % pointdArray3.Length) + 2;
                        numArray1[3] = numArray3[i] + 2;
                        mesh.Faces.AddFace(numArray1);
                    }
                }
                else
                {
                    //多线共节点的正常情况
                    KPlanktonMesh mesh2 = SpherePointsConvexHull(points).ToKPlanktonMesh();
                    //得到convex hull 为mesh2
                    mesh2.GetPositions().ToArray();
                    mesh = CircumDual(mesh2, center);//得到补形mesh
                    output[0] = new KPlanktonMesh(mesh);
                    for (int i = 0; i < mesh2.Faces.Count; i++)
                    {
                        Vector3d vectord5 = mesh.Vertices[i].ToPoint3d() - center;
                        Vector3d vectord6 = mesh2.Vertices[mesh2.Faces.GetFaceVertices(i)[0]].ToPoint3d() - center;
                        vectord5.Unitize();
                        double deg = Vector3d.VectorAngle(vectord5, vectord6);
                        mesh.Vertices.SetVertex(i, center + ((Point3d)(vectord5 * (1.0 / Math.Sin(deg)))));
                    }
                    int count = mesh.Halfedges.Count;
                    for (int i = 0; i < (count - 1); i += 2)
                    {
                        Point3d pointd2 = 0.5 * (mesh.Vertices[mesh.Halfedges[i].StartVertex].ToPoint3d() +
                            mesh.Vertices[mesh.Halfedges[i + 1].StartVertex].ToPoint3d());
                        Vector3d vectord7 = pointd2 - center;
                        Vector3d vectord8 = mesh2.Vertices[mesh.Halfedges[i].AdjacentFace].ToPoint3d() - center;
                        vectord7.Unitize();
                        double deg = Vector3d.VectorAngle(vectord7, vectord8);
                        pointd2 = center + ((Point3d)(vectord7 * (1.0 / Math.Sin(deg))));//新点位置
                        int num20 = mesh.Halfedges.SplitEdge(i);
                        mesh.Vertices.SetVertex(mesh.Halfedges[num20].StartVertex, pointd2);
                    }
                    output[1] = new KPlanktonMesh(mesh);
                    //每条边插入一个点，形成两条边， 然后再从小到大开始删面   
                    count = mesh.Halfedges.Count;
                    for (int i = 0; i < (count - 1); i += 2)
                    {
                        int length = mesh.Faces.GetHalfedges(mesh.Halfedges[i + 1].AdjacentFace).Length;
                        int length2 = mesh.Faces.GetHalfedges(mesh.Halfedges[i].AdjacentFace).Length;
                        if ((length2 > sides) && (length > sides))
                        {
                            int startVertex = mesh.Halfedges[i].StartVertex;
                            int endVertex = mesh.Halfedges[i + 1].StartVertex;
                            Point3d[] pointdArray6 = new Point3d[] { mesh.Vertices[startVertex].ToPoint3d(), mesh.Vertices[endVertex].ToPoint3d() };
                            if (pointdArray6[0].DistanceTo(pointdArray6[1]) < minEdge)
                            {
                                int valence1 = mesh.Vertices.GetValence(startVertex);
                                int valence2 = mesh.Vertices.GetValence(endVertex);
                                if ((valence1 == 2) && (valence2 != 2))
                                {
                                    mesh.Halfedges.CollapseEdge(i + 1);
                                }
                                else if ((valence2 == 2) && (valence1 != 2))
                                {
                                    mesh.Halfedges.CollapseEdge(i);
                                }
                                else
                                {
                                    Point3d pointd3 = 0.5 * (pointdArray6[0] + pointdArray6[1]);
                                    int num27 = mesh.Halfedges.CollapseEdge(i);
                                    mesh.Vertices.SetVertex(mesh.Halfedges[num27].StartVertex, pointd3);
                                }
                            }
                        }
                    }
                }
                mesh.Compact();
                output[2] = new KPlanktonMesh(mesh);
            }
            double num = radius;
            if (points.Count < 3)
            {
                num = endRadius;
            }

            for (int ii = 0; ii < 3; ii++)
            {
                mesh = output[ii];
                for (int i = 0; i < mesh.Vertices.Count; i++)
                {
                    mesh.Vertices.SetVertex(i, center + ((Point3d)((mesh.Vertices[i].ToPoint3d() - center) * num)));
                }
                output[ii] = mesh;
            }
            return output;
        }
        public KPlanktonMesh VoroBall(List<Point3d> points, int sides, double minEdge, Point3d center, double radius, double endRadius)
        {
            //side 多边形截面管子的边数 minEdge 补形最短边
            //该函数做核心节点的，做之前需要处理掉convex hull内部的点
            KPlanktonMesh mesh = new KPlanktonMesh();
            if (points.Count == 1)
            {
                Point3d[] pointdArray = points.ToArray();
                Polyline polyline = Polyline.CreateCircumscribedPolygon(new Circle(new Plane(center, (Vector3d)(center - pointdArray[0])), 0.66666666), sides);
                List<int> indices = new List<int>();
                for (int j = 0; j < sides; j++)
                {
                    mesh.Vertices.Add(polyline[j]);
                    indices.Add(j);
                }
                mesh.Faces.AddFace(indices);
            }
            else if (points.Count == 2)
            {
                Plane plane;
                Point3d[] pointdArray2 = points.ToArray();
                Vector3d a = pointdArray2[0] - pointdArray2[1];
                Vector3d b = ((pointdArray2[0] - center) + pointdArray2[1]) - center;
                Vector3d yDirection = Vector3d.CrossProduct(a, b);
                if (yDirection.IsValid && !yDirection.IsTiny(1E-06))
                {
                    plane = new Plane(center, Vector3d.CrossProduct(a, yDirection), yDirection);
                }
                else
                {
                    plane = new Plane(center, pointdArray2[1] - pointdArray2[0]);
                }
                double num3 = Vector3d.VectorAngle(pointdArray2[0] - center, pointdArray2[1] - center);
                mesh.Vertices.Add(plane.PointAt(1.0 / Math.Sin(num3 / 2.0), 0.0));
                mesh.Vertices.Add(plane.PointAt(0.0, 1.0));
                mesh.Vertices.Add(plane.PointAt(-1.0 / Math.Sin(num3 / 2.0), 0.0));
                mesh.Vertices.Add(plane.PointAt(0.0, -1.0));
                mesh.Faces.AddFace(0, 1, 2, 3);
                mesh.Faces.AddFace(0, 3, 2, 1);
            }
            else
            {
                if (Point3d.ArePointsCoplanar(points, 0.001))
                {
                    Point3d[] pointdArray3 = points.ToArray();
                    Plane plane2 = new Plane();
                    Circle circle = new Circle();
                    Circle.TryFitCircleToPoints(pointdArray3, out circle);
                    plane2 = circle.Plane;
                    plane2.Origin = center;
                    double angle = Vector3d.VectorAngle(plane2.XAxis, (Vector3d)(pointdArray3[0] - center), plane2);
                    plane2.Rotate(angle, plane2.ZAxis);
                    double[] source = new double[pointdArray3.Length];
                    int[] items = new int[pointdArray3.Length];
                    double num5 = 0.0;
                    for (int k = 0; k < pointdArray3.Length; k++)
                    {
                        source[k] = Vector3d.VectorAngle(plane2.XAxis, pointdArray3[k] - center, plane2);
                        items[k] = k;
                        num5 += Vector3d.VectorAngle(plane2.ZAxis, pointdArray3[k] - center);
                    }
                    source[0] = 0.0;
                    num5 *= 1.0 / pointdArray3.Length;
                    double w = 1.0 / Math.Sin(num5);
                    Array.Sort(source.ToArray(), pointdArray3);
                    Array.Sort(source, items);
                    mesh.Vertices.Add(plane2.PointAt(0.0, 0.0, w));
                    mesh.Vertices.Add(plane2.PointAt(0.0, 0.0, -w));
                    for (int i = 0; i < pointdArray3.Length; i++)
                    {
                        int index = (i + 1) % pointdArray3.Length;
                        Point3d vertex = pointdArray3[i];
                        double num10 = 0.5 * (source[index] - source[i]);
                        if (index == 0)
                        {
                            num10 = 0.5 * (6.2831853071795862 - source[i]);
                        }
                        vertex.Transform(Transform.Rotation((double)(1.0 * num10), plane2.ZAxis, center));
                        Vector3d vectord4 = vertex - center;
                        vectord4.Unitize();
                        vertex = center + vectord4;
                        double num11 = 0.5 * (Vector3d.VectorAngle(vectord4, (Vector3d)(pointdArray3[index] - center)) + Vector3d.VectorAngle(vectord4, (Vector3d)(pointdArray3[i] - center)));
                        vertex.Transform(Transform.Scale(center, 1.0 / Math.Sin(num11)));
                        mesh.Vertices.Add(vertex);
                    }
                    int[] numArray3 = new int[pointdArray3.Length];
                    for (int i = 0; i < items.Length; i++)
                    {
                        numArray3[i] = Array.IndexOf<int>(items, i);
                    }
                    for (int i = 0; i < items.Length; i++)
                    {
                        int[] numArray1 = new int[4];
                        numArray1[0] = 1;
                        numArray1[1] = (((numArray3[i] - 1) + pointdArray3.Length) % pointdArray3.Length) + 2;
                        numArray1[3] = numArray3[i] + 2;
                        mesh.Faces.AddFace(numArray1);
                    }
                }
                else
                {
                    //多线共节点的正常情况
                    KPlanktonMesh mesh2 = SpherePointsConvexHull(points).ToKPlanktonMesh();
                    //得到convex hull 为mesh2
                    mesh2.GetPositions().ToArray();
                    mesh = CircumDual(mesh2, center);//得到补形mesh
                    for (int i = 0; i < mesh2.Faces.Count; i++)
                    {
                        Vector3d vectord5 = mesh.Vertices[i].ToPoint3d() - center;
                        Vector3d vectord6 = mesh2.Vertices[mesh2.Faces.GetFaceVertices(i)[0]].ToPoint3d() - center;
                        vectord5.Unitize();
                        double deg = Vector3d.VectorAngle(vectord5, vectord6);
                        mesh.Vertices.SetVertex(i, center + ((Point3d)(vectord5 * (1.0 / Math.Sin(deg)))));
                    }
                    int count = mesh.Halfedges.Count;
                    for (int i = 0; i < (count - 1); i += 2)
                    {
                        Point3d pointd2 = 0.5 * (mesh.Vertices[mesh.Halfedges[i].StartVertex].ToPoint3d() +
                            mesh.Vertices[mesh.Halfedges[i + 1].StartVertex].ToPoint3d());
                        Vector3d vectord7 = pointd2 - center;
                        Vector3d vectord8 = mesh2.Vertices[mesh.Halfedges[i].AdjacentFace].ToPoint3d() - center;
                        vectord7.Unitize();
                        double deg = Vector3d.VectorAngle(vectord7, vectord8);
                        pointd2 = center + ((Point3d)(vectord7 * (1.0 / Math.Sin(deg))));//新点位置
                        int num20 = mesh.Halfedges.SplitEdge(i);
                        mesh.Vertices.SetVertex(mesh.Halfedges[num20].StartVertex, pointd2);
                    }
                    //每条边插入一个点，形成两条边， 然后再从小到大开始删面   
                    count = mesh.Halfedges.Count;
                    for (int i = 0; i < (count - 1); i += 2)
                    {
                        int length = mesh.Faces.GetHalfedges(mesh.Halfedges[i + 1].AdjacentFace).Length;
                        int length2 = mesh.Faces.GetHalfedges(mesh.Halfedges[i].AdjacentFace).Length;
                        if ((length2 > sides) && (length > sides))
                        {
                            int startVertex = mesh.Halfedges[i].StartVertex;
                            int endVertex = mesh.Halfedges[i + 1].StartVertex;
                            Point3d[] pointdArray6 = new Point3d[] { mesh.Vertices[startVertex].ToPoint3d(), mesh.Vertices[endVertex].ToPoint3d() };
                            if (pointdArray6[0].DistanceTo(pointdArray6[1]) < minEdge)
                            {
                                int valence1 = mesh.Vertices.GetValence(startVertex);
                                int valence2 = mesh.Vertices.GetValence(endVertex);
                                if ((valence1 == 2) && (valence2 != 2))
                                {
                                    mesh.Halfedges.CollapseEdge(i + 1);
                                }
                                else if ((valence2 == 2) && (valence1 != 2))
                                {
                                    mesh.Halfedges.CollapseEdge(i);
                                }
                                else
                                {
                                    Point3d pointd3 = 0.5 * (pointdArray6[0] + pointdArray6[1]);
                                    int num27 = mesh.Halfedges.CollapseEdge(i);
                                    mesh.Vertices.SetVertex(mesh.Halfedges[num27].StartVertex, pointd3);
                                }
                            }
                        }
                    }
                }
                mesh.Compact();
            }
            double num = radius;
            if (points.Count < 3)
            {
                num = endRadius;
            }
            for (int i = 0; i < mesh.Vertices.Count; i++)
            {
                mesh.Vertices.SetVertex(i, center + ((Point3d)((mesh.Vertices[i].ToPoint3d() - center) * num)));
            }
            return mesh;
        }
        public KPlanktonMesh CircumDual(KPlanktonMesh P, Point3d center)
        {
            //获得补形
            KPlanktonMesh mesh = new KPlanktonMesh();
            for (int i = 0; i < P.Faces.Count; i++)
            {
                //获取每个面的中心点投影
                int[] faceVertices = P.Faces.GetFaceVertices(i);
                Point3d[] pointdArray = new Point3d[3];
                for (int k = 0; k < 3; k++)
                {
                    pointdArray[k] = P.Vertices[faceVertices[k]].ToPoint3d();
                }
                Vector3d b = pointdArray[0] - pointdArray[1];
                Vector3d vectord2 = Vector3d.CrossProduct(pointdArray[2] - pointdArray[1], b);
                vectord2.Unitize();
                mesh.Vertices.Add(center + vectord2);
            }
            for (int j = 0; j < P.Vertices.Count; j++)
            {
                mesh.Faces.AddFace(P.Vertices.GetVertexFaces(j));
            }
            return mesh;
        }
        public List<Line> CleanInput(List<Curve> list, double KinkAngle)
        {
            List<Line> lines = new List<Line>();
            foreach (Curve curve in list)
            {
                if (curve.IsPolyline())
                {
                    Polyline polyline = new Polyline();
                    if (curve.TryGetPolyline(out polyline))
                    {
                        foreach (Line line in polyline.GetSegments())
                        {
                            lines.Add(line);
                        }
                    }
                }
                else
                {
                    foreach (Line line2 in curve.ToPolyline(0.0, KinkAngle, 0.0, 0.0).ToPolyline().GetSegments())
                    {
                        lines.Add(line2);
                    }
                }
            }
            return KangarooSolver.Util.RemoveDupLn2(lines, 1E-05);
        }
        public Mesh Fatten(List<Line> lines, List<Point3d> nodePts, List<double> radius, double strutSizeMultiplier, double offset, double segmentLength, bool capped, double cubeFit, double tolerance, double KinkAngle)
        {
            HalfEdgeGraph skeleton = new HalfEdgeGraph();
            foreach (Line line in lines)
            {
                skeleton.AddEdge(line, tolerance);
            }
            List<KPlanktonMesh> hubs = new List<KPlanktonMesh>();
            List<KPlanktonMesh> list = new List<KPlanktonMesh>();
            List<double> nodeRadius = new List<double>();
            if (radius.Count == 1)
            {
                for (int m = 0; m < skeleton.Nodes.Count; m++)
                {
                    nodeRadius.Add(radius[0]);
                }
            }
            else
            {
                nodeRadius = smoothValues(skeleton, nodePts, radius);
            }
            for (int i = 0; i < skeleton.Nodes.Count; i++)
            {
                Node node = skeleton.Nodes[i];
                int valence = node.GetValence();
                Vector3d[] v = new Vector3d[valence];
                List<Point3d> neighbourPositions = node.GetNeighbourPositions();
                for (int j = 0; j < valence; j++)
                {
                    v[j] = neighbourPositions[j] - node.Position;
                }
                List<Point3d> points = new List<Point3d>();
                for (int j = 0; j < valence; j++)
                {
                    if (v[j].Length > 1E-05)
                    {
                        v[j].Unitize();
                    }
                    points.Add(node.Position + v[j]);
                }
                double radius = nodeRadius[i] * 1.5;
                KPlanktonMesh[] meshArray = new KPlanktonMesh[2];
                if (((cubeFit > 0.0) && (valence < 7)) && (valence > 1))
                {
                    meshArray = FitCube(v, node.Position, cubeFit, radius);
                }
                if ((meshArray != null) && (meshArray[0] != null))
                {
                    hubs.Add(meshArray[0]);
                    if (meshArray[1] != null)
                    {
                        list.Add(meshArray[1]);
                    }
                }
                else
                {
                    hubs.Add(VoroBall(points, 4, 1.0, node.Position, radius, radius * strutSizeMultiplier));
                }
            }
            Mesh[] struts = new Mesh[skeleton.HalfEdges.Count / 2];
            Parallel.For(0, skeleton.HalfEdges.Count / 2, delegate (int index)
            {
                int num = 2 * index;
                int num2 = skeleton.HalfEdges[num].StartNode.Index;
                int num3 = skeleton.HalfEdges[num + 1].StartNode.Index;
                int indexInStartOutgoing = skeleton.HalfEdges[num].GetIndexInStartOutgoing();
                int f = skeleton.HalfEdges[num + 1].GetIndexInStartOutgoing();
                int[] faceVertices = hubs[num2].Faces.GetFaceVertices(indexInStartOutgoing);
                Point3d[] startLoop = new Point3d[faceVertices.Length];
                for (int i = 0; i < startLoop.Length; i++)
                {
                    startLoop[i] = hubs[num2].Vertices[faceVertices[i]].ToPoint3d();
                }
                int[] numArray2 = hubs[num3].Faces.GetFaceVertices(f);
                Point3d[] endLoop = new Point3d[numArray2.Length];
                for (int j = 0; j < endLoop.Length; j++)
                {
                    endLoop[j] = hubs[num3].Vertices[numArray2[j]].ToPoint3d();
                }
                bool startNaked = skeleton.Nodes[num2].GetValence() == 1;
                bool endNaked = skeleton.Nodes[num3].GetValence() == 1;
                double startRadius = (strutSizeMultiplier * 0.5) * (nodeRadius[num2] + nodeRadius[num3]);
                double startOffset = offset * nodeRadius[num2];
                double endOffset = offset * nodeRadius[num3];
                double num9 = 3.1415926535897931 - KinkAngle;
                if (skeleton.Nodes[num2].GetValence() < 3)
                {
                    List<Vector3d> neighbourVectors = skeleton.Nodes[num2].GetNeighbourVectors(true);
                    bool flag3 = true;
                    if (neighbourVectors.Count == 2)
                    {
                        flag3 = Vector3d.VectorAngle(neighbourVectors[0], neighbourVectors[1]) > num9;
                    }
                    if (flag3)
                    {
                        startOffset = 0.0;
                    }
                }
                if (skeleton.Nodes[num3].GetValence() < 3)
                {
                    List<Vector3d> list2 = skeleton.Nodes[num3].GetNeighbourVectors(true);
                    bool flag4 = true;
                    if (list2.Count == 2)
                    {
                        flag4 = Vector3d.VectorAngle(list2[0], list2[1]) > num9;
                    }
                    if (flag4)
                    {
                        endOffset = 0.0;
                    }
                }
                struts[index] = QuadBridge(startLoop, endLoop, skeleton.HalfEdges[num].GetLine(), startOffset, endOffset, segmentLength, startRadius, startRadius, startNaked, endNaked);
            });
            Mesh mesh = new Mesh();
            foreach (Mesh mesh2 in struts)
            {
                mesh.Append(mesh2);
            }
            foreach (KPlanktonMesh mesh3 in list)
            {
                mesh.Append(mesh3.ToRhinoMesh());
            }
            mesh.Vertices.CullUnused();
            mesh.Vertices.CombineIdentical(true, true);
            mesh.UnifyNormals();
            mesh.RebuildNormals();
            return mesh;
        }
        /// /////////通用函数//////////////////        
        public Mesh SpherePointsConvexHull(List<Point3d> points)
        {
            //绘制凸包
            bool flag;
            Point3d[] vertices = points.ToArray();
            Mesh mesh = new Mesh();
            mesh.Vertices.AddVertices(vertices);
            if (vertices.Length == 3)
            {
                mesh.Faces.AddFace(0, 1, 2);
                return mesh;
            }
            bool[] flagArray = new bool[vertices.Length];
            mesh.Faces.AddFace(0, 1, 2);
            mesh.FaceNormals.ComputeFaceNormals();
            double[] keys = new double[vertices.Length - 3];
            int[] items = new int[vertices.Length - 3];
            for (int i = 0; i < keys.Length; i++)
            {
                Vector3f vectorf = (Vector3f)(vertices[i + 3] - vertices[0]);
                keys[i] = Math.Abs((double)(mesh.FaceNormals[0] * vectorf));
                items[i] = i;
            }
            Array.Sort(keys, items);
            int index = items[items.Length - 1] + 3;
            if ((mesh.FaceNormals[0] * ((Vector3f)(vertices[index] - vertices[0]))) > 0.0)
            {
                mesh.Flip(true, true, true);
            }
            MeshFace face = mesh.Faces[0];
            face = mesh.Faces[0];
            mesh.Faces.AddFace(face[1], face[0], index);
            face = mesh.Faces[0];
            face = mesh.Faces[0];
            mesh.Faces.AddFace(face[2], face[1], index);
            face = mesh.Faces[0];
            face = mesh.Faces[0];
            mesh.Faces.AddFace(face[0], face[2], index);
            flagArray[index] = flag = true;
            flagArray[2] = flag;
            flagArray[0] = flagArray[1] = flag;
            mesh.FaceNormals.ComputeFaceNormals();
            for (int j = 3; j < vertices.Length; j++)
            {
                if (!flagArray[j])
                {
                    bool[] flagArray2 = new bool[mesh.Faces.Count];
                    List<int> list = new List<int>();
                    for (int k = 0; k < mesh.Faces.Count; k++)
                    {
                        face = mesh.Faces[k];
                        flagArray2[k] = (mesh.FaceNormals[k] * ((Vector3f)(vertices[j] - mesh.Vertices[face[0]]))) < 0.0001;
                        if (!flagArray2[k])
                        {
                            list.Add(k);
                        }
                    }
                    List<MeshFace> faces = new List<MeshFace>();
                    for (int m = 0; m < list.Count; m++)
                    {
                        bool[] sameOrientation = new bool[3];
                        int[] edgesForFace = mesh.TopologyEdges.GetEdgesForFace(list[m], out sameOrientation);
                        for (int num6 = 0; num6 < 3; num6++)
                        {
                            int[] connectedFaces = mesh.TopologyEdges.GetConnectedFaces(edgesForFace[num6]);
                            if (flagArray2[connectedFaces[0]] || flagArray2[connectedFaces[1]])
                            {
                                IndexPair topologyVertices = mesh.TopologyEdges.GetTopologyVertices(edgesForFace[num6]);
                                faces.Add(new MeshFace(topologyVertices[sameOrientation[num6] ? 0 : 1], topologyVertices[sameOrientation[num6] ? 1 : 0], j));
                            }
                        }
                    }
                    for (int n = 0; n < mesh.Faces.Count; n++)
                    {
                        if (flagArray2[n])
                        {
                            faces.Add(mesh.Faces[n]);
                        }
                    }
                    mesh.Faces.Clear();
                    mesh.Faces.AddFaces(faces);
                    mesh.FaceNormals.ComputeFaceNormals();
                    flagArray[j] = true;
                }
            }
            mesh.Normals.ComputeNormals();
            return mesh;
        }
        public Mesh joinLoops(Polyline startLoop, Polyline endLoop)
        {
            if (startLoop.IsClosed) startLoop.RemoveAt(0);
            if (endLoop.IsClosed) endLoop.RemoveAt(0);
            Point3d cen1 = new Point3d(), cen2 = new Point3d();
            for (int i = 0; i < startLoop.Count; i++)
            {
                cen1 += startLoop[i];
            }
            for (int i = 0; i < endLoop.Count; i++)
            {
                cen2 += endLoop[i];
            }
            cen1 /= startLoop.Count;
            cen2 /= endLoop.Count;        
            return joinLoops(startLoop.ToArray(), endLoop.ToArray(), new Plane(cen1, cen2 - cen1));
        }
        public Mesh joinLoops(Point3d[] x, Point3d[] y, Plane pln)
        {
            //桥接生成管子
            Transform xform = Transform.PlanarProjection(pln);
            Point3d[] pointdArray = new Point3d[x.Length];
            Point3d[] initialPoints = new Point3d[y.Length];
            for (int i = 0; i < x.Length; i++)
            {
                pointdArray[i] = x[i];
                pointdArray[i].Transform(xform);
            }
            for (int j = 0; j < y.Length; j++)
            {
                initialPoints[j] = y[j];
                initialPoints[j].Transform(xform);
            }
            Mesh mesh = new Mesh();
            mesh.Vertices.AddVertices(x);
            mesh.Vertices.AddVertices(y);
            int length1 = x.Length;
            int length2 = y.Length;
            int index1 = 0;
            int index2 = new Point3dList(initialPoints).ClosestIndex(pointdArray[index1]);
            int index3 = 0;
            for (int k = 0; (mesh.Faces.Count < (length1 + length2)) && (k < 200); k++)
            {
                int num10 = (index2 + index3) % length2;
                int num11 = (num10 + 1) % length2;
                int num12 = (index1 + 1) % length1;
                if (((pointdArray[index1].DistanceTo(initialPoints[num11]) > pointdArray[num12].DistanceTo(initialPoints[num10])) && (index1 != (length1 - 1))) || (index3 == length2))
                {
                    mesh.Faces.AddFace(index1, num10 + length1, num12);
                    index1 = num12;
                }
                else
                {
                    mesh.Faces.AddFace(index1, num10 + length1, num11 + length1);
                    index3++;
                }
            }
            return mesh;
        }    
        public Mesh QuadBridge(Polyline startLoop, Polyline endLoop ,double Offset,double radius)
        {
            // 桥接管子
            if (startLoop.IsClosed) startLoop.RemoveAt(0);
            if (endLoop.IsClosed) endLoop.RemoveAt(0);
            Point3d cen1 = new Point3d(), cen2 = new Point3d();
            for (int i = 0; i < startLoop.Count; i++)
            {
                cen1 += startLoop[i];
            }
            for (int i = 0; i < endLoop.Count; i++)
            {
                cen2 += endLoop[i];
            }
            cen1 /= startLoop.Count;
            cen2 /= endLoop.Count;
            Line l1 = new Line(cen1, cen2);

            return QuadBridge(startLoop.ToArray(), endLoop.ToArray(), l1, Offset, Offset, 0, radius, radius, false, false);
            //  endOffset startOffset 两端环切线距离
            //segmentLength中央等分换切线距离 endOffset startOffset为0时不起作用
        }
        public Mesh QuadBridge(Point3d[] startLoop, Point3d[] endLoop, Line centerLine, double startOffset, double endOffset, double segmentLength, double startRadius, double endRadius, bool startNaked, bool endNaked)
        {
            Mesh mesh = new Mesh();
            double radius = 0.5 * (startRadius + endRadius);
            Array.Reverse(endLoop);
            Vector3d direction = centerLine.Direction;
            direction.Unitize();
            Plane plane = new Plane(centerLine.From, direction);
            double length = centerLine.Length;
            Point3d[] pointdArray = new Point3d[startLoop.Length];
            for (int i = 0; i < startLoop.Length; i++)
            {
                pointdArray[i] = plane.ClosestPoint(startLoop[i]);
            }
            Point3d[] initialPoints = new Point3d[endLoop.Length];
            for (int i = 0; i < endLoop.Length; i++)
            {
                initialPoints[i] = plane.ClosestPoint(endLoop[i]);
            }
            Point3dList list = new Point3dList(initialPoints);
            if (startNaked)
            {
                double minvalue = double.MaxValue;
                for (int i = 0; i < initialPoints.Length; i++)
                {
                    double num8 = SignedAngle(Vector3d.VectorAngle(initialPoints[i] - plane.Origin, plane.XAxis, plane));
                    if (Math.Abs(num8) < Math.Abs(minvalue))
                    {
                        minvalue = num8;
                    }
                }
                plane.Rotate(-minvalue, plane.ZAxis);
                double num6 = Vector3d.VectorAngle(plane.XAxis, pointdArray[0] - plane.Origin, plane);
                for (int num9 = 0; num9 < startLoop.Length; num9++)
                {
                    startLoop[num9].Transform(Transform.Rotation(-num6, direction, plane.Origin));
                }
            }
            else if (endNaked)
            {
                double minvalue = double.MaxValue;
                for (int i = 0; i < pointdArray.Length; i++)
                {
                    double num13 = SignedAngle(Vector3d.VectorAngle(pointdArray[i] - plane.Origin, plane.XAxis, plane));
                    if (Math.Abs(num13) < Math.Abs(minvalue))
                    {
                        minvalue = num13;
                    }
                }
                plane.Rotate(-minvalue, plane.ZAxis);
                double num11 = Vector3d.VectorAngle(plane.XAxis, initialPoints[0] - plane.Origin, plane);
                for (int num14 = 0; num14 < endLoop.Length; num14++)
                {
                    endLoop[num14].Transform(Transform.Rotation(-num11, direction, plane.Origin));
                }
            }
            else
            {
                double minvalue = double.MaxValue;
                int index = -1;
                int num17 = -1;
                for (int num19 = 0; num19 < pointdArray.Length; num19++)
                {
                    int num20 = list.ClosestIndex(pointdArray[num19]);
                    double num21 = initialPoints[num20].DistanceToSquared(pointdArray[num19]);
                    if (num21 < minvalue)
                    {
                        index = num19;
                        num17 = num20;
                        minvalue = num21;
                    }
                }
                Point3d pointd = 0.5 * (pointdArray[index] + initialPoints[num17]);
                double angle = Vector3d.VectorAngle(plane.XAxis, pointd - plane.Origin, plane);
                plane.Rotate(angle, plane.ZAxis);
            }
            if ((endOffset == 0.0) && (startOffset == 0.0))
            {
                Mesh introduced61 = joinLoops(startLoop, endLoop, plane);
                return Quadrangulate3(introduced61, plane.ZAxis);
            }
            double maxvalue = 0.0;
            for (int k = 0; k < startLoop.Length; k++)
            {
                double num26 = plane.DistanceTo(startLoop[k]);
                if ((num26 > maxvalue) && (startOffset != 0.0))
                {
                    maxvalue = num26;
                }
            }
            startOffset += maxvalue;
            double num23 = 0.0;
            for (int m = 0; m < endLoop.Length; m++)
            {
                double num28 = length - plane.DistanceTo(endLoop[m]);
                if ((num28 > num23) && (endOffset != 0.0))
                {
                    num23 = num28;
                }
            }
            endOffset += num23;
            double num24 = 0.5 * (endOffset + startOffset);
            if ((length - (maxvalue + num23)) < (0.5 * num24))
            {
                Mesh introduced62 = joinLoops(startLoop, endLoop, plane);
                return Quadrangulate3(introduced62, plane.ZAxis);
            }
            int sideCount = Math.Min(startLoop.Length, endLoop.Length);
            double num30 = centerLine.Length - (endOffset + startOffset);
            Point3d pointd2 = plane.PointAt(radius * polygonRadius(sideCount), 0.0);
            double num31 = 6.2831853071795862 / ((double)sideCount);
            Point3d[] pointdArray3 = new Point3d[sideCount];
            for (int n = 0; n < sideCount; n++)
            {
                pointd2.Transform(Transform.Rotation(-num31, direction, plane.Origin));
                pointdArray3[n] = pointd2;
            }
            if (num30 < (num24 * 1.0))
            {
                Vector3d vectord2 = (Vector3d)((direction * length) * 0.5);
                Point3d[] pointdArray4 = new Point3d[sideCount];
                for (int num33 = 0; num33 < sideCount; num33++)
                {
                    pointdArray4[num33] = pointdArray3[num33] + vectord2;
                }
                Mesh introduced63 = joinLoops(startLoop, pointdArray4, plane);
                mesh.Append(Quadrangulate3(introduced63, plane.ZAxis));
                Mesh introduced64 = joinLoops(pointdArray4, endLoop, plane);
                mesh.Append(Quadrangulate3(introduced64, plane.ZAxis));
                return mesh;
            }
            int num34 = (segmentLength > 0.0) ? (((int)(num30 / segmentLength)) + 1) : 1;
            Point3d[] y = new Point3d[sideCount];
            Point3d[] x = new Point3d[sideCount];
            for (int num35 = 0; num35 < sideCount; num35++)
            {
                y[num35] = pointdArray3[num35] + ((Point3d)(direction * startOffset));
                x[num35] = pointdArray3[num35] + ((Point3d)(direction * (length - endOffset)));
            }
            if (startOffset == 0.0)
            {
                Mesh mesh2 = joinLoops(x, endLoop, plane);
                mesh.Append(Quadrangulate3(mesh2, plane.ZAxis));
                Mesh mesh3 = joinLoops(startLoop, x, plane);
                mesh.Append(Quadrangulate3(mesh3, plane.ZAxis));
                return mesh;
            }
            if (endOffset == 0.0)
            {
                Mesh mesh4 = joinLoops(startLoop, y, plane);
                mesh.Append(Quadrangulate3(mesh4, plane.ZAxis));
                Mesh mesh5 = joinLoops(y, endLoop, plane);
                mesh.Append(Quadrangulate3(mesh5, plane.ZAxis));
                return mesh;
            }
            Mesh introduced65 = joinLoops(startLoop, y, plane);
            Mesh other = Quadrangulate3(introduced65, plane.ZAxis);
            Mesh introduced66 = joinLoops(x, endLoop, plane);
            Mesh mesh7 = Quadrangulate3(introduced66, plane.ZAxis);
            mesh.Append(other);
            mesh.Append(mesh7);
            Mesh mesh8 = new Mesh();
            Vector3d vectord3 = (Vector3d)(direction * (num30 / ((double)num34)));
            for (int i = 0; i < (num34 + 1); i++)
            {
                for (int num37 = 0; num37 < sideCount; num37++)
                {
                    mesh8.Vertices.Add((pointdArray3[num37] + (i * vectord3)) + ((Point3d)(direction * startOffset)));
                }
                if (i > 0)
                {
                    int num38 = (i - 1) * sideCount;
                    int num39 = i * sideCount;
                    for (int num40 = 0; num40 < sideCount; num40++)
                    {
                        mesh8.Faces.AddFace(num38 + num40, num38 + ((num40 + 1) % sideCount), num39 + ((num40 + 1) % sideCount), num39 + num40);
                    }
                }
            }
            mesh.Append(mesh8);
            return mesh;
        }
        public double SignedAngle(double positiveAngle)
        {
            positiveAngle = positiveAngle % 6.2831853071795862;
            if (positiveAngle <= 3.1415926535897931)
            {
                return positiveAngle;
            }
            return (positiveAngle - 6.2831853071795862);
        }
        public Mesh Quadrangulate3(Mesh M, Vector3d axis)
        {
            //网格局部三角化
            bool flag = false;
            while (!flag)
            {
                MeshFace face;
                flag = true;
                int count = M.TopologyEdges.Count;
                double[] numArray = new double[count];
                for (int i = 0; i < count; i++)
                {
                    int[] connectedFaces = M.TopologyEdges.GetConnectedFaces(i);
                    if (connectedFaces.Length == 1)
                    {
                        numArray[i] = 0.0;
                    }
                    else
                    {
                        face = M.Faces[connectedFaces[0]];
                        if (face.IsTriangle)
                        {
                            face = M.Faces[connectedFaces[1]];
                            if (face.IsTriangle)
                            {
                                goto Label_008E;
                            }
                        }
                        numArray[i] = 0.0;
                    }
                    continue;
                Label_008E:
                    M.TopologyVertices.IndicesFromFace(connectedFaces[0]);
                    M.TopologyVertices.IndicesFromFace(connectedFaces[1]);
                    IndexPair topologyVertices = M.TopologyEdges.GetTopologyVertices(i);
                    int[] numArray3 = M.TopologyVertices.ConnectedTopologyVertices(topologyVertices.I, true);
                    int[] numArray4 = M.TopologyVertices.ConnectedTopologyVertices(topologyVertices.J, true);
                    if ((numArray3.Length == 3) || (numArray4.Length == 3))
                    {
                        numArray[i] = 0.0;
                    }
                    else
                    {
                        for (int j = 0; j < numArray3.Length; j++)
                        {
                            if (numArray3[j] == topologyVertices.J)
                            {
                                int topologyVertexIndex = numArray3[((j - 1) + numArray3.Length) % numArray3.Length];
                                int num5 = numArray3[((j + 1) + numArray3.Length) % numArray3.Length];
                                int[] numArray5 = M.TopologyVertices.MeshVertexIndices(topologyVertexIndex);
                                int[] numArray6 = M.TopologyVertices.MeshVertexIndices(num5);
                                Point3d[] pointdArray = new Point3d[] { M.TopologyEdges.EdgeLine(i).From, M.Vertices.Point3dAt(numArray5[0]), M.TopologyEdges.EdgeLine(i).To, M.Vertices.Point3dAt(numArray6[0]) };
                                Vector3d vectord = (Vector3d)(pointdArray[2] - pointdArray[0]);
                                Vector3d vectord2 = (Vector3d)(pointdArray[3] - pointdArray[1]);
                                vectord -= axis * (vectord * axis);
                                vectord2 -= axis * (vectord2 * axis);
                                double num6 = vectord.Length / vectord2.Length;
                                if (num6 > 1.0)
                                {
                                    num6 = 1.0 / num6;
                                }
                                numArray[i] = num6;
                                break;
                            }
                        }
                    }
                    if (numArray[i] > 0.0)
                    {
                        flag = false;
                    }
                }
                if (!flag)
                {
                    double num7 = 0.0;
                    int topologyEdgeIndex = -1;
                    for (int k = 0; k < numArray.Length; k++)
                    {
                        if (numArray[k] > num7)
                        {
                            num7 = numArray[k];
                            topologyEdgeIndex = k;
                        }
                    }
                    int[] faceIndexes = M.TopologyEdges.GetConnectedFaces(topologyEdgeIndex);
                    List<int> list = new List<int>();
                    int num9 = M.TopologyEdges.GetTopologyVertices(topologyEdgeIndex)[0];
                    int num10 = M.TopologyEdges.GetTopologyVertices(topologyEdgeIndex)[1];
                    num9 = M.TopologyVertices.MeshVertexIndices(num9)[0];
                    num10 = M.TopologyVertices.MeshVertexIndices(num10)[0];
                    for (int m = 0; m < 3; m++)
                    {
                        face = M.Faces[faceIndexes[0]];
                        if (face[m] == num9)
                        {
                            face = M.Faces[faceIndexes[0]];
                            int num13 = (face[(m + 1) % 3] == num10) ? 1 : 0;
                            face = M.Faces[faceIndexes[0]];
                            list.Add(face[(m + num13) % 3]);
                            face = M.Faces[faceIndexes[0]];
                            list.Add(face[((m + num13) + 1) % 3]);
                            face = M.Faces[faceIndexes[0]];
                            list.Add(face[((m + num13) + 2) % 3]);
                            break;
                        }
                    }
                    for (int n = 0; n < 3; n++)
                    {
                        face = M.Faces[faceIndexes[1]];
                        int item = face[n];
                        if (!list.Contains(item))
                        {
                            list.Add(item);
                        }
                    }
                    M.Faces.AddFace(list[0], list[1], list[2], list[3]);
                    M.Faces.DeleteFaces(faceIndexes);
                }
            }
            return M;
        }
        public Mesh Quadrangulate2(Mesh M, double maxRatio)
        {
            bool flag = false;
            while (!flag)
            {
                MeshFace face;
                M.FaceNormals.ComputeFaceNormals();
                flag = true;
                int count = M.TopologyEdges.Count;
                double[] numArray = new double[count];
                for (int i = 0; i < count; i++)
                {
                    int[] connectedFaces = M.TopologyEdges.GetConnectedFaces(i);
                    if (connectedFaces.Length == 1)
                    {
                        numArray[i] = 0.0;
                    }
                    else
                    {
                        face = M.Faces[connectedFaces[0]];
                        if (face.IsTriangle)
                        {
                            face = M.Faces[connectedFaces[1]];
                            if (face.IsTriangle)
                            {
                                goto Label_009A;
                            }
                        }
                        numArray[i] = 0.0;
                    }
                    continue;
                Label_009A:
                    numArray[i] = Math.Abs((double)(M.FaceNormals[connectedFaces[0]] * M.FaceNormals[connectedFaces[1]]));
                    IndexPair topologyVertices = M.TopologyEdges.GetTopologyVertices(i);
                    int[] numArray3 = M.TopologyVertices.ConnectedTopologyVertices(topologyVertices.J, true);
                    if ((M.TopologyVertices.ConnectedTopologyVertices(topologyVertices.I, true).Length == 3) || (numArray3.Length == 3))
                    {
                        numArray[i] = 0.0;
                    }
                    if (numArray[i] > 0.0)
                    {
                        flag = false;
                    }
                }
                if (!flag)
                {
                    bool[] flagArray;
                    int num5;
                    double num3 = 0.0;
                    int topologyEdgeIndex = -1;
                    for (int j = 0; j < numArray.Length; j++)
                    {
                        if (numArray[j] > num3)
                        {
                            num3 = numArray[j];
                            topologyEdgeIndex = j;
                        }
                    }
                    int[] faceIndexes = M.TopologyEdges.GetConnectedFaces(topologyEdgeIndex, out flagArray);
                    List<int> list = new List<int>();
                    if (!flagArray[0])
                    {
                        num5 = M.TopologyEdges.GetTopologyVertices(topologyEdgeIndex)[0];
                    }
                    else
                    {
                        num5 = M.TopologyEdges.GetTopologyVertices(topologyEdgeIndex)[1];
                    }
                    num5 = M.TopologyVertices.MeshVertexIndices(num5)[0];
                    for (int k = 0; k < 3; k++)
                    {
                        face = M.Faces[faceIndexes[0]];
                        if (face[k] == num5)
                        {
                            face = M.Faces[faceIndexes[0]];
                            list.Add(face[k]);
                            face = M.Faces[faceIndexes[0]];
                            list.Add(face[(k + 1) % 3]);
                            face = M.Faces[faceIndexes[0]];
                            list.Add(face[(k + 2) % 3]);
                            break;
                        }
                    }
                    for (int m = 0; m < 3; m++)
                    {
                        face = M.Faces[faceIndexes[1]];
                        int item = face[m];
                        if (!list.Contains(item))
                        {
                            list.Add(item);
                        }
                    }
                    M.Faces.AddFace(list[0], list[1], list[2], list[3]);
                    M.Faces.DeleteFaces(faceIndexes);
                }
            }
            return M;
        }
        public double polygonRadius(int sideCount)
        {
            double[] output = new double[] {
        1.5, 1.299254, 1.2, 1.143515, 1.108194, 1.084581, 1.067989, 1.055872, 1.046746,
                1.039697, 1.034137, 1.029673, 1.026034, 1.023028, 1.020515, 1.018393,
        1.016585, 1.015032, 1.013687, 1.012516};
            return output[Math.Min(sideCount - 4, 0x13)];
        }
    }
    public class HalfEdge
    {
        // Fields
        public int Index;
        public Node StartNode;
        public HalfEdge Twin;
        // Methods
        public Vector3d GetDirection(bool unitize)
        {
            Vector3d vectord = GetEndPoint() - GetStartPoint();
            if (unitize)
            {
                vectord.Unitize();
            }
            return vectord;
        }
        public Node GetEnd()
        {
            return Twin.StartNode;
        }
        public Point3d GetEndPoint()
        {
            return GetEnd().Position;
        }
        public int GetIndexInStartOutgoing()
        {
            return StartNode.OutgoingHalfEdges.IndexOf(this);
        }
        public Line GetLine()
        {
            return new Line(GetStartPoint(), GetEndPoint());
        }
        public Point3d GetStartPoint()
        {
            return StartNode.Position;
        }
    }
    public class HalfEdgeGraph
    {
        // Fields
        public List<HalfEdge> HalfEdges = new List<HalfEdge>();
        public List<Node> Nodes = new List<Node>();
        // Methods
        public void AddEdge(Line L, double tolerance)
        {
            Node node = AddNode(L.From, tolerance);
            Node node2 = AddNode(L.To, tolerance);
            if (node.Index != node2.Index)
            {
                HalfEdge item = new HalfEdge();
                HalfEdge edge2 = new HalfEdge();
                item.StartNode = node;
                edge2.StartNode = node2;
                item.Twin = edge2;
                edge2.Twin = item;
                item.Index = HalfEdges.Count;
                HalfEdges.Add(item);
                edge2.Index = HalfEdges.Count;
                HalfEdges.Add(edge2);
                node.OutgoingHalfEdges.Add(item);
                node2.OutgoingHalfEdges.Add(edge2);
            }
        }
        public Node AddNode(Point3d P, double tolerance)
        {
            int i = -1;
            for (int j = 0; j < Nodes.Count; j++)
            {
                if (Nodes[j].Position.DistanceTo(P) < tolerance)
                {
                    i = j;
                    break;
                }
            }
            if (i == -1)
            {
                i = Nodes.Count;
                Nodes.Add(new Node(P, i));
            }
            return Nodes[i];
        }
        public int GetNodeIndex(Point3d P, double tolerance)
        {
            for (int i = 0; i < Nodes.Count; i++)
            {
                if (Nodes[i].Position.DistanceTo(P) < tolerance)
                {
                    return i;
                }
            }
            return -1;
        }
        public List<Vector3d> GetNormals()
        {
            List<Vector3d> list = new List<Vector3d>();
            foreach (Node node in Nodes)
            {
                list.Add(node.Normal);
            }
            return list;
        }
        public List<Point3d> GetPositions()
        {
            List<Point3d> list = new List<Point3d>();
            foreach (Node node in Nodes)
            {
                list.Add(node.Position);
            }
            return list;
        }
    }
    public class Node
    {
        // Fields
        public int Index;
        public Vector3d Normal;
        public List<HalfEdge> OutgoingHalfEdges;
        public Point3d Position;
        // Methods
        public Node(Point3d P, int i)
        {
            Position = P;
            Index = i;
            OutgoingHalfEdges = new List<HalfEdge>();
            Normal = Vector3d.Zero;
        }
        public List<Point3d> GetNeighbourPositions()
        {
            List<Point3d> list = new List<Point3d>();
            foreach (Node node in GetNeighbours())
            {
                list.Add(node.Position);
            }
            return list;
        }
        public List<Node> GetNeighbours()
        {
            List<Node> list = new List<Node>();
            foreach (HalfEdge edge in OutgoingHalfEdges)
            {
                list.Add(edge.GetEnd());
            }
            return list;
        }
        public List<Vector3d> GetNeighbourVectors(bool unitize)
        {
            List<Vector3d> list = new List<Vector3d>();
            using (List<Point3d>.Enumerator enumerator = GetNeighbourPositions().GetEnumerator())
            {
                while (enumerator.MoveNext())
                {
                    Vector3d item = enumerator.Current - Position;
                    if (unitize)
                    {
                        item.Unitize();
                    }
                    list.Add(item);
                }
            }
            return list;
        }
        public int GetValence() { return OutgoingHalfEdges.Count; }
    }
}
