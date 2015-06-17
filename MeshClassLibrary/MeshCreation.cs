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
    public class MeshCreation
    {
        public MeshCreation() { }
        ///// MeshCreation     
        #region offset
        public Polyline QuadFaceOffset(Point3d p1, Point3d p2, Point3d p3, Point3d p4, Vector3d N, double distance)
        {
            Point3d cen = (p1 + p2 + p3 + p4) / 4;
            Line lcen = new Line(cen, cen + N);
            double u, v;
            Line l1 = new Line(p1, p2);
            Rhino.Geometry.Intersect.Intersection.LineLine(lcen, l1, out u, out v);
            Vector3d v1 = lcen.PointAt(u) - l1.PointAt(v);
            v1.Unitize(); v1 *= distance;
            l1.Transform(Transform.Translation(v1));
            Line l2 = new Line(p2, p3);
            Rhino.Geometry.Intersect.Intersection.LineLine(lcen, l2, out u, out v);
            v1 = lcen.PointAt(u) - l2.PointAt(v);
            v1.Unitize(); v1 *= distance;
            l2.Transform(Transform.Translation(v1));
            Line l3 = new Line(p3, p4);
            Rhino.Geometry.Intersect.Intersection.LineLine(lcen, l3, out u, out v);
            v1 = lcen.PointAt(u) - l3.PointAt(v);
            v1.Unitize(); v1 *= distance;
            l3.Transform(Transform.Translation(v1));
            Line l4 = new Line(p4, p1);
            Rhino.Geometry.Intersect.Intersection.LineLine(lcen, l4, out u, out v);
            v1 = lcen.PointAt(u) - l4.PointAt(v);
            v1.Unitize(); v1 *= distance;
            l4.Transform(Transform.Translation(v1));
            Polyline output = new Polyline();
            Rhino.Geometry.Intersect.Intersection.LineLine(l1, l4, out u, out v);
            output.Add((l1.PointAt(u) + l4.PointAt(v)) / 2);
            Rhino.Geometry.Intersect.Intersection.LineLine(l2, l1, out u, out v);
            output.Add((l2.PointAt(u) + l1.PointAt(v)) / 2);
            Rhino.Geometry.Intersect.Intersection.LineLine(l3, l2, out u, out v);
            output.Add((l3.PointAt(u) + l2.PointAt(v)) / 2);
            Rhino.Geometry.Intersect.Intersection.LineLine(l4, l3, out u, out v);
            output.Add((l4.PointAt(u) + l3.PointAt(v)) / 2);
            return output;
        }
        public Mesh QuadFaceOffset(Point3d p1, Point3d p2, Point3d p3, Point3d p4, double distance, bool type)
        {
            Mesh mesh = MeshFromPoints(p1, p2, p3, p4);
            Polyline pl2 = new Polyline();
            pl2.Add(p1);
            pl2.Add(p2);
            pl2.Add(p3);
            pl2.Add(p4);
            mesh.FaceNormals.ComputeFaceNormals();
            Polyline pl = QuadFaceOffset(p1, p2, p3, p4, mesh.FaceNormals[0], distance);
            if (type) return ClosedBridge(pl, pl2);
            else return MeshFromPoints(pl[0], pl[1], pl[2], pl[3]);
        }
        public Mesh MeshQuadFaceOffset(Mesh mesh, double Distance, bool type)
        {
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            Mesh meshoutput = new Mesh();
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                int[] index = vs.IndicesFromFace(i);
                if (index.Length == 4)
                {
                    meshoutput.Append(QuadFaceOffset(
                        vs[index[0]], vs[index[1]], vs[index[2]], vs[index[3]],
                        Distance, type));
                }
            }
            return meshoutput;
        }
        public Mesh MeshOffset(Polyline pl, double t, int n)
        {
            Polyline pl2; Mesh mesh = new Mesh();
            if (n < 1) return mesh;
            for (int i = 1; i < n; i++)
            {
                mesh.Append(MeshOffset(pl, t / n, out pl2));
                pl = pl2;
            }
            return mesh;
        }
        public Mesh MeshOffset(Polyline pl, double t, out Polyline pl2)
        {
            Mesh mesh = new Mesh(); pl2 = new Polyline();
            if (pl.Count < 3) return mesh;
            Point3d cen = new Point3d();
            for (int i = 1; i < pl.Count; i++)
            {
                cen += pl[i];
            }
            cen /= pl.Count - 1;
            for (int i = 0; i < pl.Count; i++)
            {
                Vector3d v = cen - pl[i]; v *= t;
                pl2.Add(pl[i] + v);
            }
            return MeshLoft(pl, pl2, false, false);
        }
        #endregion
        #region Topo
        public Mesh Topo1(Mesh x)
        {
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = x.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = x.TopologyVertices;
            List<Point3d> FaceC = new List<Point3d>();
            for (int i = 0; i < x.Faces.Count; i++)
            {
                Point3d f = new Point3d();
                if (x.Faces[i].IsQuad)
                {
                    f += x.Vertices[x.Faces[i].A];
                    f += x.Vertices[x.Faces[i].B];
                    f += x.Vertices[x.Faces[i].C];
                    f += x.Vertices[x.Faces[i].D];
                    f /= 4;
                }
                else if (x.Faces[i].IsTriangle)
                {
                    f += x.Vertices[x.Faces[i].A];
                    f += x.Vertices[x.Faces[i].B];
                    f += x.Vertices[x.Faces[i].C];
                    f /= 3;
                }
                FaceC.Add(f);
            }
            Mesh mesh = new Mesh();
            for (int i = 0; i < el.Count; i++)
            {
                if (el.GetConnectedFaces(i).Length == 2)
                {
                    int C = mesh.Vertices.Count;
                    mesh.Vertices.Add(vs[el.GetTopologyVertices(i).I]);
                    mesh.Vertices.Add(FaceC[el.GetConnectedFaces(i)[0]]);
                    mesh.Vertices.Add(vs[el.GetTopologyVertices(i).J]);
                    mesh.Vertices.Add(FaceC[el.GetConnectedFaces(i)[1]]);

                    mesh.Faces.AddFace(C, C + 1, C + 2, C + 3);
                }
                else if (el.GetConnectedFaces(i).Length == 1)
                {
                    int C = mesh.Vertices.Count;
                    mesh.Vertices.Add(vs[el.GetTopologyVertices(i).I]);
                    mesh.Vertices.Add(FaceC[el.GetConnectedFaces(i)[0]]);
                    mesh.Vertices.Add(vs[el.GetTopologyVertices(i).J]);
                    mesh.Faces.AddFace(C, C + 1, C + 2);
                }
            }
            mesh.UnifyNormals();
            return mesh;
        }
        // Doo-Sabin
        public List<Polyline> Topo2(Mesh mesh)
        {
            List<Polyline> pls = new List<Polyline>();
            List<List<Line>> pll = new List<List<Line>>();
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            List<Point3d> FaceC = new List<Point3d>();
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                Point3d f = new Point3d();
                if (mesh.Faces[i].IsQuad)
                {
                    f += mesh.Vertices[mesh.Faces[i].A];
                    f += mesh.Vertices[mesh.Faces[i].B];
                    f += mesh.Vertices[mesh.Faces[i].C];
                    f += mesh.Vertices[mesh.Faces[i].D];
                    f /= 4;
                }
                else if (mesh.Faces[i].IsTriangle)
                {
                    f += mesh.Vertices[mesh.Faces[i].A];
                    f += mesh.Vertices[mesh.Faces[i].B];
                    f += mesh.Vertices[mesh.Faces[i].C];
                    f /= 3;
                }
                FaceC.Add(f);
            }
            for (int i = 0; i < vs.Count; i++)
            {
                List<Line> ls = new List<Line>();
                pll.Add(ls);
            }
            for (int i = 0; i < el.Count; i++)
            {
                int[] index = el.GetConnectedFaces(i);
                if (index.Length == 2)
                {
                    Line l = new Line(FaceC[index[0]], FaceC[index[1]]);
                    int a = el.GetTopologyVertices(i).I;
                    int b = el.GetTopologyVertices(i).J;
                    pll[a].Add(l); pll[b].Add(l);
                }
            }
            for (int i = 0; i < vs.Count; i++)
            {
                //Print(pll[i].Count.ToString());
                Polyline cs = new Polyline();
                try
                {
                    cs = RoundCombine(pll[i]);
                }
                catch { }
                if (cs != null && cs.Count > 3) pls.Add(cs);
            }
            return pls;
        }
        private Polyline RoundCombine(List<Line> x)
        {
            List<IndexPair> id; List<BasicVertice> vs;
            BasicVertice.CreateCollection(x, out id, out vs);
            Polyline pl = new Polyline();
            for (int i = 0; i < vs.Count; i++)
            {
                if (vs[i].refer.Count != 2) return pl;
            }
            if (vs.Count == 2)
            {
                pl.Add(vs[0].pos);
                pl.Add(vs[1].pos);
            }
            if (vs.Count == 3)
            {
                pl.Add(vs[0].pos);
                pl.Add(vs[1].pos);
                pl.Add(vs[2].pos);
            }
            if (vs.Count > 3)
            {
                int start = 0; int second = 1;
                pl.Add(vs[start].pos); vs[start].energy = 1;
                pl.Add(vs[second].pos); vs[second].energy = 1;
                start = second;
                for (int i = 2; i < vs.Count; i++)
                {
                    if (vs[vs[second].refer[0]].energy == 0)
                    {
                        second = vs[second].refer[0];
                    }
                    else if (vs[vs[second].refer[1]].energy == 0)
                    {
                        second = vs[second].refer[1];
                    }
                    else { break; }
                    pl.Add(vs[second].pos); vs[second].energy = 1;
                    start = second;
                }
            }
            pl.Add(pl[0]);
            return pl;
        }
        //separate the poly in >2 degree vertice
        public List<Polyline> Topo3(List<Line> x)
        {
            List<IndexPair> id; List<BasicVertice> vs;
            BasicVertice.CreateCollection(x, out id, out vs);
            List<Polyline> pls = new List<Polyline>();
            for (int ii = 0; ii < vs.Count; ii++)
            {
                BasicVertice v1 = vs[ii];
                if (v1.refer.Count > 2)
                {
                    v1.energy = 1;
                    for (int i = 0; i < v1.refer.Count; i++)
                    {
                        int a = v1.refer[i];
                        if (vs[a].refer.Count != 2 && vs[a].energy == 0)
                        {
                            Polyline pl = new Polyline();
                            pl.Add(v1.pos); pl.Add(vs[a].pos);
                            pls.Add(pl);
                            //  Print(pl.Count.ToString() + "/0001");
                        }
                        else if (vs[a].energy == 0)
                        {
                            vs[a].energy = 1;
                            Polyline pl = new Polyline();
                            pl.Add(v1.pos); pl.Add(vs[a].pos);

                            while (vs[a].refer.Count == 2)
                            {
                                BasicVertice v2 = vs[vs[a].refer[0]];
                                if (v2.energy == 0) { a = vs[a].refer[0]; }
                                else { a = vs[a].refer[1]; }
                                if (vs[a].refer.Count == 2) vs[a].energy = 1;
                                pl.Add(vs[a].pos);
                            }
                            pls.Add(pl);
                            //   Print(pl.Count.ToString());
                        }
                    }
                }
            } return pls;
        }
        #endregion
        #region Analysis
        public double areaTri(Point3d p1, Point3d p2, Point3d p3)
        {
            double a = p1.DistanceTo(p2);
            double b = p1.DistanceTo(p3);
            double c = p2.DistanceTo(p3);
            double p = (a + b + c) / 2;
            return Math.Sqrt(p * (p - a) * (p - b) * (p - c));
        }

        public double MeshFaceArea(Mesh mesh)
        {        
            if (mesh.Faces.Count < 0 || mesh.Vertices.Count < 3) return 0;
            mesh.Faces.ConvertQuadsToTriangles();
            double t = 0;          
                for (int i = 0; i < mesh.Faces.Count; i++)
                {
                    Point3d p1 = mesh.Vertices[mesh.Faces[i].A];
                    Point3d p2 = mesh.Vertices[mesh.Faces[i].B];
                    Point3d p3 = mesh.Vertices[mesh.Faces[i].C];
                    t += areaTri(p1, p2, p3);
                }
            return t;
        }
        public List<Mesh> MeshExplode(Mesh mesh)
        {
            List<Mesh> meshes = new List<Mesh>();
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                Mesh mesh1 = new Mesh();
                mesh1.Vertices.Add(mesh.Vertices[mesh.Faces[i].A]);
                mesh1.Vertices.Add(mesh.Vertices[mesh.Faces[i].B]);
                mesh1.Vertices.Add(mesh.Vertices[mesh.Faces[i].C]);
                mesh1.Vertices.Add(mesh.Vertices[mesh.Faces[i].D]);
                mesh1.Faces.AddFace(0, 1, 2, 3);
                mesh1.Normals.ComputeNormals();
                meshes.Add(mesh1);
            }
            return meshes;
        }
        List<Point3d> MeshEdgeVertice(Mesh mesh)
        {
            List<Point3d> ls = new List<Point3d>();
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            List<bool> sign = new List<bool>();
            for (int i = 0; i < vs.Count; i++)
            {
                sign.Add(false);
            }
            for (int i = 0; i < el.Count; i++)
            {
                if (el.GetConnectedFaces(i).Length != 2)
                {
                    sign[el.GetTopologyVertices(i).I] = true;
                    sign[el.GetTopologyVertices(i).J] = true;
                }
            }
            for (int i = 0; i < vs.Count; i++)
            {
                if (sign[i]) ls.Add(vs[i]);
            }
            return ls;
        }
        public List<Line> MeshEdge(Mesh mesh)
        {
            List<Line> ls = new List<Line>();
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            for (int i = 0; i < el.Count; i++)
            {
                if (el.GetConnectedFaces(i).Length != 2)
                    ls.Add(el.EdgeLine(i));
            }
            return ls;
        }
        public List<Line> MeshEdge(Mesh mesh, double t)
        {
            List<Line> ls = new List<Line>();
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            Rhino.Geometry.Collections.MeshFaceNormalList ns = mesh.FaceNormals;
            if (!ns.ComputeFaceNormals()) return null;
            for (int i = 0; i < el.Count; i++)
            {
                if (el.GetConnectedFaces(i).Length == 1)
                {
                    ls.Add(el.EdgeLine(i));
                }
                else if (el.GetConnectedFaces(i).Length == 2)
                {
                    Vector3d v1 = ns[el.GetConnectedFaces(i)[0]];
                    Vector3d v2 = ns[el.GetConnectedFaces(i)[1]];
                    double y = Vector3d.VectorAngle(v1, v2);
                    if (y > t) { ls.Add(el.EdgeLine(i)); }
                }
            }
            return ls;
        }
        public List<Point3d> MeshProfileVertex(Mesh mesh)
        {
            return MeshProfileVertex(mesh, 0.001, 0.001);
        }
        public List<Point3d> MeshProfileVertex(Mesh mesh, double t1, double t2)
        {
            List<Line> ls = new List<Line>();
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            Rhino.Geometry.Collections.MeshFaceNormalList ns = mesh.FaceNormals;
            if (!ns.ComputeFaceNormals()) return null;
            for (int i = 0; i < el.Count; i++)
            {
                if (el.GetConnectedFaces(i).Length == 1)
                {
                    ls.Add(el.EdgeLine(i));
                }
                else if (el.GetConnectedFaces(i).Length == 2)
                {
                    Vector3d v1 = ns[el.GetConnectedFaces(i)[0]];
                    Vector3d v2 = ns[el.GetConnectedFaces(i)[1]];
                    double y = Vector3d.VectorAngle(v1, v2);
                    if (y > t1) { ls.Add(el.EdgeLine(i)); }
                }
            }
            List<IndexPair> id;
            List<BasicVertice> vs;
            List<Point3d> output = new List<Point3d>();
            BasicVertice.CreateCollection(ls, out id, out vs);
            for (int i = 0; i < vs.Count; i++)
            {
                if (vs[i].refer.Count != 2) { output.Add(vs[i].pos); }
                else
                {
                    Point3d pos1 = vs[vs[i].refer[0]].pos;
                    Point3d pos3 = vs[vs[i].refer[1]].pos;
                    Point3d pos2 = vs[i].pos;
                    Vector3d v1 = pos1 - pos2; Vector3d v2 = pos3 - pos2;
                    double y = Vector3d.VectorAngle(v1, v2);
                    if (y >= 3.1415926) y = 0;
                    if (y <= -3.1415926) y = 0;
                    // Print(y.ToString());
                    if (y > t2)
                    {
                        output.Add(vs[i].pos);
                    }
                }
            }
            return output;
        }
        public double Planeness(Point3d p1, Point3d p2, Point3d p3, Point3d p4)
        {
            Line l1 = new Line(p1, p3);
            Line l2 = new Line(p2, p4);
            return l1.MinimumDistanceTo(l2);
        }
        public List<Mesh> MeshPlaneness(Mesh mesh, out double Max, out List<double> values)
        {
            //min value=0
            List<Mesh> meshes = new List<Mesh>();
            List<double> t = new List<double>();
            double max = double.MinValue;
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                MeshFace f = mesh.Faces[i];
                if (f.IsTriangle) { t.Add(0); }
                else if (f.IsQuad)
                {
                    double value1 = Planeness(mesh.Vertices[f.A], mesh.Vertices[f.B], mesh.Vertices[f.C], mesh.Vertices[f.D]);
                    meshes.Add(MeshFromPoints(mesh.Vertices[f.A], mesh.Vertices[f.B], mesh.Vertices[f.C], mesh.Vertices[f.D]));
                    t.Add(value1);
                    if (value1 > max) max = value1;
                }
            }
            Max = max;
            values = t;
            for (int i = 0; i < t.Count; i++)
            {
                for (int j = 0; j < meshes[i].Vertices.Count; j++)
                {
                    double T = t[i] / max; double R; double G;
                    if (T >= 0.5)
                    {
                        R = 255.0;
                        G = 510.0 * (1 - T);
                    }
                    else
                    {
                        R = 510.0 * (T);
                        G = 255.0;
                    }
                    //R = 255 * T;G = 255 * (1 - T);
                    if (R > 255) R = 255; if (G > 255) G = 255;
                    if (G < 0) G = 0; if (R < 0) R = 0;
                    meshes[i].VertexColors.Add((int)R, (int)G, 0);
                }
            }
            return meshes;
        }
        public List<Mesh> MeshPlaneness(Mesh mesh)
        {
            double max; List<double> t;
            return MeshPlaneness(mesh, out max, out t);
        }
        public List<Mesh> MeshFaceDisplay(Mesh mesh, List<double> data)
        {
            //min value=0
            List<Mesh> meshes = new List<Mesh>();
            List<double> t = data;
            double max = double.MinValue;
            double min = double.MaxValue;
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                MeshFace f = mesh.Faces[i];
                if (f.IsTriangle)
                {
                    meshes.Add(MeshFromPoints(mesh.Vertices[f.A], mesh.Vertices[f.B], mesh.Vertices[f.C]));
                }
                else if (f.IsQuad)
                {
                    meshes.Add(MeshFromPoints(mesh.Vertices[f.A], mesh.Vertices[f.B], mesh.Vertices[f.C], mesh.Vertices[f.D]));
                }
                if (t[i] > max) max = t[i];
                if (t[i] < min) min = t[i];
            }

            for (int i = 0; i < t.Count; i++)
            {
                for (int j = 0; j < meshes[i].Vertices.Count; j++)
                {
                    double T = (t[i] - min) / (max - min); double R; double G;
                    if (T >= 0.5)
                    {
                        R = 255.0;
                        G = 510.0 * (1 - T);
                    }
                    else
                    {
                        R = 510.0 * (T);
                        G = 255.0;
                    }
                    //R = 255 * T;G = 255 * (1 - T);
                    if (R > 255) R = 255; if (G > 255) G = 255;
                    if (G < 0) G = 0; if (R < 0) R = 0;
                    meshes[i].VertexColors.Add((int)R, (int)G, 0);
                }
            }
            return meshes;
        }
        public List<Point3d> MeshFaceCenter(Mesh x)
        {
            List<Point3d> FaceC = new List<Point3d>();
            for (int i = 0; i < x.Faces.Count; i++)
            {
                Point3d f = new Point3d();
                if (x.Faces[i].IsQuad)
                {
                    f += x.Vertices[x.Faces[i].A];
                    f += x.Vertices[x.Faces[i].B];
                    f += x.Vertices[x.Faces[i].C];
                    f += x.Vertices[x.Faces[i].D];
                    f /= 4;
                }
                else if (x.Faces[i].IsTriangle)
                {
                    f += x.Vertices[x.Faces[i].A];
                    f += x.Vertices[x.Faces[i].B];
                    f += x.Vertices[x.Faces[i].C];
                    f /= 3;
                }
                FaceC.Add(f);
            }
            return FaceC;
        }
        public List<Rhino.Display.Text3d> GeoID(List<Object> x)
        {
            List<Rhino.Display.Text3d> output = new List<Rhino.Display.Text3d>();
            for (int i = 0; i < x.Count; i++)
            {
                Point3d p;
                // index.Add(i);
                if (x[i].GetType() == typeof(Line))
                {
                    Line l = (Line)(x[i]);
                    p = (l.From + l.To) / 2;
                }
                else
                {
                    Rhino.Geometry.GeometryBase base1 = (Rhino.Geometry.GeometryBase)x[i];
                    // pos.Add(base1.GetBoundingBox(true).Center);

                    p = base1.GetBoundingBox(true).Center;
                }

                Rhino.Display.Text3d text = new Rhino.Display.Text3d(i.ToString(), new Plane(p, Vector3d.ZAxis), 1);
                output.Add(text);
            }
            return output;

        }
        public List<Rhino.Display.Text3d> MeshID(List<Mesh> x)
        {
            List<Rhino.Display.Text3d> output = new List<Rhino.Display.Text3d>();
            for (int i = 0; i < x.Count; i++)
            {
                Point3d p;
                p = x[i].GetBoundingBox(true).Center;
                Rhino.Display.Text3d text = new Rhino.Display.Text3d(i.ToString(), new Plane(p, Vector3d.ZAxis), 1);
                output.Add(text);
            }
            return output;
        }
        public List<Rhino.Display.Text3d> MeshFaceID(Mesh x)
        {
            List<Rhino.Display.Text3d> output = new List<Rhino.Display.Text3d>();
            for (int i = 0; i < x.Faces.Count; i++)
            {
                Point3d f = new Point3d();
                if (x.Faces[i].IsQuad)
                {
                    f += x.Vertices[x.Faces[i].A];
                    f += x.Vertices[x.Faces[i].B];
                    f += x.Vertices[x.Faces[i].C];
                    f += x.Vertices[x.Faces[i].D];
                    f /= 4;
                }
                else if (x.Faces[i].IsTriangle)
                {
                    f += x.Vertices[x.Faces[i].A];
                    f += x.Vertices[x.Faces[i].B];
                    f += x.Vertices[x.Faces[i].C];
                    f /= 3;
                }
                Rhino.Display.Text3d te = new Rhino.Display.Text3d(i.ToString(), new Plane(f, Vector3d.ZAxis), 1);
                output.Add(te);
            }
            return output;
        }
        public List<Rhino.Display.Text3d> MeshEdgeID(Mesh x)
        {
            List<Rhino.Display.Text3d> output = new List<Rhino.Display.Text3d>();
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = x.TopologyEdges;
            for (int i = 0; i < el.Count; i++)
            {
                Point3d f = new Point3d();
                Line l1 = el.EdgeLine(i);
                f = l1.PointAt(0.5);
                Rhino.Display.Text3d te = new Rhino.Display.Text3d(i.ToString(), new Plane(f, Vector3d.ZAxis), 1);
                output.Add(te);
            }
            return output;
        }
        public List<Rhino.Display.Text3d> MeshVerticeID(Mesh x)
        {
            List<Rhino.Display.Text3d> output = new List<Rhino.Display.Text3d>();
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = x.TopologyVertices;
            for (int i = 0; i < vs.Count; i++)
            {
                Point3d f = vs[i];
                Rhino.Display.Text3d te = new Rhino.Display.Text3d(i.ToString(), new Plane(f, Vector3d.ZAxis), 1);
                output.Add(te);
            }
            return output;
        }
        public List<Rhino.Display.Text3d> MeshVerticeData(Mesh x, List<double> data)
        {
            List<Rhino.Display.Text3d> output = new List<Rhino.Display.Text3d>();
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = x.TopologyVertices;
            if (data.Count < vs.Count) return output;
            for (int i = 0; i < vs.Count; i++)
            {
                Point3d f = vs[i];
                Rhino.Display.Text3d te = new Rhino.Display.Text3d(data[i].ToString(), new Plane(f, Vector3d.ZAxis), 1);
                output.Add(te);
            }
            return output;
        }
        public Mesh MeshVerticeDisplay(Mesh mesh, List<double> data)
        {
            Mesh output = new Mesh();
            output.Append(mesh);
            double max = double.MinValue;
            double min = double.MaxValue;
            List<double> t = data;
            Rhino.Geometry.Collections.MeshVertexList vs = mesh.Vertices;

            if (data.Count < vs.Count) return output;

            for (int i = 0; i < vs.Count; i++)
            {
                if (t[i] > max) max = t[i];
                if (t[i] < min) min = t[i];
            }
            for (int i = 0; i < t.Count; i++)
            {
                double T;
                if (max == min) { T = 0; }
                else { T = (t[i] - min) / (max - min); }
                double R; double G;
                if (T >= 0.5)
                {
                    R = 255.0;
                    G = 510.0 * (1 - T);
                }
                else
                {
                    R = 510.0 * (T);
                    G = 255.0;
                }
                //R = 255 * T;G = 255 * (1 - T);
                if (R > 255) R = 255;
                if (G > 255) G = 255;
                if (G < 0) G = 0;
                if (R < 0) R = 0;
                output.VertexColors.Add((int)R, (int)G, 0);
            }
            return output;
        }
        #endregion
        #region clean
        public void MeshClean(ref Mesh mesh, double tolerance)
        {
            try
            {
                List<int> mapping = new List<int>();
                List<Point3d> ps = new List<Point3d>();
                List<double> MapCount = new List<double>();
                ps.Add(mesh.Vertices[0]); mapping.Add(0); MapCount.Add(1);
                for (int j = 1; j < mesh.Vertices.Count; j++)
                {
                    double min = double.MaxValue;
                    int sign = 0;

                    for (int i = 0; i < ps.Count; i++)
                    {
                        double tempt = ps[i].DistanceTo((Point3d)mesh.Vertices[j]);
                        if (tempt < min)
                        {
                            min = tempt; sign = i;
                        }
                    }
                    if (min < tolerance)
                    {
                        mapping.Add(sign);
                        MapCount[sign]++;
                        ps[sign] *= (MapCount[sign] - 1) / MapCount[sign];
                        Point3d tempp = (Point3d)(mesh.Vertices[j]); tempp *= (1 / MapCount[sign]);
                        ps[sign] += tempp;
                    }
                    else { mapping.Add(ps.Count); ps.Add(mesh.Vertices[j]); MapCount.Add(1); }
                }
                Mesh mesh2 = new Mesh();
                for (int i = 0; i < ps.Count; i++)
                {
                    mesh2.Vertices.Add(ps[i]);
                }
                for (int i = 0; i < mesh.Faces.Count; i++)
                {
                    if (mesh.Faces[i].IsQuad)
                    {
                        int p1 = mapping[mesh.Faces[i].A];
                        int p2 = mapping[mesh.Faces[i].B];
                        int p3 = mapping[mesh.Faces[i].C];
                        int p4 = mapping[mesh.Faces[i].D];
                        if (noRepeat(p1, p2, p3, p4))
                        {
                            mesh2.Faces.AddFace(p1, p2, p3, p4);
                        }
                    }
                    if (mesh.Faces[i].IsTriangle)
                    {
                        int p1 = mapping[mesh.Faces[i].A];
                        int p2 = mapping[mesh.Faces[i].B];
                        int p3 = mapping[mesh.Faces[i].C];
                        if (noRepeat(p1, p2, p3))
                        {
                            mesh2.Faces.AddFace(p1, p2, p3);
                        }
                    }
                }
                mesh2.Compact();
                mesh2.UnifyNormals();
                mesh = mesh2;
            }
            catch (Exception ex)
            {
                System.Windows.Forms.MessageBox.Show(ex.ToString());
            }
        }
        private bool noRepeat(int a1, int a2, int a3, int a4)
        {
            if (a1 == a2) return false;
            else if (a1 == a3) return false;
            else if (a1 == a4) return false;
            else if (a2 == a3) return false;
            else if (a2 == a4) return false;
            else if (a3 == a4) return false;
            else return true;
        }
        private bool noRepeat(int a1, int a2, int a3)
        {
            if (a1 == a2) return false;
            else if (a1 == a3) return false;
            else if (a2 == a3) return false;
            else return true;
        }
        #endregion
        #region shape
        public Mesh MeshPlannar(Polyline pl)
        {
            Point3d cen = new Point3d();
            Mesh mesh = new Mesh();
            mesh.Vertices.Add(pl[0]);
            for (int i = 1; i < pl.Count; i++)
            {
                cen += pl[i];
                mesh.Vertices.Add(pl[i]);
                mesh.Faces.AddFace(pl.Count, i, i - 1);
            }
            cen /= pl.Count - 1;
            mesh.Vertices.Add(cen);
            mesh.Normals.ComputeNormals();
            return mesh;
        }
        public Mesh MeshTorus(Circle c, double t)
        {
            double cut = 64;
            List<Polyline> ls2 = new List<Polyline>();
            for (int i = 0; i < cut; i++)
            {
                List<Point3d> ls = new List<Point3d>();
                ls.Add(new Point3d(1, 0, 0));
                ls.Add(new Point3d(0.707, 0, 0.707));
                ls.Add(new Point3d(0, 0, 1));
                ls.Add(new Point3d(-0.707, 0, 0.707));
                ls.Add(new Point3d(-1, 0, 0));
                ls.Add(new Point3d(-0.707, 0, -0.707));
                ls.Add(new Point3d(0, 0, -1));
                ls.Add(new Point3d(0.707, 0, -0.707));
                Polyline l1 = new Polyline();
                l1.AddRange(ls);
                l1.Reverse();
                l1.Transform(Transform.Scale(Point3d.Origin, t));
                //l1.Transform(Transform.Translation(new Vector3d(c.Radius - t, 0, 0)));
                l1.Transform(Transform.Translation(new Vector3d(c.Radius, 0, 0)));
                l1.Transform(Transform.Rotation(Math.PI / (double)cut * 2 * i, Point3d.Origin));
                ls2.Add(l1);
            }
            Mesh mesh1 = MeshLoft(ls2, true, true);
            mesh1.Transform(Transform.PlaneToPlane(Plane.WorldXY, c.Plane));
            return mesh1;
        }
        public Mesh MeshSweep1(Curve l, Polyline ls, Plane SectionPos, int Count)
        {
            List<Polyline> ps = new List<Polyline>();
            Mesh mesh = new Mesh();
            double[] div = l.DivideByCount(Count, true);

            for (int i = 0; i < div.Length; i++)
            {
                Polyline l1 = new Polyline(ls);
                Plane plane;
                if (l.PerpendicularFrameAt(div[i], out plane))
                {
                    l1.Transform(Transform.PlaneToPlane(SectionPos, plane));
                    ps.Add(l1);
                }
            }
            mesh.Append(MeshLoft(ps, false, false));
            return mesh;
        }
        public Mesh MeshSweep1(Curve l, Polyline ls, int Count)
        {
            Mesh mesh = new Mesh(); if (ls.Count < 1) return mesh;
            Plane plane1;
            Plane.FitPlaneToPoints(ls.ToArray(), out plane1);
            Point3d Origin = ls[0];
            Vector3d v1 = ls[ls.Count - 1] - ls[0];
            v1.Unitize();
            Vector3d v2 = Vector3d.CrossProduct(v1, plane1.Normal);
            plane1 = new Plane(Origin, v2, v1);
            return MeshSweep1(l, ls, plane1, Count);
        }
        public Mesh MeshSweep1(Curve l, Curve l2, int Count, int Count2)
        {
            Polyline ls = new Polyline();
            double[] div = l2.DivideByCount(Count2, true);
            for (int i = 0; i < div.Length; i++)
            {
                ls.Add(l2.PointAt(div[i]));
            }
            return MeshSweep1(l, ls, Count);
        }
        public Mesh MeshPipe(Curve l, double t, int Count)
        {
            Mesh mesh = new Mesh();
            List<Point3d> ls = new List<Point3d>();
            ls.Add(new Point3d(1, 0, 0));
            ls.Add(new Point3d(0.707, 0.707, 0));
            ls.Add(new Point3d(0, 1, 0));
            ls.Add(new Point3d(-0.707, 0.707, 0));
            ls.Add(new Point3d(-1, 0, 0));
            ls.Add(new Point3d(-0.707, -0.707, 0));
            ls.Add(new Point3d(0, -1, 0));
            ls.Add(new Point3d(0.707, -0.707, 0));
            ls.Add(new Point3d(1, 0, 0));
            Polyline l1 = new Polyline(ls);
            l1.Transform(Transform.Scale(Point3d.Origin, t));
            return MeshSweep1(l, l1, Plane.WorldXY, Count);
        }
        public Mesh MeshPipe(Line l, double t)
        {
            Mesh mesh = new Mesh();
            List<Point3d> ls = new List<Point3d>();
            ls.Add(new Point3d(1, 0, 0));
            ls.Add(new Point3d(0.707, 0.707, 0));
            ls.Add(new Point3d(0, 1, 0));
            ls.Add(new Point3d(-0.707, 0.707, 0));
            ls.Add(new Point3d(-1, 0, 0));
            ls.Add(new Point3d(-0.707, -0.707, 0));
            ls.Add(new Point3d(0, -1, 0));
            ls.Add(new Point3d(0.707, -0.707, 0));
            Polyline l1 = new Polyline(ls);
            Polyline l2 = new Polyline(ls);
            l1.Transform(Transform.Scale(new Point3d(0, 0, 0), t));
            l2.Transform(Transform.Scale(new Point3d(0, 0, 0), t));
            Vector3d v = l.To - l.From;
            l1.Transform(Transform.PlaneToPlane(Plane.WorldXY, new Plane(l.From, v)));
            l2.Transform(Transform.PlaneToPlane(Plane.WorldXY, new Plane(l.To, v)));
            mesh.Append(MeshLoft(l1, l2, true, false));
            return mesh;
        }
        public Mesh MeshPipe(Line l, double t1, double t2)
        {
            Mesh mesh = new Mesh();
            List<Point3d> ls = new List<Point3d>();
            ls.Add(new Point3d(1, 0, 0));
            ls.Add(new Point3d(0.707, 0.707, 0));
            ls.Add(new Point3d(0, 1, 0));
            ls.Add(new Point3d(-0.707, 0.707, 0));
            ls.Add(new Point3d(-1, 0, 0));
            ls.Add(new Point3d(-0.707, -0.707, 0));
            ls.Add(new Point3d(0, -1, 0));
            ls.Add(new Point3d(0.707, -0.707, 0));
            Polyline l1 = new Polyline(ls);
            Polyline l2 = new Polyline(ls);
            l1.Transform(Transform.Scale(new Point3d(0, 0, 0), t1));
            l2.Transform(Transform.Scale(new Point3d(0, 0, 0), t2));
            Vector3d v = l.To - l.From;
            l1.Transform(Transform.PlaneToPlane(Plane.WorldXY, new Plane(l.From, v)));
            l2.Transform(Transform.PlaneToPlane(Plane.WorldXY, new Plane(l.To, v)));
            mesh.Append(MeshLoft(l1, l2, true, false));
            return mesh;
        }
        public Mesh MeshBeam(Line Line, double height, double width, Vector3d N)
        {
            List<Point3d> ps = new List<Point3d>();
            ps.Add(new Point3d(height, width, 0));
            ps.Add(new Point3d(-height, width, 0));
            ps.Add(new Point3d(-height, -width, 0));
            ps.Add(new Point3d(height, -width, 0));
            Polyline l1 = new Polyline(ps);
            Polyline l2 = new Polyline(ps);
            Plane ori = Plane.WorldXY;
            Vector3d vy = Line.To - Line.From;
            Vector3d vx = Vector3d.CrossProduct(vy, N);

            Plane plane1 = new Plane(Line.From, N, vx);
            Plane plane2 = new Plane(Line.To, N, vx);
            l1.Transform(Transform.PlaneToPlane(ori, plane1));
            l2.Transform(Transform.PlaneToPlane(ori, plane2));
            Mesh mesh = new Mesh();
            mesh.Vertices.AddVertices(l1);
            mesh.Vertices.AddVertices(l2);
            mesh.Faces.AddFace(0, 1, 2, 3);
            mesh.Faces.AddFace(4, 5, 6, 7);
            mesh.Append(MeshLoft(l1, l2, true, false));
            return mesh;
        }
        #endregion
        #region loft
        public Mesh ClosedBridge(Polyline pl1, Polyline pl2)
        {
            if (pl1.Count != pl2.Count) return null;
            if (pl1.Count < 2 || pl2.Count < 2) return null;
            if (pl1[0].DistanceTo(pl1[pl1.Count - 1]) < 0.01) { pl1.RemoveAt(pl1.Count - 1); }
            if (pl2[0].DistanceTo(pl2[pl2.Count - 1]) < 0.01) { pl2.RemoveAt(pl2.Count - 1); }
            Polyline pl3 = new Polyline();
            int sign = 0; double min = double.MaxValue;
            for (int i = 0; i < pl2.Count; i++)
            {
                double Length = 0;
                for (int j = 0; j < pl2.Count; j++)
                {
                    int index = i + j;
                    if (index > pl2.Count - 1) index -= pl2.Count;
                    Length += pl2[index].DistanceTo(pl1[j]);
                }
                if (Length < min) { min = Length; sign = i; }
            }
            for (int j = 0; j < pl2.Count; j++)
            {
                int index = sign + j;
                if (index > pl2.Count - 1) index -= pl2.Count;
                pl3.Add(pl2[index]);
            }
            return MeshLoft(pl1, pl3, true, false);
        }
        public Mesh MeshLoft(List<Curve> cs, int t, bool isClosed)
        {
            if (t < 1) return null;
            List<Polyline> pls = new List<Polyline>();
            for (int i = 0; i < cs.Count; i++)
            {
                Polyline pl = new Polyline();
                double[] nodes = cs[i].DivideByCount(t, true);
                for (int j = 0; j < nodes.Length; j++)
                {
                    pl.Add(cs[i].PointAt(nodes[j]));
                }
                pls.Add(pl);
            }
            return MeshLoft(pls, false, isClosed);
        }
        public Mesh MeshLoft(Curve c1, Curve c2, int t)
        {
            return MeshLoft(c1, c2, t, t);
        }
        public Mesh MeshLoft(Curve c1, Curve c2, int u, int v)
        {
            double[] uList1 = c1.DivideByCount(u, true);
            double[] uList2 = c2.DivideByCount(u, true);
            List<Polyline> input1 = new List<Polyline>();
            for (int i = 0; i < uList1.Length; i++)
            {
                Point3d p1 = c1.PointAt(uList1[i]);
                Point3d p2 = c2.PointAt(uList2[i]);
                Vector3d V = p2 - p1;
                V /= v;
                Polyline pl = new Polyline();
                for (int k = 0; k < v + 1; k++)
                {
                    pl.Add(p1 + (V * k));
                }
                input1.Add(pl);
            }
            return MeshLoft(input1, false, false);
        }
        public Mesh MeshLoft(Polyline pl1, Polyline pl2, bool isPolyClosed, bool isClosed)
        {

            List<Polyline> pls = new List<Polyline>();
            pls.Add(pl1);
            pls.Add(pl2);
            return MeshLoft(pls, isPolyClosed, isClosed);
        }
        public Mesh MeshLoft(List<Polyline> pl, bool isPolyClosed, bool isClosed)
        {
            int U = pl[0].Count;
            int V = pl.Count;
            if (isPolyClosed) { U++; }
            if (isClosed) { pl.Add(pl[0]); V++; }
            List<Point3d> pls = new List<Point3d>();

            for (int i = 0; i < pl.Count; i++)
            {
                for (int j = 0; j < pl[i].Count; j++)
                {
                    pls.Add(pl[i][j]);
                }
                if (isPolyClosed) { pls.Add(pl[i][0]); }
            }
            return MeshFromPoints(pls, U, V);
        }
        public Mesh MeshLoft(Line l1, Line l2)
        {
            return MeshFromPoints(l1.From, l1.To, l2.To, l2.From);
        }
        public Mesh MeshLoft(List<Line> ls, bool isClosed)
        {
            Polyline l1 = new Polyline(), l2 = new Polyline();
            for (int i = 0; i < ls.Count; i++)
            {
                l1.Add(ls[i].From); l2.Add(ls[i].To);
            }
            return MeshLoft(l1, l2, isClosed, false);
        }
        public Mesh TriangleMeshLoft(Polyline pl1, Point3d c2)
        {
            Mesh mesh = new Mesh();
            mesh.Vertices.Add(c2); mesh.Vertices.AddVertices(pl1);
            for (int i = 1; i < pl1.Count; i++)
            {
                mesh.Faces.AddFace(0, i, i + 1);
            }
            return mesh;
        }
        public Mesh TriangleMeshLoft(Curve c1, Point3d c2, int t)
        {
            List<Point3d> ps = new List<Point3d>();
            ps.Add(c2);
            for (int i = 2; i < t + 2; i++)
            {
                double[] t1 = c1.DivideByCount(i - 1, true);
                for (int k = 0; k < t1.Length; k++)
                {
                    Vector3d v = c1.PointAt(t1[k]) - c2;
                    v *= ((double)i - 1) / (double)t;
                    ps.Add(c2 + v);
                }
            }
            return TriangleMeshFromPoints(ps, t + 1);
        }
        public Mesh TriangleMeshLoft(Curve c1, Curve c2, int t1, int t2)
        {//t2>t1;
            Vector3d v1 = c1.PointAtEnd - c1.PointAtStart;
            Vector3d v2 = c2.PointAtEnd - c2.PointAtStart;
            if (Vector3d.VectorAngle(v1, v2) > Math.PI / 2) { c2.Reverse(); }

            List<Point3d> ps = new List<Point3d>();
            double[] t0 = c2.DivideByCount(t2 - 1, true);
            for (int i = 0; i < t0.Length; i++)
            { ps.Add(c2.PointAt(t0[i])); }

            for (int i = t2; i < t1; i++)
            {
                t0 = c2.DivideByCount(i, true);
                double[] t01 = c1.DivideByCount(i, true);
                for (int k = 0; k < t01.Length; k++)
                {
                    Vector3d v = c1.PointAt(t01[k]) - c2.PointAt(t0[k]);
                    v *= (double)((i - t2 + 1)) / (double)(t1 - t2);
                    ps.Add(c2.PointAt(t0[k]) + v);
                }
            }
            return TriangleMeshFromPoints(ps, t2, t1);
        }
        public Mesh TriangleMeshLoft2(Curve c1, Curve c2, int t1, int t2)
        {
            Vector3d v1 = c1.PointAtEnd - c1.PointAtStart;
            Vector3d v2 = c2.PointAtEnd - c2.PointAtStart;
            if (Vector3d.VectorAngle(v1, v2) > Math.PI / 2) { c2.Reverse(); }
            Mesh mesh = new Mesh();
            List<Point3d> ps = new List<Point3d>();
            int Count = t1 - t2;
            double[] t01 = c1.DivideByCount(Count, true);
            double[] t02 = c2.DivideByCount(Count, true);
            for (int i = t2; i <= t1; i++)
            {
                Point3d p1 = c1.PointAt(t01[i - t2]);
                Point3d p2 = c2.PointAt(t02[i - t2]);
                Vector3d v = p2 - p1;
                for (int k = 0; k < i; k++)
                {
                    double t3 = 0;
                    if (i > 1) { t3 = ((double)k / (double)(i - 1)); }
                    ps.Add(p1 + v * t3);
                }
            }
            if (t2 > 1)
            {
                return TriangleMeshFromPoints(ps, t2, t1);
            }
            if (t2 == 1) { return TriangleMeshFromPoints(ps); }
            else return mesh;
        }
        public Mesh TriangleMeshLoft2(Curve c1, Curve c2, int t)
        {
            return TriangleMeshLoft2(c1, c2, t, 1);
        }
        public Mesh TriangleMeshLoft(Polyline pl1, Polyline pl2)
        {// pl1.Count=pl2.Count+1;
            Mesh mesh = new Mesh();
            Polyline poly1, poly2;
            if (pl1.Count == pl2.Count + 1) { poly2 = new Polyline(pl1); poly1 = new Polyline(pl2); }
            else if (pl1.Count == pl2.Count - 1) { poly1 = new Polyline(pl1); poly2 = new Polyline(pl2); }
            else { return mesh; }
            for (int i = 0; i < poly1.Count; i++)
            {
                mesh.Vertices.Add(poly2[i]);
                mesh.Vertices.Add(poly1[i]);
            }
            mesh.Vertices.Add(poly2[poly2.Count - 1]);
            for (int i = 0; i < mesh.Vertices.Count - 2; i++)
            {
                mesh.Faces.AddFace(i, i + 1, i + 2);
            }
            return mesh;
        }
        #endregion
        #region extrute
        public Mesh MeshExtrute(Mesh meshOral, Vector3d v)
        {
            Mesh mesh = new Mesh();
            mesh.Append(MeshExtruteEdge(meshOral, v));
            mesh.Append(MeshOffset(meshOral, v));
            return mesh;
        }
        public Mesh MeshExtrute(Mesh meshOral, double v)
        {
            Mesh mesh = new Mesh();
            mesh.Append(MeshExtruteEdge(meshOral, v));
            mesh.Append(MeshOffset(meshOral, v));
            return mesh;
        }
        public Mesh MeshExtruteEdge(Mesh meshOral, double v)
        {
            Mesh mesh = new Mesh();
            meshOral.UnifyNormals();
            meshOral.Compact();
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = meshOral.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = meshOral.TopologyVertices;
            for (int i = 0; i < el.Count; i++)
            {
                if (el.GetConnectedFaces(i).Length == 1)
                {
                    int p1 = el.GetTopologyVertices(i).I;
                    int p2 = el.GetTopologyVertices(i).J;
                    Vector3d v1 = new Vector3d(meshOral.Normals[vs.MeshVertexIndices(p1)[0]]);
                    Vector3d v2 = new Vector3d(meshOral.Normals[vs.MeshVertexIndices(p2)[0]]);
                    int VS = mesh.Vertices.Count;
                    mesh.Vertices.Add(vs[p1]);
                    mesh.Vertices.Add(vs[p2]);
                    mesh.Vertices.Add(new Point3d(vs[p2]) + v2 * v);
                    mesh.Vertices.Add(new Point3d(vs[p1]) + v1 * v);
                    mesh.Faces.AddFace(VS, VS + 1, VS + 2, VS + 3);
                }
            }
            mesh.UnifyNormals();
            return mesh;
        }
        public Mesh MeshExtruteEdge(Mesh meshOral, Vector3d v)
        {
            Mesh mesh = new Mesh();
            meshOral.UnifyNormals();
            meshOral.Compact();
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = meshOral.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = meshOral.TopologyVertices;
            for (int i = 0; i < el.Count; i++)
            {
                if (el.GetConnectedFaces(i).Length == 1)
                {
                    int p1 = el.GetTopologyVertices(i).I;
                    int p2 = el.GetTopologyVertices(i).J;
                    int VS = mesh.Vertices.Count;
                    mesh.Vertices.Add(vs[p1]);
                    mesh.Vertices.Add(vs[p2]);
                    mesh.Vertices.Add(new Point3d(vs[p2]) + v);
                    mesh.Vertices.Add(new Point3d(vs[p1]) + v);
                    mesh.Faces.AddFace(VS, VS + 1, VS + 2, VS + 3);
                }
            }
            mesh.UnifyNormals();
            return mesh;
        }
        public Mesh MeshOffset(Mesh meshOral, double v)
        {
            Mesh mesh = new Mesh();
            meshOral.UnifyNormals();
            meshOral.Compact();
            for (int i = 0; i < meshOral.Vertices.Count; i++)
            {
                mesh.Vertices.Add(new Point3d(meshOral.Vertices[i]) + new Vector3d(meshOral.Normals[i]) * v);
            }
            for (int i = 0; i < meshOral.Faces.Count; i++)
            {
                mesh.Faces.AddFace(meshOral.Faces[i]);
            }
            mesh.UnifyNormals();
            return mesh;
        }
        public Mesh MeshOffset(Mesh meshOral, Vector3d v)
        {
            Mesh mesh = new Mesh();
            meshOral.Compact();
            for (int i = 0; i < meshOral.Vertices.Count; i++)
            {
                mesh.Vertices.Add(new Point3d(meshOral.Vertices[i]) + v);
            }
            for (int i = 0; i < meshOral.Faces.Count; i++)
            {
                mesh.Faces.AddFace(meshOral.Faces[i]);
            }
            mesh.UnifyNormals();
            return mesh;
        }
        #endregion
        #region basic
        public Mesh TriangleMeshFromPoints(List<Point3d> pl, int t1, int t2)
        {
            Mesh mesh = new Mesh();
            if (t1 < 1) return mesh;
            if (t2 <= t1) return mesh;
            int n = ((t1 + t2) * (t2 - t1 + 1)) / 2;
            if (n > pl.Count) return mesh;
            mesh.Vertices.AddVertices(pl);
            List<int> layer1; List<int> layer2 = new List<int>();
            for (int i = 0; i < t1; i++)
            {
                layer2.Add(i);
            }
            for (int i = t1 - 1; i < t2; i++)
            {
                layer1 = new List<int>(layer2);
                for (int j = 0; j < layer2.Count; j++)
                {
                    layer2[j] += i + 1;
                }
                layer2.Add(layer2[layer2.Count - 1] + 1);

                if (layer1.Count > 1)
                {
                    for (int j = 0; j < layer1.Count - 1; j++)
                    {
                        mesh.Faces.AddFace(layer1[j], layer1[j + 1], layer2[j + 1]);
                    }
                }
                for (int j = 0; j < layer1.Count; j++)
                {
                    mesh.Faces.AddFace(layer2[j], layer1[j], layer2[j + 1]);
                }
            }
            mesh.Compact();
            mesh.Normals.ComputeNormals();
            return mesh;
        }
        public Mesh TriangleMeshFromPoints(List<Point3d> pl)
        {
            //triangle MeshTopo Points From topo of the pyramid to the base
            double t = pl.Count;
            double l = Math.Sqrt(t * 8 + 1) - 1;
            l /= 2;
            return TriangleMeshFromPoints(pl, (int)l);
        }
        public Mesh TriangleMeshFromPoints(List<Point3d> pl, int t)
        {
            //triangle MeshTopo Points From topo of the pyramid to the base
            Mesh mesh = new Mesh();
            if (t < 2) return mesh;
            int n = ((1 + t) * t) / 2;
            if (n > pl.Count) return mesh;

            mesh.Vertices.AddVertices(pl);
            List<int> layer1; List<int> layer2 = new List<int>();
            layer2.Add(0);
            for (int i = 0; i < t - 1; i++)
            {
                layer1 = new List<int>(layer2);
                for (int j = 0; j < layer2.Count; j++)
                {
                    layer2[j] += i + 1;
                }
                layer2.Add(layer2[layer2.Count - 1] + 1);


                if (layer1.Count > 1)
                {
                    for (int j = 0; j < layer1.Count - 1; j++)
                    {
                        mesh.Faces.AddFace(layer1[j], layer1[j + 1], layer2[j + 1]);
                    }
                }
                for (int j = 0; j < layer1.Count; j++)
                {
                    mesh.Faces.AddFace(layer2[j], layer1[j], layer2[j + 1]);
                }
            }
            mesh.Compact();
            mesh.Normals.ComputeNormals();
            return mesh;
        }
        public Mesh MeshFromPoints(List<Point3d> pl, int u, int v)
        {
            if (u * v > pl.Count || u < 2 || v < 2) return null;
            Mesh mesh = new Mesh();
            for (int i = 0; i < pl.Count; i++)
            {
                mesh.Vertices.Add(pl[i]);
            }
            for (int i = 1; i < u; i++)
            {
                for (int j = 1; j < v; j++)
                {
                    mesh.Faces.AddFace(new MeshFace(
                    (j - 1) * u + i - 1,
                    (j - 1) * u + i,
                    (j) * u + i,
                    (j) * u + i - 1));
                }
            }
            mesh.Normals.ComputeNormals();
            return mesh;
        }
        public Mesh MeshFromPoints(Point3d p1, Point3d p2, Point3d p3, Point3d p4)
        {
            Mesh mesh = new Mesh();
            mesh.Vertices.Add(p1);
            mesh.Vertices.Add(p2);
            mesh.Vertices.Add(p3);
            mesh.Vertices.Add(p4);
            mesh.Faces.AddFace(0, 1, 2, 3);
            mesh.Normals.ComputeNormals();
            return mesh;
        }
        public Mesh MeshFromPoints(Point3d p1, Point3d p2, Point3d p3)
        {
            Mesh mesh = new Mesh();
            mesh.Vertices.Add(p1);
            mesh.Vertices.Add(p2);
            mesh.Vertices.Add(p3);
            mesh.Faces.AddFace(0, 1, 2);
            mesh.Normals.ComputeNormals();
            return mesh;
        }
        #endregion
    }
}