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
    public class MeshSimplify
    {
        private class MeshSimplify_Point
        {
            public Point3d pos;
            public List<MeshSimplify_Point> refpoints;
            private double tol = 0.01;
            public Vector3d N;
            public double order = 0;
            public List<MeshSimplify_Point> children = new List<MeshSimplify_Point>();
            /*for sort and lay functions
             in sort method  order means clockwise  angel
             in  lay order equals boolean */
            public MeshSimplify_Point(Point3d Pos)
            {
                refpoints = new List<MeshSimplify_Point>();
                pos = Pos;
            }
            public bool isDump(MeshSimplify_Point other)
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
            private static int CompareDinosByLength(MeshSimplify_Point x, MeshSimplify_Point y)
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
                    MeshSimplify_Point p = new MeshSimplify_Point(this.pos);
                    p.refpoints.AddRange(this.refpoints);
                    p.refpoints.Add(this);
                    children.Add(p); children.Add(p);
                }
                for (int i = 0; i < this.refpoints.Count; i++)
                {
                    MeshSimplify_Point p = new MeshSimplify_Point(this.pos);
                    int after = i + 1; if (after > refpoints.Count - 1) after = 0;
                    p.refpoints.Add(refpoints[i]);
                    p.refpoints.Add(refpoints[after]);
                    children.Add(p);
                }

            }
        }
        private class MeshSimplify_Face
        {
            public Rhino.Geometry.MeshFace face = new MeshFace();
            public List<MeshSimplify_Face> reffaces = new List<MeshSimplify_Face>();
            public double order = 0;
            public MeshSimplify_Face() { }
            public MeshSimplify_Face(MeshFace f)
            {
                this.face = f;
            }
            public MeshSimplify_Face(int a, int b, int c, int d)
            {
                this.face = new MeshFace(a, b, c, d);
            }
            public MeshSimplify_Face(int a, int b, int c)
            {
                this.face = new MeshFace(a, b, c);
            }
            public void Add(MeshSimplify_Face f)
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
            public static List<List<MeshSimplify_Face>> Group(List<MeshSimplify_Face> faces)
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
                List<List<MeshSimplify_Face>> output = new List<List<MeshSimplify_Face>>();
                for (int i = 1; i < level; i++)
                {
                    List<MeshSimplify_Face> out1 = new List<MeshSimplify_Face>();
                    output.Add(out1);
                }

                for (int i = 0; i < faces.Count; i++)
                {
                    output[(int)faces[i].order - 1].Add(faces[i]);
                }
                return output;
            }
            public static bool FindNext(MeshSimplify_Face face, int level)
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
        public MeshSimplify() { }
        private List<MeshSimplify_Point> preDivide(Mesh mesh)
        {
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            List<MeshSimplify_Point> PointList = new List<MeshSimplify_Point>();
            for (int i = 0; i < vs.Count; i++)
            {
                MeshSimplify_Point pt = new MeshSimplify_Point(vs[i]);
                if (vs.MeshVertexIndices(i).Length > 0)
                {
                    pt.N = mesh.Normals[vs.MeshVertexIndices(i)[0]];
                }
                else
                {
                    pt.computeNormal(mesh);
                }
                PointList.Add(pt);
            }
            for (int i = 0; i < vs.Count; i++)
            {
                int[] index = vs.ConnectedTopologyVertices(i);
                for (int j = 0; j < index.Length; j++)
                {
                    PointList[i].refpoints.Add(PointList[index[j]]);
                }
                PointList[i].Sort();
            }

            /////////////////////////////////////////////////////
            for (int i = 0; i < vs.Count; i++)
            {
                PointList[i].order = 0;
            }
            for (int i = 0; i < el.Count; i++)
            {
                if (el.GetConnectedFaces(i).Length == 1)
                {
                    PointList[el.GetTopologyVertices(i).I].order = 5;
                    PointList[el.GetTopologyVertices(i).J].order = 5;
                }
            }
            for (int i = 0; i < vs.Count; i++)
            {
                if (PointList[i].order == 5)
                {
                    if (PointList[i].refpoints.Count != 3)
                    {
                        PointList[i].order = 4;
                    }
                }
                else
                {
                    if (PointList[i].refpoints.Count != 4) PointList[i].order = 4;
                }
            }
            //////////////////////////////////////////////////////////
            for (int k = 0; k < PointList.Count; k++)
            {
                bool sign = true;
                for (int i = 0; i < PointList.Count; i++)
                {
                    if (PointList[i].order == 4)
                    {
                        sign = false;
                        PointList[i].order++;
                        if (PointList[i].refpoints.Count == 4)
                        {
                            if (PointList[i].refpoints[0].order == 5 && PointList[i].refpoints[2].order != 5) { PointList[i].refpoints[2].order = 1; }
                            if (PointList[i].refpoints[1].order == 5 && PointList[i].refpoints[3].order != 5) { PointList[i].refpoints[3].order = 1; }
                            if (PointList[i].refpoints[2].order == 5 && PointList[i].refpoints[0].order != 5) { PointList[i].refpoints[0].order = 1; }
                            if (PointList[i].refpoints[3].order == 5 && PointList[i].refpoints[1].order != 5) { PointList[i].refpoints[1].order = 1; }
                        }
                        else
                        {
                            for (int j = 0; j < PointList[i].refpoints.Count; j++)
                            {
                                PointList[i].refpoints[j].order++;
                            }
                        }
                    }
                }
                for (int i = 0; i < PointList.Count; i++)
                {
                    if (PointList[i].order > 0 && PointList[i].order < 4) PointList[i].order = 4;
                }
                if (sign) { break; }
            }
            return PointList;
        }
        public List<Line> MeshProfile(Mesh mesh)
        {
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            List<MeshSimplify_Point> PointList = preDivide(mesh);
            List<Line> output = new List<Line>();
            for (int i = 0; i < el.Count; i++)
            {
                int a = el.GetTopologyVertices(i).I;
                int b = el.GetTopologyVertices(i).J;
                if (PointList[a].order > 0 && PointList[b].order > 0)
                {
                    output.Add(el.EdgeLine(i));
                }
            }
            /*   List<Point3d> output1 = new List<Point3d>();
               List<string> output2 = new List<string>();
               for(int i = 0;i < PointList.Count;i++){
                 output1.Add(PointList[i].pos);
                 string str = PointList[i].order.ToString();
                 if(PointList[i].order == 0)str = "";
                 output2.Add(str);
               }*/
            return output;
        }
        public List<Mesh> MeshSeperate(Mesh mesh)
        {
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            List<MeshSimplify_Point> PointList = preDivide(mesh);
            List<MeshSimplify_Face> faces = new List<MeshSimplify_Face>();
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                faces.Add(new MeshSimplify_Face(mesh.Faces[i]));
            }
            for (int i = 0; i < el.Count; i++)
            {
                int a = el.GetTopologyVertices(i).I;
                int b = el.GetTopologyVertices(i).J;
                if (PointList[a].order == 0 || PointList[b].order == 0)
                {
                    int[] index = el.GetConnectedFaces(i);
                    if (index.Length == 2)
                    {
                        faces[index[0]].reffaces.Add(faces[index[1]]);
                        faces[index[1]].reffaces.Add(faces[index[0]]);
                    }
                }
            }
            List<Mesh> outputMesh = new List<Mesh>();
            List<List<MeshSimplify_Face>> outtemp = MeshSimplify_Face.Group(faces);
            for (int i = 0; i < outtemp.Count; i++)
            {
                Mesh meshout = new Mesh();
                meshout.Vertices.AddVertices(mesh.Vertices);
                for (int j = 0; j < outtemp[i].Count; j++)
                {
                    meshout.Faces.AddFace(outtemp[i][j].face);
                }
                meshout.Compact();
                outputMesh.Add(meshout);
            }
            return outputMesh;

            /*
    List<Point3d> outpos = new List<Point3d>();
    List<string> outsign = new  List<string>();
    for(int i = 0;i < faces.Count;i++){
      outpos.Add(faces[i].pos(mesh.Vertices));
      string str = i.ToString() + "/" + faces[i].reffaces.Count.ToString();
      outsign.Add(str);
    }
    B = outpos;C = outsign;
    */
        }
        public List<NurbsSurface> MeshSeperate2Nurbs(Mesh mesh)
        {
            MeshConvert conv = new MeshConvert();
            List<NurbsSurface> output = new List<NurbsSurface>();
            List<Mesh> input = MeshSeperate(mesh);

            input.ForEach(delegate(Mesh mesh1)
            {
                NurbsSurface Surf;
                conv.Mesh2Nurbs(mesh1, out Surf);
                output.Add(Surf);
            });
            return output;
        }
    }

}