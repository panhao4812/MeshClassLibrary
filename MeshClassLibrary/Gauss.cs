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
   public class Gauss
    {
        public Gauss(Mesh input_mesh)
        {
            mesh = input_mesh;
            mesh.Compact();
            mesh.UnifyNormals();

            el = mesh.TopologyEdges;
            vs = mesh.TopologyVertices;
            ps = new List<VertexProperties>();

            for (int i = 0; i < vs.Count; i++)
            {
                ps.Add(new VertexProperties
                  (mesh.Normals[vs.MeshVertexIndices(i)[0]]));
                // outputs2.Add(new Vector3d());
            }
        }
        public void caculate(out List<double> v1, out List<double> v2, out List<double> v3)
        {
            CaculateAm();
            CaculateK();
            v1 = new List<double>();
            v2 = new List<double>();
            v3 = new List<double>();
            for (int i = 0; i < mesh.Vertices.Count; i++)
            {
                v1.Add(0); v2.Add(0); v3.Add(0);
            }
            for (int i = 0; i < vs.Count; i++)
            {
                double sq = ps[i].KH * ps[i].KH - ps[i].KG;
                if (sq > 0)
                {
                    sq = Math.Sqrt(sq);
                }
                else { sq = 0; };

                int[] indexV = vs.MeshVertexIndices(i);
                for (int j = 0; j < indexV.Length; j++)
                {
                    v1[indexV[j]] = ps[i].KH + sq;
                    v2[indexV[j]] = ps[i].KH - sq;
                    v3[indexV[j]] = ps[i].KG;
                }
            }
        }

        Mesh mesh;
        List<VertexProperties> ps;
        Rhino.Geometry.Collections.MeshTopologyEdgeList el;
        Rhino.Geometry.Collections.MeshTopologyVertexList vs;
        List<Point3d> outputs1 = new List<Point3d>();
        List<Vector3d> outputs2 = new List<Vector3d>();
        List<Line> outputs3 = new List<Line>();

        #region Members
        class VertexProperties
        {
            public double v1 = 0;
            public double v2 = 0;
            public double Am = 0;
            public double KG = 0;
            public double KH = 0;
            public Vector3d n = new Vector3d();
            public VertexProperties(Vector3d _n)
            {
                _n.Unitize();
                this.n = _n;
            }
        }
        public void CaculateAm()
        {
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                int[] f = vs.IndicesFromFace(i);
                Point3d p1 = vs[f[0]];
                Point3d p2 = vs[f[1]];
                Point3d p3 = vs[f[2]];
                Circle circle = new Circle(p1, p2, p3);
                double a = p1.DistanceTo(p2);
                double b = p1.DistanceTo(p3);
                double c = p2.DistanceTo(p3);
                // Print(f.Length.ToString());
                // Print(a.ToString() + "/" + b.ToString() + "/" + c.ToString());
                Point3d ci = new Point3d();
                int sign = 0;
                if (c >= a && c >= b)
                {
                    if ((c * c) < (a * a + b * b))
                    {
                        ci = circle.Center;
                    }
                    else { ci = (p2 + p3) / 2; sign = 1; }
                }
                else if (a >= c && a >= b)
                {
                    if ((a * a) < (c * c + b * b))
                    {
                        ci = circle.Center;
                    }
                    else { ci = (p2 + p1) / 2; sign = 2; }
                }
                else if (b >= a && b >= c)
                {
                    if ((b * b) < (a * a + c * c))
                    {
                        ci = circle.Center;
                    }
                    else { ci = (p1 + p3) / 2; sign = 3; }
                }
                else
                { //Print("error");
                }

                Point3d p1p2 = (p1 + p2) / 2;
                Point3d p1p3 = (p1 + p3) / 2;
                Point3d p2p3 = (p3 + p2) / 2;

                if (sign == 0)
                {
                    ps[f[0]].Am += areaTri(ci, p1, p1p2) + areaTri(ci, p1, p1p3);
                    ps[f[1]].Am += areaTri(ci, p2, p1p2) + areaTri(ci, p2, p2p3);
                    ps[f[2]].Am += areaTri(ci, p3, p1p3) + areaTri(ci, p3, p2p3);
                }
                else if (sign == 1)
                {
                    ps[f[1]].Am += areaTri(ci, p2, p1p2);
                    ps[f[2]].Am += areaTri(ci, p3, p1p3);
                    ps[f[0]].Am += areaTri(ci, p1, p1p2) + areaTri(ci, p1, p1p3);
                }
                else if (sign == 2)
                {
                    ps[f[1]].Am += areaTri(ci, p2, p2p3);
                    ps[f[0]].Am += areaTri(ci, p1, p1p3);
                    ps[f[2]].Am += areaTri(ci, p3, p1p3) + areaTri(ci, p3, p2p3);
                }
                else if (sign == 3)
                {
                    ps[f[0]].Am += areaTri(ci, p1, p1p2);
                    ps[f[2]].Am += areaTri(ci, p3, p2p3);
                    ps[f[1]].Am += areaTri(ci, p2, p1p2) + areaTri(ci, p2, p2p3);
                }
                else
                {// Print("error");
                }
                //////////////////
                ps[f[0]].KG += Vector3d.VectorAngle(p2 - p1, p3 - p1);
                ps[f[1]].KG += Vector3d.VectorAngle(p1 - p2, p3 - p2);
                ps[f[2]].KG += Vector3d.VectorAngle(p2 - p3, p1 - p3);
                /////////////////
            }
            for (int i = 0; i < el.Count; i++)
            {
                int[] f = el.GetConnectedFaces(i);
                if (f.Length == 2)
                {
                    ///////////////
                    int pi1 = el.GetTopologyVertices(i).I;
                    int pi2 = el.GetTopologyVertices(i).J;
                    int pf1 = 0; int pf2 = 0;

                    int[] vi1 = vs.IndicesFromFace(f[0]);
                    for (int j = 0; j < 3; j++)
                    {
                        if (vi1[j] != pi1 && vi1[j] != pi2)
                        {
                            pf1 = vi1[j]; break;
                        }
                    }
                    int[] vi2 = vs.IndicesFromFace(f[1]);
                    for (int j = 0; j < 3; j++)
                    {
                        if (vi2[j] != pi1 && vi2[j] != pi2)
                        {
                            pf2 = vi2[j]; break;
                        }
                    }
                    double ang1 = Vector3d.VectorAngle(vs[pi1] - vs[pf1], vs[pi2] - vs[pf1]);
                    double ang2 = Vector3d.VectorAngle(vs[pi1] - vs[pf2], vs[pi2] - vs[pf2]);
                    if (ang1 == Math.PI / 2) { ang1 = 0; } else { ang1 = 1 / Math.Tan(ang1); }
                    if (ang2 == Math.PI / 2) { ang2 = 0; } else { ang2 = 1 / Math.Tan(ang2); }
                    double total = ang1 + ang2;
                    double t1 = Vector3d.Multiply(ps[pi1].n, (vs[pi1] - vs[pi2]));
                    double t2 = Vector3d.Multiply(ps[pi2].n, (vs[pi2] - vs[pi1]));
                    ps[pi1].KH += t1 * total;
                    ps[pi2].KH += t2 * total;
                    ////////////////
                }
            }

        }
        public void CaculateK()
        {
            for (int i = 0; i < vs.Count; i++)
            {
                //  Print(ps[i].KG.ToString());
                ps[i].KG = (Math.PI * 2 - ps[i].KG) / ps[i].Am;
                ps[i].KH = ps[i].KH / (ps[i].Am * 4);
            }
        }
        #endregion

        public double areaTri(Point3d p1, Point3d p2, Point3d p3)
        {
            double a = p1.DistanceTo(p2);
            double b = p1.DistanceTo(p3);
            double c = p2.DistanceTo(p3);
            double p = (a + b + c) / 2;
            return Math.Sqrt(p * (p - a) * (p - b) * (p - c));
        }


    }
}
