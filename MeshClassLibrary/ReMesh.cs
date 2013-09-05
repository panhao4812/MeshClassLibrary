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
    public class ReMesh
    {
        public ReMesh() { }
        /// RemeshFunctions
        public Mesh MeshFromClosedPoly(List<Polyline> x)
        {
            Mesh mesh = new Mesh();
            for (int i = 0; i < x.Count; i++)
            {
                if (x[i].Count == 4)
                {
                    int n = mesh.Vertices.Count;
                    mesh.Vertices.Add(x[i][0]);
                    mesh.Vertices.Add(x[i][1]);
                    mesh.Vertices.Add(x[i][2]);
                    mesh.Faces.AddFace(new MeshFace(n, n + 1, n + 2));
                }
                else if (x[i].Count == 5)
                {
                    int n = mesh.Vertices.Count;
                    mesh.Vertices.Add(x[i][0]);
                    mesh.Vertices.Add(x[i][1]);
                    mesh.Vertices.Add(x[i][2]);
                    mesh.Vertices.Add(x[i][3]);
                    mesh.Faces.AddFace(new MeshFace(n, n + 1, n + 2, n + 3));
                }
            }
            mesh.Normals.ComputeNormals();
            return mesh;
        }
        public List<Polyline> Remesh(List<Line> x, Vector3d y)
        {
            List<M_Point> PointList = remesh_preparations(x);
            for (int i = 0; i < PointList.Count; i++)
            {
                PointList[i].computeNormal(y);
                PointList[i].lay();
            }
            return remesh_subsequent(PointList);
        }
        public List<Polyline> Remesh(List<Line> x, Surface y)
        {
            List<M_Point> PointList = remesh_preparations(x);
            for (int i = 0; i < PointList.Count; i++)
            {
                PointList[i].computeNormal(y);
                PointList[i].lay();
            }
            return remesh_subsequent(PointList);
        }
        public List<Polyline> Remesh(List<Line> x, Mesh y)
        {
            List<M_Point> PointList = remesh_preparations(x);
            for (int i = 0; i < PointList.Count; i++)
            {
                PointList[i].computeNormal(y);
                PointList[i].lay();
            }
            return remesh_subsequent(PointList);
        }
        public List<Polyline> Remesh(List<Line> x)
        {
            List<M_Point> PointList = remesh_preparations(x);
            for (int i = 0; i < PointList.Count; i++)
            {
                PointList[i].computeNormal();
                PointList[i].lay();
            }
            return remesh_subsequent(PointList);
        }
        public List<Polyline> Remesh(List<Line> x, List<Vector3d> y)
        {
            List<M_Point> PointList = remesh_preparations(x);
            for (int i = 0; i < PointList.Count; i++)
            {
                int iy = i;
                if (iy > y.Count - 1) iy = y.Count - 1;
                PointList[i].computeNormal(y[iy]);
                PointList[i].lay();
            }
            return remesh_subsequent(PointList);
        }
        private List<M_Point> remesh_preparations(List<Line> x)
        {
            List<M_Point> PointList = new List<M_Point>();
            PointList.Add(new M_Point(x[0].From));
            PointList.Add(new M_Point(x[0].To));
            PointList[0].refpoints.Add(PointList[1]);
            PointList[1].refpoints.Add(PointList[0]);
            for (int i = 1; i < x.Count; i++)
            {
                bool sign1 = true;
                bool sign2 = true;
                int C1 = -1; int C2 = -1;
                M_Point P1 = new M_Point(x[i].From);
                M_Point P2 = new M_Point(x[i].To);
                for (int j = 0; j < PointList.Count; j++)
                {
                    if (sign1)
                    {
                        if (PointList[j].isDump(P1)) { sign1 = false; C1 = j; }
                    }
                    if (sign2)
                    {
                        if (PointList[j].isDump(P2)) { sign2 = false; C2 = j; }
                    }
                    if (sign1 == false && sign2 == false) break;
                }
                if (sign1) { PointList.Add(P1); C1 = PointList.Count - 1; }
                if (sign2) { PointList.Add(P2); C2 = PointList.Count - 1; }
                PointList[C1].refpoints.Add(PointList[C2]);
                PointList[C2].refpoints.Add(PointList[C1]);
            }
            return PointList;
        }
        private List<Polyline> remesh_subsequent(List<M_Point> PointList)
        {
            List<Polyline> pls = new List<Polyline>();
            for (int i = 0; i < PointList.Count; i++)
            {
                for (int j = 0; j < PointList[i].children.Count; j++)
                {
                    M_Point Pt = PointList[i].children[j];
                    if (Pt.order == 0)
                    {
                        Polyline pl = new Polyline();
                        int or = 0;
                        pl.Add(Pt.pos); Pt.order = 1;
                        for (int ii = 0; ii < PointList.Count; ii++)
                        {
                            bool signii = true;
                            M_Point pt2 = Pt.refpoints[or]; if (or == 0) { or = 1; } else { or = 0; }
                            for (int jj = 0; jj < pt2.children.Count; jj++)
                            {
                                if (pt2.children[jj].order == 0)
                                {
                                    if (pt2.children[jj].refpoints[or].isDump(Pt))
                                    {
                                        Pt = pt2.children[jj];
                                        pl.Add(Pt.pos); Pt.order = 1;
                                        if (or == 0) { or = 1; } else { or = 0; }
                                        signii = false;
                                        break;
                                    }
                                }
                            }
                            if (signii) { break; }
                        }
                        if (pl.Count > 2) { pl.Add(pl[0]); pls.Add(pl); }
                    }
                }
            }
            return pls;
        }

    }
}
