﻿using Rhino;
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
    class MeshUnfold
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
            f.AddTransform(Transform.PlaneToPlane(new Plane(f.Center(), f.Normal), Plane.WorldXY));

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
                } return cen;
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

    }
}
