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
    class MeshUnfold
    {
    public class face{
        public Point3d[] Pts;
        public int sign = -1;
        public Vector3d Normal;
        public Transform Xform = Rhino.Geometry.Transform.Identity;
        public void AddTransform(Transform xform){
            this.Xform= Rhino.Geometry.Transform.Multiply(this.Xform, xform);
        }
        public face(List<Point3d> pts)
        {
            this.Pts = pts.ToArray();
            Normal = new Plane(this.Pts[0], this.Pts[1], this.Pts[2]).Normal;
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
        public void transform(Transform xform){
            for(int i=0;i<this.Pts.Length;i++){
                this.Pts[i].Transform(xform);
            }
        }
        public Point3d Center()
        {
            Point3d cen = new Point3d();
            if (this.Pts.Length > 3)
            {
                cen += this.Pts[0];
                cen += this.Pts[1];
                cen += this.Pts[2];
                cen /= 3;
            } return cen;
        }
        public Mesh  DrawFace()
        {
            Mesh mesh = new Mesh();
            if (this.Pts.Length > 3)
            {
                mesh.Vertices.Add(this.Pts[0]);
                mesh.Vertices.Add(this.Pts[1]);
                mesh.Vertices.Add(this.Pts[2]);
                mesh.Faces.AddFace(0, 1, 2);
                mesh.Normals.ComputeNormals();
            }
            return mesh;
        }

}
    public class edge{
        public Transform fold(face f1, edge e1, face f2)
        {
            Plane p1 ,p2;
            Point3d cen = this.From;
            Vector3d v = this.To - this.From;
            p1 = new Plane(cen, v, f1.Normal);
            p2 = new Plane(cen, v, f2.Normal);
            return Transform.PlaneToPlane(p1, p2);
        }
        public List< face> Faces = new  List< face>();
        public Point3d From;
        public Point3d To;
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
            if (this.From.Equals(el.From) && this.To.Equals(el.To)) return true;
            if (this.From.Equals(el.To) && this.To.Equals(el.From)) return true;
            return false;
        }
        public bool IsValid()
        {
            return this.From.Equals(this.To);
        }
}
}
}
