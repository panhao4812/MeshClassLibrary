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
    public class KdTree3D
    {
        public bool die = false;
        //1 means X 2 means Y 3 meansZ
        public List<Point3d> points;
        public BoundingBox Counter;
        public KdTree3D(List<Point3d> Points, BoundingBox box)
        {
            points = Points; Counter = box;
        }
        public KdTree3D(List<Point3d> Points)
        {
            this.points = new List<Point3d>();
            double x1, x2, y1, y2, z1, z2;
            Points.Sort(CompareDinos_X);
            x1 = Points[0].X; x2 = Points[Points.Count - 1].X;
            Points.Sort(CompareDinos_Y);
            y1 = Points[0].Y; y2 = Points[Points.Count - 1].Y;
            Points.Sort(CompareDinos_Z);
            z1 = Points[0].Z; z2 = Points[Points.Count - 1].Z;

            for (int i = 0; i < Points.Count; i++)
            {
                Point3d p = Points[i];
                if (p.X != x1 && p.X != x2 && p.Y != y1 && p.Y != y2 && p.Z != z1 && p.Z != z2)
                {
                    this.points.Add(p);
                }
            }
            this.Counter = new BoundingBox(x1, y1, z1, x2, y2, z2);
        }
        public List<KdTree3D> Split()
        {
            if (this.points.Count < 1) { return null; }
            List<KdTree3D> trees = new List<KdTree3D>();
            this.die = true;
            BoundingBox box1 = new BoundingBox();
            BoundingBox box2 = new BoundingBox();
            List<Point3d> ps1 = new List<Point3d>();
            List<Point3d> ps2 = new List<Point3d>();

            double x1 = Counter.Min.X;
            double y1 = Counter.Min.Y;
            double z1 = Counter.Min.Z;
            double x2 = Counter.Max.X;
            double y2 = Counter.Max.Y;
            double z2 = Counter.Max.Z;
            double t1 = x2 - x1;
            double t2 = y2 - y1;
            double t3 = z2 - z1;

            if (t1 >= t2 && t1 >= t3)
            {
                this.points.Sort(CompareDinos_X);
                int count = cut(this.points.Count);
                double t4 = this.points[count].X;
                box1 = new BoundingBox(x1, y1, z1, t4, y2, z2);
                box2 = new BoundingBox(t4, y1, z1, x2, y2, z2);
            }
            else if (t2 >= t1 && t2 >= t3)
            {
                this.points.Sort(CompareDinos_Y);
                int count = cut(this.points.Count);
                double t4 = this.points[count].Y;
                box1 = new BoundingBox(x1, y1, z1, x2, t4, z2);
                box2 = new BoundingBox(x1, t4, z1, x2, y2, z2);
            }
            else if (t3 >= t2 && t3 >= t1)
            {
                this.points.Sort(CompareDinos_Z);
                int count = cut(this.points.Count);
                double t4 = this.points[count].Z;
                box1 = new BoundingBox(x1, y1, z1, x2, y2, t4);
                box2 = new BoundingBox(x1, y1, t4, x2, y2, z2);
            }

            for (int i = 0; i < this.points.Count; i++)
            {
                Point3d p = this.points[i];
                if (isPointin(p, box1)) ps1.Add(p);
                if (isPointin(p, box2)) ps2.Add(p);
            }
            trees.Add(new KdTree3D(ps1, box1));
            trees.Add(new KdTree3D(ps2, box2));
            return trees;
        }
        private static int CompareDinos_X(Point3d x, Point3d y)
        {
            if (x == null) { if (y == null) { return 0; } else { return -1; } }
            else
            {
                if (y == null) { return 1; }
                else
                {
                    if (x.X > y.X) return 1;
                    if (x.X == y.X) return 0;
                    if (x.X < y.X) return -1;
                    else return 0;
                }
            }
        }
        private static int CompareDinos_Y(Point3d x, Point3d y)
        {
            if (x == null) { if (y == null) { return 0; } else { return -1; } }
            else
            {
                if (y == null) { return 1; }
                else
                {
                    if (x.Y > y.Y) return 1;
                    if (x.Y == y.Y) return 0;
                    if (x.Y < y.Y) return -1;
                    else return 0;
                }
            }
        }
        private static int CompareDinos_Z(Point3d x, Point3d y)
        {
            if (x == null) { if (y == null) { return 0; } else { return -1; } }
            else
            {
                if (y == null) { return 1; }
                else
                {
                    if (x.Z > y.Z) return 1;
                    if (x.Z == y.Z) return 0;
                    if (x.Z < y.Z) return -1;
                    else return 0;
                }
            }
        }
        public int cut(double i)
        {
            if (i < 1) { return -1; }
            return (int)Math.Floor(i / 2);
        }
        public bool isPointin(Point3d p, BoundingBox pl)
        {
            double x1 = pl.Min.X;
            double y1 = pl.Min.Y;
            double z1 = pl.Min.Z;
            double x2 = pl.Max.X;
            double y2 = pl.Max.Y;
            double z2 = pl.Max.Z;
            if (p.X > x1 && p.X < x2 &&
              p.Y > y1 && p.Y < y2 &&
              p.Z > z1 && p.Z < z2
            ) { return true; }
            else { return false; }
        }
        public static bool isPointinBox(Point3d p, BoundingBox pl)
        {
            double x1 = pl.Min.X;
            double y1 = pl.Min.Y;
            double z1 = pl.Min.Z;
            double x2 = pl.Max.X;
            double y2 = pl.Max.Y;
            double z2 = pl.Max.Z;
            if (p.X > x1 && p.X < x2 &&
              p.Y > y1 && p.Y < y2 &&
              p.Z > z1 && p.Z < z2
            ) { return true; }
            else { return false; }
        }
        public static List<BoundingBox> SolveKdTree3D(List<Point3d> points, BoundingBox box)
        {

            List<BoundingBox> PL = new List<BoundingBox>();
            if (points.Count < 1) return PL;
            List<KdTree3D> trees = new List<KdTree3D>();

            List<Point3d> pts = new List<Point3d>();
            points.ForEach(delegate(Point3d p)
            {
                if (isPointinBox(p, box)) pts.Add(p);
            });
            trees.Add(new KdTree3D(pts, box));

            bool toggle = true;
            for (int i = 0; i < points.Count; i++)
            {
                if (toggle == false) break;
                toggle = false;
                for (int j = 0; j < trees.Count; j++)
                {
                    if (trees[j].die == false && trees[j].points.Count > 0)
                    {
                        trees.AddRange(trees[j].Split()); toggle = true;
                    }
                }
            }
            for (int i = 0; i < trees.Count; i++)
            { if (trees[i].die == false)PL.Add(trees[i].Counter); }
            return PL;
        }
        public static List<BoundingBox> SolveKdTree3D(List<Point3d> points)
        {
            List<BoundingBox> PL = new List<BoundingBox>();
            if (points.Count < 1) return PL;
            List<KdTree3D> trees = new List<KdTree3D>();
            trees.Add(new KdTree3D(points));
            bool toggle = true;
            for (int i = 0; i < points.Count; i++)
            {
                if (toggle == false) break;
                toggle = false;
                for (int j = 0; j < trees.Count; j++)
                {
                    if (trees[j].die == false && trees[j].points.Count > 0)
                    {
                        trees.AddRange(trees[j].Split()); toggle = true;
                    }
                }
            }
            for (int i = 0; i < trees.Count; i++)
            { if (trees[i].die == false)PL.Add(trees[i].Counter); }
            return PL;
        }
        public void Dispose()
        {
            this.points = default(List<Point3d>);
            this.Counter = default(BoundingBox);
        }
    }

    public class KdTree
    {
        public bool die = false;
        public int side = 1;
        //1 means X 2 means Y 3 meansZ
        public List<Point3d> points;
        public Polyline Counter;
        public KdTree(List<Point3d> Points, Polyline polyline, int Side)
        {
            points = Points; Counter = polyline; side = Side;
        }
        public KdTree(List<Point3d> Points)
        {
            this.points = new List<Point3d>();
            double x1, x2, y1, y2, z1, z2;
            Points.Sort(CompareDinos_X);
            x1 = Points[0].X; x2 = Points[Points.Count - 1].X;
            Points.Sort(CompareDinos_Y);
            y1 = Points[0].Y; y2 = Points[Points.Count - 1].Y;
            Points.Sort(CompareDinos_Z);
            z1 = Points[0].Z; z2 = Points[Points.Count - 1].Z;
            if ((x2 - x1) > (y2 - y1)) { side = 1; } else { side = 2; }
            for (int i = 0; i < Points.Count; i++)
            {
                Point3d p = Points[i];
                if (p.X != x1 && p.X != x2 && p.Y != y1 && p.Y != y2)
                {
                    this.points.Add(p);
                }
            }
            Point3d p1, p2, p3, p4;
            p1 = new Point3d(x1, y1, 0);
            p2 = new Point3d(x2, y1, 0);
            p3 = new Point3d(x2, y2, 0);
            p4 = new Point3d(x1, y2, 0);
            this.Counter = new Polyline();
            this.Counter.Add(p1);
            this.Counter.Add(p2);
            this.Counter.Add(p3);
            this.Counter.Add(p4);
            this.Counter.Add(p1);
        }
        public List<KdTree> Split()
        {
            if (this.points.Count < 1) { return null; }
            List<KdTree> trees = new List<KdTree>();
            this.die = true;
            Point3d p1, p2, p3, p4;
            p1 = this.Counter[0];
            p2 = this.Counter[1];
            p3 = this.Counter[2];
            p4 = this.Counter[3];
            Polyline pl1 = new Polyline();
            Polyline pl2 = new Polyline();
            List<Point3d> ps1 = new List<Point3d>();
            List<Point3d> ps2 = new List<Point3d>();
            int side2 = 1;
            if (side == 1)
            {
                this.points.Sort(CompareDinos_X); side2 = 2;
                int count = cut(this.points.Count);
                Point3d p5 = new Point3d(this.points[count].X, p1.Y, 0);
                Point3d p6 = new Point3d(this.points[count].X, p3.Y, 0);
                pl1.Add(p1); pl1.Add(p5); pl1.Add(p6); pl1.Add(p4); pl1.Add(p1);
                pl2.Add(p5); pl2.Add(p2); pl2.Add(p3); pl2.Add(p6); pl2.Add(p5);
            }
            else if (side == 2)
            {
                this.points.Sort(CompareDinos_Y);
                int count = cut(this.points.Count);
                Point3d p5 = new Point3d(p1.X, this.points[count].Y, 0);
                Point3d p6 = new Point3d(p2.X, this.points[count].Y, 0);
                pl1.Add(p1); pl1.Add(p2); pl1.Add(p6); pl1.Add(p5); pl1.Add(p1);
                pl2.Add(p5); pl2.Add(p6); pl2.Add(p3); pl2.Add(p4); pl2.Add(p5);
            }
            //   if(side==3){this.points.Sort(CompareDinos_Z);}
            for (int i = 0; i < this.points.Count; i++)
            {
                Point3d p = this.points[i];
                if (isPointin(p, pl1)) ps1.Add(p);
                if (isPointin(p, pl2)) ps2.Add(p);
            }
            trees.Add(new KdTree(ps1, pl1, side2));
            trees.Add(new KdTree(ps2, pl2, side2));
            return trees;
        }
        public List<KdTree> Split2()
        {
            if (this.points.Count < 1) { return null; }
            List<KdTree> trees = new List<KdTree>();
            this.die = true;
            Point3d p1, p2, p3, p4;
            p1 = this.Counter[0];
            p2 = this.Counter[1];
            p3 = this.Counter[2];
            p4 = this.Counter[3];
            Polyline pl1 = new Polyline();
            Polyline pl2 = new Polyline();
            List<Point3d> ps1 = new List<Point3d>();
            List<Point3d> ps2 = new List<Point3d>();
            int side2 = 1;

            if ((p2.X - p1.X) > (p3.Y - p2.Y)) { side = 1; } else { side = 2; }


            if (side == 1)
            {
                this.points.Sort(CompareDinos_X);
                int count = cut(this.points.Count);
                Point3d p5 = new Point3d(this.points[count].X, p1.Y, 0);
                Point3d p6 = new Point3d(this.points[count].X, p3.Y, 0);
                pl1.Add(p1); pl1.Add(p5); pl1.Add(p6); pl1.Add(p4); pl1.Add(p1);
                pl2.Add(p5); pl2.Add(p2); pl2.Add(p3); pl2.Add(p6); pl2.Add(p5);
            }
            else if (side == 2)
            {
                this.points.Sort(CompareDinos_Y);
                int count = cut(this.points.Count);
                Point3d p5 = new Point3d(p1.X, this.points[count].Y, 0);
                Point3d p6 = new Point3d(p2.X, this.points[count].Y, 0);
                pl1.Add(p1); pl1.Add(p2); pl1.Add(p6); pl1.Add(p5); pl1.Add(p1);
                pl2.Add(p5); pl2.Add(p6); pl2.Add(p3); pl2.Add(p4); pl2.Add(p5);
            }
            //   if(side==3){this.points.Sort(CompareDinos_Z);}
            for (int i = 0; i < this.points.Count; i++)
            {
                Point3d p = this.points[i];
                if (isPointin(p, pl1)) ps1.Add(p);
                if (isPointin(p, pl2)) ps2.Add(p);
            }
            trees.Add(new KdTree(ps1, pl1, side2));
            trees.Add(new KdTree(ps2, pl2, side2));
            return trees;
        }
        private static int CompareDinos_X(Point3d x, Point3d y)
        {
            if (x == null) { if (y == null) { return 0; } else { return -1; } }
            else
            {
                if (y == null) { return 1; }
                else
                {
                    if (x.X > y.X) return 1;
                    if (x.X == y.X) return 0;
                    if (x.X < y.X) return -1;
                    else return 0;
                }
            }
        }
        private static int CompareDinos_Y(Point3d x, Point3d y)
        {
            if (x == null) { if (y == null) { return 0; } else { return -1; } }
            else
            {
                if (y == null) { return 1; }
                else
                {
                    if (x.Y > y.Y) return 1;
                    if (x.Y == y.Y) return 0;
                    if (x.Y < y.Y) return -1;
                    else return 0;
                }
            }
        }
        private static int CompareDinos_Z(Point3d x, Point3d y)
        {
            if (x == null) { if (y == null) { return 0; } else { return -1; } }
            else
            {
                if (y == null) { return 1; }
                else
                {
                    if (x.Z > y.Z) return 1;
                    if (x.Z == y.Z) return 0;
                    if (x.Z < y.Z) return -1;
                    else return 0;
                }
            }
        }
        public int cut(double i)
        {
            if (i < 1) { return -1; }
            return (int)Math.Floor(i / 2);
        }
        public bool isPointin(Point3d p, Polyline pl)
        {
            double x1 = pl[0].X; double y1 = pl[0].Y;
            double x2 = pl[2].X; double y2 = pl[2].Y;
            if (p.X > x1 && p.X < x2 && p.Y > y1 && p.Y < y2) { return true; }
            else { return false; }
        }
        public static List<Polyline> SolveKdTree(List<Point3d> points)
        {
            List<Polyline> PL = new List<Polyline>();
            if (points.Count < 1) return PL;
            List<KdTree> trees = new List<KdTree>();
            trees.Add(new KdTree(points));
            bool toggle = true;
            for (int i = 0; i < points.Count; i++)
            {
                if (toggle == false) break;
                toggle = false;
                for (int j = 0; j < trees.Count; j++)
                {
                    if (trees[j].die == false && trees[j].points.Count > 0)
                    {
                        trees.AddRange(trees[j].Split()); toggle = true;
                    }
                }
            }

            for (int i = 0; i < trees.Count; i++)
            {
                if (trees[i].die == false)
                {
                    PL.Add(trees[i].Counter);
                }
            }
            return PL;
        }
        public static List<Polyline> SolveKdTree2(List<Point3d> points)
        {
            List<Polyline> PL = new List<Polyline>();
            if (points.Count < 1) return PL;
            List<KdTree> trees = new List<KdTree>();
            trees.Add(new KdTree(points));
            bool toggle = true;
            for (int i = 0; i < points.Count; i++)
            {
                if (toggle == false) break;
                toggle = false;
                for (int j = 0; j < trees.Count; j++)
                {
                    if (trees[j].die == false && trees[j].points.Count > 0)
                    {
                        trees.AddRange(trees[j].Split2()); toggle = true;
                    }
                }
            }

            for (int i = 0; i < trees.Count; i++)
            {
                if (trees[i].die == false)
                {
                    PL.Add(trees[i].Counter);
                }
            }
            return PL;
        }
        public static List<Polyline> SolveKdTree2(List<Point3d> points, Polyline pl)
        {
            List<Point3d> points2 = new List<Point3d>();
            List<Point3d> Points = new List<Point3d>();
            for (int i = 0; i < pl.Count; i++) { Points.Add(pl[i]); }
            int side = 1;
            double x1, x2, y1, y2;
            Points.Sort(CompareDinos_X);
            x1 = Points[0].X; x2 = Points[Points.Count - 1].X;
            Points.Sort(CompareDinos_Y);
            y1 = Points[0].Y; y2 = Points[Points.Count - 1].Y;
            if ((x2 - x1) > (y2 - y1)) { side = 1; } else { side = 2; }
            for (int i = 0; i < points.Count; i++)
            {
                Point3d p = points[i];
                if (p.X > x1 && p.X < x2 && p.Y > y1 && p.Y < y2)
                {
                    points2.Add(p);
                }
            }
            Point3d p1, p2, p3, p4;
            p1 = new Point3d(x1, y1, 0);
            p2 = new Point3d(x2, y1, 0);
            p3 = new Point3d(x2, y2, 0);
            p4 = new Point3d(x1, y2, 0);
            pl = new Polyline();
            pl.Add(p1);
            pl.Add(p2);
            pl.Add(p3);
            pl.Add(p4);
            pl.Add(p1);

            points = points2;

            List<Polyline> PL = new List<Polyline>();
            if (points.Count < 1) return PL;
            List<KdTree> trees = new List<KdTree>();
            trees.Add(new KdTree(points, pl, side));
            bool toggle = true;
            for (int i = 0; i < points.Count; i++)
            {
                if (toggle == false) break;
                toggle = false;
                for (int j = 0; j < trees.Count; j++)
                {
                    if (trees[j].die == false && trees[j].points.Count > 0)
                    {
                        trees.AddRange(trees[j].Split2()); toggle = true;
                    }
                }
            }

            for (int i = 0; i < trees.Count; i++)
            {
                if (trees[i].die == false)
                {
                    PL.Add(trees[i].Counter);
                }
            }
            return PL;
        }
        public static List<Polyline> SolveKdTree2(List<Point3d> points, Polyline pl, Plane p)
        {
            //Plane p = new Plane(pl[0], pl[1], pl[2]);
            Transform xf = Transform.PlaneToPlane(p, Plane.WorldXY);
            pl.Transform(xf);
            for (int i = 0; i < points.Count; i++)
            {
                Point3d pt = new Point3d(points[i]);
                pt.Transform(Transform.PlanarProjection(p));
                pt.Transform(xf);
                points[i] = pt;
            }
            List<Polyline> pls = KdTree.SolveKdTree2(points, pl);
            xf = Transform.PlaneToPlane(Plane.WorldXY, p);
            for (int i = 0; i < pls.Count; i++)
            {
                pls[i].Transform(xf);
            }
            return pls;
        }
    }
}
