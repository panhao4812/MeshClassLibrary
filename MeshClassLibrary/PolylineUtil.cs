using Rhino.Geometry;
using System.Collections.Generic;

namespace MeshClassLibrary
{
    public class PolyPlannar
    {
        public PolyPlannar() { }
        private MeshCreation mc = new MeshCreation();
        public Mesh Plannar(Polyline pl)
        {
            if (pl[0].DistanceTo(pl[pl.Count - 1]) > Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance)
            {
                pl.Add(pl[0]);
            }
            Mesh mesh = new Mesh();
            if (pl.Count > 100)
            {
                Brep[] b = Brep.CreatePlanarBreps(pl.ToNurbsCurve(),Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
                if (b.Length == 0) return mesh;
                return Mesh.CreateFromBrep(b[0],MeshingParameters.Default)[0];
            }
            return Mesh.CreateFromClosedPolyline(pl);
            //When number of edges is larger than 100 use Brep2Mesh instead
        }
        public Mesh Plannar2D(Polyline pl)
        {
            Mesh mesh = new Mesh();
            if (pl[0].DistanceTo(pl[pl.Count - 1]) < Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance)
            {
                pl.RemoveAt(pl.Count - 1);
            }
            Polyline pl2 = new Polyline(pl);
            if (pl.Count < 3) return mesh;
            if (pl.Count == 3) { mesh.Append(mc.MeshFromPoints(pl2[0], pl2[1], pl2[2])); return mesh; }
            ///////////////////////////////////
            double[] polyX = new double[pl.Count];//  =  horizontalcoordinates of corners
            double[] polyY = new double[pl.Count];  // =  verticalcoordinates of corners
            for (int ii = 0; ii < pl.Count; ii++)
            {
                polyX[ii] = pl[ii].X;
                polyY[ii] = pl[ii].Y;
            }
            ////////////////////////////////
            while (pl2.Count >= 3)
            {
                int sign = -1; int before = -1; int after = -1;
                for (int i = 0; i < pl2.Count; i++)
                {
                    if (i == 0) before = pl2.Count - 1; else before = i - 1;
                    if (i == pl2.Count - 1) after = 0; else after = i + 1;

                    Point3d cen = (pl2[before] + pl2[after]) / 2;
                    if (!LinePolyline(new Line(pl2[before], pl2[after]), pl))
                    {
                        if (pointInPolygon2D_1(polyX, polyY, cen.X, cen.Y))
                        {
                            sign = i;
                            if (sign == 0) before = pl2.Count - 1; else before = sign - 1;
                            if (sign == pl2.Count - 1) after = 0; else after = sign + 1;
                            mesh.Append(mc.MeshFromPoints(pl2[before], pl2[sign], pl2[after]));
                            pl2.RemoveAt(sign);
                            break;
                        }
                    }
                }
            }
            mesh.Weld(Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance);
            return mesh;
        }
        private bool LinePolyline(Line l, Polyline pl)
        {
            for (int i = 0; i < pl.Count; i++)
            {
                int before = -1; if (i == 0) before = pl.Count - 1; else before = i - 1;
                Line l1 = new Line(pl[i], pl[before]);
                double a, b;
                Rhino.Geometry.Intersect.Intersection.LineLine(l, l1, out a, out b, 0.001, false);
                if (a > 0 && a < 1 && b > 0 && b < 1)
                {
                    return true;
                }
            }
            return false;
        }

        bool pointInPolygon2D_1(Point3d p, Polyline pl)
        {
            int polySides = pl.Count;// how many cornersthe polygon has
            double[] polyX = new double[pl.Count];//  =  horizontalcoordinates of corners
            double[] polyY = new double[pl.Count];  // =  verticalcoordinates of corners
            double x = p.X, y = p.Y;//=  point to be tested
            for (int ii = 0; ii < pl.Count; ii++)
            {
                polyX[ii] = pl[ii].X;
                polyY[ii] = pl[ii].Y;
            }
            bool oddNodes = true;
            int i, j = polySides - 1;
            for (i = 0; i < polySides - 1; i++)
            {
                if ((polyY[i] < y && polyY[j] >= y
                  || polyY[j] < y && polyY[i] >= y)
                  && (polyX[i] <= x || polyX[j] <= x))
                {
                    oddNodes ^= (polyX[i] + (y - polyY[i]) / (polyY[j] - polyY[i]) * (polyX[j] - polyX[i]) < x);
                }
                j = i;
            }
            return !oddNodes;
        }
        bool pointInPolygon2D_1(double[] polyX, double[] polyY, double x, double y)
        {
            int polySides = polyX.Length;
            bool oddNodes = true;
            int i, j = polySides - 1;
            for (i = 0; i < polySides - 1; i++)
            {
                if ((polyY[i] < y && polyY[j] >= y
                  || polyY[j] < y && polyY[i] >= y)
                  && (polyX[i] <= x || polyX[j] <= x))
                {
                    oddNodes ^= (polyX[i] + (y - polyY[i]) / (polyY[j] - polyY[i]) * (polyX[j] - polyX[i]) < x);
                }
                j = i;
            }
            return !oddNodes;


        }
        bool pointInPolygon2D_2(double[] polyX, double[] polyY, double x, double y)
        {
            int polySides = polyX.Length;
            bool oddNodes = true;
            int i, j = polySides - 1;
            for (i = 0; i < polySides-1; i++)
            {
                if (polyY[i] < y && polyY[j] >= y
                || polyY[j] < y && polyY[i] >= y)
                {
                    if (polyX[i] + (y - polyY[i]) / (polyY[j] - polyY[i]) * (polyX[j] - polyX[i]) < x)
                    {
                        oddNodes = !oddNodes;
                    }
                }
                j = i;
            }
            return !oddNodes;
        }
        bool pointInPolygon2D_3(double[] polyX, double[] polyY, double x, double y)
        {
            int polySides = polyX.Length;
            bool oddNodes = true;
            int i, j = polySides - 1;
            for (i = 0; i < polySides-1; i++)
            {
                if ((polyY[i] < y && polyY[j] >= y
                || polyY[j] < y && polyY[i] >= y)
                && (polyX[i] <= x || polyX[j] <= x))
                {
                    if (polyX[i] + (y - polyY[i]) / (polyY[j] - polyY[i]) * (polyX[j] - polyX[i]) < x)
                    {
                        oddNodes = !oddNodes;
                    }
                }
                j = i;
            }
            return !oddNodes;
        }
        #region UsedFUnc
        public double MinimalAngle(Point3d p1, Point3d p2, Point3d p3)
        {
            double output = double.MaxValue;
            Vector3d v1 = p2 - p1; Vector3d v2 = p3 - p1;
            double t = Vector3d.VectorAngle(v1, v2);
            if (t < output) output = t;
            v1 = p1 - p2; v2 = p3 - p2;
            t = Vector3d.VectorAngle(v1, v2);
            if (t < output) output = t;
            v1 = p2 - p3; v2 = p1 - p3;
            t = Vector3d.VectorAngle(v1, v2);
            if (t < output) output = t;
            return output;
        }
        public bool isHullVertice(int t, Polyline pl)
        {
            Grasshopper.Kernel.Geometry.Node2List list2 = new Grasshopper.Kernel.Geometry.Node2List();
            List<int> hull = new List<int>();
            for (int i = 0; i < pl.Count; i++)
            {
                Grasshopper.Kernel.Geometry.Node2 node = new Grasshopper.Kernel.Geometry.Node2(pl[i].X, pl[i].Y);
                list2.Append(node);
            }
            Grasshopper.Kernel.Geometry.ConvexHull.Solver.Compute(list2, hull);
            return hull.Contains(t);
        }
        #endregion
    }
    public class PolylineSmooth
    {
        public static Polyline cc_Subdivide(Polyline ptlist)
        {
            List<Point3d> ps2 = new List<Point3d>();
            if (ptlist.Count < 3)
            {
                return ptlist;
            }
            if (ptlist[0].DistanceTo(ptlist[ptlist.Count - 1]) > 0.001)
            {
                ps2.Add(ptlist[0]);
                Point3d pt = (ptlist[0] + ptlist[1]) / 2; ps2.Add(pt);
                for (int i = 1; i < ptlist.Count - 1; i++)
                {
                    Point3d p1 = new Point3d(ptlist[i - 1]);
                    Point3d p2 = new Point3d(ptlist[i]);
                    Point3d p3 = new Point3d(ptlist[i + 1]);
                    Point3d p4 = (p2 + p3) / 2;
                    ps2.Add(p1 * (1 / 6.0) + p2 * (2 / 3.0) + p3 * (1 / 6.0));
                    ps2.Add(p4);
                }
                ps2.Add(ptlist[ptlist.Count - 1]);
            }
            else
            {
                if (ptlist.Count < 4)
                {
                    return ptlist;
                }
                Point3d p1 = new Point3d(ptlist[ptlist.Count - 2]);
                Point3d p2 = new Point3d(ptlist[0]);
                Point3d p3 = new Point3d(ptlist[0 + 1]);
                Point3d p4 = (p2 + p3) / 2;
                ps2.Add(p1 * (1 / 6.0) + p2 * (2 / 3.0) + p3 * (1 / 6.0));
                ps2.Add(p4);
                for (int i = 1; i < ptlist.Count - 1; i++)
                {
                    p1 = new Point3d(ptlist[i - 1]);
                    p2 = new Point3d(ptlist[i]);
                    p3 = new Point3d(ptlist[i + 1]);
                    p4 = (p2 + p3) / 2;
                    ps2.Add(p1 * (1 / 6.0) + p2 * (2 / 3.0) + p3 * (1 / 6.0));
                    ps2.Add(p4);
                }
                ps2.Add(ps2[0]);
            }

            Polyline pl2 = new Polyline(ps2);
            return pl2;
        }
        public static Polyline cc_Subdivide(Polyline ptlist, int level)
        {
            if (level >= 1)
            {
                for (int i = 0; i < level; i++)
                {
                    ptlist = cc_Subdivide(ptlist);
                }
            }
            return ptlist;
        }
    }
}
