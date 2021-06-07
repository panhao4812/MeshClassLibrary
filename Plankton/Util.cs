

using System;
using System.Collections.Generic;
using System.Linq;
using Rhino.Geometry;
using Rhino;

namespace Kangaroo
{
    public class Util
    {
        public static Point3d AveragePoint(List<Point3d> Pts)
        {
            Point3d pointd = new Point3d();
            foreach (Point3d pointd2 in Pts)
            {
                pointd += pointd2;
            }
            return (Point3d)(pointd * (1.0 / ((double)Pts.Count)));
        }
        public static List<Point3d>[] GetHingePoints(Mesh M)
        {
            M.Vertices.CombineIdentical(true, true);
            M.Faces.ConvertQuadsToTriangles();
            M.Weld(3.1415926535897931);
            List<Point3d> list = new List<Point3d>();
            List<Point3d> list2 = new List<Point3d>();
            List<Point3d> list3 = new List<Point3d>();
            List<Point3d> list4 = new List<Point3d>();
            M.Faces.ConvertQuadsToTriangles();
            Point3d[] pointdArray = M.Vertices.ToPoint3dArray();
            for (int i = 0; i < M.TopologyEdges.Count; i++)
            {
                int[] connectedFaces = M.TopologyEdges.GetConnectedFaces(i);
                if (connectedFaces.Length == 2)
                {
                    list.Add(M.TopologyEdges.EdgeLine(i).From);
                    list2.Add(M.TopologyEdges.EdgeLine(i).To);
                    M.TopologyVertices.IndicesFromFace(connectedFaces[0]);
                    M.TopologyVertices.IndicesFromFace(connectedFaces[1]);
                    IndexPair topologyVertices = M.TopologyEdges.GetTopologyVertices(i);
                    int[] numArray2 = M.TopologyVertices.ConnectedTopologyVertices(topologyVertices.I, true);
                    for (int j = 0; j < numArray2.Length; j++)
                    {
                        if (numArray2[j] == topologyVertices.J)
                        {
                            int topologyVertexIndex = numArray2[((j - 1) + numArray2.Length) % numArray2.Length];
                            int num4 = numArray2[((j + 1) + numArray2.Length) % numArray2.Length];
                            int[] numArray3 = M.TopologyVertices.MeshVertexIndices(topologyVertexIndex);
                            int[] numArray4 = M.TopologyVertices.MeshVertexIndices(num4);
                            list3.Add(pointdArray[numArray3[0]]);
                            list4.Add(pointdArray[numArray4[0]]);
                        }
                    }
                }
            }
            return new List<Point3d>[] { list, list2, list3, list4 };
        }
        public static List<Line> RemoveDupLn2(List<Line> lines, double tolerance)
        { return lines.Distinct<Line>(new LineEqualityComparer(tolerance)).ToList<Line>(); }
        public static List<Point3d> RemoveDupPts2(List<Point3d> mypoints, double tolerance)
        {
            return mypoints.Distinct<Point3d>(new PointEqualityComparer(tolerance)).ToList<Point3d>();
        }
        public static bool OrthoClose(Point3d Point1, Point3d Point2, double t)
        {
            return (((Math.Abs((double)(Point1.X - Point2.X)) < t) &&
                (Math.Abs((double)(Point1.Y - Point2.Y)) < t)) &&
                (Math.Abs((double)(Point1.Z - Point2.Z)) < t));
        }
        public static double MeshVol(Mesh M)
        {
            M.Faces.ConvertQuadsToTriangles();
            M.Vertices.CombineIdentical(true, true);
            M.Vertices.CullUnused();
            M.Weld(3.1415926535897931);
            Point3d[] pointdArray = M.Vertices.ToPoint3dArray();
            double num = 0.0;
            for (int i = 0; i < M.Faces.Count; i++)
            {
                MeshFace face = M.Faces[i];
                Point3d pointd = pointdArray[face.A];
                face = M.Faces[i];
                Point3d pointd2 = pointdArray[face.B];
                face = M.Faces[i];
                Point3d pointd3 = pointdArray[face.C];
                Vector3d vectord = Vector3d.CrossProduct((Vector3d)(pointd2 - pointd), (Vector3d)(pointd3 - pointd));
                num += vectord * ((Vector3d)pointd);
            }
            return (num / 6.0);
        }
        public static List<Curve> InterConnect(List<Point3d> pts)
        {
            List<Curve> list = new List<Curve>();
            for (int i = 0; i <= (pts.Count - 2); i++)
            {
                for (int j = i + 1; j <= (pts.Count - 1); j++)
                {
                    Line line = new Line(pts[i], pts[j]);
                    list.Add(NurbsCurve.CreateFromLine(line));
                }
            }
            return list;
        }
        /// //////////////////////////////////////
        public class PointEqualityComparer : IEqualityComparer<Point3d>
        {
            // Fields
            private double t;
            // Methods
            public PointEqualityComparer()
            {
            }
            public PointEqualityComparer(double _t)
            {
                this.t = _t;
            }
            public bool Equals(Point3d P1, Point3d P2)
            {
                return Util.OrthoClose(P1, P2, this.t);
            }
            public int GetHashCode(Point3d obj)
            {
                return 0;
            }
        }
        public class LineEqualityComparer : IEqualityComparer<Line>
        {
            // Fields
            private double t;
            // Methods
            public LineEqualityComparer()
            {
            }
            public LineEqualityComparer(double _t)
            {
                this.t = _t;
            }
            public bool Equals(Line L1, Line L2)
            {
                if ((!Util.OrthoClose(L1.From, L2.From, this.t) || !Util.OrthoClose(L1.To, L2.To, this.t)) && (!Util.OrthoClose(L1.To, L2.From, this.t) || !Util.OrthoClose(L1.From, L2.To, this.t)))
                {
                    return false;
                }
                return true;
            }
            public int GetHashCode(Line obj)
            {
                return 0;
            }
        }
    }
}
