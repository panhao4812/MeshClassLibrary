using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeoTools
{
    public class Polyline : List<Point3d>
    {
        public Polyline() { }
        public Polyline(IEnumerable<Point3d> collection) : base(collection) { }
        public Point3d PointAt(double t)
        {
            int count = Count;
            if (count < 1) { return Point3d.Origin; }
            if (count == 1) { return this[0]; }

            double floor = Math.Floor(t);

            int idx = (int)floor;
            if (idx < 0) { return this[0]; }
            if (idx >= count - 1) { return this[count - 1]; }

            t -= floor;
            if (t <= 0.0)
            {
                return this[idx];
            }
            if (t >= 1.0)
            {
                return this[idx + 1];
            }

            Point3d A = this[idx];
            Point3d B = this[idx + 1];

            double s = 1.0 - t;
            return new Point3d((A.m_x == B.m_x) ? A.m_x : s * A.m_x + t * B.m_x,
                               (A.m_y == B.m_y) ? A.m_y : s * A.m_y + t * B.m_y,
                               (A.m_z == B.m_z) ? A.m_z : s * A.m_z + t * B.m_z);
        }
        public Vector3d TangentAt(double t)
        {
            int count = Count;
            if (count < 2) { return Vector3d.Zero; }

            int segment_index = (int)Math.Floor(t);
            if (segment_index < 0)
            {
                segment_index = 0;
            }
            else if (segment_index > count - 2)
            {
                segment_index = count - 2;
            }

            Vector3d tangent = this[segment_index + 1] - this[segment_index];
            tangent.Unitize();

            return tangent;
        }
        public Polyline Trim(Interval domain)
        {
            int count = Count;
            int N = count - 1;

            // Polyline parameters
            double t0 = domain.Min;
            double t1 = domain.Max;

            // Segment indices
            int si0 = (int)Math.Floor(t0);
            int si1 = (int)Math.Floor(t1);

            // Segment parameters
            double st0 = t0 - si0;
            double st1 = t1 - si1;
            if (st0 < 0.0) { st0 = 0.0; }
            if (st0 >= 1.0) { si0++; st0 = 0.0; }
            if (st1 < 0.0) { st1 = 0.0; }
            if (st1 >= 1.0) { si1++; st1 = 0.0; }

            // Limit to polyline domain.
            if (si0 < 0) { si0 = 0; st0 = 0.0; }
            if (si0 >= N) { si0 = N; st0 = 0.0; }
            if (si1 < 0) { si1 = 0; st1 = 0.0; }
            if (si1 >= N) { si1 = N; st1 = 0.0; }

            // Build trimmed polyline.
            Polyline rc = new Polyline { PointAt(t0) };
            for (int i = si0 + 1; i <= si1; i++)
            {
                rc.Add(this[i]);
            }
            if (st1 > 0.0) { rc.Add(PointAt(t1)); }
            return rc;
        }
    }
    public struct Line
    {
        #region members
        internal Point3d m_from;
        internal Point3d m_to;
        public Point3d From
        {
            get { return m_from; }
            set { m_from = value; }
        }
        public double FromX
        {
            get { return m_from.X; }
            set { m_from.X = value; }
        }
        public double FromY
        {
            get { return m_from.Y; }
            set { m_from.Y = value; }
        }
        public double FromZ
        {
            get { return m_from.Z; }
            set { m_from.Z = value; }
        }
        public Point3d To
        {
            get { return m_to; }
            set { m_to = value; }
        }
        public double ToX
        {
            get { return m_to.X; }
            set { m_to.X = value; }
        }
        public double ToY
        {
            get { return m_to.Y; }
            set { m_to.Y = value; }
        }
        public double ToZ
        {
            get { return m_to.Z; }
            set { m_to.Z = value; }
        }
        #endregion
        #region constructors
        public Line(Point3d from, Point3d to)
        {
            m_from = from;
            m_to = to;
        }
        public Line(Point3d start, Vector3d span)
        {
            m_from = start;
            m_to = start + span;
        }
        public Line(Point3d start, Vector3d direction, double length)
        {
            Vector3d dir = direction;
            if (!dir.Unitize())
                dir = new Vector3d(0, 0, 1);

            m_from = start;
            m_to = start + dir * length;
        }
        public Line(double x0, double y0, double z0, double x1, double y1, double z1)
        {
            m_from = new Point3d(x0, y0, z0);
            m_to = new Point3d(x1, y1, z1);
        }
        #endregion
        #region operators
        public static bool operator ==(Line a, Line b)
        {
            return a.From == b.From && a.To == b.To;
        }
        public static bool operator !=(Line a, Line b)
        {
            return a.From != b.From || a.To != b.To;
        }
        public override bool Equals(object obj)
        {
            return obj is Line && this == (Line)obj;
        }
        public override int GetHashCode()
        {
            return From.GetHashCode() ^ To.GetHashCode();
        }
        #endregion
        public bool IsValid
        {
            get
            {
                return From.IsValid && To.IsValid;
            }
        }
        public double Length
        {
            get { return From.DistanceTo(To); }
            set
            {
                Vector3d dir = To - From;
                if (!dir.Unitize())
                    dir = new Vector3d(0, 0, 1);

                To = From + dir * value;
            }
        }
        public Vector3d Direction
        {
            get { return To - From; }
        }
        public Vector3d UnitTangent
        {
            get
            {
                Vector3d v = To - From;
                v.Unitize();
                return v;
            }
        }
        public bool Equals(Line other)
        {
            return this == other;
        }
        public void Flip()
        {
            Point3d temp = From;
            From = To;
            To = temp;
        }
        public Point3d PointAt(double t)
        {
            double s = 1.0 - t;
            return new Point3d((From.m_x == To.m_x) ? From.m_x : s * From.m_x + t * To.m_x,
                               (From.m_y == To.m_y) ? From.m_y : s * From.m_y + t * To.m_y,
                               (From.m_z == To.m_z) ? From.m_z : s * From.m_z + t * To.m_z);
        }
        public double ClosestParameter(Point3d point)
        {
            double t = 0;
            Vector3d D = Direction;
            double DoD = D.LengthSquared();
            if (DoD > 0.0)
            {
                if (point.DistanceTo(From) <= point.DistanceTo(To))
                {
                    t = ((point - From) * D) / DoD;
                }
                else
                {
                    t = 1.0 + ((point - To) * D) / DoD;
                }
            }
            return t;
        }
        public Point3d ClosestPoint(Point3d testPoint, bool limitToFiniteSegment)
        {
            double t = ClosestParameter(testPoint);
            if (limitToFiniteSegment)
            {
                t = Math.Max(t, 0.0);
                t = Math.Min(t, 1.0);
            }
            return PointAt(t);
        }
        public double DistanceTo(Point3d testPoint, bool limitToFiniteSegment)
        {
            Point3d pp = ClosestPoint(testPoint, limitToFiniteSegment);
            return pp.DistanceTo(testPoint);
        }
        public void Transform(Transform xform)
        {
            this.m_from.Transform(xform);
            this.m_to.Transform(xform);
        }
        public bool Extend(double startLength, double endLength)
        {
            if (Length == 0.0) { return false; }
            Point3d A = m_from;
            Point3d B = m_to;

            Vector3d tan = UnitTangent;

            if (startLength != 0.0) { A = m_from - startLength * tan; }
            if (endLength != 0.0) { B = m_to + endLength * tan; }

            m_from = A;
            m_to = B;

            return true;
        }
        public override string ToString()
        {
            return string.Format("{0},{1}", From.ToString(), To.ToString());
        }
        public static bool TryFitLineToPoints(List<Point3d> points, out Line fitLine)
        {
            Plane p;
            GeoSolver.LinePlaneEstimate(points, out fitLine, out p);
            if (fitLine.IsValid) return true;
            return false;  
        }
        
        public static bool TryFitLineToPoints(Point3d[] ptArray, out Line fitLine)
        {
            Plane p;
            GeoSolver.LinePlaneEstimate(new List < Point3d > (ptArray), out fitLine, out p);
            if (fitLine.IsValid) return true;
            return false;
        }
        public bool TryGetPlane(out Plane plane, double tolerance)
        {
            plane = new Plane();
            Vector3d v = To - From;
            bool bTinyX = (Math.Abs(v.X) <= tolerance);
            bool bTinyY = (Math.Abs(v.Y) <= tolerance);
            bool bTinyZ = (Math.Abs(v.Z) <= tolerance);
            bool rc = true;
            Vector3d X = new Vector3d();
            Vector3d Y = new Vector3d();
            if (bTinyZ && (!bTinyX || !bTinyY))
            {
                X = Vector3d.XAxis;
                Y = Vector3d.YAxis;
            }
            else if (bTinyX && (!bTinyY || !bTinyZ))
            {
                X = Vector3d.YAxis;
                Y = Vector3d.ZAxis;
            }
            else if (bTinyY && (!bTinyZ || !bTinyX))
            {
                X = Vector3d.ZAxis;
                Y = Vector3d.XAxis;
            }
            else
            {
                X = v;
                X.Unitize();
                Y.PerpendicularTo(X);
                if (bTinyX && bTinyY && bTinyZ)
                {
                    rc = false;
                    if (X.IsZero)
                    {
                        X = Vector3d.XAxis;
                        Y = Vector3d.YAxis;
                    }
                }
            }
            plane.CreateFromFrame(From, X, Y);
            return rc;
        }
      
    }
}
