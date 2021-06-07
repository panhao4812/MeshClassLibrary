using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Kangaroo
{
    public struct Interval
    {
        #region Members
        private double m_t0;
        private double m_t1;
        public double T0 { get { return m_t0; } set { m_t0 = value; } }

        /// <summary>
        /// Gets or sets the upper bound of the Interval.
        /// </summary>
        public double T1 { get { return m_t1; } set { m_t1 = value; } }

        /// <summary>
        /// Gets or sets the indexed bound of this Interval.
        /// </summary>
        /// <param name="index">Bound index (0 = lower; 1 = upper).</param>
        public double this[int index]
        {
            get
            {
                if (0 == index) { return m_t0; }
                if (1 == index) { return m_t1; }

                // IronPython works with indexing is we thrown an IndexOutOfRangeException
                throw new IndexOutOfRangeException();
            }
            set
            {
                if (0 == index) { m_t0 = value; }
                else if (1 == index) { m_t1 = value; }
                else { throw new IndexOutOfRangeException(); }
            }
        }
        #endregion
        public Interval(double t0, double t1)
        {
            m_t0 = t0;
            m_t1 = t1;
        }
        public Interval(Interval other)
        {
            m_t0 = other.m_t0;
            m_t1 = other.m_t1;
        }
        #region Operators
        public static bool operator ==(Interval a, Interval b)
        {
            return a.CompareTo(b) == 0;
        }
        public static bool operator !=(Interval a, Interval b)
        {
            return a.CompareTo(b) != 0;
        }
        public static Interval operator +(Interval interval, double number)
        {
            return new Interval(interval.m_t0 + number, interval.m_t1 + number);
        }
        public static Interval operator +(double number, Interval interval)
        {
            return new Interval(interval.m_t0 + number, interval.m_t1 + number);
        }
        public static Interval operator -(Interval interval, double number)
        {
            return new Interval(interval.m_t0 - number, interval.m_t1 - number);
        }
        public static Interval operator -(double number, Interval interval)
        {
            return new Interval(number - interval.m_t0, number - interval.m_t1);
        }
        public static bool operator <(Interval a, Interval b)
        {
            return a.CompareTo(b) < 0;
        }
        public static bool operator <=(Interval a, Interval b)
        {
            return a.CompareTo(b) <= 0;
        }
        public static bool operator >(Interval a, Interval b)
        {
            return a.CompareTo(b) > 0;
        }
        public static bool operator >=(Interval a, Interval b)
        {
            return a.CompareTo(b) >= 0;
        }
        public override int GetHashCode()
        {
            // MSDN docs recommend XOR'ing the internal values to get a hash code
            return m_t0.GetHashCode() ^ m_t1.GetHashCode();
        }
        public override bool Equals(object obj)
        {
            return (obj is Interval && this == (Interval)obj);
        }
        #endregion
        public double Min
        {
            get { return (m_t0 <= m_t1 ? m_t0 : m_t1); }
        }
        public double Max
        {
            get { return (m_t0 <= m_t1 ? m_t1 : m_t0); }
        }
        public double Mid
        {
            get { return (m_t0 == m_t1 ? m_t0 : 0.5 * (m_t0 + m_t1)); }
        }
        public double Length
        {
            get { return m_t1 - m_t0; }
        }
        public bool Equals(Interval other)
        {
            return this == other;
        }
        public int CompareTo(Interval other)
        {
            if (m_t0 < other.m_t0)
                return -1;
            if (m_t0 > other.m_t0)
                return 1;
            if (m_t1 < other.m_t1)
                return -1;
            if (m_t1 > other.m_t1)
                return 1;
            return 0;
        }
    }
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
    }
    public struct Point3d
    {
        #region operators
        public static Point3d operator +(Point3d point1, Point3d point2)
        {
            return new Point3d(point1.m_x + point2.m_x, point1.m_y + point2.m_y, point1.m_z + point2.m_z);
        }
        public static Point3d operator +(Point3d point1, Vector3d point2)
        {
            return new Point3d(point1.m_x + point2.m_x, point1.m_y + point2.m_y, point1.m_z + point2.m_z);
        }
        public static Vector3d operator -(Point3d point1, Point3d point2)
        {
            return new Vector3d(point1.m_x - point2.m_x, point1.m_y - point2.m_y, point1.m_z - point2.m_z);
        }
        public static Point3d operator -(Point3d point1, Vector3d point2)
        {
            return new Point3d(point1.m_x - point2.m_x, point1.m_y - point2.m_y, point1.m_z - point2.m_z);
        }
        public static Vector3d Subtract(Point3d point1, Point3d point2)
        {
            return new Vector3d(point1.m_x - point2.m_x, point1.m_y - point2.m_y, point1.m_z - point2.m_z);
        }
        public static Point3d operator /(Point3d point, double t)
        {
            return new Point3d(point.m_x / t, point.m_y / t, point.m_z / t);
        }
        public static Point3d operator *(Point3d point, double t)
        {
            return new Point3d(point.m_x * t, point.m_y * t, point.m_z * t);
        }
        public static Point3d operator *(double t, Point3d point)
        {
            return new Point3d(point.m_x * t, point.m_y * t, point.m_z * t);
        }
        public static bool operator ==(Point3d a, Point3d b)
        {
            return a.m_x == b.m_x && a.m_y == b.m_y && a.m_z == b.m_z;
        }
        public static bool operator !=(Point3d a, Point3d b)
        {
            return a.m_x != b.m_x || a.m_y != b.m_y || a.m_z != b.m_z;
        }
        public override bool Equals(object obj)
        {
            return (obj is Point3d && this == (Point3d)obj);
        }
        public override int GetHashCode()
        {
            // MSDN docs recommend XOR'ing the internal values to get a hash code
            return m_x.GetHashCode() ^ m_y.GetHashCode() ^ m_z.GetHashCode();
        }
        #endregion
        #region members
        internal double m_x;
        internal double m_y;
        internal double m_z;
        public double X { get { return m_x; } set { m_x = value; } }

        public double Y { get { return m_y; } set { m_y = value; } }

        public double Z { get { return m_z; } set { m_z = value; } }
        #endregion
        public Point3d(double x, double y, double z)
        {
            m_x = x;
            m_y = y;
            m_z = z;
        }
        public Point3d(Point3d p)
        {
            m_x = p.m_x;
            m_y = p.m_y;
            m_z = p.m_z;
        }
        public Point3d(Vector3d v)
        {
            m_x = v.m_x;
            m_y = v.m_y;
            m_z = v.m_z;
        }
        public static Point3d Origin
        {
            get { return new Point3d(); }
        }
        public override string ToString()
        {
            return m_x.ToString() + "," + m_y.ToString() + "," + m_z.ToString();
        }
        public void Transform(Transform xform)
        {
            //David: this method doesn't test for validity. Should it?
            double ww = xform.m_30 * m_x + xform.m_31 * m_y + xform.m_32 * m_z + xform.m_33;
            if (ww != 0.0) { ww = 1.0 / ww; }

            double tx = ww * (xform.m_00 * m_x + xform.m_01 * m_y + xform.m_02 * m_z + xform.m_03);
            double ty = ww * (xform.m_10 * m_x + xform.m_11 * m_y + xform.m_12 * m_z + xform.m_13);
            double tz = ww * (xform.m_20 * m_x + xform.m_21 * m_y + xform.m_22 * m_z + xform.m_23);
            m_x = tx;
            m_y = ty;
            m_z = tz;
        }
        public static Point3d GetCenter(List<Point3d> pts)
        {
            Point3d cen = new Point3d();
            if (pts.Count < 1) return cen;
            for (int i = 0; i < pts.Count; i++)
            {
                cen += pts[i];
            }
            cen /= pts.Count;
            return cen;
        }
        public double DistanceTo(Point3d other)
        {
            return Math.Sqrt(Math.Pow(X - other.X, 2) + Math.Pow(Y - other.Y, 2) + Math.Pow(Z - other.Z, 2));
        }
        public double DistanceTo(Point3f other)
        {
            return Math.Sqrt(
                Math.Pow(X - other.X, 2) +
                Math.Pow(Y - other.Y, 2) +
                Math.Pow(Z - other.Z, 2));
        }
    }
    public struct Vector3d
    {
        #region operators
        public static Vector3d operator *(Vector3d vector, double t)
        {
            return new Vector3d(vector.m_x * t, vector.m_y * t, vector.m_z * t);
        }
        public static Vector3d Multiply(Vector3d vector, double t)
        {
            return new Vector3d(vector.m_x * t, vector.m_y * t, vector.m_z * t);
        }
        public static Vector3d operator *(double t, Vector3d vector)
        {
            return new Vector3d(vector.m_x * t, vector.m_y * t, vector.m_z * t);
        }
        public static Vector3d Multiply(double t, Vector3d vector)
        {
            return new Vector3d(vector.m_x * t, vector.m_y * t, vector.m_z * t);
        }
        public static Vector3d operator /(Vector3d vector, double t)
        {
            return new Vector3d(vector.m_x / t, vector.m_y / t, vector.m_z / t);
        }
        public static Vector3d Divide(Vector3d vector, double t)
        {
            return new Vector3d(vector.m_x / t, vector.m_y / t, vector.m_z / t);
        }
        public static Vector3d operator +(Vector3d vector1, Vector3d vector2)
        {
            return new Vector3d(vector1.m_x + vector2.m_x, vector1.m_y + vector2.m_y, vector1.m_z + vector2.m_z);
        }
        public static Vector3d Add(Vector3d vector1, Vector3d vector2)
        {
            return new Vector3d(vector1.m_x + vector2.m_x, vector1.m_y + vector2.m_y, vector1.m_z + vector2.m_z);
        }
        public static Vector3d operator -(Vector3d vector1, Vector3d vector2)
        {
            return new Vector3d(vector1.m_x - vector2.m_x, vector1.m_y - vector2.m_y, vector1.m_z - vector2.m_z);
        }
        public static Vector3d Subtract(Vector3d vector1, Vector3d vector2)
        {
            return new Vector3d(vector1.m_x - vector2.m_x, vector1.m_y - vector2.m_y, vector1.m_z - vector2.m_z);
        }
        public static double operator *(Vector3d vector1, Vector3d vector2)
        {
            return (vector1.m_x * vector2.m_x + vector1.m_y * vector2.m_y + vector1.m_z * vector2.m_z);
        }
        public static double Multiply(Vector3d vector1, Vector3d vector2)
        {
            return (vector1.m_x * vector2.m_x + vector1.m_y * vector2.m_y + vector1.m_z * vector2.m_z);
        }
        public static double DotProduct(Vector3d vector1, Vector3d vector2)
        {
            return vector1 * vector2;
        }
        public static Vector3d operator -(Vector3d vector)
        {
            return new Vector3d(-vector.m_x, -vector.m_y, -vector.m_z);
        }
        public static Vector3d Negate(Vector3d vector)
        {
            return new Vector3d(-vector.m_x, -vector.m_y, -vector.m_z);
        }
        public static bool operator ==(Vector3d a, Vector3d b)
        {
            return a.m_x == b.m_x && a.m_y == b.m_y && a.m_z == b.m_z;
        }
        public static bool operator !=(Vector3d a, Vector3d b)
        {
            return a.m_x != b.m_x || a.m_y != b.m_y || a.m_z != b.m_z;
        }
        public static Vector3d CrossProduct(Vector3d a, Vector3d b)
        {
            return new Vector3d(a.m_y * b.m_z - b.m_y * a.m_z, a.m_z * b.m_x - b.m_z * a.m_x, a.m_x * b.m_y - b.m_x * a.m_y);
        }
        public static double VectorAngle(Vector3d a, Vector3d b)
        {
            if (!a.Unitize() || !b.Unitize())
                return double.NaN;

            //compute dot product
            double dot = a.m_x * b.m_x + a.m_y * b.m_y + a.m_z * b.m_z;
            // remove any "noise"
            if (dot > 1.0) dot = 1.0;
            if (dot < -1.0) dot = -1.0;
            double radians = Math.Acos(dot);
            return radians;
        }
        public override bool Equals(object obj)
        {
            return (obj is Vector3d && this == (Vector3d)obj);
        }
        public override int GetHashCode()
        {
            // MSDN docs recommend XOR'ing the internal values to get a hash code
            return m_x.GetHashCode() ^ m_y.GetHashCode() ^ m_z.GetHashCode();
        }
        #endregion
        #region members
        internal double m_x;
        internal double m_y;
        internal double m_z;
        public double X { get { return m_x; } set { m_x = value; } }

        public double Y { get { return m_y; } set { m_y = value; } }

        public double Z { get { return m_z; } set { m_z = value; } }
        #endregion
        public Vector3d(double x, double y, double z)
        {
            m_x = x;
            m_y = y;
            m_z = z;
        }
        public Vector3d(Point3d p)
        {
            m_x = p.m_x;
            m_y = p.m_y;
            m_z = p.m_z;
        }
        public Vector3d(Vector3d v)
        {
            m_x = v.m_x;
            m_y = v.m_y;
            m_z = v.m_z;
        }
        public static Vector3d Zero
        {
            get { return new Vector3d(); }
        }
        public static Vector3d XAxis
        {
            get { return new Vector3d(1.0, 0.0, 0.0); }
        }
        public static Vector3d YAxis
        {
            get { return new Vector3d(0.0, 1.0, 0.0); }
        }
        public static Vector3d ZAxis
        {
            get { return new Vector3d(0.0, 0.0, 1.0); }
        }
    public double Length
        {
            get {
                if (this == Vector3d.Zero) return 0;
                return Math.Sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
            }
        }    
        public double LengthSquared()
        {
            return (m_x * m_x + m_y * m_y + m_z * m_z);
        }
        public bool Unitize()
        {
            bool rc = false;
            double d = Length();
            if (d > 0.0)
            {
                m_x = (double)(m_x / d);
                m_y = (double)(m_y / d);
                m_z = (double)(m_z / d);
                rc = true;
            }
            return rc;
        }
        public void Transform(Transform transformation)
        {
            double xx = transformation.m_00 * m_x + transformation.m_01 * m_y + transformation.m_02 * m_z;
            double yy = transformation.m_10 * m_x + transformation.m_11 * m_y + transformation.m_12 * m_z;
            double zz = transformation.m_20 * m_x + transformation.m_21 * m_y + transformation.m_22 * m_z;

            m_x = xx;
            m_y = yy;
            m_z = zz;
        }
        public override string ToString()
        {
            return m_x.ToString() + "," + m_y.ToString() + "," + m_z.ToString();
        }
        public bool PerpendicularTo(Vector3d v)
        {
            //bool rc = false;
            int i, j, k;
            double a, b;
            k = 2;
            if (Math.Abs(v.Y) > Math.Abs(v.X))
            {
                if (Math.Abs(v.Z) > Math.Abs(v.Y))
                {
                    // |v.Z| > |v.Y| > |v.X|
                    i = 2;
                    j = 1;
                    k = 0;
                    a = v.Z;
                    b = -v.Y;
                }
                else if (Math.Abs(v.Z) >= Math.Abs(v.X))
                {
                    // |v.Y| >= |v.Z| >= |v.X|
                    i = 1;
                    j = 2;
                    k = 0;
                    a = v.Y;
                    b = -v.Z;
                }
                else
                {
                    // |v.Y| > |v.X| > |v.Z|
                    i = 1;
                    j = 0;
                    k = 2;
                    a = v.Y;
                    b = -v.X;
                }
            }
            else if (Math.Abs(v.Z) > Math.Abs(v.X))
            {
                // |v.Z| > |v.X| >= |v.Y|
                i = 2;
                j = 0;
                k = 1;
                a = v.Z;
                b = -v.X;
            }
            else if (Math.Abs(v.Z) > Math.Abs(v.Y))
            {
                // |v.X| >= |v.Z| > |v.Y|
                i = 0;
                j = 2;
                k = 1;
                a = v.X;
                b = -v.Z;
            }
            else
            {
                // |v.X| >= |v.Y| >= |v.Z|
                i = 0;
                j = 1;
                k = 2;
                a = v.X;
                b = -v.Y;
            }
            double[] this_v = new double[3];
            this_v[i] = b;
            this_v[j] = a;
            this_v[k] = 0.0;
            m_x = this_v[0];
            m_y = this_v[1];
            m_z = this_v[2];
            return (a != 0.0) ? true : false;
        }
        public bool isPerpendicularTo(Vector3d v, double angle_tolerance)
        {
            bool rc = false;
            double ll = Length() * v.Length();
            if (ll > 0.0)
            {
                if (Math.Abs((X * v.X + Y * v.Y + Z * v.Z) / ll) <= Math.Sin(angle_tolerance))
                    rc = true;
            }
            return rc;
        }
        public int IsParallelTo(Vector3d v, double angle_tolerance)
        {
            int rc = 0;
            double ll = Length() * v.Length();
            if (ll > 0.0)
            {
                double cos_angle = (X * v.X + Y * v.Y + Z * v.Z) / ll;
                double cos_tol = Math.Cos(angle_tolerance);
                if (cos_angle >= cos_tol)
                    rc = 1;
                else if (cos_angle <= -cos_tol)
                    rc = -1;
            }
            return rc;
        }
        public void Reverse()
        {
            m_x = -m_x;
            m_y = -m_y;
            m_z = -m_z;
        }
        public bool PerpendicularTo(Point3d P0, Point3d P1, Point3d P2)
        {
            Vector3d V0, V1, V2, N0, N1, N2;
            V0 = P2 - P1;
            V1 = P0 - P2;
            V2 = P1 - P0;
            N0 = Vector3d.CrossProduct(V1, V2);
            if (!N0.Unitize())
                return false;
            N1 = Vector3d.CrossProduct(V2, V0);
            if (!N1.Unitize())
                return false;
            N2 = Vector3d.CrossProduct(V0, V1);
            if (!N2.Unitize())
                return false;
            double s0 = 1.0 / V0.Length();
            double s1 = 1.0 / V1.Length();
            double s2 = 1.0 / V2.Length();
            double e0 = s0 * Math.Abs(Vector3d.DotProduct(N0, V0)) + s1 * Math.Abs(Vector3d.DotProduct(N0, V1)) + s2 * Math.Abs(Vector3d.DotProduct(N0, V2));
            double e1 = s0 * Math.Abs(Vector3d.DotProduct(N1, V0)) + s1 * Math.Abs(Vector3d.DotProduct(N1, V1)) + s2 * Math.Abs(Vector3d.DotProduct(N1, V2));
            double e2 = s0 * Math.Abs(Vector3d.DotProduct(N2, V0)) + s1 * Math.Abs(Vector3d.DotProduct(N2, V1)) + s2 * Math.Abs(Vector3d.DotProduct(N2, V2));
            if (e0 <= e1)
            {
                if (e0 <= e2) { this = N0; }
                else { this = N2; }
            }
            else if (e1 <= e2) { this = N1; }
            else { this = N2; }
            return true;
        }
    }
    public struct Transform
    {
        #region operators
        public static bool operator ==(Transform a, Transform b)
        {
            return a.m_00 == b.m_00 && a.m_01 == b.m_01 && a.m_02 == b.m_02 && a.m_03 == b.m_03 &&
              a.m_10 == b.m_10 && a.m_11 == b.m_11 && a.m_12 == b.m_12 && a.m_13 == b.m_13 &&
              a.m_20 == b.m_20 && a.m_21 == b.m_21 && a.m_22 == b.m_22 && a.m_23 == b.m_23 &&
              a.m_30 == b.m_30 && a.m_31 == b.m_31 && a.m_32 == b.m_32 && a.m_33 == b.m_33;
        }
        public static bool operator !=(Transform a, Transform b)
        {
            return a.m_00 != b.m_00 || a.m_01 != b.m_01 || a.m_02 != b.m_02 || a.m_03 != b.m_03 ||
              a.m_10 != b.m_10 || a.m_11 != b.m_11 || a.m_12 != b.m_12 || a.m_13 != b.m_13 ||
              a.m_20 != b.m_20 || a.m_21 != b.m_21 || a.m_22 != b.m_22 || a.m_23 != b.m_23 ||
              a.m_30 != b.m_30 || a.m_31 != b.m_31 || a.m_32 != b.m_32 || a.m_33 != b.m_33;
        }
        public static Transform operator *(Transform a, Transform b)
        {
            Transform xf = new Transform();
            xf.m_00 = a.m_00 * b.m_00 + a.m_01 * b.m_10 + a.m_02 * b.m_20 + a.m_03 * b.m_30;
            xf.m_01 = a.m_00 * b.m_01 + a.m_01 * b.m_11 + a.m_02 * b.m_21 + a.m_03 * b.m_31;
            xf.m_02 = a.m_00 * b.m_02 + a.m_01 * b.m_12 + a.m_02 * b.m_22 + a.m_03 * b.m_32;
            xf.m_03 = a.m_00 * b.m_03 + a.m_01 * b.m_13 + a.m_02 * b.m_23 + a.m_03 * b.m_33;

            xf.m_10 = a.m_10 * b.m_00 + a.m_11 * b.m_10 + a.m_12 * b.m_20 + a.m_13 * b.m_30;
            xf.m_11 = a.m_10 * b.m_01 + a.m_11 * b.m_11 + a.m_12 * b.m_21 + a.m_13 * b.m_31;
            xf.m_12 = a.m_10 * b.m_02 + a.m_11 * b.m_12 + a.m_12 * b.m_22 + a.m_13 * b.m_32;
            xf.m_13 = a.m_10 * b.m_03 + a.m_11 * b.m_13 + a.m_12 * b.m_23 + a.m_13 * b.m_33;

            xf.m_20 = a.m_20 * b.m_00 + a.m_21 * b.m_10 + a.m_22 * b.m_20 + a.m_23 * b.m_30;
            xf.m_21 = a.m_20 * b.m_01 + a.m_21 * b.m_11 + a.m_22 * b.m_21 + a.m_23 * b.m_31;
            xf.m_22 = a.m_20 * b.m_02 + a.m_21 * b.m_12 + a.m_22 * b.m_22 + a.m_23 * b.m_32;
            xf.m_23 = a.m_20 * b.m_03 + a.m_21 * b.m_13 + a.m_22 * b.m_23 + a.m_23 * b.m_33;

            xf.m_30 = a.m_30 * b.m_00 + a.m_31 * b.m_10 + a.m_32 * b.m_20 + a.m_33 * b.m_30;
            xf.m_31 = a.m_30 * b.m_01 + a.m_31 * b.m_11 + a.m_32 * b.m_21 + a.m_33 * b.m_31;
            xf.m_32 = a.m_30 * b.m_02 + a.m_31 * b.m_12 + a.m_32 * b.m_22 + a.m_33 * b.m_32;
            xf.m_33 = a.m_30 * b.m_03 + a.m_31 * b.m_13 + a.m_32 * b.m_23 + a.m_33 * b.m_33;
            return xf;
        }
        public static Transform Multiply(Transform a, Transform b)
        {
            return a * b;
        }
        public static Point3d operator *(Transform m, Point3d p)
        {
            double x = p.m_x; // optimizer should put x,y,z in registers
            double y = p.m_y;
            double z = p.m_z;
            Point3d rc = new Point3d();
            rc.m_x = m.m_00 * x + m.m_01 * y + m.m_02 * z + m.m_03;
            rc.m_y = m.m_10 * x + m.m_11 * y + m.m_12 * z + m.m_13;
            rc.m_z = m.m_20 * x + m.m_21 * y + m.m_22 * z + m.m_23;
            double w = m.m_30 * x + m.m_31 * y + m.m_32 * z + m.m_33;
            if (w != 0.0)
            {
                w = 1.0 / w;
                rc.m_x *= w;
                rc.m_y *= w;
                rc.m_z *= w;
            }
            return rc;
        }
        public static Point3d Multiply(Transform m, Point3d p)
        {
            return m * p;
        }
        public static Vector3d operator *(Transform m, Vector3d v)
        {
            double x = v.m_x; // optimizer should put x,y,z in registers
            double y = v.m_y;
            double z = v.m_z;
            Vector3d rc = new Vector3d();
            rc.m_x = m.m_00 * x + m.m_01 * y + m.m_02 * z;
            rc.m_y = m.m_10 * x + m.m_11 * y + m.m_12 * z;
            rc.m_z = m.m_20 * x + m.m_21 * y + m.m_22 * z;
            return rc;
        }
        public static Vector3d Multiply(Transform m, Vector3d v)
        {
            return m * v;
        }
        public static Transform operator +(Transform a, Transform b)
        {
            Transform xf = new Transform();
            xf.m_00 = a.m_00 + b.m_00;
            xf.m_01 = a.m_01 + b.m_01;
            xf.m_02 = a.m_02 + b.m_02;
            xf.m_03 = a.m_03 + b.m_03;

            xf.m_10 = a.m_10 + b.m_10;
            xf.m_11 = a.m_11 + b.m_11;
            xf.m_12 = a.m_12 + b.m_12;
            xf.m_13 = a.m_13 + b.m_13;

            xf.m_20 = a.m_20 + b.m_20;
            xf.m_21 = a.m_21 + b.m_21;
            xf.m_22 = a.m_22 + b.m_22;
            xf.m_23 = a.m_23 + b.m_23;

            xf.m_30 = a.m_30 + b.m_30;
            xf.m_31 = a.m_31 + b.m_31;
            xf.m_32 = a.m_32 + b.m_32;
            xf.m_33 = a.m_33 + b.m_33;
            return xf;
        }
        public static Transform operator -(Transform a, Transform b)
        {
            Transform xf = new Transform();
            xf.m_00 = a.m_00 - b.m_00;
            xf.m_01 = a.m_01 - b.m_01;
            xf.m_02 = a.m_02 - b.m_02;
            xf.m_03 = a.m_03 - b.m_03;

            xf.m_10 = a.m_10 - b.m_10;
            xf.m_11 = a.m_11 - b.m_11;
            xf.m_12 = a.m_12 - b.m_12;
            xf.m_13 = a.m_13 - b.m_13;

            xf.m_20 = a.m_20 - b.m_20;
            xf.m_21 = a.m_21 - b.m_21;
            xf.m_22 = a.m_22 - b.m_22;
            xf.m_23 = a.m_23 - b.m_23;

            xf.m_30 = a.m_30 - b.m_30;
            xf.m_31 = a.m_31 - b.m_31;
            xf.m_32 = a.m_32 - b.m_32;
            xf.m_33 = a.m_33 - b.m_33;
            return xf;
        }
        public override bool Equals(object obj)
        {
            return obj is Transform && Equals((Transform)obj);
        }
        public override int GetHashCode()
        {
            // MSDN docs recommend XOR'ing the internal values to get a hash code
            return m_00.GetHashCode() ^ m_01.GetHashCode() ^ m_02.GetHashCode() ^ m_03.GetHashCode() ^
                   m_10.GetHashCode() ^ m_11.GetHashCode() ^ m_12.GetHashCode() ^ m_13.GetHashCode() ^
                   m_20.GetHashCode() ^ m_21.GetHashCode() ^ m_22.GetHashCode() ^ m_23.GetHashCode() ^
                   m_30.GetHashCode() ^ m_31.GetHashCode() ^ m_32.GetHashCode() ^ m_33.GetHashCode();
        }
        #endregion
        #region members
        internal double m_00, m_01, m_02, m_03;
        internal double m_10, m_11, m_12, m_13;
        internal double m_20, m_21, m_22, m_23;
        internal double m_30, m_31, m_32, m_33;
        public double this[int row, int column]
        {
            get
            {
                if (row < 0) { row = 0; }
                if (row > 3) { row = 3; }
                if (column < 0) { column = 0; }
                if (column > 3) { column = 3; }
                if (row == 0)
                {
                    if (column == 0) { return m_00; }
                    if (column == 1) { return m_01; }
                    if (column == 2) { return m_02; }
                    if (column == 3) { return m_03; }
                }
                else if (row == 1)
                {
                    if (column == 0) { return m_10; }
                    if (column == 1) { return m_11; }
                    if (column == 2) { return m_12; }
                    if (column == 3) { return m_13; }
                }
                else if (row == 2)
                {
                    if (column == 0) { return m_20; }
                    if (column == 1) { return m_21; }
                    if (column == 2) { return m_22; }
                    if (column == 3) { return m_23; }
                }
                else if (row == 3)
                {
                    if (column == 0) { return m_30; }
                    if (column == 1) { return m_31; }
                    if (column == 2) { return m_32; }
                    if (column == 3) { return m_33; }
                }
                return double.NaN;
            }
            set
            {
                if (row < 0) { row = 0; }
                if (row > 3) { row = 3; }
                if (column < 0) { column = 0; }
                if (column > 3) { column = 3; }
                if (row == 0)
                {
                    if (column == 0)
                    { m_00 = value; }
                    else if (column == 1)
                    { m_01 = value; }
                    else if (column == 2)
                    { m_02 = value; }
                    else if (column == 3)
                    { m_03 = value; }
                }
                else if (row == 1)
                {
                    if (column == 0)
                    { m_10 = value; }
                    else if (column == 1)
                    { m_11 = value; }
                    else if (column == 2)
                    { m_12 = value; }
                    else if (column == 3)
                    { m_13 = value; }
                }
                else if (row == 2)
                {
                    if (column == 0)
                    { m_20 = value; }
                    else if (column == 1)
                    { m_21 = value; }
                    else if (column == 2)
                    { m_22 = value; }
                    else if (column == 3)
                    { m_23 = value; }
                }
                else if (row == 3)
                {
                    if (column == 0)
                    { m_30 = value; }
                    else if (column == 1)
                    { m_31 = value; }
                    else if (column == 2)
                    { m_32 = value; }
                    else if (column == 3)
                    { m_33 = value; }
                }
            }
        }
        public double M00 { get { return m_00; } set { m_00 = value; } }
        /// <summary>Gets or sets this[0,1].</summary>
        public double M01 { get { return m_01; } set { m_01 = value; } }
        /// <summary>Gets or sets this[0,2].</summary>
        public double M02 { get { return m_02; } set { m_02 = value; } }
        /// <summary>Gets or sets this[0,3].</summary>
        public double M03 { get { return m_03; } set { m_03 = value; } }
        /// <summary>Gets or sets this[1,0].</summary>
        public double M10 { get { return m_10; } set { m_10 = value; } }
        /// <summary>Gets or sets this[1,1].</summary>
        public double M11 { get { return m_11; } set { m_11 = value; } }
        /// <summary>Gets or sets this[1,2].</summary>
        public double M12 { get { return m_12; } set { m_12 = value; } }
        /// <summary>Gets or sets this[1,3].</summary>
        public double M13 { get { return m_13; } set { m_13 = value; } }

        /// <summary>Gets or sets this[2,0].</summary>
        public double M20 { get { return m_20; } set { m_20 = value; } }
        /// <summary>Gets or sets this[2,1].</summary>
        public double M21 { get { return m_21; } set { m_21 = value; } }
        /// <summary>Gets or sets this[2,2].</summary>
        public double M22 { get { return m_22; } set { m_22 = value; } }
        /// <summary>Gets or sets this[2,3].</summary>
        public double M23 { get { return m_23; } set { m_23 = value; } }

        /// <summary>Gets or sets this[3,0].</summary>
        public double M30 { get { return m_30; } set { m_30 = value; } }
        /// <summary>Gets or sets this[3,1].</summary>
        public double M31 { get { return m_31; } set { m_31 = value; } }
        /// <summary>Gets or sets this[3,2].</summary>
        public double M32 { get { return m_32; } set { m_32 = value; } }
        /// <summary>Gets or sets this[3,3].</summary>
        public double M33 { get { return m_33; } set { m_33 = value; } }
        #endregion
        public static Transform Zero
        {
            get
            {
                Transform xf = new Transform();
                return xf;
            }
        }
        public static Transform Identity
        {
            get
            {
                Transform xf = new Transform();
                xf.m_00 = 1.0;
                xf.m_11 = 1.0;
                xf.m_22 = 1.0;
                xf.m_33 = 1.0;
                return xf;
            }
        }
        public Transform(Point3d P, Vector3d X, Vector3d Y, Vector3d Z)
        {
            m_00 = X.X;
            m_01 = X.Y;
            m_02 = X.Z;
            m_03 = 0;

            m_10 = Y.X;
            m_11 = Y.Y;
            m_12 = Y.Z;
            m_13 = 0;

            m_20 = Z.X;
            m_21 = Z.Y;
            m_22 = Z.Z;
            m_23 = 0;

            m_30 = P.X;
            m_31 = P.Y;
            m_32 = P.Z;
            m_33 = 1;
        }
        public static Transform Translation(double dx, double dy, double dz)
        {
            Transform xf = Identity;
            xf.m_03 = dx;
            xf.m_13 = dy;
            xf.m_23 = dz;
            xf.m_33 = 1.0;
            return xf;
        }
        public static Transform Translation(Vector3d v)
        {
            return Translation(v.X, v.Y, v.Z);
        }
        public static Transform Mirror(Plane P)
        {
            return Mirror(P.Origin, P.Normal);
        }
        public static Transform Mirror(Point3d P, Vector3d N)
        {
            Transform output = Transform.Identity;
            N.Unitize();
            Vector3d V = (2.0 * (N.X * P.X + N.Y * P.Y + N.Z * P.Z)) * N;
            output.m_00 = 1 - 2.0 * N.X * N.X;
            output.m_01 = -2.0 * N.X * N.Y;
            output.m_02 = -2.0 * N.X * N.Z;
            output.m_03 = V.X;

            output.m_10 = -2.0 * N.Y * N.X;
            output.m_11 = 1.0 - 2.0 * N.Y * N.Y;
            output.m_12 = -2.0 * N.Y * N.Z;
            output.m_13 = V.Y;

            output.m_20 = -2.0 * N.Z * N.X;
            output.m_21 = -2.0 * N.Z * N.Y;
            output.m_22 = 1.0 - 2.0 * N.Z * N.Z;
            output.m_23 = V.Z;

            output.m_30 = 0.0;
            output.m_31 = 0.0;
            output.m_32 = 0.0;
            output.m_33 = 1.0;
            return output;
        }
        public static Transform Scale(Vector3d v)
        {
            return Scale(v.X, v.Y, v.Z);
        }
        public static Transform Scale(double x, double y, double z)
        {
            Transform xf = new Transform();
            xf.m_00 = x;
            xf.m_11 = y;
            xf.m_22 = z;
            xf.m_33 = 1.0;
            return xf;
        }
        public static Transform Scale(Point3d fixed_point, double x, double y, double z)
        {
            if (fixed_point.X == 0.0 && fixed_point.Y == 0.0 && fixed_point.Z == 0.0)
            {
                return Scale(x, y, z);
            }
            else
            {
                Transform tr0 = Transform.Translation(Point3d.Origin - fixed_point);
                Transform tr1 = Transform.Scale(x, y, z);
                Transform tr2 = Transform.Translation(fixed_point - Point3d.Origin);
                return (tr0 * tr1 * tr2);
            }
        }
        public static Transform Scale(Plane plane, Vector3d x1, Vector3d y1, Vector3d z1)
        {
            return Shear(plane, x1, y1, z1);
        }
        public static Transform Shear(Plane plane, Vector3d x1, Vector3d y1, Vector3d z1)
        {

            Transform t0 = Transform.Translation(Point3d.Origin - plane.Origin);
            Transform s0 = new Transform();
            Transform s1 = new Transform();
            s0.m_00 = plane.XAxis.X;
            s0.m_01 = plane.XAxis.Y;
            s0.m_02 = plane.XAxis.Z;
            s0.m_10 = plane.YAxis.X;
            s0.m_11 = plane.YAxis.Y;
            s0.m_12 = plane.YAxis.Z;
            s0.m_20 = plane.ZAxis.X;
            s0.m_21 = plane.ZAxis.Y;
            s0.m_22 = plane.ZAxis.Z;
            s1.m_00 = x1.X;
            s1.m_01 = x1.Y;
            s1.m_02 = x1.Z;
            s1.m_10 = y1.X;
            s1.m_11 = y1.Y;
            s1.m_12 = y1.Z;
            s1.m_20 = z1.X;
            s1.m_21 = z1.Y;
            s1.m_22 = z1.Z;
            Transform t1 = Transform.Translation(plane.Origin - Point3d.Origin);
            return (t1 * s1 * s0 * t0);
        }
        public static Transform ChangeBasis(Plane plane0, Plane plane1)
        {
            return Transform.ChangeBasis(
            plane0.Origin, plane0.XAxis, plane0.YAxis, plane0.ZAxis,
            plane1.Origin, plane1.XAxis, plane1.YAxis, plane1.ZAxis
            );
        }
        public static Transform ChangeBasis(Point3d P0, Vector3d X0, Vector3d Y0, Vector3d Z0, Point3d P1, Vector3d X1, Vector3d Y1, Vector3d Z1)
        {
            Transform F0 = new Transform(P0, X0, Y0, Z0);		// Frame 0
            Transform T1 = Transform.Translation(-P1.X, -P1.Y, -P1.Z);
            Transform CB = Transform.ChangeBasis(Vector3d.XAxis, Vector3d.YAxis, Vector3d.ZAxis, X1, Y1, Z1);
            return CB * T1 * F0;
        }
        public static Transform ChangeBasis(Vector3d X0, Vector3d Y0, Vector3d Z0, Vector3d X1, Vector3d Y1, Vector3d Z1)
        {
            Transform xform = Transform.Identity;
            double a, b, c, d;
            a = X1 * Y1;
            b = X1 * Z1;
            c = Y1 * Z1;
            double[,] R = new double[3,6]{
                {X1*X1,a,b,X1*X0,X1*Y0,X1*Z0},
                {a,Y1*Y1,c,Y1*X0,Y1*Y0,Y1*Z0},
                {b,c,Z1*Z1,Z1*X0,Z1*Y0,Z1*Z0}};
            //double R[3,6] = {{X1*X1,      a,      b,       X0*X1, X0*Y1, X0*Z1},
            //                 {    a,  Y1*Y1,      c,       Y0*X1, Y0*Y1, Y0*Z1},
            //                 {    b,      c,  Z1*Z1,       Z0*X1, Z0*Y1, Z0*Z1}};
            // row reduce R
            int i0 = (R[0, 0] >= R[1, 1]) ? 0 : 1;
            if (R[2, 2] > R[i0, i0])
                i0 = 2;
            int i1 = (i0 + 1) % 3;
            int i2 = (i1 + 1) % 3;
            if (R[i0, i0] == 0.0)
                return xform;
            d = 1.0 / R[i0, i0];
            R[i0, 0] *= d;
            R[i0, 1] *= d;
            R[i0, 2] *= d;
            R[i0, 3] *= d;
            R[i0, 4] *= d;
            R[i0, 5] *= d;
            R[i0, i0] = 1.0;
            if (R[i1, i0] != 0.0)
            {
                d = -R[i1, i0];
                R[i1, 0] += d * R[i0, 0];
                R[i1, 1] += d * R[i0, 1];
                R[i1, 2] += d * R[i0, 2];
                R[i1, 3] += d * R[i0, 3];
                R[i1, 4] += d * R[i0, 4];
                R[i1, 5] += d * R[i0, 5];
                R[i1, i0] = 0.0;
            }
            if (R[i2, i0] != 0.0)
            {
                d = -R[i2, i0];
                R[i2, 0] += d * R[i0, 0];
                R[i2, 1] += d * R[i0, 1];
                R[i2, 2] += d * R[i0, 2];
                R[i2, 3] += d * R[i0, 3];
                R[i2, 4] += d * R[i0, 4];
                R[i2, 5] += d * R[i0, 5];
                R[i2, i0] = 0.0;
            }
            if (Math.Abs(R[i1, i1]) < Math.Abs(R[i2, i2]))
            {
                int i = i1; i1 = i2; i2 = i;
            }
            if (R[i1, i1] == 0.0)
                return xform;
            d = 1.0 / R[i1, i1];
            R[i1, 0] *= d;
            R[i1, 1] *= d;
            R[i1, 2] *= d;
            R[i1, 3] *= d;
            R[i1, 4] *= d;
            R[i1, 5] *= d;
            R[i1, i1] = 1.0;
            if (R[i0, i1] != 0.0)
            {
                d = -R[i0, i1];
                R[i0, 0] += d * R[i1, 0];
                R[i0, 1] += d * R[i1, 1];
                R[i0, 2] += d * R[i1, 2];
                R[i0, 3] += d * R[i1, 3];
                R[i0, 4] += d * R[i1, 4];
                R[i0, 5] += d * R[i1, 5];
                R[i0, i1] = 0.0;
            }
            if (R[i2, i1] != 0.0)
            {
                d = -R[i2, i1];
                R[i2, 0] += d * R[i1, 0];
                R[i2, 1] += d * R[i1, 1];
                R[i2, 2] += d * R[i1, 2];
                R[i2, 3] += d * R[i1, 3];
                R[i2, 4] += d * R[i1, 4];
                R[i2, 5] += d * R[i1, 5];
                R[i2, i1] = 0.0;
            }

            if (R[i2, i2] == 0.0)
                return xform;
            d = 1.0 / R[i2, i2];
            R[i2, 0] *= d;
            R[i2, 1] *= d;
            R[i2, 2] *= d;
            R[i2, 3] *= d;
            R[i2, 4] *= d;
            R[i2, 5] *= d;
            R[i2, i2] = 1.0;
            if (R[i0, i2] != 0.0)
            {
                d = -R[i0, i2];
                R[i0, 0] += d * R[i2, 0];
                R[i0, 1] += d * R[i2, 1];
                R[i0, 2] += d * R[i2, 2];
                R[i0, 3] += d * R[i2, 3];
                R[i0, 4] += d * R[i2, 4];
                R[i0, 5] += d * R[i2, 5];
                R[i0, i2] = 0.0;
            }
            if (R[i1, i2] != 0.0)
            {
                d = -R[i1, i2];
                R[i1, 0] += d * R[i2, 0];
                R[i1, 1] += d * R[i2, 1];
                R[i1, 2] += d * R[i2, 2];
                R[i1, 3] += d * R[i2, 3];
                R[i1, 4] += d * R[i2, 4];
                R[i1, 5] += d * R[i2, 5];
                R[i1, i2] = 0.0;
            }

            xform.M00 = R[0, 3];
            xform.M01 = R[0, 4];
            xform.M02 = R[0, 5];
            xform.M10 = R[1, 3];
            xform.M11 = R[1, 4];
            xform.M12 = R[1, 5];
            xform.M20 = R[2, 3];
            xform.M21 = R[2, 4];
            xform.M22 = R[2, 5];
            return xform;
        }
        public static Transform PlaneToPlane(Plane plane0, Plane plane1)
        {
            Transform rc = Transform.Identity;
            rc = ChangeBasis(plane0, plane1);
            return rc;
        }
        public void Transpose()
        {
            double t;
            t = M01; M01 = M10; M10 = t;
            t = M02; M02 = M20; M20 = t;
            t = M03; M03 = M30; M30 = t;
            t = M12; M12 = M21; M21 = t;
            t = M13; M13 = M31; M31 = t;
            t = M23; M23 = M32; M32 = t;
        }
        public static Transform PlanarProjection(Plane plane)
        {
            Transform xf = Identity;
            int i, j;
            double[] x = new double[3] { plane.XAxis.X, plane.XAxis.Y, plane.XAxis.Z };
            double[] y = new double[3] { plane.YAxis.X, plane.YAxis.Y, plane.YAxis.Z };
            double[] p = new double[3] { plane.OriginX, plane.OriginY, plane.OriginZ };
            double[] q = new double[3];
            for (i = 0; i < 3; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    xf[i, j] = x[i] * x[j] + y[i] * y[j];
                }
                q[i] = xf[i, 0] * p[0] + xf[i, 1] * p[1] + xf[i, 2] * p[2];
            }
            for (i = 0; i < 3; i++)
            {
                xf[3, i] = 0.0;
                xf[i, 3] = p[i] - q[i];
            }
            xf[3, 3] = 1.0;
            return xf;
        }
        public Transform(double[,] xf)
        {
            this = Transform.Identity;
            int ii = xf.GetLength(0);
            int jj = xf.GetLength(1);
            if (ii > 4) ii = 4;
            if (jj > 4) jj = 4;         
                for (int i = 0; i < ii; i++)
                {
                    for (int j = 0; j < jj; j++)
                    {
                        this[i, j] = xf[i, j];
                    }
                }           
        }
        public double[,] ToArray()
        {
            double[,] M = new double[4, 4];
            M[0, 0] = M00; M[1, 0] = M10; M[2, 0] = M20; M[3, 0] = M30;
            M[0, 1] = M01; M[1, 1] = M11; M[2, 1] = M21; M[3, 1] = M31;
            M[0, 2] = M02; M[1, 2] = M12; M[2, 2] = M22; M[3, 2] = M32;
            M[0, 3] = M03; M[1, 3] = M13; M[2, 3] = M23; M[3, 3] = M33;
            return M;
        }
        public int Invert()
        {
            // return (rank == 4) ? true : false;
            // returns rank (0, 1, 2, 3, or 4), inverse, and smallest pivot
            double[,] M = ToArray(), I = new double[4, 4];
            int[] col = new int[4] { 0, 1, 2, 3 };
            int swapcount = 0, rank = 0;
            double pivot = 0.0, determinant = 0.0;
            int ix, jx; double x, c, d;
            I[0, 0] = I[1, 1] = I[2, 2] = I[3, 3] = 1.0;
            ix = jx = 0;
            x = Math.Abs(M[0, 0]);
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    if (Math.Abs(M[i, j]) > x)
                    {
                        ix = i;
                        jx = j;
                        x = Math.Abs(M[i, j]);
                    }
                }
            }
            pivot = x;
            if (ix != 0)
            {
                SwapRow(ref M, 0, ix);
                SwapRow(ref I, 0, ix);
                swapcount++;
            }
            if (jx != 0)
            {
                SwapCol(ref M, 0, jx);
                col[0] = jx;
                swapcount++;
            }
            if (x > 0.0)
            {
                rank++;
                c = M[0, 0];
                M[0, 1] /= c; M[0, 2] /= c; M[0, 3] /= c;
                I[0, 0] /= c; I[0, 1] /= c; I[0, 2] /= c; I[0, 3] /= c;
                d = 1.0 / c;
                x *= double.Epsilon;
                if (Math.Abs(M[1, 0]) > x)
                {
                    c = -M[1, 0];
                    M[1, 1] += c * M[0, 1]; M[1, 2] += c * M[0, 2]; M[1, 3] += c * M[0, 3];
                    AddCxRow(ref I, c, 0, 1);
                }
                if (Math.Abs(M[2, 0]) > x)
                {
                    c = -M[2, 0];
                    M[2, 1] += c * M[0, 1]; M[2, 2] += c * M[0, 2]; M[2, 3] += c * M[0, 3];
                    AddCxRow(ref I, c, 0, 2);
                }
                if (Math.Abs(M[3, 0]) > x)
                {
                    c = -M[3, 0];
                    M[3, 1] += c * M[0, 1]; M[3, 2] += c * M[0, 2]; M[3, 3] += c * M[0, 3];
                    AddCxRow(ref I, c, 0, 3);
                }
                ix = jx = 1;
                x = Math.Abs(M[1, 1]);
                for (int i = 1; i < 4; i++)
                {
                    for (int j = 1; j < 4; j++)
                    {
                        if (Math.Abs(M[i, j]) > x)
                        {
                            ix = i;
                            jx = j;
                            x = Math.Abs(M[i, j]);
                        }
                    }
                }
                if (x < pivot)
                    pivot = x;
                if (ix != 1)
                {
                    SwapRow(ref M, 1, ix);
                    SwapRow(ref I, 1, ix);
                    swapcount++;
                }
                if (jx != 1)
                {
                    SwapCol(ref M, 1, jx);
                    col[1] = jx;
                    swapcount++;
                }
                if (x > 0.0)
                {
                    rank++;
                    c = M[1, 1];
                    M[1, 2] /= c; M[1, 3] /= c;
                    I[1, 0] /= c; I[1, 1] /= c; I[1, 2] /= c; I[1, 3] /= c;
                    d /= c;
                    x *= double.Epsilon;
                    if (Math.Abs(M[0, 1]) > x)
                    {
                        c = -M[0, 1];
                        M[0, 2] += c * M[1, 2]; M[0, 3] += c * M[1, 3];
                        AddCxRow(ref I, c, 1, 0);
                    }
                    if (Math.Abs(M[2, 1]) > x)
                    {
                        c = -M[2, 1];
                        M[2, 2] += c * M[1, 2]; M[2, 3] += c * M[1, 3];
                        AddCxRow(ref I, c, 1, 2);
                    }
                    if (Math.Abs(M[3, 1]) > x)
                    {
                        c = -M[3, 1];
                        M[3, 2] += c * M[1, 2]; M[3, 3] += c * M[1, 3];
                        AddCxRow(ref I, c, 1, 3);
                    }
                    ix = jx = 2;
                    x = Math.Abs(M[2, 2]);
                    for (int i = 2; i < 4; i++) for (int j = 2; j < 4; j++)
                        {
                            if (Math.Abs(M[i, j]) > x)
                            {
                                ix = i;
                                jx = j;
                                x = Math.Abs(M[i, j]);
                            }
                        }
                    if (x < pivot)
                        pivot = x;
                    if (ix != 2)
                    {
                        SwapRow(ref M, 2, ix);
                        SwapRow(ref I, 2, ix);
                        swapcount++;
                    }
                    if (jx != 2)
                    {
                        SwapCol(ref M, 2, jx);
                        col[2] = jx;
                        swapcount++;
                    }
                    if (x > 0.0)
                    {
                        rank++;
                        c = M[2, 2];
                        M[2, 3] /= c;
                        I[2, 0] /= c; I[2, 1] /= c; I[2, 2] /= c; I[2, 3] /= c;
                        d /= c;
                        x *= double.Epsilon;
                        if (Math.Abs(M[0, 2]) > x)
                        {
                            c = -M[0, 2];
                            M[0, 3] += c * M[2, 3];
                            AddCxRow(ref I, c, 2, 0);
                        }
                        if (Math.Abs(M[1, 2]) > x)
                        {
                            c = -M[1, 2];
                            M[1, 3] += c * M[2, 3];
                            AddCxRow(ref I, c, 2, 1);
                        }
                        if (Math.Abs(M[3, 2]) > x)
                        {
                            c = -M[3, 2];
                            M[3, 3] += c * M[2, 3];
                            AddCxRow(ref I, c, 2, 3);
                        }
                        x = Math.Abs(M[3, 3]);
                        if (x < pivot) pivot = x;
                        if (x > 0.0)
                        {
                            rank++;
                            c = M[3, 3];
                            I[3, 0] /= c; I[3, 1] /= c; I[3, 2] /= c; I[3, 3] /= c;
                            d /= c;
                            x *= double.Epsilon;
                            if (Math.Abs(M[0, 3]) > x)
                            {
                                AddCxRow(ref I, -M[0, 3], 3, 0);
                            }
                            if (Math.Abs(M[1, 3]) > x)
                            {
                                AddCxRow(ref I, -M[1, 3], 3, 1);
                            }
                            if (Math.Abs(M[2, 3]) > x)
                            {
                                AddCxRow(ref I, -M[2, 3], 3, 2);
                            }
                            determinant = ((swapcount % 2) != 0) ? -d : d;
                        }
                    }
                }
            }
            if (col[3] != 3)
                SwapRow(ref I, 3, col[3]);
            if (col[2] != 2)
                SwapRow(ref I, 2, col[2]);
            if (col[1] != 1)
                SwapRow(ref I, 1, col[1]);
            if (col[0] != 0)
                SwapRow(ref I, 0, col[0]);
            // memcpy(dst, I, sizeof(I) );
            this = new Transform(I);
            return rank;
        }
        void AddCxRow(ref double[,] matrix, double c, int i0, int i1)
        {
            matrix[i1, 0] += c * matrix[i0, 0];
            matrix[i1, 1] += c * matrix[i0, 1];
            matrix[i1, 2] += c * matrix[i0, 2];
            matrix[i1, 3] += c * matrix[i0, 3];
        }
        void SwapRow(ref double[,] matrix, int i0, int i1)
        {
            double t;
            t = matrix[i0, 0]; matrix[i0, 0] = matrix[i1, 0]; matrix[i1, 0] = t;
            t = matrix[i0, 1]; matrix[i0, 1] = matrix[i1, 1]; matrix[i1, 1] = t;
            t = matrix[i0, 2]; matrix[i0, 2] = matrix[i1, 2]; matrix[i1, 2] = t;
            t = matrix[i0, 3]; matrix[i0, 3] = matrix[i1, 3]; matrix[i1, 3] = t;
        }
        void SwapCol(ref double[,] matrix, int i0, int i1)
        {
            double t;
            t = matrix[0, i0]; matrix[0, i0] = matrix[0, i1]; matrix[0, i1] = t;
            t = matrix[1, i0]; matrix[1, i0] = matrix[1, i1]; matrix[1, i1] = t;
            t = matrix[2, i0]; matrix[2, i0] = matrix[2, i1]; matrix[2, i1] = t;
            t = matrix[3, i0]; matrix[3, i0] = matrix[3, i1]; matrix[3, i1] = t;
        }
        public override string ToString()
        {
            return M00.ToString() + "," + M01.ToString() + "," + M02.ToString() + "," + M03.ToString() + "\r\n"
            + M10.ToString() + "," + M11.ToString() + "," + M12.ToString() + "," + M13.ToString() + "\r\n"
            + M20.ToString() + "," + M21.ToString() + "," + M22.ToString() + "," + M23.ToString() + "\r\n"
            + M30.ToString() + "," + M31.ToString() + "," + M32.ToString() + "," + M33.ToString();
        }
        public static Transform Rotation(double angle, Vector3d axis, Point3d center)
        {
            return Rotation(Math.Sin(angle), Math.Cos(angle), axis, center);
        }
        public static Transform Rotation(double sin_angle, double cos_angle, Vector3d axis, Point3d center)
        {
            Transform m_xform = Transform.Identity;
            for (;;)
            {
                if (Math.Abs(sin_angle) >= 1.0 && Math.Abs(cos_angle) <= 0)
                {
                    cos_angle = 0.0;
                    sin_angle = (sin_angle < 0.0) ? -1.0 : 1.0;
                    break;
                }
                if (Math.Abs(cos_angle) >= 1.0 && Math.Abs(sin_angle) <= 0)
                {
                    cos_angle = (cos_angle < 0.0) ? -1.0 : 1.0;
                    sin_angle = 0.0;
                    break;
                }
                if (Math.Abs(cos_angle * cos_angle + sin_angle * sin_angle - 1.0) > 0)
                {
                    Vector3d cs = new Vector3d(cos_angle, sin_angle, 0);
                    if (cs.Unitize())
                    {
                        cos_angle = cs.X;
                        sin_angle = cs.Y;
                        // no break here
                    }
                    else
                    {
                        cos_angle = 1.0;
                        sin_angle = 0.0;
                        break;
                    }
                }
                if (Math.Abs(cos_angle) > 1.0 || Math.Abs(sin_angle) < 0)
                {
                    cos_angle = (cos_angle < 0.0) ? -1.0 : 1.0;
                    sin_angle = 0.0;
                    break;
                }
                if (Math.Abs(sin_angle) > 1.0 || Math.Abs(cos_angle) < 0)
                {
                    cos_angle = 0.0;
                    sin_angle = (sin_angle < 0.0) ? -1.0 : 1.0;
                    break;
                }
                break;
            }
            if (sin_angle != 0.0 || cos_angle != 1.0)
            {
                double one_minus_cos_angle = 1.0 - cos_angle;
                Vector3d a = axis;
                if (Math.Abs(a.LengthSquared() - 1.0) > 0)
                    a.Unitize();
                m_xform[0, 0] = a.X * a.X * one_minus_cos_angle + cos_angle;
                m_xform[0, 1] = a.X * a.Y * one_minus_cos_angle - a.Z * sin_angle;
                m_xform[0, 2] = a.X * a.Z * one_minus_cos_angle + a.Y * sin_angle;

                m_xform[1, 0] = a.Y * a.X * one_minus_cos_angle + a.Z * sin_angle;
                m_xform[1, 1] = a.Y * a.Y * one_minus_cos_angle + cos_angle;
                m_xform[1, 2] = a.Y * a.Z * one_minus_cos_angle - a.X * sin_angle;

                m_xform[2, 0] = a.Z * a.X * one_minus_cos_angle - a.Y * sin_angle;
                m_xform[2, 1] = a.Z * a.Y * one_minus_cos_angle + a.X * sin_angle;
                m_xform[2, 2] = a.Z * a.Z * one_minus_cos_angle + cos_angle;

                if (center.X != 0.0 || center.Y != 0.0 || center.Z != 0.0)
                {
                    m_xform[0, 3] = -((m_xform[0, 0] - 1.0) * center.X + m_xform[0, 1] * center.Y + m_xform[0, 2] * center.Z);
                    m_xform[1, 3] = -(m_xform[1, 0] * center.X + (m_xform[1, 1] - 1.0) * center.Y + m_xform[1, 2] * center.Z);
                    m_xform[2, 3] = -(m_xform[2, 0] * center.X + m_xform[2, 1] * center.Y + (m_xform[2, 2] - 1.0) * center.Z);
                }
                m_xform[3, 0] = m_xform[3, 1] = m_xform[3, 2] = 0.0;
                m_xform[3, 3] = 1.0;
            }
            return m_xform;
        }
        public static Transform Rotation(Vector3d start_dir, Vector3d end_dir, Point3d rotation_center)
        {
            if (Math.Abs(start_dir.Length() - 1.0) > 0)
                start_dir.Unitize();
            if (Math.Abs(end_dir.Length() - 1.0) > 0)
                end_dir.Unitize();
            double cos_angle = start_dir * end_dir;
            Vector3d axis = Vector3d.CrossProduct(start_dir, end_dir);
            double sin_angle = axis.Length();
            if (0.0 == sin_angle || !axis.Unitize())
            {
                axis.PerpendicularTo(start_dir);
                axis.Unitize();
                sin_angle = 0.0;
                cos_angle = (cos_angle < 0.0) ? -1.0 : 1.0;
            }
            return Rotation(sin_angle, cos_angle, axis, rotation_center);
        }
    }
    public struct Plane
    {
        #region members
        internal Point3d m_origin;
        internal Vector3d m_xaxis;
        internal Vector3d m_yaxis;
        internal Vector3d m_zaxis;
        public Point3d Origin
        {
            get { return m_origin; }
            set { m_origin = value; }
        }
        public double OriginX
        {
            get { return m_origin.X; }
            set { m_origin.X = value; }
        }
        public double OriginY
        {
            get { return m_origin.Y; }
            set { m_origin.Y = value; }
        }
        public double OriginZ
        {
            get { return m_origin.Z; }
            set { m_origin.Z = value; }
        }
        public Vector3d XAxis
        {
            get { return m_xaxis; }
            set { m_xaxis = value; }
        }
        public Vector3d YAxis
        {
            get { return m_yaxis; }
            set { m_yaxis = value; }
        }
        public Vector3d ZAxis
        {
            get { return m_zaxis; }
            set { m_zaxis = value; }
        }
        public static Plane WorldXY
        {
            get
            {
                return new Plane { XAxis = new Vector3d(1, 0, 0), YAxis = new Vector3d(0, 1, 0), ZAxis = new Vector3d(0, 0, 1) };
            }
        }
        public static Plane WorldYZ
        {
            get
            {
                return new Plane { XAxis = new Vector3d(0, 1, 0), YAxis = new Vector3d(0, 0, 1), ZAxis = new Vector3d(1, 0, 0) };
            }
        }
        public static Plane WorldZX
        {
            get
            {
                return new Plane { XAxis = new Vector3d(0, 0, 1), YAxis = new Vector3d(1, 0, 0), ZAxis = new Vector3d(0, 1, 0) };
            }
        }
        public Vector3d Normal
        {
            get { return ZAxis; }
        }
        #endregion
        #region constructors
        public Plane(Plane other)
        {
            this = other;
        }
        public Plane(Point3d P, Vector3d N) : this()
        {
            CreateFromNormal(P, N);
        }
        public Plane(double x, double y, double z, double a, double b, double c) : this()
        {
            CreateFromNormal(new Point3d(x, y, z), new Vector3d(a, b, c));
        }
        public Plane(Point3d origin, Vector3d xDirection, Vector3d yDirection) : this()
        {
            CreateFromFrame(origin, xDirection, yDirection);
        }
        public Plane(Point3d origin, Point3d xPoint, Point3d yPoint) : this()
        {
            CreateFromPoints(origin, xPoint, yPoint);
        }
        public Plane(double a, double b, double c, double d) : this()
        {
            CreateFromEquation(a, b, c, d);
            this.m_zaxis.Unitize();
        }
        public bool CreateFromNormal(Point3d P, Vector3d N)
        {
            m_origin = P;
            m_zaxis = N;
            bool b = m_zaxis.Unitize();
            m_xaxis.PerpendicularTo(m_zaxis);
            m_xaxis.Unitize();
            m_yaxis = Vector3d.CrossProduct(m_zaxis, m_xaxis);
            m_yaxis.Unitize();
            return b;
        }
        public bool CreateFromFrame(Point3d P, Vector3d X, Vector3d Y)
        {
            m_origin = P;
            m_xaxis = X;
            m_xaxis.Unitize();
            m_yaxis = Y - Vector3d.DotProduct(Y, m_xaxis) * m_xaxis;
            m_yaxis.Unitize();
            m_zaxis = Vector3d.CrossProduct(m_xaxis, m_yaxis);
            bool b = m_zaxis.Unitize();
            return b;
        }
        public bool CreateFromPoints(Point3d P, Point3d Q, Point3d R)
        {
            m_origin = P;
            bool rc = m_zaxis.PerpendicularTo(P, Q, R);
            m_xaxis = Q - P;
            m_xaxis.Unitize();
            m_yaxis = Vector3d.CrossProduct(m_zaxis, m_xaxis);
            m_yaxis.Unitize();
            return rc;
        }
        public bool CreateFromEquation(double X, double Y, double Z, double D)
        {
            bool b = false;
            m_zaxis.X = X;
            m_zaxis.Y = Y;
            m_zaxis.Z = Z;
            double d = m_zaxis.Length();
            if (d > 0.0)
            {
                d = 1.0 / d;
                m_origin = new Point3d(m_zaxis * (-d * D));
                b = true;
            }
            m_xaxis.PerpendicularTo(m_zaxis);
            m_xaxis.Unitize();
            m_yaxis = Vector3d.CrossProduct(m_zaxis, m_xaxis);
            m_yaxis.Unitize();
            return b;
        }
        #endregion
        #region operators
        public static bool operator ==(Plane a, Plane b)
        {
            return a.Equals(b);
        }
        public static bool operator !=(Plane a, Plane b)
        {
            return (a.m_origin != b.m_origin) ||
                   (a.m_xaxis != b.m_xaxis) ||
                   (a.m_yaxis != b.m_yaxis) ||
                   (a.m_zaxis != b.m_zaxis);
        }
        public override bool Equals(object obj)
        {
            return ((obj is Plane) && (this == (Plane)obj));
        }
        public bool Equals(Plane plane)
        {
            return (m_origin == plane.m_origin) &&
                   (m_xaxis == plane.m_xaxis) &&
                   (m_yaxis == plane.m_yaxis) &&
                   (m_zaxis == plane.m_zaxis);
        }
        public override int GetHashCode()
        {
            return m_origin.GetHashCode() ^ m_xaxis.GetHashCode() ^ m_yaxis.GetHashCode() ^ m_zaxis.GetHashCode();
        }
        public override string ToString()
        {
            return "";
        }
        #endregion
        public static double[] PlaneEquationCreate(Point3d P, Vector3d N)
        {
            double[] output = new double[4];
            output[0] = N.X;
            output[1] = N.Y;
            output[2] = N.Z;
            output[3] = -(new Vector3d(P) * N);
            return output;
        }
        public double[] GetPlaneEquation()
        {
            return Plane.PlaneEquationCreate(this.Origin, this.Normal);
        }
        public Point3d ClosestPoint(Point3d P)
        {
            double x = m_zaxis.X; double y = m_zaxis.Y; double z = m_zaxis.Z;
            double d = -(new Vector3d(m_origin) * m_zaxis);
            double t = -(x * P.X + y * P.Y + z * P.Z + d) / (x * x + y * y + z * z);
            return new Point3d(P.X + t * x, P.Y + t * y, P.Z + t * z);
        }
        public double DistanceTo(Point3d testPoint)
        {
            return (testPoint - m_origin) * m_zaxis;
        }
        public void Flip()
        {
            Vector3d v = m_xaxis;
            m_xaxis = m_yaxis;
            m_yaxis = v;
            m_zaxis = -m_zaxis;
        }
        public bool Transform(Transform xform)
        {
            m_origin.Transform(xform);
            m_xaxis.Transform(xform);
            m_yaxis.Transform(xform);
            return CreateFromFrame(m_origin, m_xaxis, m_yaxis);
        }
    }
    public struct Point3f
    {
        #region members
        internal float m_x;
        internal float m_y;
        internal float m_z;
        public float X { get { return m_x; } set { m_x = value; } }

        public float Y { get { return m_y; } set { m_y = value; } }

        public float Z { get { return m_z; } set { m_z = value; } }
        #endregion
        #region operators
        public static Point3f operator +(Point3f point1, Point3f point2)
        {
            return new Point3f(point1.m_x + point2.m_x, point1.m_y + point2.m_y, point1.m_z + point2.m_z);
        }
        public static Point3f operator +(Point3f point1, Vector3d point2)
        {
            return new Point3f(point1.m_x + point2.m_x, point1.m_y + point2.m_y, point1.m_z + point2.m_z);
        }
        public static Vector3d operator -(Point3f point1, Point3f point2)
        {
            return new Vector3d(point1.m_x - point2.m_x, point1.m_y - point2.m_y, point1.m_z - point2.m_z);
        }
        public static Point3f operator -(Point3f point1, Vector3d point2)
        {
            return new Point3f(point1.m_x - point2.m_x, point1.m_y - point2.m_y, point1.m_z - point2.m_z);
        }
        public static Vector3d Subtract(Point3f point1, Point3f point2)
        {
            return new Vector3d(point1.m_x - point2.m_x, point1.m_y - point2.m_y, point1.m_z - point2.m_z);
        }
        public static Point3f operator /(Point3f point, float t)
        {
            return new Point3f(point.m_x / t, point.m_y / t, point.m_z / t);
        }
        public static Point3f operator *(Point3f point, float t)
        {
            return new Point3f(point.m_x * t, point.m_y * t, point.m_z * t);
        }
        public static bool operator ==(Point3f a, Point3f b)
        {
            return a.m_x == b.m_x && a.m_y == b.m_y && a.m_z == b.m_z;
        }
        public static bool operator !=(Point3f a, Point3f b)
        {
            return a.m_x != b.m_x || a.m_y != b.m_y || a.m_z != b.m_z;
        }
        public override bool Equals(object obj)
        {
            return (obj is Point3f && this == (Point3f)obj);
        }
        public override int GetHashCode()
        {
            // MSDN docs recommend XOR'ing the internal values to get a hash code
            return m_x.GetHashCode() ^ m_y.GetHashCode() ^ m_z.GetHashCode();
        }
        #endregion
        public Point3f(float x, float y, float z)
        {
            m_x = x;
            m_y = y;
            m_z = z;
        }
        public Point3f(double x, double y, double z)
        {
            m_x = (float)x;
            m_y = (float)y;
            m_z = (float)z;
        }
        public Point3f(Point3d p)
        {
            m_x = (float)p.m_x;
            m_y = (float)p.m_y;
            m_z = (float)p.m_z;
        }
        public float DistanceTo(Point3f other)
        {
            return (float)Math.Sqrt(
                Math.Pow(X - other.X, 2) +
                Math.Pow(Y - other.Y, 2) +
                Math.Pow(Z - other.Z, 2));
        }
        public float DistanceTo(Point3d other)
        {
            return (float)Math.Sqrt(
                Math.Pow(X - (float)other.X, 2) +
                Math.Pow(Y - (float)other.Y, 2) +
                Math.Pow(Z - (float)other.Z, 2));
        }
    }
    public struct Vector3f
    {
        #region members
        internal float m_x;
        internal float m_y;
        internal float m_z;
        public float X { get { return m_x; } set { m_x = value; } }

        public float Y { get { return m_y; } set { m_y = value; } }

        public float Z { get { return m_z; } set { m_z = value; } }
        #endregion
        #region operators
        public static Vector3f operator *(Vector3f vector, float t)
        {
            return new Vector3f(vector.m_x * t, vector.m_y * t, vector.m_z * t);
        }
        public static Vector3f Multiply(Vector3f vector, float t)
        {
            return new Vector3f(vector.m_x * t, vector.m_y * t, vector.m_z * t);
        }
        public static Vector3f operator *(float t, Vector3f vector)
        {
            return new Vector3f(vector.m_x * t, vector.m_y * t, vector.m_z * t);
        }
        public static Vector3f Multiply(float t, Vector3f vector)
        {
            return new Vector3f(vector.m_x * t, vector.m_y * t, vector.m_z * t);
        }
        public static Vector3f operator /(Vector3f vector, float t)
        {
            return new Vector3f(vector.m_x / t, vector.m_y / t, vector.m_z / t);
        }
        public static Vector3f Divide(Vector3f vector, float t)
        {
            return new Vector3f(vector.m_x / t, vector.m_y / t, vector.m_z / t);
        }
        public static Vector3f operator +(Vector3f vector1, Vector3f vector2)
        {
            return new Vector3f(vector1.m_x + vector2.m_x, vector1.m_y + vector2.m_y, vector1.m_z + vector2.m_z);
        }
        public static Vector3f Add(Vector3f vector1, Vector3f vector2)
        {
            return new Vector3f(vector1.m_x + vector2.m_x, vector1.m_y + vector2.m_y, vector1.m_z + vector2.m_z);
        }
        public static Vector3f operator -(Vector3f vector1, Vector3f vector2)
        {
            return new Vector3f(vector1.m_x - vector2.m_x, vector1.m_y - vector2.m_y, vector1.m_z - vector2.m_z);
        }
        public static Vector3f Subtract(Vector3f vector1, Vector3f vector2)
        {
            return new Vector3f(vector1.m_x - vector2.m_x, vector1.m_y - vector2.m_y, vector1.m_z - vector2.m_z);
        }
        public static float operator *(Vector3f vector1, Vector3f vector2)
        {
            return (vector1.m_x * vector2.m_x + vector1.m_y * vector2.m_y + vector1.m_z * vector2.m_z);
        }
        public static float Multiply(Vector3f vector1, Vector3f vector2)
        {
            return (vector1.m_x * vector2.m_x + vector1.m_y * vector2.m_y + vector1.m_z * vector2.m_z);
        }
        public static float DotProduct(Vector3f vector1, Vector3f vector2)
        {
            return vector1 * vector2;
        }
        public static Vector3f operator -(Vector3f vector)
        {
            return new Vector3f(-vector.m_x, -vector.m_y, -vector.m_z);
        }
        public static Vector3f Negate(Vector3f vector)
        {
            return new Vector3f(-vector.m_x, -vector.m_y, -vector.m_z);
        }
        public static bool operator ==(Vector3f a, Vector3f b)
        {
            return a.m_x == b.m_x && a.m_y == b.m_y && a.m_z == b.m_z;
        }
        public static bool operator !=(Vector3f a, Vector3f b)
        {
            return a.m_x != b.m_x || a.m_y != b.m_y || a.m_z != b.m_z;
        }
        public static Vector3f CrossProduct(Vector3f a, Vector3f b)
        {
            return new Vector3f(a.m_y * b.m_z - b.m_y * a.m_z, a.m_z * b.m_x - b.m_z * a.m_x, a.m_x * b.m_y - b.m_x * a.m_y);
        }
        public static float VectorAngle(Vector3f a, Vector3f b)
        {
            if (!a.Unitize() || !b.Unitize())
                return float.NaN;

            //compute dot product
            float dot = a.m_x * b.m_x + a.m_y * b.m_y + a.m_z * b.m_z;
            // remove any "noise"
            if (dot > 1.0) dot = 1.0f;
            if (dot < -1.0) dot = -1.0f;
            float radians = (float)Math.Acos(dot);
            return radians;
        }
        public override bool Equals(object obj)
        {
            return (obj is Vector3f && this == (Vector3f)obj);
        }
        public override int GetHashCode()
        {
            // MSDN docs recommend XOR'ing the internal values to get a hash code
            return m_x.GetHashCode() ^ m_y.GetHashCode() ^ m_z.GetHashCode();
        }
        #endregion
        public Vector3f(float x, float y, float z)
        {
            m_x = x;
            m_y = y;
            m_z = z;
        }
        public Vector3f(double x, double y, double z)
        {
            m_x = (float)x;
            m_y = (float)y;
            m_z = (float)z;
        }
        public Vector3f(Point3d p)
        {
            m_x = (float)p.m_x;
            m_y = (float)p.m_y;
            m_z = (float)p.m_z;
        }
        public Vector3f(Vector3d v)
        {
            m_x = (float)v.m_x;
            m_y = (float)v.m_y;
            m_z = (float)v.m_z;
        }
        public static Vector3f Zero
        {
            get { return new Vector3f(); }
        }
        public static Vector3f XAxis
        {
            get { return new Vector3f(1.0f, 0.0f, 0.0f); }
        }
        public static Vector3f YAxis
        {
            get { return new Vector3f(0.0f, 1.0f, 0.0f); }
        }
        public static Vector3f ZAxis
        {
            get { return new Vector3f(0.0, 0.0, 1.0); }
        }
        public float Length()
        {
            if (this == Vector3f.Zero) return 0;
            return (float)Math.Sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
        }
        public float LengthSquared()
        {
            return (float)(m_x * m_x + m_y * m_y + m_z * m_z);
        }
        public bool Unitize()
        {
            bool rc = false;
            float d = Length();
            if (d > 0.0)
            {
                m_x = (float)(m_x / d);
                m_y = (float)(m_y / d);
                m_z = (float)(m_z / d);
                rc = true;
            }
            return rc;
        }
    }
    public struct Circle
    {
        #region members
        internal Plane m_plane;
        internal double m_radius;
        #endregion
        #region constructors
        public Circle(double radius) : this(Plane.WorldXY, radius) { }
        public Circle(Plane plane, double radius)
        {
            m_plane = plane;
            m_radius = radius;
        }
        public Circle(Point3d center, double radius)
        {
            m_plane = Plane.WorldXY;
            m_plane.Origin = center;
            m_radius = radius;
        }
        public Circle(Point3d P, Point3d Q, Point3d R)
        {
            m_plane = Plane.WorldXY;
            m_radius = 0.0;
            double a = P.DistanceTo(Q); double b = P.DistanceTo(R); double c = R.DistanceTo(Q);
            if (a < RhinoMath.tol || b < RhinoMath.tol || c < RhinoMath.tol) return;
            if (a >= (b + c) || b >= (a + c) || c >= (b + a)) return;
            m_plane = new Plane(P, Q, R);
            m_radius = a * b * c / Math.Sqrt((a + b + c) * (a + b - c) * (a + c - b) * (b + c - a));
            m_plane.Origin = RhinoMath.SolveCenterPointOfCircle(P, Q, R);
        }
        public static Circle CircleInside(Point3d P, Point3d Q, Point3d R)
        {
            Plane plane = Plane.WorldXY;
            double radius = 0.0;
            double a = R.DistanceTo(Q); double b = P.DistanceTo(R); double c = P.DistanceTo(Q);
            if (a < RhinoMath.tol || b < RhinoMath.tol || c < RhinoMath.tol) return new Circle(plane, radius);
            if (a >= (b + c) || b >= (a + c) || c >= (b + a)) return new Circle(plane, radius);
            plane = new Plane(P, Q, R);
            radius = Math.Sqrt((a + b + c) * (a + b - c) * (a + c - b) * (b + c - a))
                / (2 * (a + b + c));
            plane.Origin = RhinoMath.SolveCenterPointInCircle(P, Q, R, a, b, c);
            return new Circle(plane, radius);
        }
        public Circle(Plane plane, Point3d center, double radius)
        {
            m_plane = plane;
            m_radius = radius;
            m_plane.Origin = center;
        }
        #endregion
        #region properties
        public double Radius
        {
            get { return m_radius; }
            set { m_radius = value; }
        }
        public double Diameter
        {
            get { return m_radius * 2.0; }
            set { m_radius = 0.5 * value; }
        }
        public Plane Plane
        {
            get { return m_plane; }
            set { m_plane = value; }
        }
        public Point3d Center
        {
            get { return m_plane.Origin; }
            set { m_plane.Origin = value; }
        }
        public Vector3d Normal
        {
            get { return m_plane.ZAxis; }
        }
        public double Circumference
        {
            get
            {
                return Math.Abs(2.0 * Math.PI * m_radius);
            }
            set
            {
                m_radius = value / (2.0 * Math.PI);
            }
        }
        private static double Length2d(double x, double y)
        {
            double len;
            x = Math.Abs(x);
            y = Math.Abs(y);
            if (y > x)
            {
                len = x;
                x = y;
                y = len;
            }
            if (x > double.Epsilon)
            {
                len = 1.0 / x;
                y *= len;
                len = x * Math.Sqrt(1.0 + y * y);
            }
            else if (x > 0.0 && !double.IsInfinity(x))
            {
                len = x;
            }
            else
            {
                len = 0.0;
            }
            return len;
        }
        #endregion
        public void Reverse()
        {
            m_plane.YAxis = -m_plane.YAxis;
            m_plane.ZAxis = -m_plane.ZAxis;
        }
        public bool Transform(Transform xform)
        {
            return m_plane.Transform(xform);
        }
    }
}
