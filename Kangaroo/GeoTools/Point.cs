using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeoTools
{
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
        public bool IsValid
        {
            get { return RhinoMath.IsValidDouble(m_x) && RhinoMath.IsValidDouble(m_y) && RhinoMath.IsValidDouble(m_z); }
        }
        public static Point3d Unset
        {
            get { return new Point3d(RhinoMath.UnsetValue, RhinoMath.UnsetValue, RhinoMath.UnsetValue); }
        }
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
        public static implicit operator Point3d(Point3f point)
        {
            return new Point3d(point);
        }
        public static implicit operator Point3d(Vector3d point)
        {
            return new Point3d(point);
        }
    }
    public struct Vector3d
    {
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
        public static double VectorAngle(Vector3d a, Vector3d b, Plane plane)
        {
            { // Project vectors onto plane.
                Point3d pA = plane.Origin + a;
                Point3d pB = plane.Origin + b;

                pA = plane.ClosestPoint(pA);
                pB = plane.ClosestPoint(pB);

                a = pA - plane.Origin;
                b = pB - plane.Origin;
            }

            // Abort on invalid cases.
            if (!a.Unitize()) { return RhinoMath.UnsetValue; }
            if (!b.Unitize()) { return RhinoMath.UnsetValue; }

            double dot = a * b;
            { // Limit dot product to valid range.
                if (dot >= 1.0)
                { dot = 1.0; }
                else if (dot < -1.0)
                { dot = -1.0; }
            }

            double angle = Math.Acos(dot);
            { // Special case (anti)parallel vectors.
                if (Math.Abs(angle) < 1e-64) { return 0.0; }
                if (Math.Abs(angle - Math.PI) < 1e-64) { return Math.PI; }
            }

            Vector3d cross = Vector3d.CrossProduct(a, b);
            if (plane.ZAxis.IsParallelTo(cross) == +1)
                return angle;
            return 2.0 * Math.PI - angle;
        }
        public bool IsValid
        {
            get { return RhinoMath.IsValidDouble(m_x) && RhinoMath.IsValidDouble(m_y) && RhinoMath.IsValidDouble(m_z); }
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
        public static Vector3d Unset
        {
            get { return new Vector3d(RhinoMath.UnsetValue, RhinoMath.UnsetValue, RhinoMath.UnsetValue); }
        }
        public static implicit operator Vector3d(Point3d vector)
        {
            return new Vector3d(vector);
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
        public double SquareLength
        {
            get { return (m_x * m_x) + (m_y * m_y) + (m_z * m_z); }
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
            get
            {
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
            double d = Length;
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
            double ll = Length * v.Length;
            if (ll > 0.0)
            {
                if (Math.Abs((X * v.X + Y * v.Y + Z * v.Z) / ll) <= Math.Sin(angle_tolerance))
                    rc = true;
            }
            return rc;
        }
        public int IsParallelTo(Vector3d other)
        {
            return IsParallelTo(other, RhinoMath.DefaultAngleTolerance);
        }
        public int IsParallelTo(Vector3d v, double angle_tolerance)
        {
            int rc = 0;
            double ll = Length * v.Length;
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
            double s0 = 1.0 / V0.Length;
            double s1 = 1.0 / V1.Length;
            double s2 = 1.0 / V2.Length;
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
        public bool IsZero
        {
            get
            {
                return (m_x == 0.0 && m_y == 0.0 && m_z == 0.0);
            }
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
        public Point3f(Vector3d p)
        {
            m_x = (float)p.m_x;
            m_y = (float)p.m_y;
            m_z = (float)p.m_z;
        }
        public bool IsValid
        {
            get { return RhinoMath.IsValidDouble(m_x) && RhinoMath.IsValidDouble(m_y) && RhinoMath.IsValidDouble(m_z); }
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
        public static implicit operator Point3f(Point3d point)
        {
            return new Point3f(point);
        }
        public static implicit operator Point3f(Vector3d point)
        {
            return new Point3f(point);
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
}
