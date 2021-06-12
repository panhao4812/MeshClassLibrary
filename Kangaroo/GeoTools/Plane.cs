using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeoTools
{
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
            double d = m_zaxis.Length;
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
        public static bool IsOrthogonalFrame(Vector3d X, Vector3d Y, Vector3d Z)
        {
            // returns true if X, Y, Z is an orthogonal frame
            if (!X.IsValid || !Y.IsValid || !Z.IsValid)
                return false;

            double lx = X.Length;
            double ly = Y.Length;
            double lz = Z.Length;
            if (lx <= RhinoMath.SqrtEpsilon)
                return false;
            if (ly <= RhinoMath.SqrtEpsilon)
                return false;
            if (lz <= RhinoMath.SqrtEpsilon)
                return false;
            lx = 1.0 / lx;
            ly = 1.0 / ly;
            lz = 1.0 / lz;
            double xy = (X.X * Y.X + X.Y * Y.Y + X.Z * Y.Z) * lx * ly;
            double yz = (Y.X * Z.X + Y.Y * Z.Y + Y.Z * Z.Z) * ly * lz;
            double zx = (Z.X * X.X + Z.Y * X.Y + Z.Z * X.Z) * lz * lx;
            if (Math.Abs(xy) > RhinoMath.SqrtEpsilon
                 || Math.Abs(yz) > RhinoMath.SqrtEpsilon
                 || Math.Abs(zx) > RhinoMath.SqrtEpsilon
               )
            {
                double t = 0.0000152587890625;
                if (Math.Abs(xy) >= t || Math.Abs(yz) >= t || Math.Abs(zx) >= t)
                    return false;
                Vector3d V;
                V = (lx * ly) * Vector3d.CrossProduct(X, Y);
                t = Math.Abs((V.X * Z.X + V.Y * Z.Y + V.Z * Z.Z) * lz);
                if (Math.Abs(t - 1.0) > RhinoMath.SqrtEpsilon)
                    return false;

                V = (ly * lz) * Vector3d.CrossProduct(Y, Z);
                t = Math.Abs((V.X * X.X + V.Y * X.Y + V.Z * X.Z) * lx);
                if (Math.Abs(t - 1.0) > RhinoMath.SqrtEpsilon)
                    return false;

                V = (lz * lx) * Vector3d.CrossProduct(Z, X);
                t = Math.Abs((V.X * Y.X + V.Y * Y.Y + V.Z * Y.Z) * ly);
                if (Math.Abs(t - 1.0) > RhinoMath.SqrtEpsilon)
                    return false;
            }
            return true;
        }
        public static bool IsOrthonormalFrame(Vector3d X, Vector3d Y, Vector3d Z)
        {
            // returns true if X, Y, Z is an orthonormal frame
            if (!IsOrthogonalFrame(X, Y, Z))
                return false;
            double x = X.Length;
            if (Math.Abs(x - 1.0) > RhinoMath.SqrtEpsilon)
                return false;
            x = Y.Length;
            if (Math.Abs(x - 1.0) > RhinoMath.SqrtEpsilon)
                return false;
            x = Z.Length;
            if (Math.Abs(x - 1.0) > RhinoMath.SqrtEpsilon)
                return false;

            return true;
        }
        public static bool IsRightHandFrame(Vector3d X, Vector3d Y, Vector3d Z)
        {
            // returns true if X, Y, Z is an orthonormal right hand frame
            if (!IsOrthonormalFrame(X, Y, Z))
                return false;
            double x = Vector3d.DotProduct(Vector3d.CrossProduct(X, Y), Z);
            if (x <= RhinoMath.SqrtEpsilon)
                return false;
            return true;
        }
        public bool IsValid
        {
            get
            {
                double[] equation = GetPlaneEquation();

                if (!
                    (RhinoMath.IsValidDouble(equation[0]) && RhinoMath.IsValidDouble(equation[1]) &&
                    RhinoMath.IsValidDouble(equation[2]) && RhinoMath.IsValidDouble(equation[3]))
                    ) return false;
                if (!IsRightHandFrame(XAxis, YAxis, ZAxis))
                    return false;
                return true;
            }
        }
        public bool Rotate(double angle, Vector3d axis)
        {
            return Rotate(Math.Sin(angle), Math.Cos(angle), axis);
        }
        public bool Rotate(double sinAngle, double cosAngle, Vector3d axis)
        {
            bool rc = true;
            if (axis == ZAxis)
            {
                Vector3d x = cosAngle * XAxis + sinAngle * YAxis;
                Vector3d y = cosAngle * YAxis - sinAngle * XAxis;
                XAxis = x;
                YAxis = y;
            }
            else
            {
                Point3d origin_pt = Origin;
                rc = Rotate(sinAngle, cosAngle, axis, Origin);
                Origin = origin_pt; // to kill any fuzz
            }
            return rc;
        }
        public bool Rotate(double angle, Vector3d axis, Point3d centerOfRotation)
        {
            return Rotate(Math.Sin(angle), Math.Cos(angle), axis, centerOfRotation);
        }
        public bool Rotate(double sinAngle, double cosAngle, Vector3d axis, Point3d centerOfRotation)
        {
            if (centerOfRotation == Origin)
            {
                Transform rot = GeoTools.Transform.Rotation(sinAngle, cosAngle, axis, Point3d.Origin);
                XAxis = rot * XAxis;
                YAxis = rot * YAxis;
                ZAxis = rot * ZAxis;
                return true;
            }
            Transform rot2 = GeoTools.Transform.Rotation(sinAngle, cosAngle, axis, centerOfRotation);
            return Transform(rot2);
        }
        public static bool FitPlaneToPoints(List<Point3d> points, out Plane plane)
        {
            Line fitLine;
            GeoSolver.LinePlaneEstimate(points, out fitLine, out plane);
            if (plane.IsValid) return true;
            return false;
        }
        public static bool FitPlaneToPoints(Point3d[] points, out Plane plane)
        {
            return FitPlaneToPoints(new List<Point3d>(points), out plane);       
        }

    }
}
