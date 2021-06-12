using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeoTools
{
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
