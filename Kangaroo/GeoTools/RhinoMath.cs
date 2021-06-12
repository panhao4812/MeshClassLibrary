using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeoTools
{
    public class RhinoMath
    {
        public static bool OrthoClose(Point3d Point1, Point3d Point2, double t)
        {
            return (((Math.Abs((double)(Point1.X - Point2.X)) < t) &&
                (Math.Abs((double)(Point1.Y - Point2.Y)) < t)) &&
                (Math.Abs((double)(Point1.Z - Point2.Z)) < t));
        }
        public const double DefaultAngleTolerance = PI / 180.0;
        public const double SqrtEpsilon = 1.490116119385000000e-8;
        public const double UnsetValue = -1.23432101234321e+308;
        public const float UnsetSingle = -1.234321e+38f;
        public const double ZeroTolerance = 1.0e-12;
        const double PI = 3.141592653589793238462643;
        public static double ToRadians(double degrees)
        {
            return degrees * PI / 180.0;
        }
        public static double ToDegrees(double radians)
        {
            return radians * 180.0 / PI;
        }
        public static bool IsValidDouble(double x)
        {
            return (x != UnsetValue) && (!double.IsInfinity(x)) && (!double.IsNaN(x));
        }
        public static double tol = 0.00001;
        public static int Solve2x2(double m00, double m01, double m10, double m11, double d0, double d1,
                      ref double x_addr, ref double y_addr, ref double pivot_ratio)
        {
            int i = 0;
            double maxpiv, minpiv;
            double x = Math.Abs(m00);
            double y = Math.Abs(m01); if (y > x) { x = y; i = 1; }
            y = Math.Abs(m10); if (y > x) { x = y; i = 2; }
            y = Math.Abs(m11); if (y > x) { x = y; i = 3; }
            pivot_ratio = x_addr = y_addr = 0.0;
            if (x == 0.0) { return 0; }//rank=0
            minpiv = maxpiv = x;
            if ((i % 2) != 0)
            {
                { double tmp = x_addr; x_addr = y_addr; y_addr = tmp; }
                x = m00; m00 = m01; m01 = x;
                x = m10; m10 = m11; m11 = x;
            }
            if (i > 1)
            {
                x = d0; d0 = d1; d1 = x;
                x = m00; m00 = m10; m10 = x;
                x = m01; m01 = m11; m11 = x;
            }
            x = 1.0 / m00;
            m01 *= x; d0 *= x;
            if (m10 != 0.0) { m11 -= m10 * m01; d1 -= m10 * d0; }
            if (m11 == 0.0) { return 1; }// rank = 1}
            y = Math.Abs(m11);
            if (y > maxpiv) { maxpiv = y; }
            else if (y < minpiv) { minpiv = y; }
            d1 /= m11;
            if (m01 != 0.0) { d0 -= m01 * d1; }
            x_addr = d0;
            y_addr = d1;
            pivot_ratio = minpiv / maxpiv;
            return 2;
        }
        public static Point3d SolveCenterPointOfCircle(Point3d p1, Point3d p2, Point3d p3)
        {
            double a1, b1, c1, d1;
            double a2, b2, c2, d2;
            double a3, b3, c3, d3;
            double x1 = p1.X, y1 = p1.Y, z1 = p1.Z;
            double x2 = p2.X, y2 = p2.Y, z2 = p2.Z;
            double x3 = p3.X, y3 = p3.Y, z3 = p3.Z;
            a1 = (y1 * z2 - y2 * z1 - (y1 * z3 - y3 * z1) + y2 * z3 - y3 * z2);
            b1 = -(x1 * z2 - x2 * z1 - (x1 * z3 - x3 * z1) + x2 * z3 - x3 * z2);
            c1 = (x1 * y2 - x2 * y1 - (x1 * y3 - x3 * y1) + x2 * y3 - x3 * y2);
            d1 = -(x1 * y2 * z3 - x1 * y3 * z2 - (x2 * y1 * z3 - x2 * y3 * z1) + x3 * y1 * z2 - x3 * y2 * z1);
            a2 = 2 * (x2 - x1);
            b2 = 2 * (y2 - y1);
            c2 = 2 * (z2 - z1);
            d2 = x1 * x1 + y1 * y1 + z1 * z1 - (x2 * x2 + y2 * y2 + z2 * z2);
            a3 = 2 * (x3 - x1);
            b3 = 2 * (y3 - y1);
            c3 = 2 * (z3 - z1);
            d3 = x1 * x1 + y1 * y1 + z1 * z1 - (x3 * x3 + y3 * y3 + z3 * z3);
            double t = (a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 + a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1);
            double x = -(b1 * c2 * d3 - b1 * c3 * d2 - b2 * c1 * d3 + b2 * c3 * d1 + b3 * c1 * d2 - b3 * c2 * d1) / t;
            double y = (a1 * c2 * d3 - a1 * c3 * d2 - a2 * c1 * d3 + a2 * c3 * d1 + a3 * c1 * d2 - a3 * c2 * d1) / t;
            double z = -(a1 * b2 * d3 - a1 * b3 * d2 - a2 * b1 * d3 + a2 * b3 * d1 + a3 * b1 * d2 - a3 * b2 * d1) / t;
            return new Point3d(x, y, z);
        }
        public static Point3d SolveCenterPointInCircle(Point3d p1, Point3d p2, Point3d p3,
            double a, double b, double c)
        {
            double x = (a * p1.X + b * p2.X + c * p3.X) / (a + b + c);
            double y = (a * p1.Y + b * p2.Y + c * p3.Y) / (a + b + c);
            double z = (a * p1.Z + b * p2.Z + c * p3.Z) / (a + b + c);
            return new Point3d(x, y, z);
        }
        public static Point3d SolveCenterPointInCircle(Point3d p1, Point3d p2, Point3d p3)
        {
            double a = p2.DistanceTo(p3); double b = p1.DistanceTo(p3); double c = p1.DistanceTo(p2);
            return SolveCenterPointInCircle(p1, p2, p3, a, b, c);
        }

    }
}
