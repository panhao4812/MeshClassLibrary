using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Kangaroo
{
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
            double[,] R = new double[3, 6]{
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
            if (Math.Abs(start_dir.Length - 1.0) > 0)
                start_dir.Unitize();
            if (Math.Abs(end_dir.Length - 1.0) > 0)
                end_dir.Unitize();
            double cos_angle = start_dir * end_dir;
            Vector3d axis = Vector3d.CrossProduct(start_dir, end_dir);
            double sin_angle = axis.Length;
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
}
