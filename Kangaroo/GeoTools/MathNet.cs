using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Double.Factorization;
using MathNet.Numerics.LinearAlgebra.Generic;
using System;
using System.Collections.Generic;

namespace Kangaroo
{
    public class GeoSolver
    {
        public GeoSolver() { }
        public static void LinePlaneEstimate(List<Point3d> datas, out Line line, out Plane plane)
        {
            Point3d cen = Point3d.GetCenter(datas);
            DenseMatrix jacobian = new DenseMatrix(datas.Count, 3);
            foreach (Point3d temp in datas)
            {
                Vector<double> gradient = new DenseVector(3);
                gradient[0] = temp.X - cen.X;
                gradient[1] = temp.Y - cen.Y;
                gradient[2] = temp.Z - cen.Z;
                jacobian.SetRow(datas.IndexOf(temp), gradient);
            }
            Svd svd = jacobian.Svd(true);
            Matrix<double> V = svd.VT().Transpose();
            Vector<double> parameters = V.Column(2);
            Vector<double> parameters2 = V.Column(0);
            plane = new Plane(cen, new Vector3d(parameters[0], parameters[1], parameters[2]));
            line = new Line(cen, cen + new Vector3d(parameters2[0], parameters2[1], parameters2[2]));
        }
        public static Transform Kabsch(List<Point3d> P, List<Point3d> Q)
        {
            Transform xf = Transform.Identity;
            if (P.Count != Q.Count) return xf;
            Point3d CenP = new Point3d(), CenQ = new Point3d();
            for (int i = 0; i < P.Count; i++)
            {
                CenP += P[i];
                CenQ += Q[i];
            }
            CenP /= P.Count; CenQ /= Q.Count;
            DenseMatrix MX = new DenseMatrix(P.Count, 3);
            DenseMatrix MY = new DenseMatrix(P.Count, 3);
            for (int i = 0; i < P.Count; i++)
            {
                DenseVector v1 = new DenseVector(3);
                DenseVector v2 = new DenseVector(3);
                v1[0] = P[i].X - CenP.X; v2[0] = Q[i].X - CenQ.X;
                v1[1] = P[i].Y - CenP.Y; v2[1] = Q[i].Y - CenQ.Y;
                v1[2] = P[i].Z - CenP.Z; v2[2] = Q[i].Z - CenQ.Z;
                MX.SetRow(i, v1); MY.SetRow(i, v2);
            }
            DenseMatrix H = DenseMatrix.OfMatrix(MX.TransposeThisAndMultiply(MY));
            Svd svd = H.Svd(true);
            Matrix<double> UT = svd.U().Transpose();
            Matrix<double> V = svd.VT().Transpose();
            Matrix<double> R = V.Multiply(UT);
            double d = R.Determinant();
            if (d > 0) { d = 1; } else { d = -1; }
            DenseMatrix I = new DenseMatrix(3, 3, new double[] { 1, 0, 0, 0, 1, 0, 0, 0, d });
            R = V.Multiply(I).Multiply(UT);
            xf.M00 = R[0, 0]; xf.M01 = R[0, 1]; xf.M02 = R[0, 2];
            xf.M10 = R[1, 0]; xf.M11 = R[1, 1]; xf.M12 = R[1, 2];
            xf.M20 = R[2, 0]; xf.M21 = R[2, 1]; xf.M22 = R[2, 2];
            CenP.Transform(xf);
            Transform tf = Transform.Translation(CenQ - CenP);
            return Transform.Multiply(tf, xf);
        }
        public static Point3d bezier3func(double uu, List<Point3d> ControlPoint)
        {
            //...p1,p1_right,p2_left,p2...
            if (ControlPoint.Count != 4)
            {
                return new Point3d();
            }
            Point3d part1 = ControlPoint[0] * uu * uu * uu;
            Point3d part2 = ControlPoint[1] * uu * uu * (1 - uu);
            Point3d part3 = ControlPoint[2] * uu * (1 - uu) * (1 - uu);
            Point3d part4 = ControlPoint[3] * (1 - uu) * (1 - uu) * (1 - uu);
            return part1 + part2 + part3 + part4;
        }
        public static Point3d bezier3func(double uu, Point3d[] ControlPoint)
        {
            //...p1,p1_right,p2_left,p2...
            if (ControlPoint.Length != 4)
            {
                return new Point3d();
            }
            Point3d part1 = ControlPoint[0] * uu * uu * uu;
            Point3d part2 = ControlPoint[1] * uu * uu * (1 - uu);
            Point3d part3 = ControlPoint[2] * uu * (1 - uu) * (1 - uu);
            Point3d part4 = ControlPoint[3] * (1 - uu) * (1 - uu) * (1 - uu);
            return part1 + part2 + part3 + part4;
        }
        public static List<Point3d> Bezier(List<Point3d> OriginPoint, bool IsClosed)
        {
            return Bezier(OriginPoint, IsClosed, 1, 0.025);
        }
        public static List<Point3d> Bezier(List<Point3d> OriginPoint, bool IsClosed, double d, double uu)
        {
            List<Point3d> output = new List<Point3d>();
            if (OriginPoint.Count < 3) return output;
            Point3d FirstPoint; Point3d LastPoint;
            if (IsClosed)
            {
                FirstPoint = OriginPoint[OriginPoint.Count - 1];
                LastPoint = OriginPoint[0];
            }
            else
            {
                Vector3d v = OriginPoint[0] - OriginPoint[1];
                FirstPoint = OriginPoint[0] + v;
                v = OriginPoint[OriginPoint.Count - 1] - OriginPoint[OriginPoint.Count - 2];
                LastPoint = OriginPoint[OriginPoint.Count - 1] + v;
            }
            OriginPoint.Add(LastPoint); OriginPoint.Insert(0, FirstPoint);
            List<Point3d> LeftPoint = new List<Point3d>();
            List<Point3d> RightPoint = new List<Point3d>();
            LeftPoint.Add(new Point3d()); RightPoint.Add(new Point3d());
            for (int i = 1; i < OriginPoint.Count - 1; i++)
            {
                Point3d leftpt = (OriginPoint[i] + OriginPoint[i - 1]) / 2;
                Point3d rightpt = (OriginPoint[i] + OriginPoint[i + 1]) / 2;
                Vector3d v1 = OriginPoint[i] - OriginPoint[i - 1];
                Vector3d v2 = OriginPoint[i] - OriginPoint[i + 1];
                Vector3d v3 = rightpt - leftpt;
                v3 *= v1.Length() / (v1.Length() + v2.Length());
                Point3d cenpt = leftpt + v3;
                v1 = leftpt - cenpt; v2 = rightpt - cenpt; v1 *= d; v2 *= d;
                leftpt = cenpt + v1; rightpt = cenpt + v2;
                Vector3d v4 = OriginPoint[i] - cenpt;
                leftpt += v4; rightpt += v4;
                LeftPoint.Add(leftpt); RightPoint.Add(rightpt);
            }
            double indexlast = 2;
            if (IsClosed) { indexlast = 1; LeftPoint.Add(LeftPoint[1]); }

            for (int i = 1; i < OriginPoint.Count - indexlast; i++)
            {
                Point3d[] controlP = new Point3d[4];
                controlP[0] = OriginPoint[i];
                controlP[1] = RightPoint[i];
                controlP[2] = LeftPoint[i + 1];
                controlP[3] = OriginPoint[i + 1];
                double u = 1;
                while (u >= 0)
                {
                    output.Add(bezier3func(u, controlP));
                    u -= uu;
                }
            }
            return output;
        }
    }
    public class CylinderSolver
    {
        public CylinderSolver() { }
        public static double IterationEndValue = 0.00001;
        public Plane CylinderPlane;
        public double CylinderRadius;
        public Line GetLine(double length)
        {
            Point3d pt = CylinderPlane.Origin;
            Vector3d v = CylinderPlane.Normal;
            return new Line(pt, pt + v * length);
        }
        public void init_paras(List<Point3d> pts, Line L, out Plane plane, out double R)
        {
            List<double> paras = new List<double>();
            double x0, y0, z0, a, b, c;
            Vector3d v1 = L.To - L.From; v1.Unitize();
            x0 = L.FromX; y0 = L.FromY; z0 = L.FromZ;
            a = v1.X; b = v1.Y; c = v1.Z; R = 0;
            foreach (Point3d temp in pts)
            {
                double x = temp.X;
                double y = temp.Y;
                double z = temp.Z;
                double u = c * (y - y0) - b * (z - z0);
                double v = a * (z - z0) - c * (x - x0);
                double w = b * (x - x0) - a * (y - y0);
                R += Math.Sqrt(u * u + v * v + w * w);
            }
            R /= pts.Count;
            plane = new Plane(new Point3d(x0, y0, z0), new Vector3d(a, b, c));
        }
        private double getDeltaValue(Plane p, double radius, List<Point3d> datas)
        {
            double x0 = p.OriginX;
            double y0 = p.OriginY;
            double z0 = p.OriginZ;
            double a = p.Normal.X;
            double b = p.Normal.Y;
            double c = p.Normal.Z;
            double result = 0;
            foreach (Point3d temp in datas)
            {
                double u, v, w, x, y, z, r;
                x = temp.X; y = temp.Y; z = temp.Z;
                u = c * (y - y0) - b * (z - z0);
                v = a * (z - z0) - c * (x - x0);
                w = b * (x - x0) - a * (y - y0);
                r = Math.Sqrt(u * u + v * v + w * w) / Math.Sqrt(a * a + b * b + c * c);
                result += r - radius;
            }
            return result / datas.Count;
        }
        public void Estimate(List<Point3d> Rhinodatas, Line L)
        {
            if (Rhinodatas.Count < 6) return;
            Plane plane_old, plane_new;
            double radius_old, radius_new;
            init_paras(Rhinodatas, L, out plane_old, out radius_old);
            double deltavalue_old = getDeltaValue(plane_old, radius_old, Rhinodatas);
            double deltavalue_new = 0;
            int step = 0;
            while (step < 10000)
            {
                step++;
                //先将圆柱平面移到worldXY（将7个参数简化到5个-->P矩阵）
                List<Point3d> pos = new List<Point3d>();
                foreach (Point3d tempPoint in Rhinodatas)
                {
                    tempPoint.Transform(Transform.PlaneToPlane(plane_old, Plane.WorldXY));
                    pos.Add(tempPoint);
                }
                //构造jacobian矩阵 
                DenseMatrix J = new DenseMatrix(pos.Count, 5);
                Vector<double> d = new DenseVector(pos.Count);
                for (int i = 0; i < pos.Count; i++)
                {
                    double x, y, z, r;
                    x = pos[i].X;
                    y = pos[i].Y;
                    z = pos[i].Z;
                    r = Math.Sqrt(x * x + y * y);
                    Vector<double> single = new DenseVector(
                        new double[] { -x / r, -y / r, -x * z / r, -y * z / r, -1 });
                    J.SetRow(i, single);
                    d[i] = -(r - radius_old);
                }
                //解J*P=-d
                Vector<double> P = J.Transpose().Multiply(J).Cholesky().Solve(J.Transpose().Multiply(d));
                //还原位置
                plane_new = new Plane(
                    new Point3d(P[0], P[1], -P[0] * P[2] - P[1] * P[3]),
                    new Vector3d(P[2], P[3], 1));
                plane_new.Transform(Transform.PlaneToPlane(Plane.WorldXY, plane_old));
                radius_new = radius_old + P[4];
                //终止判断
                deltavalue_new = getDeltaValue(plane_new, radius_new, Rhinodatas);
                if (Math.Abs(deltavalue_new - deltavalue_old) < IterationEndValue) { break; }
                else
                {
                    plane_old = plane_new;
                    radius_old = radius_new;
                    deltavalue_old = deltavalue_new;
                }
            }
            this.CylinderPlane = plane_old;
            this.CylinderRadius = radius_old;
        }
    }
    public class NumSolver{
    public NumSolver() { }
        public static double[] LinearSolver(double[,] A, double[]B)
        {
            //double[,] A = { { 8, -3, 2 },{4,11,-1 },{2,1,4 } };
           // double[] B = new double[] { 20, 33, 12 };
            Matrix<double> ma = DenseMatrix.OfArray(A); 
            Vector<double> d = new DenseVector(B);
            Vector<double>output= ma.LU().Solve(d);
            return output.ToArray();
    }
        public static double[] NewtonSolver(double[,] A, double[] B)
        {
            ////解J*P=-d
            //double[,] A = { { -10,0 },{1,-10 } };
            // double[] B = new double[] {8,8};
            Matrix<double> J = DenseMatrix.OfArray(A);
            Vector<double> d = new DenseVector(B);
            Vector<double> P=J.Inverse().Multiply(d).Multiply(-1);
            return P.ToArray();
        }

    }


}
