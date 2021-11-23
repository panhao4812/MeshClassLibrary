using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Generic;
using MathNet.Numerics.LinearAlgebra.Double.Factorization;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace FitAndInterpolation
{
    public class RandomPointEstimate
    {
        public RandomPointEstimate() { }
        public List<Point3d> GetRandomVertex(Mesh mesh, int count)
        {
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            List<Point3d> output = new List<Point3d>();
            if (count < 7) return output;
            BoundingBox box = mesh.GetBoundingBox(true);
            double distance = box.Max.DistanceTo(box.Min);
            distance /= Math.Sqrt(count * 2); int stepCount = 0;
            while (output.Count < count && stepCount < count * 100)
            {
                bool sign = true;
                Point3d pt2 = new Point3d(vs[rnd.Next(0, vs.Count)]);
                for (int i = 0; i < output.Count; i++)
                {
                    Point3d pt = output[i];
                    if (pt.DistanceTo(pt2) < distance)
                    {
                        sign = false;
                        break;
                    }
                }
                if (sign) { output.Add(pt2); stepCount = 0; } else { stepCount++; }
            }
            return output;
        }
        public List<Point3d> ShakePoints(List<Point3d> Pts, double radius)
        {
            List<Point3d> pts = new List<Point3d>();
            for (int i = 0; i < Pts.Count; i++)
            {
                Point3d pt = new Point3d(Pts[i].X + radius * (rnd.NextDouble() * 2 - 1),
                    Pts[i].Y + radius * (rnd.NextDouble() * 2 - 1),
                    Pts[i].Z + radius * (rnd.NextDouble() * 2 - 1));
                pts.Add(pt);
            }
            return pts;
        }
        Random rnd = new Random();
    }
    public class CylinderSolver
    {
        public CylinderSolver()
        {
            CylinderPlane = Plane.WorldXY;
            CylinderRadius = 0;
        }
        public Plane CylinderPlane;
        public double CylinderRadius;
        public Line GetLine(double length)
        {
            Point3d pt = CylinderPlane.Origin;
            Vector3d v = CylinderPlane.Normal;
            v.Unitize();
            return new Line(pt, pt + v * length);
        }
        private double getDeltaValue(Plane p, double radius, List<Point3d> datas)
        {
            double x0 = p.OriginX, y0 = p.OriginY, z0 = p.OriginZ;
            double a = p.Normal.X, b = p.Normal.Y, c = p.Normal.Z;
            double result = 0;
            foreach (Point3d temp in datas)
            {
                double u, v, w, x, y, z, r;
                x = temp.X; y = temp.Y; z = temp.Z;
                u = c * (y - y0) - b * (z - z0);
                v = a * (z - z0) - c * (x - x0);
                w = b * (x - x0) - a * (y - y0);
                r = Math.Sqrt(Math.Pow(u, 2) + Math.Pow(v, 2) + Math.Pow(w, 2))
                    /
                    Math.Sqrt(Math.Pow(a, 2) + Math.Pow(b, 2) + Math.Pow(c, 2));
                result += r - radius;
            }
            return result / datas.Count;
        }
        public void init_paras(List<Point3d> pts, Line L)
        {
            double x0, y0, z0, a, b, c;
            Vector3d v1 = L.To - L.From; v1.Unitize();
            x0 = L.FromX; y0 = L.FromY; z0 = L.FromZ;
            a = v1.X; b = v1.Y; c = v1.Z;
            CylinderRadius = 0;
            foreach (Point3d temp in pts)
            {
                double x, y, z, u, v, w;
                x = temp.X; y = temp.Y; z = temp.Z;
                u = c * (y - y0) - b * (z - z0);
                v = a * (z - z0) - c * (x - x0);
                w = b * (x - x0) - a * (y - y0);
                CylinderRadius += Math.Sqrt(u * u + v * v + w * w);
            }
            CylinderRadius /= pts.Count;
            CylinderPlane = new Plane(new Point3d(x0, y0, z0), new Vector3d(a, b, c));
        }
        public void Estimate(List<Point3d> datas, Line L)
        {
            if (datas.Count < 6) return;
            Plane plane_new; double radius_new;
            init_paras(datas, L);
            double delt_old = getDeltaValue(CylinderPlane, CylinderRadius, datas);
            double delt_new = 0;
            int step = 0;
            while (step < 10000)
            {
                step++;
                List<Point3d> pos = new List<Point3d>();
                foreach (Point3d tempPoint in datas)
                {
                    tempPoint.Transform(Transform.PlaneToPlane(CylinderPlane, Plane.WorldXY));
                    pos.Add(tempPoint);
                }
                DenseMatrix J = new DenseMatrix(pos.Count, 5);
                Vector<double> dev = new DenseVector(pos.Count);
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
                    dev[i] = -(r - CylinderRadius);
                }
                Vector<double> P = J.Transpose().Multiply(J).Cholesky().Solve(J.Transpose().Multiply(dev));
                plane_new = new Plane(new Point3d(P[0], P[1], -P[0] * P[2] - P[1] * P[3]),
                    new Vector3d(P[2], P[3], 1));
                plane_new.Transform(Transform.PlaneToPlane(Plane.WorldXY, CylinderPlane));
                radius_new = CylinderRadius + P[4];
                delt_new = getDeltaValue(plane_new, radius_new, datas);
                if (Math.Abs(delt_new - delt_old) < 0.00001)
                {
                    break;
                }
                else
                {
                    CylinderPlane = plane_new;
                    CylinderRadius = radius_new;
                }
            }
        }
    }
    public class GeoSolver
    {
        public Transform KabschEstimate(List<Plane> P, List<Plane> Q, double NormalLength)
        {
            Transform output = Transform.Identity;
            if (P.Count != Q.Count) return output;
            List<Point3d> P1 = new List<Point3d>();
            List<Point3d> Q1 = new List<Point3d>();
            for (int i = 0; i < P.Count; i++)
            {
                P1.Add(new Point3d(P[i].Origin));
                Vector3d v1 = P[i].Normal;
                v1.Unitize(); v1 *= NormalLength;
                P1.Add(P[i].Origin + v1); P1.Add(P[i].Origin - v1);
                Q1.Add(new Point3d(Q[i].Origin));
                Vector3d v2 = Q[i].Normal;
                v2.Unitize(); v2 *= NormalLength;
                Q1.Add(Q[i].Origin + v1); Q1.Add(Q[i].Origin - v1);
            }
            return KabschEstimate(P1, Q1);
        }
        public Transform KabschEstimate(List<Plane> P, List<Plane> Q, List<Vector3d> V1,List<Vector3d> V2)
        {
            Transform output = Transform.Identity;
            if (P.Count != Q.Count) return output;
            List<Point3d> P1 = new List<Point3d>();
            Vector3d vx,vy,vz;
            List<Point3d> Q1 = new List<Point3d>();
            for (int i = 0; i < P.Count; i++)
            {
                P1.Add(new Point3d(P[i].Origin));
                if (V1[i].X != 0) { vx = P[i].XAxis; vx.Unitize(); vx *= V1[i].X; P1.Add(P[i].Origin + vx); P1.Add(P[i].Origin - vx); }
                if (V1[i].Y != 0) { vy = P[i].YAxis; vy.Unitize(); vy *= V1[i].Y; P1.Add(P[i].Origin + vy); P1.Add(P[i].Origin - vy); }
                if (V1[i].Z != 0) { vz = P[i].ZAxis; vz.Unitize(); vz *= V1[i].Z; P1.Add(P[i].Origin + vz); P1.Add(P[i].Origin - vz); }
                Q1.Add(new Point3d(Q[i].Origin));
                if (V2[i].X != 0) { vx = Q[i].XAxis; vx.Unitize(); vx *= V2[i].X; Q1.Add(Q[i].Origin + vx); Q1.Add(Q[i].Origin - vx); }
                if (V2[i].Y != 0) { vy = Q[i].YAxis; vy.Unitize(); vy *= V2[i].Y; Q1.Add(Q[i].Origin + vy); Q1.Add(Q[i].Origin - vy); }
                if (V2[i].Z != 0) { vz = Q[i].ZAxis; vz.Unitize(); vz *= V2[i].Z; Q1.Add(Q[i].Origin + vz); Q1.Add(Q[i].Origin - vz); }
            }
            if (P1.Count != Q1.Count) return output;
            return KabschEstimate(P1, Q1);
        }
        public double RMS(List<Point3d> P, List<Point3d> Q)
        {
            double t = 0;
            for (int i = 0; i < P.Count; i++)
            {
                t += P[i].DistanceTo(Q[i]);
            }
            t /= P.Count;
            return Math.Sqrt(t);
        }
        public Transform KabschEstimate(List<Point3d> x, List<Line> y, double RMStol, bool limit)
        {
            List<Point3d> x1 = new List<Point3d>();
            List<Point3d> x2 = new List<Point3d>();
            Transform xform = Transform.Identity;
            for (int j = 0; j < y.Count; j++)
            {
                Line l = y[j];
                Point3d pt2 = l.ClosestPoint(x[j], limit);
                x2.Add(pt2); x1.Add(x[j]);
            }
            for (int i = 0; i < 10000; i++)
            {
                xform = KabschEstimate(x1, x2);

                for (int j = 0; j < x1.Count; j++)
                {
                    Point3d pt = x1[j];
                    pt.Transform(xform);
                    x1[j] = pt;
                }
                for (int j = 0; j < x1.Count; j++)
                {
                    Line l = y[j];
                    Point3d pt2 = l.ClosestPoint(x1[j], limit);
                    x2[j] = pt2;
                }
                // Print(i.ToString());
                if (RMS(x1, x2) < RMStol) break;
            }
            return KabschEstimate(x, x1);
        }
        public Transform KabschEstimate(List<Line> x, List<Line> y, double RMStol, bool limit)
        {
            List<Point3d> x1 = new List<Point3d>();
            List<Point3d> x2 = new List<Point3d>();
            List<Point3d> x0 = new List<Point3d>();
            List<Line> line0 = new List<Line>();
            for (int j = 0; j < x.Count; j++)
            {
                Line l = y[j];
                Point3d pt2 = l.ClosestPoint(x[j].From, limit);
                x1.Add(x[j].From); x2.Add(pt2); line0.Add(l);
                Point3d pt3 = l.ClosestPoint(x[j].To, limit);
                x1.Add(x[j].To); x2.Add(pt3); line0.Add(l);
                if (x.Count < 3)
                {
                    Point3d pt4 = l.ClosestPoint((x[j].To + x[j].From) / 2, limit);
                    x1.Add((x[j].To + x[j].From) / 2); x2.Add(pt4); line0.Add(l);
                }
            }
            x0.AddRange(x1);
            Transform xform = Transform.Identity;
            for (int i = 0; i < 10000; i++)
            {
                xform = KabschEstimate(x1, x2);
                for (int j = 0; j < x1.Count; j++)
                {
                    Point3d pt = x1[j];
                    pt.Transform(xform);
                    x1[j] = pt;
                }
                for (int j = 0; j < x1.Count; j++)
                {
                    Line l = line0[j];
                    Point3d pt2 = l.ClosestPoint(x1[j], limit);
                    x2[j] = pt2;
                }
                if (RMS(x1, x2) < RMStol) break;
            }
            return KabschEstimate(x0, x1);
        }
        public Transform KabschEstimate(List<Plane> P, List<Plane> Q, double RMStol, bool limit, double NormalLength)
        {
            List<Line> ls1 = new List<Line>();
            List<Line> ls2 = new List<Line>();
            for (int i = 0; i < P.Count; i++)
            {
                ls1.Add(new Line(P[i].Origin - P[i].Normal * NormalLength, P[i].Origin + P[i].Normal * NormalLength));
                ls2.Add(new Line(Q[i].Origin - Q[i].Normal * NormalLength, Q[i].Origin + Q[i].Normal * NormalLength));
            }
            return KabschEstimate(ls1, ls2, RMStol, limit);
        }
        public static Point3d GetCenter(List<Point3d> pts)
        {
            Point3d cen = new Point3d();
            foreach (Point3d pt in pts)
            {
                cen += pt;
            }
            cen /= pts.Count;
            return cen;
        }
        public GeoSolver() { }
        public void LinePlaneEstimate(List<Point3d> datas, out Line line, out Plane plane)
        {
            Point3d cen = GetCenter(datas);
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
            Vector<double> para1 = new DenseVector(3);
            Vector<double> para2 = new DenseVector(3);
            para1 = V.Column(0); para2 = V.Column(2);
            plane = new Plane(cen, new Vector3d(para2[0], para2[1], para2[2]));
            line = new Line(cen, cen + new Vector3d(para1[0], para1[1], para1[2]));
        }
        public Transform KabschEstimate(List<Point3d> P, List<Point3d> Q)
        {
            /* 
            一组空间点匹配另一组空间点 kabsch算法
            将两组点的中心算出来 ，然后平移中点到原点
            协方差矩阵 H = X * Y.transpose();
            [U,S,VT]=SVD(H)
            I=(V*UT).Determinant判断方向
            旋转R = V*I*UT;
            平移T = -R * cenA + cenB;
            新点坐标P=P*T*R
            */
            Transform xf = Transform.Identity;
            if (P.Count != Q.Count) return xf;
            Point3d CenP = new Point3d(); Point3d CenQ = new Point3d();
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
        public Point3d bezier3func(double uu, Point3d[] controlP)
        {
            //p1..p1(right)..p2(left)..p2
            Point3d part0 = controlP[0] * uu * uu * uu;
            Point3d part1 = 3 * controlP[1] * uu * uu * (1 - uu);
            Point3d part2 = 3 * controlP[2] * uu * (1 - uu) * (1 - uu);
            Point3d part3 = controlP[3] * (1 - uu) * (1 - uu) * (1 - uu);
            return part0 + part1 + part2 + part3;
        }
        public List<Point3d> BezierEstimate(List<Point3d> OriginPoint, bool IsClosed)
        {
            return BezierEstimate(OriginPoint, IsClosed, 1, 0.025);
        }
        public List<Point3d> BezierEstimate(List<Point3d> OriginPoint, bool IsClosed, double d, double uu)
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
                v3 *= v1.Length / (v1.Length + v2.Length);
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

}