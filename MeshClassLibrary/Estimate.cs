using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Generic;
using MathNet.Numerics.LinearAlgebra.Double.Factorization;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace MeshClassLibrary
{
    public class RandomPointEstimate
    {
        public RandomPointEstimate() { }
        public List<Point3d> GetRandomPoint(Mesh mesh, int count)
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
        public List<Point3d> GetRandomPos(List<Point3d> Pts, double radius)
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
        public CylinderSolver() { }
        // public String output = "";
        public static double IterationEndValue = 0.00001;
        public static string Output_Format = "F6";
        public double x, y, z, i, j, k, radius_result;
        public Line GetLine(double length)
        {
            Point3d pt = new Point3d(x, y, z);
            Vector3d v = new Vector3d(i, j, k);
            return new Line(pt, pt + v * length);
        }
        public List<double> init_paras(List<Point3d> pts, Line L)
        {
            List<double> paras = new List<double>();
            double x0, y0, z0, a, b, c, radius;
            Vector3d v1 = L.To - L.From; v1.Unitize();
            x0 = L.FromX; y0 = L.FromY; z0 = L.FromZ;
            a = v1.X; b = v1.Y; c = v1.Z;
            radius = 0;
            foreach (Point3d temp in pts)
            {
                double x, y, z, u, v, w;
                x = temp.X; y = temp.Y; z = temp.Z;
                u = c * (y - y0) - b * (z - z0);
                v = a * (z - z0) - c * (x - x0);
                w = b * (x - x0) - a * (y - y0);
                radius += Math.Sqrt(u * u + v * v + w * w);
            }
            radius /= pts.Count;
            paras.Add(x0);
            paras.Add(y0);
            paras.Add(z0);
            paras.Add(a);
            paras.Add(b);
            paras.Add(c);
            paras.Add(radius);
            return paras;
        }
        public void Estimate(List<Point3d> Rhinodatas, Line L)
        {
            if (Rhinodatas.Count < 6) return;
            List<MathNetPoint> datas = new List<MathNetPoint>();
            for (int ii = 0; ii < Rhinodatas.Count; ii++)
            {
                datas.Add(new MathNetPoint(Rhinodatas[ii]));
                // output += datas[ii].ToString() + "\r\n";
            }
            //output += datas.Count.ToString() + "\r\n";
            double c2, s2, c1, s1;
            double x0, y0, z0, a, b, c, radius;
            double deltavalue_old;
            double deltavalue_new;
            #region try to fit by gauss-newton method
            List<double> current_para = init_paras(Rhinodatas, L);
            deltavalue_old = getDeltaValue(current_para, datas);
            int step = 0;
            while (step < 10000)
            {
                step++;
                x0 = current_para[0];
                y0 = current_para[1];
                z0 = current_para[2];
                a = current_para[3];
                b = current_para[4];
                c = current_para[5];
                radius = current_para[6];

                Vector<double> pos_old = new DenseVector(new double[]{
                        x0,y0,z0
                        });
                Vector<double> vec_old = new DenseVector(new double[]{
                        a,b,c
                        });
                //step I
                #region step 1
                MathNetPoint zero_new = new MathNetPoint(pos_old);
                List<MathNetPoint> transfer_data = new List<MathNetPoint>();
                #endregion
                //step II
                #region step 2
                if (Math.Abs(a - 1) < IterationEndValue)
                {
                    s1 = 0;
                    c1 = 1;
                    s2 = -1;
                    c2 = 0;
                }
                else
                {
                    c1 = c / Math.Sqrt(Math.Pow(b, 2) + Math.Pow(c, 2));
                    s1 = -b / Math.Sqrt(Math.Pow(b, 2) + Math.Pow(c, 2));
                    c2 = (c * c1 - b * s1) / Math.Sqrt(
                        Math.Pow(a, 2) + Math.Pow(c * c1 - b * s1, 2)
                        );
                    s2 = -a / Math.Sqrt(
                        Math.Pow(a, 2) + Math.Pow(c * c1 - b * s1, 2)
                        );
                }
                DenseMatrix U1 = new DenseMatrix(3, 3, new double[]
                       {
                           c2,0,s2,
                           0,1,0,
                           -s2,0,c2,
                        }
                          );
                DenseMatrix U2 = new DenseMatrix(3, 3, new double[]
                       {
                           1,0,0,
                           0,c1,s1,
                           0,-s1,c1,
                       }
                );

                DenseMatrix U = U2 * U1;

                //move all the points to the new positions.
                foreach (MathNetPoint temp in datas)
                {
                    transfer_data.Add((temp - zero_new) * U);
                }
                x0 = 0;
                y0 = 0;
                z0 = 0;
                a = 0;
                b = 0;
                c = 1;
                #endregion
                //step III
                //calculate the jacobian matrix & deviations
                #region step 3
                DenseMatrix J = new DenseMatrix(transfer_data.Count, 5);
                Vector<double> dev = new DenseVector(transfer_data.Count);

                foreach (MathNetPoint temp in transfer_data)
                {
                    double u, v, x, y, z, r; //w;
                    x = temp.x;
                    y = temp.y;
                    z = temp.z;

                    //in theory, but x0 = y0 = z0 =0
                    //          and  a = b = 0 , c=1
                    //u = c * (y - y0) - b * (z - z0);
                    //v = a * (z - z0) - c * (x - x0);
                    //w = b * (x - x0) - a * (y - y0);
                    u = y;
                    v = -x;
                    // w = 0;
                    //double rr = Math.Sqrt(Math.Pow(u, 2) + Math.Pow(v, 2) + Math.Pow(w, 2))
                    r = Math.Sqrt(u * u + v * v);

                    Vector<double> single = new DenseVector(
                        new double[]{
                            -x/r, -y/r,-x*z/r,-y*z/r,-1
                        });
                    J.SetRow(transfer_data.IndexOf(temp), single);
                    dev[transfer_data.IndexOf(temp)] = -(r - radius);
                }
                #endregion
                //step IV
                #region step 4
                Vector<double> P =
               J.Transpose().Multiply(J).Cholesky().Solve(
               J.Transpose().Multiply(dev)
               );
                #endregion
                //step V
                #region step 5
                Vector<double> pos_off =
                    new DenseVector(
                        new double[]{
                            P[0],P[1],-P[0]*P[2]-P[1]*P[3]
                        }
                        );
                Vector<double> pos_new = pos_old + pos_off * U.Transpose();

                Vector<double> vec_off = new DenseVector(
                        new double[]{
                            P[2],P[3],1
                        }
                        );
                Vector<double> vec_new = vec_off * U.Transpose();

                double radius_new = radius + P[4];
                #endregion
                List<double> paras_new = new List<double>();
                paras_new.AddRange(pos_new);
                paras_new.AddRange(vec_new);
                paras_new.Add(radius_new);
                deltavalue_new = getDeltaValue(paras_new, datas);

                double value = Judge_out(paras_new, current_para);
                double value2 = Math.Abs(deltavalue_new - deltavalue_old);
                if (value < IterationEndValue || value2 < IterationEndValue)
                {
                    break;
                }
                else
                {
                    current_para = paras_new;
                    deltavalue_old = deltavalue_new;
                }
            }
            this.x = current_para[0];
            this.y = current_para[1];
            this.z = current_para[2];
            this.i = current_para[3];
            this.j = current_para[4];
            this.k = current_para[5];
            this.radius_result = current_para[6];
            #endregion
        }
        private double getDeltaValue(List<double> para, List<MathNetPoint> datas)
        {
            double x0 = para[0];
            double y0 = para[1];
            double z0 = para[2];
            double a = para[3];
            double b = para[4];
            double c = para[5];
            double radius = para[6];
            double result = 0;
            foreach (MathNetPoint temp in datas)
            {
                double u, v, w;
                double x, y, z;
                double r;
                x = temp.x;
                y = temp.y;
                z = temp.z;
                u = c * (y - y0) - b * (z - z0);
                v = a * (z - z0) - c * (x - x0);
                w = b * (x - x0) - a * (y - y0);
                //u = c * y - b * z;
                //v = a * z - c * x;
                //w = b * x - a * y;
                // r = Math.Sqrt(y * y + x * x);
                r = Math.Sqrt(Math.Pow(u, 2) + Math.Pow(v, 2) + Math.Pow(w, 2))
                    /
                    Math.Sqrt(Math.Pow(a, 2) + Math.Pow(b, 2) + Math.Pow(c, 2));
                result += r - radius;
            }
            return result / datas.Count;
        }
        private double Judge_out(List<double> paras_new, List<double> current_para)
        {
            DenseVector a = DenseVector.OfEnumerable(paras_new);
            DenseVector b = DenseVector.OfEnumerable(current_para);
            return (a - b).Norm(2);
        }
        public override string ToString()
        {
            string temp = "x:  " + x.ToString(Output_Format) + "\r\n" +
                           "y:  " + y.ToString(Output_Format) + "\r\n" +
                           "z:  " + z.ToString(Output_Format) + "\r\n" +
                           "a:  " + i.ToString(Output_Format) + "\r\n" +
                           "b:  " + j.ToString(Output_Format) + "\r\n" +
                           "c:  " + k.ToString(Output_Format) + "\r\n" +
                           "diameter:  " + (radius_result * 2).ToString(Output_Format);
            return temp;
        }
    }
    public class MathNetPoint
    {
        public double x
        {
            get
            {
                return pos[0];
            }
        }
        public double y
        {
            get
            {
                return pos[1];
            }
        }
        public double z
        {
            get
            {
                return pos[2];
            }
        }
        internal Vector<double> _pos;
        public Vector<double> pos
        {
            get
            {
                return _pos.SubVector(0, 3);
            }
        }
        public MathNetPoint(Point3d pt)
        {
            _pos = new DenseVector(new double[] { (double)pt.X, (double)pt.Y, (double)pt.Z, 1 });
        }
        public MathNetPoint(Point3f pt)
        {
            _pos = new DenseVector(new double[] { (double)pt.X, (double)pt.Y, (double)pt.Z, 1 });
        }
        public MathNetPoint(Vector3d pt)
        {
            _pos = new DenseVector(new double[] { (double)pt.X, (double)pt.Y, (double)pt.Z, 1 });
        }
        public MathNetPoint(Vector3f pt)
        {
            _pos = new DenseVector(new double[] { (double)pt.X, (double)pt.Y, (double)pt.Z, 1 });
        }
        public MathNetPoint(double[] input)
        {
            _pos = new DenseVector(new double[] { input[0], input[1], input[2], 1 });
        }
        public MathNetPoint(Vector<double> v)
            : this(v.ToArray())
        {
        }
        public static MathNetPoint operator -(MathNetPoint a, MathNetPoint b)
        {
            return new MathNetPoint(a.pos.Subtract(b.pos));
        }
        public static MathNetPoint operator *(MathNetPoint a, DenseMatrix b)
        {
            return new MathNetPoint(a.pos * b);
        }
        public void Transform(DenseMatrix trans)
        {
            this._pos = (this.pos * trans);
        }
        public static MathNetPoint operator *(DenseMatrix a, MathNetPoint b)
        {
            return new MathNetPoint(a * b.pos);
        }
        public override string ToString()
        {
            return this.x.ToString() + " " + this.y.ToString() + " " + this.z.ToString();
        }
    }
    public class LineSolver
    {
        public LineSolver() { }
        // public String output = "";
        public static double IterationEndValue = 0.00001;
        public static string Output_Format = "F6";
        public double x, y, z, i, j, k, radius_result;
        public void Estimate(List<MathNetPoint> datas)
        {
            double sum_x = 0;
            double sum_y = 0;
            double sum_z = 0;
            foreach (MathNetPoint temp in datas)
            {
                sum_x += temp.x;
                sum_y += temp.y;
                sum_z += temp.z;
            }
            sum_x /= datas.Count;
            sum_y /= datas.Count;
            sum_z /= datas.Count;

            DenseMatrix jacobian = new DenseMatrix(datas.Count, 3);
            foreach (MathNetPoint temp in datas)
            {
                Vector<double> gradient = new DenseVector(3);
                gradient[0] = temp.x - sum_x;
                gradient[1] = temp.y - sum_y;
                gradient[2] = temp.z - sum_z;
                jacobian.SetRow(datas.IndexOf(temp), gradient);
            }
            Svd svd = jacobian.Svd(true);
            // get matrix of left singular vectors with first n columns of U
            Matrix<double> U1 = svd.U().SubMatrix(0, datas.Count, 0, 3);
            // get matrix of singular values
            Matrix<double> S = new DiagonalMatrix(3, 3, svd.S().ToArray());
            // get matrix of right singular vectors
            Matrix<double> V = svd.VT().Transpose();

            Vector<double> parameters = new DenseVector(3);
            parameters = V.Column(0);
            x = sum_x;
            y = sum_y;
            z = sum_z;
            i = parameters[0];
            j = parameters[1];
            k = parameters[2];
        }
        public override string ToString()
        {
            string temp = "x:  " + x.ToString(Output_Format) + "\r\n" +
                "y:  " + y.ToString(Output_Format) + "\r\n" +
                "z:  " + z.ToString(Output_Format) + "\r\n" +
                "a:  " + i.ToString(Output_Format) + "\r\n" +
                "b:  " + j.ToString(Output_Format) + "\r\n" +
                "c:  " + k.ToString(Output_Format);
            return temp;
        }
    }
    public class PlaneSolver
    {
        public PlaneSolver() { }
        // public String output = "";
        public static double IterationEndValue = 0.00001;
        public static string Output_Format = "F6";
        public double x, y, z, i, j, k, radius_result;
        public void Estimate(List<MathNetPoint> datas)
        {
            double sum_x = 0;
            double sum_y = 0;
            double sum_z = 0;
            foreach (MathNetPoint temp in datas)
            {
                sum_x += temp.x;
                sum_y += temp.y;
                sum_z += temp.z;
            }
            sum_x /= datas.Count;
            sum_y /= datas.Count;
            sum_z /= datas.Count;

            DenseMatrix jacobian = new DenseMatrix(datas.Count, 3);
            foreach (MathNetPoint temp in datas)
            {
                Vector<double> gradient = new DenseVector(3);
                gradient[0] = temp.x - sum_x;
                gradient[1] = temp.y - sum_y;
                gradient[2] = temp.z - sum_z;
                jacobian.SetRow(datas.IndexOf(temp), gradient);
            }
            MathNet.Numerics.LinearAlgebra.Double.Factorization.Svd
                svd = jacobian.Svd(true);
            // get matrix of left singular vectors with first n columns of U
            Matrix<double> U1 = svd.U().SubMatrix(0, datas.Count, 0, 3);
            // get matrix of singular values
            Matrix<double> S = new DiagonalMatrix(3, 3, svd.S().ToArray());
            // get matrix of right singular vectors
            Matrix<double> V = svd.VT().Transpose();

            Vector<double> parameters = new DenseVector(3);
            parameters = V.Column(2);
            x = sum_x;
            y = sum_y;
            z = sum_z;
            i = parameters[0];
            j = parameters[1];
            k = parameters[2];
        }
        public override string ToString()
        {
            string temp = "x:  " + x.ToString(Output_Format) + "\r\n" +
                "y:  " + y.ToString(Output_Format) + "\r\n" +
                "z:  " + z.ToString(Output_Format) + "\r\n" +
                "a:  " + i.ToString(Output_Format) + "\r\n" +
                "b:  " + j.ToString(Output_Format) + "\r\n" +
                "c:  " + k.ToString(Output_Format);
            return temp;
        }
    }
}
