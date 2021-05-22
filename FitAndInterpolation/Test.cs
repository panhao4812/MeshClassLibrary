using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Generic;
using MathNet.Numerics.LinearAlgebra.Double.Factorization;
using Rhino.Geometry;
using System.Collections.Generic;
using System;

namespace FitAndInterpolation
{
    class PlaneFitSolver
    {
        public Plane Slerp(Plane starting, Plane ending, double t)
        {

            double[] p1 = starting.GetPlaneEquation();
            double[] p2 = starting.GetPlaneEquation();
            //double temp = 0;
           // temp = p1[0]; p1[0] = p1[3]; p1[3] = temp;
           // temp = p2[0]; p2[0] = p2[3]; p2[3] = temp;
            double[] sq = Slerp(p1, p2, t);
            return new Plane(sq[0], sq[1], sq[2], sq[3]);
        }
        double[] Slerp(double[] starting, double[] ending, double t)
        {

            double[] result = new double[4];
            double cosa = starting[0] * ending[0] + starting[1] * ending[1] + starting[2] * ending[2] + starting[3] * ending[3];

            // If the dot product is negative, the quaternions have opposite handed-ness and slerp won't take
            // the shorter path. Fix by reversing one quaternion.
            if (cosa < 0.0)
            {
                ending[0] = -ending[0];
                ending[1] = -ending[1];
                ending[2] = -ending[2];
                ending[3] = -ending[3];
                cosa = -cosa;
            }

            double k0, k1;

            // If the inputs are too close for comfort, linearly interpolate
            if (cosa > 0.9995)
            {
                k0 = 1.0f - t;
                k1 = t;
            }
            else
            {
                double sina = Math.Sqrt(1.0 - cosa * cosa);
                double a = Math.Atan2(sina, cosa);
                k0 = Math.Sin((1.0 - t) * a) / sina;
                k1 = Math.Sin(t * a) / sina;
            }
            result[0] = starting[0] * k0 + ending[0] * k1;
            result[1] = starting[1] * k0 + ending[1] * k1;
            result[2] = starting[2] * k0 + ending[2] * k1;
            result[3] = starting[3] * k0 + ending[3] * k1;
            return result;
        }
    }
}
