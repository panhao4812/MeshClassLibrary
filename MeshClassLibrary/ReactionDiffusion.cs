using Rhino;
using Rhino.Geometry;
using Rhino.DocObjects;
using Rhino.Collections;

using GH_IO;
using GH_IO.Serialization;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

using System;
using System.IO;
using System.Xml;
using System.Xml.Linq;
using System.Linq;
using System.Data;
using System.Drawing;
using System.Reflection;
using System.Collections;
using System.Windows.Forms;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using Rhino.Geometry.Collections;
namespace MeshClassLibrary
{
    class ReactionDiffusion
    {
        List<Vertice3> vs;
        public ReactionDiffusion(Mesh x, Curve y)
        {
            Vertice3.CreateCollection(x, out vs);
            Random rnd = new Random();
            for (int i = 0; i < vs.Count; i++)
            {

                Point3d P1 = new Point3d(vs[i].pos.X, vs[i].pos.Y, 0);
                if (y.Contains(P1) == Rhino.Geometry.PointContainment.Inside)
                {
                    vs[i].U = 0.5 * (rnd.NextDouble() * 2);
                    vs[i].V = 0.25 * (rnd.NextDouble() * 2);
                }
            }

        }
        public List<double> RunReactionDiffusion(Mesh x, Curve y, bool z, double iso)
        {
            for (int i = 0; i < vs.Count; i++)
            {
                vs[i].ComputeLaplation3(vs);
            }
            for (int i = 0; i < vs.Count; i++)
            {
                vs[i].ComputeUV1();
            }
            List<double> U2 = new List<double>();
            for (int i = 0; i < vs.Count; i++)
            {
                U2.Add(1 - vs[i].U);
            }
            return U2;
            /*
          A = mc.MeshTopoVerticeDisplay(x, U2);
          if(z) B = tmf.followlines3(x,
              MeshClassLibrary.MeshTopoVerticeConvert.Data_TopoVertices2Vertice(x, U2), iso);
        }
             * */
        }
    }
    public class Vertice3
    {
        double K = 0.062;
        double F = 0.062;
        public int[] refer;
        public double U = 1, V = 0;
        public double dU = 0, dV = 0;
        //lapU means laplace equation
        //u and v means two different kinds of chemical solution.
        //We always use a energy number to define the density.
        public Point3d pos;
        public Vertice3(Point3d p)
        {
            pos = p;
        }
        public static void CreateCollection(Mesh mesh, out  List<Vertice3> vs)
        {
            Rhino.Geometry.Collections.MeshTopologyVertexList vs1 = mesh.TopologyVertices;
            vs = new List<Vertice3>();
            for (int i = 0; i < vs1.Count; i++)
            {
                vs.Add(new Vertice3(new Point3d(vs1[i].X, vs1[i].Y, vs1[i].Z)));
            }
            for (int i = 0; i < vs1.Count; i++)
            {
                vs[i].refer = vs1.ConnectedTopologyVertices(i);
            }
        }
        public void ComputeLaplation1(List<Vertice3> vs)
        {
            double lapU = 0, lapV = 0;
            for (int i = 0; i < this.refer.Length; i++)
            {
                lapU += vs[this.refer[i]].U;
                lapV += vs[this.refer[i]].V;
            }
            lapU -= U * this.refer.Length;
            lapV -= V * this.refer.Length;
            lapU *= 0.19; lapV *= 0.08;
            dU = lapU - U * V * V + F * (1 - U);
            dV = lapV + U * V * V - (K + F) * V;
        }
        public void ComputeLaplation2(List<Vertice3> vs)
        {
            double lapU = 0, lapV = 0;
            double tot = 0;
            for (int i = 0; i < this.refer.Length; i++)
            {
                double t1 = vs[refer[i]].pos.DistanceTo(this.pos);
                lapU += vs[this.refer[i]].U * t1;
                lapV += vs[this.refer[i]].V * t1;
                tot += t1;
            }
            lapU /= tot; lapU -= U;
            lapV /= tot; lapV -= V;
            lapU *= 0.19 * 2; lapV *= 0.08 * 2;
            dU = lapU - U * V * V + F * (1 - U);
            dV = lapV + U * V * V - (K + F) * V;
        }
        public void ComputeLaplation3(List<Vertice3> vs)
        {
            double lapU = 0, lapV = 0;
            double tot = 0;
            for (int i = 0; i < this.refer.Length; i++)
            {
                double t1 = vs[refer[i]].pos.DistanceTo(this.pos);
                lapU += vs[this.refer[i]].U * 0.1 / t1;
                lapV += vs[this.refer[i]].V * 0.1 / t1;
                tot += 0.1 / t1;
            }
            lapU /= tot; lapU -= U;
            lapV /= tot; lapV -= V;
            lapU *= 0.19 * 2; lapV *= 0.085 * 2;
            dU = lapU - U * V * V + F * (1 - U);
            dV = lapV + U * V * V - (K + F) * V;
        }
        public void ComputeUV1()
        {
            this.U += dU;
            this.V += dV;
        }

    }
}