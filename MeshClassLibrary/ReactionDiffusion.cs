using Rhino.Geometry;

using System;
using System.Collections.Generic;
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
                Plane p;         
                if (y.TryGetPlane(out p)) { 
                if (y.Contains(P1,p,Rhino.RhinoDoc.ActiveDoc.ModelAbsoluteTolerance) == Rhino.Geometry.PointContainment.Inside)
                {
                    vs[i].U = 0.5 * (rnd.NextDouble() * 2);
                    vs[i].V = 0.25 * (rnd.NextDouble() * 2);
                }
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
}