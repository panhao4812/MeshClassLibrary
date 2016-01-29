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
}