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
   public class MeshTree
    {
        public MeshTree() { }
        public MeshCreation mc = new MeshCreation();
        List<IndexPair> id;
        List<Vertice> vs;
        /// <summary>
        /// unstable method, the level is 10
        /// </summary>
   
        public  Mesh ComputeMeshTree(List<Line> x, Point3d y,double firstEnergy,double EnergyDecrease) {
            Vertice.CreateCollection(x, out this.id, out this.vs);
            for (int i = 0; i < vs.Count; i++)
            {
                if (vs[i].equalTo(y)) { vs[i].energe = firstEnergy; break; }
            }
            for (int i = 0; i < 10; i++)
            {
                vs.ForEach(delegate(Vertice v) { v.transferEnerge(EnergyDecrease, ref vs); });
            }

            for (int i = 0; i < vs.Count; i++)
            {
                vs[i].CrateEdges(vs);
                //Print(vs[i].edges.Count.ToString());
            }
            ////////////////

            Mesh mesh = new Mesh();
            for (int i = 0; i < id.Count; i++)
            {
                Polyline pl1 = new Polyline(); Polyline pl2 = new Polyline();
                if (vs[id[i].J].refer.Count == 3)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        if (vs[id[i].J].refer[j] == id[i].I)
                        {
                            pl1 = vs[id[i].J].edges[j]; break;
                        }
                    }
                }
                if (vs[id[i].I].refer.Count == 3)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        if (vs[id[i].I].refer[j] == id[i].J)
                        {
                            pl2 = vs[id[i].I].edges[j]; break;
                        }
                    }
                }
                //Print(pl1.Count.ToString());
                if (pl1.Count == 4 && pl2.Count == 0)
                {
                    Plane p = new Plane(vs[id[i].I].pos, vs[vs[id[i].I].refer[0]].pos - vs[id[i].I].pos);
                    pl2.AddRange(pl1);
                    pl2.Transform(Transform.PlanarProjection(p));

                }
                if (pl1.Count == 0 && pl2.Count == 4)
                {
                    Plane p = new Plane(vs[id[i].J].pos, vs[vs[id[i].J].refer[0]].pos - vs[id[i].J].pos);
                    pl1.AddRange(pl2);
                    pl1.Transform(Transform.PlanarProjection(p));

                }
                if (pl1.Count == 4 && pl2.Count == 4)
                {

                    Plane p1 = new Plane(pl1[0], pl1[1], pl1[2]);
                    Plane p2 = new Plane(pl2[0], pl2[1], pl2[2]);
                    if (Vector3d.VectorAngle(p1.Normal, p2.Normal) > Math.PI / 2) pl2.Reverse();
                    mesh.Append(mc.ClosedBridge(pl1, pl2));

                }
            }

          return mesh;


        }      
    }
}
