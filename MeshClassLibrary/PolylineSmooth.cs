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
    class PolylineSmooth
    {
        public Polyline cc_Subdivide(Polyline ptlist)
        {
            List<Point3d> ps2 = new List<Point3d>();
            if (ptlist.Count < 3)
            {
                return ptlist;
            }
            if (ptlist[0].DistanceTo(ptlist[ptlist.Count - 1]) > 0.001)
            {
                ps2.Add(ptlist[0]);
                Point3d pt = (ptlist[0] + ptlist[1]) / 2; ps2.Add(pt);
                for (int i = 1; i < ptlist.Count - 1; i++)
                {
                    Point3d p1 = new Point3d(ptlist[i - 1]);
                    Point3d p2 = new Point3d(ptlist[i]);
                    Point3d p3 = new Point3d(ptlist[i + 1]);
                    Point3d p4 = (p2 + p3) / 2;
                    ps2.Add(p1 * (1 / 6.0) + p2 * (2 / 3.0) + p3 * (1 / 6.0));
                    ps2.Add(p4);
                }
                ps2.Add(ptlist[ptlist.Count - 1]);
            }
            else
            {
                if (ptlist.Count < 4)
                {
                    return ptlist;
                }
                Point3d p1 = new Point3d(ptlist[ptlist.Count - 2]);
                Point3d p2 = new Point3d(ptlist[0]);
                Point3d p3 = new Point3d(ptlist[0 + 1]);
                Point3d p4 = (p2 + p3) / 2;
                ps2.Add(p1 * (1 / 6.0) + p2 * (2 / 3.0) + p3 * (1 / 6.0));
                ps2.Add(p4);
                for (int i = 1; i < ptlist.Count - 1; i++)
                {
                    p1 = new Point3d(ptlist[i - 1]);
                    p2 = new Point3d(ptlist[i]);
                    p3 = new Point3d(ptlist[i + 1]);
                    p4 = (p2 + p3) / 2;
                    ps2.Add(p1 * (1 / 6.0) + p2 * (2 / 3.0) + p3 * (1 / 6.0));
                    ps2.Add(p4);
                }
                ps2.Add(ps2[0]);
            }

            Polyline pl2 = new Polyline(ps2);
            return pl2;
        }
        public Polyline cc_Subdivide(Polyline ptlist, int level)
        {
            if (level >= 1)
            {
                for (int i = 0; i < level; i++)
                {
                    ptlist = cc_Subdivide(ptlist);
                }
            }
            return ptlist;
        }
    }
}
