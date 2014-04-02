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
    public class BasicFace
    {

        public static List<Line> MeshMaze(Mesh x)
        {
            List<bool> sign;
            List<BasicFace> fs;
            Random rnd = new Random();
            fs = new List<BasicFace>();
            sign = new List<bool>();
            for (int i = 0; i < x.Faces.Count; i++)
            {
                fs.Add(new BasicFace(i));
            }
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = x.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = x.TopologyVertices;
            for (int i = 0; i < el.Count; i++)
            {
                sign.Add(true);
                int[] parel = el.GetConnectedFaces(i);
                if (parel.Length == 2)
                {
                    fs[parel[0]].EdgeIndex.Add(i); fs[parel[0]].FaceIndex.Add(parel[1]);
                    fs[parel[1]].EdgeIndex.Add(i); fs[parel[1]].FaceIndex.Add(parel[0]);
                }
            }
            for (int i = 0; i < x.Faces.Count; i++)
            {
                fs[i].WaveList(rnd);
            }
            int step = 0;
            for (int i = 0; i < fs.Count * 2; i++)
            {
                step = fs[step].FindNext(ref fs, ref sign);
                if (step == -1) break;
                //Print(step.ToString());
            }
            List<Line> output = new List<Line>();
            for (int i = 0; i < el.Count; i++)
            {
                if (sign[i]) output.Add(el.EdgeLine(i));
            }
            return output;
        }


        public List<int> EdgeIndex = new List<int>();
        public List<int> FaceIndex = new List<int>();
        public int ID = -1;
        public int parent = -1;
        public double energy = 0;
        public BasicFace(int i)
        {
            this.ID = i;
        }
        public bool WaveList(Random Rnd)
        {
            int dt = 0;
            if (EdgeIndex.Count != FaceIndex.Count) return false;
            if (this.EdgeIndex.Count <= 1) return false;
            if (this.FaceIndex.Count <= 1) return false;
            for (int i = 0; i <= this.EdgeIndex.Count; i++)
            {
                dt = Rnd.Next(FaceIndex.Count);
                int rep = EdgeIndex[0]; int rep2 = EdgeIndex[dt];
                EdgeIndex[0] = rep2;
                EdgeIndex[dt] = rep;
                rep = FaceIndex[0]; rep2 = FaceIndex[dt];
                FaceIndex[0] = rep2;
                FaceIndex[dt] = rep;
            }
            return true;
        }
        public int FindNext(ref List<BasicFace> fs, ref  List<bool> sign)
        {
            this.energy = 1;
            for (int i = 0; i < this.FaceIndex.Count; i++)
            {
                if (fs[this.FaceIndex[i]].energy == 0)
                {
                    fs[this.FaceIndex[i]].parent = this.ID;
                    sign[this.EdgeIndex[i]] = false;
                    return this.FaceIndex[i];
                }
            }
            return this.parent;
        }
    }
}