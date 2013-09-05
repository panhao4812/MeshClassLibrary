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
    public class MeshConvert
    {
        public MeshConvert() { }
        /// 2Nurbs 
        public string Mesh2Nurbs(Mesh mesh, out  NurbsSurface s)
        {
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            string str = "";
            List<int> firstLoop1;
            str += FirstEdge(mesh, out firstLoop1);
            double column = (double)vs.Count / (double)firstLoop1.Count;
            int Column = (int)column;
            if (column - Column != 0) str += "Points Count error,Please confirm the topo is quad";
            int[] energy = new int[vs.Count];
            List<int> indexPt = new List<int>();
            indexPt.AddRange(firstLoop1);
            for (int i = 0; i < firstLoop1.Count; i++) { energy[firstLoop1[i]] = 1; }
            for (int i = 0; i < Column - 1; i++)
            {
                bool sign = true;
                for (int j = 0; j < firstLoop1.Count; j++)
                {
                    int[] index = vs.ConnectedTopologyVertices(firstLoop1[j]);
                    for (int k = 0; k < index.Length; k++)
                    {
                        //Print("j:" + j.ToString() + " k:" + k.ToString() + " energy: " + energy[index[k]].ToString());////////
                        energy[index[k]]++;
                    }
                }
                //Print("///");
                for (int j = 0; j < firstLoop1.Count; j++)
                {

                    int[] index = vs.ConnectedTopologyVertices(firstLoop1[j]);
                    for (int k = 0; k < index.Length; k++)
                    {
                        // Print("j:" + j.ToString() + " k:" + k.ToString() + " energy: " + energy[index[k]].ToString());////////
                        if (energy[index[k]] == 1)
                        {
                            firstLoop1[j] = index[k]; sign = false; break;
                        }
                    }
                }
                if (sign) { str += " Loop false,Not a quad topo Or To the end"; }
                else { indexPt.AddRange(firstLoop1); }
            }
            List<Point3d> output1 = new List<Point3d>();
            for (int i = 0; i < indexPt.Count; i++)
            {
                output1.Add(vs[indexPt[i]]);
            }
            s = NurbsSurface.CreateFromPoints(output1, Column, firstLoop1.Count, 3, 3);
            return str;
        }
        public string FirstEdge(Mesh mesh, out  List<int> firstLoop1)
        {
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            firstLoop1 = new List<int>();
            int firstPoint = -1;
            for (int i = 0; i < vs.Count; i++)
            {
                if (vs.ConnectedTopologyVertices(i).Length == 2)
                {
                    firstPoint = i; break;
                }
            }
            if (firstPoint == -1) { return "can not find a 2 degree Vertex,Please check the mesh boundary"; }
            int SecondPoint = vs.ConnectedTopologyVertices(firstPoint)[0];
            int ThirdPoint = vs.ConnectedTopologyVertices(firstPoint)[1];
            int firstPoint1 = firstPoint;
            firstLoop1.Add(firstPoint);



            if (vs.ConnectedTopologyVertices(ThirdPoint).Length == 2)
            {
                firstLoop1.Add(ThirdPoint);
            }
            else if (vs.ConnectedTopologyVertices(SecondPoint).Length == 2)
            {
                firstLoop1.Add(SecondPoint);
            }
            else
            {
                firstLoop1.Add(SecondPoint);
                for (int i = 0; i < vs.Count; i++)
                {
                    bool stop = false;
                    int[] index = vs.ConnectedTopologyVertices(SecondPoint);
                    for (int j = 0; j < index.Length; j++)
                    {
                        if (index[j] != firstPoint1)
                        {
                            if (vs.ConnectedTopologyVertices(index[j]).Length == 3)
                            {
                                firstPoint1 = SecondPoint;
                                SecondPoint = index[j];
                                firstLoop1.Add(SecondPoint);
                                break;
                            }//if
                            if (vs.ConnectedTopologyVertices(index[j]).Length == 2)
                            {
                                firstPoint1 = SecondPoint;
                                SecondPoint = index[j];
                                firstLoop1.Add(SecondPoint);
                                stop = true;
                                break;
                            }//if
                        }
                    }
                    if (stop) break;
                }
            }
            return "Done";
        }
    }
}
