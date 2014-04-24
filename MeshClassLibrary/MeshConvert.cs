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
        MeshCreation mc;
        public MeshConvert() { mc = new MeshCreation(); }
        /// 2Nurbs 
        /// 
        public string Mesh2Topo(ref Mesh mesh, out Mesh mesh2, double length, double width)
        {
            return Mesh2Topo(ref mesh, out mesh2, length, width, 0, 0, length, width);
        }
        public string Mesh2Topo(ref Mesh mesh, out Mesh mesh2, double length, double width, double mX, double mY, double maxX, double maxY)
        {
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            string str = "";
            List<int> firstLoop1;
            str += FirstEdge(mesh, out firstLoop1);
            double column = (double)vs.Count / (double)firstLoop1.Count;
            int Column = (int)column;
            if (column - Column != 0) str += "Points Count error,Please confirm the topo to be quad style";
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
                if (sign) { str += " Loop false,Not quad topo Or To the end"; }
                else { indexPt.AddRange(firstLoop1); }
            }
            ///*******
            double MX = (double)Column - 1;
            double MY = (double)firstLoop1.Count - 1;
            List<Point3d> output1 = new List<Point3d>();
            List<Point3d> output2 = new List<Point3d>();
            int iCount = 0;

            Color[] cl = new Color[mesh.Vertices.Count];

            for (int i = 0; i < Column; i++)
            {
                for (int j = 0; j < firstLoop1.Count; j++)
                {
                    output1.Add(vs[indexPt[iCount]]);
                    int indexV = vs.MeshVertexIndices(indexPt[iCount])[0];
                    //Print(mesh.VertexColors[0].ToString());
                    // Print(indexV.ToString());
                    cl[iCount] = mesh.VertexColors[indexV];

                    iCount++;
                    Point3d pt = new Point3d(length / MX * i, width / MY * (double)j, 0);
                    output2.Add(pt);
                }
            }
            mesh2 = mc.MeshFromPoints(output2, firstLoop1.Count, Column);
            mesh = mc.MeshFromPoints(output1, firstLoop1.Count, Column);
            mesh2.Transform(Transform.Translation(mX, mY, 0));

            mesh.VertexColors.Clear();
            mesh.VertexColors.AppendColors(cl);

            mesh2.VertexColors.Clear();
            mesh2.VertexColors.AppendColors(cl);

            iCount = 0;
            //Print(mesh.TextureCoordinates.Count.ToString());//unit=100
            for (double i = 0; i < Column; i++)
            {
                for (double j = 0; j < firstLoop1.Count; j++)
                {

                    mesh.TextureCoordinates.Add(new Point2f((Single)((i * 100.0 / Column + mX) / maxX),
                      (Single)(j * 100.0 / firstLoop1.Count + mY) / maxY));
                    mesh2.TextureCoordinates.Add(new Point2f((Single)((i * 100.0 / Column + mX) / maxX),
                      (Single)(j * 100.0 / firstLoop1.Count + mY) / maxY));

                    iCount++;
                }
            }
            return str;
        }
        public string Mesh2Nurbs(Mesh mesh, out  NurbsSurface s)
        {
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            string str = "";
            List<int> firstLoop1;
            str += FirstEdge(mesh, out firstLoop1);
            double column = (double)vs.Count / (double)firstLoop1.Count;
            int Column = (int)column;
            if (column - Column != 0) str += "Points Count error,Please confirm the topo to be quad style";
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
                if (sign) { str += " Loop false,Not quad topo Or To the end"; }
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
