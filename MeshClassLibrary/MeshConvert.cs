using Rhino;
using Rhino.Geometry;

using System;
using System.Drawing;
using System.Collections.Generic;

namespace MeshClassLibrary
{
    public class MeshConvert
    {
        MeshCreation mc;
        public MeshConvert() { mc = new MeshCreation(); }
        /// 2Nurbs 
        /// 
        public string Mesh2Topo(ref Mesh mesh, out Mesh mesh2, int edge, out IndexPair data)
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
                        energy[index[k]]++;
                    }
                }
                for (int j = 0; j < firstLoop1.Count; j++)
                {
                    int[] index = vs.ConnectedTopologyVertices(firstLoop1[j]);
                    for (int k = 0; k < index.Length; k++)
                    {
                        if (energy[index[k]] == 1)
                        {
                            firstLoop1[j] = index[k]; sign = false; break;
                        }
                    }
                }
                if (sign) { str += " Loop false,Not quad topo Or To the end"; }
                else { indexPt.AddRange(firstLoop1); }
            }
            ///
            ///*******
            double MX = (double)Column - 1;
            double MY = (double)firstLoop1.Count - 1;
            List<Point3d> output1 = new List<Point3d>();
            List<Point3d> output2 = new List<Point3d>();
            int iCount = 0;
            Color[] cl = new Color[mesh.Vertices.Count];
            //edge
            List<int> pl1 = new List<int>();
            List<int> pl2 = new List<int>();
            List<int> pl3 = new List<int>();
            List<int> pl4 = new List<int>();
            int edge1 = 0, edge2 = 0, edge3 = 0, edge4 = 0;
            //////////
            for (int i = 0; i < Column; i++)
            {
                for (int j = 0; j < firstLoop1.Count; j++)
                {
                    if (i == 0) pl1.Add(iCount);
                    if (i == column - 1) pl3.Add(iCount);
                    if (j == 0) pl2.Add(iCount);
                    if (j == firstLoop1.Count - 1) pl4.Add(iCount);
                    if (i == 0 && j == 0) edge1 = iCount;
                    if (i == 0 && j == firstLoop1.Count - 1) edge2 = iCount;
                    if (i == column - 1 && j == firstLoop1.Count - 1) edge3 = iCount;
                    if (i == column - 1 && j == 0) edge4 = iCount;
                    output1.Add(vs[indexPt[iCount]]);
                    int indexV = vs.MeshVertexIndices(indexPt[iCount])[0];
                    cl[iCount] = mesh.VertexColors[indexV];
                    iCount++;
                    Point3d pt = new Point3d(edge * i, edge * j, 0);
                    output2.Add(pt);
                }
            }
            mesh2 = mc.MeshFromPoints(output2, firstLoop1.Count, Column);
            mesh = mc.MeshFromPoints(output1, firstLoop1.Count, Column);
            mesh.VertexColors.Clear();
            mesh.VertexColors.AppendColors(cl);
            mesh2.VertexColors.Clear();
            mesh2.VertexColors.AppendColors(cl);
            ///edge
            List<Point3d> pts;


            /////////////////////////////////////////////////////////
            Color[] cl1 = new Color[pl1.Count * 2]; pts = new List<Point3d>();
            for (int i = 0; i < pl1.Count; i++)
            {
                pts.Add(new Point3d(mesh2.Vertices[pl1[i]]));
                pts.Add((Point3d)mesh2.Vertices[pl1[i]] + new Vector3d(-edge, 0, 0));
                cl1[i * 2] = mesh2.VertexColors[pl1[i]]; cl1[i * 2 + 1] = mesh2.VertexColors[pl1[i]];
            }
            Mesh mesh3 = mc.MeshFromPoints(pts, 2, pl1.Count);
            mesh3.VertexColors.AppendColors(cl1);
            mesh2.Append(mesh3);
            /////////////////////////////////////////////////////////
            Color[] cl2 = new Color[pl2.Count * 2]; pts = new List<Point3d>();
            for (int i = 0; i < pl2.Count; i++)
            {
                pts.Add(new Point3d(mesh2.Vertices[pl2[i]]));
                pts.Add((Point3d)mesh2.Vertices[pl2[i]] + new Vector3d(0, -edge, 0));
                cl2[i * 2] = mesh2.VertexColors[pl2[i]]; cl2[i * 2 + 1] = mesh2.VertexColors[pl2[i]];
            }
            mesh3 = mc.MeshFromPoints(pts, 2, pl2.Count);
            mesh3.VertexColors.AppendColors(cl2);
            mesh2.Append(mesh3);
            /////////////////////////////////////////////////////////
            Color[] cl3 = new Color[pl3.Count * 2]; pts = new List<Point3d>();
            for (int i = 0; i < pl3.Count; i++)
            {
                pts.Add(new Point3d(mesh2.Vertices[pl3[i]]));
                pts.Add((Point3d)mesh2.Vertices[pl3[i]] + new Vector3d(edge, 0, 0));
                cl3[i * 2] = mesh2.VertexColors[pl3[i]]; cl3[i * 2 + 1] = mesh2.VertexColors[pl3[i]];
            }
            mesh3 = mc.MeshFromPoints(pts, 2, pl3.Count);
            mesh3.VertexColors.AppendColors(cl3);
            mesh2.Append(mesh3);
            /////////////////////////////////////////////////////////
            Color[] cl4 = new Color[pl4.Count * 2]; pts = new List<Point3d>();
            for (int i = 0; i < pl4.Count; i++)
            {
                pts.Add(new Point3d(mesh2.Vertices[pl4[i]]));
                pts.Add((Point3d)mesh2.Vertices[pl4[i]] + new Vector3d(0, edge, 0));
                cl4[i * 2] = mesh2.VertexColors[pl4[i]]; cl4[i * 2 + 1] = mesh2.VertexColors[pl4[i]];
            }
            mesh3 = mc.MeshFromPoints(pts, 2, pl4.Count);
            mesh3.VertexColors.AppendColors(cl4);
            mesh2.Append(mesh3);
            /////////////////////////////////////////////////////////

            Point3d ept = (Point3d)mesh2.Vertices[edge1];
            mesh3 = mc.MeshFromPoints(ept, ept + new Vector3d(-edge, 0, 0), ept + new Vector3d(-edge, -edge, 0), ept + new Vector3d(0, -edge, 0));
            mesh3.VertexColors.Add(mesh2.VertexColors[edge1]);
            mesh3.VertexColors.Add(mesh2.VertexColors[edge1]);
            mesh3.VertexColors.Add(mesh2.VertexColors[edge1]);
            mesh3.VertexColors.Add(mesh2.VertexColors[edge1]);
            mesh2.Append(mesh3);
            ept = (Point3d)mesh2.Vertices[edge2];
            mesh3 = mc.MeshFromPoints(ept, ept + new Vector3d(-edge, 0, 0), ept + new Vector3d(-edge, edge, 0), ept + new Vector3d(0, edge, 0));
            mesh3.VertexColors.Add(mesh2.VertexColors[edge2]);
            mesh3.VertexColors.Add(mesh2.VertexColors[edge2]);
            mesh3.VertexColors.Add(mesh2.VertexColors[edge2]);
            mesh3.VertexColors.Add(mesh2.VertexColors[edge2]);
            mesh2.Append(mesh3);
            ept = (Point3d)mesh2.Vertices[edge3];
            mesh3 = mc.MeshFromPoints(ept, ept + new Vector3d(edge, 0, 0), ept + new Vector3d(edge, edge, 0), ept + new Vector3d(0, edge, 0));
            mesh3.VertexColors.Add(mesh2.VertexColors[edge3]);
            mesh3.VertexColors.Add(mesh2.VertexColors[edge3]);
            mesh3.VertexColors.Add(mesh2.VertexColors[edge3]);
            mesh3.VertexColors.Add(mesh2.VertexColors[edge3]);
            mesh2.Append(mesh3);
            ept = (Point3d)mesh2.Vertices[edge4];
            mesh3 = mc.MeshFromPoints(ept, ept + new Vector3d(edge, 0, 0), ept + new Vector3d(edge, -edge, 0), ept + new Vector3d(0, -edge, 0));
            mesh3.VertexColors.Add(mesh2.VertexColors[edge4]);
            mesh3.VertexColors.Add(mesh2.VertexColors[edge4]);
            mesh3.VertexColors.Add(mesh2.VertexColors[edge4]);
            mesh3.VertexColors.Add(mesh2.VertexColors[edge4]);
            mesh2.Append(mesh3);


            //////////////////////////////////////////////////////
            for (double i = 1; i <= Column; i++)
            {
                for (double j = 1; j <= firstLoop1.Count; j++)
                {
                    mesh.TextureCoordinates.Add(i / (Column + 1), j / (firstLoop1.Count + 1));
                }
            }
            data = new IndexPair((Column + 1) * edge, (firstLoop1.Count + 1) * edge);
            return str;
        }
        public string Mesh2Topo(ref Mesh mesh, out Mesh mesh2, double length, double width, double mX, double mY, double maxX, double maxY)
        {
            //
            //construct a single bitmap but the edge is not fit well.The Vt has not been tested yet.
            //
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
                        energy[index[k]]++;
                    }
                }
                for (int j = 0; j < firstLoop1.Count; j++)
                {
                    int[] index = vs.ConnectedTopologyVertices(firstLoop1[j]);
                    for (int k = 0; k < index.Length; k++)
                    {
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
                    //not tested
                    mesh.TextureCoordinates.Add(new Point2f((Single)((i * length / MX + mX) / maxX),
                      (Single)(j * width / MY + mY) / maxY));
                    mesh2.TextureCoordinates.Add(new Point2f((Single)((i * length / MX + mX) / maxX),
                      (Single)(j * width / MY + mY) / maxY));
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
        public string Mesh2Mesh(Mesh mesh, out  Mesh s)
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

            //   s = NurbsSurface.CreateFromPoints(output1, Column, firstLoop1.Count, 3, 3);

            s = mc.MeshFromPoints(output1, firstLoop1.Count, Column);
            return str;
        }
        public string Mesh2Mesh(Mesh mesh, out  Mesh s, double uscale, double vscale)
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

            //   s = NurbsSurface.CreateFromPoints(output1, Column, firstLoop1.Count, 3, 3);

            s = mc.MeshFromPoints(output1, firstLoop1.Count, Column, uscale, vscale);
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
