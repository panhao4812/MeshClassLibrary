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
    public class MeshSimplify
    {
        public MeshSimplify() { }
        private List<M_Point> preDivide(Mesh mesh)
        {
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            List<M_Point> PointList = new List<M_Point>();
            for (int i = 0; i < vs.Count; i++)
            {
                M_Point pt = new M_Point(vs[i]);
                if (vs.MeshVertexIndices(i).Length > 0)
                {
                    pt.N = mesh.Normals[vs.MeshVertexIndices(i)[0]];
                }
                else
                {
                    pt.computeNormal(mesh);
                }
                PointList.Add(pt);
            }
            for (int i = 0; i < vs.Count; i++)
            {
                int[] index = vs.ConnectedTopologyVertices(i);
                for (int j = 0; j < index.Length; j++)
                {
                    PointList[i].refpoints.Add(PointList[index[j]]);
                }
                PointList[i].Sort();
            }

            /////////////////////////////////////////////////////
            for (int i = 0; i < vs.Count; i++)
            {
                PointList[i].order = 0;
            }
            for (int i = 0; i < el.Count; i++)
            {
                if (el.GetConnectedFaces(i).Length == 1)
                {
                    PointList[el.GetTopologyVertices(i).I].order = 5;
                    PointList[el.GetTopologyVertices(i).J].order = 5;
                }
            }
            for (int i = 0; i < vs.Count; i++)
            {
                if (PointList[i].order == 5)
                {
                    if (PointList[i].refpoints.Count != 3)
                    {
                        PointList[i].order = 4;
                    }
                }
                else
                {
                    if (PointList[i].refpoints.Count != 4) PointList[i].order = 4;
                }
            }
            //////////////////////////////////////////////////////////
            for (int k = 0; k < PointList.Count; k++)
            {
                bool sign = true;
                for (int i = 0; i < PointList.Count; i++)
                {
                    if (PointList[i].order == 4)
                    {
                        sign = false;
                        PointList[i].order++;
                        if (PointList[i].refpoints.Count == 4)
                        {
                            if (PointList[i].refpoints[0].order == 5 && PointList[i].refpoints[2].order != 5) { PointList[i].refpoints[2].order = 1; }
                            if (PointList[i].refpoints[1].order == 5 && PointList[i].refpoints[3].order != 5) { PointList[i].refpoints[3].order = 1; }
                            if (PointList[i].refpoints[2].order == 5 && PointList[i].refpoints[0].order != 5) { PointList[i].refpoints[0].order = 1; }
                            if (PointList[i].refpoints[3].order == 5 && PointList[i].refpoints[1].order != 5) { PointList[i].refpoints[1].order = 1; }
                        }
                        else
                        {
                            for (int j = 0; j < PointList[i].refpoints.Count; j++)
                            {
                                PointList[i].refpoints[j].order++;
                            }
                        }
                    }
                }
                for (int i = 0; i < PointList.Count; i++)
                {
                    if (PointList[i].order > 0 && PointList[i].order < 4) PointList[i].order = 4;
                }
                if (sign) { break; }
            }
            return PointList;

        }
        public List<Line> MeshProfile(Mesh mesh)
        {
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            List<M_Point> PointList = preDivide(mesh);
            List<Line> output = new List<Line>();
            for (int i = 0; i < el.Count; i++)
            {
                int a = el.GetTopologyVertices(i).I;
                int b = el.GetTopologyVertices(i).J;
                if (PointList[a].order > 0 && PointList[b].order > 0)
                {
                    output.Add(el.EdgeLine(i));
                }
            }
            /*   List<Point3d> output1 = new List<Point3d>();
               List<string> output2 = new List<string>();
               for(int i = 0;i < PointList.Count;i++){
                 output1.Add(PointList[i].pos);
                 string str = PointList[i].order.ToString();
                 if(PointList[i].order == 0)str = "";
                 output2.Add(str);
               }*/
            return output;
        }
        public List<Mesh> MeshSeperate(Mesh mesh)
        {
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            List<M_Point> PointList = preDivide(mesh);
            List<M_Face> faces = new List<M_Face>();
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                faces.Add(new M_Face(mesh.Faces[i]));
            }
            for (int i = 0; i < el.Count; i++)
            {
                int a = el.GetTopologyVertices(i).I;
                int b = el.GetTopologyVertices(i).J;
                if (PointList[a].order == 0 || PointList[b].order == 0)
                {
                    int[] index = el.GetConnectedFaces(i);
                    if (index.Length == 2)
                    {
                        faces[index[0]].reffaces.Add(faces[index[1]]);
                        faces[index[1]].reffaces.Add(faces[index[0]]);
                    }
                }
            }
            List<Mesh> outputMesh = new List<Mesh>();
            List<List<M_Face>> outtemp = M_Face.Group(faces);
            for (int i = 0; i < outtemp.Count; i++)
            {
                Mesh meshout = new Mesh();
                meshout.Vertices.AddVertices(mesh.Vertices);
                for (int j = 0; j < outtemp[i].Count; j++)
                {
                    meshout.Faces.AddFace(outtemp[i][j].face);
                }
                meshout.Compact();
                outputMesh.Add(meshout);
            }
            return outputMesh;

            /*
    List<Point3d> outpos = new List<Point3d>();
    List<string> outsign = new  List<string>();
    for(int i = 0;i < faces.Count;i++){
      outpos.Add(faces[i].pos(mesh.Vertices));
      string str = i.ToString() + "/" + faces[i].reffaces.Count.ToString();
      outsign.Add(str);
    }
    B = outpos;C = outsign;
    */
        }
        public List<NurbsSurface> MeshSeperate2Nurbs(Mesh mesh)
        {
            MeshConvert conv = new MeshConvert();
            List<NurbsSurface> output = new List<NurbsSurface>();
            List<Mesh> input = MeshSeperate(mesh);

            input.ForEach(delegate(Mesh mesh1)
            {
                NurbsSurface Surf;
                conv.Mesh2Nurbs(mesh1, out Surf);
                output.Add(Surf);
            });
            return output;
        }
    }

}