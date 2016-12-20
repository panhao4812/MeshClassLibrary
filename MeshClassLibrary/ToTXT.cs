using Rhino.Geometry;

using System;
using System.Collections.Generic;
using System.IO;

namespace MeshClassLibrary
{
    class ToTXT
    {
        public ToTXT() { }
        public bool Write(List<Polyline> pls, string path)
        {
            if (pls.Count < 1) return false;
            FileStream fs = new FileStream(path, FileMode.OpenOrCreate);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("count: " + pls.Count.ToString());
            for (int i = 0; i < pls.Count; i++)
            {
                Polyline pl = pls[i];
                sw.WriteLine("length: " + pl.Count.ToString());
                for (int j = 0; j < pl.Count; j++)
                {
                    sw.WriteLine(pl[j].X.ToString() + "," + pl[j].Y.ToString() + "," + pl[j].Z.ToString());
                }
            }
            //清空缓冲区、关闭流
            sw.Flush();
            sw.Close();
            return true;
        }
        public bool Write(Polyline pl, string path)
        {
            if (!pl.IsValid) return false;
            FileStream fs = new FileStream(path, FileMode.OpenOrCreate);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("length: " + pl.Count.ToString());
            for (int j = 0; j < pl.Count; j++)
            {
                sw.WriteLine(pl[j].X.ToString() + "," + pl[j].Y.ToString() + "," + pl[j].Z.ToString());
            }
            //清空缓冲区、关闭流
            sw.Flush();
            sw.Close();
            return true;
        }
        public bool Write(List<Point3d> pts, string path)
        {
            if (pts.Count < 1) return false;
            FileStream fs = new FileStream(path, FileMode.OpenOrCreate);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("count: " + pts.Count.ToString());
            for (int i = 0; i < pts.Count; i++)
            {
                Point3d pt = pts[i];
                sw.WriteLine(pt.X.ToString() + "," + pt.Y.ToString() + "," + pt.Z.ToString());
            }
            //清空缓冲区、关闭流
            sw.Flush();
            sw.Close();
            return true;
        }
        public string WriteScene(int index, string Name, string Texture)
        {
            string str = "";
            str += "Mesh," + index.ToString() + ",wmesh," + Name + "\n";
            str += "Mesh," + index.ToString() + ",wmeshcollider," + Name + "\n";
            str += "Mesh," + index.ToString() + ",texture," + Texture + "\n";
            str += "Mesh," + index.ToString() + ",location,0,0,0" + "\n";
            str += "Mesh," + index.ToString() + ",scale,1,1,1" + "\n";
            return str;
        }
        public bool Write(Mesh mesh, string path)
        {
            FileStream fs = new FileStream(path, FileMode.OpenOrCreate);
            StreamWriter sw = new StreamWriter(fs);
            sw.Write(MeshToObj(mesh));
            //清空缓冲区、关闭流
            sw.Flush();
            sw.Close();
            return true;
        }
        public bool Write(List<Mesh> mesh, string folderpath, string Name)
        {
            bool hr = false;
            if (mesh.Count < 1) return hr;
            for (int i = 0; i < mesh.Count; i++)
            {
                string path = folderpath + "\\" + Name + i.ToString() + ".obj";
                hr = Write(mesh[i], path);
            }
            return hr;
        }
        public bool Write(List<Mesh> mesh, string folderpath, List<string> Name)
        {
            bool hr = false;
            if (mesh.Count != Name.Count) return hr;
            for (int i = 0; i < mesh.Count; i++)
            {
                string path = folderpath + "\\" + Name[i] + ".obj";
                hr = Write(mesh[i], path);
            }
            return hr;
        }
        public string MeshToObj(Mesh mesh)
        {
            System.Text.StringBuilder sb = new System.Text.StringBuilder();
            Rhino.Geometry.Collections.MeshVertexList vs = mesh.Vertices;
            sb.Append("#vs:" + vs.Count.ToString() + "\n");
            for (int i = 0; i < vs.Count; i++)
            {
                sb.Append(string.Format("v {0} {1} {2}\n", vs[i].X, vs[i].Z, vs[i].Y));
            }

            Rhino.Geometry.Collections.MeshVertexNormalList ns = mesh.Normals;
            sb.Append("#ns:" + ns.Count.ToString() + "\n");
            for (int i = 0; i < ns.Count; i++)
            {
                sb.Append(string.Format("vn {0} {1} {2}\n", ns[i].X, ns[i].Z, ns[i].Y));
            }

            Rhino.Geometry.Collections.MeshTextureCoordinateList ts = mesh.TextureCoordinates;
            sb.Append("#ts:" + ts.Count.ToString() + "\n");
            for (int i = 0; i < ts.Count; i++)
            {
                sb.Append(string.Format("vt {0} {1}\n", ts[i].X, ts[i].Y));
            }
            Rhino.Geometry.Collections.MeshFaceList fs = mesh.Faces;
            sb.Append("#fs:" + fs.Count.ToString() + "\n");
            for (int i = 0; i < fs.Count; i++)
            {
                if (fs[i].IsTriangle)
                {
                }
                if (fs[i].IsQuad)
                {
                }
            }
            return sb.ToString();
        }
    }
}
