using Rhino.Geometry;

using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;

namespace MeshClassLibrary
{
    class HeightMapTools
    {
        public HeightMapTools() { }
        public Grasshopper.Kernel.Types.GH_Material GetGHMaterial(string Path)
        {
            Rhino.Display.DisplayMaterial mat = new Rhino.Display.DisplayMaterial();
            mat.SetBitmapTexture(Path, true);
            Grasshopper.Kernel.Types.GH_Material mater = new Grasshopper.Kernel.Types.GH_Material();
            mater.Value = mat;
            return mater;
        }
        public Rhino.Display.DisplayMaterial GetRHMaterial(string Path)
        {
            Rhino.Display.DisplayMaterial mat = new Rhino.Display.DisplayMaterial();
            mat.SetBitmapTexture(Path, true);
            return mat;
        }
        public Mesh LoadHeightMap(double uscale, double vscale, int u, int v, double heightscale, double texturescale, string HeightMapPath)
        {
            if (uscale < 0.001) uscale = 0.001;
            if (vscale < 0.001) vscale = 0.001;
            if (heightscale < 0.001) heightscale = 0.001;
            if (uscale > 10000) uscale = 10000;
            if (vscale > 10000) vscale = 10000;
            if (heightscale > 10000) heightscale = 10000;
            if (u < 4) u = 4;
            if (v < 4) v = 4;
            if (u > 4096) u = 4096;
            if (v > 4096) v = 4096;
            if (texturescale > 4096) texturescale = 4096;
            if (texturescale < 0.001) texturescale = 0.001;
            Mesh mesh = new Mesh();
            try
            {
                Bitmap te2 = new Bitmap(HeightMapPath);
                if (te2.Width < u || te2.Height < v)
                {
                    u = te2.Width; v = te2.Height;
                }
                for (int j = 0; j < v; j++)
                {
                    for (int i = 0; i < u; i++)
                    {
                        mesh.Vertices.Add(new Point3d(i * uscale, te2.GetPixel(i, j).GetBrightness() * heightscale, j * vscale));
                        mesh.TextureCoordinates.Add((double)i / (double)u * texturescale, (double)j / (double)v * texturescale);
                        if (i > 0 && j > 0)
                        {
                            mesh.Faces.AddFace(new MeshFace((j - 1) * u + i - 1, j * u + i - 1, j * u + i, (j - 1) * u + i));
                        }
                    }
                }
                mesh.Compact();
                mesh.UnifyNormals();
                return mesh;
            }
            catch { return null; }
        }
    }
    class BesidgeTools
    {
        public BesidgeTools() { }
        public enum TextureFileType
        {
            JPG,
            PNG,
        }
        private string WorldMachine_RhinoToBesiege(List<Mesh> x, Point3d y)
        {
            string output = "";
            {
                output += "MSnow,size,0" + "\n";
                output += "Wind,size,0" + "\n";
                output += "Cloud,location," + ((int)y.X).ToString() + "," + ((int)y.Z).ToString() + "," + ((int)y.Y).ToString() + "\n";
                output += "Cloud,floorScale,2500,150,2500" + "\n";
                output += "Cloud,size,80" + "\n";
                output += "Camera,farClipPlane,5000" + "\n";
                output += "MWater,size,0" + "\n";
                output += "Triggers,size,0" + "\n";
                output += "Meshes,NoShadow," + x.Count.ToString() + "\n";
                string folderpath = "C:\\Users\\Administrator\\Desktop\\MEP";
                for (int i = 0; i < x.Count; i++)
                {
                    string Name = "L01_Land" + i.ToString();
                    output += "Mesh," + i.ToString() + ",wmesh," + Name + "\n";
                    output += "Mesh," + i.ToString() + ",wmeshcollider," + Name + "\n";
                    output += "Mesh," + i.ToString() + ",materialPropcopy,0" + "\n";
                    output += "Mesh," + i.ToString() + ",location,0,0,0" + "\n";
                    output += "Mesh," + i.ToString() + ",scale,1,1,1" + "\n";
                    string path = folderpath + "\\" + Name + ".obj";
                    //    Write(x[i], path);
                }
            }
            return output;
        }
        private static int SortByNameAscending(string name1, string name2)
        {
            name1 = Path.GetFileNameWithoutExtension(name1);
            name2 = Path.GetFileNameWithoutExtension(name2);
            return name1.CompareTo(name2);
        }
        public string WorldMachineToBesiege(string Folder, Transform xform, string ProgramName)
        {
            List<string> objfile = new List<string>();
            List<string> pngFile = new List<string>();

            _SearchFiles(Folder, "mesh", ref objfile);
            _SearchFiles(Folder, "image", ref pngFile);
            objfile.Sort(SortByNameAscending);
            pngFile.Sort(SortByNameAscending);
            string output = "";
            {
                output += "Snow,size,0" + "\n";
                output += "Wind,size,0" + "\n";
                output += "Cloud,size,0" + "\n";
                output += "Camera,farClipPlane,3000" + "\n";
                output += "Water,size,0" + "\n";
                output += "Triggers,size,0" + "\n";
                output += "Meshes,size," + objfile.Count.ToString() + "\n";
                for (int i = 0; i < objfile.Count; i++)
                {
                    Mesh mesh = ObjToMesh(objfile[i]);
                    mesh.Transform(xform);
                    string name = Path.GetFileNameWithoutExtension(objfile[i]);
                    name = name.Replace("mesh", ProgramName);
                    string str = Folder + "/"
                      + name + ".obj";
                    if (mesh.Normals.Count == 0) mesh.Normals.ComputeNormals();
                    Write(mesh, str);
                    output += eWriteScene(i, str, pngFile[i]);
                }
            }
            return output;
        }
        private void _SearchFiles(string dir, string type, ref List<string> output)
        {
            //通过文件名或者或者扩展名
            if (Directory.Exists(dir))
            {
                foreach (string d in Directory.GetFileSystemEntries(dir))
                {
                    if (File.Exists(d))
                    {
                        string name = Path.GetFileNameWithoutExtension(d);
                        string extension = Path.GetExtension(d);

                        if (name.Contains(type) || extension.Contains(type)) { output.Add(d); }
                    }
                    else
                        _SearchFiles(d, type, ref output);
                }
            }
        }
        public string eSceneToBesiege(List<Mesh> z, string Folder, string ProgramName, TextureFileType filetype)
        {
            string output = "";
            output += "Snow,size,0" + "\n";
            output += "Wind,size,0" + "\n";
            output += "Cloud,size,0" + "\n";
            output += "Camera,farClipPlane,3000" + "\n";
            output += "Water,size,0" + "\n";
            output += "Triggers,size,0" + "\n";
            output += "Meshes,size," + z.Count.ToString() + "\n";
            int co = 0;
            Writebrackets(z, Folder, ProgramName);
            for (int i = 0; i < z.Count; i++)
            {
                if (filetype == TextureFileType.PNG)
                {
                    output += eWriteScene(co, Folder + "/" + ProgramName + " (" + (i + 1).ToString() + ").obj",
                      Folder + "/" + ProgramName + " (" + (i + 1).ToString() + ").png"); co++;
                }
                else if (filetype == TextureFileType.JPG)
                {
                    output += eWriteScene(co, Folder + "/" + ProgramName + " (" + (i + 1).ToString() + ").obj",
                      Folder + "/" + ProgramName + " (" + (i + 1).ToString() + ").jpg"); co++;
                }
            }
            return output;
        }
        public string wSceneToBesiege(List<Mesh> z, string Folder, string ProgramName)
        {
            string output = "";
            output += "Snow,size,0" + "\n";
            output += "Wind,size,0" + "\n";
            output += "Cloud,size,0" + "\n";
            output += "Camera,farClipPlane,3000" + "\n";
            output += "Water,size,0" + "\n";
            output += "Triggers,size,0" + "\n";
            output += "Meshes,size," + z.Count.ToString() + "\n";
            int co = 0;
            Writebrackets(z, Folder, ProgramName);
            for (int i = 0; i < z.Count; i++)
            {
                output += wWriteScene(co, ProgramName, ProgramName);
                co++;
            }
            return output;
        }
        private bool Writebrackets(List<Mesh> mesh, string folderpath, string Name)
        {
            bool hr = false;
            if (mesh.Count < 1) return hr;
            for (int i = 0; i < mesh.Count; i++)
            {
                string path = folderpath + "/" + Name + " (" + (i + 1).ToString() + ").obj";
                hr = Write(mesh[i], path);
            }
            return hr;
        }
        private string eWriteScene(int index, string Name, string Texture)
        {
            string str = "";
            str += "Mesh," + index.ToString() + ",emesh," + Name + "\n";
            str += "Mesh," + index.ToString() + ",emeshcollider," + Name + "\n";
            str += "Mesh," + index.ToString() + ",shader,Standard" + "\n";
            str += "Mesh," + index.ToString() + ",setfloat,_Glossiness,1" + "\n";
            str += "Mesh," + index.ToString() + ",etexture," + Texture + "\n";
            str += "Mesh," + index.ToString() + ",location,0,0,0" + "\n";
            str += "Mesh," + index.ToString() + ",scale,1,1,1" + "\n";
            return str;
        }
        private string wWriteScene(int index, string Name, string Texture)
        {
            string str = "";
            str += "Mesh," + index.ToString() + ",wmesh," + Name + "\n";
            str += "Mesh," + index.ToString() + ",wmeshcollider," + Name + "\n";
            str += "Mesh," + index.ToString() + ",shader,Standard" + "\n";
            str += "Mesh," + index.ToString() + ",setfloat,_Glossiness,1" + "\n";
            str += "Mesh," + index.ToString() + ",stexture," + Texture + "\n";
            str += "Mesh," + index.ToString() + ",location,0,0,0" + "\n";
            str += "Mesh," + index.ToString() + ",scale,1,1,1" + "\n";
            return str;
        }
        private string sWriteScene(int index, string Name, string Name2)
        {
            string str = "";
            str += "Mesh," + index.ToString() + ",wmesh," + Name + "\n";
            str += "Mesh," + index.ToString() + ",wmeshcollider," + Name + "\n";
            str += "Mesh," + index.ToString() + ",smaterialcopy," + Name2 + "\n";
            str += "Mesh," + index.ToString() + ",location,0,0,0" + "\n";
            str += "Mesh," + index.ToString() + ",scale,1,1,1" + "\n";
            return str;
        }
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
                    sb.Append(string.Format("f {0}/{0}/{0} {2}/{2}/{2} {1}/{1}/{1}\n", fs[i].A + 1, fs[i].B + 1, fs[i].C + 1));
                }
                if (fs[i].IsQuad)
                {
                    sb.Append(string.Format("f {3}/{3}/{3} {2}/{2}/{2} {1}/{1}/{1} {0}/{0}/{0}\n", fs[i].A + 1, fs[i].B + 1, fs[i].C + 1, fs[i].D + 1));
                }
            }
            return sb.ToString();
        }
        public Mesh ObjToMesh(string Objpath)
        {
            Mesh mesh = new Mesh();
            StreamReader srd;
            try
            {
                srd = File.OpenText(Objpath);
            }
            catch
            {
                return null;
            }
            try
            {
                while (srd.Peek() != -1)
                {
                    string str = srd.ReadLine();
                    string[] chara = str.Split(new string[] { " " }, StringSplitOptions.RemoveEmptyEntries);
                    if (chara.Length > 2)
                    {
                        if (chara[0] == "v")
                        {
                            mesh.Vertices.Add(
                              Convert.ToSingle(chara[1]),
                              -Convert.ToSingle(chara[3]),
                              Convert.ToSingle(chara[2]));

                        }
                        else if (chara[0] == "vt")
                        {
                            mesh.TextureCoordinates.Add(
                              Convert.ToSingle(chara[1]),
                              Convert.ToSingle(chara[2]));
                        }
                        else if (chara[0] == "vn")
                        {
                            mesh.Normals.Add(
                              Convert.ToSingle(chara[1]),
                              -Convert.ToSingle(chara[3]),
                              Convert.ToSingle(chara[2]));

                        }
                        else if (chara[0] == "f")
                        {
                            if (chara.Length == 4)
                            {
                                int a = Convert.ToInt32(chara[1].Split('/')[0]) - 1;
                                int c = Convert.ToInt32(chara[2].Split('/')[0]) - 1;
                                int b = Convert.ToInt32(chara[3].Split('/')[0]) - 1;
                                mesh.Faces.AddFace(a, b, c);
                            }
                            if (chara.Length == 5)
                            {
                                int a = Convert.ToInt32(chara[1].Split('/')[0]) - 1;
                                int b = Convert.ToInt32(chara[2].Split('/')[0]) - 1;
                                int c = Convert.ToInt32(chara[3].Split('/')[0]) - 1;
                                int d = Convert.ToInt32(chara[4].Split('/')[0]) - 1;
                                mesh.Faces.AddFace(a, b, c, d);
                            }
                        }
                    }
                }
                srd.Close();
                mesh.Compact();
                mesh.UnifyNormals();
            }
            catch
            {
                return null;
            }
            return mesh;
        }
    }
}
