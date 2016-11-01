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
                for (int j = 0; j < pls.Count; j++)
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
    }
}
