using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Parameters;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace MeshClassLibrary
{

   public  class Global
    {/// 启动控制台 
            public Global()
        {
            initConsole();
        }
        [DllImport("kernel32.dll")]
        public static extern bool AllocConsole();
        /// 释放控制台 
        [DllImport("kernel32.dll")]
        public static extern bool FreeConsole();
        public IGH_ActiveObject GetCapsure(string Name,string NickName)
        {
            GH_Document ghdoc = Grasshopper.Instances.ActiveCanvas.Document;
            List<IGH_ActiveObject> aos = ghdoc.ActiveObjects();
            for (int i = 0; i < aos.Count; i++)
            {
                IGH_ActiveObject ao = aos[i];
                Print("==> " + ao.Name);
                if (ao.Name == "Point")
                {            
                    Print("==> " + ao.NickName);
                    if (ao.NickName == NickName)
                    {
                        return ao;
                    }
                }
            }
            return null;
        }
        public void SetParamPointCapsure(string NickName , IEnumerable<Point3d> pts)
        {
            GH_Document ghdoc = Grasshopper.Instances.ActiveCanvas.Document;
            List<IGH_ActiveObject> aos = ghdoc.ActiveObjects();
            for (int i = 0; i < aos.Count; i++)
            {
                IGH_ActiveObject ao = aos[i];
                Print("==> "+ao.Name);
                if (ao.Name == "Point")
                {
                    Param_Point comp = (Param_Point)ao;
                    Print("==> " + comp.NickName);
                    if (ao.NickName == NickName)
                    {
                        GH_Path path = new GH_Path(1);        
                        comp.AddVolatileDataList(path, pts);
                    }
                }
            }        
        }
        public virtual void initConsole()
        {
            AllocConsole();
            Console.WindowWidth = 136;
            Console.BackgroundColor = ConsoleColor.White;
            Console.ForegroundColor = ConsoleColor.Black;
            Console.Clear();
        }
        public void Print(string str)
        {
            Console.WriteLine(str);
        }
        public void Print(double str)
        {
            Console.WriteLine(str.ToString());
        }
        public void Print(int str)
        {
            Console.WriteLine(str.ToString());
        }
        public void Print(float str)
        {
            Console.WriteLine(str.ToString());
        }
        public void Print(IEnumerable<string> collection)
        {
            foreach(string str in collection){
                Console.WriteLine(str.ToString());
            }
            
        }
        public void Print(IEnumerable<double> collection)
        {
            foreach (double str in collection)
            {
                Console.WriteLine(str.ToString());
            }

        }
        public void Print(IEnumerable<float> collection)
        {
            foreach (float str in collection)
            {
                Console.WriteLine(str.ToString());
            }

        }
        public void Print(IEnumerable<int> collection)
        {
            foreach (int str in collection)
            {
                Console.WriteLine(str.ToString());
            }

        }
    }
}
