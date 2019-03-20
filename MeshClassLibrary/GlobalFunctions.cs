using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Parameters;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.InteropServices;
using System.Text;
//using System.Management;
namespace MeshClassLibrary
{
    public class GH_Global
    {
        public GH_Global()
        {
            initConsole();
        }
        [DllImport("kernel32.dll")]
        public static extern bool AllocConsole();
        [DllImport("kernel32.dll")]
        public static extern bool FreeConsole();
        public IGH_ActiveObject GetCapsure(string Name, string NickName)
        {
            GH_Document ghdoc = Grasshopper.Instances.ActiveCanvas.Document;
            List<IGH_ActiveObject> aos = ghdoc.ActiveObjects();
            for (int i = 0; i < aos.Count; i++)
            {
                IGH_ActiveObject ao = aos[i];
                Print("Type==> " + ao.Name);
                if (ao.Name == Name)
                {
                    Print("Name==> " + ao.NickName);
                    if (ao.NickName == NickName)
                    {
                        return ao;
                    }
                }
            }
            return null;
        }
        public virtual void SetParamPointCapsure(string NickName, IEnumerable<Point3d> pts)
        {
            GH_Document ghdoc = Grasshopper.Instances.ActiveCanvas.Document;
            List<IGH_ActiveObject> aos = ghdoc.ActiveObjects();
            for (int i = 0; i < aos.Count; i++)
            {
                IGH_ActiveObject ao = aos[i];
                Print("Type==> " + ao.Name);
                if (ao.Name == "Point")
                {
                    Param_Point comp = (Param_Point)ao;
                    Print("Name==> " + comp.NickName);
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
            foreach (string str in collection)
            {
                Console.WriteLine(str);
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
        public void Print(Point3d P)
        {
            Console.WriteLine(P.X.ToString() + "," + P.Y.ToString() + "," + P.Z.ToString());
        }
        public void Print(Point3f P)
        {
            Console.WriteLine(P.X.ToString() + "," + P.Y.ToString() + "," + P.Z.ToString());
        }
        public void Print(Vector3d P)
        {
            Console.WriteLine(P.X.ToString() + "," + P.Y.ToString() + "," + P.Z.ToString());
        }
        public void Print(Vector3f P)
        {
            Console.WriteLine(P.X.ToString() + "," + P.Y.ToString() + "," + P.Z.ToString());
        }
        public void Print(Plane P)
        {
            double[] PlaneEquation = P.GetPlaneEquation();
            if (PlaneEquation.Length != 4)
            {
                Console.WriteLine("InValid Plane");
            }
            else
            {
                Console.WriteLine(
                    PlaneEquation[0].ToString() + "," +
                    PlaneEquation[1].ToString() + "," +
                    PlaneEquation[2].ToString() + "," +
                    PlaneEquation[3].ToString());
            }
        }
        public void Print(Line L)
        {
            Point3d P = L.From; Point3d P2 = L.To;
            Console.WriteLine(L.Length.ToString()
                + "/" + P.X.ToString() + "," + P.Y.ToString() + "," + P.Z.ToString()
                + "/" + P2.X.ToString() + "," + P2.Y.ToString() + "," + P2.Z.ToString()
                );
        }
    }
    public class Global
    {
        private bool _enable = false;
        public bool Enable
        {
            get
            {
                return _enable;
            }
            set
            {
                _enable = value;
                initConsole();
            }
        }
        public Global()
        {
            this.Enable = true;
        }
        public static string GetMachineName()
        {
            try
            {
                return System.Environment.MachineName;
            }
            catch 
            {
                return "uMnNk";
            }
        }
        /*
        public static string GetLocalMac()
        {
            string mac = null;
            ManagementObjectSearcher query = new ManagementObjectSearcher("SELECT * FROM Win32_NetworkAdapterConfiguration");
            ManagementObjectCollection queryCollection = query.Get();
            foreach (ManagementObject mo in queryCollection)
            {
                if (mo["IPEnabled"].ToString() == "True")
                    mac = mo["MacAddress"].ToString();
            }
            return (mac);
        }
        */
        public static string GetComputerName()
        {
            return Environment.GetEnvironmentVariable("COMPUTERNAME");
        }
        public Global(bool enable)
        {
            this.Enable = enable;
        }
        public virtual void initConsole()
        {
            if (!_enable) return;
            AllocConsole();
            Console.WindowWidth = 136;
            Console.Clear();
        }
        public void Print(string str)
        {
            if (_enable) Console.WriteLine(str);
        }
        public void Print(double str)
        {
            if (_enable) Console.WriteLine(str.ToString());
        }
        public void Print(int str)
        {
            if (_enable) Console.WriteLine(str.ToString());
        }
        public void Print(float str)
        {
            if (_enable) Console.WriteLine(str.ToString());
        }
        public void Print(IEnumerable<string> collection)
        {
            if (!_enable) return;
            foreach (string str in collection)
            {
                Console.WriteLine(str);
            }

        }
        public void Print(IEnumerable<double> collection)
        {
            if (!_enable) return;
            foreach (double str in collection)
            {
                Console.WriteLine(str.ToString());
            }

        }
        public void Print(IEnumerable<float> collection)
        {
            if (!_enable) return;
            foreach (float str in collection)
            {
                Console.WriteLine(str.ToString());
            }

        }
        public void Print(IEnumerable<int> collection)
        {
            if (!_enable) return;
            foreach (int str in collection)
            {
                Console.WriteLine(str.ToString());
            }
        }
        [DllImport("kernel32.dll")]
        public static extern bool AllocConsole();
        [DllImport("kernel32.dll")]
        public static extern bool FreeConsole();
        public static bool _deletefolder_withoutroot(string dir)
        {
            if (Directory.Exists(dir))
            {
                foreach (string d in Directory.GetFileSystemEntries(dir))
                {
                    if (File.Exists(d))
                        File.Delete(d);
                    else
                        _deletefolder(d);
                }
            }
            else { return false; }
            return true;
        }
        public static bool _deletefolder(string dir)
        {
            if (Directory.Exists(dir))
            {
                foreach (string d in Directory.GetFileSystemEntries(dir))
                {
                    if (File.Exists(d))
                        File.Delete(d);
                    else
                        _deletefolder(d);
                }
                Directory.Delete(dir);
                return true;
            }
            return false;
        }
        public static void _deletefiles_withWhitelist(string dir, ref List<string> _whitelist)
        {
            if (_whitelist.Contains(dir)) return;
            if (Directory.Exists(dir))
            {
                foreach (string d in Directory.GetFileSystemEntries(dir))
                {
                    if (!_whitelist.Contains(d))
                    {
                        if (File.Exists(d)) { File.Delete(d); }
                        else { _deletefiles_withWhitelist(d, ref _whitelist); }
                    }
                }
            }
        }
        public static bool _deletefiles(string dir)
        {
            if (Directory.Exists(dir))
            {
                foreach (string d in Directory.GetFileSystemEntries(dir))
                {
                    if (File.Exists(d))
                        File.Delete(d);
                    else
                        _deletefiles(d);
                }
                return true;
            }
            return false;
        }
        public static void _SearchFiles(string dir, ref List<string> type, ref List<string> output)
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
                        for (int i = 0; i < type.Count; i++)
                        {
                            if (name == type[i]) { output.Add(d); break; }
                            if (extension == type[i]) { output.Add(d); break; }
                        }
                    }
                    else
                        _SearchFiles(d, ref type, ref output);
                }
            }
        }
        public static void _SearchFiles(string dir, string type, ref List<string> output)
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

                        if (name == type || extension == type) { output.Add(d); }
                    }
                    else
                        _SearchFiles(d, type, ref output);
                }
            }
        }
        public static void _SearchFolders(string dir, ref List<string> type, ref List<string> output)
        {//文件夹名称
            if (Directory.Exists(dir))
            {
                string basefolder = Path.GetFileName(dir);
                for (int i = 0; i < type.Count; i++)
                {
                    if (basefolder == type[i]) { output.Add(dir); break; }
                }
                foreach (string d in Directory.GetFileSystemEntries(dir))
                {
                    _SearchFolders(d, ref type, ref output);
                }
            }
        }
        public static void _SearchFolders(string dir, string type, ref List<string> output)
        {//文件夹名称
            if (Directory.Exists(dir))
            {
                string basefolder = Path.GetFileName(dir);
                if (basefolder == type) { output.Add(dir); }
                foreach (string d in Directory.GetFileSystemEntries(dir))
                {
                    _SearchFolders(d, type, ref output);
                }
            }
        }
        public static List<string> From0Find1by2Name(string intput, string folder)
        {
            // For example "D:\0\1\2\..." "D:\0\4\2\..."
            // "D:\0"=intput "D:\0\1\2"=folder
            //{ "D:\0\1","D:\0\4"}=output
            List<string> output = new List<string>();
            foreach (string d1 in Directory.GetFileSystemEntries(intput))
            {
                if (Directory.Exists(d1))
                {
                    bool sign = false;
                    foreach (string d2 in Directory.GetFileSystemEntries(d1))
                    {
                        if (Directory.Exists(d2) && d2.Contains(folder))
                        {
                            sign = true;
                            break;
                        }
                    }

                    if (sign)
                    {
                        output.Add(d1);
                    }
                }
            }
            return output;
        }
        public static bool _OpenFolder(string dir)
        {
            if (Directory.Exists(dir))
            {
                System.Diagnostics.Process.Start("Explorer.exe", dir);
            }
            else { return false; }
            return true;
        }
        public static void ListAllFiles(string dir, ref List<string> output)
        {
            if (Directory.Exists(dir))
            {
                foreach (string d in Directory.GetFileSystemEntries(dir))
                {                
                    if (File.Exists(d))
                    {
                        output.Add(d);
                    }
                    else
                        ListAllFiles(d, ref output);
                }
            }
        }
        public string GetMD5HashFromFile(string fileName)
        {
            try
            {
                FileStream file = new FileStream(fileName, FileMode.Open);
                System.Security.Cryptography.MD5 md5 = new System.Security.Cryptography.MD5CryptoServiceProvider();
                byte[] retVal = md5.ComputeHash(file);
                file.Close();
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < retVal.Length; i++)
                {
                    sb.Append(retVal[i].ToString("x2"));
                }
                return sb.ToString();
            }
            catch (Exception ex)
            {
                Print("GetMD5HashFromFile() fail, error:" +ex.Message);
                return null;                    
            }
        }
    }
}
