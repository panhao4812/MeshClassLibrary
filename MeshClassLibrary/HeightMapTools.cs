using Rhino.Geometry;
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
            catch{return null;}
        }
    }
}
