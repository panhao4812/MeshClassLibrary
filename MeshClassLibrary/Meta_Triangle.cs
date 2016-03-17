using Rhino.Geometry;

using System;
using System.Collections.Generic;
namespace MeshClassLibrary
{/*//If you want to have a go at defining your own fields, 
      * you’ll also need to include the gradient function.
      * For this I find Wolfram Alpha comes in very handy – you can just plug in Grad
      * (some function of xyz) and get the result. 
      * For example: grad(sin x * cos y + sin y * cos z + sin z * cos x)
       //*/
    public abstract class ScalarField
    {
        public double Q = 1;
        public void MeshPush(ref Mesh mesh)
        {

            List<Point3d> Pt = new List<Point3d>();
            for (int i = 0; i < mesh.Vertices.Count; i++)
            {
                Pt.Add(new Point3d(mesh.Vertices[i]));
            }
            System.Threading.Tasks.Parallel.For(0, Pt.Count, j =>
            {
                if (Pt[j].IsValid)
                {
                    {
                        double FValue = 1.0;
                        int counter = 0;
                        while (Math.Abs(FValue) > Acc && counter < Iter) //accuracy and max iterations
                        {
                            FValue = this.ValueAt(Pt[j]);
                            Vector3d Gradient = this.GradientAt(Pt[j]);
                            Gradient *= 1.0 / Gradient.SquareLength;
                            Pt[j] = Pt[j] - FValue * Gradient;
                            counter++;
                        }
                    }
                }
            });

            for (int i = 0; i < mesh.Vertices.Count; i++)
            {
                mesh.Vertices[i] = new Point3f((float)Pt[i].X, (float)Pt[i].Y, (float)Pt[i].Z);
            }
        }
        private double Acc = 0.001;
        public int Iter = 15;
        public abstract double ValueAt(Point3d p);
        public abstract Vector3d GradientAt(Point3d p);
    }
    public class ScalarField_Gyroid : ScalarField
    {
        public ScalarField_Gyroid() { }
        public ScalarField_Gyroid(double _Q)
        {
            Q = _Q;
        }
        public override double ValueAt(Point3d P)
        {
            return
              (Math.Sin(P.X) * Math.Cos(P.Y) +
              Math.Sin(P.Y) * Math.Cos(P.Z) +
              Math.Sin(P.Z) * Math.Cos(P.X)) + Q;
        }
        public override Vector3d GradientAt(Point3d P)
        {
            double x = P.X;
            double y = P.Y;
            double z = P.Z;
            return
              new Vector3d(
              Math.Cos(x) * Math.Cos(y) - Math.Sin(x) * Math.Sin(z),
              Math.Cos(y) * Math.Cos(z) - Math.Sin(x) * Math.Sin(y),
              Math.Cos(x) * Math.Cos(z) - Math.Sin(y) * Math.Sin(z));
        }
    }
    public class ScalarField_Blobs : ScalarField
    {
        public ScalarField_Blobs() { }
        public ScalarField_Blobs(double _Q, List<Point3d> P)
        {
            Q = _Q;
            Pts = P;
        }
        public List<Point3d> Pts = new List<Point3d>();
        public override double ValueAt(Point3d P)
        {
            double Sum = 0;
            foreach (Point3d Pt in Pts)
            {
                Sum += (1.0 / (P - Pt).SquareLength);
            }
            return Sum-Q;
        }
        public override Vector3d GradientAt(Point3d P)
        {
            var VSum = new Vector3d();
            foreach (Point3d Pt in Pts)
            {
                Vector3d V = P - Pt;
                double Value = Math.Pow((V.SquareLength), 2);
                VSum -= 2 * V / Value;
            }
            return VSum;
        }
    }

    public class ScalarField_Curve : ScalarField
    {
        public ScalarField_Curve() { }
        public List<Curve> Crv;
        public ScalarField_Curve(double _Q, List<Curve> C)
        {
            Q = _Q;
            Crv = C;
        }
        public override double ValueAt(Point3d P)
        {
            double Sum = 0;
            double t;
            foreach (Curve C in Crv)
            {
                C.ClosestPoint(P, out t);
                Point3d Pt = C.PointAt(t);
                Sum += (1.0 / (P - Pt).Length);
            }
            return Sum-Q;
        }
        public override Vector3d GradientAt(Point3d P)
        {
            double x = P.X;
            double y = P.Y;
            double z = P.Z;
            double t;
            var VSum = new Vector3d();
            foreach (Curve C in Crv)
            {
                C.ClosestPoint(P, out t);
                Point3d Pt = C.PointAt(t);
                Vector3d V = P - Pt;
                double Value = Math.Pow((V.SquareLength), 3 / 2);
                VSum -= new Vector3d(V.X / Value, V.Y / Value, V.Z / Value);
            }
            return VSum;
        }
    }
    public class Meta_Triangle
    {

        public class Matrix_Octahedra
        {
            public List<Mesh> mesh1 = new List<Mesh>();
            public List<Mesh> mesh2 = new List<Mesh>();
            public Matrix_Octahedra(int XCells, int YCells, int ZCells, double CellSize)
            {
                double CS = CellSize;
                double CH = 0.5 * CS;
                double CQ = 0.25 * CS;
                var Grids = new Point3d[2][, ,];
                Grids[0] = new Point3d[XCells, YCells, ZCells];
                Grids[1] = new Point3d[XCells, YCells, ZCells];
                for (int i = 0; i < XCells; i++)
                {
                    for (int j = 0; j < YCells; j++)
                    {
                        for (int k = 0; k < ZCells; k++)
                        {
                            for (int g = 0; g < 2; g++)
                            {
                                Point3d P = new Point3d(i * CS, j * CS, k * CS) + g * new Vector3d(CH, CH, CH);
                                Mesh M = new Mesh();
                                Grids[g][i, j, k] = P;
                                #region mesh construct
                                int VCount;
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CH, P.Y + CQ, P.Z + 00);
                                M.Vertices.Add(P.X + CH, P.Y + 00, P.Z + CQ);
                                M.Vertices.Add(P.X + CH, P.Y - CQ, P.Z + 00);
                                M.Vertices.Add(P.X + CH, P.Y + 00, P.Z - CQ);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 3);
                                M.Faces.AddFace(VCount + 1, VCount + 2, VCount + 3);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y + CH, P.Z + 00);
                                M.Vertices.Add(P.X + 00, P.Y + CH, P.Z + CQ);
                                M.Vertices.Add(P.X - CQ, P.Y + CH, P.Z + 00);
                                M.Vertices.Add(P.X + 00, P.Y + CH, P.Z - CQ);
                                M.Faces.AddFace(VCount + 3, VCount + 1, VCount);
                                M.Faces.AddFace(VCount + 3, VCount + 2, VCount + 1);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y + 00, P.Z + CH);
                                M.Vertices.Add(P.X + 00, P.Y + CQ, P.Z + CH);
                                M.Vertices.Add(P.X - CQ, P.Y + 00, P.Z + CH);
                                M.Vertices.Add(P.X + 00, P.Y - CQ, P.Z + CH);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 3);
                                M.Faces.AddFace(VCount + 1, VCount + 2, VCount + 3);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X - CH, P.Y + CQ, P.Z + 00);
                                M.Vertices.Add(P.X - CH, P.Y + 00, P.Z + CQ);
                                M.Vertices.Add(P.X - CH, P.Y - CQ, P.Z + 00);
                                M.Vertices.Add(P.X - CH, P.Y + 00, P.Z - CQ);
                                M.Faces.AddFace(VCount + 3, VCount + 1, VCount);
                                M.Faces.AddFace(VCount + 3, VCount + 2, VCount + 1);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y - CH, P.Z + 00);
                                M.Vertices.Add(P.X + 00, P.Y - CH, P.Z + CQ);
                                M.Vertices.Add(P.X - CQ, P.Y - CH, P.Z + 00);
                                M.Vertices.Add(P.X + 00, P.Y - CH, P.Z - CQ);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 3);
                                M.Faces.AddFace(VCount + 1, VCount + 2, VCount + 3);

                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y + 00, P.Z - CH);
                                M.Vertices.Add(P.X + 00, P.Y + CQ, P.Z - CH);
                                M.Vertices.Add(P.X - CQ, P.Y + 00, P.Z - CH);
                                M.Vertices.Add(P.X + 00, P.Y - CQ, P.Z - CH);
                                M.Faces.AddFace(VCount + 3, VCount + 1, VCount);
                                M.Faces.AddFace(VCount + 3, VCount + 2, VCount + 1);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y + CQ, P.Z + CQ);
                                M.Vertices.Add(P.X + CH, P.Y + CQ, P.Z + 00);
                                M.Vertices.Add(P.X + CQ, P.Y + CH, P.Z + 00);
                                M.Vertices.Add(P.X + 00, P.Y + CH, P.Z + CQ);
                                M.Vertices.Add(P.X + 00, P.Y + CQ, P.Z + CH);
                                M.Vertices.Add(P.X + CQ, P.Y + 00, P.Z + CH);
                                M.Vertices.Add(P.X + CH, P.Y + 00, P.Z + CQ);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 6);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 1);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X - CQ, P.Y + CQ, P.Z + CQ);
                                M.Vertices.Add(P.X - CQ, P.Y + CH, P.Z + 00);
                                M.Vertices.Add(P.X - CH, P.Y + CQ, P.Z + 00);
                                M.Vertices.Add(P.X - CH, P.Y + 00, P.Z + CQ);
                                M.Vertices.Add(P.X - CQ, P.Y + 00, P.Z + CH);
                                M.Vertices.Add(P.X + 00, P.Y + CQ, P.Z + CH);
                                M.Vertices.Add(P.X + 00, P.Y + CH, P.Z + CQ);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 6);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 1);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X - CQ, P.Y - CQ, P.Z + CQ);
                                M.Vertices.Add(P.X - CH, P.Y - CQ, P.Z + 00);
                                M.Vertices.Add(P.X - CQ, P.Y - CH, P.Z + 00);
                                M.Vertices.Add(P.X + 00, P.Y - CH, P.Z + CQ);
                                M.Vertices.Add(P.X + 00, P.Y - CQ, P.Z + CH);
                                M.Vertices.Add(P.X - CQ, P.Y + 00, P.Z + CH);
                                M.Vertices.Add(P.X - CH, P.Y + 00, P.Z + CQ);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 6);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 1);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y - CQ, P.Z + CQ);
                                M.Vertices.Add(P.X + CQ, P.Y - CH, P.Z + 00);
                                M.Vertices.Add(P.X + CH, P.Y - CQ, P.Z + 00);
                                M.Vertices.Add(P.X + CH, P.Y + 00, P.Z + CQ);
                                M.Vertices.Add(P.X + CQ, P.Y + 00, P.Z + CH);
                                M.Vertices.Add(P.X + 00, P.Y - CQ, P.Z + CH);
                                M.Vertices.Add(P.X + 00, P.Y - CH, P.Z + CQ);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 6);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 1);
                                //******************
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y + CQ, P.Z - CQ);
                                M.Vertices.Add(P.X + CH, P.Y + CQ, P.Z - 00);
                                M.Vertices.Add(P.X + CQ, P.Y + CH, P.Z - 00);
                                M.Vertices.Add(P.X + 00, P.Y + CH, P.Z - CQ);
                                M.Vertices.Add(P.X + 00, P.Y + CQ, P.Z - CH);
                                M.Vertices.Add(P.X + CQ, P.Y + 00, P.Z - CH);
                                M.Vertices.Add(P.X + CH, P.Y + 00, P.Z - CQ);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 1);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 6);

                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X - CQ, P.Y + CQ, P.Z - CQ);
                                M.Vertices.Add(P.X - CQ, P.Y + CH, P.Z - 00);
                                M.Vertices.Add(P.X - CH, P.Y + CQ, P.Z - 00);
                                M.Vertices.Add(P.X - CH, P.Y + 00, P.Z - CQ);
                                M.Vertices.Add(P.X - CQ, P.Y + 00, P.Z - CH);
                                M.Vertices.Add(P.X + 00, P.Y + CQ, P.Z - CH);
                                M.Vertices.Add(P.X + 00, P.Y + CH, P.Z - CQ);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 1);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 6);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X - CQ, P.Y - CQ, P.Z - CQ);
                                M.Vertices.Add(P.X - CH, P.Y - CQ, P.Z - 00);
                                M.Vertices.Add(P.X - CQ, P.Y - CH, P.Z - 00);
                                M.Vertices.Add(P.X + 00, P.Y - CH, P.Z - CQ);
                                M.Vertices.Add(P.X + 00, P.Y - CQ, P.Z - CH);
                                M.Vertices.Add(P.X - CQ, P.Y + 00, P.Z - CH);
                                M.Vertices.Add(P.X - CH, P.Y + 00, P.Z - CQ);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 1);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 6);
                                //----------------------------------------------
                                VCount = M.Vertices.Count;
                                M.Vertices.Add(P.X + CQ, P.Y - CQ, P.Z - CQ);
                                M.Vertices.Add(P.X + CQ, P.Y - CH, P.Z - 00);
                                M.Vertices.Add(P.X + CH, P.Y - CQ, P.Z - 00);
                                M.Vertices.Add(P.X + CH, P.Y + 00, P.Z - CQ);
                                M.Vertices.Add(P.X + CQ, P.Y + 00, P.Z - CH);
                                M.Vertices.Add(P.X + 00, P.Y - CQ, P.Z - CH);
                                M.Vertices.Add(P.X + 00, P.Y - CH, P.Z - CQ);
                                M.Faces.AddFace(VCount, VCount + 2, VCount + 1);
                                M.Faces.AddFace(VCount, VCount + 3, VCount + 2);
                                M.Faces.AddFace(VCount, VCount + 4, VCount + 3);
                                M.Faces.AddFace(VCount, VCount + 5, VCount + 4);
                                M.Faces.AddFace(VCount, VCount + 6, VCount + 5);
                                M.Faces.AddFace(VCount, VCount + 1, VCount + 6);

                                #endregion

                                if (g == 0) mesh1.Add(M);
                                if (g == 1) mesh2.Add(M);
                            }
                        }
                    }
                }
            }
        }
        public class tribox
        {
            #region table
            readonly int[,] triTable = {
      {-1,-1,-1,-1},
      {0,2,1,-1},
      {3,5,0,-1},
      {1,3,5,2},
      {1,4,3,-1},
      {0,2,4,3},
      {4,5,0,1},
      {2,4,5,-1},
      //***********************
      {5,4,2,-1},
      {1,0,5,4},
      {3,4,2,0},
      {3,4,1,-1},
      {2,5,3,1},
      {0,5,3,-1},
      {1,2,0,-1},
      {-1,-1,-1,-1} };
            #endregion
            public edge[] Cedges;//= new edge[6];
            public point[] Points;//= new Point3d[4];
            public tribox(point[] pts, edge[] els)
            {
                Cedges = new edge[6];
                Points = new point[4];
                Cedges[0] = els[0]; Cedges[1] = els[1]; Cedges[2] = els[2];
                Cedges[3] = els[3]; Cedges[4] = els[4]; Cedges[5] = els[5];
                Points[0] = pts[0]; Points[1] = pts[1]; Points[2] = pts[2]; Points[3] = pts[3];
            }
            public tribox(point[] pts)
            {
                Cedges = new edge[6];
                Points = new point[4];
                Points[0] = pts[0]; Points[1] = pts[1]; Points[2] = pts[2]; Points[3] = pts[3];
                Cedges[0] = new edge(pts[0], pts[1]);
                Cedges[1] = new edge(pts[0], pts[2]);
                Cedges[2] = new edge(pts[0], pts[3]);
                Cedges[3] = new edge(pts[1], pts[2]);
                Cedges[4] = new edge(pts[2], pts[3]);
                Cedges[5] = new edge(pts[3], pts[1]);
            }
            public tribox(List<Point3d> pts)
            {
                Cedges = new edge[6];
                Points = new point[4];
                Points[0] = new point(pts[0]); Points[1] = new point(pts[1]);
                Points[2] = new point(pts[2]); Points[3] = new point(pts[3]);
                Cedges[0] = new edge(Points[0], Points[1]);
                Cedges[1] = new edge(Points[0], Points[2]);
                Cedges[2] = new edge(Points[0], Points[3]);
                Cedges[3] = new edge(Points[1], Points[2]);
                Cedges[4] = new edge(Points[2], Points[3]);
                Cedges[5] = new edge(Points[3], Points[1]);
            }
            public tribox(Point3f[] pts)
            {
                Cedges = new edge[6];
                Points = new point[4];
                Points[0] = new point(pts[0]); Points[1] = new point(pts[1]);
                Points[2] = new point(pts[2]); Points[3] = new point(pts[3]);
                Cedges[0] = new edge(Points[0], Points[1]);
                Cedges[1] = new edge(Points[0], Points[2]);
                Cedges[2] = new edge(Points[0], Points[3]);
                Cedges[3] = new edge(Points[1], Points[2]);
                Cedges[4] = new edge(Points[2], Points[3]);
                Cedges[5] = new edge(Points[3], Points[1]);
            }
            public void getValue(double val1, double val2, double val3, double val4)
            {
                this.Points[0].Value = val1;
                this.Points[1].Value = val2;
                this.Points[2].Value = val3;
                this.Points[3].Value = val4;

            }
            public bool Compute(out Mesh mesh, double isolevel)
            {
                for (int i = 0; i < this.Cedges.Length; i++)
                {
                    this.Cedges[i].VertexInterp(isolevel);
                }
                mesh = new Mesh();
                int sign = 0;
                if (Points[0].Value > isolevel) sign |= 1;
                if (Points[1].Value > isolevel) sign |= 2;
                if (Points[2].Value > isolevel) sign |= 4;
                if (Points[3].Value > isolevel) sign |= 8;
                if (sign == 0 || sign == 15) { return false; }
                else
                {
                    int a = triTable[sign, 0];
                    int b = triTable[sign, 1];
                    int c = triTable[sign, 2];
                    int d = triTable[sign, 3];
                    if (d == -1)
                    {
                        mesh.Vertices.Add(this.Cedges[a].Cpoint);
                        mesh.Vertices.Add(this.Cedges[b].Cpoint);
                        mesh.Vertices.Add(this.Cedges[c].Cpoint);
                        mesh.Faces.AddFace(0, 1, 2);
                    }
                    else
                    {
                        mesh.Vertices.Add(this.Cedges[a].Cpoint);
                        mesh.Vertices.Add(this.Cedges[b].Cpoint);
                        mesh.Vertices.Add(this.Cedges[c].Cpoint);
                        mesh.Vertices.Add(this.Cedges[d].Cpoint);
                        mesh.Faces.AddFace(0, 1, 2, 3);
                    }
                    return true;
                }
            }
            public static List<tribox> FromMesh(Mesh mesh)
            {
                List<tribox> output = new List<tribox>();
                if (mesh.Vertices.Count == 0) return output;
                Point3d cen = new Point3d();
                for (int i = 0; i < mesh.Vertices.Count; i++)
                {
                    cen += mesh.Vertices[i];
                }
                cen /= mesh.Vertices.Count;
                for (int i = 0; i < mesh.Faces.Count; i++)
                {
                    Point3f[] pts = new Point3f[4];
                    pts[0] = new Point3f((float)cen.X, (float)cen.Y, (float)cen.Z);
                    pts[1] = mesh.Vertices[mesh.Faces[i].A];
                    pts[2] = mesh.Vertices[mesh.Faces[i].B];
                    pts[3] = mesh.Vertices[mesh.Faces[i].C];
                    output.Add(new tribox(pts));
                }
                return output;
            }
        }
        public class point
        {
            public Point3f loc = new Point3f();
            public point(Point3d pt)
            {
                this.loc = new Point3f((float)pt.X, (float)pt.Y, (float)pt.Z);
            }
            public point(Point3f pt)
            {
                this.loc = pt;
            }
            public Point3d DrawPoint()
            {
                return new Point3d(this.loc);
            }
            public double Value = 0;
            public bool EqualTo(point pt)
            {
                return this.loc.Equals(pt.loc);
            }
        }
        public class edge
        {
            public Point3f Cpoint = new Point3f();
            public point From;
            public point To;
            public edge(int a, int b, List<point> pts)
            {
                From = pts[a]; To = pts[b];
            }
            public edge(point p1, point p2)
            {
                From = p1; To = p2;
            }
            public void VertexInterp(double isolevel)
            {
                if ((From.Value > isolevel && To.Value < isolevel) || (From.Value < isolevel && To.Value > isolevel))
                    VertexInterp(isolevel, From.loc, To.loc, From.Value, To.Value);
            }
            private void VertexInterp(double isolevel, Point3f p1, Point3f p2, double valp1, double valp2)
            {
                double mu;
                if (Math.Abs(isolevel - valp1) < 0.00001) Cpoint = p1;
                if (Math.Abs(isolevel - valp2) < 0.00001) Cpoint = p2;
                if (Math.Abs(valp1 - valp2) < 0.00001) Cpoint = p1;
                Point3f p = new Point3f();
                mu = (isolevel - valp1) / (valp2 - valp1);
                p.X = (float)(p1.X + (p2.X - p1.X) * mu);
                p.Y = (float)(p1.Y + (p2.Y - p1.Y) * mu);
                p.Z = (float)(p1.Z + (p2.Z - p1.Z) * mu);
                Cpoint = p;
            }
            private void VertexInterp(double isolevel, Point3d p1, Point3d p2, double valp1, double valp2)
            {
                double mu;
                if (Math.Abs(isolevel - valp1) < 0.00001)
                    Cpoint = new Point3f((float)p1.X, (float)p1.Y, (float)p1.Z);
                if (Math.Abs(isolevel - valp2) < 0.00001)
                    Cpoint = new Point3f((float)p2.X, (float)p2.Y, (float)p2.Z);
                if (Math.Abs(valp1 - valp2) < 0.00001)
                    Cpoint = new Point3f((float)p1.X, (float)p1.Y, (float)p1.Z);
                Point3f p = new Point3f();
                mu = (isolevel - valp1) / (valp2 - valp1);
                p.X = (float)(p1.X + (p2.X - p1.X) * mu);
                p.Y = (float)(p1.Y + (p2.Y - p1.Y) * mu);
                p.Z = (float)(p1.Z + (p2.Z - p1.Z) * mu);
                Cpoint = p;
            }
            public Line DrawLine(List<Point3f> Points)
            {
                Point3d p1 = new Point3d(From.loc);
                Point3d p2 = new Point3d(To.loc);
                return new Line(p1, p2);
            }
            public Line DrawLine(List<point> Points)
            {
                Point3d p1 = new Point3d(From.loc);
                Point3d p2 = new Point3d(To.loc);
                return new Line(p1, p2);
            }
            public bool EqualTo(edge el)
            {
                if (this.From.Equals(el.From) && this.To.Equals(el.To)) return true;
                if (this.From.Equals(el.To) && this.To.Equals(el.From)) return true;
                return false;
            }
            public bool IsValid()
            {
                return this.From.EqualTo(this.To);
            }
        }
        public Meta_Triangle() { }
        public Mesh Compute(double isolevel, int length, int width, int height, double cellsize)
        {
            Mesh mesh = new Mesh();
            Matrix_Octahedra mo = new Matrix_Octahedra(length, width, height, cellsize);
            List<Mesh> x = new List<Mesh>();
            x.AddRange(mo.mesh1); x.AddRange(mo.mesh2);
            List<tribox> boxes = new List<tribox>();
            for (int i = 0; i < x.Count; i++)
            {
                boxes.AddRange(tribox.FromMesh(x[i]));
            }
            //**************************************************
            Mesh[] meshes = new Mesh[boxes.Count];
            bool[] signs = new bool[boxes.Count];
            //for(int i = 0;i < boxes.Count;i++){
            System.Threading.Tasks.Parallel.For(0, boxes.Count, i =>
            {

                for (int j = 0; j < boxes[i].Points.Length; j++)
                {
                    Point3f P = boxes[i].Points[j].loc;
                    boxes[i].Points[j].Value =
                      (Math.Sin(P.X) * Math.Cos(P.Y) +
                      Math.Sin(P.Y) * Math.Cos(P.Z) +
                      Math.Sin(P.Z) * Math.Cos(P.X)) + 0.8;
                }
                Mesh mesh1;
                signs[i] = boxes[i].Compute(out mesh1, isolevel);
                meshes[i] = mesh1;

            });

            //**************************************************

            for (int i = 0; i < boxes.Count; i++)
            {
                if (signs[i]) mesh.Append(meshes[i]);
            }

            return mesh;
        }
    }

}

