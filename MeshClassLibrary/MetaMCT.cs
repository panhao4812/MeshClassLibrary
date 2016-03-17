using Rhino.Geometry;
using Grasshopper;
using Grasshopper.Kernel.Data;

using System;
using System.Linq;
using System.Collections.Generic;
namespace MeshClassLibrary
{
    public class MetaMCT
    {
        public MetaMCT(List<Point3d> c, List<double> rad, double Res)
        {
            this.use_charges = c;
            this.use_radii = rad;
            this.res = Res;
        }
        public MetaMCT(List<Point3d> c)
        {
            this.use_charges = c;
            use_radii = new List<double>();
            for (int i = 0; i < c.Count; i++)
            {
                use_radii.Add(100);
            }
        }
        private List<Point3d> use_charges;
        private List<double> use_radii;
        private double cA = -0.444444;
        private double cB = 1.888889;
        private double cC = -2.444444;
        public double res = 50;
        public double use_iso = 0.1;
        public Mesh Metaball(double iso)
        {
            this.use_iso = iso;
            Plane plane = new Plane();
            Plane.FitPlaneToPoints(use_charges, out plane);
            plane.Origin = use_charges[0];
            Interval xSize = new Interval(-use_radii[0], use_radii[0]);
            Box box = new Box(plane, xSize, xSize, xSize);
            int num27 = use_charges.Count - 1;
            for (int i = 1; i <= num27; i++)
            {
                plane.Origin = use_charges[i];
                box.Union(plane.PointAt(-use_radii[i], -use_radii[i], -use_radii[i]));
                box.Union(plane.PointAt(use_radii[i], use_radii[i], use_radii[i]));
            }
            box.Inflate(res);
            int xSet = (int)Math.Round((double)(box.X.Length / res));
            int ySet = (int)Math.Round((double)(box.Y.Length / res));
            int zSet = (int)Math.Round((double)(box.Z.Length / res));
            double xLength = box.X.Length / ((double)xSet);
            double yLength = box.Y.Length / ((double)ySet);
            double zLength = box.Z.Length / ((double)zSet);
            double xBase = xLength / 2.0;
            double yBase = yLength / 2.0;
            double zBase = zLength / 2.0;
            plane.Origin = box.GetCorners()[0];
            List<Point3d> list = new List<Point3d>();
            Mesh mesh = new Mesh();

            for (int j = 0; j <= xSet - 1; j++)
            {
                for (int m = 0; m <= ySet - 1; m++)
                {
                    for (int n = 0; n <= zSet - 1; n++)
                    {
                        Point3d item = plane.PointAt(xBase + (xLength * j), yBase + (yLength * m), zBase + (zLength * n));
                        int xSet1 = use_charges.Count - 1;
                        for (int k = 0; k <= xSet1; k++)
                        {
                            if (item.DistanceTo(use_charges[k]) < (use_radii[k] + res))
                            {
                                Plane pl = plane;
                                pl.Origin = item;
                                Mesh other = this.local_tet(pl, xBase, yBase, zBase);
                                mesh.Append(other);
                                list.Add(item);
                                break;
                            }
                        }
                    }
                }
            }
            mesh.Vertices.CombineIdentical(true, true);
            mesh.UnifyNormals();
            return mesh;
        }
        public Mesh[] Smooth(Mesh mesh, int smooth)
        {
            ////////////////////////////
            Mesh[] msh_smooth = mesh.SplitDisjointPieces();

            for (int k = 0; k < msh_smooth.Length; k++)
            {
                msh_smooth[k].Vertices.CombineIdentical(true, true);
                msh_smooth[k].FaceNormals.ComputeFaceNormals();
                msh_smooth[k].Normals.ComputeNormals();
                msh_smooth[k].UnifyNormals();
                for (int j = 0; j < smooth; j++)
                {
                    #region LaplaceFunction
                    Point3f[] new_v = msh_smooth[k].Vertices.ToArray<Point3f>();

                    for (int i = 0; i < msh_smooth[k].TopologyVertices.Count; i++)
                    {
                        int[] mv = msh_smooth[k].TopologyVertices.MeshVertexIndices(i);
                        int[] ctv = msh_smooth[k].TopologyVertices.ConnectedTopologyVertices(i);
                        Single new_x = 0;
                        Single new_y = 0;
                        Single new_z = 0;
                        foreach (int ctv_i in ctv)
                        {
                            int[] vtc = msh_smooth[k].TopologyVertices.MeshVertexIndices(ctv_i);

                            new_x += msh_smooth[k].Vertices[vtc[0]].X / (Single)ctv.Length;
                            new_y += msh_smooth[k].Vertices[vtc[0]].Y / (Single)ctv.Length;
                            new_z += msh_smooth[k].Vertices[vtc[0]].Z / (Single)ctv.Length;
                        }
                        new_v[mv[0]] = new Point3f(new_x, new_y, new_z);
                    }
                    for (int i = 0; i < msh_smooth[k].Vertices.Count; i++)
                    {
                        msh_smooth[k].Vertices[i] = new_v[i];
                    }

                    #endregion
                    msh_smooth[k].FaceNormals.ComputeFaceNormals();
                    msh_smooth[k].Normals.ComputeNormals();

                    for (int i = 0; i < new_v.Length; i++)
                    {
                        Point3d pointd3 = new_v[i];
                        Point3d pointd4 = pointd3 + (Vector3d)msh_smooth[k].Normals[i] * (res / 10);
                        double v1 = this.calculate_field(pointd3);
                        double v2 = this.calculate_field(pointd4);
                        Point3d pointd2 = new Point3d();
                        if (v2 != v1)
                        {
                            pointd2 = this.interp_vertex(pointd3, pointd4, v1, v2);
                        }
                        else
                        {
                            pointd2 = pointd3;
                        }
                        if (pointd3.DistanceTo(pointd2) > this.res)
                        {
                            new_v[i] = new Point3f((float)pointd3.X, (float)pointd3.Y, (float)pointd3.Z);
                        }
                        else
                        {
                            new_v[i] = new Point3f((float)pointd2.X, (float)pointd2.Y, (float)pointd2.Z);
                        }
                    }

                    for (int i = 0; i < msh_smooth[k].Vertices.Count; i++)
                    {
                        msh_smooth[k].Vertices[i] = new_v[i];
                    }
                }
            }
            //*/
            return msh_smooth;
        }
        private Point3d interp_vertex(Point3d p1, Point3d p2, double v1, double v2)
        {
            return new Point3d(p1 + ((Point3d)(((this.use_iso - v1) / (v2 - v1)) * (p2 - p1))));
        }
        private Mesh local_tet(Plane pl, double xBase, double yBase, double zBase)
        {
            List<Point3d> list = new List<Point3d> {
          pl.PointAt(-xBase, yBase, -zBase),
          pl.PointAt(xBase, yBase, -zBase),
          pl.PointAt(xBase, -yBase, -zBase),
          pl.PointAt(-xBase, -yBase, -zBase),
          pl.PointAt(-xBase, yBase, zBase),
          pl.PointAt(xBase, yBase, zBase),
          pl.PointAt(xBase, -yBase, zBase),
          pl.PointAt(-xBase, -yBase, zBase)
          };

            List<double> list2 = new List<double>();
            foreach (Point3d pointd in list)
            {
                double item = this.calculate_field(pointd);
                list2.Add(item);
            }
            DataTree<Point3d> tree = new DataTree<Point3d>();
            Point3d[] data = new Point3d[] { list[0], list[2], list[3], list[7] };
            tree.AddRange(data, new GH_Path(0));
            data = new Point3d[] { list[0], list[2], list[6], list[7] };
            tree.AddRange(data, new GH_Path(1));
            data = new Point3d[] { list[0], list[4], list[6], list[7] };
            tree.AddRange(data, new GH_Path(2));
            data = new Point3d[] { list[0], list[6], list[1], list[2] };
            tree.AddRange(data, new GH_Path(3));
            data = new Point3d[] { list[0], list[6], list[1], list[4] };
            tree.AddRange(data, new GH_Path(4));
            data = new Point3d[] { list[5], list[6], list[1], list[4] };
            tree.AddRange(data, new GH_Path(5));
            DataTree<double> tree2 = new DataTree<double>();
            tree2.AddRange(new double[] { list2[0], list2[2], list2[3], list2[7] }, new GH_Path(0));
            tree2.AddRange(new double[] { list2[0], list2[2], list2[6], list2[7] }, new GH_Path(1));
            tree2.AddRange(new double[] { list2[0], list2[4], list2[6], list2[7] }, new GH_Path(2));
            tree2.AddRange(new double[] { list2[0], list2[6], list2[1], list2[2] }, new GH_Path(3));
            tree2.AddRange(new double[] { list2[0], list2[6], list2[1], list2[4] }, new GH_Path(4));
            tree2.AddRange(new double[] { list2[5], list2[6], list2[1], list2[4] }, new GH_Path(5));
            Mesh mesh = new Mesh();
            foreach (GH_Path path in tree2.Paths)
            {
                List<Point3d> list3 = tree.Branch(path);
                List<double> list4 = tree2.Branch(path);
                int num2 = 0;
                if (list4[0] < this.use_iso)
                {
                    num2 |= 1;
                }
                if (list4[1] < this.use_iso)
                {
                    num2 |= 2;
                }
                if (list4[2] < this.use_iso)
                {
                    num2 |= 4;
                }
                if (list4[3] < this.use_iso)
                {
                    num2 |= 8;
                }
                Mesh other = new Mesh();
                switch (num2)
                {
                    case 1:
                    case 14:
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[1], list4[0], list4[1]));
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[2], list4[0], list4[2]));
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[3], list4[0], list4[3]));
                        other.Faces.AddFace(0, 1, 2);
                        break;

                    case 2:
                    case 13:
                        other.Vertices.Add(this.interp_vertex(list3[1], list3[0], list4[1], list4[0]));
                        other.Vertices.Add(this.interp_vertex(list3[1], list3[3], list4[1], list4[3]));
                        other.Vertices.Add(this.interp_vertex(list3[1], list3[2], list4[1], list4[2]));
                        other.Faces.AddFace(0, 1, 2);
                        break;

                    case 3:
                    case 12:
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[3], list4[0], list4[3]));
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[2], list4[0], list4[2]));
                        other.Vertices.Add(this.interp_vertex(list3[1], list3[3], list4[1], list4[3]));
                        other.Faces.AddFace(0, 1, 2);
                        other.Vertices.Add(this.interp_vertex(list3[1], list3[2], list4[1], list4[2]));
                        other.Faces.AddFace(2, 3, 1);
                        break;

                    case 4:
                    case 11:
                        other.Vertices.Add(this.interp_vertex(list3[2], list3[0], list4[2], list4[0]));
                        other.Vertices.Add(this.interp_vertex(list3[2], list3[1], list4[2], list4[1]));
                        other.Vertices.Add(this.interp_vertex(list3[2], list3[3], list4[2], list4[3]));
                        other.Faces.AddFace(0, 1, 2);
                        break;

                    case 5:
                    case 10:
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[1], list4[0], list4[1]));
                        other.Vertices.Add(this.interp_vertex(list3[2], list3[3], list4[2], list4[3]));
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[3], list4[0], list4[3]));
                        other.Faces.AddFace(0, 1, 2);
                        other.Vertices.Add(this.interp_vertex(list3[1], list3[2], list4[1], list4[2]));
                        other.Faces.AddFace(0, 3, 1);
                        break;

                    case 6:
                    case 9:
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[1], list4[0], list4[1]));
                        other.Vertices.Add(this.interp_vertex(list3[1], list3[3], list4[1], list4[3]));
                        other.Vertices.Add(this.interp_vertex(list3[2], list3[3], list4[2], list4[3]));
                        other.Faces.AddFace(0, 1, 2);
                        other.Vertices.Add(this.interp_vertex(list3[0], list3[2], list4[0], list4[2]));
                        other.Faces.AddFace(0, 3, 2);
                        break;

                    case 7:
                    case 8:
                        other.Vertices.Add(this.interp_vertex(list3[3], list3[0], list4[3], list4[0]));
                        other.Vertices.Add(this.interp_vertex(list3[3], list3[2], list4[3], list4[2]));
                        other.Vertices.Add(this.interp_vertex(list3[3], list3[1], list4[3], list4[1]));
                        other.Faces.AddFace(0, 1, 2);
                        break;
                }
                mesh.Append(other);
            }
            mesh.Vertices.CombineIdentical(true, true);
            return mesh;
        }
        private double calculate_field(Point3d test_pt)
        {
            double num2 = 0.0;
            int num5 = this.use_charges.Count - 1;
            for (int i = 0; i <= num5; i++)
            {
                double x = test_pt.DistanceTo(this.use_charges[i]);
                if (x < this.use_radii[i])
                {
                    num2 += ((((this.cA * Math.Pow(x, 3.0)) / Math.Pow(this.use_radii[i], 3.0)) + ((this.cB * Math.Pow(x, 2.0)) / Math.Pow(this.use_radii[i], 2.0))) + ((this.cC * x) / this.use_radii[i])) + 1.0;
                }
            }
            if (num2 == 0.0)
            {
                num2 = -1.0;
            }
            return num2;
        }
    }
}
