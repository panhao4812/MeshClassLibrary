using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GeoTools;

namespace Kangaroo
{
    public class Pressure : GoalObject
    {
        public class LaplacianSmooth : GoalObject
        {
            public double Strength;

            public LaplacianSmooth()
            {
            }

            public LaplacianSmooth(int[] P, double k)
            {
                PIndex = P;
                Move = new Vector3d[P.Length];
                Weighting = new double[P.Length];
                Strength = k;
            }

            public LaplacianSmooth(Point3d[] P, double k)
            {
                PPos = P;
                Move = new Vector3d[P.Length];
                Weighting = new double[P.Length];
                Strength = k;
            }

            public override void Calculate(List<Particle> p)
            {
                Point3d Avg = new Point3d();
                for (int i = 1; i < PIndex.Length; i++)
                {
                    Avg = Avg + p[PIndex[i]].Position;
                }
                double Inv = 1.0 / (PIndex.Length - 1);
                Avg = Avg * Inv;
                Vector3d Smooth = 0.5 * (Avg - p[PIndex[0]].Position);
                Move[0] = Smooth;
                Weighting[0] = Strength;
                Smooth *= -Inv;
                for (int i = 1; i < PIndex.Length; i++)
                {
                    Move[i] = Smooth;
                    Weighting[i] = Strength;
                }
            }

        }
        public class PolygonArea : GoalObject
        {
            public double Strength;
            public double TargetArea;
            public PolygonArea()
            {
            }
            public PolygonArea(List<int> V, double Area, double k)
            {
                int L = V.Count;
                PIndex = V.ToArray();
                Move = new Vector3d[L];
                Weighting = new double[L];
                for (int i = 0; i < L; i++)
                {
                    Weighting[i] = k;
                }
                TargetArea = Area;
                Strength = k;
            }
            public PolygonArea(List<Point3d> V, double Area, double k)
            {
                int L = V.Count;
                PPos = V.ToArray();
                Move = new Vector3d[L];
                Weighting = new double[L];
                for (int i = 0; i < L; i++)
                {
                    Weighting[i] = k;
                }
                TargetArea = Area;
                Strength = k;
            }
            public override void Calculate(List<Particle> p)
            {
                double CurrentAreaDoubled = 0;
                double TotalLength = 0;
                int L = PIndex.Length;

                for (int i = 0; i < L; i++)
                {
                    CurrentAreaDoubled += Vector3d.CrossProduct(
                      (Vector3d)p[PIndex[i]].Position,
                      (Vector3d)p[PIndex[(i + 1) % L]].Position).Z;
                    //note - points must be ordered CCW

                    TotalLength += p[PIndex[i]].Position.DistanceTo
                      (p[PIndex[(i + 1) % L]].Position);

                    Move[i] = Vector3d.Zero;
                }

                double AreaShortage = TargetArea - 0.5 * CurrentAreaDoubled;
                double Offset = AreaShortage / TotalLength;

                for (int i = 0; i < L; i++)
                {
                    int NextVert = (i + 1) % L;
                    Vector3d Edge = (p[PIndex[(i + 1) % L]].Position - p[PIndex[i]].Position);
                    Edge.Unitize();
                    Vector3d Pressure = Offset * Vector3d.CrossProduct(Edge, Vector3d.ZAxis);
                    Move[i] += Pressure;
                    Move[NextVert] += Pressure;
                    Weighting[i] = Strength;
                }
            }
        }
        public Pressure(Point3d[] Pts, double k)
        {
            PPos = Pts;
            Move = new Vector3d[3];
            Weighting = new double[3] { k, k, k };
        }
        public override void Calculate(List<Particle> p)
        {
            Point3d PA = p[PIndex[0]].Position;
            Point3d PB = p[PIndex[1]].Position;
            Point3d PC = p[PIndex[2]].Position;

            Vector3d AB = PB - PA;
            Vector3d BC = PC - PB;
            Vector3d CA = PA - PC;

            Vector3d Normal = Vector3d.CrossProduct(AB, BC); //this gives us a vector normal to the triangle, and of length twice its area            
            //halve this and divide it evenly over the 3 vertices. TODO - add Cotan weighting option here
            Vector3d PressureVector = (1.0 / 6.0) * Normal;

            Move[0] = Move[1] = Move[2] = PressureVector;
        }
    }
    public class PlasticHinge : GoalObject
    {
        public double RestAngle;
        public double PlasticLimit;
        public PlasticHinge(Point3d P0, Point3d P1, Point3d P2, Point3d P3, double Angle, double Plastic, double K)
        {
            PPos = new Point3d[4] { P0, P1, P2, P3 };
            Move = new Vector3d[4];
            Weighting = new double[4] { K, K, K, K };
            RestAngle = Angle;
            PlasticLimit = Plastic;
        }
        public override void Calculate(List<Particle> p)
        {
            Point3d P0 = p[PIndex[0]].Position;
            Point3d P1 = p[PIndex[1]].Position;
            Point3d P2 = p[PIndex[2]].Position;
            Point3d P3 = p[PIndex[3]].Position;

            Vector3d V01 = P1 - P0;
            Vector3d V02 = P2 - P0;
            Vector3d V03 = P3 - P0;
            Vector3d V21 = P1 - P2;
            Vector3d V31 = P1 - P3;

            double L01 = V01.Length;
            double invL01 = 1.0 / L01;

            double H0 = (V02 - V02 * invL01 * invL01 * V01 * V01).Length;
            double H1 = (V03 - V03 * invL01 * invL01 * V01 * V01).Length;
            double H = 0.5 / (H0 + H1);

            double Dot0201 = V02 * V01;
            double Dot0301 = V03 * V01;
            double Dot2101 = V21 * V01;
            double Dot3101 = V31 * V01;

            Vector3d Cross0 = Vector3d.CrossProduct(V02, V01);
            Vector3d Cross1 = Vector3d.CrossProduct(V01, V03);

            double CurrentAngle = Vector3d.VectorAngle(Cross0, Cross1, new Plane(P0, V01));
            if (CurrentAngle > Math.PI) { CurrentAngle = CurrentAngle - 2 * Math.PI; }

            double AngleError = CurrentAngle - RestAngle;

            double OverFold = Math.Abs(AngleError) - PlasticLimit; // the amount of folding beyond the elastic/plastic threshold
            if (OverFold > 0)
            {
                if (AngleError > 0)
                {
                    RestAngle += OverFold;
                    AngleError -= OverFold;
                }
                else
                {
                    RestAngle -= OverFold;
                    AngleError += OverFold;
                }
            }

            double InvL = 1.0 / Cross0.Length;
            double Cot0u = Dot0201 * InvL;
            double Cot0v = Dot2101 * InvL;
            InvL = 1.0 / Cross1.Length;
            double Cot1u = Dot0301 * InvL;
            double Cot1v = Dot3101 * InvL;

            double D = AngleError * H * 0.5;
            Vector3d A = Cross0 * D;
            Move[0] = Cot0v * A;
            Move[1] = Cot0u * A;
            Move[2] = -(Cot0v + Cot0u) * A;
            A = Cross1 * D;
            Move[0] += Cot1v * A;
            Move[1] += Cot1u * A;
            Move[3] = -(Cot1v + Cot1u) * A;
        }
    }
    /// <summary>
    /// Tries to make the 4 vertices of a quad lie on a common circle.
    /// Should be used in conjunction with planarize.
    /// Can be used for generating circular meshes, which have useful offset properties for beam structures.
    /// </summary>
    public class CyclicQuad : GoalObject
    {
        public double Strength;
        public CyclicQuad()
        {
        }
        public CyclicQuad(int[] P, double k)
        {
            PIndex = P;
            Move = new Vector3d[4];
            Weighting = new double[4];
            Strength = k;
        }
        public CyclicQuad(Point3d[] P, double k)
        {
            PPos = P;
            Move = new Vector3d[4];
            Weighting = new double[4];
            Strength = k;
        }
        public override void Calculate(List<Particle> p)
        {
            Point3d P1 = p[PIndex[0]].Position;
            Point3d P2 = p[PIndex[1]].Position;
            Point3d P3 = p[PIndex[2]].Position;
            Point3d P4 = p[PIndex[3]].Position;

            Point3d Center =
              (
              (new Circle(P1, P2, P3)).Center +
              (new Circle(P2, P3, P4)).Center +
              (new Circle(P3, P4, P1)).Center +
              (new Circle(P4, P1, P2)).Center
              ) * 0.25;
            double D1 = Center.DistanceTo(P1);
            double D2 = Center.DistanceTo(P2);
            double D3 = Center.DistanceTo(P3);
            double D4 = Center.DistanceTo(P4);
            double AvgDist = 0.25 * (D1 + D2 + D3 + D4);

            double stretchfactor;
            stretchfactor = 1.0 - AvgDist / D1;
            Move[0] = (Center - P1) * stretchfactor;
            stretchfactor = 1.0 - AvgDist / D2;
            Move[1] = (Center - P2) * stretchfactor;
            stretchfactor = 1.0 - AvgDist / D3;
            Move[2] = (Center - P3) * stretchfactor;
            stretchfactor = 1.0 - AvgDist / D4;
            Move[3] = (Center - P4) * stretchfactor;

            Weighting[0] = Strength;
            Weighting[1] = Strength;
            Weighting[2] = Strength;
            Weighting[3] = Strength;
        }
    }
    /// <summary>
    /// Bending resistance between a pair of triangles
    /// Based on ideas from
    /// "Discrete Shells" by Grinspun et al
    /// http://www.cs.columbia.edu/cg/pdfs/10_ds.pdf
    /// and
    /// "Interactive Form-Finding of Elastic Origami" by Tachi
    /// http://www.tsg.ne.jp/TT/cg/ElasticOrigami_Tachi_IASS2013.pdf
    /// </summary>
    public class Hinge : GoalObject
    {
        public double Strength;
        public double RestAngle;
        public Hinge()
        {
        }
        public Hinge(int P0, int P1, int P2, int P3, double RA, double k)
        {
            PIndex = new int[4] { P0, P1, P2, P3 };
            Move = new Vector3d[4];
            Weighting = new double[4];
            Strength = k;
            RestAngle = RA;
        }
        public Hinge(Point3d P0, Point3d P1, Point3d P2, Point3d P3, double RA, double k)
        {
            PPos = new Point3d[4] { P0, P1, P2, P3 };
            PIndex = new int[4];
            Move = new Vector3d[4];
            Weighting = new double[4];
            Strength = k;
            RestAngle = RA;
        }
        public override void Calculate(List<Particle> p)
        {
            Point3d P0 = p[PIndex[0]].Position;
            Point3d P1 = p[PIndex[1]].Position;
            Point3d P2 = p[PIndex[2]].Position;
            Point3d P3 = p[PIndex[3]].Position;

            Vector3d V01 = P1 - P0;
            Vector3d V02 = P2 - P0;
            Vector3d V03 = P3 - P0;
            Vector3d V21 = P1 - P2;
            Vector3d V31 = P1 - P3;

            double L01 = V01.Length;
            double invL01 = 1.0 / L01;
            //get heights

            //there is some re-use possible here to speed up
            double H0 = (V02 - V02 * invL01 * invL01 * V01 * V01).Length;
            double H1 = (V03 - V03 * invL01 * invL01 * V01 * V01).Length;
            double H = 0.5 / (H0 + H1);

            double Dot0201 = V02 * V01;
            double Dot0301 = V03 * V01;
            double Dot2101 = V21 * V01;
            double Dot3101 = V31 * V01;
            //Cot is dot over mag of cross

            //Get normals
            Vector3d Cross0 = Vector3d.CrossProduct(V02, V01);
            Vector3d Cross1 = Vector3d.CrossProduct(V01, V03);

            double CurrentAngle = Vector3d.VectorAngle(Cross0, Cross1, new Plane(P0, V01));
            if (CurrentAngle > Math.PI) { CurrentAngle = CurrentAngle - 2 * Math.PI; }
            double AngleError = CurrentAngle - RestAngle;
            double InvL = 1.0 / Cross0.Length;
            double Cot0u = Dot0201 * InvL;
            double Cot0v = Dot2101 * InvL;
            InvL = 1.0 / Cross1.Length;
            double Cot1u = Dot0301 * InvL;
            double Cot1v = Dot3101 * InvL;

            double D = AngleError * H * 0.5;
            Vector3d A = Cross0 * D;
            Move[0] = Cot0v * A;
            Move[1] = Cot0u * A;
            Move[2] = -(Cot0v + Cot0u) * A;
            A = Cross1 * D;
            Move[0] += Cot1v * A;
            Move[1] += Cot1u * A;
            Move[3] = -(Cot1v + Cot1u) * A;

            Weighting[0] = Strength;
            Weighting[1] = Strength;
            Weighting[2] = Strength;
            Weighting[3] = Strength;
        }
    }
    public class TangentialSmooth : GoalObject
    {
        public double Strength;

        public TangentialSmooth(Point3d[] P, double k)
        {
            PPos = P;
            Move = new Vector3d[P.Length];
            Weighting = new double[P.Length];
            Strength = k;
        }

        public override void Calculate(List<Particle> p)
        {
            Point3d Avg = new Point3d();
            var Vecs = new Vector3d[PIndex.Length - 1];
            for (int i = 1; i < PIndex.Length; i++)
            {
                Avg = Avg + p[PIndex[i]].Position;
                Vecs[i - 1] = p[PIndex[i]].Position - p[PIndex[0]].Position;
            }
            double Inv = 1.0 / (PIndex.Length - 1);
            Avg = Avg * Inv;
            Vector3d Smooth = 0.5 * (Avg - p[PIndex[0]].Position);

            Vector3d Normal = new Vector3d();
            for (int i = 0; i < Vecs.Length; i++)
            {
                Normal += Vector3d.CrossProduct(Vecs[i], Vecs[(i + 1) % Vecs.Length]);
            }
            Normal.Unitize();
            Smooth -= Normal * (Normal * Smooth);

            Move[0] = Smooth;
            Weighting[0] = Strength;
            Smooth *= -Inv;
            for (int i = 1; i < PIndex.Length; i++)
            {
                Move[i] = Smooth;
                Weighting[i] = Strength;
            }
        }
    }
    public class TangentIncircles : GoalObject
    {
        public double K;

        public TangentIncircles()
        {
        }

        public TangentIncircles(int P0, int P1, int P2, int P3, double Strength)
        {
            PIndex = new int[4] { P0, P1, P2, P3 };
            Move = new Vector3d[4];
            Weighting = new double[4];
            K = Strength;
        }

        public TangentIncircles(Point3d P0, Point3d P1, Point3d P2, Point3d P3, double Strength)
        {
            PPos = new Point3d[4] { P0, P1, P2, P3 };
            PIndex = new int[4];
            Move = new Vector3d[4];
            Weighting = new double[4];
            K = Strength;
        }

        public override void Calculate(List<Particle> p)
        {
            Point3d P0 = p[PIndex[0]].Position;
            Point3d P1 = p[PIndex[1]].Position;
            Point3d P2 = p[PIndex[2]].Position;
            Point3d P3 = p[PIndex[3]].Position;

            Vector3d V0 = P1 - P0;
            Vector3d V1 = P2 - P1;
            Vector3d V2 = P3 - P2;
            Vector3d V3 = P0 - P3;

            double L0 = V0.Length;
            double L1 = V1.Length;
            double L2 = V2.Length;
            double L3 = V3.Length;

            double L0L2 = L0 + L2;
            double L1L3 = L1 + L3;

            double MeanSum = 0.5 * (L0L2 + L1L3);

            double Stretch02 = 0.5 * (L0L2 - MeanSum);
            double Stretch13 = 0.5 * (L1L3 - MeanSum);

            V0.Unitize();
            V1.Unitize();
            V2.Unitize();
            V3.Unitize();

            Vector3d M0 = V0 * Stretch02;
            Vector3d M1 = V1 * Stretch13;
            Vector3d M2 = V2 * Stretch02;
            Vector3d M3 = V3 * Stretch13;

            Move[0] = M0 - M3;
            Move[1] = M1 - M0;
            Move[2] = M2 - M1;
            Move[3] = M3 - M2;

            Weighting[0] = Weighting[1] = Weighting[2] = Weighting[3] = K;
        }
    }
    /// <summary>
    /// Simple isotropic soap film element, as used for finding minimal surfaces
    /// Equivalent to cotan weighted Laplacian smoothing
    /// http://people.bath.ac.uk/abscjkw/LectureNotes/LightweightStructures/OtherMaterial/SoapFilmElement.pdf
    /// http://www.cs.cmu.edu/~kmcrane/Projects/DGPDEC/paper.pdf
    /// </summary>
    public class SoapFilm : GoalObject
    {
        public double Strength;

        public SoapFilm(Point3d[] Pts, double k)
        {
            PPos = Pts;
            Move = new Vector3d[3];
            Weighting = new double[3] { k, k, k };
            Strength = k;
        }

        public override void Calculate(List<Particle> p)
        {
            Point3d PA = p[PIndex[0]].Position;
            Point3d PB = p[PIndex[1]].Position;
            Point3d PC = p[PIndex[2]].Position;

            Vector3d AB = PB - PA;
            Vector3d BC = PC - PB;
            Vector3d CA = PA - PC;

            Vector3d Normal = Vector3d.CrossProduct(AB, BC);
            Normal.Unitize();

            Vector3d V0 = 0.5 * Vector3d.CrossProduct(BC, Normal);
            Vector3d V1 = 0.5 * Vector3d.CrossProduct(CA, Normal);
            Move[0] = V0;
            Move[1] = V1;
            Move[2] = -V0 - V1;

            Weighting[0] = Strength;
            Weighting[1] = Strength;
            Weighting[2] = Strength;
        }
    }
}
