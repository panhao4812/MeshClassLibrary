using Rhino;
using Rhino.Geometry;

using System;
using System.Collections.Generic;

namespace MeshClassLibrary
{
    public class IntPolyline
    {
        public int From = -1;
        public int To = -1;
        private int[] _MidPoints;
        private int MidLength = 0;
        public int[] ToArray()
        {
            int[] output = new int[MidLength + 2];
            output[0] = From;
            if (this.MidLength > 0)
            {
                for (int i = 0; i < MidLength; i++)
                {
                    output[i + 1] = this._MidPoints[i];
                }
            }
            output[output.Length - 1] = To;
            return output;
        }
        public IntPolyline(IntPolyline pl2)
        {
            this.From = pl2.From;
            this.To = pl2.To;
            this._MidPoints = pl2.MidPoints;
            this.MidLength = pl2.GetMidLength();
        }
        public IntPolyline(int MidLength)
        {
            _MidPoints = new int[MidLength];
            for (int i = 0; i < _MidPoints.Length; i++)
            {
                _MidPoints[i] = -1;
            }
        }
        public IntPolyline(int _from, int _to, int MidLength)
        {
            From = _from; To = _to;
            _MidPoints = new int[MidLength];///{0,0,0,0...}
            for (int i = 0; i < _MidPoints.Length; i++)
            {
                _MidPoints[i] = -1;
            }
        }
        public void Add(int pt)
        {
            _MidPoints[MidLength] = pt;
            MidLength++;
        }
        public void AddRange(int[] pt)
        {
            int sign = 0;
            for (int i = MidLength; i < _MidPoints.Length; i++)
            {
                if (sign < pt.Length)
                {
                    _MidPoints[i] = pt[sign];
                    sign++; MidLength++;
                }
                else
                {
                    return;
                }
            }
        }
        public void AddRange(IntPolyline pt)
        {
            int sign = 0;
            for (int i = MidLength; i < _MidPoints.Length; i++)
            {
                if (sign < pt.GetMidLength())
                {
                    _MidPoints[i] = pt.MidPoints[sign];
                    sign++; MidLength++;
                }
                else
                {
                    return;
                }
            }
        }
        public int[] MidPoints
        {
            get { return _MidPoints; }
            set
            {
                MidLength = 0;
                for (int i = 0; i < _MidPoints.Length; i++)
                {
                    if (i < value.Length)
                    {
                        _MidPoints[i] = value[i];
                        MidLength++;
                    }
                    else
                    {
                        _MidPoints[i] = -1;
                    }
                }
            }
        }
        public void AddRange(List<int> pt)
        {
            AddRange(pt.ToArray());
        }
        public void Reverse()
        {
            int sign = From;
            this.From = this.To; this.To = sign;
            int[] _MidPoints2 = new int[this._MidPoints.Length];
            for (int i = 0; i < this._MidPoints.Length; i++)
            {
                if (i < this.MidLength)
                {
                    _MidPoints2[i] = this._MidPoints[this.MidLength - 1 - i];
                }
                else
                {
                    _MidPoints2[i] = -1;
                }
            }
            this._MidPoints = _MidPoints2;

        }
        public int GetMidLength() { return this.MidLength; }
        public int Connet(IntPolyline pl2)
        {//this func does nothing for colosed polyline(return 4 or 5)
            if (pl2.GetMidLength() + this.GetMidLength() > this._MidPoints.Length) return 0;
            if (pl2.From == this.To)
            {
                if (pl2.To == this.From) return 5;
                this.To = pl2.To;
                this.Add(pl2.From);
                this.AddRange(pl2);

                return 1;
            }
            if (pl2.To == this.From)
            {
                if (pl2.From == this.To) return 5;
                pl2.To = this.To;
                pl2.Add(this.From);
                pl2.AddRange(this);
                this.From = pl2.From;
                this.To = pl2.To;
                this._MidPoints = pl2.MidPoints;
                this.MidLength = pl2.GetMidLength();
                return 2;
            }
            if (pl2.To == this.To)
            {
                if (pl2.From == this.From) return 6;
                pl2.Reverse();
                this.To = pl2.To;
                this.Add(pl2.From);
                this.AddRange(pl2);

                return 3;
            }
            if (pl2.From == this.From)
            {
                if (pl2.To == this.To) return 6;
                pl2.Reverse();
                pl2.To = this.To;
                pl2.Add(this.From);
                pl2.AddRange(this);
                this.From = pl2.From;
                this.To = pl2.To;
                this._MidPoints = pl2.MidPoints;
                this.MidLength = pl2.GetMidLength();
                return 4;
            }

            return 7;
        }
        public static string ConnectOpenIndex(List<IntPolyline> refer4, out int[] refer)
        {
            string _out = ""; refer = null;
            if (refer4.Count == 0) { _out += "input is empty"; return _out; }
            if (refer4.Count == 1) { refer = refer4[0].ToArray(); return _out; }
            for (int i = 0; i < refer4.Count - 1; i++)
            {
                for (int j = i + 1; j < refer4.Count; j++)
                {
                    int sign = refer4[j].Connet(refer4[i]);
                    if (sign == 1 || sign == 2 || sign == 3 || sign == 4) { break; }
                    if (sign == 5 || sign == 6) { _out = "closed"; return _out; }
                }
            }
            refer = refer4[refer4.Count - 1].ToArray();
            return _out;
        }
    }
    public class BasicVertice
    {
        /////////////////////basic
        public Point3d pos;
        public bool dead = false;
        public List<int> refer = new List<int>();
        public double energy = 0;
        public BasicVertice(Point3d p)
        {
            pos = new Point3d(p);
        }
        public BasicVertice(Point3d p, int index)
        {
            pos = new Point3d(p);
            this.refer.Add(index);
        }
        public void Add(int i)
        {
            this.refer.Add(i);
        }
        public bool equalTo(Point3d pt)
        {
            if (pos.DistanceTo(pt) < 0.01) { return true; }
            return false;
        }
        /// //////////////////static
        public static void CreateCollection(List<Line> x, out List<IndexPair> id, out List<BasicVertice> vs)
        {
            id = new List<IndexPair>(); vs = new List<BasicVertice>();
            id.Add(new IndexPair(0, 1));
            vs.Add(new BasicVertice(x[0].From, 1));
            vs.Add(new BasicVertice(x[0].To, 0));
            for (int i = 1; i < x.Count; i++)
            {
                bool sign1 = true;
                bool sign2 = true;
                int a = 0, b = 0;
                for (int j = 0; j < vs.Count; j++)
                {
                    if (vs[j].equalTo(x[i].From)) { sign1 = false; a = j; }
                    if (vs[j].equalTo(x[i].To)) { sign2 = false; b = j; }
                    if (!sign1 && !sign2) { break; }
                }
                if (sign1) { vs.Add(new BasicVertice(x[i].From)); a = vs.Count - 1; }
                if (sign2) { vs.Add(new BasicVertice(x[i].To)); b = vs.Count - 1; }
                vs[a].Add(b); vs[b].Add(a);
                id.Add(new IndexPair(a, b));
            }
        }
        public static List<Point3d> DisplayPos(List<BasicVertice> vs)
        {
            List<Point3d> output = new List<Point3d>();
            vs.ForEach(delegate (BasicVertice v) { output.Add(v.pos); });
            return output;
        }
        public static List<string> Displayenergy(List<BasicVertice> vs)
        {
            List<string> output = new List<string>();
            vs.ForEach(delegate (BasicVertice v) { output.Add(v.energy.ToString()); });
            return output;
        }
        public static List<string> DisplayLife(List<BasicVertice> vs)
        {
            List<string> output = new List<string>();
            vs.ForEach(delegate (BasicVertice v) { output.Add(v.dead.ToString()); });
            return output;
        }
    }
    public class Vertice1 : BasicVertice
    {
        /// static
        public Vertice1(Point3d p) : base(p) { }
        public static List<Point3d> DisplayPos(List<Vertice1> vs)
        {
            List<Point3d> output = new List<Point3d>();
            vs.ForEach(delegate (Vertice1 v) { output.Add(v.pos); });
            return output;
        }
        public static List<string> Displayenergy(List<Vertice1> vs)
        {
            List<string> output = new List<string>();
            vs.ForEach(delegate (Vertice1 v) { output.Add(v.energy.ToString()); });
            return output;
        }
        public static List<string> DisplayLife(List<Vertice1> vs)
        {
            List<string> output = new List<string>();
            vs.ForEach(delegate (Vertice1 v) { output.Add(v.dead.ToString()); });
            return output;
        }
        public static void CreateCollection(List<Line> x, out List<IndexPair> id, out List<Vertice1> vs)
        {
            id = new List<IndexPair>(); vs = new List<Vertice1>();
            id.Add(new IndexPair(0, 1));
            vs.Add(new Vertice1(x[0].From, 1));
            vs.Add(new Vertice1(x[0].To, 0));
            for (int i = 1; i < x.Count; i++)
            {
                bool sign1 = true;
                bool sign2 = true;
                int a = 0, b = 0;
                for (int j = 0; j < vs.Count; j++)
                {
                    if (vs[j].equalTo(x[i].From)) { sign1 = false; a = j; }
                    if (vs[j].equalTo(x[i].To)) { sign2 = false; b = j; }
                    if (!sign1 && !sign2) { break; }
                }
                if (sign1) { vs.Add(new Vertice1(x[i].From)); a = vs.Count - 1; }
                if (sign2) { vs.Add(new Vertice1(x[i].To)); b = vs.Count - 1; }
                vs[a].Add(b); vs[b].Add(a);
                id.Add(new IndexPair(a, b));
            }
        }
        /// ////////////////
        public Vertice1(Point3d p, int index) : base(p, index) { }
        public List<Polyline> edges = new List<Polyline>();
        public bool transferenergy(double percentage, ref List<Vertice1> vs)
        {
            bool sign = false;
            if (!this.dead && this.energy != 0)
            {
                this.dead = true;
                for (int i = 0; i < this.refer.Count; i++)
                {
                    if (vs[this.refer[i]].energy == 0)
                    {
                        vs[this.refer[i]].energy = this.energy * percentage;
                        sign = true;
                    }
                }
            }
            return sign;
        }
        public void CrateEdges(List<Vertice1> vs)
        {
            if (this.refer.Count == 3)
            {
                Point3d p1 = vs[this.refer[0]].pos; Vector3d v1 = p1 - this.pos; v1.Unitize();
                Point3d p2 = vs[this.refer[1]].pos; Vector3d v2 = p2 - this.pos; v2.Unitize();
                Point3d p3 = vs[this.refer[2]].pos; Vector3d v3 = p3 - this.pos; v3.Unitize();
                Plane p = new Plane(p1, p2, p3);
                Vector3d n = p.Normal; n.Unitize();
                Point3d N1 = this.pos + n * energy;
                Point3d N2 = this.pos - n * energy;
                Vector3d v0 = v1 + v2; v0.Unitize(); v0 *= this.energy;
                Point3d p12 = this.pos + v0;
                v0 = v2 + v3; v0.Unitize(); v0 *= this.energy;
                Point3d p23 = this.pos + v0;
                v0 = v3 + v1; v0.Unitize(); v0 *= this.energy;
                Point3d p31 = this.pos + v0;
                Polyline pl1 = new Polyline(); pl1.Add(N1); pl1.Add(p12); pl1.Add(N2); pl1.Add(p31);
                Polyline pl2 = new Polyline(); pl2.Add(N1); pl2.Add(p23); pl2.Add(N2); pl2.Add(p12);
                Polyline pl3 = new Polyline(); pl3.Add(N1); pl3.Add(p23); pl3.Add(N2); pl3.Add(p31);
                edges.Add(pl1); edges.Add(pl2); edges.Add(pl3);
            }
        }
    }
    public class Vertice2 : BasicVertice
    {
        /// static
        public Vertice2(Point3d p) : base(p) { }
        public static List<Point3d> DisplayPos(List<Vertice2> vs)
        {
            List<Point3d> output = new List<Point3d>();
            vs.ForEach(delegate (Vertice2 v) { output.Add(v.pos); });
            return output;
        }
        public static List<string> Displayenergy(List<Vertice2> vs)
        {
            List<string> output = new List<string>();
            vs.ForEach(delegate (Vertice2 v) { output.Add(v.energy.ToString()); });
            return output;
        }
        public static List<string> DisplayLife(List<Vertice2> vs)
        {
            List<string> output = new List<string>();
            vs.ForEach(delegate (Vertice2 v) { output.Add(v.dead.ToString()); });
            return output;
        }
        public static List<string> DisplayRef(List<Vertice2> vs)
        {
            List<string> output = new List<string>();
            vs.ForEach(delegate (Vertice2 v)
            {
                string str = "";
                for (int i = 0; i < v.refer.Count; i++)
                {
                    str += v.refer[i].ToString() + "/";
                }
                output.Add(str);
            });
            return output;
        }
        public static void CreateCollection(List<Line> x, out List<IndexPair> id, out List<Vertice2> vs)
        {
            id = new List<IndexPair>(); vs = new List<Vertice2>();
            id.Add(new IndexPair(0, 1));
            vs.Add(new Vertice2(x[0].From, 1));
            vs.Add(new Vertice2(x[0].To, 0));
            for (int i = 1; i < x.Count; i++)
            {
                bool sign1 = true;
                bool sign2 = true;
                int a = 0, b = 0;
                for (int j = 0; j < vs.Count; j++)
                {
                    if (vs[j].equalTo(x[i].From)) { sign1 = false; a = j; }
                    if (vs[j].equalTo(x[i].To)) { sign2 = false; b = j; }
                    if (!sign1 && !sign2) { break; }
                }
                if (sign1) { vs.Add(new Vertice2(x[i].From)); a = vs.Count - 1; }
                if (sign2) { vs.Add(new Vertice2(x[i].To)); b = vs.Count - 1; }
                vs[a].Add(b); vs[b].Add(a);
                id.Add(new IndexPair(a, b));
            }
        }
        /// ////////////////   
        public Vertice2(Point3d p, int index) : base(p, index) { }
        Vector3d N = Vector3d.ZAxis;
        public void computeNormal(Vector3d v)
        {
            N = v;
        }
        public void computeNormal(Surface s)
        {
            double u, v;
            s.ClosestPoint(this.pos, out u, out v);
            N = s.NormalAt(u, v);
        }
        public void computeNormal(Mesh s)
        {
            Point3d outpt; Vector3d outNormal;
            int output = s.ClosestPoint(this.pos, out outpt, out outNormal, double.MaxValue);
            if (output == -1) { N = new Vector3d(0, 0, 1); }
            else
            {
                N = outNormal;
            }
        }
        public static List<Polyline> Remesh(List<Vertice2> vs)
        {
            List<Polyline> output = new List<Polyline>();
            List<List<IndexPair>> index = new List<List<IndexPair>>();
            List<List<bool>> sign = new List<List<bool>>();
            int TCount = 0;
            for (int i = 0; i < vs.Count; i++)
            {
                List<IndexPair> children = new List<IndexPair>();
                List<bool> sign2 = new List<bool>();
                for (int j = 0; j < vs[i].refer.Count; j++)
                {
                    TCount++;
                    int after = j + 1;
                    if (after == vs[i].refer.Count) after = 0;
                    children.Add(new IndexPair(vs[i].refer[j], vs[i].refer[after]));

                    //   Print(vs[i].refer[j].ToString() + "~" + i.ToString() + "~" + vs[i].refer[after].ToString());
                    sign2.Add(true);
                }
                index.Add(children);
                sign.Add(sign2);
            }
            ////////////////////////////////////////////////////////
            for (int i = 0; i < vs.Count; i++)
            {
                for (int j = 0; j < vs[i].refer.Count; j++)
                {
                    if (sign[i][j])
                    {
                        sign[i][j] = false;
                        Polyline pl = new Polyline();
                        pl.Add(vs[i].pos);
                        ///to find a start vertice to construct polyline
                        bool signST = true;
                        int before = i;
                        int next = index[before][j].J;
                        //string error = "";
                        for (int loop = 0; loop < TCount; loop++)
                        {
                            signST = false;
                            // error += next.ToString() + "-";
                            for (int k = 0; k < index[next].Count; k++)
                            {
                                if (index[next][k].I == before)
                                {
                                    if (sign[next][k])
                                    {
                                        sign[next][k] = false;
                                        pl.Add(vs[next].pos);
                                        before = next;
                                        next = index[before][k].J;
                                        signST = true;
                                        break;
                                    }
                                }
                            }
                            if (!signST) break;
                        }
                        //    Print(error);
                        /////////////////
                        if (pl.Count > 2) pl.Add(pl[0]);
                        output.Add(pl);
                    }//if
                }//j
            }//i
            return output;
        }
        public void Sort(List<Vertice2> vs)
        { //sort the refer points in clockwise order
            List<IndexPair> refpoints = new List<IndexPair>();
            if (this.refer.Count < 3) return;
            Plane p1 = new Plane(this.pos, this.N);
            Plane p2 = new Plane(new Point3d(0, 0, 0), new Vector3d(0, 0, 1));
            for (int i = 0; i < this.refer.Count; i++)
            {
                Point3d P = new Point3d(vs[this.refer[i]].pos);
                P.Transform(Transform.PlaneToPlane(p1, p2));
                Vector3d v = new Vector3d(P.X, P.Y, 0);
                double t = 0;
                if (P.Y >= 0) { t = Vector3d.VectorAngle(new Vector3d(1, 0, 0), v); }
                else { t = Math.PI * 2 - Vector3d.VectorAngle(new Vector3d(1, 0, 0), v); }
                refpoints.Add(new IndexPair(this.refer[i], (int)(t * 1000000)));
            }
            refpoints.Sort(CompareDinosByLength);
            this.refer.Clear();
            for (int i = 0; i < refpoints.Count; i++)
            {
                this.refer.Add(refpoints[i].I);
            }
        }
        private static int CompareDinosByLength(IndexPair x, IndexPair y)
        {

            if (x.J > y.J) return 1;
            if (x.J == y.J) return 0;
            if (x.J < y.J) return -1;
            else return 0;
        }
        private static int CompareDinosByLength2(Polyline x, Polyline y)
        {
            if (x.Count > y.Count) return 1;
            if (x.Count == y.Count) return 0;
            if (x.Count < y.Count) return -1;
            else return 0;
        }
        public static List<Vertice2> CleanEdge(List<Vertice2> vs)
        {
            for (int k = 0; k < vs.Count; k++)
            {
                for (int i = 0; i < vs.Count; i++)
                {
                    int iCount = 0;
                    for (int j = 0; j < vs[i].refer.Count; j++)
                    {
                        if (vs[vs[i].refer[j]].dead == false) iCount++;
                    }
                    if (iCount < 2) { vs[i].dead = true; continue; }
                }
            }
            return Clean(vs);
        }
        public static List<Vertice2> Clean(List<Vertice2> vs)
        {
            List<Vertice2> vs2 = new List<Vertice2>();
            List<int> index = new List<int>();
            for (int i = 0; i < vs.Count; i++)
            {
                if (vs[i].dead == false) { index.Add(vs2.Count); vs2.Add(vs[i]); }
                else { index.Add(-1); }
            }
            for (int i = 0; i < vs2.Count; i++)
            {
                List<int> refer2 = new List<int>();
                for (int j = 0; j < vs2[i].refer.Count; j++)
                {
                    int indexn = index[vs2[i].refer[j]];
                    if (indexn != -1) refer2.Add(indexn);
                }
                vs2[i].refer = refer2;
            }
            return vs2;
        }
    }
    public class Vertice3 : BasicVertice
    {
        double K = 0.062;
        double F = 0.062;
        public double U = 1, V = 0;
        public double dU = 0, dV = 0;
        //lapU means laplace equation
        //u and v means two different kinds of chemical solution.
        //We always use a energy number to define the density.
        public Vertice3(Point3d p) : base(p) { }
        public static void CreateCollection(Mesh mesh, out List<Vertice3> vs)
        {
            Rhino.Geometry.Collections.MeshTopologyVertexList vs1 = mesh.TopologyVertices;
            vs = new List<Vertice3>();
            for (int i = 0; i < vs1.Count; i++)
            {
                vs.Add(new Vertice3(new Point3d(vs1[i].X, vs1[i].Y, vs1[i].Z)));
            }
            for (int i = 0; i < vs1.Count; i++)
            {
                vs[i].refer = new List<int>(vs1.ConnectedTopologyVertices(i));
            }
        }
        public void ComputeLaplation1(List<Vertice3> vs)
        {
            double lapU = 0, lapV = 0;
            for (int i = 0; i < this.refer.Count; i++)
            {
                lapU += vs[this.refer[i]].U;
                lapV += vs[this.refer[i]].V;
            }
            lapU -= U * this.refer.Count;
            lapV -= V * this.refer.Count;
            lapU *= 0.19; lapV *= 0.08;
            dU = lapU - U * V * V + F * (1 - U);
            dV = lapV + U * V * V - (K + F) * V;
        }
        public void ComputeLaplation2(List<Vertice3> vs)
        {
            double lapU = 0, lapV = 0;
            double tot = 0;
            for (int i = 0; i < this.refer.Count; i++)
            {
                double t1 = vs[refer[i]].pos.DistanceTo(this.pos);
                lapU += vs[this.refer[i]].U * t1;
                lapV += vs[this.refer[i]].V * t1;
                tot += t1;
            }
            lapU /= tot; lapU -= U;
            lapV /= tot; lapV -= V;
            lapU *= 0.19 * 2; lapV *= 0.08 * 2;
            dU = lapU - U * V * V + F * (1 - U);
            dV = lapV + U * V * V - (K + F) * V;
        }
        public void ComputeLaplation3(List<Vertice3> vs)
        {
            double lapU = 0, lapV = 0;
            double tot = 0;
            for (int i = 0; i < this.refer.Count; i++)
            {
                double t1 = vs[refer[i]].pos.DistanceTo(this.pos);
                lapU += vs[this.refer[i]].U * 0.1 / t1;
                lapV += vs[this.refer[i]].V * 0.1 / t1;
                tot += 0.1 / t1;
            }
            lapU /= tot; lapU -= U;
            lapV /= tot; lapV -= V;
            lapU *= 0.19 * 2; lapV *= 0.085 * 2;
            dU = lapU - U * V * V + F * (1 - U);
            dV = lapV + U * V * V - (K + F) * V;
        }
        public void ComputeUV1()
        {
            this.U += dU;
            this.V += dV;
        }
    }
    public class Vertice4 : BasicVertice
    {//simplify mesh
        public bool isNaked = true;
        public Vertice4(Point3d p) : base(p) { }
        private List<int> refer2 = new List<int>();
        private List<int> refer3 = new List<int>();
        private void SortRefer()
        {
            // string _out = "";
            if (!ConfirmRefer(refer2)) return;//_out += "error in refer2" + "\r\n"; return _out; }
            if (!ConfirmRefer(refer3)) return;//_out += "error in refer3" + "\r\n"; return _out; }
            List<IntPolyline> refer4 = new List<IntPolyline>();
            for (int j = 0; j < refer2.Count; j++)
            {
                refer4.Add(new IntPolyline(refer2[j], refer3[j], refer2.Count));
            }
            if (refer4.Count == 0) return;// _out = "error in refer4"; return _out; }
            if (refer4.Count == 1)
            {
                this.refer = new List<int>(refer4[0].ToArray());
                return;//_out += "Vertice4 is OK" + "\r\n"; return _out;
            }
            for (int i = 0; i < refer4.Count - 1; i++)
            {
                for (int j = i + 1; j < refer4.Count; j++)
                {
                    int sign = refer4[j].Connet(refer4[i]);
                    //mesh normals must be simplified.
                    if (sign == 5) { refer4[j] = refer4[i]; this.isNaked = false; }
                    if (sign == 1 || sign == 2) break;
                }
            }
            this.refer = new List<int>(refer4[refer4.Count - 1].ToArray());
            //_out += "Vertice4 is OK" + "\r\n";
            //////////////////////////////////////////////
            return;// _out;
        }
        private bool ConfirmRefer(List<int> refer4)
        {
            List<int> p1 = new List<int>(refer4);
            if (p1.Count <= 0) return false;
            if (p1.Count == 1) return true;
            p1.Sort(); int j = p1.Count - 1;
            for (int i = 0; i < p1.Count; i++)
            {
                if (p1[i] == p1[j]) return false;
                j = i;
            }
            return true;
        }
        //static function MeshSimplify has bugs that are not fixed.
        public static List<Line> MeshProfile(Mesh mesh)
        {
            List<Vertice4> v4s; Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            CreateCollection(mesh, out v4s);
            List<bool> edge = new List<bool>();
            for (int i = 0; i < el.Count; i++) { edge.Add(false); }
            for (int step = 0; step < v4s.Count; step++)
            {
                bool sign = false;
                for (int i = 0; i < v4s.Count; i++)
                {
                    //  v4s[i].p.Print(i.ToString() + "==>");
                    bool sign2 = Vertice4.TransEnergy(i, ref v4s, ref edge);
                    sign = sign || sign2;
                }
                // Print(sign.ToString());
                if (!sign && step > 0) { break; }
            }
            List<Line> output = new List<Line>();
            for (int i = 0; i < el.Count; i++)
            {
                if (edge[i]) { output.Add(el.EdgeLine(i)); }
                else
                {
                    int a = el.GetTopologyVertices(i).I;
                    int b = el.GetTopologyVertices(i).J;
                    if (v4s[a].isNaked && v4s[b].isNaked)
                    {
                        output.Add(el.EdgeLine(i));
                    }
                }
            }
            return output;
        }
        //public static List<Mesh> MeshSeperate(Mesh mesh){ }
        static bool[,] VerticeTable =   {
            {false,false,false,false,false},
            {true,true,false,true,false},
            {true,false,true,false,true},
            {true,true,true,true,true},
            {true,true,false,true,false},
            {false,true,false,true,false},
            {true,true,true,true,true},
            {true,true,true,true,true},
            {true,false,true,false,true},
            {true,true,true,true,true},
            {false,false,true,false,true},
            {true,true,true,true,true},
            {true,true,true,true,true},
            {true,true,true,true,true},
            {true,true,true,true,true},
            {false,true,true,true,true}
            };
        public static bool TransEnergy(int I, ref List<Vertice4> v4s, ref List<bool> signlist)
        {
            /*
             *There are two kings of vertice,the generator and transfer.
             * the generator is (naked && degree!=3) and  (!naked && degree!=4)
             * the transfer is (naked && degree!=4)
             * transfer activate with opposite edges.
             */
            Vertice4 V = v4s[I];
            if (V.isNaked)
            {
                if (V.refer.Count == 3) return false;
                for (int i = 0; i < V.refer.Count; i++)
                {
                    signlist[V.refer[i]] = true;
                }
                return false;
            }
            else
            {
                if (V.refer.Count != 4)
                {
                    for (int i = 0; i < V.refer.Count; i++)
                    {
                        signlist[V.refer[i]] = true;
                    }
                    return false;
                }
                if (V.refer.Count == 4)
                {
                    int sign = 0x00;
                    if (signlist[V.refer[0]])
                    {
                        sign |= 1;
                    }
                    if (signlist[V.refer[1]])
                    {
                        sign |= 2;
                    }
                    if (signlist[V.refer[2]])
                    {
                        sign |= 4;
                    }
                    if (signlist[V.refer[3]])
                    {
                        sign |= 8;
                    }
                    signlist[V.refer[0]] = VerticeTable[sign, 1];
                    signlist[V.refer[1]] = VerticeTable[sign, 2];
                    signlist[V.refer[2]] = VerticeTable[sign, 3];
                    signlist[V.refer[3]] = VerticeTable[sign, 4];

                    return VerticeTable[sign, 0];
                }
            }
            return false;
        }
        public static void CreateCollection(Mesh mesh, out List<Vertice4> v4s)
        {
            v4s = new List<Vertice4>();
            Rhino.Geometry.Collections.MeshTopologyVertexList vs = mesh.TopologyVertices;
            for (int i = 0; i < vs.Count; i++)
            {
                v4s.Add(new Vertice4(new Point3d(vs[i].X, vs[i].Y, vs[i].Z)));
            }
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                int[] index = vs.IndicesFromFace(i);
                if (index.Length == 4)
                {
                    v4s[index[0]].refer3.Add(index[3]); v4s[index[0]].refer2.Add(index[1]);
                    v4s[index[1]].refer3.Add(index[0]); v4s[index[1]].refer2.Add(index[2]);
                    v4s[index[2]].refer3.Add(index[1]); v4s[index[2]].refer2.Add(index[3]);
                    v4s[index[3]].refer3.Add(index[2]); v4s[index[3]].refer2.Add(index[0]);
                }
                if (index.Length == 3)
                {
                    v4s[index[0]].refer3.Add(index[2]); v4s[index[0]].refer2.Add(index[1]);
                    v4s[index[1]].refer3.Add(index[0]); v4s[index[1]].refer2.Add(index[2]);
                    v4s[index[2]].refer3.Add(index[1]); v4s[index[2]].refer2.Add(index[0]);
                }
            }
            Rhino.Geometry.Collections.MeshTopologyEdgeList el = mesh.TopologyEdges;
            for (int i = 0; i < v4s.Count; i++)
            {
                v4s[i].SortRefer();
                v4s[i].refer2.Clear();
                for (int j = 0; j < v4s[i].refer.Count; j++)
                {
                    v4s[i].refer2.Add(el.GetEdgeIndex(i, v4s[i].refer[j]));
                }
                v4s[i].refer = v4s[i].refer2;
            }
            return;
        }
        public static List<Rhino.Display.Text3d> DisplayRef(List<Vertice4> vs)
        {
            List<Rhino.Display.Text3d> output = new List<Rhino.Display.Text3d>();

            for (int i = 0; i < vs.Count; i++)
            {
                Point3d f = vs[i].pos;
                List<int> fs = vs[i].refer;
                List<int> fs2 = vs[i].refer2;
                List<int> fs3 = vs[i].refer3;
                string str = "";
                for (int j = 0; j < fs2.Count; j++)
                {
                    str += fs2[j].ToString() + "-";
                }
                Rhino.Display.Text3d te = new Rhino.Display.Text3d(str.ToString(), new Plane(f, Vector3d.ZAxis), 1);
                output.Add(te);
            }
            return output;
        }
    }
    #region random
    class M_Random
    {
        public int seed = -1;
        public M_Random() { }
        public double NextDouble()
        {
            seed++;
            if (seed == NL.Length) seed = 0;
            return NL[seed];
        }
        double[] NL = {0.84,0.13,0.14,0.23,0.76,0.54,0.67,0.59,0.91,0.74,
0.19,0.71,0.88,0.46,0.78,0.15,0.22,0.53,0.21,0.14,0.56,0.14,0.09,
0.75,0.33,0.67,0.26,0.24,0.22,0.04,0.03,0.96,0.71,0.13,0.43,0.94,
0.11,0.79,0.16,0.42,0.67,0.17,0.11,0.09,0.02,0.32,0.9,0.31,0.57,
0.05,0.94,0.6,0.6,0.82,0.04,0.7,0.05,0.39,0.9,0.09,0.28,0.43,0.37,
0.87,0.71,0.23,0,0.75,0.03,0.85,0.04,0.43,0.37,0.79,0.47,0.39,0.03,
1,0.73,0.02,0.77,0.95,0.67,0.17,0.11,0.43,0.36,0.89,0.09,0.73,0.89,
0.72,0.89,0.07,0.14,0.23,0.81,0.24,0.38,0.79,0.31,0.15,0.28,0.73,
0.01,0.51,0.23,0.81,0.35,0.65,0.67,0.05,0.67,0.88,0.32,0.33,0.77,
0.2,0.77,0.28,0.86,0.12,0.66,0.3,0.96,0.32,0.53,0.3,0.66,0.23,
0.58,0.79,0.62,0.93,0.7,0.62,0.67,0.94,0.16,0.6,0.2,0.56,0.53,
0.43,0.06,0.84,0.06,0.02,0.74,0.8,0.47,0.61,0.47,0.1,0.93,0.2,
0.49,0.99,0.77,0.7,0.97,0.92,0.15,0.12,0.07,0.88,0.43,0.73,0.17,
0.71,0.66,0.83,0.03,0.17,0.08,0.31,0.58,0.23,0.24,0.12,0.26,0.52,
0.56,0.86,0.77,0.97,0.32,0.52,1,0.51,0.12,0.68,0.17,0.46,0.62,
0.28,0.41,0.42,0.36,0.18,0.42,0.32,0.84,0.04,0.14,0.64,0.58,0.3,
0.02,0.63,0.62,0.26,0.74,0.65,0.44,0.46,0.37,0.29,0.35,0.1,0.56,
0.9,0.74,0.67,0.58,0.98,0.66,0.57,0.55,0.8,0.89,0.17,0.87,0.06,
0.7,0.94,0.68,0.52,0.72,0.13,0.4,0.02,0.5,0.37,0.89,0.86,0.94,
0.52,0.03,0.17,0.91,0.12,0.07,0.27,0.62,0.51,0.59,0.17,0.46,
0.16,0.98,0.01,0.75,0.23,0.73,0.45,0.4,0.68,0.95,0.5,0.78,0.85,
0.57,0.22,0.7,0.53,0.4,0.36,0.78,0.72,0.04,0.14,0.54,0.38,0.88,
0.77,0.1,0.6,0.44,0.19,0.35,0.5,0.06,0.56,0.15,0.39,0.27,0.28,
0.64,0.44,0.46,0.26,0.57,0.53,0.39,0.07,0.55,0.85,0.57,0.09,0.11,
0.23,0.39,0.74,0.12,0.85,0.46,0.37,0.35,0.96,0.35,0.8,0.24,0.77,
0.15,0.28,0.79,0.01,0.07,0.3,0.26,0.13,0.72,0.34,0.26,0.79,0.57,
0.02,0.99,0.81,0.22,0.25,0.03,0.35,0.07,0.12,0.11,0.32,0.45,0.3,
0.93,0.9,0.93,0.68,0.09,0.66,0.02,0.8,0.37,0.11,0.28,0.54,0.79,
0.26,0.83,0.99,0.51,0.05,0.48,0.33,0.28,0.45,0.38,0.53,0.74,0.1,
0.77,0.89,0.69,0.03,0.16,0.46,0.57,0.76,0.38,0.36,0.2,0.04,0.25,
0.6,0.77,0.77,0.64,0.88,0.53,0.68,0.46,0.76,0.52,0.08,0.62,0.06,
0.84,0.12,0.02,0.48,0.52,0.4,0.94,0.98,0.89,0.12,0.11,0.35,0.95,
0.82,0.97,0.02,0.89,0.47,0.79,0.47,0.8,0.88,0.56,0.51,0.8,0.5,0,
0.06,0.64,0.01,0.38,0.61,0.41,0.11,0.62,0.45,0.74,0.89,0.84,0.8,
0.1,0.27,0.71,0.65,0.66,0.29,0.93,0.72,0.7,0.44,0.87,0.04,0.29,
0.15,0.25,0.97,0.56,0.51,0.68,0.02,0.4,0.88,0.34,0.88,0.75,0.5,
0.94,0.84,0.2,0.52,0.29,1,0.64,0.99,0.37,0.54,0.16,0.91,0.85,
0.51,0.56,0.28,0.36,0.21,0.13,0.33,0.32,0.26,0.85,0.65,0.89,
0.23,0.21,0.81,0.4,0.22,0.93,0.83,0.9,0.16,0.35,0.09,0.52,0.18,
0.15,0.88,0.41,0.3,0.78,0.72,0.8,0.65,0.66,0.17,0.46,0.11,0.52,
0.13,0.75,0.41,0.18,0.68,0.99,0.55,0.63,0.06,0.78,0.82,0.59,0.15,
0.61,0.33,0.01,0.69,0.16,0.47,0.77,0.18,0.06,0.25,0.93,0.02,0.48,
0.13,0.85,0.25,0.57,0.05,0.35,0.29,0.71,0.79,0.09,0.49,0.98,0.67,
0.1,0.97,0.55,0.09,0.1,0.59,0.71,0.63,0.11,0.48,0.63,0.97,0.01,
0.99,0.35,0.34,0.08,0.49,0.49,0.16,0.2,0.86,0.7,0.38,0.5,0.73,
0.47,0.3,0.44,0.82,0.24,0.52,0.71,0.49,0.37,0.79,0.63,0.97,0.16,
0.34,0.31,0.85,0.03,0.37,0.61,0.6,0.04,0.36,0.94,0.37,0.72,0.6,0,
0.83,0.47,0.24,0.27,0.17,0.59,0.36,0.12,0.41,0.19,0.29,0.23,0.11,
0.26,0.52,0.62,0.55,0.71,0.11,0.33,0.14,0.84,0.35,0.83,0.33,0.77,
0.9,0.69,0.11,0.36,0.07,0.1,0.65,0.52,0.88,0.01,0.12,0.52,0.46,
0.38,0.79,0.22,0.91,0.65,0.74,0.14,0.5,0.34,0.52,0.74,0.39,0.66,
0.61,0.26,0.86,0.98,0.12,0.42,0.95,0.4,0.7,0.67,0.01,0.05,0.12,
0.18,0.59,0.59,0.38,0.51,0.5,0.03,0.25,0.73,0.54,0.92,0.94,0.7,
0.09,0.19,0.27,0.56,0.17,0.37,0.96,0.41,0.49,0.38,0.67,0.9,0.89,
0.71,0.58,0.07,0.68,0.12,0.21,0.86,0.53,0.56,0.55,0.92,0.96,0.01,
0.23,0.36,0.41,0.88,0.72,0.94,0.05,0.42,0.33,0.76,0.13,0.14,0.5,
0.64,0.09,0.71,0.69,0.21,0.92,0.48,0.61,0.79,0.45,0.18,0.05,0.41,
0.71,0.08,0.17,0.53,0.64,0.35,0.6,0.16,0.14,0.6,0.01,0.61,0.66,0.73,
0.86,0.47,0.38,0.82,0.94,0.54,0.63,0.57,0.77,0.82,0.87,0.34,1,0.49,
0.4,0.44,0.91,0.22,0.83,0.31,0.23,0.97,0.25,0.8,0.11,0.78,0.54,0.34,
0.5,0.49,0.7,0.09,0.55,0.19,0.62,0.14,0.42,0.63,0.25,0.51,0.79,0.13,
0.31,0.35,0.66,0.3,0.35,0.11,0.76,0.7,0.69,0.79,0.77,0.35,0.5,0.89,
0.22,0.57,0.71,0.16,1,0.29,0.07,0.28,0.12,0.78,0.8,0.8,0.87,0.26,
0.02,0.28,0.97,0.33,0.52,0.09,0.66,0.91,0.14,0.81,0.42,0.43,0.58,0.8,
0.8,0.91,0.31,0.2,0.69,0.73,0.92,0.84,0.93,0.09,0.51,0.5,0.06,0.03,
0.22,0.88,0.5,0.55,0.25,0.5,0.67,0.42,0.81,0.45,0.83,0.41,0.23,0.31,
0.43,0.89,0.74,0.57,0.71,0.27,0.49,0.21,0.46,0.6,0.12,0.14,0.34,0.18,
0.35,0.89,0.81,0.02,0.03,0.63,0.69,0.26,0.31,0.87,0.19,0.07,0.13,0.38,
0.1,0.87,0.37,0.28,0.5,0.61,0.42,0.03,0.35,0.94,0.78,0.79,0.54,0.01,
0.42,0.9,0.43,0.11,0.17,0.49,0.07,0.93,0.63,0.81,0.38,0.6,0.62,0.17,
0.58,0.86,0.38,0.64,0.15,0.1,0.12,0.6,0.23,0.84,0.64,0.72,0.76,0.32,
0.54,0.87,0.24,0.24,0.09,0.68,0.84,0.41,0.44,0.08,0.91,0.63,0.31,0.17,
0.23,0.56,0.9,0.9,0.99,0.25,0.45,0.48,0.56,0.14,0.64,0.44,0.89,0.82,
0.67,0.59,0.47,0.44,0.73,0.75,0.39,0.76,0.57,0.14,0.51,0.95,0.32,0.17,
0.42,0.3,0.73,0.51,0.79,0.95,0.37,0.68,0.94,0.74,0.73,0.51,0.87,0.06,
0.31,0.1,0.6,0.65,0.78,0.02,0.74};
    }
    #endregion
}