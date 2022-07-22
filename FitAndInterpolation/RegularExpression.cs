using System;
using System.Collections.Generic;
using System.Linq;

namespace FitAndInterpolation
{
    public struct IndexPair
    {
        public int I;
        public int J;
        public IndexPair(int i, int j)
        {
            I = i;
            J = j;
        }
    }
    public class RegularExpression
    {
        private string str;
        private void Print(string text)
        {
            str += text + "\r\n";
        }
        public void PrintPairList(List<IndexPair> pair)
        {
            for (int i = 0; i < pair.Count; i++)
            {
                Print(pair[i].I.ToString() + " " + pair[i].J.ToString());
            }
        }
        public void PrintList(List<int> pair)
        {
            for (int i = 0; i < pair.Count; i++)
            {
                Print(i.ToString() + " " + pair[i].ToString() + " ");
            }
        }
        public List<IndexPair> CombinePair(List<IndexPair> PairList)
        {
            List<int> list_int = new List<int>();
            for (int i = 0; i < PairList.Count; i++)
            {
                list_int.Add(PairList[i].I);
                list_int.Add(PairList[i].J);
            }
            list_int = list_int.Distinct<int>().ToList();

            List<IndexPair> ShortPairs = new List<IndexPair>();
            for (int i = 0; i < list_int.Count - 1; i++)
            {
                if (list_int[i + 1] == list_int[i] + 1) ShortPairs.Add(new IndexPair(list_int[i], list_int[i + 1]));
            }
            PairList.AddRange(ShortPairs);
            List<bool> sign = new List<bool>();
            for (int i = 0; i < list_int.Count; i++)
            {
                sign.Add(EdgeTest(PairList, list_int[i]));
            }
            List<int> list_output = new List<int>();
            for (int i = 0; i < list_int.Count; i++)
            {
                if (!sign[i]) list_output.Add(list_int[i]);
            }
            return FromPairArray(list_output);
        }
        private bool EdgeTest(List<IndexPair> pairs, int input)
        {
            double a = input - 0.1; double b = input + 0.1;
            bool signa = false; bool signb = false;
            foreach (IndexPair pair in pairs)
            {
                double I = (double)pair.I;
                double J = (double)pair.J;
                if (a <= J && a >= I) signa = true;
                if (b <= J && b >= I) signb = true;
                if (signa && signb) break;
            }
            return (signa && signb);
        }
        public int CompairPair(IndexPair x, IndexPair y)
        {
            if (x.I < y.I) return -1;
            if (x.I > y.I) return 1;
            return 0;
        }
        public List<IndexPair> FromIntArray(List<int> input)
        {
            List<IndexPair> output = new List<IndexPair>();
            input = input.Distinct<int>().ToList();
            for (int i = 0; i < input.Count - 1; i++)
            {
                output.Add(new IndexPair(input[i], input[i + 1]));
            }
            return output;
        }
        public List<IndexPair> FromPairArray(int[] input)
        {
            List<IndexPair> output = new List<IndexPair>();
            for (int i = 0; i < input.Length / 2.0; i++)
            {
                int a1 = input[i * 2], a2 = input[i * 2 + 1];
                if (a1 > a2)
                {
                    output.Add(new IndexPair(a2, a1));
                }
                else if (a2 > a1)
                {
                    output.Add(new IndexPair(a1, a2));
                }
            }
            return output;
        }
        public List<IndexPair> FromPIArray(int[] input)
        {
            List<IndexPair> output = new List<IndexPair>();
            for (int i = 0; i < input.Length / 2.0; i++)
            {
                int a1 = input[i * 2], a2 = input[i * 2 + 1];
                if (a1 > a2)
                {
                    if (a1 < 359)
                    {
                        output.Add(new IndexPair(a1, 359));
                    }
                    if (a2 > 0)
                    {
                        output.Add(new IndexPair(0, a2));
                    }
                }
                else if (a2 > a1)
                {
                    output.Add(new IndexPair(a1, a2));
                }
            }
            return output;
        }
        public List<IndexPair> FromPairArray(List<int> input)
        {
            return FromPairArray(input.ToArray());
        }
        //必须是I<J
        public string ExpressionPair(IndexPair pair)
        {
            List<IndexPair> pairs = SplitPair(pair);
            //string str = "^(";
            string str = "";
            for (int i = 0; i < pairs.Count; i++)
            {
                str += ExpressionShortPair(pairs[i]);
                if (i != (pairs.Count - 1)) str += "|";
            }
            //str += ")$";
            return str;
        }
        public string ExpressionShortPair(IndexPair pair)
        {
            char[] Char_a = pair.I.ToString().ToCharArray();
            char[] Char_b = pair.J.ToString().ToCharArray();
            string output = "";
            if (pair.I == pair.J) return pair.I.ToString();
            if (Char_a[Char_a.Length - 1] == '9')
            {
                output = pair.I.ToString() + "|";
                Char_a = (pair.I + 1).ToString().ToCharArray();
            }
            if (Char_a.Length != Char_b.Length) { return ""; }

            for (int i = 0; i < Char_a.Length; i++)
            {
                if (Char_a[i] == Char_b[i]) { output += Char_a[i]; }
                else { output += "[" + Char_a[i] + "-" + Char_b[i] + "]"; }
            }
            return output;
        }
        public List<IndexPair> SplitPair(IndexPair pair)
        {
            List<IndexPair> output = new List<IndexPair>();
            char[] Char_a = pair.I.ToString().ToCharArray();
            char[] Char_b = pair.J.ToString().ToCharArray();

            List<int> up_list = new List<int>();
            up_list.Add(pair.I);
            /*
            if (Char_a[Char_a.Length - 1] == '9')
            {
                Char_a = (pair.I + 1).ToString().ToCharArray();
            }
            */
            for (int i = 0; i < Char_b.Length; i++)
            {
                string db = "";
                for (int j = 0; j < Char_a.Length - i - 1; j++)
                {
                    db += Char_a[j];
                }
                for (int j = 0; j <= i; j++)
                {
                    db += "9";
                }
                up_list.Add(Convert.ToInt32(db));
            }
            List<int> up_list2 = new List<int>();
            for (int i = 0; i < up_list.Count; i++)
            {
                if (up_list[i] < pair.J) up_list2.Add(up_list[i]);
            }
            up_list2.Add(pair.J);
            up_list2 = up_list2.Distinct().ToList();
            for (int i = 0; i < up_list2.Count - 1; i++)
            {
                if (i == 0) { output.Add(new IndexPair(up_list2[i], up_list2[i + 1])); }
                else
                {
                    output.Add(new IndexPair(up_list2[i] + 1, up_list2[i + 1]));
                }
            }
            //  if (output.Count == 1) return output;
            IndexPair lastpair = output[output.Count - 1];
            output.RemoveAt(output.Count - 1);
            up_list.Clear();
            up_list.Add(pair.J);
            for (int i = 0; i < Char_b.Length; i++)
            {
                if (i == 0 && Char_b[Char_b.Length - 1] == '9') continue;
                string db = "";
                for (int j = 0; j < Char_b.Length - i - 1; j++)
                {
                    db += Char_b[j];
                }
                for (int j = 0; j <= i; j++)
                {
                    db += "0";
                }
                up_list.Add(Convert.ToInt32(db));
            }
            up_list2.Clear();
            for (int i = 0; i < up_list.Count; i++)
            {
                if (up_list[i] >= lastpair.I)
                {
                    up_list2.Add(up_list[i]);
                }
            }
            up_list2.Reverse();
            if (up_list2.Count == 1) { output.Add(new IndexPair(lastpair.I, up_list2[0])); return output; }
            if (up_list2[0] > lastpair.I) output.Add(new IndexPair(lastpair.I, up_list2[0] - 1));

            for (int i = 0; i < up_list2.Count - 1; i++)
            {
                if (i == up_list2.Count - 2) { output.Add(new IndexPair(up_list2[i], up_list2[i + 1])); }
                else
                {
                    output.Add(new IndexPair(up_list2[i], up_list2[i + 1] - 1));
                }
            }
            return output;
        }
    }
    public class PI2
    {
        public int _a, _b;
        public override string ToString()
        {
            string str = _a.ToString() + " ";
            str += _b.ToString();
            return str;
        }
        public int Filt(int a)
        {
            a = a % 360;
            if (a < 0) a = a + 360;
            return a;
        }
        public PI2(int a, int b)
        {
            _a = Filt(a);
            _b = Filt(b);
        }
        public List<string> ToExpression()
        {
            RegularExpression RE = new RegularExpression();
            List<string> output = new List<string>();
            List<IndexPair> pa = new List<IndexPair>();
            List<IndexPair> pb = new List<IndexPair>();
            int a = Filt(_a);
            int b = Filt(_b);
            List<int> list_int = new List<int>();
            list_int.Add(a); list_int.Add(Filt(b - 1));
            list_int.Add(b); list_int.Add(Filt(a - 1));

            list_int = list_int.Distinct<int>().ToList();
            if (list_int.Count != 4) return output;
            _a = a;
            _b = b;
            pa = RE.FromPIArray(new int[] { a, b - 1 });
            pb = RE.FromPIArray(new int[] { b, a - 1 });
            string str = "";
            for (int i = 0; i < pa.Count; i++)
            {
                if (i > 0) str += "|";
                str += RE.ExpressionPair(pa[i]);
            }
            output.Add(str); str = "";
            for (int i = 0; i < pb.Count; i++)
            {
                if (i > 0) str += "|";
                str += RE.ExpressionPair(pb[i]);
            }
            output.Add(str); str = "";
            return output;
        }
    }
    public class PI4
    {
        public int _a, _b, _c, _d;
        public override string ToString()
        {
            string str = _a.ToString() + " ";
            str += _b.ToString() + " ";
            str += _c.ToString() + " ";
            str += _d.ToString();
            return str;
        }
        public int Filt(int a)
        {
            a = a % 360;
            if (a < 0) a = a + 360;
            return a;
        }
        public PI4(int a, int b, int c, int d)
        {
            _a = Filt(a);
            _b = Filt(b);
            _c = Filt(c);
            _d = Filt(d);
        }
        public List<string> ToExpression()
        {
            RegularExpression RE = new RegularExpression();
            List<string> output = new List<string>();
            List<IndexPair> pa = new List<IndexPair>();
            List<IndexPair> pb = new List<IndexPair>();
            List<IndexPair> pc = new List<IndexPair>();
            List<IndexPair> pd = new List<IndexPair>();
            int a = Filt(_a);
            int b = Filt(_b);
            int c = Filt(_c);
            int d = Filt(_d);
            List<int> list_int = new List<int>();
            list_int.Add(a); list_int.Add(Filt(b - 1));
            list_int.Add(b); list_int.Add(Filt(c - 1));
            list_int.Add(c); list_int.Add(Filt(d - 1));
            list_int.Add(d); list_int.Add(Filt(a - 1));
            list_int = list_int.Distinct<int>().ToList();
            if (list_int.Count != 8) return output;
            _a = a;
            _b = b;
            _c = c;
            _d = d;
            pa = RE.FromPIArray(new int[] { a, b - 1 });
            pb = RE.FromPIArray(new int[] { b, c - 1 });
            pc = RE.FromPIArray(new int[] { c, d - 1 });
            pd = RE.FromPIArray(new int[] { d, a - 1 });
            string str = "";
            for (int i = 0; i < pa.Count; i++)
            {
                if (i > 0) str += "|";
                str += RE.ExpressionPair(pa[i]);
            }
            output.Add(str); str = "";
            for (int i = 0; i < pb.Count; i++)
            {
                if (i > 0) str += "|";
                str += RE.ExpressionPair(pb[i]);
            }
            output.Add(str); str = "";
            for (int i = 0; i < pc.Count; i++)
            {
                if (i > 0) str += "|";
                str += RE.ExpressionPair(pc[i]);
            }
            output.Add(str); str = "";
            for (int i = 0; i < pd.Count; i++)
            {
                if (i > 0) str += "|";
                str += RE.ExpressionPair(pd[i]);
            }
            output.Add(str); str = "";
            return output;
        }
    }
}
