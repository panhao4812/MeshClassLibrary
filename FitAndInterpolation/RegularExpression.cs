﻿using System;
using System.Collections.Generic;
using System.Linq;

namespace FitAndInterpolation
{
    public struct IndexDomain
    {
        public int I;
        public int J;
        public IndexDomain(int i, int j)
        {
            I = i;
            J = j;
        }
    }
    public class RegularExpression
    {
        public int Filt(int a)
        {
            a = a % 360;
            if (a < 0) a = a + 360;
            return a;
        }
        private string str;
        private void Print(string text)
        {
            str += text + "\r\n";
        }
        public void PrintPairList(List<IndexDomain> pair)
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
        public List<IndexDomain> CombinePair(List<IndexDomain> PairList)
        {
            List<int> list_int = new List<int>();
            for (int i = 0; i < PairList.Count; i++)
            {
                list_int.Add(PairList[i].I);
                list_int.Add(PairList[i].J);
            }
            list_int = list_int.Distinct<int>().ToList();

            List<IndexDomain> ShortPairs = new List<IndexDomain>();
            for (int i = 0; i < list_int.Count - 1; i++)
            {
                if (list_int[i + 1] == list_int[i] + 1) ShortPairs.Add(new IndexDomain(list_int[i], list_int[i + 1]));
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
        private bool EdgeTest(List<IndexDomain> pairs, int input)
        {
            double a = input - 0.1; double b = input + 0.1;
            bool signa = false; bool signb = false;
            foreach (IndexDomain pair in pairs)
            {
                double I = (double)pair.I;
                double J = (double)pair.J;
                if (a <= J && a >= I) signa = true;
                if (b <= J && b >= I) signb = true;
                if (signa && signb) break;
            }
            return (signa && signb);
        }
        public int CompairPair(IndexDomain x, IndexDomain y)
        {
            if (x.I < y.I) return -1;
            if (x.I > y.I) return 1;
            return 0;
        }
        public List<IndexDomain> FromIntArray(List<int> input)
        {
            List<IndexDomain> output = new List<IndexDomain>();
            input = input.Distinct<int>().ToList();
            for (int i = 0; i < input.Count - 1; i++)
            {
                output.Add(new IndexDomain(input[i], input[i + 1]));
            }
            return output;
        }
        public List<IndexDomain> FromPairArray(int[] input)
        {
            List<IndexDomain> output = new List<IndexDomain>();
            for (int i = 0; i < input.Length / 2.0; i++)
            {
                int a1 = input[i * 2], a2 = input[i * 2 + 1];
                if (a1 > a2)
                {
                    output.Add(new IndexDomain(a2, a1));
                }
                else if (a2 > a1)
                {
                    output.Add(new IndexDomain(a1, a2));
                }
            }
            return output;
        }
        public List<IndexDomain> FromPIArray(int[] input, int index)
        {
            List<IndexDomain> output = new List<IndexDomain>();
            for (int i = 0; i < input.Length / 2.0; i++)
            {
                int a1 = input[i * 2], a2 = input[i * 2 + 1];
                if (a1 > a2)
                {
                    if (a1 < 359)
                    {
                        output.Add(new IndexDomain(a1 + index * 360, 359 + index * 360));
                    }
                    if (a2 > 0)
                    {
                        output.Add(new IndexDomain(0 + index * 360, a2 + index * 360));
                    }
                }
                else if (a2 > a1)
                {
                    output.Add(new IndexDomain(a1 + index * 360, a2 + index * 360));
                }
            }
            return output;
        }
        public List<IndexDomain> FromPairArray(List<int> input)
        {
            return FromPairArray(input.ToArray());
        }
        //必须是I<J
        public string ExpressionPair(IndexDomain pair)
        {
            List<IndexDomain> pairs = SplitPair(pair);
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
        public string ExpressionShortPair(IndexDomain pair)
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
            if (Char_b[Char_b.Length - 1] == '0')
            {
                output = pair.J.ToString() + "|";
                Char_b = (pair.J - 1).ToString().ToCharArray();
            }
            if (Char_a.Length != Char_b.Length) { return ""; }
            for (int i = 0; i < Char_a.Length; i++)
            {
                if (Char_a[i] == Char_b[i]) { output += Char_a[i]; }
                else { output += "[" + Char_a[i] + "-" + Char_b[i] + "]"; }
            }
            return output;
        }
        public List<IndexDomain> SplitPair(IndexDomain pair)
        {
            List<IndexDomain> output = new List<IndexDomain>();
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
                if (i == 0 && Char_a[Char_a.Length - 1] == '0' && Char_a.Length != 1) continue;
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
                if (i == 0) { output.Add(new IndexDomain(up_list2[i], up_list2[i + 1])); }
                else
                {
                    output.Add(new IndexDomain(up_list2[i] + 1, up_list2[i + 1]));
                }
            }
            //  if (output.Count == 1) return output;
            IndexDomain lastpair = output[output.Count - 1];
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
            up_list2 = up_list2.Distinct().ToList();
            if (up_list2.Count == 1) { output.Add(new IndexDomain(lastpair.I, up_list2[0])); return output; }
            if (up_list2[0] > lastpair.I) output.Add(new IndexDomain(lastpair.I, up_list2[0] - 1));

            for (int i = 0; i < up_list2.Count - 1; i++)
            {
                if (i == up_list2.Count - 2) { output.Add(new IndexDomain(up_list2[i], up_list2[i + 1])); }
                else
                {
                    output.Add(new IndexDomain(up_list2[i], up_list2[i + 1] - 1));
                }
            }
            return output;
        }
        public string SolvePI2(IndexDomain pos1, IndexDomain pos2)
        {
            string str = "";
            int int_a = Filt(pos1.I);
            int int_b = Filt(pos1.J);
            str += int_a.ToString() + " " + int_b.ToString() + "\r\n";
            PI2 pi1 = new PI2(int_a, int_b, 1);
            int_a = Filt(pos2.J);
            int_b = Filt(pos2.I);
            str += int_a.ToString() + " " + int_b.ToString() + "\r\n";
            PI2 pi2 = new PI2(int_a, int_b, 2);
            List<string> output1 = pi1.ToExpression();
            List<string> output2 = pi2.ToExpression();
            str += "\r\n";
            for (int i = 0; i < 2; i++)
            {
                str += "^(" + output1[i] + "|" + output2[i] + ")$";
                str += "\r\n"; str += "\r\n";
            }
            return str;
        }
        public string SolvePI3(IndexDomain pos1, IndexDomain pos2, IndexDomain pos3)
        {
            string str = "";
            int Length1 = Convert.ToInt32(180.0 - Math.Abs(pos1.J - pos1.I) / 2.0);
            int int_a = Filt(pos1.I);
            int int_b = Filt(pos1.J);
            int int_c = Filt(pos1.J + Length1);
            str += int_a.ToString() + " " + int_b.ToString() + " " + int_c.ToString() + "\r\n";
            PI3 pi1 = new PI3(int_a, int_b, int_c, 1);
            Length1 = Convert.ToInt32(180.0 - Math.Abs(pos2.J - pos2.I) / 2.0);
            int_a = Filt(pos2.I - Length1);
            int_b = Filt(pos2.I);
            int_c = Filt(pos2.J);
            str += int_a.ToString() + " " + int_b.ToString() + " " + int_c.ToString() + "\r\n";
            PI3 pi2 = new PI3(int_a, int_b, int_c, 2);
            Length1 = Convert.ToInt32(180.0 - Math.Abs(pos3.J - pos3.I) / 2.0);
            int_a = Filt(pos3.J);
            int_b = Filt(pos3.J + Length1);
            int_c = Filt(pos3.I);
            str += int_a.ToString() + " " + int_b.ToString() + " " + int_c.ToString() + "\r\n";
            PI3 pi3 = new PI3(int_a, int_b, int_c, 3);
            List<string> output1 = pi1.ToExpression();
            List<string> output2 = pi2.ToExpression();
            List<string> output3 = pi3.ToExpression();
            str += "\r\n";
            for (int i = 0; i < 3; i++)
            {
                str += "^(" + output1[i] + "|" + output2[i] + "|"
                      + output3[i] + ")$";
                str += "\r\n"; str += "\r\n";
            }
            return str;
        }
        public string SolvePI4(IndexDomain pos1, IndexDomain pos2, IndexDomain pos3, IndexDomain pos4)
        {
            string str = "";
            int Length1 = Convert.ToInt32(120.0 - Math.Abs(pos1.J - pos1.I) / 3.0);
            int int_a = Filt(pos1.I);
            int int_b = Filt(pos1.J);
            int int_c = Filt(pos1.J + Length1);
            int int_d = Filt(pos1.I - Length1);
            str += int_a.ToString() + " " + int_b.ToString() + " " +
                int_c.ToString() + " " + int_d.ToString() + "\r\n";
            PI4 pi1 = new PI4(int_a, int_b, int_c, int_d, 1);
            Length1 = Convert.ToInt32(120.0 - Math.Abs(pos2.J - pos2.I) / 3.0);
            int_a = Filt(pos2.I - Length1);
            int_b = Filt(pos2.I);
            int_c = Filt(pos2.J);
            int_d = Filt(pos2.J + Length1);
            str += int_a.ToString() + " " + int_b.ToString() + " " +
               int_c.ToString() + " " + int_d.ToString() + "\r\n";
            PI4 pi2 = new PI4(int_a, int_b, int_c, int_d, 2);
            Length1 = Convert.ToInt32(120.0 - Math.Abs(pos3.J - pos3.I) / 3.0);
            int_a = Filt(pos3.J + Length1);
            int_b = Filt(pos3.I - Length1);
            int_c = Filt(pos3.I);
            int_d = Filt(pos3.J);
            str += int_a.ToString() + " " + int_b.ToString() + " " +
               int_c.ToString() + " " + int_d.ToString() + "\r\n";
            PI4 pi3 = new PI4(int_a, int_b, int_c, int_d, 3);
            Length1 = Convert.ToInt32(120.0 - Math.Abs(pos4.J - pos4.I) / 3.0);
            int_a = Filt(pos4.J);
            int_b = Filt(pos4.J + Length1);
            int_c = Filt(pos4.I - Length1);
            int_d = Filt(pos4.I);
            str += int_a.ToString() + " " + int_b.ToString() + " " +
               int_c.ToString() + " " + int_d.ToString() + "\r\n";
            PI4 pi4 = new PI4(int_a, int_b, int_c, int_d, 4);
            List<string> output1 = pi1.ToExpression();
            List<string> output2 = pi2.ToExpression();
            List<string> output3 = pi3.ToExpression();
            List<string> output4 = pi4.ToExpression();
            str += "\r\n";
            for (int i = 0; i < 4; i++)
            {
                str += "^(" + output1[i] + "|" + output2[i] + "|"
                      + output3[i] + "|" + output4[i] + ")$";
                str += "\r\n"; str += "\r\n";
            }
            return str;
        }
        public string SolvePI8(IndexDomain pos1, IndexDomain pos2, IndexDomain pos3, IndexDomain pos4,
            IndexDomain pos5, IndexDomain pos6, IndexDomain pos7, IndexDomain pos8)
        {
            string str = "";

            int Length1 = Convert.ToInt32(120.0 - Math.Abs(pos1.J - pos1.I) / 3.0);
            int int_a = Filt(pos1.I);
            int int_b = Filt(pos1.J);
            int int_c = Filt(pos1.J + Length1);
            int int_d = Filt(pos1.I - Length1);
            str += int_a.ToString() + " " + int_b.ToString() + " " +
                int_c.ToString() + " " + int_d.ToString() + "\r\n";
            PI8 pi1 = new PI8(int_a, int_b, int_c, int_d, 1);

            Length1 = Convert.ToInt32(120.0 - Math.Abs(pos2.J - pos2.I) / 3.0);
            int_a = Filt(pos2.I - Length1);
            int_b = Filt(pos2.I);
            int_c = Filt(pos2.J);
            int_d = Filt(pos2.J + Length1);
            str += int_a.ToString() + " " + int_b.ToString() + " " +
               int_c.ToString() + " " + int_d.ToString() + "\r\n";
            PI8 pi2 = new PI8(int_a, int_b, int_c, int_d, 2);

            Length1 = Convert.ToInt32(120.0 - Math.Abs(pos3.J - pos3.I) / 3.0);
            int_a = Filt(pos3.J + Length1);
            int_b = Filt(pos3.I - Length1);
            int_c = Filt(pos3.I);
            int_d = Filt(pos3.J);
            str += int_a.ToString() + " " + int_b.ToString() + " " +
               int_c.ToString() + " " + int_d.ToString() + "\r\n";
            PI8 pi3 = new PI8(int_a, int_b, int_c, int_d, 3);

            Length1 = Convert.ToInt32(120.0 - Math.Abs(pos4.J - pos4.I) / 3.0);
            int_a = Filt(pos4.J);
            int_b = Filt(pos4.J + Length1);
            int_c = Filt(pos4.I - Length1);
            int_d = Filt(pos4.I);
            str += int_a.ToString() + " " + int_b.ToString() + " " +
               int_c.ToString() + " " + int_d.ToString() + "\r\n";
            PI8 pi4 = new PI8(int_a, int_b, int_c, int_d, 4);

            Length1 = Convert.ToInt32(120.0 - Math.Abs(pos5.J - pos5.I) / 3.0);
            int_a = Filt(pos5.I);
            int_b = Filt(pos5.J);
            int_c = Filt(pos5.J + Length1);
            int_d = Filt(pos5.I - Length1);
            str += int_a.ToString() + " " + int_b.ToString() + " " +
               int_c.ToString() + " " + int_d.ToString() + "\r\n";
            PI8 pi5 = new PI8(int_a, int_b, int_c, int_d, 5);

            Length1 = Convert.ToInt32(120.0 - Math.Abs(pos6.J - pos6.I) / 3.0);
            int_a = Filt(pos6.I - Length1);
            int_b = Filt(pos6.I);
            int_c = Filt(pos6.J);
            int_d = Filt(pos6.J + Length1);
            str += int_a.ToString() + " " + int_b.ToString() + " " +
               int_c.ToString() + " " + int_d.ToString() + "\r\n";
            PI8 pi6 = new PI8(int_a, int_b, int_c, int_d, 6);

            Length1 = Convert.ToInt32(120.0 - Math.Abs(pos7.J - pos7.I) / 3.0);
            int_a = Filt(pos7.J + Length1);
            int_b = Filt(pos7.I - Length1);
            int_c = Filt(pos7.I);
            int_d = Filt(pos7.J);
            str += int_a.ToString() + " " + int_b.ToString() + " " +
               int_c.ToString() + " " + int_d.ToString() + "\r\n";
            PI8 pi7 = new PI8(int_a, int_b, int_c, int_d, 7);

            Length1 = Convert.ToInt32(120.0 - Math.Abs(pos8.J - pos8.I) / 3.0);
            int_a = Filt(pos8.J);
            int_b = Filt(pos8.J + Length1);
            int_c = Filt(pos8.I - Length1);
            int_d = Filt(pos8.I);
            str += int_a.ToString() + " " + int_b.ToString() + " " +
               int_c.ToString() + " " + int_d.ToString() + "\r\n";
            PI8 pi8 = new PI8(int_a, int_b, int_c, int_d, 8);

            List<string> output1 = pi1.ToExpression();
            List<string> output2 = pi2.ToExpression();
            List<string> output3 = pi3.ToExpression();
            List<string> output4 = pi4.ToExpression();
            List<string> output5 = pi5.ToExpression();
            List<string> output6 = pi6.ToExpression();
            List<string> output7 = pi7.ToExpression();
            List<string> output8 = pi8.ToExpression();
            str += "\r\n";
            str += "^(" + output5[0] + "|" + output2[1] + "|" + output3[2] + "|" + output4[3] + ")$";
            str += "\r\n"; str += "\r\n";
            str += "^(" + output1[0] + "|" + output6[1] + "|" + output3[2] + "|" + output4[3] + ")$";
            str += "\r\n"; str += "\r\n";
            str += "^(" + output1[0] + "|" + output2[1] + "|" + output7[2] + "|" + output4[3] + ")$";
            str += "\r\n"; str += "\r\n";
            str += "^(" + output1[0] + "|" + output2[1] + "|" + output3[2] + "|" + output8[3] + ")$";
            str += "\r\n"; str += "\r\n";
            str += "^(" + output1[0] + "|" + output6[1] + "|" + output7[2] + "|" + output8[3] + ")$";
            str += "\r\n"; str += "\r\n";
            str += "^(" + output5[0] + "|" + output2[1] + "|" + output7[2] + "|" + output8[3] + ")$";
            str += "\r\n"; str += "\r\n";
            str += "^(" + output5[0] + "|" + output6[1] + "|" + output3[2] + "|" + output8[3] + ")$";
            str += "\r\n"; str += "\r\n";
            str += "^(" + output5[0] + "|" + output6[1] + "|" + output7[2] + "|" + output4[3] + ")$";
            str += "\r\n"; str += "\r\n";
            return str;
        }
    }
    public class PI2
    {
        RegularExpression RE = new RegularExpression();
        public int _a, _b, _index;
        public override string ToString()
        {
            string str = _a.ToString() + " ";
            str += _b.ToString();
            return str;
        }
        public PI2(int a, int b, int index)
        {
            _a = RE.Filt(a);
            _b = RE.Filt(b);
            _index = index;
        }
        public List<string> ToExpression()
        {

            List<string> output = new List<string>();
            List<IndexDomain> pa = new List<IndexDomain>();
            List<IndexDomain> pb = new List<IndexDomain>();
            int a = RE.Filt(_a);
            int b = RE.Filt(_b);
            List<int> list_int = new List<int>();
            list_int.Add(a); list_int.Add(RE.Filt(b - 1));
            list_int.Add(b); list_int.Add(RE.Filt(a - 1));
            list_int = list_int.Distinct<int>().ToList();
            if (list_int.Count != 4) return output;
            pa = RE.FromPIArray(new int[] { a, RE.Filt(b - 1) }, _index);
            pb = RE.FromPIArray(new int[] { b, RE.Filt(a - 1) }, _index);
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
        public int _a, _b, _c, _d, _index;
        public override string ToString()
        {
            string str = _a.ToString() + " ";
            str += _b.ToString() + " ";
            str += _c.ToString() + " ";
            str += _d.ToString();
            return str;
        }
        RegularExpression RE = new RegularExpression();
        public PI4(int a, int b, int c, int d, int index)
        {
            _a = RE.Filt(a);
            _b = RE.Filt(b);
            _c = RE.Filt(c);
            _d = RE.Filt(d);
            _index = index;
        }
        public List<string> ToExpression()
        {

            List<string> output = new List<string>();
            List<IndexDomain> pa = new List<IndexDomain>();
            List<IndexDomain> pb = new List<IndexDomain>();
            List<IndexDomain> pc = new List<IndexDomain>();
            List<IndexDomain> pd = new List<IndexDomain>();
            int a = RE.Filt(_a);
            int b = RE.Filt(_b);
            int c = RE.Filt(_c);
            int d = RE.Filt(_d);
            List<int> list_int = new List<int>();
            list_int.Add(a); list_int.Add(RE.Filt(b - 1));
            list_int.Add(b); list_int.Add(RE.Filt(c - 1));
            list_int.Add(c); list_int.Add(RE.Filt(d - 1));
            list_int.Add(d); list_int.Add(RE.Filt(a - 1));
            list_int = list_int.Distinct<int>().ToList();
            if (list_int.Count != 8) return output;
            pa = RE.FromPIArray(new int[] { a, RE.Filt(b - 1) }, _index);
            pb = RE.FromPIArray(new int[] { b, RE.Filt(c - 1) }, _index);
            pc = RE.FromPIArray(new int[] { c, RE.Filt(d - 1) }, _index);
            pd = RE.FromPIArray(new int[] { d, RE.Filt(a - 1) }, _index);
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
    public class PI3
    {
        public int _a, _b, _c, _index;
        RegularExpression RE = new RegularExpression();
        public override string ToString()
        {
            string str = _a.ToString() + " ";
            str += _b.ToString() + " " + _c.ToString();
            return str;
        }
        public PI3(int a, int b, int c, int index)
        {
            _a = RE.Filt(a);
            _b = RE.Filt(b);
            _c = RE.Filt(c);
            _index = index;
        }
        public List<string> ToExpression()
        {

            List<string> output = new List<string>();
            List<IndexDomain> pa = new List<IndexDomain>();
            List<IndexDomain> pb = new List<IndexDomain>();
            List<IndexDomain> pc = new List<IndexDomain>();
            int a = RE.Filt(_a);
            int b = RE.Filt(_b);
            int c = RE.Filt(_c);
            List<int> list_int = new List<int>();
            list_int.Add(a); list_int.Add(RE.Filt(b - 1));
            list_int.Add(b); list_int.Add(RE.Filt(c - 1));
            list_int.Add(c); list_int.Add(RE.Filt(a - 1));
            list_int = list_int.Distinct<int>().ToList();
            if (list_int.Count != 6) return output;
            pa = RE.FromPIArray(new int[] { a, RE.Filt(b - 1) }, _index);
            pb = RE.FromPIArray(new int[] { b, RE.Filt(c - 1) }, _index);
            pc = RE.FromPIArray(new int[] { c, RE.Filt(a - 1) }, _index);
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
            return output;
        }
    }
    public class PI8
    {
        public int _a, _b, _c, _d, _index;
        public override string ToString()
        {
            string str = _a.ToString() + " ";
            str += _b.ToString() + " ";
            str += _c.ToString() + " ";
            str += _d.ToString();
            return str;
        }
        RegularExpression RE = new RegularExpression();
        public PI8(int a, int b, int c, int d, int index)
        {
            _a = RE.Filt(a);
            _b = RE.Filt(b);
            _c = RE.Filt(c);
            _d = RE.Filt(d);
            _index = index;
        }
        public List<string> ToExpression()
        {
            List<string> output = new List<string>();
            List<IndexDomain> pa = new List<IndexDomain>();
            List<IndexDomain> pb = new List<IndexDomain>();
            List<IndexDomain> pc = new List<IndexDomain>();
            List<IndexDomain> pd = new List<IndexDomain>();
            int a = RE.Filt(_a);
            int b = RE.Filt(_b);
            int c = RE.Filt(_c);
            int d = RE.Filt(_d);
            List<int> list_int = new List<int>();
            list_int.Add(a); list_int.Add(RE.Filt(b - 1));
            list_int.Add(b); list_int.Add(RE.Filt(c - 1));
            list_int.Add(c); list_int.Add(RE.Filt(d - 1));
            list_int.Add(d); list_int.Add(RE.Filt(a - 1));
            list_int = list_int.Distinct<int>().ToList();
            if (list_int.Count != 8) return output;
            int _index2 = 0;
            if (_index >= 5)
            {
                _index2 = _index - 4;
            }
            else
            {
                _index2 = _index + 4;
            }
            pa = RE.FromPIArray(new int[] { a, RE.Filt(b - 1) }, _index);
            pb = RE.FromPIArray(new int[] { b, RE.Filt(c - 1) }, _index);
            pc = RE.FromPIArray(new int[] { c, RE.Filt(d - 1) }, _index);
            pd = RE.FromPIArray(new int[] { d, RE.Filt(a - 1) }, _index);
            if (_index == 1 || _index == 5) pa = RE.FromPIArray(new int[] { a, RE.Filt(b - 1) }, _index2);
            if (_index == 2 || _index == 6) pb = RE.FromPIArray(new int[] { b, RE.Filt(c - 1) }, _index2);
            if (_index == 3 || _index == 7) pc = RE.FromPIArray(new int[] { c, RE.Filt(d - 1) }, _index2);
            if (_index == 4 || _index == 8) pd = RE.FromPIArray(new int[] { d, RE.Filt(a - 1) }, _index2);
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
