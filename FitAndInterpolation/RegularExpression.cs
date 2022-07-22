using Rhino;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FitAndInterpolation
{
    public class RegularExpression
    {
        public List<IndexPair> FromPairArray(int[] input)
        {
            List<IndexPair> output = new List<IndexPair>();
            for (int i = 0; i < input.Length / 2.0; i++)
            {
                int a1 = input[i], a2 = input[i + 1];
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
        public List<IndexPair> FromPairArray(List<int> input)
        {
            return FromPairArray(input.ToArray());
        }
        //必须是J>I
        public string ExpressionPair(IndexPair pair)
        {
            char[] Char_a = pair.I.ToString().ToCharArray();
            char[] Char_b = pair.J.ToString().ToCharArray();
            int[] int_a = new int[Char_a.Length];
            int[] int_b = new int[Char_b.Length];
            for (int i = 0; i < Char_a.Length; i++)
            {
                int_a[i] = Convert.ToInt32(Char_a[i]);
            }
            for (int i = 0; i < Char_b.Length; i++)
            {
                int_b[i] = Convert.ToInt32(Char_b[i]);
            }

            if (int_a.Length != int_b.Length) { return ""; }
            string output = "";
            for (int i = 0; i < int_a.Length; i++)
            {
                if (int_a[i] == int_b[i]) { output += Char_a[i]; }
                else { output += "[" + Char_a[i] + "-" + Char_b[i] + "]"; }
            }
            return output;
        }
        public List<IndexPair> SplitPair(IndexPair pair)
        {
            List<IndexPair> output = new List<IndexPair>();
            char[] Char_a = pair.I.ToString().ToCharArray();
            char[] Char_b = pair.J.ToString().ToCharArray();
            int[] int_a = new int[Char_a.Length];
            int[] int_b = new int[Char_b.Length];
            for (int i = 0; i < Char_a.Length; i++)
            {
                int_a[i] = Convert.ToInt32(Char_a[i]);
            }
            for (int i = 0; i < Char_b.Length; i++)
            {
                int_b[i] = Convert.ToInt32(Char_b[i]);
            }
            List<int> up_list = new List<int>();
            up_list.Add(pair.I);
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
            for (int i = 0; i < up_list2.Count - 1; i++)
            {
                if (i == 0) { output.Add(new IndexPair(up_list2[i], up_list2[i + 1])); }
                else
                {
                    output.Add(new IndexPair(up_list2[i] + 1, up_list2[i + 1]));
                }
            }
            if (output.Count == 1) return output;
            IndexPair lastpair = output[output.Count - 1];
            output.RemoveAt(output.Count - 1);
            up_list.Clear();
            up_list.Add(pair.J);
            for (int i = 0; i < Char_b.Length; i++)
            {
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
}
