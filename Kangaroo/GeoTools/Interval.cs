using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeoTools
{
    public struct Interval
    {
        #region Members
        private double m_t0;
        private double m_t1;
        public double T0 { get { return m_t0; } set { m_t0 = value; } }

        /// <summary>
        /// Gets or sets the upper bound of the Interval.
        /// </summary>
        public double T1 { get { return m_t1; } set { m_t1 = value; } }

        /// <summary>
        /// Gets or sets the indexed bound of this Interval.
        /// </summary>
        /// <param name="index">Bound index (0 = lower; 1 = upper).</param>
        public double this[int index]
        {
            get
            {
                if (0 == index) { return m_t0; }
                if (1 == index) { return m_t1; }

                // IronPython works with indexing is we thrown an IndexOutOfRangeException
                throw new IndexOutOfRangeException();
            }
            set
            {
                if (0 == index) { m_t0 = value; }
                else if (1 == index) { m_t1 = value; }
                else { throw new IndexOutOfRangeException(); }
            }
        }
        #endregion
        public Interval(double t0, double t1)
        {
            m_t0 = t0;
            m_t1 = t1;
        }
        public Interval(Interval other)
        {
            m_t0 = other.m_t0;
            m_t1 = other.m_t1;
        }
        #region Operators
        public static bool operator ==(Interval a, Interval b)
        {
            return a.CompareTo(b) == 0;
        }
        public static bool operator !=(Interval a, Interval b)
        {
            return a.CompareTo(b) != 0;
        }
        public static Interval operator +(Interval interval, double number)
        {
            return new Interval(interval.m_t0 + number, interval.m_t1 + number);
        }
        public static Interval operator +(double number, Interval interval)
        {
            return new Interval(interval.m_t0 + number, interval.m_t1 + number);
        }
        public static Interval operator -(Interval interval, double number)
        {
            return new Interval(interval.m_t0 - number, interval.m_t1 - number);
        }
        public static Interval operator -(double number, Interval interval)
        {
            return new Interval(number - interval.m_t0, number - interval.m_t1);
        }
        public static bool operator <(Interval a, Interval b)
        {
            return a.CompareTo(b) < 0;
        }
        public static bool operator <=(Interval a, Interval b)
        {
            return a.CompareTo(b) <= 0;
        }
        public static bool operator >(Interval a, Interval b)
        {
            return a.CompareTo(b) > 0;
        }
        public static bool operator >=(Interval a, Interval b)
        {
            return a.CompareTo(b) >= 0;
        }
        public override int GetHashCode()
        {
            // MSDN docs recommend XOR'ing the internal values to get a hash code
            return m_t0.GetHashCode() ^ m_t1.GetHashCode();
        }
        public override bool Equals(object obj)
        {
            return (obj is Interval && this == (Interval)obj);
        }
        #endregion
        public double Min
        {
            get { return (m_t0 <= m_t1 ? m_t0 : m_t1); }
        }
        public double Max
        {
            get { return (m_t0 <= m_t1 ? m_t1 : m_t0); }
        }
        public double Mid
        {
            get { return (m_t0 == m_t1 ? m_t0 : 0.5 * (m_t0 + m_t1)); }
        }
        public double Length
        {
            get { return m_t1 - m_t0; }
        }
        public bool Equals(Interval other)
        {
            return this == other;
        }
        public int CompareTo(Interval other)
        {
            if (m_t0 < other.m_t0)
                return -1;
            if (m_t0 > other.m_t0)
                return 1;
            if (m_t1 < other.m_t1)
                return -1;
            if (m_t1 > other.m_t1)
                return 1;
            return 0;
        }
    }
}
