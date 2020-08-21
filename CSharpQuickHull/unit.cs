using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace qhull
{
    public class point3d : vector3d
    {

        public point3d()
        {
        }

        public point3d(vector3d v)
        {
            set(v);
        }

        public point3d(double x, double y, double z)
        {
            set(x, y, z);
        }
    }


    public class vector3d
    {
        private const double DOUBLE_PREC = 2.2204460492503131e-16;
        public double x;
        public double y;
        public double z;
        public vector3d()
        {
        }
        public vector3d(vector3d v)
        {
            set(v);
        }
        public vector3d(double x, double y, double z)
        {
            set(x, y, z);
        }
        public double get(int i)
        {
            switch (i)
            {
                case 0:
                    {
                        return x;
                    }
                case 1:
                    {
                        return y;
                    }
                case 2:
                    {
                        return z;
                    }
                default:
                    {
                        throw new Exception("ArrayIndexOutOfBounds");
                    }
            }
        }
        public void set(int i, double value)
        {
            switch (i)
            {
                case 0:
                    {
                        x = value;
                        break;
                    }
                case 1:
                    {
                        y = value;
                        break;
                    }
                case 2:
                    {
                        z = value;
                        break;
                    }
                default:
                    {
                        throw new Exception("ArrayIndexOutOfBounds");
                    }
            }
        }
        public void set(vector3d v1)
        {
            x = v1.x;
            y = v1.y;
            z = v1.z;
        }
        public void add(vector3d v1, vector3d v2)
        {
            x = v1.x + v2.x;
            y = v1.y + v2.y;
            z = v1.z + v2.z;
        }
        public void add(vector3d v1)
        {
            x += v1.x;
            y += v1.y;
            z += v1.z;
        }
        public void sub(vector3d v1, vector3d v2)
        {
            x = v1.x - v2.x;
            y = v1.y - v2.y;
            z = v1.z - v2.z;
        }
        public void sub(vector3d v1)
        {
            x -= v1.x;
            y -= v1.y;
            z -= v1.z;
        }
        public void scale(double s)
        {
            x = s * x;
            y = s * y;
            z = s * z;
        }
        public void scale(double s, vector3d v1)
        {
            x = s * v1.x;
            y = s * v1.y;
            z = s * v1.z;
        }
        public double norm()
        {
            return Math.Sqrt(x * x + y * y + z * z);
        }
        public double normSquared()
        {
            return x * x + y * y + z * z;
        }
        public double distance(vector3d v)
        {
            double dx = x - v.x;
            double dy = y - v.y;
            double dz = z - v.z;

            return Math.Sqrt(dx * dx + dy * dy + dz * dz);
        }
        public double distanceSquared(vector3d v)
        {
            double dx = x - v.x;
            double dy = y - v.y;
            double dz = z - v.z;

            return (dx * dx + dy * dy + dz * dz);
        }
        public double dot(vector3d v1)
        {
            return x * v1.x + y * v1.y + z * v1.z;
        }
        public void normalize()
        {
            double lenSqr = x * x + y * y + z * z;
            double err = lenSqr - 1;
            if (err > (2 * DOUBLE_PREC) ||
                err < -(2 * DOUBLE_PREC))
            {
                double len = Math.Sqrt(lenSqr);
                x /= len;
                y /= len;
                z /= len;
            }
        }
        public void setZero()
        {
            x = 0;
            y = 0;
            z = 0;
        }
        public void set(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }
        public void cross(vector3d v1, vector3d v2)
        {
            double tmpx = v1.y * v2.z - v1.z * v2.y;
            double tmpy = v1.z * v2.x - v1.x * v2.z;
            double tmpz = v1.x * v2.y - v1.y * v2.x;

            x = tmpx;
            y = tmpy;
            z = tmpz;
        }
        protected void setRandom(double lower, double upper, Random generator)
        {
            double range = upper - lower;

            x = generator.NextDouble() * range + lower;
            y = generator.NextDouble() * range + lower;
            z = generator.NextDouble() * range + lower;
        }
        public String toString()
        {
            return x + " " + y + " " + z;
        }
    }


}
