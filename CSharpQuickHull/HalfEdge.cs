using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace qhull
{

   public class HalfEdge
    {
        public Vertex vertex;
        public Face face;
        public HalfEdge next;
        public HalfEdge prev;
        public HalfEdge opposite;
        public HalfEdge(Vertex v, Face f)
        {
            vertex = v;
            face = f;
        }
        public HalfEdge()
        {
        }
        public void setNext(HalfEdge edge)
        {
            next = edge;
        }
        public HalfEdge getNext()
        {
            return next;
        }
        public void setPrev(HalfEdge edge)
        {
            prev = edge;
        }
        public HalfEdge getPrev()
        {
            return prev;
        }
        public Face getFace()
        {
            return face;
        }
        public HalfEdge getOpposite()
        {
            return opposite;
        }
        public void setOpposite(HalfEdge edge)
        {
            opposite = edge;
            edge.opposite = this;
        }
        public Vertex head()
        {
            return vertex;
        }
        public Vertex tail()
        {
            return prev != null ? prev.vertex : null;
        }
        public Face oppositeFace()
        {
            return opposite != null ? opposite.face : null;
        }
        public String getVertexString()
        {
            if (tail() != null)
            {
                return "" +
                    tail().index + "-" +
                    head().index;
            }
            else
            {
                return "?-" + head().index;
            }
        }
        public double length()
        {
            if (tail() != null)
            {
                return head().pnt.distance(tail().pnt);
            }
            else
            {
                return -1;
            }
        }
        public double lengthSquared()
        {
            if (tail() != null)
            {
                return head().pnt.distanceSquared(tail().pnt);
            }
            else
            {
                return -1;
            }
        }

    }
}
