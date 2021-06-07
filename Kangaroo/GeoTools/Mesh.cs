using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Kangaroo
{
    public class Mesh
    {
        List<Point3d> m_dV;
        List<Point3f> m_V;
        List<MeshFace> m_F;
        private MeshVertexList m_vertices;
        private MeshVertexNormalList m_normals;
        private MeshFaceList m_faces;
        private MeshFaceNormalList m_facenormals;
        public MeshVertexList Vertices
        {
            get { return m_vertices ?? (m_vertices = new MeshVertexList(this)); }
        }      
        public MeshFaceList Faces
        {
            get { return m_faces ?? (m_faces = new MeshFaceList(this)); }
        }
        public MeshVertexNormalList Normals
        {
            get { return m_normals ?? (m_normals = new MeshVertexNormalList(this)); }
        }
        public MeshFaceNormalList FaceNormals
        {
            get { return m_facenormals ?? (m_facenormals = new MeshFaceNormalList(this)); }
        }
    }
    public struct MeshFace
    {
        #region members
        internal int m_a;
        internal int m_b;
        internal int m_c;
        internal int m_d;
        #endregion
        #region constructors
        public MeshFace(int a, int b, int c)
        {
            m_a = a;
            m_b = b;
            m_c = c;
            m_d = c;
        }
        public MeshFace(int a, int b, int c, int d)
        {
            m_a = a;
            m_b = b;
            m_c = c;
            m_d = d;
        }
        public static MeshFace Unset
        {
            get { return new MeshFace(int.MinValue, int.MinValue, int.MinValue); }
        }
        #endregion
        #region properties
        internal string DebuggerDisplayUtil
        {
            get
            {
                if (IsTriangle)
                {
                    return string.Format(System.Globalization.CultureInfo.InvariantCulture, "T({0}, {1}, {2})", m_a, m_b, m_c);
                }
                return string.Format(System.Globalization.CultureInfo.InvariantCulture, "Q({0}, {1}, {2}, {3})", m_a, m_b, m_c, m_d);
            }
        }
        public int A
        {
            get { return m_a; }
            set { m_a = value; }
        }
        public int B
        {
            get { return m_b; }
            set { m_b = value; }
        }
        public int C
        {
            get { return m_c; }
            set { m_c = value; }
        }
        public int D
        {
            get { return m_d; }
            set { m_d = value; }
        }
        public int this[int index]
        {
            get
            {
                if (index == 0) return m_a;
                if (index == 1) return m_b;
                if (index == 2) return m_c;
                if (index == 3) return m_d;
                throw new IndexOutOfRangeException();
            }
            set
            {
                if (index == 0) m_a = value;
                else if (index == 1) m_b = value;
                else if (index == 2) m_c = value;
                else if (index == 3) m_d = value;
                else
                    throw new IndexOutOfRangeException();
            }
        }
        public bool IsValid()
        {
            if (m_a < 0) { return false; }
            if (m_b < 0) { return false; }
            if (m_c < 0) { return false; }
            if (m_d < 0) { return false; }
            if (m_a == m_b) { return false; }
            if (m_a == m_c) { return false; }
            if (m_a == m_d) { return false; }
            if (m_b == m_c) { return false; }
            if (m_b == m_d) { return false; }
            return true;
        }
        public bool IsValid(int vertexCount)
        {
            if (!IsValid()) { return false; }

            if (m_a >= vertexCount) { return false; }
            if (m_b >= vertexCount) { return false; }
            if (m_c >= vertexCount) { return false; }
            if (m_d >= vertexCount) { return false; }
            return true;
        }
        public bool IsTriangle { get { return m_c == m_d; } }
        public bool IsQuad { get { return m_c != m_d; } }
        #endregion
        #region methods
        public void Set(int a, int b, int c)
        {
            m_a = a;
            m_b = b;
            m_c = c;
            m_d = c;
        }
        public void Set(int a, int b, int c, int d)
        {
            m_a = a;
            m_b = b;
            m_c = c;
            m_d = d;
        }
        public MeshFace Flip()
        {
            if (m_c == m_d)
                return new MeshFace(m_a, m_c, m_b, m_b);
            return new MeshFace(m_a, m_d, m_c, m_a);
        }
        #endregion
    }
    public class MeshVertexList : List<Point3f>
    {
        private readonly Mesh m_mesh;
        internal MeshVertexList(Mesh ownerMesh)
        {
            m_mesh = ownerMesh;
        } 
        public void Add(float x, float y, float z)
        {
           Add(new Point3f(x, y, z));   
        }
        public void Add(double x, double y, double z)
        {
            Add((float)x, (float)y, (float)z);
        }
        public void Add(Point3d vertex)
        {
            Add((float)vertex.X, (float)vertex.Y, (float)vertex.Z);
        }    
    }
    public class MeshFaceList : List<MeshFace>
    {
        private readonly Mesh m_mesh;
        internal MeshFaceList(Mesh ownerMesh)
        {
            m_mesh = ownerMesh;
        }
    }
    public class MeshFaceNormalList : List<Vector3f>
    {
        private readonly Mesh m_mesh;    
        internal MeshFaceNormalList(Mesh ownerMesh)
        {
            m_mesh = ownerMesh;
        }
    }
    public class MeshVertexNormalList : List<Vector3f>
    {
        private readonly Mesh m_mesh;
        internal MeshVertexNormalList(Mesh ownerMesh)
        {
            m_mesh = ownerMesh;
        }
     
    }
    }

