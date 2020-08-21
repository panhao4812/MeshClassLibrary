using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace qhull
{
    public class QuickHull3D
    {
        public void Print(String str) { }
        private const double DOUBLE_PREC = 2.2204460492503131e-16;
        public readonly double AUTOMATIC_TOLERANCE = -1;
        private double explicitTolerance = -1;//AUTOMATIC_TOLERANCE;
        public readonly int CLOCKWISE = 0x1;
        public readonly int INDEXED_FROM_ONE = 0x2;
        public readonly int INDEXED_FROM_ZERO = 0x4;
        public readonly int POINT_RELATIVE = 0x8;
        private const int NONCONVEX_WRT_LARGER_FACE = 1;
        private const int NONCONVEX = 2;
        private int findIndex = -1;
        private bool debug = false;
        private Vertex[] pointBuffer = new Vertex[0];
        private int[] vertexPointIndices = new int[0];
        private Face[] discardedFaces = new Face[3];
        private Vertex[] maxVtxs = new Vertex[3];
        private Vertex[] minVtxs = new Vertex[3];
        private List<Face> faces = new List<Face>(16);
        private List<HalfEdge> horizon = new List<HalfEdge>(16);
        private FaceList newFaces = new FaceList();
        private VertexList unclaimed = new VertexList();
        private VertexList claimed = new VertexList();
        private int numVertices;
        private int numFaces;
        private int numPoints;
        private double tolerance;
        private double charLength;
        #region build
        public QuickHull3D()
        {
        }
        public QuickHull3D(double[] coords)
        {
            build(coords, coords.Length / 3);
        }
        public QuickHull3D(point3d[] points)
        {
            build(points, points.Length);
        }
        public void build(double[] coords)
        {
            build(coords, coords.Length / 3);
        }
        public void build(double[] coords, int nump)
        {
            if (nump < 4)
            {
                throw new Exception("Less than four input points specified");
            }
            if (coords.Length / 3 < nump)
            {
                throw new Exception("Coordinate array too small for specified number of points");
            }
            initBuffers(nump);
            setPoints(coords, nump);
            buildHull();
        }
        public void build(point3d[] points)
        {
            build(points, points.Length);
        }
        public void build(point3d[] points, int nump)
        {
            if (nump < 4)
            {
                throw new Exception("Less than four input points specified");
            }
            if (points.Length < nump)
            {
                throw new Exception("Point array too small for specified number of points");
            }
            initBuffers(nump);
            setPoints(points, nump);
            buildHull();
        }
        #endregion
        private void addPointToFace(Vertex vtx, Face face)
        {
            vtx.face = face;

            if (face.outside == null)
            {
                claimed.add(vtx);
            }
            else
            {
                claimed.insertBefore(vtx, face.outside);
            }
            face.outside = vtx;
        }
        private void removePointFromFace(Vertex vtx, Face face)
        {
            if (vtx == face.outside)
            {
                if (vtx.next != null && vtx.next.face == face)
                {
                    face.outside = vtx.next;
                }
                else
                {
                    face.outside = null;
                }
            }
            claimed.delete(vtx);
        }
        private Vertex removeAllPointsFromFace(Face face)
        {
            if (face.outside != null)
            {
                Vertex end = face.outside;
                while (end.next != null && end.next.face == face)
                {
                    end = end.next;
                }
                claimed.delete(face.outside, end);
                end.next = null;
                return face.outside;
            }
            else
            {
                return null;
            }
        }
        private HalfEdge findHalfEdge(Vertex tail, Vertex head)
        {
            for (int i = 0; i < faces.Count; i++)
            {
                HalfEdge he = faces[i].findEdge(tail, head);
                if (he != null)
                {
                    return he;
                }
            }
            return null;
        }
        private void setHull(double[] coords, int nump, int[][] faceIndices, int numf)
        {
            initBuffers(nump);
            setPoints(coords, nump);
            computeMaxAndMin();
            for (int i = 0; i < numf; i++)
            {
                Face face = Face.create(pointBuffer, faceIndices[i]);
                HalfEdge he = face.he0;
                do
                {
                    HalfEdge heOpp = findHalfEdge(he.head(), he.tail());
                    if (heOpp != null)
                    {
                        he.setOpposite(heOpp);
                    }
                    he = he.next;
                }
                while (he != face.he0);
                faces.Add(face);
            }
        }
        public void triangulate()
        {
            double minArea = 1000 * charLength * DOUBLE_PREC;
            newFaces.clear();
            for (int i = 0; i < faces.Count; i++)
            {
                if (faces[i].mark == Face.VISIBLE)
                {
                    faces[i].triangulate(newFaces, minArea);
                }
            }
            for (Face face = newFaces.first(); face != null; face = face.next)
            {
                faces.Add(face);
            }
        }
        private void initBuffers(int nump)
        {
            if (pointBuffer.Length < nump)
            {
                Vertex[] newBuffer = new Vertex[nump];
                vertexPointIndices = new int[nump];
                for (int i = 0; i < pointBuffer.Length; i++)
                {
                    newBuffer[i] = pointBuffer[i];
                }
                for (int i = pointBuffer.Length; i < nump; i++)
                {
                    newBuffer[i] = new Vertex();
                }
                pointBuffer = newBuffer;
            }
            faces.Clear();
            claimed.clear();
            numFaces = 0;
            numPoints = nump;
        }
        private void setPoints(double[] coords, int nump)
        {
            for (int i = 0; i < nump; i++)
            {
                Vertex vtx = pointBuffer[i];
                vtx.pnt.set(coords[i * 3 + 0], coords[i * 3 + 1], coords[i * 3 + 2]);
                vtx.index = i;
            }
        }
        private void setPoints(point3d[] pnts, int nump)
        {
            for (int i = 0; i < nump; i++)
            {
                Vertex vtx = pointBuffer[i];
                vtx.pnt.set(pnts[i]);
                vtx.index = i;
            }
        }
        private void computeMaxAndMin()
        {
            vector3d max = new vector3d();
            vector3d min = new vector3d();

            for (int i = 0; i < 3; i++)
            {
                maxVtxs[i] = minVtxs[i] = pointBuffer[0];
            }
            max.set(pointBuffer[0].pnt);
            min.set(pointBuffer[0].pnt);

            for (int i = 1; i < numPoints; i++)
            {
                point3d pnt = pointBuffer[i].pnt;
                if (pnt.x > max.x)
                {
                    max.x = pnt.x;
                    maxVtxs[0] = pointBuffer[i];
                }
                else if (pnt.x < min.x)
                {
                    min.x = pnt.x;
                    minVtxs[0] = pointBuffer[i];
                }
                if (pnt.y > max.y)
                {
                    max.y = pnt.y;
                    maxVtxs[1] = pointBuffer[i];
                }
                else if (pnt.y < min.y)
                {
                    min.y = pnt.y;
                    minVtxs[1] = pointBuffer[i];
                }
                if (pnt.z > max.z)
                {
                    max.z = pnt.z;
                    maxVtxs[2] = pointBuffer[i];
                }
                else if (pnt.z < min.z)
                {
                    min.z = pnt.z;
                    minVtxs[2] = pointBuffer[i];
                }
            }

            // this epsilon formula comes from QuickHull, and I'm
            // not about to quibble.
            charLength = Math.Max(max.x - min.x, max.y - min.y);
            charLength = Math.Max(max.z - min.z, charLength);
            if (explicitTolerance == AUTOMATIC_TOLERANCE)
            {
                tolerance =
                    3 * DOUBLE_PREC * (Math.Max(Math.Abs(max.x), Math.Abs(min.x)) +
                      Math.Max(Math.Abs(max.y), Math.Abs(min.y)) +
                      Math.Max(Math.Abs(max.z), Math.Abs(min.z)));
            }
            else
            {
                tolerance = explicitTolerance;
            }
        }
        private void createInitialSimplex()
        {
            double max = 0;
            int imax = 0;

            for (int i = 0; i < 3; i++)
            {
                double diff = maxVtxs[i].pnt.get(i) - minVtxs[i].pnt.get(i);
                if (diff > max)
                {
                    max = diff;
                    imax = i;
                }
            }

            if (max <= tolerance) { throw new Exception("Input points appear to be coincident"); }
            Vertex[] vtx = new Vertex[4];
            vtx[0] = maxVtxs[imax];
            vtx[1] = minVtxs[imax];
            vector3d u01 = new vector3d();
            vector3d diff02 = new vector3d();
            vector3d nrml = new vector3d();
            vector3d xprod = new vector3d();
            double maxSqr = 0;
            u01.sub(vtx[1].pnt, vtx[0].pnt);
            u01.normalize();
            for (int i = 0; i < numPoints; i++)
            {
                diff02.sub(pointBuffer[i].pnt, vtx[0].pnt);
                xprod.cross(u01, diff02);
                double lenSqr = xprod.normSquared();
                if (lenSqr > maxSqr &&
                pointBuffer[i] != vtx[0] &&  // paranoid
                pointBuffer[i] != vtx[1])
                {
                    maxSqr = lenSqr;
                    vtx[2] = pointBuffer[i];
                    nrml.set(xprod);
                }
            }
            if (Math.Sqrt(maxSqr) <= 100 * tolerance)
            { throw new Exception("Input points appear to be colinear"); }
            nrml.normalize();


            double maxDist = 0;
            double d0 = vtx[2].pnt.dot(nrml);
            for (int i = 0; i < numPoints; i++)
            {
                double dist = Math.Abs(pointBuffer[i].pnt.dot(nrml) - d0);
                if (dist > maxDist &&
                pointBuffer[i] != vtx[0] &&  // paranoid
                pointBuffer[i] != vtx[1] &&
                pointBuffer[i] != vtx[2])
                {
                    maxDist = dist;
                    vtx[3] = pointBuffer[i];
                }
            }
            if (Math.Abs(maxDist) <= 100 * tolerance)
            { throw new Exception("Input points appear to be coplanar"); }

            if (debug)
            {
                Print("initial vertices:");
                Print(vtx[0].index + ": " + vtx[0].pnt);
                Print(vtx[1].index + ": " + vtx[1].pnt);
                Print(vtx[2].index + ": " + vtx[2].pnt);
                Print(vtx[3].index + ": " + vtx[3].pnt);
            }

            Face[] tris = new Face[4];

            if (vtx[3].pnt.dot(nrml) - d0 < 0)
            {
                tris[0] = Face.createTriangle(vtx[0], vtx[1], vtx[2]);
                tris[1] = Face.createTriangle(vtx[3], vtx[1], vtx[0]);
                tris[2] = Face.createTriangle(vtx[3], vtx[2], vtx[1]);
                tris[3] = Face.createTriangle(vtx[3], vtx[0], vtx[2]);

                for (int i = 0; i < 3; i++)
                {
                    int k = (i + 1) % 3;
                    tris[i + 1].getEdge(1).setOpposite(tris[k + 1].getEdge(0));
                    tris[i + 1].getEdge(2).setOpposite(tris[0].getEdge(k));
                }
            }
            else
            {
                tris[0] = Face.createTriangle(vtx[0], vtx[2], vtx[1]);
                tris[1] = Face.createTriangle(vtx[3], vtx[0], vtx[1]);
                tris[2] = Face.createTriangle(vtx[3], vtx[1], vtx[2]);
                tris[3] = Face.createTriangle(vtx[3], vtx[2], vtx[0]);

                for (int i = 0; i < 3; i++)
                {
                    int k = (i + 1) % 3;
                    tris[i + 1].getEdge(0).setOpposite(tris[k + 1].getEdge(1));
                    tris[i + 1].getEdge(2).setOpposite(tris[0].getEdge((3 - i) % 3));
                }
            }


            for (int i = 0; i < 4; i++)
            {
                faces.Add(tris[i]);
            }

            for (int i = 0; i < numPoints; i++)
            {
                Vertex v = pointBuffer[i];

                if (v == vtx[0] || v == vtx[1] || v == vtx[2] || v == vtx[3])
                {
                    continue;
                }

                maxDist = tolerance;
                Face maxFace = null;
                for (int k = 0; k < 4; k++)
                {
                    double dist = tris[k].distanceToPlane(v.pnt);
                    if (dist > maxDist)
                    {
                        maxFace = tris[k];
                        maxDist = dist;
                    }
                }
                if (maxFace != null)
                {
                    addPointToFace(v, maxFace);
                }
            }
        }
        public int getNumVertices()
        {
            return numVertices;
        }
        public point3d[] getVertices()
        {
            point3d[] vtxs = new point3d[numVertices];
            for (int i = 0; i < numVertices; i++)
            {
                vtxs[i] = pointBuffer[vertexPointIndices[i]].pnt;
            }
            return vtxs;
        }
        public int getVertices(double[] coords)
        {
            for (int i = 0; i < numVertices; i++)
            {
                point3d pnt = pointBuffer[vertexPointIndices[i]].pnt;
                coords[i * 3 + 0] = pnt.x;
                coords[i * 3 + 1] = pnt.y;
                coords[i * 3 + 2] = pnt.z;
            }
            return numVertices;
        }
        public int[] getVertexPointIndices()
        {
            int[] indices = new int[numVertices];
            for (int i = 0; i < numVertices; i++)
            {
                indices[i] = vertexPointIndices[i];
            }
            return indices;
        }
        public int getNumFaces()
        {
            return faces.Count();
        }
        public int[][] getFaces()
        {
            return getFaces(0);
        }
        public int[][] getFaces(int indexFlags)
        {
            int[][] allFaces = new int[faces.Count][];
            for (int i = 0; i < faces.Count; i++)
            {
                allFaces[i] = new int[faces[i].numVertices()];
                getFaceIndices(allFaces[i], faces[i], indexFlags);
            }
            return allFaces;
        }
        private void getFaceIndices(int[] indices, Face face, int flags)
        {
            bool ccw = ((flags & CLOCKWISE) == 0);
            bool indexedFromOne = ((flags & INDEXED_FROM_ONE) != 0);
            bool pointRelative = ((flags & POINT_RELATIVE) != 0);

            HalfEdge hedge = face.he0;
            int k = 0;
            do
            {
                int idx = hedge.head().index;
                if (pointRelative)
                {
                    idx = vertexPointIndices[idx];
                }
                if (indexedFromOne)
                {
                    idx++;
                }
                indices[k++] = idx;
                hedge = (ccw ? hedge.next : hedge.prev);
            }
            while (hedge != face.he0);
        }
        private void resolveUnclaimedPoints(FaceList newFaces)
        {
            Vertex vtxNext = unclaimed.first();
            for (Vertex vtx = vtxNext; vtx != null; vtx = vtxNext)
            {
                vtxNext = vtx.next;

                double maxDist = tolerance;
                Face maxFace = null;
                for (Face newFace = newFaces.first(); newFace != null;
                 newFace = newFace.next)
                {
                    if (newFace.mark == Face.VISIBLE)
                    {
                        double dist = newFace.distanceToPlane(vtx.pnt);
                        if (dist > maxDist)
                        {
                            maxDist = dist;
                            maxFace = newFace;
                        }
                        if (maxDist > 1000 * tolerance)
                        {
                            break;
                        }
                    }
                }
                if (maxFace != null)
                {
                    addPointToFace(vtx, maxFace);
                    if (debug && vtx.index == findIndex)
                    {
                        Print(findIndex + " CLAIMED BY " +
                         maxFace.getVertexString());
                    }
                }
                else
                {
                    if (debug && vtx.index == findIndex)
                    {
                        Print(findIndex + " DISCARDED");
                    }
                }
            }
        }
        private void deleteFacePoints(Face face, Face absorbingFace)
        {
            Vertex faceVtxs = removeAllPointsFromFace(face);
            if (faceVtxs != null)
            {
                if (absorbingFace == null)
                {
                    unclaimed.addAll(faceVtxs);
                }
                else
                {
                    Vertex vtxNext = faceVtxs;
                    for (Vertex vtx = vtxNext; vtx != null; vtx = vtxNext)
                    {
                        vtxNext = vtx.next;
                        double dist = absorbingFace.distanceToPlane(vtx.pnt);
                        if (dist > tolerance)
                        {
                            addPointToFace(vtx, absorbingFace);
                        }
                        else
                        {
                            unclaimed.add(vtx);
                        }
                    }
                }
            }
        }
        private double oppFaceDistance(HalfEdge he)
        {
            return he.face.distanceToPlane(he.opposite.face.getCentroid());
        }
        private bool doAdjacentMerge(Face face, int mergeType)
        {
            HalfEdge hedge = face.he0;

            bool convex = true;
            do
            {
                Face oppFace = hedge.oppositeFace();
                bool merge = false;
                double dist1;

                if (mergeType == NONCONVEX)
                { // then merge faces if they are definitively non-convex
                    if (oppFaceDistance(hedge) > -tolerance ||
                        oppFaceDistance(hedge.opposite) > -tolerance)
                    {
                        merge = true;
                    }
                }
                else // mergeType == NONCONVEX_WRT_LARGER_FACE
                { // merge faces if they are parallel or non-convex
                    // wrt to the larger face; otherwise, just mark
                    // the face non-convex for the second pass.
                    if (face.area > oppFace.area)
                    {
                        if ((dist1 = oppFaceDistance(hedge)) > -tolerance)
                        {
                            merge = true;
                        }
                        else if (oppFaceDistance(hedge.opposite) > -tolerance)
                        {
                            convex = false;
                        }
                    }
                    else
                    {
                        if (oppFaceDistance(hedge.opposite) > -tolerance)
                        {
                            merge = true;
                        }
                        else if (oppFaceDistance(hedge) > -tolerance)
                        {
                            convex = false;
                        }
                    }
                }

                if (merge)
                {
                    if (debug)
                    {
                        Print(
                        "  merging " + face.getVertexString() + "  and  " +
                        oppFace.getVertexString());
                    }

                    int numd = face.mergeAdjacentFace(hedge, discardedFaces);
                    for (int i = 0; i < numd; i++)
                    {
                        deleteFacePoints(discardedFaces[i], face);
                    }
                    if (debug)
                    {
                        Print(
                           "  result: " + face.getVertexString());
                    }
                    return true;
                }
                hedge = hedge.next;
            }
            while (hedge != face.he0);
            if (!convex)
            {
                face.mark = Face.NON_CONVEX;
            }
            return false;
        }
        private void calculateHorizon(point3d eyePnt, HalfEdge edge0, Face face, List<HalfEdge> horizon)
        {
            //	   oldFaces.add (face);
            deleteFacePoints(face, null);
            face.mark = Face.DELETED;
            if (debug)
            {
                Print("  visiting face " + face.getVertexString());
            }
            HalfEdge edge;
            if (edge0 == null)
            {
                edge0 = face.getEdge(0);
                edge = edge0;
            }
            else
            {
                edge = edge0.getNext();
            }
            do
            {
                Face oppFace = edge.oppositeFace();
                if (oppFace.mark == Face.VISIBLE)
                {
                    if (oppFace.distanceToPlane(eyePnt) > tolerance)
                    {
                        calculateHorizon(eyePnt, edge.getOpposite(),
                                  oppFace, horizon);
                    }
                    else
                    {
                        horizon.Add(edge);
                        if (debug)
                        {
                            Print("  adding horizon edge " +
                                    edge.getVertexString());
                        }
                    }
                }
                edge = edge.getNext();
            }
            while (edge != edge0);
        }
        private HalfEdge addAdjoiningFace(Vertex eyeVtx, HalfEdge he)
        {
            Face face = Face.createTriangle(
               eyeVtx, he.tail(), he.head());
            faces.Add(face);
            face.getEdge(-1).setOpposite(he.getOpposite());
            return face.getEdge(0);
        }
        private void addNewFaces(FaceList newFaces, Vertex eyeVtx, List<HalfEdge> horizon)
        {
            newFaces.clear();
            HalfEdge hedgeSidePrev = null;
            HalfEdge hedgeSideBegin = null;
            for (int i = 0; i < horizon.Count; i++)
            {
                HalfEdge horizonHe = horizon[i];
                HalfEdge hedgeSide = addAdjoiningFace(eyeVtx, horizonHe);
                if (debug)
                {
                    Print("new face: " + hedgeSide.face.getVertexString());
                }
                if (hedgeSidePrev != null)
                {
                    hedgeSide.next.setOpposite(hedgeSidePrev);
                }
                else
                {
                    hedgeSideBegin = hedgeSide;
                }
                newFaces.add(hedgeSide.getFace());
                hedgeSidePrev = hedgeSide;
            }
            hedgeSideBegin.next.setOpposite(hedgeSidePrev);
        }
        private Vertex nextPointToAdd()
        {
            if (!claimed.isEmpty())
            {
                Face eyeFace = claimed.first().face;
                Vertex eyeVtx = null;
                double maxDist = 0;
                for (Vertex vtx = eyeFace.outside;
                 vtx != null && vtx.face == eyeFace;
                 vtx = vtx.next)
                {
                    double dist = eyeFace.distanceToPlane(vtx.pnt);
                    if (dist > maxDist)
                    {
                        maxDist = dist;
                        eyeVtx = vtx;
                    }
                }
                return eyeVtx;
            }
            else
            {
                return null;
            }
        }
        private void addPointToHull(Vertex eyeVtx)
        {
            horizon.Clear();
            unclaimed.clear();

            if (debug)
            {
                Print("Adding point: " + eyeVtx.index);
                Print(
                   " which is " + eyeVtx.face.distanceToPlane(eyeVtx.pnt) +
                   " above face " + eyeVtx.face.getVertexString());
            }
            removePointFromFace(eyeVtx, eyeVtx.face);
            calculateHorizon(eyeVtx.pnt, null, eyeVtx.face, horizon);
            newFaces.clear();
            addNewFaces(newFaces, eyeVtx, horizon);

            // first merge pass ... merge faces which are non-convex
            // as determined by the larger face

            for (Face face = newFaces.first(); face != null; face = face.next)
            {
                if (face.mark == Face.VISIBLE)
                {
                    while (doAdjacentMerge(face, NONCONVEX_WRT_LARGER_FACE))
                        ;
                }
            }
            // second merge pass ... merge faces which are non-convex
            // wrt either face	     
            for (Face face = newFaces.first(); face != null; face = face.next)
            {
                if (face.mark == Face.NON_CONVEX)
                {
                    face.mark = Face.VISIBLE;
                    while (doAdjacentMerge(face, NONCONVEX))
                        ;
                }
            }
            resolveUnclaimedPoints(newFaces);
        }
        private void buildHull()
        {
            int cnt = 0;
            Vertex eyeVtx;
            computeMaxAndMin();
            createInitialSimplex();
            while ((eyeVtx = nextPointToAdd()) != null)
            {
                addPointToHull(eyeVtx);
                cnt++;
                if (debug)
                {
                    Print("iteration " + cnt + " done");
                }
            }
            reindexFacesAndVertices();
            if (debug)
            {
                Print("hull done");
            }
        }
        private void markFaceVertices(Face face, int mark)
        {
            HalfEdge he0 = face.getFirstEdge();
            HalfEdge he = he0;
            do
            {
                he.head().index = mark;
                he = he.next;
            }
            while (he != he0);
        }
        private void reindexFacesAndVertices()
        {
            for (int i = 0; i < numPoints; i++) { pointBuffer[i].index = -1; }
            // remove inactive faces and mark active vertices
            numFaces = 0;
            List<Face> output = new List<Face>();
            for (int i = 0; i < faces.Count; i++)
            {
                Face face = faces[i];
                if (face.mark != Face.VISIBLE)
                {
                    // it.remove();
                }
                else
                {
                    output.Add(face);
                    markFaceVertices(face, 0);
                    numFaces++;
                }
            }
            faces = output;
            // reindex vertices
            numVertices = 0;
            for (int i = 0; i < numPoints; i++)
            {
                Vertex vtx = pointBuffer[i];
                if (vtx.index == 0)
                {
                    vertexPointIndices[numVertices] = i;
                    vtx.index = numVertices++;
                }
            }
        }
        /*	public bool getDebug()
         {
           return debug;
         }

        public void setDebug (bool enable)
         { 
           debug = enable;
         }

        public double getDistanceTolerance()
         {
           return tolerance;
         }

        public void setExplicitDistanceTolerance(double tol)
         { 
           explicitTolerance = tol;
         }

        public double getExplicitDistanceTolerance()
         {
           return explicitTolerance;
         }
        private void printQhullErrors (Process proc)
           throws IOException
         {
           bool wrote = false;
           InputStream es = proc.getErrorStream();
           while (es.available() > 0)
            { System.out.write (es.read());
              wrote = true;
            }
           if (wrote)
            { Print("");
            }
         }

        private void setFromQhull (double[] coords, int nump,
                         bool triangulate)
         {
           String commandStr = "./qhull i";
           if (triangulate)
            { commandStr += " -Qt"; 
            }
           try
            { 
              Process proc = Runtime.getRuntime().exec (commandStr);
              PrintStream ps = new PrintStream (proc.getOutputStream());
              StreamTokenizer stok =
             new StreamTokenizer (
                new InputStreamReader (proc.getInputStream()));

              ps.println ("3 " + nump);
              for (int i=0; i<nump; i++)
               { ps.println (
                coords[i*3+0] + " " +
                coords[i*3+1] + " " +  
                coords[i*3+2]);
               }
              ps.flush();
              ps.close();
              Vector indexList = new Vector(3);
              stok.eolIsSignificant(true);
              printQhullErrors (proc);
	      
              do
               { stok.nextToken();
               }
              while (stok.sval == null ||
                 !stok.sval.startsWith ("MERGEexact"));
              for (int i=0; i<4; i++)
               { stok.nextToken();
               }
              if (stok.ttype != StreamTokenizer.TT_NUMBER)
               { Print ("Expecting number of faces");
             System.exit(1); 
               }
              int numf = (int)stok.nval;
              stok.nextToken(); // clear EOL
              int[][] faceIndices = new int[numf][];
              for (int i=0; i<numf; i++)
               { indexList.clear();
             while (stok.nextToken() != StreamTokenizer.TT_EOL)
              { if (stok.ttype != StreamTokenizer.TT_NUMBER)
                 { Print ("Expecting face index");
                   System.exit(1); 
                 }
                indexList.add (0, new Integer((int)stok.nval));
              }
             faceIndices[i] = new int[indexList.size()];
             int k = 0;
             for (Iterator it=indexList.iterator(); it.hasNext(); ) 
              { faceIndices[i][k++] = ((Integer)it.next()).intValue();
              }
               }
              setHull (coords, nump, faceIndices, numf);
            }
           catch (Exception e) 
            { e.printStackTrace();
              System.exit(1); 
            }
         }

        private void printPoints ()
         {
           for (int i=0; i<numPoints; i++)
            { point3d pnt = pointBuffer[i].pnt;
              Print (pnt.x + ", " + pnt.y + ", " + pnt.z + ",");
            }
         }
       public void print (PrintStream ps)
         {
           print (ps, 0);
         }

        public void print (PrintStream ps, int indexFlags)
         {
           if ((indexFlags & INDEXED_FROM_ZERO) == 0)
            { indexFlags |= INDEXED_FROM_ONE;
            }
           for (int i=0; i<numVertices; i++)
            { point3d pnt = pointBuffer[vertexPointIndices[i]].pnt;
              ps.println ("v " + pnt.x + " " + pnt.y + " " + pnt.z);
            }
           for (Iterator fi=faces.iterator(); fi.hasNext(); )
            { Face face = (Face)fi.next();
              int[] indices = new int[face.numVertices()];
              getFaceIndices (indices, face, indexFlags);

              ps.print ("f");
              for (int k=0; k<indices.Length; k++)
               { ps.print (" " + indices[k]); 
               }
              ps.println ("");
            }
         }

        
         	private void splitFace (Face face)
        	 {
          	   Face newFace = face.split();
         	   if (newFace != null)
         	    { newFaces.add (newFace);
         	      splitFace (newFace);
         	      splitFace (face);
          	    }
         	 }	
        
        private bool checkFaceConvexity (
           Face face, double tol, PrintStream ps)
         {
           double dist;
           HalfEdge he = face.he0;
           do
            { face.checkConsistency();
              // make sure edge is convex
              dist = oppFaceDistance (he);
              if (dist > tol)
               { if (ps != null)
              { ps.println ("Edge " + he.getVertexString() +
                    " non-convex by " + dist);
              }
             return false;
               }
              dist = oppFaceDistance (he.opposite);
              if (dist > tol)
               { if (ps != null)
              { ps.println ("Opposite edge " +
                    he.opposite.getVertexString() +
                    " non-convex by " + dist);
              }
             return false;
               }
              if (he.next.oppositeFace() == he.oppositeFace())
               { if (ps != null)
              { ps.println ("Redundant vertex " + he.head().index +
                    " in face " + face.getVertexString());
              }
             return false;
               }
              he = he.next;
            }
           while (he != face.he0);	   
           return true;
         }
        private bool checkFaces(double tol, PrintStream ps)
         { 
           // check edge convexity
           bool convex = true;
           for (Iterator it=faces.iterator(); it.hasNext(); ) 
            { Face face = (Face)it.next();
              if (face.mark == Face.VISIBLE)
               { if (!checkFaceConvexity (face, tol, ps))
              { convex = false;
              }
               }
            }
           return convex;
         }
        public bool check (PrintStream ps)
         {
           return check (ps, getDistanceTolerance());
         }
        public bool check (PrintStream ps, double tol)

         {
           double dist;
           double pointTol = 10*tol;

           if (!checkFaces(tolerance, ps))
            { return false; 
            }

           for (int i=0; i<numPoints; i++)
            { point3d pnt = pointBuffer[i].pnt;
              for (Iterator it=faces.iterator(); it.hasNext(); ) 
               { Face face = (Face)it.next();
             if (face.mark == Face.VISIBLE)
              { 
                dist = face.distanceToPlane (pnt);
                if (dist > pointTol)
                 { if (ps != null)
                { ps.println (
                     "Point " + i + " " + dist + " above face " +
                     face.getVertexString());
                }
                   return false;
                 }
              }
               }
            }
           return true;
         }
         * */
    }

}
