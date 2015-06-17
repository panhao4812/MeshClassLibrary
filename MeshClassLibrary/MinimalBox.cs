using Rhino;
using Rhino.Geometry;
using Rhino.DocObjects;
using Rhino.Collections;

using GH_IO;
using GH_IO.Serialization;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

using System;
using System.IO;
using System.Xml;
using System.Xml.Linq;
using System.Linq;
using System.Data;
using System.Drawing;
using System.Reflection;
using System.Collections;
using System.Windows.Forms;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using Rhino.Geometry.Collections;
namespace MeshClassLibrary
{
    class MinimalBox
    {
public MinimalBox(){}
        public Polyline MinimalBox2D(List<Point3d> x){
         Grasshopper.Kernel.Geometry.Node2List list = new Grasshopper.Kernel.Geometry.Node2List(x);
    Polyline pl = Grasshopper.Kernel.Geometry.ConvexHull.Solver.ComputeHull(list);
    // List<Polyline> boxes = new List<Polyline>();
    Polyline output = new Polyline();
    double t = double.MaxValue;
    for(int i = 0;i < pl.Count - 1;i++){
      Vector3d Xaxis = pl[i + 1] - pl[i];
      Vector3d Yaxis = Vector3d.CrossProduct(Xaxis, Vector3d.ZAxis);
      Plane p = new Plane(pl[i], Xaxis, Yaxis);
      Polyline pl2 = new Polyline(pl);
      pl2.Transform(Transform.PlaneToPlane(p, Plane.WorldXY));
      Rhino.Geometry.BoundingBox box = pl2.BoundingBox;
      Polyline pl3 = new Polyline();
      pl3.Add(box.Corner(false, false, false));
      pl3.Add(box.Corner(false, true, false));
      pl3.Add(box.Corner(true, true, false));
      pl3.Add(box.Corner(true, false, false));
      pl3.Add(box.Corner(false, false, false));
      double area = pl3[1].DistanceTo(pl3[0]) * pl3[1].DistanceTo(pl3[2]);
      if(area < t){t = area;  pl3.Transform(Transform.PlaneToPlane(Plane.WorldXY, p));output = pl3;}
      // boxes.Add(pl3);
    }
    return output;
    }
}

    /*
     * 3D box:
      Private Sub RunScript(ByVal x As String, ByVal y As List(Of Point3d), ByRef A As Object, ByRef B As Object) 
    Dim outa As New List (Of String)
    Dim ls As New List(Of Polyline)()
    Try
      Dim qproc As New System.Diagnostics.ProcessStartInfo("D:\ph的桌面\Mesh逆向\qhull-2012.1\bin\qconvex.exe")
      qproc.UseShellExecute = False
      qproc.RedirectStandardOutput = True
      qproc.RedirectStandardError = True
      qproc.RedirectStandardInput = True
      qproc.CreateNoWindow = True
      qproc.Arguments = "i"
      Dim res As New System.Diagnostics.Process()
      res = System.Diagnostics.Process.Start(qproc)

      Dim inp As System.IO.StreamWriter = res.StandardInput
      inp.WriteLine("3 grasshopper qhull call")
      inp.WriteLine(y.Count)
      For i As int32=0 To y.count - 1
        Dim fp As New Point3d(y(i))
        inp.WriteLine(fp.x.ToString + " " + fp.y.ToString + " " + fp.z.ToString)
      Next
      inp.Close()

      Dim  stdout As System.IO.StreamReader = res.StandardOutput
      Dim  stderr As System.IO.StreamReader = res.StandardError

      'Dim outa As String = stdout.ReadToEnd()

      Dim str As String = "open"
      Dim count As Int32 = 0


      Do While ((Not str = "") And count < 10000)
        Dim pl1 As New Polyline()
        str = stdout.ReadLine()
        ' Print(str)
        If str = Nothing Then Exit Do
        Dim str2 As String() = str.Split(" ")
        outa.Add(str)
        If str2.Length > 2 Then
          For j As int32=0 To str2.Length - 2
            Dim Num  As int32 = Convert.ToInt32(str2(j))
            pl1.Add(y(Num))
          Next
          pl1.Add(pl1(0))
          ls.Add(pl1)
        End If
        count += 1
      Loop
      Print("--stderr start--")
      Print(stderr.ReadToEnd())
      Print("--stderr end--")
      Dim aexit As Boolean = res.WaitForExit(5000)
      If (aExit) Then
        Print("--qhull calculations complete--")
        Print("--qhull exits--")
      Else
        Print("qhull still on. why?")
      End If


      Dim output As New Mesh()
      Dim plane1 As  Plane
      Dim  t As Double = Double.MaxValue
      For i As int32=0 To ls.Count - 1
        plane1 = New Plane(ls(i)(0), ls(i)(1), ls(i)(2))
        Dim whole As New Polyline(y)
        whole.Transform(Transform.PlaneToPlane(plane1, Plane.WorldXY))
        Dim box As Rhino.Geometry.BoundingBox = whole.BoundingBox
        Dim p0 As Point3d = box.Corner(False, False, False)
        Dim p1 As Point3d = box.Corner(False, True, False)
        Dim p2 As Point3d = box.Corner(True, False, False)
        Dim p3 As Point3d = box.Corner(False, False, True)
        Dim area As Double = p0.DistanceTo(p2) * p0.DistanceTo(p3) * p0.DistanceTo(p1)
        If area < t Then
          Print(t)
          t = area
          Dim mesh1 As New Mesh()
          Dim p4 As Point3d = box.Corner(True, True, True)
          Dim p5 As Point3d = box.Corner(False, True, True)
          Dim p6 As Point3d = box.Corner(True, False, True)
          Dim p7 As Point3d = box.Corner(True, True, False)
          mesh1.Vertices.Add(p0)
          mesh1.Vertices.Add(p1)
          mesh1.Vertices.Add(p2)
          mesh1.Vertices.Add(p3)
          mesh1.Vertices.Add(p4)
          mesh1.Vertices.Add(p5)
          mesh1.Vertices.Add(p6)
          mesh1.Vertices.Add(p7)
          mesh1.Faces.AddFace(0, 2, 7, 1)
          mesh1.Faces.AddFace(6, 3, 5, 4)
          mesh1.Faces.AddFace(6, 2, 7, 4)
          mesh1.Faces.AddFace(0, 1, 5, 3)
          mesh1.Faces.AddFace(0, 2, 6, 3)
          mesh1.Faces.AddFace(4, 5, 1, 7)
          mesh1.Transform(Transform.PlaneToPlane(Plane.WorldXY, plane1))
          output = mesh1

        End If
      Next

      A = output
    Catch ex As exception
      print(ex.ToString())

    End Try

  End Sub 
   *////
}
