'#Language "WWB-COM"
Option Explicit

Function toStr(val as Double)
    toStr = Cstr(val)
    If Left(toStr, 1) = "." Then
        toStr = "0" + toStr
    End If
    If Left(toStr, 2) = "-." Then
        toStr = "-0" + Right(toStr, Len(toStr) - 1)
    End If
End Function

Sub PrintPorts(fn As Integer)
    Dim unit
    unit = Units.GetGeometryUnitToSI()
    Print #fn, """ports"": ["
    Dim x0 As Double
    Dim x1 As Double
    Dim y0 As Double
    Dim y1 As Double
    Dim z0 As Double
    Dim z1 As Double
    Dim portNum
    portNum = 0
    Dim portCount
    portCount = 0
    While portNum < 9
        portNum = portNum + 1
        On Error GoTo tryNext
        DiscretePort.GetLength(portNum)
        DiscretePort.GetCoordinates(portNum, x0, y0, z0, x1, y1, z1)
        If portCount > 0 Then
            Print #fn, ","
        End If
        Print #fn, "  {"
        Print #fn, "    ""src"": [" + toStr(x0 * unit) + ", " + toStr(y0 * unit) + ", " + toStr(z0 * unit) + "],"
        Print #fn, "    ""dst"": [" + toStr(x1 * unit) + ", " + toStr(y1 * unit) + ", " + toStr(z1 * unit) + "] "
        Print #fn, "  }"
        portCount = portCount + 1
        tryNext:
    Wend
    Print #fn, "],"
End Sub

Sub PrintSolids(fn As Integer)
    Dim solidNum, solidName, i, jsonPath
    solidNum = Solid.GetNumberOfShapes
    jsonPath = Environ("JSON_PATH")
    For i = 0 To solidNum - 1
        solidName = Solid.GetNameOfShapeFromIndex(i)
        Print #fn, "  {"
        Print #fn, "    ""index"": " + Cstr(i) + ","
        Print #fn, "    ""name"": """ + solidName + ""","
        If Environ("EXPORT_SOLIDS") <> "" Then
            Dim idx, comp, shape, file
            idx = InStrRev(solidName, ":")
            If idx >= 0 Then
                comp = Left(solidName, idx - 1)
                shape = Right(solidName, Len(solidName) - idx)
                file = jsonPath + "." + Cstr(i)
                With STL
                    .Reset
                    .ExportFileUnits("mm")
                    .SurfaceTolerance(0.001)
                    .NormalTolerance(30)
                    .FileName(file + ".stl")
                    .Component(comp)
                    .Name(shape)
                    .Write
                End With
                Print #fn, "    ""stl"": """ + Replace(file + ".stl", "\", "\\") + ""","
            End If
        End If
        Print #fn, "    ""material"": """ + Solid.GetMaterialNameForShape(solidName) + """"
        If i < solidNum - 1 Then
            Print #fn, "  },"
        Else
            Print #fn, "  }"
        End If
    Next
End Sub

Sub PrintGrid(fn As Integer)
    Dim n, i
    n = Mesh.GetNX()
    Print #fn, "  ""xs"": ["
    For i = 0 To n - 1
        Print #fn, "    " + toStr(Mesh.GetXPos(i))
        If i < n - 1 Then
            Print #fn, ","
        End If
    Next
    Print #fn, "  ],"
    n = Mesh.GetNY()
    Print #fn, "  ""ys"": ["
    For i = 0 To n - 1
        Print #fn, "    " + toStr(Mesh.GetYPos(i))
        If i < n - 1 Then
            Print #fn, ","
        End If
    Next
    Print #fn, "  ],"
    n = Mesh.GetNZ()
    Print #fn, "  ""zs"": ["
    For i = 0 To n - 1
        Print #fn, "    " + toStr(Mesh.GetZPos(i))
        If i < n - 1 Then
            Print #fn, ","
        End If
    Next
    Print #fn, "  ]"
End Sub

Sub Main
    OpenFile(Environ("CST_PATH"))
    Dim fn
    fn = FreeFile()
    Open Environ("JSON_PATH") for output as #fn
    Print #fn, "{"
    Dim fqUnit, geoUnit, timeUnit
    fqUnit = Units.GetFrequencyUnitToSI()
    geoUnit = Units.GetGeometryUnitToSI()
    timeUnit = Units.GetTimeUnitToSI()
    Print #fn, """freq"": [" + toStr(Solver.GetFmin() * fqUnit) + "," + toStr(Solver.GetFmax() * fqUnit) + "],"
    Print #fn, """units"": {"
    Print #fn, "  ""geometry"": " + toStr(geoUnit) + ","
    Print #fn, "  ""time"": " + toStr(timeUnit) + ","
    Print #fn, "  ""frequency"": " + toStr(fqUnit)
    Print #fn, "},"
    PrintPorts(fn)
    Print #fn, """solids"": ["
    PrintSolids(fn)
    Print #fn, "],"
    Print #fn, """grid"": {"
    PrintGrid(fn)
    Print #fn, "},"
    Print #fn, """version"": 2019"
    Print #fn, "}"
    Close #fn
    If Environ("BUILD_MATRIX") <> "" Then
        With Solver
            .SteadyStateDurationType "Time"
            .SteadyStateDurationTime "0.0001"
            .Start
        End With
    End If
	SelectTreeItem("Excitation Signals")
	With ASCIIExport
		.FileName(Environ("JSON_PATH") + ".excitation.txt")
		.Execute()
	End With
End Sub
