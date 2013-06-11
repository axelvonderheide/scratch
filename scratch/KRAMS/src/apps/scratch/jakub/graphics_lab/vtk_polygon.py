import vtk

polygonPoints = vtk.vtkPoints()
polygonPoints.SetNumberOfPoints(4)
polygonPoints.InsertPoint(0, 0, 0, 0)
polygonPoints.InsertPoint(1, 1, 0, 0)
polygonPoints.InsertPoint(2, 1, 1, 0)
polygonPoints.InsertPoint(3, 0, 1, 0)
aPolygon = vtk.vtkPolygon()
aPolygon.GetPointIds().SetNumberOfIds(4)
aPolygon.GetPointIds().SetId(0, 0)
aPolygon.GetPointIds().SetId(1, 1)
aPolygon.GetPointIds().SetId(2, 2)
aPolygon.GetPointIds().SetId(3, 3)
aPolygonGrid = vtk.vtkUnstructuredGrid()
aPolygonGrid.Allocate(1, 1)
aPolygonGrid.InsertNextCell(aPolygon.GetCellType(), aPolygon.GetPointIds())
aPolygonGrid.SetPoints(polygonPoints)
aPolygonMapper = vtk.vtkDataSetMapper()
aPolygonMapper.SetInput(aPolygonGrid)
aPolygonActor = vtk.vtkActor()
aPolygonActor.SetMapper(aPolygonMapper)
aPolygonActor.AddPosition(6, 0, 2)
aPolygonActor.GetProperty().SetDiffuseColor(1, .4, .5)

# Create the usual rendering stuff.
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(300, 150)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

ren.SetBackground(.1, .2, .4)

ren.AddActor(aPolygonActor)

ren.ResetCamera()
ren.GetActiveCamera().Azimuth(30)
ren.GetActiveCamera().Elevation(20)
ren.GetActiveCamera().Dolly(2.8)
ren.ResetCameraClippingRange()

# Render the scene and start interaction.
iren.Initialize()
renWin.Render()
iren.Start()