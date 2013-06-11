#!/usr/bin/env python

# This example demonstrates how to use a programmable filter and how
# to use the special vtkDataSetToDataSet::GetOutputPort() methods

import vtk
from numpy import *

# We create a 100 by 100 point plane to sample
plane = vtk.vtkPlaneSource()
plane.SetXResolution(200)
plane.SetYResolution(200)

# We transform the plane by a factor of 10 on X and Y
transform = vtk.vtkTransform()
transform.Scale(10, 10, 1)
transF = vtk.vtkTransformPolyDataFilter()
transF.SetInputConnection(plane.GetOutputPort())
transF.SetTransform(transform)

# Compute Bessel function and derivatives. We'll use a programmable filter
# for this. Note the unusual GetPolyDataInput() & GetOutputPort() methods.
test = vtk.vtkProgrammableFilter()
test.SetInputConnection(transF.GetOutputPort())

# The SetExecuteMethod takes a Python function as an argument
# In here is where all the processing is done.
def damage():
    input = test.GetPolyDataInput()
    numPts = input.GetNumberOfPoints()
    newPts = vtk.vtkPoints()
    isolines = vtk.vtkFloatArray()

    for i in range(0, numPts):
        x = input.GetPoint(i)
        x0, x1 = x[:2]

        r = sqrt(x0**2+x1**2)
        if r < 1:
            x2 = 1
        else:
            x2 = (1/sqrt(r))*(abs(cos(x1/r)))**2 + (1/r**3)*abs((sin(x1/r)))**2
        
        if r < 1:
            isoline = 1
        else:
            isoline = (1/sqrt(r))*(abs(cos(x1/r)))**2 + (1/r**2)*abs((sin(x1/r)))**2

        newPts.InsertPoint(i, x0, x1, x2)
        isolines.InsertValue(i, isoline) 

    test.GetPolyDataOutput().CopyStructure(input)
    test.GetPolyDataOutput().SetPoints(newPts)
    test.GetPolyDataOutput().GetPointData().SetScalars(isolines)

test.SetExecuteMethod(damage) 

# We warp the plane based on the scalar values calculated above
warp = vtk.vtkWarpScalar()
warp.SetInput(test.GetPolyDataOutput())
warp.XYPlaneOn()
warp.SetScaleFactor(0.5)


# We create a mapper and actor as usual. In the case we adjust the
# scalar range of the mapper to match that of the computed scalars
mapper = vtk.vtkPolyDataMapper()
mapper.SetInput(warp.GetPolyDataOutput())
mapper.SetScalarRange(test.GetPolyDataOutput().GetScalarRange())
carpet = vtk.vtkActor()
carpet.SetMapper(mapper)

# Create the RenderWindow, Renderer
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

ren.AddActor(carpet)
renWin.SetSize(500, 500)

ren.ResetCamera()
ren.GetActiveCamera().Zoom(1.5)

iren.Initialize()
renWin.Render()
iren.Start()
