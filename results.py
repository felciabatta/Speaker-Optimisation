import numpy as np
import vtk 


filename = 'geomeshes/experimental/enclosed speaker.vtk'



reader = vtk.vtkGenericDataObjectReader()
reader.SetFileName(filename)
reader.Update()

points = np.array( reader.GetOutput().GetPoints().GetData() )

