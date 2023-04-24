import pygmsh
import meshio
import matplotlib.pyplot as plt

# Initialize the geometry object
geom = pygmsh.Geometry()

# Define the points
p1 = geom.add_point([0, 0, 0], 0.1)
p2 = geom.add_point([1, 0, 0], 0.1)
p3 = geom.add_point([1, 1, 0], 0.1)
p4 = geom.add_point([0, 1, 0], 0.1)

# Define the lines
l1 = geom.add_line(p1, p2)
l2 = geom.add_line(p2, p3)
l3 = geom.add_line(p3, p4)
l4 = geom.add_line(p4, p1)

# Define the plane surface
plate = geom.add_polygon([
    p1,
    p2,
    p3,
    p4,
], lcar=0.1)

# Generate the mesh
mesh = pygmsh.generate_mesh(geom)

# Plot the mesh
plt.triplot(mesh.points[:, 0], mesh.points[:, 1], mesh.cells["triangle"])
plt.show()
