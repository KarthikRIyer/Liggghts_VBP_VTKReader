import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from pysph.base.utils import get_particle_array
from pysph.tools.interpolator import Interpolator

filename = '/home/karthik/Desktop/BTP-Vibrating-Packed-Bed/post-spherical/post/dump28500000.vtk'

reader = vtk.vtkGenericDataObjectReader()
reader.SetFileName(filename)
reader.Update()
point_obj = reader.GetOutput()
points = point_obj.GetPoints()
points_data = point_obj.GetPointData()
v = points_data.GetArray('v')
m = points_data.GetArray('mass')
p = points.GetData()
pts_array = vtk_to_numpy(p)
vel_array = vtk_to_numpy(v)
mass_array = vtk_to_numpy(m)
masses = np.unique(mass_array)
masses = np.sort(masses)

x_array, y_array, z_array = pts_array.T
vx_array, vy_array, vz_array = vel_array.T
# print(vx)

x_small = []
y_small = []
z_small = []
m_small = []
m_large = []
x_large = []
y_large = []
z_large = []
vx_small = []
vy_small = []
vz_small = []
vx_large = []
vy_large = []
vz_large = []
rho_small = []
rho_large = []

for index in range(len(mass_array)):
    mass = mass_array[index].item()
    vx = vx_array[index].item()
    vy = vy_array[index].item()
    vz = vz_array[index].item()
    x = x_array[index].item()
    y = y_array[index].item()
    z = z_array[index].item()
    if mass == masses[0]:
        x_small.append(x)
        y_small.append(y)
        z_small.append(z)
        vx_small.append(vx)
        vy_small.append(vy)
        vz_small.append(vz)
        m_small.append(mass)
        rho_small.append(2500)
    elif mass == masses[1]:
        x_large.append(x)
        y_large.append(y)
        z_large.append(z)
        vx_large.append(vx)
        vy_large.append(vy)
        vz_large.append(vz)
        m_large.append(mass)
        rho_large.append(2500)
    else:
        print("mass doesn't exist")

p_small = get_particle_array(name='p_small', x=x_small, y=y_small, z=z_small, u=vx_small,
                             v=vy_small,
                             w=vz_small, m=m_small, rho=rho_small)
# p_small.set_output_arrays(['u', 'v', 'w', 'm'])
# print(p_small.get_property_arrays())
# p_small
interp = Interpolator([p_small], num_points=len(mass_array)*2)
# z_grad = interp.interpolate(prop='w', gradient=True)
