#! /usr/bin/env python3

from nutils import *

# construct topology, geometry and basis
verts = numpy.linspace(-0.5**0.5, 0.5**0.5, 9)
domain, geom = mesh.rectilinear([verts, verts])
basis = domain.basis('spline', degree=1)

# construct matrix
A = domain.integrate(
    basis['i,k'] * basis['j,k'],
    geometry=geom, ischeme='gauss3')

# solve linear system
w = A.solve()
