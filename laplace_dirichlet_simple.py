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

x0, x1 = geom
f = function.sin(x0) * function.exp(x1)

# construct dirichlet boundary constraints
cons = domain.boundary.project(
    f, onto=basis, geometry=geom, ischeme='gauss3')

# solve linear system
w = A.solve(constrain=cons)
