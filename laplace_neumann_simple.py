#! /usr/bin/env python3

from nutils import *

# construct topology, geometry and basis
verts = numpy.linspace(-0.5**0.5, 0.5**0.5, 9)
domain, geom = mesh.rectilinear([verts, verts])
basis = domain.basis('spline', degree=1)

# populate a Namespace
ns = function.Namespace()
ns.x = geom
ns.φ = basis
ns.uh = 'φ_n ?w_n'

# define the weak formulation
res = domain.integral('φ_n,i uh_,i' @ ns,
                      geometry=ns.x, ischeme='gauss1')

# solve linear system
dofs = solver.solve_linear('w', res)
