#! /usr/bin/env python3

from nutils import *

# construct topology, geometry and basis
verts = numpy.linspace(-0.5**0.5, 0.5**0.5, nelems+1)
domain, geom = mesh.rectilinear([verts, verts])
basis = domain.basis('spline', degree=degree)

# populate a Namespace
ns = function.Namespace()
ns.x = geom
ns.φ = basis
ns.uh = 'φ_n ?w_n'
ns.f = 'sin(x_0) exp(x_1)'

# define the weak formulation
kwargs = dict(geometry=ns.x, degree=2*degree)
res = domain.integral('φ_n,i uh_,i - φ_n f' @ ns, **kwargs)

# construct dirichlet boundary constraints
sqr = domain.boundary.integral('(uh - f)^2' @ ns, **kwargs)
cons = solver.optimize('w', sqr, droptol=1e-1)

# solve linear system
dofs = solver.solve_linear('w', res, constrain=cons)
