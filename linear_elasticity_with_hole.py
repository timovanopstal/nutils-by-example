#! /usr/bin/env python3

from nutils import *

def main(nelems=12, stress=library.Hooke(lmbda=1,mu=1), degree=2):

    # construct topology, geometry and basis
    verts = numpy.linspace(0, 1, nelems+1)
    domain, geom = mesh.rectilinear([verts,verts])

    # cut a hole with radius .2
    levelset = function.norm2(geom - (.6,.4)) - .2
    domain = domain.trim(levelset, maxrefine=3)

    dbasis = domain.basis('spline', degree=degree).vector(2)
    ischeme = 'gauss2'

    # construct matrix
    A = domain.integrate(
        dbasis['ik,l']*stress(dbasis.symgrad(geom))['jkl'],
        geometry=geom, ischeme=ischeme)

    # construct dirichlet boundary constraints
    cons = \
        domain.boundary['left'].project(
            0.0, geometry=geom,
            onto=dbasis, ischeme=ischeme) \
        | domain.boundary['right'].project(
            0.5, geometry=geom,
            onto=dbasis.dotnorm(geom), ischeme=ischeme)

    # solve system
    w = A.solve(constrain=cons)

    # construct solution function
    disp = dbasis.dot(w)

    # plot solution
    points, colors = domain.simplex.elem_eval(
        [geom+disp, stress(disp.symgrad(geom))[0,1]],
        ischeme='bezier3', separate=True)
    with plot.PyPlot('stress') as plt:
        plt.mesh(points, colors, tight=False)
        plt.colorbar()

util.run(main)
