#! /usr/bin/env python3

from nutils import *

def main(nelems=8, degree=1):

    # construct topology, geometry and basis
    verts = numpy.linspace(-0.5**0.5, 0.5**0.5, nelems+1)
    domain, geom = mesh.rectilinear([verts, verts])
    basis = domain.basis('spline', degree=degree)
    ischeme = 'gauss3'

    # construct matrix
    A = domain.integrate(
        basis['i,k'] * basis['j,k'],
        geometry=geom, ischeme=ischeme)

    x0, x1 = geom
    f = function.sin(x0) * function.exp(x1)

    # construct dirichlet boundary constraints
    cons = domain.boundary.project(
        f, onto=basis, geometry=geom, ischeme=ischeme)

    # solve linear system
    w = A.solve(constrain=cons)

    # construct solution
    u = basis.dot(w)

    # plot
    points, colors = domain.elem_eval(
        [geom, u], ischeme='bezier3', separate=True)
    with plot.PyPlot('solution') as plt:
        plt.mesh(points, colors)
        plt.colorbar()

if __name__ == '__main__':
    util.run(main)
