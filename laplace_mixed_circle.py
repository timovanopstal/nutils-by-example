#! /usr/bin/env python3

from nutils import *

def main(nelems=8, degree=1):

    # construct topology, geometry and basis
    verts = numpy.linspace(-numpy.pi/4, numpy.pi/4, nelems+1)
    domain, geom = mesh.rectilinear([verts, verts])
    geom = function.stack([
        function.sin(geom[0]) * function.cos(geom[1]),
        function.cos(geom[0]) * function.sin(geom[1]),
    ])
    basis = domain.basis('spline', degree=degree)
    ischeme = 'gauss3'

    # construct matrix
    A = domain.integrate(
        basis['i,k'] * basis['j,k'],
        geometry=geom, ischeme=ischeme)

    x0, x1 = geom
    f = function.sin(x0) * function.exp(x1)

    # construct dirichlet boundary constraints
    cons = domain.boundary['left,bottom'].project(
        f, onto=basis, geometry=geom, ischeme=ischeme)

    # construct right hand side
    n = geom.normal()
    b = domain.boundary['right,top'].integrate(
        basis['i'] * f[',k'] * n['k'],
        geometry=geom, ischeme=ischeme)

    # solve linear system
    w = A.solve(b, constrain=cons)

    # construct solution
    u = basis.dot(w)

    # plot
    points, colors = domain.elem_eval(
        [geom, u], ischeme='bezier3', separate=True)
    with plot.PyPlot('solution') as plt:
        plt.mesh(points, colors)
        plt.colorbar()
        plt.clim(-1.318, 1.318)

if __name__ == '__main__':
    util.run(main)
