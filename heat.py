#! /usr/bin/env python3

from nutils import *

def main(
      alpha: 'thermal diffusivity' = 0.01,
      nelems: 'number of elements in one direction' = 64,
      degree: 'polynomial degree of spline basis' = 2,
      dt: 'time step' = 0.01,
      tend: 'simulation time' = 1
    ):

    # construct topology, geometry and basis
    verts = numpy.linspace(0, 1, nelems+1)
    domain, geom = mesh.rectilinear([verts, verts])
    basis = domain.basis('spline', degree=degree)

    # populate namespace
    ns = function.Namespace()
    ns.x = geom
    ns.φ = basis
    ns.uh = 'φ_i ?w_i'
    ns.dt = dt
    ns.α = alpha

    # construct initial condition
    ns.l1 = '0.15 - (abs(0.3 - x_0) + abs(0.4 - x_1))'
    ns.l2 = '0.15 - ((0.7 - x_0)^2 + (0.6 - x_1)^2)^0.5'
    ns.l = function.max(ns.l1, ns.l2)
    ns.c = nelems/2
    ns.u0 = '0.5 + 0.5 tanh(c l)'
    kwargs = dict(geometry=ns.x, degree=2*degree)
    sqr = domain.integral('(uh - u0)^2' @ ns, **kwargs)
    w0 = solver.optimize('w', sqr, droptol=1e-10)

    # define the weak formulation
    res = domain.integral('0.5 α uh_,i φ_n,i' @ ns, **kwargs)
    inertia = domain.integral('uh φ_n' @ ns, **kwargs)
    
    # construct dirichlet boundary constraints
    sqr = domain.boundary.integral('uh^2' @ ns, **kwargs)
    cons = solver.optimize('w', sqr, droptol=1e-1)

    # time integration
    stepper = solver.cranknicolson('w', residual=res, inertia=inertia,
        lhs0=w0, timestep=dt, constrain=cons)
    for n, dofs in log.enumerate('timestep', stepper):

        # break if simulation time reached
        if n*dt > tend: break

        # construct solution
        sol = 'uh' @ ns(w=dofs)

        # plot
        points, colors = domain.elem_eval(
            [geom, sol], ischeme='bezier3', separate=True)
        with plot.PyPlot('temperature') as plt:
            plt.title('t={:5.2f}'.format(n*dt))
            plt.mesh(points, colors)
            plt.colorbar()
            plt.clim(0, 1)

if __name__ == '__main__':
    cli.run(main)
