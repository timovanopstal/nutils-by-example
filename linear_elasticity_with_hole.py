#! /usr/bin/env python3

from nutils import *

def main(
      nelems: 'number of elements in one direction' = 8,
      degree: 'polynomial degree of spline basis' = 2,
      lmbda: 'first Lamé parameter' = 1,
      mu: 'second Lamé parameter' = 1
    ):

    # construct topology and geometry
    verts = numpy.linspace(0, 1, nelems+1)
    domain, geom = mesh.rectilinear([verts,verts])

    # cut a hole with radius .2
    levelset = function.norm2(geom - (.6,.4)) - .2
    domain = domain.trim(levelset, maxrefine=3)

    # construct a basis
    basis = domain.basis('spline', degree=degree).vector(2)

    ns = function.Namespace()
    ns.x = geom
    ns.φ = basis
    ns.uh_k = 'φ_nk ?w_n'
    ns.λ = lmbda
    ns.μ = mu
    ns.ε_ij = '(uh_i,j + uh_j,i) / 2'
    ns.σ_ij = 'λ ε_kk δ_ij + 2 μ ε_ij'
    
    # define the weak formulation
    kwargs = dict(geometry=ns.x, degree=2*degree)
    res = domain.integral('φ_ni,j σ_ij' @ ns, **kwargs)
    
    # construct dirichlet boundary constraints
    cons = domain.boundary['left'].project(0., onto=basis, **kwargs) \
         | domain.boundary['right'].project(0.5, onto=basis.dotnorm(geom),
         **kwargs)
    
    # solve linear system
    dofs = solver.solve_linear('w', res, constrain=cons)

    # construct solution
    x = 'x_k + uh_k' @ ns(w=dofs)
    τ = 'σ_01' @ ns(w=dofs)

    # plot
    points, colors = domain.simplex.elem_eval(
        [x, τ], ischeme='bezier3', separate=True)
    with plot.PyPlot('stress') as plt:
        plt.mesh(points, colors)
        plt.colorbar()

if __name__ == '__main__':
    cli.run(main)
