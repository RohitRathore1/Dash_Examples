using TopOpt, LinearAlgebra, StatsFuns
using TopOpt.TrussTopOptProblems
using TopOptMakie
using VTKDataIO
using PyCall
using VTKDataTypes
using DashVtk
using Dash
using DashHtmlComponents

#dash = pyimport("dash")
#dash_html = pyimport("dash_html_components")
dash_vtk = pyimport("dash_vtk")

# 2D
ndim = 2
node_points, elements, mats, crosssecs, fixities, load_cases = 
    parse_truss_json(joinpath(@__DIR__, "tim_$(ndim)d.json"));
ndim, nnodes, ncells = length(node_points[1]), length(node_points), length(elements)
loads = load_cases["0"]
problem = TrussProblem(
    Val{:Linear}, node_points, elements,
    loads, fixities, mats, crosssecs,
);

xmin = 0.0001 # minimum density
x0 = fill(1.0, ncells) # initial design
p = 4.0 # penalty
V = 0.5 # maximum volume fraction

solver = FEASolver(Displacement, Direct, problem, xmin = xmin)
comp = TopOpt.Compliance(problem, solver)

function obj(x)
    # minimize compliance
    return comp(x)
end
function constr(x)
    # volume fraction constraint
    return sum(x) / length(x) - V
end

m = Model(obj)
addvar!(m, zeros(length(x0)), ones(length(x0)))
Nonconvex.add_ineq_constraint!(m, constr)

options = MMAOptions(
    maxiter=1000, tol = Tolerance(kkt = 1e-4, f = 1e-4),
)
TopOpt.setpenalty!(solver, p)
@time r = Nonconvex.optimize(
    m, MMA87(dualoptimizer = ConjugateGradient()),
    x0, options = options,
);

@show obj(r.minimizer)
@show constr(r.minimizer)

topology = r.minimizer;

mesh = VTKUnstructuredData(problem)
topology = round.(topology)
inds = findall(isequal(0), topology)
deleteat!(mesh.cell_connectivity, inds)
deleteat!(mesh.cell_types, inds)

pymesh = PyVTK(mesh)
dash_vtk_utils = pyimport("dash_vtk.utils")
dash_vtk_utils.to_mesh_state
state = dash_vtk_utils.to_mesh_state(pymesh)

content = vtk_view([
        vtk_geometryrepresentation([
        vtk_mesh(state=state)
    ]),
])

app = dash()

app.layout = html_div() do 
    html_div(
    style = Dict(
        "width" => "100%", 
        "height" => "400px"
    ),
    children=[content],
    )
end

run_server(app, "0.0.0.0", debug = true)