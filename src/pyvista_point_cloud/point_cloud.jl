using Dash
using DashVtk
using DashHtmlComponents
using DashCoreComponents
using PyCall
using VTKDataIO
using VTKDataTypes

dash_vtk = pyimport("dash_vtk")
vtk = pyimport("vtk")
dash_vtk_utils = pyimport("dash_vtk.utils")
pv = pyimport("pyvista")
np = pyimport("numpy")
examples = pyimport("pyvista.examples")

rand(42)
dataset = examples.download_lidar()
fp = _VTKDataTypes(dataset)
#pydataset = PyVTK(dataset)
subset = 0.2
selection = np.random.randint(
    low=1, high=dataset.n_points, size=round(Int,dataset.n_points * subset)
)
#selection = [1 , 3, 5]
n, m = size(dataset.points)
#jlpoints = copy(reshape(vec(dataset.points), m, n)')

points = dataset.points[selection, :]
xyz = points'[:]
elevation = points[:, end]
min_elevation = minimum(elevation)
max_elevation = maximum(elevation)
app = dash()

content = vtk_view([
    vtk_pointcloudrepresentation(
        xyz=xyz,
        scalars=elevation,
        colorDataRange=[min_elevation, max_elevation],
        property=Dict(
            "pointSize" => 2
        ),
    )
])

app.layout = html_div() do 
        html_div(
            style=Dict(
                "height" => "calc(100vh - 16px)"
            ),
            children=[(
                html_div(
                    content,
                    style=Dict(
                        "height" => "100%",
                        "width" => "100%",
                    )
                )
            )],
        )
end

run_server(app, "0.0.0.0", debug = true)