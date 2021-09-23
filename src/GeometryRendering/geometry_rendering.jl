using VTKDataIO
using PyCall
using VTKDataTypes
using DashVtk
using Dash
using DashHtmlComponents

#dash = pyimport("dash")
#dash_html = pyimport("dash_html_components")
dash_vtk = pyimport("dash_vtk")
x = y = z = [-2, -1, 0, 1, 2];
rect = VTKRectilinearData((x, y, z));
rect.cell_data["Cell scalar"] = reshape([rand() for i in 1:num_of_cells(rect)], cell_extents(rect));
unstruct = VTKUnstructuredData(rect)
pyunstruct = PyVTK(unstruct)
dash_vtk_utils = pyimport("dash_vtk.utils")

dash_vtk_utils.to_mesh_state
state = dash_vtk_utils.to_mesh_state(pyunstruct)

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
