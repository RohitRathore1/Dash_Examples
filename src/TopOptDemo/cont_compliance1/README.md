## Demo of COnt Compliance 1D 

![A demo of the usage DashVtk-TopOPt.jl app](./demo.png)

## Goal

The goal of that example is to show you how you can use DashVtk
to get TopOpt app and render it using DashVtk on the client side.

For that specific example we rely on `GeometryRepresentation` 
and encapsulate in the following structure.

```julia
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
```