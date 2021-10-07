## Demo of Point Cloud Representation

![A demo of Point Cloud Representation](demo.gif)

```julia
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
```