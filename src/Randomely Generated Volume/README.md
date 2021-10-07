## Demo of Randomely Generated Volume

![A demo of Randomely Generated Volume](demo.gif)

```julia
volume_view = vtk_view(
   children=vtk_volumedatarepresentation(
     spacing=[1, 1, 1],
     dimensions=[10, 10, 10],
     origin=[0, 0, 0],
     scalars = [rand() for _ in 1:1000],
     rescaleColorMap=false,
   )
)
```