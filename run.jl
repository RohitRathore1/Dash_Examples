using DashExamples

app = DashExamples.make_app()

port = haskey(ENV, "PORT") ? parse(Int64, ENV["PORT"]) : 8050

DashExamples.Dash.run_server(app, "0.0.0.0", port)