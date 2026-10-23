using Pkg
Pkg.activate(@__DIR__)
using Literate
Literate.notebook(
    joinpath(@__DIR__, "fermionic_one_site_scan.jl"),
    @__DIR__;
    execute = false,
    name = "fermionic_one_site_scan",
)
println("Wrote fermionic_one_site_scan.ipynb")
