using QuadricMeshSimplification
using Test
using Meshing

@testset "QuadricMeshSimplification.jl" begin
    vts, fcs = isosurface(rand(20,20,20),MarchingCubes(0.5))
    simplify_mesh_lossless!(vts,fcs)
end
