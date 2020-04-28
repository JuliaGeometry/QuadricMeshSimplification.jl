using QuadricMeshSimplification
using Test
using Meshing
using GeometryBasics

@testset "QuadricMeshSimplification.jl" begin
    #m = GLNormalMesh(Sphere(zero(Point{3,Float32}),1.0))
    #simplify_mesh_lossless!(vertices(m),faces(m))
    @testset "Sphere" begin
        sphere(v) = sqrt(sum(v.^2)) - 1
        vts, fcs = isosurface(sphere,MarchingTetrahedra())
        simplify_mesh!(vts,fcs, length(vts)>>2)
    end
    @testset "random" begin
        vts, fcs = isosurface(rand(20,20,20),MarchingCubes(0.5))
        simplify_mesh!(vts,fcs, length(vts)>>2)
    end
end
