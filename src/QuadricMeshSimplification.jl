module QuadricMeshSimplification


# Mesh Simplification Tutorial
#
# (C) by Sven Forstmann in 2014
# (C) by Steve Kelly in 2020
#
# License : MIT
# http://opensource.org/licenses/MIT
#
# Derived from:
# https://github.com/sp4cerat/Fast-Quadric-Mesh-Simplification
#
# 5/2016: Chris Rorden created minimal version for OSX/Linux/Windows compile
# 4/2020: Steve Kelly Julia port

using StaticArrays
using LinearAlgebra

export simplify_mesh!, simplify_mesh_lossless!

function SymetricMatrix(a, b, c, d)
    SVector{10,Float64}(a*a, a*b, a*c, a*d, b*b, b*c, b*d, c*c, c*d, d*d)
end

function det(m::SVector, a11, a12, a13, a21, a22, a23, a31, a32, a33)
    m[a11]*m[a22]*m[a33] + m[a13]*m[a21]*m[a32] + m[a12]*m[a23]*m[a31]
                - m[a13]*m[a22]*m[a31] - m[a11]*m[a23]*m[a32]- m[a12]*m[a21]*m[a33];
end


mutable struct Triangle
    v::SVector{3,Int}
    err::SVector{4,Float64}
    deleted::Bool
    dirty::Bool
    n::SVector{3,Float64}
end

mutable struct Vertex{PT}
    p::PT
    tstart::Int
    tcount::Int
    q::SVector{10,Float64}
    border::Bool
end

struct TRef
    tid::Int
    tvertex::Int
end

# Helper functions

"""
Main simplification function

 target_count  : target nr. of triangles
 agressiveness : sharpness to increase the threshold.
                 5..8 are good numbers
                 more iterations yield higher quality

"""
function simplify_mesh!(invts, infcs, target_count::Int, agressiveness=7, verbose=false)
    # init
    triangles = Vector{Triangle}(undef,length(infcs))
    vertices = Vector{Vertex{eltype(invts)}}(undef,length(invts))

    for i in eachindex(infcs)
        triangles[i]=Triangle(infcs[i],zero(SVector{4,Float64}),false,false,zero(SVector{3,Float64}))
    end
    for i in eachindex(invts)
        vertices[i]=Vertex(invts[i],0,0,zero(SVector{10,Float64}),false)
    end
    # main iteration loop
    deleted_triangles=0;
    deleted0 = Bool[]
    deleted1 = Bool[]
    refs = TRef[]

    triangle_count=length(triangles)

    for iteration = 0:100
        triangle_count-deleted_triangles<=target_count && break

        # update mesh once in a while
        @show iteration
        iszero(iteration%5) && update_mesh!(triangles, vertices, refs, iteration)

        # clear dirty flag
        for t in triangles
            t.dirty = false
        end

        #
        # All triangles with edges below the threshold will be removed
        #
        # The following numbers works well for most models.
        # If it does not, try to adjust the 3 parameters
        #
        threshold = 0.000000001*(iteration+3)^agressiveness

        # remove vertices & mark deleted triangles
        for i in eachindex(triangles)
            t = triangles[i]
            t.err[4]>threshold && continue
            t.deleted && continue
            t.dirty && continue

            for j = 1:3
                if t.err[j] < threshold

                    i0 = t.v[j]
                    v0 = vertices[i0]
                    i1 = t.v[(j%3)+1]
                    v1 = vertices[i1]
                    # Border check
                    v0.border != v1.border && continue

                    # Compute vertex to collapse to
                    e, p = calculate_error(vertices,i0,i1)
                    resize!(deleted0,v0.tcount) # normals temporarily
                    resize!(deleted1,v1.tcount) # normals temporarily
                    # don't remove if flipped
                    flipped(vertices,triangles,refs,p,i0,i1,v0,v1,deleted0) && continue

                    flipped(vertices,triangles,refs,p,i1,i0,v1,v0,deleted1) && continue

                    # not flipped, so remove edge
                    v0.p = p
                    v0.q = v1.q .+ v0.q
                    tstart = length(refs)

                    deleted_triangles = update_triangles(triangles,vertices,refs,i0,v0,deleted0,deleted_triangles)
                    deleted_triangles = update_triangles(triangles,vertices,refs,i0,v1,deleted1,deleted_triangles)

                    tcount = length(refs) - tstart

                    if tcount <= v0.tcount
                        # save ram
                        s = v0.tstart
                        if !iszero(tcount)
                            for ne in 0:(tcount-1)
                                refs[s+ne] = refs[tstart+ne]
                            end
                        end
                    else
                        # append
                        v0.tstart = tstart
                    end
                    v0.tcount = tcount
                    break
                end
            end
            # done?
            triangle_count-deleted_triangles <= target_count && break
        end
    end
    # clean up mesh
    compact_mesh(vertices,triangles)
end #simplify_mesh()


"""
Check if a triangle flips when this edge is removed
"""
function flipped(vertices, triangles, refs, p, i0, i1, v0::Vertex, v1::Vertex, deleted::Vector{Bool})
    for k in 1:v0.tcount
        t=triangles[refs[v0.tstart+k].tid]

        t.deleted && continue

        s=refs[v0.tstart+k].tvertex
        id1=t.v[(s%3)+1]
        id2=t.v[((s+1)%3)+1]

        if id1==i1 || id2==i1 # delete ?
            deleted[k] = true
            continue
        end
        d1 = normalize(vertices[id1].p .- p)
        d2 = normalize(vertices[id2].p .- p)
        abs(dot(d1, d2))>0.999 && return true
        n = normalize(cross(d1,d2))
        deleted[k]=false
        dot(n,t.n) < 0.2 && return true
    end
    return false
end

# Update triangle connections and edge error after a edge is collapsed
function update_triangles(triangles, vertices, refs, i0, v::Vertex, deleted, deleted_triangles)
    dt = deleted_triangles
    for k in 1:v.tcount
        r = refs[v.tstart+k-1]
        t = triangles[r.tid]
        t.deleted && continue
        if deleted[k]
            t.deleted=true
            dt += 1
            continue
        end
        ind = r.tvertex
        tv = t.v
        if ind == 1
            triangles[r.tid].v = SVector(i0,tv[2],tv[3])
        elseif ind == 2
            triangles[r.tid].v = SVector(tv[1],i0,tv[3])
        else
            triangles[r.tid].v = SVector(tv[1],tv[2],i0)
        end
        triangles[r.tid].dirty = true
        err0, p = calculate_error(vertices, tv[1], tv[2])
        err1, p = calculate_error(vertices, tv[2], tv[3])
        err2, p = calculate_error(vertices, tv[3], tv[1])
        err3 = min(err0, err1, err2)
        t.err = SVector{4,Float64}(err0,err1,err2,err3)
        push!(refs, refs[v.tstart+k])
    end
    dt
end


function update_mesh!(triangles::Vector{Triangle}, vertices::Vector, refs::Vector{TRef}, iteration)
    if iteration > 0 # compact triangles
        dst = 1
        for i in eachindex(triangles)
            if !(triangles[i].deleted)
                triangles[dst]=triangles[i]
                dst += 1
            end
        end
        resize!(triangles,dst)
    end
    #
    # Init Quadrics by Plane & Edge Errors
    #
    # required at the beginning ( iteration == 0 )
    # recomputing during the simplification is not required,
    # but mostly improves the result for closed meshes
    #
    if iszero(iteration)
        for i in eachindex(vertices)
            vertices[i].q = zero(SVector{10,Float64})
        end

        for i in eachindex(triangles)
            t = triangles[i]
            p1 = vertices[t.v[1]].p
            p2 = vertices[t.v[2]].p
            p3 = vertices[t.v[3]].p
            n = normalize(cross(p2 .- p1, p3 .- p1))
            triangles[i].n = n
            vertices[t.v[1]].q += SymetricMatrix(n[1],n[2],n[3],dot(-n,p1))
            vertices[t.v[2]].q += SymetricMatrix(n[1],n[2],n[3],dot(-n,p1))
            vertices[t.v[3]].q += SymetricMatrix(n[1],n[2],n[3],dot(-n,p1))
        end
        for i in eachindex(triangles)
            # Calc Edge Error
            t = triangles[i]
            err1, _ = calculate_error(vertices,t.v[1],t.v[2])
            err2, _ = calculate_error(vertices,t.v[2],t.v[3])
            err3, _ = calculate_error(vertices,t.v[3],t.v[1])
            err4 = min(err1,err2,err3)
            t.err = SVector{4,Float64}(err1,err2,err3,err4)
        end
    end

    # Init Reference ID list
    for i in eachindex(vertices)
        vertices[i].tstart = 0
        vertices[i].tcount = 0
    end
    @show "hello"
    for i in eachindex(triangles)
        for j in 1:3
            #@show triangles[i].v[j], length(vertices), i, length(triangles)
            vertices[triangles[i].v[j]].tcount += 1
        end
    end
    @show "exited"
    tstart = 1
    for i in eachindex(vertices)
        vertices[i].tstart = tstart
        tstart += vertices[i].tcount
        vertices[i].tcount = 1
    end

    # Write References
    resize!(refs, length(triangles)*3)
    for i in eachindex(triangles)
        t=triangles[i]
        for j in 1:3
            v=vertices[t.v[j]]
            refs[v.tstart+v.tcount-1] = TRef(i, j)
            vertices[t.v[j]].tcount += 1
        end
    end

    # Identify boundary : vertices[].border=0,1
    if iszero(iteration)
        vcount = Int[]
        vids = Int[]

        for i in eachindex(vertices)
            vertices[i].border = false
        end

        for i in eachindex(vertices)
            v=vertices[i]
            empty!(vcount)
            empty!(vids)
            for j in 1:(v.tcount-1)
                k = refs[v.tstart+j-1].tid
                t = triangles[k]
                for k = 1:3
                    ofs = 0
                    id = t.v[k]
                    while ofs < length(vcount)
                        vids[ofs+1] == id && break
                        ofs += 1
                    end
                    if ofs == length(vcount)
                        push!(vcount,1)
                        push!(vids, id)
                    else
                        vcount[ofs+1] += 1
                    end
                end
            end
            for j in eachindex(vcount)
                if vcount[j] == 1
                    vertices[vids[j]].border = true
                end
            end
        end
    end
end

# Finally compact mesh before exiting
function compact_mesh(vertices::Vector,triangles::Vector{Triangle})
    dst = 1
    for i in eachindex(vertices)
        vertices[i].tcount=0
    end
    for i in eachindex(triangles)
        if !(triangles[i].deleted)
            t = triangles[i]
            triangles[dst]=t
            dst += 1
            for j in 1:3
                vertices[t.v[j]].tcount = 1
            end
        end
    end
    resize!(triangles, dst)
    dst = 1
    for i in eachindex(vertices)
        if !iszero(vertices[i].tcount)
            vertices[i].tstart=dst;
            vertices[dst].p = vertices[i].p
            dst += 1
        end
    end
    for i in eachindex(triangles)
        @show triangles[i]
        @show length(vertices), triangles[i].v[1],triangles[i].v[2],triangles[i].v[3]
        triangles[i].v = SVector(vertices[triangles[i].v[1]].tstart,
                                 vertices[triangles[i].v[2]].tstart,
                                 vertices[triangles[i].v[3]].tstart)
    end
    resize!(vertices,dst-1)
end

"""
Error between vertex and Quadric
"""
function vertex_error(q::SVector{10,Float64}, x::Float64, y::Float64, z::Float64)
    q[1]*x*x + 2*q[2]*x*y + 2*q[3]*x*z + 2*q[4]*x + q[5]*y*y + 2*q[6]*y*z + 2*q[7]*y + q[8]*z*z + 2*q[9]*z + q[10]
end


"""
Error for one edge
"""
function calculate_error(vertices::Vector, id_v1::Int, id_v2::Int)
    # compute interpolated vertex

    q = vertices[id_v1].q .+ vertices[id_v2].q
    border = vertices[id_v1].border & vertices[id_v2].border
    d = det(q,1, 2, 3, 2, 5, 6, 3, 6, 8)

    if !iszero(d) && !border
        # q_delta is invertible
        x = -1/d*(det(q,2, 3, 4, 5, 6, 7, 6, 8, 9)) # vx = A41/det(q_delta)
        y =  1/d*(det(q,1, 3, 4, 2, 6, 7, 3, 8, 9)) # vy = A42/det(q_delta)
        z = -1/d*(det(q,1, 2, 4, 2, 5, 7, 3, 6,  9)) # vz = A43/det(q_delta)

        return vertex_error(q, x, y, z), SVector{3,Float64}(x,y,z)
    else
        # det = 0 -> try to find best result
        p1 = vertices[id_v1].p
        p2 = vertices[id_v2].p
        p3 = (p1.+p2)/2
        error1 = vertex_error(q, p1[1],p1[2],p1[3])
        error2 = vertex_error(q, p2[1],p2[2],p2[3])
        error3 = vertex_error(q, p3[1],p3[2],p3[3])
        error4 = min(error1, error2, error3)
        if error1 == error4
            return error4, p1
        elseif error2 == error4
            return error4, p2
        elseif error3 == error4
            return error4, p3
        end
    end
    return 0.0, zero(SVector{3,Float64})
end


end # module
