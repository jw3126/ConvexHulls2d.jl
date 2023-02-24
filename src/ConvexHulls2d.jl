module ConvexHulls2d
import Base.Iterators
using LinearAlgebra

struct ConvexHull{Pts}
    points::Pts
    indices::Vector{Int}
end
function vertices(h::ConvexHull)
    view(h.points, h.indices)
end
function nvertices(h)
    length(h.indices)
end
function indices(h::ConvexHull)
    h.indices
end

struct GrahamScan end

function ConvexHull(pts::AbstractVector) 
    ConvexHull(pts, GrahamScan())
end

function ConvexHull(pts::AbstractVector, ::GrahamScan)
    indices = graham_scan(pts)
    return ConvexHull(pts, indices)
end

macro foreachedge(edgename, h, code)
    edge = esc(edgename)
    # TODO inbounds
    quote
        h = $(esc(h))
        i_first = firstindex(h.indices)
        i_last  = lastindex(h.indices)
        for i in i_first:(i_last-1)
            key1 = h.indices[i]
            key2 = h.indices[i+1]
            pt1 = h.points[key1]
            pt2 = h.points[key2]
            $edge = (pt1,pt2)
            $(esc(code))
        end
        key1 = h.indices[i_last]
        key2 = h.indices[i_first]
        pt1 = h.points[key1]
        pt2 = h.points[key2]
        $edge = (pt1, pt2)
        $(esc(code))
    end
end


function distance(p1::AbstractVector, p2::AbstractVector)
    x1,y1 = p1
    x2,y2 = p2
    return sqrt((x1-x2)^2 + (y1-y2)^2)
end

function distance(h::ConvexHull, pt::AbstractVector)
    max(0, signed_distance(h, pt))
end
function distance(pt::AbstractVector, h::ConvexHull)
    max(0, signed_distance(h, pt))
end

function numtype(h::ConvexHull)
    float(typeof(h.points[1][1]))
end
function circumference(h::ConvexHull)
    T = numtype(h)
    ret = zero(T)
    @foreachedge (p1, p2) h begin
        ret += distance(p1, p2)
    end
    return ret
end

function distance_edge_point(p1, p2, p)
    v = p2 - p1
    pt = p - p1
    t = dot(pt, v) / sum(abs2, v)
    pt_proj = clamp(t, 0, 1) * v
    distance(pt_proj, pt)
end

function signed_distance(h::ConvexHull, p::AbstractVector)
    T = numtype(h)
    if nvertices(h) == 1
        p0 = only(vertices(h))
        return T(distance(p0, p))
    end
    has_ccw = false
    has_cw = false
    outside_dist = T(Inf)
    inside_dist = T(Inf)
    @foreachedge (p1, p2) h begin
        det = edge_function(p1, p2, p)
        if det > 0
            has_ccw = true
        elseif det < 0
            has_cw = true
        end
        d = distance_edge_point(p1, p2, p)
        if det >= 0
            inside_dist = min(inside_dist, d)
        end
        if det <= 0
            outside_dist = min(outside_dist, d)
        end
    end
    if !has_cw
        # inside
        return -inside_dist
    else
        # outside
        return outside_dist
    end
    return ret
end

function edge_function(p1::AbstractVector, p2::AbstractVector, p3::AbstractVector)
    # v1 = p2 - p1
    # v2 = p3 - p2
    # v1[1]*v2[2] - v1[2]*v2[1]
    A11 = p2[1] - p1[1]
    A21 = p3[1] - p2[1]
    A12 = p2[2] - p1[2]
    A22 = p3[2] - p2[2]
    A11*A22 - A21*A12
end

function new_key!(out, key_new, pts, i0)
    pt_new = pts[key_new]
    if isempty(out)
        push!(out, key_new)
        return out
    end
    while true
        pt2 = pts[out[end]]
        if pt_new == pt2
            break
        end
        if (length(out) == 1+i0)
            push!(out, key_new)
            break
        end
        @assert length(out) > i0+1
        pt1 = pts[out[end-1]]
        is_ccw = edge_function(pt1, pt2, pt_new) > 0
        if is_ccw
            push!(out, key_new)
            break
        else
            pop!(out)
        end
    end
    return out
end

function isless_graham(pt1, pt2)
    # (pt1[1] < pt2[1]) || (pt1[2] < pt2[2])
    <(pt1, pt2)
end

function graham_scan(pts)
    keys = sortperm(pts)#, lt=isless_graham)
    out = Int[]
    for key in keys
        new_key!(out, key, pts, 0)
    end
    @assert !isempty(out)
    @assert pts[last(out)] == pts[last(keys)]
    n = length(out)-1
    for key in Iterators.reverse(keys)
        new_key!(out, key, pts, n)
    end
    @assert pts[last(out)] == pts[first(keys)]
    if length(out) > 1
        pop!(out)
    end
    return out
end

end
