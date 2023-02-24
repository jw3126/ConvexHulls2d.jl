using ConvexHulls2d
import ConvexHulls2d as CH
using Test
using LinearAlgebra
using StaticArrays
import Random: MersenneTwister

@testset "edge_function" begin
    using ConvexHulls2d: edge_function
    @test edge_function([0,0], [1,0], [0,1]) ≈ 1
    @test edge_function([0,0], [1,0], [0,0]) ≈ 0 atol = 1e-14
    @test edge_function([0,0], [1,0], [0,-1]) ≈ -1 atol = 1e-14
end

@testset "fuzz edge_function" begin
    using ConvexHulls2d: edge_function
    rng = MersenneTwister(1234)
    for _ in 1:100 
        p1 = randn(rng, 2)
        p2 = randn(rng, 2)
        p3 = randn(rng, 2)
        A = randn(rng, 2,2)
        b = randn(rng, 2)
        @test det(A)*edge_function(p1,p2,p3) ≈ edge_function(A*p1+b, A*p2+b, A*p3+b)
    end
end

function test_hull(h; atol=0)
    for pt in h.points
        @test CH.signed_distance(h, pt) <= atol
    end
    nedges = length(CH.vertices(h))
    if nedges >= 3
        keys = [h.indices; h.indices[1]; h.indices[2]]
        for i in 1:nedges
            key1 = keys[i+0]
            key2 = keys[i+1]
            key3 = keys[i+2]
            p1 = h.points[key1]
            p2 = h.points[key2]
            p3 = h.points[key3]
            @test CH.edge_function(p1,p2,p3) >= -atol
        end
    end
end

@testset "ConvexHull square" begin
    pts = [
        [0.0, 0.0],
        [0.0, 1.0],
        [1.0, 0.0],
        [1.0, 1.0],
    ]
    p1,p2,p3,p4 = pts
    h = @inferred CH.ConvexHull(pts)
    test_hull(h)
    @inferred CH.signed_distance(h, p1)
    @inferred CH.distance(h, p1)
    @test CH.signed_distance(h, p1) ≈ 0
    @test CH.signed_distance(h, p2) ≈ 0
    @test CH.signed_distance(h, p3) ≈ 0
    @test CH.signed_distance(h, p4) ≈ 0
    @test CH.signed_distance(h, [0.5,0.5]) ≈ -0.5
    @test CH.distance(h, [0.5,0.5]) == 0.0
    @test CH.signed_distance(h, [0.5,0.9]) ≈ -0.1
    @test CH.signed_distance(h, [1.5,0.5]) ≈ 0.5
    @test CH.signed_distance(h, [3,2.0]) ≈ sqrt(5)
    @test CH.distance(h, [3,2.0]) ≈ sqrt(5)
    @inferred(CH.vertices(h))
    @inferred(CH.circumference(h))
    @inferred(CH.indices(h))
    @test CH.circumference(h) ≈ 4
    @test length(CH.vertices(h)) == 4
    @test CH.vertices(h) == [pts[1], pts[3], pts[4], pts[2]]
end


@testset "ConvexHull diamond" begin
    pts = [[0,0], [1,1], [1,-1], [2,0]]
    h = CH.ConvexHull(pts)
    test_hull(h)
    @test CH.circumference(h) ≈ 4*sqrt(2)
    @test CH.distance(@SVector[0,0], h) == 0
    @test CH.distance(@SVector[1,0], h) == 0
    @test CH.distance(@SVector[-1,0], h) > 0
    @test CH.distance(@SVector[1,0.99999], h) == 0
    @test CH.distance(@SVector[1,1.000001], h) > 0
end
#
#@testset "ConvexHull duplicate points" begin
#    pts = [
#        [-9.892169405717053e-7, -1.5874977989795982],
#        [1.5874985415154759, 0.0],
#        [-9.892169405717053e-7, -1.5874977989795983],
#        [0.24833866271737112, -1.5679550911531746],
#    ]
#    h = ConvexHull(pts)
#    @test circumference(h) ≈ norm(pts[1] - pts[2]) +
#        norm(pts[4] - pts[2]) +
#        norm(pts[4] - pts[1])
#    for pt in pts
#        @test isinside(pt, h)
#    end
#    @test allunique(h.indices)
#end
#
#@testset "ConvexHull duplicate points 2" begin
#    pts = [
#        [1.5874995012976747, 0.0],
#        [-5.838217168722755e-7, 1.5874995095689919],
#        [0.03589953217296475, 1.58709324677428],
#        [-5.838217168722755e-7, 1.5874995095689919],
#        [1.5855872826713822, 0.07789421068418088],
#    ]
#
#    h = ConvexHull(pts)
#    @test circumference(h) ≈ norm(pts[2] - pts[1]) +
#        norm(pts[1] - pts[5]) +
#        norm(pts[5] - pts[3]) +
#        norm(pts[3] - pts[4]) +
#        norm(pts[4] - pts[2])
#    for pt in pts
#        @test isinside(pt, h)
#    end
#    @test allunique(h.indices)
#end
#
#@testset "ConvexHull fuzz" begin
#    rng = MersenneTwister(1337)
#    for _ in 1:100
#        npts = rand(rng, 2:50)
#        pts = [SVector(randn(rng, 2)...) for _ in 1:npts]
#        h = ConvexHull(pts)
#        for pt in pts
#            # if !isinside(pt, h)
#            #     foreach(println, pts)
#            #     error()
#            # end
#            @test isinside(pt, h)
#        end
#        @test !isinside(@SVector[0,10], h)
#        @test allunique(h.indices)
#    end
#end
