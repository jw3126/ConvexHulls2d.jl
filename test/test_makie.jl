module TestMakie
using Test
import ConvexHulls2d as CH
using Makie

@testset "smoketest Makie" begin
    ch = CH.ConvexHull([randn(2) for _ in 1:10])
    plot(ch)
end

end#module
