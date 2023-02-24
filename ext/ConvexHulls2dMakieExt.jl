module ConvexHulls2dMakieExt

import ConvexHulls2d as CH
import Makie as MK

function makie_points(ch::CH.ConvexHull)
    pts = collect(MK.Point2f, CH.vertices(ch))
    push!(pts, first(pts))
    pts
end

function MK.convert_arguments(trait::MK.PointBased, ch::CH.ConvexHull)
    MK.convert_arguments(trait, makie_points(ch))
end

function MK.plottype(::CH.ConvexHull)
    MK.Lines
end

end#module
