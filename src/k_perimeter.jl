#! /bin/env julial
using Cairo, DataStructures, Printf
import LinearAlgebra

struct Point2DAny{T}
    x::T
    y::T
end

const Point2D = Point2DAny{Int};
const Point2DF = Point2DAny{Float64};

iszero( p::Point2D ) = ( ( p.x == 0 )  &&  ( p.y == 0 ) )

function Base.show(io::IO, p::Point2D)
    print(io, "($(p.x), $(p.y))")
end

function  dY(p1::Point2D, p2::Point2D)
    dx = p1.x - p2.x
    dy = p1.y - p2.y
    return sqrt(dx^2 + dy^2)
end

function  norm( p::Point2DAny{T} )::Float64 where {T}
    return sqrt(p.x^2 + p.y^2)
end

function  normSQ( p::Point2D )
    return sqrt(p.x^2 + p.y^2)
end

# Addition and Subtraction (Vector-like)
Base.:+(a::Point2DAny{T}, b::Point2DAny{T}) where {T} = Point2DAny{T}(a.x + b.x, a.y + b.y)
Base.:-(a::Point2DAny{T}, b::Point2DAny{T}) where {T} = Point2DAny{T}(a.x - b.x, a.y - b.y)

# Multiplication by a constant (Scalar multiplication)
# We define both (Point * Constant) and (Constant * Point) for convenience
Base.:*(p::Point2DAny{T}, c::Real) where {T} = Point2DAny{T}(p.x * c, p.y * c)
Base.:*(c::Real, p::Point2DAny{T}) where {T} = p * c
Base.:/(p::Point2D, c::Real) = Point2D(p.xa / c, p.y / c)
dot(a::Point2DAny{T}, b::Point2DAny{T}) where {T} = a.x * b.x + a.y * b.y

#LinearAlgebra.norm(a::Point2D) = sqrt(a.x^2 + a.y^2)
Point2DAny{T}(t::Tuple{T,T}) where {T} =  Point2DAny{T}( convert( T, t[1] ), convert( T, t[2] ) )

Point2DAny{Float64}(p::Point2D) = Point2DAny{Float64}( p.x, p.y )

Base.convert(::Type{Point2D}, x::Tuple{Int,Int}) = Point2D(x)

vtrans( P, v ) =[p + v    for p ∈ P]

function is_lefteq_turn(a::Point2D, b::Point2D, c::Point2D)
    # The cross product of vectors (b-a) and (c-a):
    # (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
    return (b.x - a.x) * (c.y - a.y) >= (b.y - a.y) * (c.x - a.x)
end

function is_left_turn(a::Point2D, b::Point2D, c::Point2D)
    # The cross product of vectors (b-a) and (c-a):
    # (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
    return (b.x - a.x) * (c.y - a.y) > (b.y - a.y) * (c.x - a.x)
end

function is_right_turn(a::Point2D, b::Point2D, c::Point2D)
    # The cross product of vectors (b-a) and (c-a):
    # (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
    return (b.x - a.x) * (c.y - a.y) < (b.y - a.y) * (c.x - a.x)
end


const VPoint2D = Vector{Point2D};


"""
    distance_to_segment(p, a, b)

Computes the minimum Euclidean distance from point `p` to the line segment defined by endpoints `a` and `b`.
"""
function distance_to_segment(p::Point2D, a::Point2D, b::Point2D)
    ab = b - a
    ap = p - a
    
    # Calculate the squared length of the segment
    ab_len_sq = dot(ab, ab)
    
    # Handle the case where the segment is a single point (a == b)
    if ab_len_sq == 0.0
        return norm(p - a)
    end
    
    # Project vector AP onto AB to find the parameter t
    t = dot(ap, ab) / ab_len_sq
    
    # Clamp t to the segment range [0, 1]
    # t <= 0: closest point is A
    # t >= 1: closest point is B
    # 0 < t < 1: closest point is on the line between A and B
    t_clamped = clamp(t, 0.0, 1.0)
    
    # Find the closest point on the segment
    closest_point = Point2DF( a ) + (t_clamped * Point2DF(ab))
    
    return norm(Point2DF(p) - closest_point)
end

"""
    polygon_boundary_distance(poly, query_point)

Calculates the minimum distance from `query_point` to the boundary (edges) of `poly`.
`poly` is expected to be a Vector of Point2D.
"""
function polygon_boundary_distance(poly::Vector{Point2D}, query_point::Point2D)
    n = length(poly)
    if n < 2
        error("Polygon must have at least 2 vertices")
    end

    min_dist = Inf

    for i in 1:n
        # Get current vertex and the next vertex (wrapping around to the start)
        p1 = poly[i]
        p2 = poly[mod1(i + 1, n)]
        
        d = distance_to_segment(query_point, p1, p2)
        
        if d < min_dist
            min_dist = d
        end
    end

    return min_dist
end

"""
    is_point_in_polygon(poly, p)

Returns `true` if point `p` is inside the polygon `poly`, and `false` otherwise.
Uses the Ray Casting (Even-Odd) algorithm.
"""
function is_point_in_polygon(poly::Vector{Point2D}, p::Point2D)
    n = length(poly)
    inside = false
    
    # Iterate through each edge of the polygon
    j = n # The last vertex connects to the first
    for i in 1:n
        # Vertices of the current edge
        v1 = poly[i]
        v2 = poly[j]
        
        # Check if the ray intersects the edge:
        # 1. The point's Y coordinate must be between the edge's Y coordinates.
        #    (v1.y > p.y) != (v2.y > p.y) checks this.
        # 2. The intersection of the edge with the horizontal line y = p.y
        #    must lie to the right of point p.
        if ((v1.y > p.y) != (v2.y > p.y)) &&
           (p.x < (v2.x - v1.x) * (p.y - v1.y) / (v2.y - v1.y) + v1.x)
            
            # Toggle the state
            inside = !inside
        end
        
        j = i # Move to the next edge
    end
    
    return inside
end

# Interprent point  (i,j) as the fraction i / j.
function fary(vec::VPoint2D, u::Point2D, v::Point2D, max_d::Int)
    if (v.x > max_d)
        return
    end

    mid = Point2D(u.x + v.x, u.y + v.y)
    if (mid.x > max_d)
        return
    end
    fary(vec, u, mid, max_d)
    push!(vec, mid)
    fary(vec, mid, v, max_d)
end


function double_triangle_area(a::Point2D, b::Point2D, c::Point2D)::Int
    return abs(a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y))
end

# Triangle area alwasy return multiply of 1/2
function triangle_area(a::Point2D, b::Point2D, c::Point2D)::Float64
    return double_triangle_area(a, b, c) / 2.0
end

function grid_points_inside_edge(u::Point2D, v::Point2D)
    if  ( u == v )
        return  0
    end
    e = u - v
    t = gcd(e.x, e.y)
    @assert( t >= 1 )
    return t - 1;
end


"""
Computes the convex hull of a set of 2D points.
Returns a Vector{Point2D} representing the vertices of the hull in counter-clockwise order.
"""
function convex_hull(points::Vector{Point2D})
    # 1. Handle edge cases (0, 1, or 2 points are their own hull)
    n = length(points)
    if n <= 2
        return points
    end

    # 2. Sort points lexicographically (by x, then by y)
    # We create a copy to avoid modifying the input vector
    sorted_points = sort(points, by = p -> (p.x, p.y))

    # 3. Build the Lower Hull
    lower = Point2D[]
    for p in sorted_points
        # While the last 3 points do NOT make a left turn (i.e., they make a right turn or are straight),
        # remove the middle point because it cannot be on the hull.
        while length(lower) >= 2 && !is_left_turn(lower[end-1], lower[end], p)
            pop!(lower)
        end
        push!(lower, p)
    end

    # 4. Build the Upper Hull
    upper = Point2D[]
    # Iterate in reverse order for the upper hull
    for p in reverse(sorted_points)
        while length(upper) >= 2 && !is_left_turn(upper[end-1], upper[end], p)
            pop!(upper)
        end
        push!(upper, p)
    end

    # 5. Concatenate Hulls
    # The last point of 'lower' is the first point of 'upper' (and vice versa),
    # so we remove the duplicate last points before concatenating.
    pop!(lower)
    pop!(upper)

    return vcat(lower, upper)
end

function  get_y_range( x, r )
    y_max = -1
    for y ∈ 0:L
        if  norm( (x, y) ) ≤ r
            y_max = y
        else
            break
        end
    end
    return  -y_max, y_max
end

function  ch_disk_origin( k )
    r = sqrt( k / pi ) + 2

    V = VPoint2D()
    L = ceil( r ) + 1
    for x ∈ -L:L
        for y ∈ -L:L
            p = Point2D( x, y )
            if  norm( p ) > r
                continue
            end
            push!( V, p )
        end
    end
    @assert( sizeof( V ) > k )
    sort!(V, by = p -> (p.x^2 + p.y^2))
    resize!( V, k )
    println( V )
    sort!( V, by = p -> (p.x, p.y) )
    out  = VPoint2D();
    for  i ∈ eachindex(V)
        f_push = false
        if  ( i == 1  ||  i == length( V ) )
            push!( out, V[i ] )
            continue
        end
        if  V[i].x == V[i - 1].x  &&   V[i].x == V[i + 1].x
            continue            
        end
        push!( out, V[i ] )
    end
    ch = convex_hull( out )
    min_y = minimum(p -> p.y, ch )
    max_x = maximum(p.x for p in ch if p.y == min_y)

    ch_m = [(p - Point2D( max_x, min_y))  for p ∈ ch]

    return  ch_m
end


function generate_primitive_vectors(max_d::Int)
    vec = VPoint2D()
    #    push!( vec, (1, 0) )
    tvec = VPoint2D()
    fary(tvec, Point2D(1, 0), Point2D(1, 1), max_d)

    push!(vec, (1, 0))
    append!(vec, tvec)
    push!(vec, (1, 1))
    for v ∈ Iterators.reverse(tvec)
        push!(vec, (v.y, v.x))
    end
    ℓ = length(vec)
    for i ∈ 1:ℓ
        p = vec[i]
        push!(vec, (- p.y, p.x))
    end
    ℓ = length(vec)
    for i ∈ 1:ℓ
        p = vec[i]
        push!(vec, (- p.x, -p.y))
    end
    return vec
end

####################################################################################


"""
    DPStateKey(loc_prev::Point2D, loc::Point2D, n_g::Int)

Key for the Dynamic Programming state dictionary. Represents the configuration of the boundary being built.
- `loc_prev`: The previous vertex, used for angle and turn checking.
- `loc`: The current vertex.
- `n_g`: The number of internal grid points enclosed by the polygon so far.
"""
struct DPStateKey
    loc_prev::Point2D
    loc::Point2D
    n_g::Int           # Number of grid points covered so far
end

"""
    DPStateValue(perimeter_so_far::Float64, prev_key::DPStateKey, dir_index::Int)

Value stored in the Dynamic Programming state dictionary.
- `perimeter_so_far`: The total length of the boundary constructed from the origin up to `loc`.
- `prev_key`: Back-pointer to the preceding state (used to extract the optimal polygon).
- `dir_index`: The index of the incoming direction vector to `loc`, used to bound the outgoing angles.
"""
struct DPStateValue
    perimeter_so_far::Float64
    prev_key::DPStateKey
    dir_index::Int
end


import Base: hash, ==, isless
function Base.isless(a::Point2D, b::Point2D)
    return (a.x < b.x) || ((a.x == b.x) && (a.y < b.y))
end
function Base.hash(p::Point2D, h::UInt)
    hash(p.y, hash(p.x, h))
end
function ==(a::Point2D, b::Point2D)
    a.x == b.x && a.y == b.y
end


function Base.isless(a::DPStateKey, b::DPStateKey)
    return a.n_g < b.n_g || (a.n_g == b.n_g && (a.loc < b.loc || (a.loc == b.loc && a.loc_prev < b.loc_prev)))
end
function Base.hash(k::DPStateKey, h::UInt)
    hash(k.n_g, hash(k.loc, hash(k.loc_prev, h)))
end
function ==(a::DPStateKey, b::DPStateKey)
    a.loc_prev == b.loc_prev && a.loc == b.loc && a.n_g == b.n_g
end

function count_distinct(a, b, c)
    if a == b == c
        return 1
    elseif a == b || a == c || b == c
        return 2
    else
        return 3
    end
end



function is_colinear(a::Point2D, b::Point2D, c::Point2D)
    # The cross product of vectors (b-a) and (c-a):
    # (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
    return (b.x - a.x) * (c.y - a.y) == (b.y - a.y) * (c.x - a.x)
end

const DictStates = Dict{DPStateKey,DPStateValue}

function triangle_count_new_points(a::Point2D, b::Point2D, c::Point2D)
    area2 = double_triangle_area(a, b, c)
    ab_g_n = grid_points_inside_edge(a, b)
    bc_g_n = grid_points_inside_edge(b, c)
    ac_g_n = grid_points_inside_edge(a, c)

    @assert(ab_g_n >= 0)
    @assert(bc_g_n >= 0)
    @assert(ac_g_n >= 0)
    
    boundary_n = ab_g_n + bc_g_n + ac_g_n + count_distinct(a, b, c)

    if area2 > 0 
        tri_i_new = div(area2 - boundary_n + 2, 2)
    else
        tri_i_new = 0
    end

    tri_b_new = 0
    @assert(a != c)
    if a == b
        tri_b_new = ac_g_n + 1
    elseif area2 == 0
        tri_b_new = bc_g_n + 1
    else
        tri_b_new = ac_g_n + bc_g_n + 1
    end 

    return tri_i_new, tri_b_new
end


"""
    comp_next_conf(cfg::DPStateKey, val::DPStateValue, v::Point2D, k, max_d, opt_perim)

Computes the next dynamic programming state resulting from stepping by vector `v` 
from the current configuration `cfg`. 

Returns `(f_valid, next_cfg, new_perim)`:
- `f_valid`: Boolean indicating if the step is valid (e.g. valid turn, within bounds, does not exceed point constraints).
- `next_cfg`: The new `DPStateKey` containing the updated locus and point count.
- `new_perim`: The updated perimeter length including the length of the new segment.
"""
function   comp_next_conf(cfg::DPStateKey, val::DPStateValue, v::Point2D, k, max_d, opt_perim )
    a = cfg.loc_prev
    b = cfg.loc
    c = cfg.loc + v

    if  ( !is_lefteq_turn(a,b,c)  ||  ( c.y < 0 )  ||  ( a == c )  || iszero( c ) 
          ||  c.x > max_d  ||  c.y > max_d  ||  -c.x > max_d  ||  -c.y > max_d
          ||  val.perimeter_so_far > opt_perim
          ||  is_right_turn( b, c, Point2D( 0, 0 ) ) )
        return  false, cfg, -1.0 
    end

    # It is not allowed to back to the origin...
    if  is_colinear( b, c, Point2D(0 , 0) )
        if  normSQ(c) <  normSQ(b)
            return  false, cfg, -1.0 
        end
    end
   
    @assert( b != c )
    @assert( v.x != 0  ||  v.y != 0 )

    
    # New triangle: ◬( origin, b, c )
    origin = Point2D( 0, 0 )
    tri_i_new, tri_b_new = triangle_count_new_points( origin, b, c )

    @assert( tri_i_new >= 0  &&  tri_b_new >= 0 ) 
    n_g = cfg.n_g + tri_i_new + tri_b_new 
    if  ( n_g > 0  &&  iszero( c ) )
        return  false, cfg, -1.0 
    end
    if  n_g > k
        return  false, cfg, -1.0 
    end

    new_cfg = DPStateKey( b, c, n_g )

    return  true, new_cfg, val.perimeter_so_far + dY( b, c )
end



function extract_solution( d, cfg )
    out = VPoint2D()

    curr = cfg
    while  true
        val = d[ curr ]
        #println( "-- curr", curr, "  ", val.perimeter_so_far )
        push!( out, curr.loc)
        if  ( val.prev_key == curr )
            break
        end
        curr = val.prev_key
 #       ( curr.loc.x != 0  ||  curr.loc.y != 0 
#             ||curr.loc_prev.x != 0  ||  curr.loc_prev.y != 0 )
#        if  ()
    end
    return  out
end

function cr_string_pdf(cr, w, h, text_content::String)
    
    # 2. Draw a background or shape (Optional)
    set_source_rgb(cr, 1.0, 1.0, 1.0 ) # Light gray
    rectangle(cr, 0, 0, w, h )
    fill(cr)

    # 3. Setup Font Properties
    # select_font_face(context, family, slant, weight)
    select_font_face(cr, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
    set_font_size(cr, 8.0)
    fh = 16
    
    lines = split(text_content, "\n")
    current_y = fh * 1.5
    
    set_source_rgb(cr, 0.0, 0.0, 0.0) # Black text
    for line in lines
        move_to(cr, 10, current_y)
        show_text(cr, line)
        
        # Move down for the next line
        current_y += fh
    end
            # 4. Position and Render Text

end

function angle_between(p1::Point2D, p2::Point2D)
    dot_prod = p1.x * p2.x + p1.y * p2.y
    det_prod = p1.x * p2.y - p1.y * p2.x # 2D "cross product"
    
    # atan(y, x) in Julia computes atan2
    return atan(det_prod, dot_prod)
end

"""
    comp_stop_indexes(V, max_angle)

Precomputes the permissible range of outgoing direction vectors for each incoming direction vector in `V`, 
such that the angle between consecutive segments does not exceed `max_angle`.
Returns a vector of indices indicating the stop index for the inner loop.
"""
function comp_stop_indexes( V, max_angle )
    stops = Vector{Int}()
    for i ∈ eachindex( V )
        dir_end = length( V )
        prev_dir = V[ i ]
        f_break = false 
        for dir_i ∈ i:dir_end
            v = V[ dir_i ]
            α = angle_between( prev_dir, v )
            if  α > max_angle 
                push!( stops, dir_i - 1 )
                f_break = true
                break
            end
        end
        if  ! f_break
            push!( stops, dir_end )
        end
    end
    @assert( length( stops ) == length( V ))
    return  stops
end

function is_store(d_all, next_cfg, new_perim, opt_perim)
    if new_perim > opt_perim
        return false
    end

    val = get(d_all, next_cfg, nothing)
    if val === nothing
        return true
    end

    # Is the stored solution better?
    if val.perimeter_so_far < new_perim
        return false
    end
    return true
end

"""
    minimize_perimeter_dp(k, good, bad)

Executes the main Dynamic Programming algorithm to find the minimum-perimeter polygon
enclosing exactly `k` grid points, starting near the origin, while strictly enclosing a `good` set
of points and strictly excluding a `bad` set.

Returns `(sol, ub_circle)`:
- `sol`: A vector of `Point2D` vertices of the optimal polygon.
- `ub_circle`: The theoretical upper bound perimeter based on a perfect circle.
"""
function minimize_perimeter_dp( k, good, bad )
    max_d = ceil( Int, sqrt(k) ) + 1

    max_angle = 2 / k^(1/8)
    ub_circle = 2.0 * sqrt( pi * k )
    sq_perim = 4 * max_d + 6;
    opt_perim = min( ub_circle, sq_perim )
   
    #println( "Upper bound on perim:", opt_perim )
    best_sol = start_key = DPStateKey((0, 0), (0, 0), 1)

    pq = BinaryMinHeap{Tuple{Int, DPStateKey}}()
    push!( pq, (1, start_key ) )
    
    d_all = DictStates()
    start_val = DPStateValue(0.0, start_key, 0 )
    d_all[ start_key ] = start_val

    max_edge_l = round( Int, k^(1/3) ) + 1
    println( max_edge_l )

    V = generate_primitive_vectors(max_edge_l)
    stops = comp_stop_indexes( V, max_angle )
    
    handled = Set{DPStateKey}()
    
    while !isempty(pq)
        nq, cfg = pop!( pq )
        val = d_all[ cfg ]
        if  cfg ∈ handled
            continue
        end
        push!( handled, cfg )
            
        dir_start = max( val.dir_index,  1 )
        dir_end = stops[ dir_start ];
        prev_dir = V[ dir_start ]
        for dir_i ∈ dir_start:dir_end
            v = V[ dir_i ]
            if  cfg.loc.y == 0  &&  v.y <= 0
                continue
            end
            p_next = cfg.loc + v
            if  p_next ∈ bad
                continue
            end
            
            f_valid, next_cfg, new_perim = comp_next_conf(cfg, val, v, k, max_d, opt_perim )
            if ! f_valid
                continue
            end
            
            if ! is_store( d_all, next_cfg, new_perim, opt_perim )
                continue;
            end
            
            d_all[next_cfg] = DPStateValue(new_perim, cfg, dir_i )
            push!( pq, (next_cfg.n_g, next_cfg ) )
 
            if ( next_cfg.n_g == k )
                total_perim = new_perim + norm( next_cfg.loc ) # Because start locaion is (0,0)
                if total_perim < opt_perim
                    opt_perim = total_perim
                    best_sol = next_cfg
                end
            end
        end
    end

    sol = extract_solution( d_all, best_sol )
    return  sol, ub_circle;
end

function print_fractions(vec)
    for p ∈ vec
        println(p.x, ",", p.y, "  : ", rad2deg(atan(p.y, p.x)))
    end
end

function   draw_polygon( cr, poly, rgba )
    pt = first( poly )
    move_to(cr, pt... )
    
    for i in 2:length( poly )
        line_to(cr, poly[ i ]... )
    end
    close_path(cr) # Connect last point back to first
    
    # Stroke the outline
    stroke_preserve(cr)
    
    # Fill with semi-transparent color
    set_source_rgba(cr, rgba... )
    
    fill(cr)
end

function  bound( polys, expand )
    min_x = max_x = 0
    min_y = max_y = 0
    f_init = false
    
    for  P ∈ polys
        xs = [p.x for p in P]
        ys = [p.y for p in P]
        _min_x, _max_x = minimum(xs) - 1, maximum(xs) + 1
        _min_y, _max_y = minimum(ys) - 1, maximum(ys) + 1
        if  !f_init
            min_x = _min_x
            min_x = _max_x
            min_y = _min_y
            min_y = _max_y
            f_init = true
        end
        min_x = min( _min_x, min_x )
        min_y = min( _min_y, min_y )
        max_x = max( _max_x, max_x )
        max_y = max( _max_y, max_y )
    end
    return  min_x - expand, max_x + expand, min_y - expand, max_y + expand
end

function draw_polygon_with_grid( poly::VPoint2D, poly_circ::VPoint2D, k, ub_circle,
    bad:: Set{Point2D},
    good::Set{Point2D}
)
    filename = @sprintf( "output/%d.pdf", k )

    min_x, max_x, min_y, max_y = bound( [poly, poly_circ ], 1 )

    # 2. Setup Translation and Padding
    # We add a small margin (e.g., 50 units) so the points aren't on the very edge
    margin = 50
    scale = round( Int, max( 100/k, 10 ) )
    width = scale * (max_x - min_x) + 2 * margin
    height = scale * (max_y - min_y) + 2 * margin
    
    # Translation function to make all coordinates positive
    translate(p::Point2D) = (scale * (p.x - min_x) + margin,
                             scale * (p.y - min_y) + margin)
    
    #println( "W:", width, " h: ", height )
    # 3. Initialize Cairo PDF Surface
    surface = CairoPDFSurface(filename, width, height)
    cr = CairoContext(surface)
    
    
    # --- Draw the Polygon ---
    set_line_width(cr, 2.0)
    set_source_rgb(cr, 0.2, 0.5, 0.8) # "Julia Blue"
    
    trans( P ) =[translate( p )   for p ∈ P]
    draw_polygon( cr, trans( poly ), (0.2, 0.5, 0.8, 1.0 ) )

    set_line_width(cr, 1.0)
    set_source_rgb(cr, 0.9, 0.5, 0.8) # "Julia Blue"
    draw_polygon( cr, trans( poly_circ ), (0.9, 0.5, 0.2, 0.1 ) )
    
    # --- Draw the Grid Points ---
    for x in min_x:max_x
        for y in min_y:max_y
            tx, ty = translate(Point2D(x, y))
            if  ( x == 0  &&  y == 0 )
                set_source_rgb(cr, 1.0, 1.0, 0.0) 
            elseif  Point2D( x, y) ∈ bad
                set_source_rgb(cr, 0.8, 0.0, 0.0) # Light gray for grid                
            elseif  Point2D( x, y) ∈ good
                set_source_rgb(cr, 0.0, 0.0, 1.0) # Light gray for grid                
            else
                set_source_rgb(cr, 0.2, 0.2, 0.2) # Light gray for grid
            end
            
            arc(cr, tx, ty, 1.5, 0, 2π) # Draw a small dot
            fill(cr)
        end
    end

    show_page( cr )

    str = @sprintf( "k = %d\nCircle ub on perimeter: %g\n", k, ub_circle )
    str_b = @sprintf( "solution parimeter: %g\n# vertices: %d\n",
        compute_perimeter( poly ), length( poly ) )
    cr_string_pdf(cr, width, height, str * str_b)
    
    # 4. Cleanup
    finish(surface)
    println("Polygon saved to $filename")
end

function compute_perimeter( sol )
    ℓ = dY( first( sol), last( sol ) )
    for i ∈ 2:length( sol )
        ℓ += dY( sol[ i-1], sol[ i ])
    end

    return  ℓ     
end

function show_perimeter( sol )
    ℓ = 0
    for i ∈ length( sol ):-1:2
        d = dY( sol[ i-1], sol[ i ])
        ℓ += d
        println( sol[ i - 1 ], " -- ", sol[ i ]," : ", ℓ )
    end
    ℓ += dY( first( sol), last( sol ) )
    println( first( sol), " -- ", last(sol), " : ", ℓ )

    return  ℓ     
end


"""
    compute_good_bad_sets(ch_m, ℓ)

Classifies points in a bounding box around a base shape `ch_m` into `good` (points that must be enclosed)
and `bad` (points that must be excluded) based on the polygon boundary distance `ℓ`.

Returns `(good, bad)` sets of `Point2D`.
"""
function  compute_good_bad_sets( ch_m, ℓ )
    min_x, max_x, min_y, max_y = bound( [ch_m], ℓ + 2 );
    
    bad = Set{Point2D}();
    good = Set{Point2D}();
    for  y ∈ min_y:max_y
        if y < 0
            continue
        end
        for  x ∈ min_x:max_x
            p = Point2D( x, y )
            f_in =  is_point_in_polygon(ch_m, p)
            d = polygon_boundary_distance( ch_m, p )

            if  d > ℓ
                push!( bad, p )
                continue
            end

            if  f_in  ||  d ≤ ℓ 
                push!( good, p)
            end
        end
    end
    return  good, bad
end


"""
    @main(args)

Command-line entry point. Takes `k` (the number of points to enclose) as an argument.
Example: `julia --project src/k_perimeter.jl 10`
"""
function (@main)(args)
    k = 6
    if isempty(args)
        println("Usage: k_perimeter [k]")
        exit(1)
    end
    k = -1
    try
        k = parse(Int, args[1])
    catch e
        if e isa ArgumentError
            println("Error: '$(args[1])' is not a valid integer.")
        else
            rethrow(e)
        end
    end

    ch_z = ch_disk_origin( k )
    
    Δ = ceil( k^(1/4)/4 )
    if  ( k > 400)
        Δ = 0
    end
    ch_m = vtrans( ch_z, Point2D( -Δ , 0  ) )

    ℓ = max( round( Int, k^(1/4)/3 ), 3 )

    good, bad =  compute_good_bad_sets( ch_m, ℓ )
    
    sol, ub_circle = minimize_perimeter_dp( k, good, bad )

    println( "k: ", k )
    perimeter = compute_perimeter( sol )
    println( "Perimeter        : ", perimeter )
    println( "circle perimeter : ", compute_perimeter( ch_m ) )
    println( "Naive perimeter  : ", ub_circle )

    draw_polygon_with_grid( sol, ch_m, k, ub_circle, bad, good )

    
#    println( sol )
    
    return 0
end
