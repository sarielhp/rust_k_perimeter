#! /bin/env julial

"""
    medial_axis.jl
    Computes the medial axis (straight skeleton) and furthest-neighbor Voronoi diagram
    of a convex polygon, finding the optimal point for the minimum width annulus.
"""

using LinearAlgebra
using Cairo

# Ensure output directory exists
output_dir = joinpath("output", "annulus")
if !ispath(output_dir)
    mkpath(output_dir)
end

struct Point
    x::Float64
    y::Float64
end

# Basic point operations
Base.:+(p1::Point, p2::Point) = Point(p1.x + p2.x, p1.y + p2.y)
Base.:-(p1::Point, p2::Point) = Point(p1.x - p2.x, p1.y - p2.y)
Base.:*(p::Point, s::Real) = Point(p.x * s, p.y * s)
Base.:*(s::Real, p::Point) = Point(p.x * s, p.y * s)
Base.:/(p::Point, s::Real) = Point(p.x / s, p.y / s)

dist(p1::Point, p2::Point) = sqrt((p1.x - p2.x)^2 + (p1.y - p2.y)^2)

struct Line
    a::Float64
    b::Float64
    c::Float64
end

function get_line(p1::Point, p2::Point)
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    len = sqrt(dx^2 + dy^2)
    a = -dy / len
    b = dx / len
    c = -(a * p1.x + b * p1.y)
    return Line(a, b, c)
end

function intersect(l1::Line, l2::Line, d::Float64)
    det = l1.a * l2.b - l1.b * l2.a
    if abs(det) < 1e-12
        return nothing
    end
    x = (l2.b * (d - l1.c) - l1.b * (d - l2.c)) / det
    y = (l1.a * (d - l2.c) - l2.a * (d - l1.c)) / det
    return Point(x, y)
end

mutable struct EdgeNode
    line::Line
    prev::Union{Nothing,EdgeNode}
    next::Union{Nothing,EdgeNode}
    id::Int
    EdgeNode(line, id) = new(line, nothing, nothing, id)
end

function get_collapse_time(e::EdgeNode)
    l1, l2, l3 = e.prev.line, e.line, e.next.line
    det12 = l1.a * l2.b - l1.b * l2.a
    if abs(det12) < 1e-12
        return Inf
    end
    x0_12, vx_12 = (l2.b * (-l1.c) - l1.b * (-l2.c)) / det12, (l2.b - l1.b) / det12
    det23 = l2.a * l3.b - l2.b * l3.a
    if abs(det23) < 1e-12
        return Inf
    end
    x0_23, vx_23 = (l3.b * (-l2.c) - l2.b * (-l3.c)) / det23, (l3.b - l2.b) / det23
    dvx, dx0 = vx_12 - vx_23, x0_23 - x0_12
    if abs(dvx) < 1e-12
        y0_12, vy_12 = (l1.a * (-l2.c) - l2.a * (-l1.c)) / det12, (l1.a - l2.a) / det12
        y0_23, vy_23 = (l2.a * (-l3.c) - l3.a * (-l2.c)) / det23, (l2.a - l3.a) / det23
        dvy, dy0 = vy_12 - vy_23, y0_23 - y0_12
        if abs(dvy) < 1e-12
            return Inf
        end
        d = dy0 / dvy
    else
        d = dx0 / dvx
    end
    return d > -1e-8 ? d : Inf
end

function read_polygon(filename)
    points = Point[]
    for line in eachline(filename)
        line = strip(line);
        if isempty(line) || startswith(line, "#")
            continue
        end
        parts = split(line)
        if length(parts) >= 2
            push!(points, Point(parse(Float64, parts[1]), parse(Float64, parts[2])))
        end
    end
    if length(points) > 1
        filtered = [points[1]]
        for i = 2:length(points)
            if dist(points[i], points[i-1]) > 1e-9
                push!(filtered, points[i])
            end
        end
        if length(filtered) > 1 && dist(filtered[end], filtered[1]) < 1e-9
            pop!(filtered)
        end
        points = filtered
    end
    
    # Ensure CCW orientation for is_inside_lenient
    if length(points) >= 3
        area = 0.0
        for i = 1:length(points)
            p1, p2 = points[i], points[mod1(i+1, length(points))]
            area += (p1.x * p2.y - p2.x * p1.y)
        end
        if area < 0
            reverse!(points)
        end
    end

    # Find integer points inside the original polygon
    integer_points = Set{Point}()
    if length(points) >= 3
        min_ox = floor(Int, minimum(p.x for p in points))
        max_ox = ceil(Int, maximum(p.x for p in points))
        min_oy = floor(Int, minimum(p.y for p in points))
        max_oy = ceil(Int, maximum(p.y for p in points))
        for ix = min_ox:max_ox
            for iy = min_oy:max_oy
                ip = Point(Float64(ix), Float64(iy))
                if is_inside_lenient(ip, points)
                    push!(integer_points, ip)
                end
            end
        end
        # Explicitly add vertices if they are integers (with tolerance)
        for p in points
            if abs(p.x - round(p.x)) < 1e-8 && abs(p.y - round(p.y)) < 1e-8
                push!(integer_points, Point(round(p.x), round(p.y)))
            end
        end
    end
    int_pts_list = collect(integer_points)

    scale_factor = 1.0
    # Translate so all coordinates are non-negative
    if !isempty(points)
        min_x = minimum(p.x for p in points)
        min_y = minimum(p.y for p in points)
        points = [Point(p.x - min_x, p.y - min_y) for p in points]
        int_pts_list = [Point(p.x - min_x, p.y - min_y) for p in int_pts_list]
        
        # Scale so coordinates are in [0, 1000]^2 and max(max_x, max_y) == 1000
        max_val = max(maximum(p.x for p in points), maximum(p.y for p in points))
        if max_val > 1e-12
            scale_factor = 1000.0 / max_val
            points = [Point(p.x * scale_factor, p.y * scale_factor) for p in points]
            int_pts_list = [Point(p.x * scale_factor, p.y * scale_factor) for p in int_pts_list]
        end
    end
    
    points = remove_redundant_vertices(points)
    return points, scale_factor, int_pts_list
end

function remove_redundant_vertices(points)
    if length(points) < 3
        return points
    end
    changed = true
    while changed && length(points) >= 3
        changed = false;
        i = 1
        while i <= length(points) && length(points) >= 3
            p_prev, p_curr, p_next = points[mod1(i-1, length(points))],
            points[i],
            points[mod1(i+1, length(points))]
            vx, vy, wx, wy = p_next.x - p_prev.x,
            p_next.y - p_prev.y,
            p_curr.x - p_prev.x,
            p_curr.y - p_prev.y
            cp = vx * wy - vy * wx
            if abs(cp) < 1e-10
                dotp, sqlen = wx * vx + wy * vy, vx * vx + vy * vy
                if dotp >= -1e-10 && dotp <= sqlen + 1e-10
                    deleteat!(points, i);
                    changed = true;
                    continue
                end
            end
            i += 1
        end
    end
    return points
end

function run_skeleton(points)
    n = length(points);
    if n < 3
        return Tuple{Point,Point}[]
    end
    active_edges = EdgeNode[]
    for i = 1:n
        p1, p2 = points[i], points[mod1(i+1, n)]
        push!(active_edges, EdgeNode(get_line(p1, p2), i))
    end
    for i = 1:n
        active_edges[i].prev = active_edges[mod1(i-1, n)];
        active_edges[i].next = active_edges[mod1(i+1, n)]
    end
    current_vertices =
        [intersect(active_edges[i].prev.line, active_edges[i].line, 0.0) for i = 1:n]
    for i = 1:n
        if current_vertices[i] === nothing
            current_vertices[i] = points[i]
        end
    end
    medial_segments = Tuple{Point,Point}[]
    while n > 3
        min_d, idx = Inf, -1
        for i = 1:n
            d = get_collapse_time(active_edges[i])
            if d < min_d
                min_d, idx = d, i
            end
        end
        if idx == -1 || min_d == Inf
            break
        end
        e = active_edges[idx]
        p_new = intersect(e.prev.line, e.line, min_d)
        if p_new === nothing
            p_new = intersect(e.line, e.next.line, min_d)
        end
        if p_new === nothing
            break
        end
        push!(medial_segments, (current_vertices[idx], p_new))
        next_idx = mod1(idx + 1, n)
        push!(medial_segments, (current_vertices[next_idx], p_new))
        new_active_edges, new_vertices = EdgeNode[], Point[]
        for i = 1:n
            if i == idx
                continue
            end
            push!(new_active_edges, active_edges[i])
            push!(new_vertices, i == next_idx ? p_new : current_vertices[i])
        end
        active_edges, current_vertices, n =
            new_active_edges, new_vertices, length(new_active_edges)
        for i = 1:n
            active_edges[i].prev = active_edges[mod1(i-1, n)];
            active_edges[i].next = active_edges[mod1(i+1, n)]
        end
    end
    if n == 3
        d_f = get_collapse_time(active_edges[2])
        p_f = intersect(active_edges[1].line, active_edges[2].line, d_f)
        if p_f === nothing
            p_f = intersect(active_edges[2].line, active_edges[3].line, d_f)
        end
        if p_f !== nothing
            for i = 1:3
                push!(medial_segments, (current_vertices[i], p_f))
            end
        end
    end
    return medial_segments
end

function clip_segment_to_polygon(p1::Point, p2::Point, poly::Vector{Point})
    t_min, t_max, d = 0.0, 1.0, p2 - p1
    for i = 1:length(poly)
        v1, v2 = poly[i], poly[mod1(i+1, length(poly))]
        e = v2 - v1;
        n = Point(-e.y, e.x)
        num, den = n.x * (v1.x - p1.x) + n.y * (v1.y - p1.y), n.x * d.x + n.y * d.y
        if abs(den) < 1e-13
            if n.x * (p1.x - v1.x) + n.y * (p1.y - v1.y) < -1e-9
                return nothing
            end
        else
            t = num / den
            if den > 0
                t_min = max(t_min, t)
            else
                t_max = min(t_max, t)
            end
        end
    end
    return t_min <= t_max + 1e-11 ? (p1 + t_min * d, p1 + t_max * d) : nothing
end

function intersect_halfplane(poly::Vector{Point}, N::Point, c::Float64)
    out_poly = Point[]
    for i = 1:length(poly)
        p1, p2 = poly[i], poly[mod1(i+1, length(poly))]
        d1, d2 = p1.x * N.x + p1.y * N.y - c, p2.x * N.x + p2.y * N.y - c
        if d1 <= 1e-9
            push!(out_poly, p1)
            if d2 > 1e-9
                t = d1 / (d1 - d2);
                push!(out_poly, Point(p1.x + t*(p2.x-p1.x), p1.y + t*(p2.y-p1.y)))
            end
        elseif d2 <= 1e-9
            t = d1 / (d1 - d2);
            push!(out_poly, Point(p1.x + t*(p2.x-p1.x), p1.y + t*(p2.y-p1.y)))
        end
    end
    return out_poly
end

function is_inside_lenient(p::Point, poly::Vector{Point})
    for i = 1:length(poly)
        v1, v2 = poly[i], poly[mod1(i+1, length(poly))]
        len = dist(v1, v2)
        if len > 1e-9
            cp = (v2.x - v1.x) * (p.y - v1.y) - (v2.y - v1.y) * (p.x - v1.x)
            if cp / len < -1e-7
                return false
            end
        end
    end
    return true
end

function compute_fnvd_raw(points::Vector{Point})
    n, fnvd_segments = length(points), Tuple{Point,Point}[]
    for i = 1:n
        cell = copy(points)
        for j = 1:n
            if i == j
                continue
            end
            N = Point(points[i].x - points[j].x, points[i].y - points[j].y)
            c = 0.5 * (points[i].x^2 + points[i].y^2 - points[j].x^2 - points[j].y^2)
            cell = intersect_halfplane(cell, Point(-N.x, -N.y), -c)
            if isempty(cell)
                break
            end
        end
        for k = 1:length(cell)
            push!(fnvd_segments, (cell[k], cell[mod1(k+1, length(cell))]))
        end
    end
    return fnvd_segments
end

function segment_intersection(p1::Point, p2::Point, p3::Point, p4::Point)
    dx1, dy1, dx2, dy2 = p2.x-p1.x, p2.y-p1.y, p4.x-p3.x, p4.y-p3.y
    det = dx1 * dy2 - dy1 * dx2;
    if abs(det) < 1e-10
        return nothing
    end
    t, u = ((p3.x-p1.x)*dy2 - (p3.y-p1.y)*dx2)/det, ((p3.x-p1.x)*dy1 - (p3.y-p1.y)*dx1)/det
    return (-1e-9 <= t <= 1.0+1e-9 && -1e-9 <= u <= 1.0+1e-9) ?
           Point(p1.x + t*dx1, p1.y + t*dy1) : nothing
end

function point_to_segment_dist(p::Point, v1::Point, v2::Point)
    dx, dy = v2.x - v1.x, v2.y - v1.y;
    len2 = dx^2 + dy^2
    if len2 < 1e-12
        return dist(p, v1)
    end
    t = clamp(((p.x - v1.x) * dx + (p.y - v1.y) * dy) / len2, 0.0, 1.0)
    return dist(p, Point(v1.x + t * dx, v1.y + t * dy))
end

function closest_point_on_segment(p::Point, v1::Point, v2::Point)
    dx, dy = v2.x - v1.x, v2.y - v1.y
    len2 = dx^2 + dy^2
    if len2 < 1e-12
        return v1
    end
    t = clamp(((p.x - v1.x) * dx + (p.y - v1.y) * dy) / len2, 0.0, 1.0)
    return Point(v1.x + t * dx, v1.y + t * dy)
end

function get_r_in(p::Point, poly::Vector{Point})
    min_r = Inf
    for i = 1:length(poly)
        min_r = min(min_r, point_to_segment_dist(p, poly[i], poly[mod1(i+1, length(poly))]))
    end
    return min_r
end

function get_r_out(p::Point, poly::Vector{Point})
    max_r = 0.0
    for v in poly
        max_r = max(max_r, dist(p, v))
    end
    return max_r
end

function draw_output(points, medial_raw, fnvd_raw, pdf_file, scale_factor, k, final_name, integer_points)
    all_segs, all_points = [medial_raw; fnvd_raw], Set{Point}()
    for (p1, p2) in all_segs
        if is_inside_lenient(p1, points)
            push!(all_points, p1)
        end
        if is_inside_lenient(p2, points)
            push!(all_points, p2)
        end
    end
    for i = 1:length(all_segs), j = (i+1):length(all_segs)
        p_int = segment_intersection(
            all_segs[i][1],
            all_segs[i][2],
            all_segs[j][1],
            all_segs[j][2],
        )
        if p_int !== nothing && is_inside_lenient(p_int, points)
            found = false
            for p in all_points
                if dist(p, p_int) < 1e-7
                    found = true;
                    break
                end
            end
            if !found
                push!(all_points, p_int)
            end
        end
    end
    inter_points = collect(all_points)
    best_p, min_diff, best_rin, best_rout = nothing, Inf, 0.0, 0.0
    for p in inter_points
        rin = get_r_in(p, points)
        rout = get_r_out(p, points);
        diff = rout - rin
        if diff < min_diff
            min_diff, best_p, best_rin, best_rout = diff, p, rin, rout
        end
    end
    
    normalized_width = 0.0
    if best_p !== nothing
        normalized_width = (best_rout - best_rin) / scale_factor
        txt_file = joinpath(output_dir, final_name * ".txt")
        open(txt_file, "w") do f
            println(f, "# annulus width : $normalized_width")
        end
    end

    medial_clipped = [
        c for
        c in [clip_segment_to_polygon(p1, p2, points) for (p1, p2) in medial_raw] if
        c !== nothing
    ]
    fnvd_clipped = [
        c for c in [clip_segment_to_polygon(p1, p2, points) for (p1, p2) in fnvd_raw] if
        c !== nothing
    ]
    xs, ys = [p.x for p in points], [p.y for p in points]
    min_x, max_x, min_y, max_y = minimum(xs), maximum(xs), minimum(ys), maximum(ys)
    w, h, p = max_x - min_x, max_y - min_y, 0.1 * max(max_x - min_x, max_y - min_y)
    sw, sh = w + 2*p, h + 2*p
    sf = max(sw, sh) > 1000 ? 1000.0 / max(sw, sh) : 1.0
    pw, maxw = sw * sf, sw * sf / 200.0
    polyw, segw, printr, bprintr = maxw/sf, (maxw*0.25)/sf, (maxw*0.4)/sf, (maxw*0.8)/sf
    c = CairoPDFSurface(pdf_file, sw * sf, sh * sf);
    cr = CairoContext(c)
    function setup(cr)
        set_source_rgb(cr, 1, 1, 1);
        paint(cr);
        save(cr);
        scale(cr, sf, sf);
        translate(cr, p-min_x, p+h+min_y);
        scale(cr, 1, -1)
    end
    
    function draw_p(cr, wd)
        set_source_rgb(cr, 0, 0, 0);
        set_line_width(cr, wd);
        move_to(cr, points[1].x, points[1].y);
        for i = 2:length(points)
            line_to(cr, points[i].x, points[i].y)
        end;
        close_path(cr);
        stroke(cr)
    end
    function draw_s(cr, sgs, r, g, b, wd)
        set_source_rgb(cr, r, g, b);
        set_line_width(cr, wd);
        for (p1, p2) in sgs
            move_to(cr, p1.x, p1.y);
            line_to(cr, p2.x, p2.y);
            stroke(cr)
        end
    end
    setup(cr);
    draw_p(cr, polyw);
    draw_s(cr, medial_clipped, 1, 0, 0, segw);
    restore(cr);
    show_page(cr)
    
    setup(cr);
    draw_p(cr, polyw);
    draw_s(cr, fnvd_clipped, 0, 0, 1, segw);
    restore(cr);
    # Draw text on second page
    set_source_rgb(cr, 0, 0, 0)
    select_font_face(cr, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
    set_font_size(cr, 20.0)
    move_to(cr, 50.0, 50.0)
    show_text(cr, "k : $k")
    move_to(cr, 50.0, 80.0)
    show_text(cr, "$normalized_width")
    show_page(cr)
    
    setup(cr);
    draw_p(cr, polyw);
    draw_s(cr, medial_clipped, 1, 0, 0, segw);
    draw_s(cr, fnvd_clipped, 0, 0, 1, segw);
    restore(cr);
    show_page(cr)
    
    setup(cr);
    draw_p(cr, polyw);
    set_source_rgb(cr, 0, 0.5, 0);
    for pt in inter_points
        arc(cr, pt.x, pt.y, printr, 0, 2π);
        fill(cr)
    end;
    restore(cr);
    show_page(cr)
    setup(cr)
    if best_p !== nothing
        set_source_rgb(cr, 1.0, 1.0, 0.878);
        arc(cr, best_p.x, best_p.y, best_rout, 0, 2π);
        fill(cr)
        stroke(cr)
        set_source_rgb(cr, 1.0, 1.0, 1.0);
        arc(cr, best_p.x, best_p.y, best_rin, 0, 2π);
        fill(cr)
        set_source_rgb(cr, 0.5, 0.5, 0.5);
        set_line_width(cr, (maxw/10.0)/sf);
        arc(cr, best_p.x, best_p.y, best_rin, 0, 2π);
        stroke(cr);
        arc(cr, best_p.x, best_p.y, best_rout, 0, 2π);
        stroke(cr)
        set_source_rgb(cr, 1, 0, 0);
        arc(cr, best_p.x, best_p.y, bprintr, 0, 2π);
        fill(cr)
        stroke(cr)
    end
    
    # Draw grid lines first (blue, thinner)
    if !isempty(integer_points)
        set_source_rgb(cr, 0, 0, 1);
        set_line_width(cr, (segw * 0.5) / 3.0);
        
        # Vertical lines
        x_vals = sort(collect(Set(p.x for p in integer_points)))
        for xv in x_vals
            pts = filter(p -> abs(p.x - xv) < 1e-9, integer_points)
            y_min = minimum(p.y for p in pts)
            y_max = maximum(p.y for p in pts)
            move_to(cr, xv, y_min)
            line_to(cr, xv, y_max)
            stroke(cr)
        end
        
        # Horizontal lines
        y_vals = sort(collect(Set(p.y for p in integer_points)))
        for yv in y_vals
            pts = filter(p -> abs(p.y - yv) < 1e-9, integer_points)
            x_min = minimum(p.x for p in pts)
            x_max = maximum(p.x for p in pts)
            move_to(cr, x_min, yv)
            line_to(cr, x_max, yv)
            stroke(cr)
        end
    end

    # Draw polygon outline
    set_source_rgb(cr, 0, 0, 0);
    set_line_width(cr, (maxw/10.0)/sf);
    move_to(cr, points[1].x, points[1].y);
    for i = 2:length(points)
        line_to(cr, points[i].x, points[i].y)
    end;
    close_path(cr);
    stroke(cr);

    # Draw vertices as red dots last
    set_source_rgb(cr, 1, 0, 0);
    for pt in points
        arc(cr, pt.x, pt.y, printr * 0.5, 0, 2π);
        fill(cr)
    end

    restore(cr);
    show_page(cr);
    finish(c)
    
    println(pdf_file)
    txt_file = joinpath(output_dir, final_name * ".txt")
    if isfile(txt_file)
        println(txt_file)
    end
end

function main()
    if length(ARGS) < 1
        println("Usage: julia medial_axis.jl <polygon_file> [output_base]");
        return
    end
    input_file = ARGS[1];
    base_name = replace(basename(input_file), r"\.[^.]+$" => "")
    final_name = length(ARGS) >= 2 ? ARGS[2] : "a_" * replace(base_name, "_poly" => "")
    pdf_file = joinpath(output_dir, final_name * ".pdf")
    points, scale_factor, integer_points = read_polygon(input_file)
    if length(points) < 3
        println("Polygon must have at least 3 vertices.");
        return
    end
    k = length(points)
    medial_raw = run_skeleton(points)
    fnvd_raw = compute_fnvd_raw(points)
    draw_output(points, medial_raw, fnvd_raw, pdf_file, scale_factor, k, final_name, integer_points)
end

main()
