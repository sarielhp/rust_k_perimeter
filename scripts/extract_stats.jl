#! /bin/env julial

using DataFrames, CSV

# --- Configuration ---
const DEFAULT_INPUT_DIR = "output/summary"
const DEFAULT_OUTPUT_FILE = "misc/extracted_data.csv"
const TRIGGER_LINE = "# k:"

"""
Smart parser that preserves large integers as Int128 
and only uses Float64 when a decimal is present.
"""
function smart_parse(s::AbstractString)
    # Try parsing as a large Integer first to keep it explicit
    i = tryparse(Int128, s)
    !isnothing(i) && return i
    
    # Fallback to Float64 if it contains a decimal or scientific notation
    f = tryparse(Float64, s)
    return !isnothing(f) ? f : s # Return string if all else fails
end

function parse_metrics_file(file_path::String)
    content = try
        read(file_path, String)
    catch
        return nothing
    end
    
    if !occursin(TRIGGER_LINE, content)
        return nothing
    end

    data = Dict{String, Any}()
    for line in split(content, '\n')
        m = match(r"#\s*([^:]+?)\s*:\s*(.*)", line)
        if !isnothing(m)
            label = strip(m[1])
            value_str = strip(m[2])
            data[label] = smart_parse(value_str)
        end
    end
    
    return isempty(data) ? nothing : data
end

function (@main)(args)
    input_dir = length(args) >= 1 ? args[1] : DEFAULT_INPUT_DIR
    
    if !isdir(input_dir)
        println("Error: Input directory '$input_dir' not found.")
        return 1
    end

    all_data_raw = Vector{Dict{String, Any}}()

    files = filter(f -> isfile(joinpath(input_dir, f)) && !startswith(f, "."), readdir(input_dir))

    for file_name in files
        full_path = joinpath(input_dir, file_name)
        result = parse_metrics_file(full_path)

        if !isnothing(result)
            push!(all_data_raw, result)
        end
    end

    if isempty(all_data_raw)
        println("No valid data found.")
        return 0
    end

    # Collect all unique keys
    all_keys = Set{String}()
    for d in all_data_raw
        for k in keys(d)
            push!(all_keys, k)
        end
    end

    # Standardize dictionaries by filling missing keys with missing
    all_data = [Dict(k => get(d, k, missing) for k in all_keys) for d in all_data_raw]

    # Build DataFrame from vector of dicts
    df = DataFrame(all_data)

    # Reorder columns: 'k' first, then others alphabetically
    cols = names(df)
    if "k" in cols
        other_cols = sort(filter(c -> c != "k", cols))
        select!(df, vcat(["k"], other_cols))
    else
        select!(df, sort(cols))
    end

    # Sort by the first column ('k')
    sort!(df, 1)

    # CSV.write preserves the integer formatting in the output file
    CSV.write(DEFAULT_OUTPUT_FILE, df)

    println("Success! $(nrow(df)) rows written to $DEFAULT_OUTPUT_FILE")
    return 0
end



