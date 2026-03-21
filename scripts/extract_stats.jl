#! /bin/env julial

using DataFrames, CSV

# --- Configuration ---
const DEFAULT_INPUT_DIR = "output/summary"
const DEFAULT_OUTPUT_FILE = "misc/extracted_data.csv"
const TRIGGER_LINE = "# Configs computed"

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

    labels = [strip(m[1]) for m in eachmatch(r"#\s*(.*?)\s*:", content)]
    
    # Use our smart_parse instead of direct Float64 parsing
    numbers = [smart_parse(m.match) for m in eachmatch(r"[-+]?\d*\.?\d+([eE][-+]?\d+)?", content)]
    
    return (labels=labels, numbers=numbers)
end

function (@main)(args)
    input_dir = length(args) >= 1 ? args[1] : DEFAULT_INPUT_DIR
    
    if !isdir(input_dir)
        println("Error: Input directory '$input_dir' not found.")
        return 1
    end

    all_rows = Vector{Vector{Any}}() # Changed to Any to hold mixed Int/Float
    column_names = String[]

    files = filter(f -> isfile(joinpath(input_dir, f)) && !startswith(f, "."), readdir(input_dir))

    for file_name in files
        full_path = joinpath(input_dir, file_name)
        result = parse_metrics_file(full_path)

        if !isnothing(result)
            if isempty(column_names)
                column_names = result.labels
            end
            push!(all_rows, result.numbers)
        end
    end

    if isempty(all_rows)
        println("No valid data found.")
        return 0
    end

    # Build DataFrame
    df = DataFrame(reduce(hcat, all_rows)', column_names)

    # Sort by the first column (handles Int128 vs Float64 comparison automatically)
    sort!(df, 1)

    # CSV.write preserves the integer formatting in the output file
    CSV.write(DEFAULT_OUTPUT_FILE, df)

    println("Success! $(nrow(df)) rows written to $DEFAULT_OUTPUT_FILE")
    return 0
end



