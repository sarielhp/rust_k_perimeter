#! /bin/env julial

using CSV
using DataFrames
using Printf

# Formatter for integers and floats with commas
function format_commas(val)
    if val isa Integer
        s = string(val)
        return replace(s, r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
    else
        parts = split(@sprintf("%.2f", val), '.')
        int_part = replace(parts[1], r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
        return "$int_part.$(parts[2])"
    end
end

function generate_triple_parallel_latex(file_path::String)
    if !isfile(file_path)
        error("File not found: $file_path")
    end

    df = CSV.read(file_path, DataFrame)
    n_rows = nrow(df)
    
    # Calculate split points for 3 columns
    chunk_size = Int(ceil(n_rows / 3))
    
    println("% Add \\usepackage{booktabs} to your LaTeX preamble")
    println("\\begin{table}[ht]")
    println("  \\centering")
    println("  \\small")

    for col in 0:2
        start_idx = col * chunk_size + 1
        end_idx = min((col + 1) * chunk_size, n_rows)
        
        # Guard against empty chunks if data is very small
        if start_idx > n_rows
            break
        end

        # Start a separate tabular for each chunk
        println("  \\begin{tabular}[t]{rr}") # [t] aligns them by the top
        println("    \\toprule")
        println("    \$k\$ Value & Time (s) \\\\")
        println("    \\midrule")
        
        for i in start_idx:end_idx
            k_str = format_commas(df[i, :k_value])
            t_str = format_commas(df[i, :Running_Time_Sec])
            println("    $(k_str) & $(t_str) \\\\")
        end
        
        println("    \\bottomrule")
        println("  \\end{tabular}")
        
        # Add spacing between tables, but not after the last one
        if col < 2
            println("  \\hfill") 
        end
    end

    println("  \\caption{Performance summary split into three independent parallel sections.}")
    println("\\end{table}")
end

generate_triple_parallel_latex("misc/results.csv")
