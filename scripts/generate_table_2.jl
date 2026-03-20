#! /bin/env julial

using CSV
using DataFrames
using Printf

# Versatile comma formatter for both Int and Float
function format_commas(val)
    # Split into integer and decimal parts
    parts = split(@sprintf("%.2f", val), '.')
    int_part = parts[1]
    dec_part = parts[2]
    
    # Add commas to the integer part
    formatted_int = replace(int_part, r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
    
    # Only return decimal if it's a float or requested
    return val isa Integer ? formatted_int : "$formatted_int.$dec_part"
end

function generate_parallel_latex(file_path::String)
    if !isfile(file_path)
        error("File not found: $file_path")
    end

    df = CSV.read(file_path, DataFrame)
    n_rows = nrow(df)
    mid = Int(ceil(n_rows / 2)) # Find the split point

    println("% Add \\usepackage{booktabs} to your LaTeX preamble")
    println("\\begin{table}[ht]")
    println("  \\centering")
    println("  \\small")
    # Define 4 columns: k, time, spacer, k, time
    # 'p{2em}' acts as a gutter between the two table halves
    println("  \\begin{tabular}{rr @{\\hspace{3em}} rr}") 
    println("    \\toprule")
    println("    \$k\$ Value & Time (s) & \$k\$ Value & Time (s) \\\\")
    println("    \\midrule")

    for i in 1:mid
        # Left side data
        k1 = format_commas(df[i, :k_value])
        t1 = format_commas(df[i, :Running_Time_Sec])

        # Right side data (check if index exists for odd-numbered totals)
        if i + mid <= n_rows
            k2 = format_commas(df[i + mid, :k_value])
            t2 = format_commas(df[i + mid, :Running_Time_Sec])
            println("    $(k1) & $(t1) & $(k2) & $(t2) \\\\")
        else
            println("    $(k1) & $(t1) & & \\\\")
        end
    end

    println("    \\bottomrule")
    println("  \\end{tabular}")
    println("  \\caption{Parallel summary of running times for various \$k\$ values.}")
    println("\\end{table}")
end

generate_parallel_latex("misc/results.csv")
