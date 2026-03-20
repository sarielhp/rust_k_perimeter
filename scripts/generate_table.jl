#! /bin/env julial

using CSV
using DataFrames
using Printf

function format_commas(n::Integer)
    s = string(n)
    return replace(s, r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
end

function generate_professional_latex(file_path::String)
    if !isfile(file_path)
        error("File not found: $file_path")
    end

    df = CSV.read(file_path, DataFrame)

    println("% Add \\usepackage{booktabs} and \\usepackage{siunitx} to your LaTeX preamble")
    println("\\begin{table}[h]")
    println("  \\centering")
    println("  \\small")
    # S column type from siunitx aligns decimals perfectly
    println("  \\begin{tabular}{r S[table-format=4.2]}") 
    println("    \\toprule")
    println("    {k Value} & {Running Time (s)} \\\\")
    println("    \\midrule")

    for row in eachrow(df)
        k_str = format_commas(row.k_value)
        # Using { } around the comma-string so siunitx doesn't try to parse it as a number
        time_val = @sprintf("%.2f", row.Running_Time_Sec)
        
        println("    {$(k_str)} & $(time_val) \\\\")
    end

    println("    \\bottomrule")
    println("  \\end{tabular}")
    println("  \\caption{Performance metrics for increasing \$k\$ values.}")
    println("\\end{table}")
end

generate_professional_latex("misc/results.csv")
