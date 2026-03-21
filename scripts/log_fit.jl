#! /bin/env julial

using CSV
using DataFrames
using GLM
using Plots
using PGFPlotsX  # <--- Switch to this backend
#using Cairo
using LaTeXStrings

# Set the backend to gr or pstats; Cairo handles the PDF rendering
#gr() 
pgfplotsx()
push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\usepackage{amsmath}")

function analyze_and_plot(file_path::String)
    # 1. Load and clean data
    df = CSV.read(file_path, DataFrame ); #, header=[:n, :time])
#    println( names(df ));
#    exit(-1)

    rename!(df, :k_value => :n, :Running_Time_Sec => :time)
    df = filter(row -> row.n > 0 && row.time > 0, df)
    
    df.log_n = log10.(df.n)
    df.log_time = log10.(df.time)

    println( "Trying to fit the model..." )
    # 2. Linear Regression in Log10 space
    model = lm(@formula(log_time ~ log_n), df)
    slope = coef(model)[2]
    intercept = coef(model)[1]

    println( "Plotting..." )
    # 3. Create the Plot

    p = plot(df.log_n, df.log_time, 
             seriestype = :scatter, 
             label = L"\text{Runs}",
             xlabel = L"\log_{10}(n)",
             ylabel = L"\log_{10}(\text{seconds})",
             title = L"{Runtime Analysis}",
             legend = :topleft)

    plot!(df.log_n, predict(model), 
        label = L"\text{Fit } (c \approx %$(round(coef(model)[2], digits=2)))")
    
    # 4. Save to PDF
    output_dir = "plots"
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    
    output_path = joinpath(output_dir, "rt_fit.pdf")
    savefig(p, output_path)
    
    println("Analysis complete. Exponent: $(round(slope, digits=4))")
    println("Plot saved to: $output_path")
end

# Usage:
analyze_and_plot("misc/results.csv")
