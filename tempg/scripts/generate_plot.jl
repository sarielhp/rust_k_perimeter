#! /bin/env julial

using CSV
using DataFrames
using LsqFit
using Plots
using Printf

# Use GR backend for clean PDF vector graphics
gr()

function run_split_pdf_analysis(filename)
    if !isfile(filename)
        println("Error: File '$filename' not found.")
        return
    end

    # 1. Create the output directory if it doesn't exist
    plot_dir = "plots"
    if !isdir(plot_dir)
        mkpath(plot_dir)
        println("Created directory: $plot_dir/")
    end

    # 2. Load Data (Columns 2 and 3)
    df = CSV.read(filename, DataFrame)
    xdata = Float64.(df[:, 2])
    ydata = Float64.(df[:, 3])

    # 3. Grid Search (0.3 to 3.0 in steps of 0.05)
    test_degrees = 0.3:0.01:3.0
    best_mse = Inf
    best_p = [1.0, 1.0, 0.0] 
    mses = Float64[]

    for d in test_degrees
        @. fixed_model(x, p) = p[1] * x^d + p[2]
        p0 = [1.0, sum(ydata)/length(ydata)]
        
        try
            fit = curve_fit(fixed_model, xdata, ydata, p0)
            mse = sum(fit.resid.^2) / length(ydata)
            push!(mses, mse)

            if mse < best_mse
                best_mse = mse
                best_p = [coef(fit)[1], d, coef(fit)[2]]
            end
        catch
            push!(mses, NaN)
        end
    end

    a, best_d, c = best_p

    # 4. Generate Plot 1: The Fit Curve
    x_range = range(minimum(xdata), stop=maximum(xdata), length=200)
    y_fit = a .* (x_range .^ best_d) .+ c

    p1 = scatter(xdata, ydata, 
                label="Data Points", 
                color=:black, 
                markerstrokewidth=0,
                title="Optimal Fractional Fit (d = $best_d)",
                xlabel="Input Parameter",
                ylabel="Measured Output")
    plot!(p1, x_range, y_fit, label="Fit Line", lw=3, linecolor=:red)

    # 5. Generate Plot 2: The MSE Search
    p2 = plot(collect(test_degrees), mses, 
              lw=2, color=:blue, label="MSE",
              title="MSE Search (Winner: $best_d)",
              xlabel="Tested Degree",
              ylabel="Mean Squared Error")
    scatter!(p2, [best_d], [best_mse], label="Minimum MSE", markersize=6, markercolor=:gold)

    # 6. Save as individual PDF files in the plots/ directory
    path1 = joinpath(plot_dir, "fit_curve.pdf")
    path2 = joinpath(plot_dir, "mse_search.pdf")
    
    savefig(p1, path1)
    savefig(p2, path2)
    
    println("------------------------------------")
    println("Success!")
    @printf("Winning Degree: %.2f\n", best_d)
    println("Files saved to $plot_dir/:")
    println("  - fit_curve.pdf")
    println("  - mse_search.pdf")
    println("------------------------------------")
end

run_split_pdf_analysis("misc/results.csv")
