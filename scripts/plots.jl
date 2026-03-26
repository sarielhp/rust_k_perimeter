#! /bin/env julial

using CSV
using DataFrames
using Plots
using GLM
using Printf

# Load the data
df = CSV.read("misc/extracted_data.csv", DataFrame)

# Ensure the plots directory exists
if !isdir("plots")
    mkdir("plots")
end

# Define the pairs to plot: (Column name, plot filename, y-axis label)
plot_pairs = [
    ("Configs computed", "configs_computed.pdf", "Configs Computed"),
    ("Longest edge len", "longest_edge_len.pdf", "Longest Edge Length"),
    ("boundary grid points", "boundary_grid_points.pdf", "Boundary Grid Points"),
    ("vertices", "vertices.pdf", "Number of Vertices")
]

#for n in names(df)
#    println("$n: $(eltype(df[!, n]))")
#    println( df[!, n] )
#end
#describe(df, :eltype)
#println( names( df ) )
for (col_name, filename, y_label) in plot_pairs
    println( "filename:", filename)
    println( "col_name [", col_name, "]" );
    println( df[!, Symbol(col_name)] )

    # Filter out rows with zero or negative values for log-log plot
    sub_df = filter(row -> row[:k] > 0 && row[Symbol(col_name)] > 0, df)
    
    #println( "sub_df: ", sub_df )
    if nrow(sub_df) == 0
        println("Skipping $col_name as there is no valid data for log-log plot.")
        continue
    end

    # Extract data
    X = sub_df.k
    Y = sub_df[!, Symbol(col_name)]
    
    log_X = log.(X)
    log_Y = log.(Y)
    
    # Fit linear model: log(Y) = β₀ + β₁ * log(X)
    data_to_fit = DataFrame(log_X = log_X, log_Y = log_Y)
    model = lm(@formula(log_Y ~ log_X), data_to_fit)
    β = coef(model)
    
    # Extract coefficients
    intercept = β[1]
    slope = β[2]
    
    # Generate the plot
    p = scatter(X, Y, 
                xscale=:log10, yscale=:log10, 
                label="Data", 
                xlabel="k", ylabel=y_label,
                title="Log-Log Plot: $y_label vs k",
                markeralpha=0.5,
                legend=:topleft)
    
    # Generate fit line points
    x_fit = range(minimum(X), maximum(X), length=100)
    y_fit = exp.(intercept .+ slope .* log.(x_fit))
    
    # Add fit line to plot
    fit_label = @sprintf("Fit: y = %.2f * x^{%.2f}", exp(intercept), slope)
    plot!(p, x_fit, y_fit, label=fit_label, linewidth=2, color=:red)
    
    # Save the plot
    output_path = joinpath("plots", filename)
    savefig(p, output_path)
    println("Saved plot to $output_path")
end
