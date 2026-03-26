#! /bin/env julial

using DelimitedFiles

# Configuration
input_list = "misc/summary_files.txt"
output_csv = "misc/results.csv"

# Regex patterns
# Captures numbers (integers or floats) after the specific triggers
k_regex = r"k:\s*(\d+\.?\d*)"
time_regex = r"# Running time in seconds\s*:\s*(\d+\.?\d*)"

function process_files()
    if !isfile(input_list)
        println("Error: List file '$input_list' not found.")
        return
    end

    # Read the list of filenames
    files = filter(!isempty, readlines(input_list))
    
    # Prepare the results storage with a header
    results = ["File Path" "k_value" "Running_Time_Sec"]

    for file_path in files
        if isfile(file_path)
            content = read(file_path, String)
            
            # Match patterns
            k_match = match(k_regex, content)
            #println( "k_match:", k_match )
            time_match = match(time_regex, content)
            #println( "time_match:", time_match )
            
            # Extract values or use "N/A" if not found
            k_val = k_match !== nothing ? k_match.captures[1] : "N/A"
            time_val = time_match !== nothing ? time_match.captures[1] : "N/A"

            #println( "k_val: ", k_val )
            #println( "time_val: ", time_val )
            # Append as a new row (1x3 matrix for hcat-style storage)
            results = vcat(results, [file_path k_val time_val])
            #println( "vcat done" )
        else
            println("Warning: Skipping '$file_path' (File not found)")
        end
    end

    # Write to CSV using a comma delimiter
    writedlm(output_csv, results, ',')
    println("Success! Data saved to $output_csv")
end

process_files()
