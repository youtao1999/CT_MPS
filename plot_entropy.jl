using Plots
using JSON
using Statistics
using Printf

"""
Load data from a JSON file and extract p values and entanglement entropy.
"""
function load_data(filename)
    # Read JSON file
    data = JSON.parsefile(filename)
    results = data["results"]
    args = data["args"]
    
    # Extract and sort p values
    p_values = sort(parse.(Float64, collect(keys(results))))
    
    # Get entropy values and errors
    if args["ancilla"] == 0
        entropy = [results[string(p)]["avg_EE"] for p in p_values]
        entropy_std = [results[string(p)]["std_EE"] for p in p_values]
        entropy_label = "EE"
    else
        entropy = [results[string(p)]["avg_SA"] for p in p_values]
        entropy_std = [results[string(p)]["std_SA"] for p in p_values]
        entropy_label = "SA"
    end
    
    return p_values, entropy, entropy_std, args, entropy_label
end

"""
Plot entanglement entropy vs p_scan for a given file.
"""
function plot_entropy(filename; show_plot=true)
    p_values, entropy, entropy_std, args, entropy_label = load_data(filename)
    
    fixed_type = args["scan_type"] == "p_ctrl" ? "p_proj" : "p_ctrl"
    
    # Create plot
    plot = Plots.plot(
        p_values, entropy,
        ribbon=entropy_std,
        xlabel=args["scan_type"],
        ylabel="Ensemble Average $(entropy_label)",
        title="L=$(args["L"]), $(fixed_type)=$(args["p_fixed"]), $(args["n_realizations"]) realizations",
        marker=:circle,
        line=:solid,
        legend=false,
        grid=true,
        dpi=300
    )
    
    if show_plot
        display(plot)
    end
    return plot
end

"""
Plot multiple datasets on the same graph.
"""
function plot_multiple_datasets(filenames)
    plot_obj = nothing
    
    for (i, filename) in enumerate(filenames)
        p_values, entropy, entropy_std, args, entropy_label = load_data(filename)
        
        if i == 1
            plot_obj = Plots.plot(
                xlabel="p",
                ylabel="Ensemble Average $(entropy_label)",
                title="L=$(args["L"]), $(args["n_realizations"]) realizations",
                grid=true,
                dpi=300
            )
        end
        
        fixed_type = args["scan_type"] == "p_ctrl" ? "p_proj" : "p_ctrl"
        label = "$(fixed_type)=$(args["p_fixed"])"
        
        Plots.plot!(
            plot_obj,
            p_values, entropy,
            ribbon=entropy_std,
            label=label,
            marker=:circle,
            line=:solid
        )
    end
    
    display(plot_obj)
    return plot_obj
end

# Example usage:
if abspath(PROGRAM_FILE) == @__FILE__
    # List all JSON files
    json_files = filter(f -> endswith(f, ".json"), readdir())
    println("Available data files:")
    for (i, f) in enumerate(json_files)
        println("$i. $f")
    end
    
    if !isempty(json_files)
        # Plot first file
        println("\nPlotting first file...")
        plot_entropy(json_files[1])
        
        # If there are multiple files, plot them together
        if length(json_files) > 1
            println("\nPlotting multiple datasets...")
            plot_multiple_datasets(json_files[1:min(3, length(json_files))])
        end
    end
end 