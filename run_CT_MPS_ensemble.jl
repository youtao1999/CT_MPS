using ArgParse
using JSON
using Printf
using Random
using LinearAlgebra
using MKL
using Serialization
using ITensors
using ITensorMPS
using CT
using Distributed
using Statistics
using Base: values
using Profile

"""
This script computes ensemble-averaged late-time entanglement entropy with 1/3 and 2/3 global control maps.
It scans over either p_ctrl or p_proj while keeping the other fixed, averaging over multiple realizations.
"""

@everywhere begin
    using Random
    using LinearAlgebra
    using MKL
    using ITensors
    using ITensorMPS
    using CT
    using Statistics
    using Printf
    using Base: values

    """
    Class to compute and manage ensemble calculations for CT_MPS
    """
    struct CTMPSComputer
        L::Int
        p_fixed::Float64
        p_fixed_name::String
        p_scan_values::Vector{Float64}
        ensemble_size::Int  # Total number of realizations
        ancilla::Int
        maxdim::Int
        
        function CTMPSComputer(L, p_fixed, p_fixed_name, p_scan_values, ensemble_size, ancilla, maxdim)
            @assert p_fixed_name in ["p_ctrl", "p_proj"] "p_fixed_name must be either 'p_ctrl' or 'p_proj'"
            new(L, p_fixed, p_fixed_name, p_scan_values, ensemble_size, ancilla, maxdim)
        end
    end

    """
    Run a single realization of the circuit
    """
    function run_single(computer::CTMPSComputer, p_ctrl::Float64, p_proj::Float64, seed::Int)
        ct_f = CT.CT_MPS(
            L=computer.L, 
            seed=seed, 
            folded=true, 
            store_op=false, 
            store_vec=false,
            ancilla=computer.ancilla, 
            debug=false, 
            xj=Set([1//3,2//3]), 
            _maxdim=computer.maxdim
        )
        
        i = 1
        T_max = computer.ancilla == 0 ? 2*(ct_f.L^2) : div(ct_f.L^2,2)
        
        for _ in 1:T_max
            i = CT.random_control!(ct_f, i, p_ctrl, p_proj)
        end
        
        if computer.ancilla == 0
            EE = CT.von_Neumann_entropy(ct_f.mps, div(ct_f.L,2))
            return Dict("EE" => EE)
        else
            SA = CT.von_Neumann_entropy(ct_f.mps, 1)
            return Dict("SA" => SA)
        end
    end

    """
    Compute assigned realizations for given worker
    """
    function compute_chunk(computer::CTMPSComputer, worker_id::Int, n_workers::Int)
        # Ensure ensemble size is divisible by number of workers
        if computer.ensemble_size % n_workers != 0
            error("Ensemble size ($(computer.ensemble_size)) must be divisible by number of workers ($n_workers).\n" *
                  "Valid worker counts are: $(join([i for i in 1:computer.ensemble_size if computer.ensemble_size % i == 0], ", "))")
        end
        
        # Calculate chunk size (number of realizations per worker)
        chunk_size = computer.ensemble_size ÷ n_workers  # Using integer division since we know it's divisible
        start_real = (worker_id - 1) * chunk_size + 1
        end_real = worker_id * chunk_size  # No need for min() since we ensure even division
        
        println("Worker $(worker_id): handling realizations $(start_real) to $(end_real)")
        
        chunk_results = Dict()
        
        for p_scan in computer.p_scan_values
            # Set parameters based on which one is fixed
            if computer.p_fixed_name == "p_ctrl"
                p_ctrl = computer.p_fixed
                p_proj = p_scan
            else  # p_fixed_name == "p_proj"
                p_ctrl = p_scan
                p_proj = computer.p_fixed
            end
            
            t_start = time()
            last_update = t_start
            
            # Run assigned realizations
            realizations = Dict{Int, Dict{String, Float64}}()
            for i in start_real:end_real
                realization_start = time()
                realizations[i] = run_single(computer, p_ctrl, p_proj, i)
                
                # Print progress every 60 seconds
                if time() - last_update > 60
                    progress = (i - start_real + 1) / chunk_size * 100
                    elapsed = time() - t_start
                    avg_time = elapsed / (i - start_real + 1)
                    remaining = avg_time * (end_real - i)
                    println("Worker $(worker_id): $(round(progress, digits=1))% complete. Avg time per realization: $(round(avg_time, digits=2))s. Est. remaining: $(round(remaining/60, digits=1)) minutes")
                    last_update = time()
                end
            end
            
            total_time = time() - t_start
            println("Worker $(worker_id) completed p=$(round(p_scan, digits=3)) in $(round(total_time/60, digits=1)) minutes ($(round(total_time/chunk_size, digits=2))s per realization)")
            
            chunk_results[string(p_scan)] = realizations
        end
        
        # Debug print for first p_scan value only
        first_p = string(first(computer.p_scan_values))
        println("Worker $(worker_id) sample data for p=$(first_p):")
        println("First realization: ", chunk_results[first_p][start_real])
        println("Last realization: ", chunk_results[first_p][end_real])
        
        return chunk_results
    end
end

"""
Parse the p_range argument into a list of values
"""
function parse_p_range(p_range_str::String)
    if contains(p_range_str, ":")
        # Format: "start:stop:num"
        parts = split(p_range_str, ":")
        if length(parts) != 3
            error("Invalid range format. Expected start:stop:num")
        end
        start = parse(Float64, parts[1])
        stop = parse(Float64, parts[2])
        num = parse(Int, parts[3])
        return range(start, stop, length=num) |> collect
    else
        # Format: "0.1,0.2,0.3"
        return parse.(Float64, split(p_range_str, ","))
    end
end

function parse_my_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--scan_type"
            arg_type = String
            required = true
            help = "Which parameter to scan: 'p_ctrl' or 'p_proj'"
        "--p_fixed"
            arg_type = Float64
            required = true
            help = "Fixed value for the non-scanned parameter"
        "--p_range"
            arg_type = String
            required = true
            help = "Range of p values to scan. Formats: '0.1,0.2,0.3' or '0.1:0.5:5'"
        "--L", "-L"
            arg_type = Int
            default = 8
            help = "system size"
        "--ancilla", "-a"
            arg_type = Int
            default = 0
            help = "number of ancilla"
        "--maxdim", "-m"
            arg_type = Int
            default = 10
            help = "set the maximal bond dim"
        "--n_realizations"
            arg_type = Int
            default = 2000
            help = "total number of realizations in the ensemble"
    end
    return parse_args(s)
end

"""
Run a single realization with profiling (for single-process mode)
"""
function run_single_with_profiling(computer::CTMPSComputer, p_ctrl::Float64, p_proj::Float64, seed::Int)
    @profile begin
        ct_f = CT.CT_MPS(
            L=computer.L, 
            seed=seed, 
            folded=true, 
            store_op=false, 
            store_vec=false,
            ancilla=computer.ancilla, 
            debug=false, 
            xj=Set([1//3,2//3]), 
            _maxdim=computer.maxdim
        )
        
        i = 1
        T_max = computer.ancilla == 0 ? 2*(ct_f.L^2) : div(ct_f.L^2,2)
        
        for _ in 1:T_max
            i = CT.random_control!(ct_f, i, p_ctrl, p_proj)
        end
        
        if computer.ancilla == 0
            EE = CT.von_Neumann_entropy(ct_f.mps, div(ct_f.L,2))
            return Dict("EE" => EE)
        else
            SA = CT.von_Neumann_entropy(ct_f.mps, 1)
            return Dict("SA" => SA)
        end
    end
end

"""
Compute all realizations in single-process mode with profiling
"""
function compute_single_process(computer::CTMPSComputer)
    all_results = Dict()
    
    for p_scan in computer.p_scan_values
        # Set parameters based on which one is fixed
        if computer.p_fixed_name == "p_ctrl"
            p_ctrl = computer.p_fixed
            p_proj = p_scan
        else  # p_fixed_name == "p_proj"
            p_ctrl = p_scan
            p_proj = computer.p_fixed
        end
        
        println("Running $(computer.ensemble_size) realizations for p=$(round(p_scan, digits=3))...")
        t_start = time()
        last_update = t_start
        
        # Run all realizations
        realizations = Dict{Int, Dict{String, Float64}}()
        for i in 1:computer.ensemble_size
            realizations[i] = run_single_with_profiling(computer, p_ctrl, p_proj, i)
            
            # Print progress every 60 seconds
            if time() - last_update > 60
                progress = i / computer.ensemble_size * 100
                elapsed = time() - t_start
                avg_time = elapsed / i
                remaining = avg_time * (computer.ensemble_size - i)
                println("$(round(progress, digits=1))% complete. Avg time per realization: $(round(avg_time, digits=2))s. Est. remaining: $(round(remaining/60, digits=1)) minutes")
                last_update = time()
            end
        end
        
        total_time = time() - t_start
        println("Completed p=$(round(p_scan, digits=3)) in $(round(total_time/60, digits=1)) minutes ($(round(total_time/computer.ensemble_size, digits=2))s per realization)")
        
        all_results[string(p_scan)] = realizations
    end
    
    return all_results
end

function main()
    # Ensure we're using 1 thread per process
    BLAS.set_num_threads(1)
    
    # Get number of workers
    n_workers = nworkers()
    println("Number of workers: ", n_workers)
    println("BLAS threads per worker: ", BLAS.get_num_threads())
    
    args = parse_my_args()
    p_values = parse_p_range(args["p_range"])
    
    # Create computer instance
    computer = CTMPSComputer(
        args["L"],
        args["p_fixed"],
        args["scan_type"],
        p_values,
        args["n_realizations"],  # This is the ensemble_size
        args["ancilla"],
        args["maxdim"]
    )
    
    # Detect if we're running in single-process mode (for profiling)
    if n_workers == 1
        println("Running in single-process mode with profiling enabled...")
        println("Total ensemble size: $(args["n_realizations"])")
        
        # Clear previous profiling data
        Profile.clear()
        
        # Run computation with profiling
        worker_results = compute_single_process(computer)
        
        # Print profiling results
        println("\n" * "="^50)
        println("PROFILING RESULTS")
        println("="^50)
        Profile.print(mincount=10)  # Only show functions called at least 10 times
        
        # Also save profiling data to file
        Profile.print(IOContext(open("profile_results.txt", "w"), :displaysize => (24, 200)), mincount=5)
        println("\nDetailed profiling results saved to profile_results.txt")
        
    else
        # Calculate and display work distribution
        chunk_size = ceil(Int, args["n_realizations"] / n_workers)
        println("Total ensemble size: $(args["n_realizations"])")
        println("Chunk size per worker: $chunk_size")
        println("Starting distributed computation...")
        
        # Use @distributed to collect results from all workers
        worker_results = @distributed (merge) for worker_id in 1:n_workers
            compute_chunk(computer, worker_id, n_workers)
        end
        
        # Debug print merged results structure
        println("\nMerged results structure:")
        first_p = string(first(p_values))
        println("Number of realizations for p=$(first_p): ", length(keys(worker_results[first_p])))
        println("Sample realizations: ", join(["$(k)=>$(v)" for (k,v) in first(worker_results[first_p], 3)], ", "))
    end
    
    # Process results for each p_scan value (same for both modes)
    all_results = Dict()
    for p_scan in p_values
        # Collect all realizations for this p_scan value
        all_realizations = worker_results[string(p_scan)]
        
        # Calculate statistics
        if args["ancilla"] == 0
            entropy_values = [r["EE"] for r in values(all_realizations)]
            avg_key, std_key = "avg_EE", "std_EE"
        else
            entropy_values = [r["SA"] for r in values(all_realizations)]
            avg_key, std_key = "avg_SA", "std_SA"
        end
        
        # Store results with proper parameter values
        if args["scan_type"] == "p_ctrl"
            p_ctrl = p_scan
            p_proj = args["p_fixed"]
        else
            p_ctrl = args["p_fixed"]
            p_proj = p_scan
        end
        
        all_results[string(p_scan)] = Dict(
            avg_key => mean(entropy_values),
            std_key => std(entropy_values),
            "p_ctrl" => p_ctrl,
            "p_proj" => p_proj
        )
    end
    
    # Debug print final statistics
    println("\nFinal statistics structure:")
    first_p = string(first(p_values))
    println("Sample for p=$(first_p): ", all_results[first_p])
    
    # Save results
    filename = "MPS_ensemble_(1-3,2-3)_L$(args["L"])_$(args["scan_type"])_scan_fixed$(@sprintf("%.3f", args["p_fixed"]))_a$(args["ancilla"]).json"
    data_to_serialize = merge(Dict("results" => all_results), Dict("args" => args))
    
    println("\nSaving to JSON...")
    open(filename, "w") do f
        write(f, JSON.json(data_to_serialize))
    end
    
    println("Computation completed. Results saved to $filename")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

# Example usage:
# # For scanning p_ctrl with fixed p_proj=0.0:
# singularity exec julia_itensor.sif julia --threads=auto run_CT_MPS_ensemble.jl --scan_type p_ctrl --p_fixed 0.0 --p_range "0.45:0.55:10" --L 8 --ancilla 0

# # For scanning p_proj with fixed p_ctrl=0.5:
# singularity exec julia_itensor.sif julia --threads=auto run_CT_MPS_ensemble.jl --scan_type p_proj --p_fixed 0.5 --p_range "0.0,0.1,0.2,0.3,0.4,0.5" --L 8 --ancilla 0