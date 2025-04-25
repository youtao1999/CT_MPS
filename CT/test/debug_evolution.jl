using Pkg
Pkg.activate("CT")
using CT

# Initialize system parameters
L = 8
steps = 20  # Total evolution time
p_ctrl = 0.5   # Probability of control operation
p_proj = 0.0   # Probability of projection
x0 = 1//2

# Create CT_MPS object with debugging and history tracking enabled
ct = CT.CT_MPS(
    L=L,
    seed=42,           # Fixed seed for reproducibility
    store_op=true,     # Store operation history
    store_vec=false,   # Don't store state history (to save memory)
    debug=true,        # Enable debug output
    xj=Set([0]),      # Using the FM fixed point setting
    _maxdim0=10,       # Initial bond dimension
    x0=x0
)

# Evolve the system
i = 1  # Starting site
println("Starting evolution for $steps steps...")
println("System size L = $L")
println("Control probability = $p_ctrl")
println("Projection probability = $p_proj")
println("Initial state = $x0 which in binary is $(CT.dec2bin(x0, L))")
println("ram_phy = $(ct.ram_phy)")
println("phy_list = $(ct.phy_list)")
println("phy_ram = $(ct.phy_ram)")
println("\nEvolution log:")
println("-------------")

for step in 1:steps
    global i  # Need to declare i as global since we modify it
    i = CT.random_control!(ct, i, p_ctrl, p_proj)
end

# Print operation history statistics
println("\nOperation Statistics:")
println("-------------------")
n_control = count(op -> op[1]["Type"] == "Control", ct.op_history)
n_bernoulli = count(op -> op[1]["Type"] == "Bernoulli", ct.op_history)
n_projection = count(op -> op[1]["Type"] == "Projection", ct.op_history)

println("Total steps: $steps")
println("Control operations: $n_control")
println("Bernoulli operations: $n_bernoulli")
println("Projection operations: $n_projection")

# Print detailed operation history
println("\nDetailed Operation History:")
println("-------------------------")
for (idx, op_list) in enumerate(ct.op_history)
    println("Step $idx:")
    for op in op_list
        if op["Type"] == "Control"
            println("  Control at site $(op["Site"][1]) with outcome $(op["Outcome"][1])")
        elseif op["Type"] == "Bernoulli"
            println("  Bernoulli map between sites $(op["Site"][1]) and $(op["Site"][2])")
        elseif op["Type"] == "Projection"
            println("  Projection at site $(op["Site"][1]) with outcome $(op["Outcome"][1])")
        end
    end
end 