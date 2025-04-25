using PackageCompiler
using Pkg

# Make sure ITensors is installed
Pkg.add("ITensors")

# Create the system image
create_sysimage(
    ["ITensors"];
    sysimage_path="/home/coding/.julia/sysimages/sys_itensors.so"
) 