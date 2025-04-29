using PackageCompiler
using Pkg

# Create a temporary environment to resolve dependencies
temp_dir = mktempdir()
Pkg.activate(temp_dir)

# Add the packages and let Julia resolve dependencies
Pkg.add(["ITensors", "ITensorMPS", "TCIITensorConversion"])

# Create the system image
create_sysimage(
    ["ITensors", "ITensorMPS", "TCIITensorConversion"];
    sysimage_path=joinpath(homedir(), ".julia/sysimages/sys_itensors.dylib")
)

# Clean up
Pkg.activate() 