# --------------------------------------------------------------
# Customizations for kernels of example application
# --------------------------------------------------------------

# CMake initialization: No customizations
# Build parameters: No customizations
# Build targets: No customizations
# Build requirements: No customizations
# Installation: No customizations

# --------------------------------------------------------------
# Generate Fortran files
# --------------------------------------------------------------

# Floating-point precision
if (SINGLE_PRECISION)

    set(REAL "real*4")
    set(TINY_NONZERO_NUMBER 1.d-35)

else (SINGLE_PRECISION)

    set(REAL "real*8")
    set(TINY_NONZERO_NUMBER 1.d-99)

endif (SINGLE_PRECISION)

# Generate fortran files
configure_file(interface_kernels_2d.f.in interface_kernels_2d.f)
configure_file(interface_kernels_3d.f.in interface_kernels_3d.f)
configure_file(pore_kernels_2d.f.in pore_kernels_2d.f)
configure_file(pore_kernels_3d.f.in pore_kernels_3d.f)
