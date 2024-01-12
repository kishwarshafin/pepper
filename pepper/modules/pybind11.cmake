# Include the ExternalProject module
include(ExternalProject)

# Set the version and URL for the Pybind11 library
set(PYBIND11_VERSION "2.6.2")
set(PYBIND11_URL "https://github.com/pybind/pybind11/archive/v${PYBIND11_VERSION}.tar.gz")

# Use ExternalProject_Add to download and extract Pybind11
ExternalProject_Add(
    pybind11
    URL ${PYBIND11_URL}
    PREFIX ${CMAKE_BINARY_DIR}/pybind11-${PYBIND11_VERSION}
    DOWNLOAD_EXTRACT_TIMESTAMP true
    # Specify where to install Pybind11 (if necessary)
    # INSTALL_DIR ${CMAKE_BINARY_DIR}/pybind11-install
    # Add any additional configuration options for Pybind11 here
)

# Set the include directory for Pybind11 (useful for other CMake files to locate Pybind11 headers)
set(PYBIND11_INCLUDE_DIR "${CMAKE_BINARY_DIR}/pybind11-${PYBIND11_VERSION}/src/pybind11/include")


# Optional: Add a custom target to build Pybind11 (if needed)
# This can be useful if you want to explicitly control when Pybind11 is built
# add_custom_target(build_pybind11 DEPENDS pybind11)

