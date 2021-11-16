CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

# name of the project is PEPPER
PROJECT(pepper)
set (CMAKE_CXX_FLAGS "-fPIC -O3 -pipe")
set (CMAKE_C_FLAGS "-fPIC -O3 -pipe -shared -rdynamic")

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

FIND_PACKAGE(PythonInterp 3 REQUIRED)
FIND_PACKAGE(PythonLibs 3 REQUIRED)

if (PYTHONINTERP_FOUND)
    message("Python found")
else()
    message("Python not found")
endif()

# enable installing dependencies
option(INSTALL_DEPENDENCIES
        "Install project dependencies"
        ON)

INCLUDE(pepper/modules/htslib.cmake)
INCLUDE(pepper/modules/pybind11.cmake)

# pybind11 to interface
pybind11_add_module(PEPPER pepper/modules/src/pybind_api.cpp)

add_dependencies(PEPPER htslib)

# add all the external libraries
target_link_libraries(PEPPER PRIVATE z)
target_link_libraries(PEPPER PRIVATE ${HTSLIB_SRC_DIR}/libhts.a)

# PEPPER_VARIANT
pybind11_add_module(PEPPER_VARIANT pepper_variant/modules/cpp/pybind_api.cpp)

add_dependencies(PEPPER_VARIANT htslib)

# add all the external libraries
target_link_libraries(PEPPER_VARIANT PRIVATE z)
target_link_libraries(PEPPER_VARIANT PRIVATE ${HTSLIB_SRC_DIR}/libhts.a)

#-------------------
# Test
#-------------------
#ENABLE_TESTING()
#ADD_SUBDIRECTORY(modules/unit_tests)

