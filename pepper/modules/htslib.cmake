# build htslib
set(htslib_PREFIX ${CMAKE_BINARY_DIR}/htslib)
set (FLAGS "-fPIC")
# Enable ExternalProject CMake module
include(ExternalProject)

ExternalProject_Add(htslib
        BUILD_IN_SOURCE 1
        URL https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
        PREFIX ${htslib_PREFIX}
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND autoconf && ./configure --prefix=${CMAKE_BINARY_DIR}/htslib --disable-lzma --disable-s3 --disable-plugins --disable-bz2 --disable-libcurl
        BUILD_COMMAND make CFLAGS=${CMAKE_C_FLAGS}
        INSTALL_COMMAND "")
ExternalProject_Get_Property(htslib SOURCE_DIR)
set(HTSLIB_SRC_DIR ${SOURCE_DIR})

include_directories("${HTSLIB_SRC_DIR}/htslib")