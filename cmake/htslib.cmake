set(htslib_PREFIX ${CMAKE_BINARY_DIR}/contrib/htslib-prefix)
set(htslib_INSTALL ${CMAKE_BINARY_DIR}/contrib/htslib-install)

set (CMAKE_STATIC_LINKER_FLAGS "-Wl,--as-needed.-lcurl")

ExternalProject_Add(htslib
        PREFIX ${htslib_PREFIX}
        GIT_REPOSITORY "https://github.com/samtools/htslib.git"
        GIT_TAG ""
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND autoheader COMMAND autoconf COMMAND ${htslib_PREFIX}/src/htslib/configure
        BUILD_COMMAND make -C ${htslib_PREFIX}/src/htslib "CFLAGS=-g -Wall -O3" libhts.a
        INSTALL_COMMAND ""
        )

add_dependencies(htslib zlib)
include_directories(${htslib_PREFIX}/src/htslib)
set(htslib_LIB ${htslib_PREFIX}/src/htslib/libhts.a)