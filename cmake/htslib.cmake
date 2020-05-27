set(htslib_PREFIX ${CMAKE_BINARY_DIR}/contrib/htslib-prefix)
set(htslib_INSTALL ${CMAKE_BINARY_DIR}/contrib/htslib-install)

set (CMAKE_STATIC_LINKER_FLAGS "-Wl,--as-needed.-lcurl")

ExternalProject_Add(htslib
        PREFIX ${htslib_PREFIX}
        GIT_REPOSITORY "https://github.com/samtools/htslib.git"
        GIT_TAG ""
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND autoreconf && ${htslib_PREFIX}/src/htslib/configure --prefix=${htslib_INSTALL}
        BUILD_COMMAND make lib-static
        INSTALL_COMMAND make install prefix=${htslib_INSTALL}
        )

add_dependencies(htslib zlib)
include_directories(${htslib_INSTALL}/include)
set(htslib_LIB ${htslib_INSTALL}/lib/libhts.a)