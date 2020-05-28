# build zlib
set(zlib_PREFIX ${CMAKE_BINARY_DIR}/contrib/zlib-prefix)
set(zlib_INSTALL ${CMAKE_BINARY_DIR}/contrib/zlib-install)

ExternalProject_Add(zlib
        PREFIX ${zlib_PREFIX}
        GIT_REPOSITORY "https://github.com/madler/zlib.git"
        GIT_TAG "v1.2.11"
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ${zlib_PREFIX}/src/zlib/configure --prefix=${zlib_INSTALL} --static
        INSTALL_DIR ${zlib_INSTALL}
        LOG_DOWNLOAD 1
        LOG_INSTALL 1
        )

include_directories(${zlib_INSTALL}/include)
set(zlib_LIB ${zlib_INSTALL}/lib/libz.a)