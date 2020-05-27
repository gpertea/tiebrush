# build gclib
set(gclib_PREFIX ${CMAKE_BINARY_DIR}/contrib/gclib)

ExternalProject_Add(gclib
                    PREFIX ${gclib_PREFIX}
                    GIT_REPOSITORY "https://github.com/gpertea/gclib.git"
                    CONFIGURE_COMMAND ""
                    BUILD_COMMAND ""
                    INSTALL_COMMAND ""
)

set(gclib_LIB ${gclib_PREFIX}/src/gclib/)