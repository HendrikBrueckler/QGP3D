set(CLI_NAME "qgp3d_cli")

if (NOT TARGET CLI11::CLI11)
    add_subdirectory(${PROJECT_SOURCE_DIR}/extern/MC3D/extern/CLI11 extern/CLI11 EXCLUDE_FROM_ALL)
endif()

add_executable(${CLI_NAME} main.cpp)
target_link_libraries(${CLI_NAME} ${QGP3D_LIB_NAME} CLI11::CLI11)
set_target_properties(${CLI_NAME} PROPERTIES
                        CXX_STANDARD 17
                        CXX_STANDARD_REQUIRED ON
                        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/cli)
