### Better safe than sorry...
list(APPEND ${quant3d_LIB_NAME}_COMPILE_OPTIONS "-Wall" "-Wextra" "-Werror" "-Wshadow")

### Add all source files
list(APPEND ${quant3d_LIB_NAME}_SOURCE_LIST
     "MCQuantizer.cpp"
     "ConstraintExtractor.cpp"
     "ConstraintWriter.cpp"
     )
add_library(${quant3d_LIB_NAME} STATIC ${${quant3d_LIB_NAME}_SOURCE_LIST})

### Find dependencies
# mc3d
if (NOT TARGET mc3d)
    add_subdirectory(${PROJECT_SOURCE_DIR}/extern/MC3D extern/MC3D EXCLUDE_FROM_ALL)
endif()
list(APPEND ${quant3d_LIB_NAME}_LIB_LIST mc3d)

#gurobi #TODO free software alternative at some point
find_package(GUROBI REQUIRED)
list(APPEND ${quant3d_LIB_NAME}_LIB_LIST Gurobi::GurobiC)
list(APPEND ${quant3d_LIB_NAME}_LIB_LIST Gurobi::GurobiCXX)

### Link with dependencies
target_link_libraries_system(${quant3d_LIB_NAME} PUBLIC OpenVolumeMesh::OpenVolumeMesh ${${quant3d_LIB_NAME}_LIB_LIST})

### Set compile options
set_target_properties(${quant3d_LIB_NAME} PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON)
target_compile_options(${quant3d_LIB_NAME} PRIVATE ${${quant3d_LIB_NAME}_COMPILE_OPTIONS})
target_compile_definitions(${quant3d_LIB_NAME} PRIVATE ${${quant3d_LIB_NAME}_COMPILE_DEFINITIONS})

### Include own headers with public access
target_include_directories(${quant3d_LIB_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)

### Include headers as system headers to prevent warnings/errors from strict compile definitions
target_include_directories(${quant3d_LIB_NAME} SYSTEM PUBLIC ${GMPXX_INCLUDE_DIR})
target_include_directories(${quant3d_LIB_NAME} SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/extern/MC3D/include)
target_include_directories(${quant3d_LIB_NAME} SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/extern/MC3D/extern/OpenVolumeMesh/src)
target_include_directories(${quant3d_LIB_NAME} SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/extern/MC3D/extern/glog/src)
