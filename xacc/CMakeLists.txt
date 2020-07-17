
find_package(XACC REQUIRED)


add_library(FunctorExecutor STATIC Executor/Functor.cpp Executor/Executor.cpp)
target_compile_options(FunctorExecutor PRIVATE -Wall -Wextra -Wpedantic -Wno-unused-parameter)
target_include_directories(FunctorExecutor PUBLIC  ./Utils ./Executor ${QUAC_INCLUDE_DIRS} ${XACC_QUAC_BASE_INCLUDE_DIR} ./tpls/cereal/include)
target_link_libraries(FunctorExecutor PUBLIC QuaC base lapack)
target_compile_options(FunctorExecutor PRIVATE -fPIC)
target_compile_definitions(FunctorExecutor PRIVATE NODE_SERVICE_EXE_DIR=${CMAKE_CURRENT_BINARY_DIR})

set (LIBRARY_NAME xacc_quac)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file (GLOB HEADERS *.hpp)
file (GLOB SRC *.cpp Accelerator/*.cpp Utils/*.cpp Backend/*.cpp)

usFunctionGetResourceSource(TARGET ${LIBRARY_NAME} OUT SRC)
usFunctionGenerateBundleInit(TARGET ${LIBRARY_NAME} OUT SRC)

add_library(${LIBRARY_NAME} SHARED ${SRC})
add_dependencies(${LIBRARY_NAME} base FunctorExecutor)
# Enable more warnings for our C++ code 
target_compile_options(${LIBRARY_NAME} PRIVATE -Wall -Wextra -Wpedantic -Wno-unused-parameter)

set(_bundle_name xacc_quac_accelerator)
set_target_properties(${LIBRARY_NAME} PROPERTIES
    # This is required for every bundle
    COMPILE_DEFINITIONS US_BUNDLE_NAME=${_bundle_name}
    # This is for convenience, used by other CMake functions
    US_BUNDLE_NAME ${_bundle_name}
    )

# Embed meta-data from a manifest.json file
usFunctionEmbedResources(TARGET ${LIBRARY_NAME}
    WORKING_DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}
    FILES
    manifest.json
    )

message ("XACC_QUAC_BASE_INCLUDE_DIR = ${XACC_QUAC_BASE_INCLUDE_DIR}")

target_include_directories(${LIBRARY_NAME} PUBLIC   . 
                                                    ./Accelerator
                                                    ./Utils
                                                    ./Executor
                                                    ./Backend
                                                    ${XACC_ROOT}/include/armadillo
                                                    ${QUAC_INCLUDE_DIRS}
                                                    ${XACC_QUAC_BASE_INCLUDE_DIR}
                                                    ${XACC_ROOT}/include
                                                    ${XACC_ROOT}/include/eigen)


target_link_libraries(${LIBRARY_NAME} PUBLIC xacc::xacc xacc::quantum_gate FunctorExecutor)
set_target_properties(${LIBRARY_NAME} PROPERTIES LINK_FLAGS "-shared")

set_target_properties(${LIBRARY_NAME} PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")

xacc_configure_plugin_rpath(${LIBRARY_NAME})
set (CMAKE_INSTALL_PREFIX "${XACC_ROOT}" CACHE PATH "default install path" FORCE )

install(FILES $<TARGET_FILE:${LIBRARY_NAME}> DESTINATION ${CMAKE_INSTALL_PREFIX}/plugins) 
add_subdirectory(tests)

# Build MPI PlatformNode
set (MPI_EXECUTOR PlatformNode)
add_executable(${MPI_EXECUTOR} Executor/MPI/PlatformNode.cpp)
target_include_directories(${MPI_EXECUTOR} PRIVATE  ./Executor ./Executor/MPI ./Utils ${QUAC_INCLUDE_DIRS} ${XACC_QUAC_BASE_INCLUDE_DIR} ./tpls/cereal/include)
target_link_libraries(${MPI_EXECUTOR} PRIVATE FunctorExecutor)

find_package(Python COMPONENTS Interpreter Development)
if(Python_FOUND)
    add_subdirectory(python)
else()
    message(WARNING "Cannot find Python. Skip building Python bindings.")
endif()