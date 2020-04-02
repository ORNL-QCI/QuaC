include_directories(${CMAKE_BINARY_DIR})
add_xacc_test(Simple)
target_link_libraries(SimpleTester xacc xacc_quac xacc::quantum_gate)

add_xacc_test(OpenPulse)
target_link_libraries(OpenPulseTester xacc xacc_quac xacc::quantum_gate)
target_compile_definitions(OpenPulseTester PRIVATE GATEIR_TEST_FILE_DIR="${CMAKE_SOURCE_DIR}/xacc_examples/resources")

add_xacc_test(Circuit)
target_link_libraries(CircuitTester xacc xacc_quac xacc::quantum_gate)

add_xacc_test(DigitalGate)
target_link_libraries(DigitalGateTester xacc xacc_quac xacc::quantum_gate)

if (ENABLE_MPI_TESTS)
    # Macro to add an MPI-enabled gtest: just change the driver command to starting with 1 process.
    macro(add_mpi_test _TEST_NAME)
       add_executable(${_TEST_NAME}Tester ${_TEST_NAME}Tester.cpp)
       target_include_directories(${_TEST_NAME}Tester PRIVATE ${GTEST_INCLUDE_DIRS})
       add_test(NAME xacc_${_TEST_NAME}Tester COMMAND sh -c "${MPIEXEC_PROGRAM} -np 1 ./${_TEST_NAME}Tester")
       target_link_libraries(${_TEST_NAME}Tester gtest)
    endmacro()
   
   # Add the a simple MPI test
   add_mpi_test(QuaCMPI)
   target_link_libraries(QuaCMPITester xacc xacc_quac xacc::quantum_gate)
   target_compile_definitions(QuaCMPITester PRIVATE GATEIR_TEST_FILE_DIR="${CMAKE_SOURCE_DIR}/xacc_examples/resources")
endif()

add_xacc_test(MultiLevel)
target_link_libraries(MultiLevelTester xacc xacc_quac xacc::quantum_gate)

add_xacc_test(IrTransform)
target_link_libraries(IrTransformTester xacc xacc_quac xacc::quantum_gate)


add_xacc_test(DefaultBackend)
target_link_libraries(DefaultBackendTester xacc xacc_quac xacc::quantum_gate)