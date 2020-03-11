find_package(XACC REQUIRED)
include_directories(${XACC_INCLUDE_DIRS})
link_libraries(-Wl,-rpath,${XACC_ROOT}/lib -L${XACC_ROOT}/lib)

add_executable(one_qubit one_qubit.cpp)
target_link_libraries(one_qubit PRIVATE xacc_quac xacc::quantum_gate)

add_executable(cross_resonant_interaction cross_resonant_interaction.cpp)
target_link_libraries(cross_resonant_interaction PRIVATE xacc_quac xacc::quantum_gate)

add_executable(pulse_optimization pulse_optimization.cpp)
target_link_libraries(pulse_optimization PRIVATE xacc_quac xacc::quantum_gate)

add_executable(quantum_process_tomography quantum_process_tomography.cpp)
target_link_libraries(quantum_process_tomography PRIVATE xacc_quac xacc::quantum_gate)

add_executable(pulse_command_def pulse_command_def.cpp)
target_link_libraries(pulse_command_def PRIVATE xacc_quac xacc::quantum_gate)
target_compile_definitions(pulse_command_def PRIVATE GATEIR_TEST_FILE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/resources")

add_executable(hamiltonian_json_input hamiltonian_json_input.cpp)
target_link_libraries(hamiltonian_json_input PRIVATE xacc_quac xacc::quantum_gate)
target_compile_definitions(hamiltonian_json_input PRIVATE GATEIR_TEST_FILE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/resources")


add_executable(goat_optimized_pulse goat_optimized_pulse.cpp)
target_link_libraries(goat_optimized_pulse PRIVATE xacc_quac xacc::quantum_gate)

add_executable(freq_sweep freq_sweep.cpp)
target_link_libraries(freq_sweep PRIVATE xacc_quac xacc::quantum_gate)
target_compile_definitions(freq_sweep PRIVATE GATEIR_TEST_FILE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/resources")

add_executable(fourier_pulse fourier_pulse.cpp)
target_link_libraries(fourier_pulse PRIVATE xacc_quac xacc::quantum_gate)

add_executable(mpi_demo mpi_demo.cpp)
target_link_libraries(mpi_demo PRIVATE xacc_quac xacc::quantum_gate)
target_compile_definitions(mpi_demo PRIVATE GATEIR_TEST_FILE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/resources")
