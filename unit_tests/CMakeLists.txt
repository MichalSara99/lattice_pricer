find_package(Boost REQUIRED COMPONENTS unit_test_framework)

include_directories(${PROJECT_SOURCE_DIR}/include)

set(TEST_SOURCES ${PROJECT_SOURCE_DIR}/unit_tests)
set(TEST_LIBS lattice_pricer_lib)



add_executable(lattice_printing_test ${TEST_SOURCES}/lattice_printing_t)
add_executable(lattice_test ${TEST_SOURCES}/lattice_t)
add_executable(lattice_forward_test ${TEST_SOURCES}/lattice_forward_t)

target_link_libraries(lattice_printing_test ${TEST_LIBS} ${Boost_LIBRARIES})
target_link_libraries(lattice_test ${TEST_LIBS} ${Boost_LIBRARIES})
target_link_libraries(lattice_forward_test ${TEST_LIBS} ${Boost_LIBRARIES})


add_test(lattice_printing_t lattice_printing_test COMMAND lattice_printing_test)
add_test(lattice_t lattice_test COMMAND lattice_test)
add_test(lattice_forward_t lattice_forward_test COMMAND lattice_forward_test)