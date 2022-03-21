# CMake generated Testfile for 
# Source directory: /home/xsarm06/Development/cpp_projects/lattice_pricer/unit_tests
# Build directory: /home/xsarm06/Development/cpp_projects/lattice_pricer/build/unit_tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(lattice_printing_t "lattice_printing_test" "COMMAND" "lattice_printing_test")
set_tests_properties(lattice_printing_t PROPERTIES  _BACKTRACE_TRIPLES "/home/xsarm06/Development/cpp_projects/lattice_pricer/unit_tests/CMakeLists.txt;19;add_test;/home/xsarm06/Development/cpp_projects/lattice_pricer/unit_tests/CMakeLists.txt;0;")
add_test(lattice_t "lattice_test" "COMMAND" "lattice_test")
set_tests_properties(lattice_t PROPERTIES  _BACKTRACE_TRIPLES "/home/xsarm06/Development/cpp_projects/lattice_pricer/unit_tests/CMakeLists.txt;20;add_test;/home/xsarm06/Development/cpp_projects/lattice_pricer/unit_tests/CMakeLists.txt;0;")
add_test(lattice_forward_t "lattice_forward_test" "COMMAND" "lattice_forward_test")
set_tests_properties(lattice_forward_t PROPERTIES  _BACKTRACE_TRIPLES "/home/xsarm06/Development/cpp_projects/lattice_pricer/unit_tests/CMakeLists.txt;21;add_test;/home/xsarm06/Development/cpp_projects/lattice_pricer/unit_tests/CMakeLists.txt;0;")
