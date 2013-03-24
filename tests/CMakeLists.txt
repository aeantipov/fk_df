include_directories (. ../jsoncpp/include ../src)

link_libraries(
fk
)

enable_testing()

# Here all the tests are set. The source is file is assumed to be ${test}.cpp
set (tests 
SolverTest
DMFTCubicTest
DMFTBetheTest
DMFTCubicSuscTest
DMFTCubicSuscTest02
DMFT_DFCubicTest01
DMFT_DFCubicTest02
DMFT_DFCubicTest03
DiagramsTest
DFCubicTest01
)

foreach (test ${tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    add_test(${test} ${test})
endforeach(test)
