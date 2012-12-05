include_directories (. ../jsoncpp/include ../src)

link_libraries(
fk
)

enable_testing()

# Here all the tests are set. The source is file is assumed to be ${test}.cpp
set (tests 
GridTest
ContainerTest
GridObjectTest
SolverTest
)

foreach (test ${tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    add_test(${test} ${test})
endforeach(test)
