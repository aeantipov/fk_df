include_directories (. ../src)

link_libraries(
    ${Boost_LIBRARIES}
    ${PYTHON_LIBRARIES}
    ${HDF5_LIBRARIES}
    fk
)

set (execs 
fk_DF
#fk_DMFT_cubic
#fk_DMFT_bethe
#fk_DMFT_cubicinf
)

foreach (ex ${execs})
    set(ex_src ${ex}.cpp)
    add_executable(${ex} ${ex_src})
    install ( TARGETS ${ex} DESTINATION bin )
endforeach(ex)

