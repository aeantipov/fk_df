include_directories (. ../src)

#message(STATUS "Will compile DMFT for the following lattices:" ${dmft_lattices} )
#message(STATUS "Will compile DF for the following lattices:" ${df_lattices} )

set (execs 
${prog_targets}
)

file(GLOB_RECURSE other_headers ../src/*.h ../src/*.hpp)

#DMFT
foreach (dmft_lattice ${dmft_lattices})
    message(STATUS "DMFT: ${dmft_lattice}")
    set (ex_src fk_DMFT.cpp)
    set (ex fk_DMFT_${dmft_lattice}${FKPOSTFIX})
    add_executable(${ex} ${ex_src} ${other_headers})
    set_target_properties(${ex} PROPERTIES COMPILE_DEFINITIONS "LATTICE_${dmft_lattice}") 
    target_link_libraries(${ex} fk${FKPOSTFIX} ${FFTW_LIBRARIES})
    install ( TARGETS ${ex} DESTINATION bin )
endforeach(dmft_lattice)

#DF
foreach (df_lattice ${df_lattices})
    message(STATUS "DF: ${df_lattice}")
    set (ex_src fk_DF.cpp)
    set (ex fk_DF_${df_lattice}${FKPOSTFIX})
    add_executable(${ex} ${ex_src} ${other_headers})
    set_target_properties(${ex} PROPERTIES COMPILE_DEFINITIONS "LATTICE_${df_lattice}") 
    target_link_libraries(${ex} fk${FKPOSTFIX} ${FFTW_LIBRARIES})
    install ( TARGETS ${ex} DESTINATION bin )
endforeach(df_lattice)


foreach (ex_name ${execs})
    set(ex_src ${ex_name}.cpp)
    set(ex ${ex_name}${FKPOSTFIX})
    add_executable(${ex} ${ex_src} ${other_headers})
    target_link_libraries(${ex} fk${FKPOSTFIX} ${FFTW_LIBRARIES})
    install ( TARGETS ${ex} DESTINATION bin )
endforeach(ex_name)

add_executable(dmft_cubic_path EXCLUDE_FROM_ALL DMFTCubicSuscPath.cpp) 
target_link_libraries(dmft_cubic_path fk)

