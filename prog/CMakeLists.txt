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
    set (ex fk_DMFT_${dmft_lattice})
    add_executable(${ex} ${ex_src} ${other_headers})
    set_target_properties(${ex} PROPERTIES COMPILE_DEFINITIONS "LATTICE_${dmft_lattice}") 
    target_link_libraries(${ex} fk ${FFTW_LIBRARIES})
    install ( TARGETS ${ex} DESTINATION bin )
endforeach(dmft_lattice)

#DF
foreach (df_lattice ${df_lattices})
    message(STATUS "DF: ${df_lattice}")
    set (ex_src fk_DF.cpp)
    set (ex fk_DF_${df_lattice})
    add_executable(${ex} ${ex_src} ${other_headers})
    set_target_properties(${ex} PROPERTIES COMPILE_DEFINITIONS "LATTICE_${df_lattice}") 
    target_link_libraries(${ex} fk ${FFTW_LIBRARIES})
    install ( TARGETS ${ex} DESTINATION bin )
endforeach(df_lattice)


foreach (ex ${execs})
    set(ex_src ${ex}.cpp)
    add_executable(${ex} ${ex_src} ${other_headers})
    target_link_libraries(${ex} fk ${FFTW_LIBRARIES})
    set_target_properties(${ex} PROPERTIES COTIRE_PREFIX_HEADER_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/src" )
    install ( TARGETS ${ex} DESTINATION bin )
endforeach(ex)

add_executable(dmft_cubic_path EXCLUDE_FROM_ALL DMFTCubicSuscPath.cpp) 
target_link_libraries(dmft_cubic_path fk)

