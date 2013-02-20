include_directories (. ../src)

set (execs 
fk_DF
fk_DMFT
)

file(GLOB_RECURSE other_headers ../src/*.h ../src/*.hpp)

foreach (ex ${execs})
    set(ex_src ${ex}.cpp)
    add_executable(${ex} ${ex_src} ${other_headers})
    target_link_libraries(${ex} fk)
    set_target_properties(${ex} PROPERTIES COTIRE_PREFIX_HEADER_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/src" )
    install ( TARGETS ${ex} DESTINATION bin )
endforeach(ex)

add_executable(dmft_cubic_path EXCLUDE_FROM_ALL DMFTCubicSuscPath.cpp) 
target_link_libraries(dmft_cubic_path fk)

