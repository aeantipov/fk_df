include_directories (.)

set (fk_src 
    fk_common.h
    gf_wrap.h gf_wrap.cpp
    )
add_library(fk SHARED ${fk_src})

target_link_libraries(fk
    ${Boost_LIBRARIES}
#    ${PYTHON_LIBRARIES}
#    ${HDF5_LIBRARIES}
)

install(TARGETS fk
    LIBRARY DESTINATION lib
    )

install(DIRECTORY . DESTINATION include/fk
          FILES_MATCHING PATTERN "*.h")
