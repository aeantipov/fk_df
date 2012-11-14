include_directories (.)

set (FK_src 
    FK.h
    gf_wrap.h gf_wrap.cpp
    )
add_library(FK SHARED ${FK_src})

target_link_libraries(FK
    ${Boost_LIBRARIES}
#    ${PYTHON_LIBRARIES}
#    ${HDF5_LIBRARIES}
)

install(TARGETS FK
    LIBRARY DESTINATION lib
    )

install(DIRECTORY . DESTINATION include/FK
          FILES_MATCHING PATTERN "*.h")
