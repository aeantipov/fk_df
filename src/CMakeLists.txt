include_directories (.)

set (fk_src 
#    Logger.h Logger.cpp
    GFWrap.h GFWrap.cpp
    Solver.h Solver.cpp
    DMFT.h DMFT.hpp DMFT.cpp
    Diagrams.h Diagrams.hpp
    DF.h DF.hpp 
    #DF.cpp
    )
add_library(fk SHARED ${fk_src})
set_target_properties(fk PROPERTIES LINKER_LANGUAGE CXX)

#target_link_libraries(fk
#    ${Boost_LIBRARIES}
#    ${PYTHON_LIBRARIES}
#    ${HDF5_LIBRARIES}
#)

install(TARGETS fk
    LIBRARY DESTINATION lib
    )

install(DIRECTORY . DESTINATION include/fk
          FILES_MATCHING PATTERN "*.h *.hpp")
