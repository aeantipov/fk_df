include_directories (.)

set (fk_src 
    EigenIterator.h
    FKCommon.h
    Logger.h Logger.cpp
    Grid.h 
    Container.h Container.hpp
    GridObject.h GridObject.hpp
    GFWrap.h 
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
          FILES_MATCHING PATTERN "*.h")
