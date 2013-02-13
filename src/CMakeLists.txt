include_directories (.)

set (base_src
    FKCommon.h
    Grid.h Grid.hpp 
    Container.h Container.hpp
    GridObject.h GridObject.hpp
    )

add_library(fk_base OBJECT ${base_src})

set (fk_src 
    Logger.h Logger.cpp
    GFWrap.h GFWrap.cpp
    Solver.h Solver.cpp
    SelfConsistency.h SelfConsistency.hpp SelfConsistency.cpp
    Diagrams.h Diagrams.hpp
    DF.h DF.hpp DF.cpp
    )
add_library(fk SHARED ${fk_src} $<TARGET_OBJECTS:fk_base>)
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
