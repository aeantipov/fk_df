include_directories (. ../src)

link_libraries(
    ${Boost_LIBRARIES}
    ${PYTHON_LIBRARIES}
    ${HDF5_LIBRARIES}
    FK
)

add_executable(fk_DF fk_DF.cpp)
install ( TARGETS fk_DF DESTINATION bin )
