add_library(iceicle_io pvd_writer.cpp)
target_link_libraries(iceicle_io PUBLIC iceicle_fe)
target_include_directories(iceicle_io PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
