
set(RT_SOURCES ${CMAKE_CURRENT_LIST_DIR}/src/stb_write_impl.cpp)
set(RT_INCLUDES ${CMAKE_CURRENT_LIST_DIR}/include)
message(${CMAKE_CURRENT_LIST_DIR})
add_executable(RT RT.cpp)
target_include_directories(RT PUBLIC ${RT_INCLUDES})
target_sources(RT PUBLIC ${RT_SOURCES})