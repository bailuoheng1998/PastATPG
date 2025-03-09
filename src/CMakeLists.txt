#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
#set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

# add ictest src to a static library
file(GLOB_RECURSE SRC CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_library(ictest STATIC ${SRC})

#add_executable(ictest-shell ${SRC})
target_include_directories(ictest PUBLIC ${PROJECT_SOURCE_DIR}/src/include)
target_link_libraries(ictest PRIVATE ZLIB::ZLIB)
#target_include_directories(ictest-shell PUBLIC ${PROJECT_SOURCE_DIR}/src/include)


# ##########################################
# set(OpenGL_GL_PREFERENCE LEGACY)
# # FIND GLUT
# find_package(GLUT REQUIRED)
# include_directories(${GLUT_INCLUDE_DIRS})
# link_directories(${GLUT_LIBRARY_DIRS})
# add_definitions(${GLUT_DEFINITIONS})
# if(NOT GLUT_FOUND)
#     message(ERROR " GLUT not found!")
# endif(NOT GLUT_FOUND)
# # FIND OPENGL
# find_package(OpenGL REQUIRED)
# include_directories(${OpenGL_INCLUDE_DIRS})
# link_directories(${OpenGL_LIBRARY_DIRS})
# add_definitions(${OpenGL_DEFINITIONS})
# if(NOT OPENGL_FOUND)
#     message(ERROR " OPENGL not found!")
# endif(NOT OPENGL_FOUND)
# ##########################################
 target_link_libraries(ictest PUBLIC OpenMP::OpenMP_CXX)
 target_link_libraries(ictest PUBLIC spdlog_header_only)

# target_link_libraries(ictest PUBLIC ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} -ldl -lglut)


target_include_directories(ictest PUBLIC ${PROJECT_SOURCE_DIR}/thirdparty/boost/linux/include)
target_include_directories(ictest PUBLIC ${PROJECT_SOURCE_DIR}/thirdparty/CLI11)
target_include_directories(ictest PUBLIC ${PROJECT_SOURCE_DIR}/thirdparty/prompt)
target_include_directories(ictest PUBLIC ${PROJECT_SOURCE_DIR}/thirdparty/nlohmann)

#target_link_libraries(ictest-shell PUBLIC spdlog_header_only)
#target_include_directories(ictest-shell PUBLIC ${PROJECT_SOURCE_DIR}/thirdparty/boost/linux/include)
#target_include_directories(ictest-shell PUBLIC ${PROJECT_SOURCE_DIR}/thirdparty/CLI11)
#target_include_directories(ictest-shell PUBLIC ${PROJECT_SOURCE_DIR}/thirdparty/prompt)
