# 1. Reference pure .h lib
# ---------------------------------------------------------------------------------------
# add spdlog library
# Header only version without .cpp
# ---------------------------------------------------------------------------------------
find_package(Threads REQUIRED)
add_library(spdlog_header_only INTERFACE)
add_library(spdlog::spdlog_header_only ALIAS spdlog_header_only)
target_include_directories(spdlog_header_only INTERFACE ${PROJECT_SOURCE_DIR}/thirdparty)
target_link_libraries(spdlog_header_only INTERFACE Threads::Threads)

# 2. Reference thirdparty lib by src
# ---------------------------------------------------------------------------------------
# add googletest library
# ---------------------------------------------------------------------------------------
# adds the targers: gtest, gtest_main, gmock, gmock_main
add_subdirectory(${PROJECT_SOURCE_DIR}/thirdparty/googletest-release-1.12.1)


