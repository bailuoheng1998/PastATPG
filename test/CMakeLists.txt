# include gtest(no need)
#include_directories(${PROJECT_SOURCE_DIR}/thirdparty/googletest-release-1.12.1/googletest/include)

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

# add uint test function
function(ictest_test test_file)
    # get file name
    get_filename_component(test_target_name ${test_file} NAME_WE)
    add_executable(${test_target_name} ${test_file})
    target_link_libraries(${test_target_name} PUBLIC ictest gtest gtest_main)
    add_test(NAME ${test_target_name} COMMAND ${test_target_name})
endfunction(ictest_test)

ictest_test(TestSATATPG_ND_backtrack.cpp)
