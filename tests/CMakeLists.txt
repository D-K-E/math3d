# declare test executable 


include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.14.0
)
FetchContent_MakeAvailable(googletest)
add_library(GTest::GTest INTERFACE IMPORTED)
target_link_libraries(GTest::GTest INTERFACE gtest_main)

add_executable(
    ${PROJECT_NAME}_Test
    "test_matn.cpp"
    "test_vecn.cpp"
    "test_quaternion.cpp"
    "test_lu.cpp"
    )

target_link_libraries(${PROJECT_NAME}_Test PRIVATE GTest::GTest ${PROJECT_NAME})
