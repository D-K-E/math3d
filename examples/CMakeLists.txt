
add_executable(vectorExample "vector.cpp")
target_link_libraries(vectorExample PRIVATE ${PROJECT_NAME})


add_executable(matrixExample "matrix.cpp")
target_link_libraries(matrixExample PRIVATE ${PROJECT_NAME})
