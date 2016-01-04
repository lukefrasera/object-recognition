CMAKE_MINIMUM_REQUIRED(VERSION 2.8.3)
project(object-recognition)

INCLUDE_DIRECTORIES(include)
option(LOG_INFO_PRINTING "Add logging to builds executable for debugging")
if( ${LOG_INFO_PRINTING} )
  set(LOG_INFO_PRINTING 1)
  add_definitions(-DLOG_INFO_PRINTING=${LOG_INFO_PRINTING})
endif( ${LOG_INFO_PRINTING} )

find_package(OpenCV REQUIRED)
option(BUILD_TESTS "Add test build functionality")

add_library(quad_tree
  src/quad_tree.cc
)

if( ${BUILD_TESTS} )
  add_executable(qr_test
    src/qr_test.cc
    src/qr_detect.cc
  )
  target_link_libraries(qr_test
    ${OpenCV_LIBS}
    quad_tree
  )
endif( ${BUILD_TESTS} )

option(BUILD_UTILS "Build utility programs" True)
if( ${BUILD_UTILS} )
  add_subdirectory(utils/)
endif( ${BUILD_UTILS} )
