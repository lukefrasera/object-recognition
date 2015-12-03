CMAKE_MINIMUM_REQUIRED(VERSION 2.8.3)
project(object-recognition)

INCLUDE_DIRECTORIES(include)
option(LOG_INFO_PRINTING "Add logging to builds executable for debugging")
if( ${LOG_INFO_PRINTING} )
  set(LOG_INFO_PRINTING 1)
  add_definitions(-DLOG_INFO_PRINTING=${LOG_INFO_PRINTING})
endif( ${LOG_INFO_PRINTING} )


add_subdirectory(utils/)
