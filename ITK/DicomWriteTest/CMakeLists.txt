cmake_minimum_required(VERSION 2.8)
 
project(TestDicom)
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(DicomSeriesReadSeriesWrite DicomSeriesReadSeriesWrite.cxx )
target_link_libraries(DicomSeriesReadSeriesWrite  ${ITK_LIBRARIES})

