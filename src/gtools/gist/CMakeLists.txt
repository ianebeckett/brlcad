set(gist_ignore_files
  CMakeLists.txt
  CONTRIBUTING.md
  FactsHandler.h
  IFPainter.h
  InformationGatherer.h
  LICENSE
  Options.h
  PerspectiveGatherer.h
  Position.h
  README.md
  RenderHandler.h
  TODO
  pch.h
  picohash.h
)

set(GIST_SRCS
  FactsHandler.cpp
  IFPainter.cpp
  InformationGatherer.cpp
  main.cpp
  Options.cpp
  PerspectiveGatherer.cpp
  Position.cpp
  RenderHandler.cpp
)

CMAKEFILES(${gist_ignore_files} ${GIST_SRCS})

find_package(OpenCV COMPONENTS core features2d imgproc highgui)
# FIND_PACKAGE(OpenCV)
if(NOT OpenCV_FOUND)
  message("OpenCV is DISABLED")
  return()
endif(NOT OpenCV_FOUND)
message("Found OpenCV: ${OpenCV_LIBRARIES}")

set(GIST_INCLUDE_DIRS
  ${BU_INCLUDE_DIRS}
  ${RT_INCLUDE_DIRS}
  ${GED_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
)
list(REMOVE_DUPLICATES GIST_INCLUDE_DIRS)
include_directories(${GIST_INCLUDE_DIRS})

BRLCAD_ADDEXEC(gist "${GIST_SRCS}"  "libbu;librt;libged;${OpenCV_LIBRARIES};")

# FIXME: should be cxx_std_14, but uses std::filesystem.  should use libbu.
target_compile_features(gist PRIVATE cxx_std_17)

if (CMAKE_SHARED_LIBRARY_SUFFIX STREQUAL ".dll")
  add_custom_command(TARGET gist POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:gist> $<TARGET_FILE_DIR:gist>
    COMMAND_EXPAND_LISTS
  )
endif (CMAKE_SHARED_LIBRARY_SUFFIX STREQUAL ".dll")


# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8