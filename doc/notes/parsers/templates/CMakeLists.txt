set(READER foo)

PERPLEX_TARGET(${READER}_scanner ${READER}_scanner.perplex
    ${CMAKE_CURRENT_BINARY_DIR}/${READER}_scanner.c
    ${CMAKE_CURRENT_BINARY_DIR}/${READER}_scanner.h)
LEMON_TARGET(${READER}_parser ${READER}_parser.lemon
    ${CMAKE_CURRENT_BINARY_DIR}/${READER}_parser.c
    ${CMAKE_CURRENT_BINARY_DIR}/${READER}_parser.h)
ADD_PERPLEX_LEMON_DEPENDENCY(${READER}_scanner ${READER}_parser)

set (SOURCES
  ${PERPLEX_${PARSER}_scanner_OUTPUTS}
  ${LEMON_${PARSER}_parser_OUTPUT_SOURCE}
  )

set (LOCAL_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
)

BRLCAD_ADDLIB(lib${READER} "${SOURCES}" "librt;libbn;libbu" "" "${LOCAL_INCLUDE_DIRS}")

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
