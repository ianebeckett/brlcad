# Define a macro for building lists of files.  Distcheck needs
# to know what files are "supposed" to be present in order to make
# sure the source tree is clean prior to building a distribution
# tarball, hence this macro stores its results in files and not
# variables  It's a no-op in a SUBBUILD.
MACRO(CMAKEFILES)
	IF(NOT BRLCAD_IS_SUBBUILD)
		FOREACH(ITEM ${ARGN})
			IF(NOT ${ITEM} MATCHES "^SHARED$" AND NOT ${ITEM} MATCHES "^STATIC$" AND NOT x${ITEM} MATCHES "^xWIN32$")
				GET_FILENAME_COMPONENT(ITEM_PATH ${ITEM} PATH)
				GET_FILENAME_COMPONENT(ITEM_NAME ${ITEM} NAME)
				IF(NOT ITEM_PATH STREQUAL "")
					GET_FILENAME_COMPONENT(ITEM_ABS_PATH ${ITEM_PATH} ABSOLUTE)
					IF(NOT ${ITEM_PATH} MATCHES "^${ITEM_ABS_PATH}$")
						IF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ITEM})
							MESSAGE(FATAL_ERROR "Attempting to ignore non-existent file ${ITEM}, in directory ${CMAKE_CURRENT_SOURCE_DIR}")
						ENDIF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ITEM})
						IF(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ITEM})
							FILE(APPEND	${CMAKE_BINARY_DIR}/cmakedirs.cmake	"${ITEM_ABS_PATH}/${ITEM_NAME}\n")
						ELSE(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ITEM})
							FILE(APPEND	${CMAKE_BINARY_DIR}/cmakefiles.cmake "${ITEM_ABS_PATH}/${ITEM_NAME}\n")
						ENDIF(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ITEM})
						FILE(APPEND	${CMAKE_BINARY_DIR}/cmakefiles.cmake "${ITEM_ABS_PATH}\n")
						WHILE(NOT ITEM_PATH STREQUAL "")
							get_filename_component(ITEM_NAME	${ITEM_PATH} NAME)
							get_filename_component(ITEM_PATH	${ITEM_PATH} PATH)
							IF(NOT ITEM_PATH STREQUAL "")
								GET_FILENAME_COMPONENT(ITEM_ABS_PATH ${ITEM_PATH} ABSOLUTE)
								IF(NOT ${ITEM_NAME} MATCHES "\\.\\.")
									FILE(APPEND	${CMAKE_BINARY_DIR}/cmakefiles.cmake "${ITEM_ABS_PATH}\n")
								ENDIF(NOT ${ITEM_NAME} MATCHES "\\.\\.")
							ENDIF(NOT ITEM_PATH STREQUAL "")
						ENDWHILE(NOT ITEM_PATH STREQUAL "")
					ENDIF(NOT ${ITEM_PATH} MATCHES "^${ITEM_ABS_PATH}$")
				ELSE(NOT ITEM_PATH STREQUAL "")
					IF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ITEM})
						MESSAGE(FATAL_ERROR "Attempting to ignore non-existent file ${ITEM}, in directory ${CMAKE_CURRENT_SOURCE_DIR}")
					ENDIF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ITEM})
					IF(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ITEM})
						FILE(APPEND	${CMAKE_BINARY_DIR}/cmakedirs.cmake "${CMAKE_CURRENT_SOURCE_DIR}/${ITEM_NAME}\n")
					ELSE(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ITEM})
						FILE(APPEND	${CMAKE_BINARY_DIR}/cmakefiles.cmake "${CMAKE_CURRENT_SOURCE_DIR}/${ITEM_NAME}\n")
					ENDIF(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ITEM})
				ENDIF(NOT ITEM_PATH STREQUAL "")
			ENDIF()
		ENDFOREACH(ITEM ${ARGN})
	ENDIF(NOT BRLCAD_IS_SUBBUILD)
ENDMACRO(CMAKEFILES FILESLIST)

