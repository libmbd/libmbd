if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    execute_process(
        COMMAND git describe --tags --dirty=.dirty
        OUTPUT_VARIABLE VERSION_TAG
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "Setting version tag to ${VERSION_TAG} from Git")
else()
    include(LibmbdVersionTag)
endif()

set(PROJECT_VERSION ${VERSION_TAG})
string(REGEX MATCH "^([0-9]+)\.([0-9]+)\.([0-9]+)-?(.*)?$" VERSION_TAG ${VERSION_TAG})

set(PROJECT_VERSION_MAJOR ${CMAKE_MATCH_1})
set(PROJECT_VERSION_MINOR ${CMAKE_MATCH_2})
set(PROJECT_VERSION_PATCH ${CMAKE_MATCH_3})
set(PROJECT_VERSION_SUFFIX ${CMAKE_MATCH_4})
