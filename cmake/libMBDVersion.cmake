if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    find_package(Git REQUIRED)

    execute_process(
        COMMAND git describe --tags --dirty=.dirty
        RESULTS_VARIABLE RUN_TAG_EXTRACTION_FAIL
        OUTPUT_VARIABLE VERSION_TAG
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    # when tag based description fails, use hash
    if(RUN_TAG_EXTRACTION_FAIL)
        execute_process(
            COMMAND git describe --tags --dirty=.dirty --always
            OUTPUT_VARIABLE VERSION_TAG
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()

    message(STATUS "Setting version tag to ${VERSION_TAG} from Git")
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cmake/libMBDVersionTag.cmake")
    include(libMBDVersionTag)
else()
    message(FATAL_ERROR
        "Not in a Git repository and version tag is missing, you most likely "
        "attempt to install from a copy of the source tree. Obtain the source "
        "distribution (libmbd-<version>.tar.gz) from a Github release page "
        "instead.")
endif()

set(PROJECT_VERSION ${VERSION_TAG})

if(RUN_TAG_EXTRACTION_FAIL)
    set(PROJECT_VERSION_MAJOR 0)
    set(PROJECT_VERSION_MINOR 0)
    set(PROJECT_VERSION_PATCH 0)
    set(PROJECT_VERSION_SUFFIX ${VERSION_TAG})
else()
    string(REGEX MATCH "^([0-9]+)\.([0-9]+)\.([0-9]+)-?(.*)?$" VERSION_TAG ${VERSION_TAG})

    set(PROJECT_VERSION_MAJOR ${CMAKE_MATCH_1})
    set(PROJECT_VERSION_MINOR ${CMAKE_MATCH_2})
    set(PROJECT_VERSION_PATCH ${CMAKE_MATCH_3})
    set(PROJECT_VERSION_SUFFIX ${CMAKE_MATCH_4})
endif()
