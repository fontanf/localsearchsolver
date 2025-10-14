# Enable FetchContent.
include(FetchContent)

# Fetch fontanf/optimizationtools.
set(OPTIMIZATIONTOOLS_BUILD_TEST OFF)
FetchContent_Declare(
    optimizationtools
    GIT_REPOSITORY https://github.com/fontanf/optimizationtools.git
    GIT_TAG e8c379203792fcad764cce239986b325619475fd
    #SOURCE_DIR "${PROJECT_SOURCE_DIR}/../optimizationtools/"
    EXCLUDE_FROM_ALL)
FetchContent_MakeAvailable(optimizationtools)

# Fetch fontanf/orproblems.
if(LOCALSEARCHSOLVER_BUILD_EXAMPLES)
    FetchContent_Declare(
        orproblems
        GIT_REPOSITORY https://github.com/fontanf/orproblems.git
        GIT_TAG ba099d6bd0ac186e4113fe340789a1800fa72852)
        #SOURCE_DIR "${PROJECT_SOURCE_DIR}/../orproblems/")
    FetchContent_MakeAvailable(orproblems)
endif()
