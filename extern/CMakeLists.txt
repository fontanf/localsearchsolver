# Enable FetchContent.
include(FetchContent)

# Fetch fontanf/orproblems.
FetchContent_Declare(
    orproblems
    GIT_REPOSITORY https://github.com/fontanf/orproblems.git
    GIT_TAG b9e5b32ea612fd2cf23e71de0a969b6dd56ae5d3)
    #SOURCE_DIR "${PROJECT_SOURCE_DIR}/../orproblems/")
FetchContent_MakeAvailable(orproblems)
