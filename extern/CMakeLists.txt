# Enable FetchContent.
include(FetchContent)

# Fetch fontanf/orproblems.
FetchContent_Declare(
    orproblems
    GIT_REPOSITORY https://github.com/fontanf/orproblems.git
    GIT_TAG ba099d6bd0ac186e4113fe340789a1800fa72852)
    #SOURCE_DIR "${PROJECT_SOURCE_DIR}/../orproblems/")
FetchContent_MakeAvailable(orproblems)
