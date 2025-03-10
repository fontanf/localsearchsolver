add_library(LocalSearchSolver_localsearchsolver)
target_sources(LocalSearchSolver_localsearchsolver PRIVATE
    common.cpp)
target_include_directories(LocalSearchSolver_localsearchsolver PUBLIC
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(LocalSearchSolver_localsearchsolver PUBLIC
    OptimizationTools::containers
    OptimizationTools::utils)
add_library(LocalSearchSolver::localsearchsolver ALIAS LocalSearchSolver_localsearchsolver)
set_target_properties(LocalSearchSolver_localsearchsolver PROPERTIES OUTPUT_NAME "localsearchsolver")
install(TARGETS LocalSearchSolver_localsearchsolver)
