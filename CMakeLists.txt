# stop if cmake version is below 3.0
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# project name and enable C++ support
project(nlofit CXX)

set(CMAKE_BUILD_TYPE Release)

# Git module
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_describe(GIT_DESCRIBE)
git_local_changes(GIT_LOCAL_CHANGES)

# Headers and libraries for AmplitudeLib v2
# and Minuit and Cuba (if not available on system)
link_directories("../amplitudelib_v2/build/lib/")
link_directories("../libs/lib/")
include_directories("../amplitudelib_v2/")
include_directories("../libs/include")

add_subdirectory(src)

set_target_properties(nlofit PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )
