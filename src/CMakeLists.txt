add_subdirectory(ompl)
include_directories(${PYTHON_INCLUDE_DIRS} ${TRAC_IK_INCLUDE_DIRS})

find_package(trac_ik_lib REQUIRED)
include_directories(${trac_ik_lib_INCLUDE_DIRS})

find_package(nlopt REQUIRED)
include_directories(${nlopt_INCLUDE_DIRS})

install(DIRECTORY ompl/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ompl"
    COMPONENT ompl
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN ".DS_Store" EXCLUDE
    PATTERN "sample_contrib" EXCLUDE
    REGEX "/src$" EXCLUDE
    REGEX "/doc$" EXCLUDE
    REGEX "/tests$" EXCLUDE)     