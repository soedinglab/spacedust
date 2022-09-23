include_directories(commons)
add_subdirectory(commons)
add_subdirectory(util)
add_subdirectory(version)
add_subdirectory(workflow)

add_executable(clustersearch
        ${commons_source_files}
        ${util_header_files}
        ${util_source_files}
        ${workflow_source_files}
        LocalCommandDeclarations.h
        clustersearchBase.cpp
        )
mmseqs_setup_derived_target(clustersearch)
target_link_libraries(clustersearch version)

add_dependencies(clustersearch local-generated)

install(TARGETS clustersearch DESTINATION bin)
