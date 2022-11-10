include_directories(commons)
add_subdirectory(commons)
add_subdirectory(util)
add_subdirectory(version)
add_subdirectory(workflow)

add_executable(spacedust
        ${commons_source_files}
        ${util_header_files}
        ${util_source_files}
        ${workflow_source_files}
        LocalCommandDeclarations.h
        spacedust.cpp
        )
mmseqs_setup_derived_target(spacedust)
target_link_libraries(spacedust version)

add_dependencies(spacedust local-generated)

install(TARGETS spacedust DESTINATION bin)
