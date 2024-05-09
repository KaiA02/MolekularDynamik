# make doc_doxygen optional if someone does not have / like doxygen

option(BUILD_DOCS "Build Doxygen documentation" ON)


# Doxygen custom target
if(BUILD_DOCS)
    find_package(Doxygen REQUIRED)

    # Set input and output directories
    set(DOXYGEN_INPUT_DIR "${CMAKE_SOURCE_DIR}" CACHE PATH "Path to input directory for Doxygen")
    set(DOXYGEN_OUTPUT_DIR "${CMAKE_BINARY_DIR}/docs" CACHE PATH "Path to output directory for Doxygen")

    # Specify existing Doxyfile
    set(DOXYGEN_CONFIG_FILE "${CMAKE_BINARY_DIR}/Doxyfile")
    configure_file(${CMAKE_SOURCE_DIR}/Doxyfile ${DOXYGEN_CONFIG_FILE} @ONLY)

    # Add custom target for building the documentation
    add_custom_target(doc
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONFIG_FILE}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM
    )


    # Add custom target for cleaning documentation
    add_custom_target(clean-doc
            COMMAND ${CMAKE_COMMAND} -E remove_directory ${DOXYGEN_OUTPUT_DIR}
            COMMENT "Removing generated documentation"
            VERBATIM
    )

    # Exclude the Doxygen target from the all/default target
    set_target_properties(doc PROPERTIES EXCLUDE_FROM_ALL TRUE)

    # Add option to disable creating the Doxygen target
    option(SKIP_DOXYGEN_TARGET "Skip creating the Doxygen target" OFF)
    if(SKIP_DOXYGEN_TARGET)
        set_target_properties(doc PROPERTIES EXCLUDE_FROM_DEFAULT_BUILD TRUE)
    endif()
endif()