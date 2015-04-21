# ================================================ +
# PyroTools -- A toolkit of computational programs +
# for sequencing data analysis                     +
#                                                  +
# src/                                             +
# ================================================ +

# ================================================ +
# Link to BamTools library                         +
# ================================================ +
# set Bamtools include and library path
include_directories(${PyroTools_SOURCE_DIR}/bamtools/include)
link_directories(${PyroTools_SOURCE_DIR}/bamtools/lib)

# ================================================ +
# Link to NlOpt library                            +
# ================================================ +
# set Nlopt include and library path
include_directories(${PyroTools_SOURCE_DIR}/nlopt/include)
link_directories(${PyroTools_SOURCE_DIR}/nlopt/lib)

# ================================================ +
# Compile setting for PyroTools                    +
# ================================================ +

# set external link libraries
set(LIBS m z pthread bamtools bamtools-utils nlopt gomp)

# source file
aux_source_directory(. SRC_LIST)

# compile main PyroTools application
add_executable(
               PyroTools_cmd
               ${SRC_LIST}
              )


# set PyroTools application properties
set_target_properties(PyroTools_cmd PROPERTIES
                      OUTPUT_NAME "PyroTools")


# define libraries to link
target_link_libraries(PyroTools_cmd ${LIBS})

# set applification install destination
install(TARGETS PyroTools_cmd DESTINATION "bin")
