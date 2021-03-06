# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/fengzeng/workspace/PyroTools/bamtools

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/fengzeng/workspace/PyroTools/bamtools/build

# Utility rule file for UtilsHeaders.

# Include the progress variables for this target.
include src/utils/CMakeFiles/UtilsHeaders.dir/progress.make

src/utils/CMakeFiles/UtilsHeaders:
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/fengzeng/workspace/PyroTools/bamtools/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Exporting UtilsHeaders"

UtilsHeaders: src/utils/CMakeFiles/UtilsHeaders
UtilsHeaders: src/utils/CMakeFiles/UtilsHeaders.dir/build.make
	cd /Users/fengzeng/workspace/PyroTools/bamtools/build/src/utils && /opt/local/bin/cmake -E copy_if_different /Users/fengzeng/workspace/PyroTools/bamtools/src/utils/bamtools_fasta.h /Users/fengzeng/workspace/PyroTools/bamtools/include/utils/bamtools_fasta.h
	cd /Users/fengzeng/workspace/PyroTools/bamtools/build/src/utils && /opt/local/bin/cmake -E copy_if_different /Users/fengzeng/workspace/PyroTools/bamtools/src/utils/bamtools_filter_engine.h /Users/fengzeng/workspace/PyroTools/bamtools/include/utils/bamtools_filter_engine.h
	cd /Users/fengzeng/workspace/PyroTools/bamtools/build/src/utils && /opt/local/bin/cmake -E copy_if_different /Users/fengzeng/workspace/PyroTools/bamtools/src/utils/bamtools_filter_properties.h /Users/fengzeng/workspace/PyroTools/bamtools/include/utils/bamtools_filter_properties.h
	cd /Users/fengzeng/workspace/PyroTools/bamtools/build/src/utils && /opt/local/bin/cmake -E copy_if_different /Users/fengzeng/workspace/PyroTools/bamtools/src/utils/bamtools_filter_ruleparser.h /Users/fengzeng/workspace/PyroTools/bamtools/include/utils/bamtools_filter_ruleparser.h
	cd /Users/fengzeng/workspace/PyroTools/bamtools/build/src/utils && /opt/local/bin/cmake -E copy_if_different /Users/fengzeng/workspace/PyroTools/bamtools/src/utils/bamtools_options.h /Users/fengzeng/workspace/PyroTools/bamtools/include/utils/bamtools_options.h
	cd /Users/fengzeng/workspace/PyroTools/bamtools/build/src/utils && /opt/local/bin/cmake -E copy_if_different /Users/fengzeng/workspace/PyroTools/bamtools/src/utils/bamtools_pileup_engine.h /Users/fengzeng/workspace/PyroTools/bamtools/include/utils/bamtools_pileup_engine.h
	cd /Users/fengzeng/workspace/PyroTools/bamtools/build/src/utils && /opt/local/bin/cmake -E copy_if_different /Users/fengzeng/workspace/PyroTools/bamtools/src/utils/bamtools_utilities.h /Users/fengzeng/workspace/PyroTools/bamtools/include/utils/bamtools_utilities.h
	cd /Users/fengzeng/workspace/PyroTools/bamtools/build/src/utils && /opt/local/bin/cmake -E copy_if_different /Users/fengzeng/workspace/PyroTools/bamtools/src/utils/bamtools_variant.h /Users/fengzeng/workspace/PyroTools/bamtools/include/utils/bamtools_variant.h
	cd /Users/fengzeng/workspace/PyroTools/bamtools/build/src/utils && /opt/local/bin/cmake -E copy_if_different /Users/fengzeng/workspace/PyroTools/bamtools/src/utils/utils_global.h /Users/fengzeng/workspace/PyroTools/bamtools/include/utils/utils_global.h
.PHONY : UtilsHeaders

# Rule to build all files generated by this target.
src/utils/CMakeFiles/UtilsHeaders.dir/build: UtilsHeaders
.PHONY : src/utils/CMakeFiles/UtilsHeaders.dir/build

src/utils/CMakeFiles/UtilsHeaders.dir/clean:
	cd /Users/fengzeng/workspace/PyroTools/bamtools/build/src/utils && $(CMAKE_COMMAND) -P CMakeFiles/UtilsHeaders.dir/cmake_clean.cmake
.PHONY : src/utils/CMakeFiles/UtilsHeaders.dir/clean

src/utils/CMakeFiles/UtilsHeaders.dir/depend:
	cd /Users/fengzeng/workspace/PyroTools/bamtools/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/fengzeng/workspace/PyroTools/bamtools /Users/fengzeng/workspace/PyroTools/bamtools/src/utils /Users/fengzeng/workspace/PyroTools/bamtools/build /Users/fengzeng/workspace/PyroTools/bamtools/build/src/utils /Users/fengzeng/workspace/PyroTools/bamtools/build/src/utils/CMakeFiles/UtilsHeaders.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/utils/CMakeFiles/UtilsHeaders.dir/depend

