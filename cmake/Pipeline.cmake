# SPDX-License-Identifier: GPL-3.0-or-later
add_executable(ttauto_pipeline_export "${CMAKE_SOURCE_DIR}/native/export.cpp")
target_link_libraries(ttauto_pipeline_export PRIVATE ttauto)
ttauto_add_compile_flags(ttauto_pipeline_export)
ttauto_set_inplace_output(ttauto_pipeline_export "${CMAKE_BINARY_DIR}/bin")
if(MINGW)
  target_link_options(ttauto_pipeline_export PRIVATE -static)
endif()
add_executable(ttauto_trace_prefix "${CMAKE_SOURCE_DIR}/native/trace_prefix.cpp")
target_compile_options(ttauto_trace_prefix PRIVATE -O3 -fno-fast-math)
ttauto_set_inplace_output(ttauto_trace_prefix "${CMAKE_BINARY_DIR}/bin")
if(MINGW)
  target_link_options(ttauto_trace_prefix PRIVATE -static)
endif()
if(BUILD_TESTING)
  find_package(Python3 3.10 COMPONENTS Interpreter REQUIRED)
  add_test(NAME pipeline_python COMMAND "${Python3_EXECUTABLE}" -m unittest discover
    -s "${CMAKE_SOURCE_DIR}/python/tests" -p "test_*.py")
  set_tests_properties(pipeline_python PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/python;TTAUTO_BUILD=${CMAKE_BINARY_DIR}"
    LABELS "pipeline;deterministic" TIMEOUT 240)
endif()
