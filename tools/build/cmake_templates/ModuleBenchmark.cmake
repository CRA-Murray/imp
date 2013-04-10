# Autogenerated file, run tools/build/setup_cmake.py to regenerate

include_directories(%(includepath)s)
link_directories(%(libpath)s)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${%(NAME)s_CXX_FLAGS}")

math(EXPR timeout "${IMP_TIMEOUT_FACTOR} * 300")

if (NOT ${IMP_MAX_CHECKS} MATCHES "NONE")
set(testarg "--run_quick_test")
else()
set(testarg "")
endif()

set(cppbenchmarks %(cppbenchmarks)s)

foreach (bin ${cppbenchmarks})
   GET_FILENAME_COMPONENT(name ${bin} NAME_WE)
   add_executable(%(name)s.${name} ${bin})
   target_link_libraries(%(name)s.${name}     imp_%(name)s
    %(modules)s
    ${IMP_BENCHMARK_LIBRARY}
    %(dependencies)s)
   set_target_properties(%(name)s.${name} PROPERTIES
                         RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/benchmark/%(name)s"
                         OUTPUT_NAME ${name})
   add_test(%(name)s.${name} ${IMP_TEST_SETUP}
            "${PROJECT_BINARY_DIR}/benchmark/%(name)s/${name}${CMAKE_EXECUTABLE_SUFFIX}" "${testarg}")
   set_tests_properties("%(name)s.${name}" PROPERTIES LABELS "IMP.%(name)s;benchmark")
   set_tests_properties("%(name)s.${name}" PROPERTIES TIMEOUT ${timeout})
   set_tests_properties("%(name)s.${name}" PROPERTIES COST 1)
   set(executables ${executables} %(name)s.${name})
endforeach(bin)

add_custom_target("imp_%(name)s_benchmarks" ALL DEPENDS ${executables}
  # add dummy dep as empty targets seem to go away
  imp_%(name)s imp_base)

set(pybenchmarks %(pybenchmarks)s)
foreach (test ${pybenchmarks})
 GET_FILENAME_COMPONENT(name ${test} NAME_WE)
 add_test("%(name)s.${name}" ${IMP_TEST_SETUP} python ${test})
 set_tests_properties("%(name)s.${name}" PROPERTIES LABELS "IMP.%(name)s;benchmark")
 set_tests_properties("%(name)s.${name}" PROPERTIES TIMEOUT ${timeout})
 set_tests_properties("%(name)s.${name}" PROPERTIES COST 4)
endforeach(test)
