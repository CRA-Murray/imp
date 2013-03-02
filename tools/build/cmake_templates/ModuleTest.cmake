# Autogenerated file, run tools/build/setup_cmake.py to regenerate

include_directories(%(includepath)s)
link_directories(%(libpath)s)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${%(NAME)s_CXX_FLAGS}")


File(GLOB runtimepytests "${PROJECT_BINARY_DIR}/test/%(name)s/test_*.py")
set(pytests %(pytests)s %(expytests)s)

foreach (test ${runtimepyttests} ${pytests})
  GET_FILENAME_COMPONENT(name ${test} NAME_WE)
  if(EXISTS "${PROJECT_BINARY_DIR}/test/%(name)s/${name}.pytests")
    FILE(READ "${PROJECT_BINARY_DIR}/test/%(name)s/${name}.pytests" contents)
    STRING(REGEX REPLACE ";" "\\\\;" contents "${contents}")
    STRING(REGEX REPLACE "\n" ";" contents "${contents}")
    foreach(testline ${contents})
      string(REGEX REPLACE "([A-Za-z0-9_]+\\.[A-Za-z0-9_]+) (.*)" 
                           "\\1;\\2" split "${testline}")
      list(GET split 0 methname)
      list(GET split 1 docstring)
      add_test("%(name)s.${name}.${methname}" ${PROJECT_BINARY_DIR}/setup_environment.sh ${IMP_PYTHON} ${test} "${methname}")
      set_tests_properties("%(name)s.${name}.${methname}" PROPERTIES LABELS "%(name)s;test")
      set_tests_properties("%(name)s.${name}.${methname}" PROPERTIES MEASUREMENT "docstring=${docstring}")
    endforeach()
  else()
    add_test("%(name)s.${name}" ${PROJECT_BINARY_DIR}/setup_environment.sh ${IMP_PYTHON} ${test})
    set_tests_properties("%(name)s.${name}" PROPERTIES LABELS "%(name)s;test")
  endif()
endforeach(test)

set(cpp_tests %(cpptests)s %(excpptests)s)

foreach (test ${cpp_tests})
   GET_FILENAME_COMPONENT(name ${test} NAME_WE)
   add_executable("%(name)s.${name}" ${test})
   target_link_libraries("%(name)s.${name}"     imp_%(name)s
    %(modules)s
    %(dependencies)s)
   set_target_properties("%(name)s.${name}" PROPERTIES
                         RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/test/%(name)s/"
                         OUTPUT_NAME ${name})
   add_test("%(name)s.${name}" ${PROJECT_BINARY_DIR}/setup_environment.sh
            "${PROJECT_BINARY_DIR}/test/%(name)s/${name}")
   set_tests_properties("%(name)s.${name}" PROPERTIES LABELS "%(name)s;test")
   set(executables ${executables} "%(name)s.${name}")
endforeach(test)

add_custom_target("imp_%(name)s_tests" ALL DEPENDS ${executables})
