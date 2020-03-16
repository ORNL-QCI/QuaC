set (LIBRARY_NAME _pyquaC)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")

if(APPLE)
   set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
endif(APPLE)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing -O2 -g -pipe -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -Wformat -fexceptions --param=ssp-buffer-size=4 -grecord-gcc-switches -m64 -mtune=generic -D_GNU_SOURCE -fPIC -fwrapv")
add_library(${LIBRARY_NAME} SHARED py_systemModel.cpp)
target_include_directories(${LIBRARY_NAME} PUBLIC ${Python_INCLUDE_DIRS} ${XACC_ROOT}/include/pybind11/include)

set_target_properties(${LIBRARY_NAME} PROPERTIES PREFIX "")
target_link_libraries(${LIBRARY_NAME} PUBLIC xacc_quac)

if(APPLE)
    set_target_properties(${LIBRARY_NAME} PROPERTIES INSTALL_RPATH "@loader_path/lib")
    set_target_properties(${LIBRARY_NAME} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
else()
    set_target_properties(${LIBRARY_NAME} PROPERTIES INSTALL_RPATH "$ORIGIN/lib")
    set_target_properties(${LIBRARY_NAME} PROPERTIES LINK_FLAGS "-shared")
endif()

install(FILES $<TARGET_FILE:${LIBRARY_NAME}> DESTINATION ${CMAKE_INSTALL_PREFIX}) 


