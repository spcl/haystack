# define source files
set (SOURCES 
   model_options.c
   CacheModel.cpp
   Program.cpp
   Access.cpp
   Timer.cpp
   isl-helpers.cpp)

add_library(CacheModel OBJECT ${SOURCES})
add_executable(StackDistance main.cpp $<TARGET_OBJECTS:CacheModel>)

llvm_map_components_to_libnames(LLVM_LIBRARIES core option scalaropts)
set(CLANG_LIBRARIES 
   clangFrontendTool 
   clangFrontend 
   clangDriver 
   clangFrontend 
   clangParse 
   clangSema 
   clangEdit 
   clangAnalysis 
   clangAST 
   clangLex 
   clangBasic 
   clangDriver 
   clangSerialization)

target_link_libraries(StackDistance 
   barvinok.a pet.a isl.a polylibgmp.a ntl.a m.a gmp.a
   ${Boost_LIBRARIES} 
   ${LLVM_LIBRARIES} 
   ${CLANG_LIBRARIES} 
   -lpthread -lm) 
