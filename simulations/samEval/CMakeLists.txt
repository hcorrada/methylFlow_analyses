include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/samEvaluation
  ${CMAKE_BINARY_DIR}/samEvaluation
)

add_executable(samEvaluation
  main.cpp
)

TARGET_LINK_LIBRARIES(samEvaluation 
 mflib ${LEMON_LIBRARIES} glpk
) 

INSTALL(
  TARGETS samEvaluation
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)


