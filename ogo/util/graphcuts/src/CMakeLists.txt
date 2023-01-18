# Sources and headers --------------------------------------------
set (SHEET_SRCS sheetness.cxx)

# Build, link, install
add_executable(Sheetness ${SHEET_SRCS})
target_link_libraries(Sheetness ${ITK_LIBRARIES})
install (TARGETS Sheetness RUNTIME DESTINATION bin)

# Sources and headers --------------------------------------------
set (PERI_SRCS periosteal_segmentation.cxx)

# Build, link, install
add_executable(PeriostealSegmentation ${PERI_SRCS})
target_link_libraries(PeriostealSegmentation ${ITK_LIBRARIES})
install (TARGETS PeriostealSegmentation RUNTIME DESTINATION bin)

# Sources and headers --------------------------------------------
set (ENDO_SRCS endosteal_segmentation.cxx)

# Build, link, install
add_executable(EndostealSegmentation ${ENDO_SRCS})
target_link_libraries(EndostealSegmentation ${ITK_LIBRARIES})
install (TARGETS EndostealSegmentation RUNTIME DESTINATION bin)

# Sources and headers --------------------------------------------
set (HU_SRCS hu_periosteal_segmentation.cxx)

# Build, link, install
add_executable(HU ${HU_SRCS})
target_link_libraries(HU ${ITK_LIBRARIES})
install (TARGETS HU RUNTIME DESTINATION bin)

# Sources and headers --------------------------------------------
set (SHEET2_SRCS sheetness2.cxx)

# Build, link, install
add_executable(Sheetness2 ${SHEET2_SRCS})
target_link_libraries(Sheetness2 ${ITK_LIBRARIES})
install (TARGETS Sheetness2 RUNTIME DESTINATION bin)

