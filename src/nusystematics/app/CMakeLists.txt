LIST(APPEND TARGETS_TO_BUILD 
GenerateSystProviderConfigNuSyst
DumpConfiguredTweaksNuSyst
PlotSystVariationsNuSyst
)

foreach(targ ${TARGETS_TO_BUILD})
  add_executable(${targ} ${targ}.cxx)
  target_link_libraries(${targ} 
    nusyst::systproviders
    nusyst::utility
    ${ROOT_LIBRARIES}
    ROOT::Geom
    ROOT::MathMore
  )

endforeach()

install(TARGETS ${TARGETS_TO_BUILD} DESTINATION bin/)
