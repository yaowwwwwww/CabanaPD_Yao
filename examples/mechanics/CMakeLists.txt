add_executable(ElasticWave elastic_wave.cpp)
target_link_libraries(ElasticWave LINK_PUBLIC CabanaPD)

add_executable(KalthoffWinkler kalthoff_winkler.cpp)
target_link_libraries(KalthoffWinkler LINK_PUBLIC CabanaPD)

add_executable(CrackBranching crack_branching.cpp)
target_link_libraries(CrackBranching LINK_PUBLIC CabanaPD)

add_executable(FragmentingCylinder fragmenting_cylinder.cpp)
target_link_libraries(FragmentingCylinder LINK_PUBLIC CabanaPD)

add_executable(RandomCracks random_cracks.cpp)
target_link_libraries(RandomCracks LINK_PUBLIC CabanaPD)

add_executable(DogboneTensileTest dogbone_tensile_test.cpp)
target_link_libraries(DogboneTensileTest LINK_PUBLIC CabanaPD)

add_executable(PlateWithHole plate_with_hole.cpp)
target_link_libraries(PlateWithHole LINK_PUBLIC CabanaPD)

add_executable(ColdSprayImpact simple_impact.cpp)  
target_link_libraries(ColdSprayImpact LINK_PUBLIC CabanaPD)

install(TARGETS ElasticWave KalthoffWinkler CrackBranching FragmentingCylinder RandomCracks DogboneTensileTest PlateWithHole ColdSprayImpact DESTINATION ${CMAKE_INSTALL_BINDIR})
