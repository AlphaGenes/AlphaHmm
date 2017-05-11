SET(ALPHAHMM_src
AlphaHmmInMod.f90
CMakeLists.txt
hmm.f90
hmm.f90.bak
hmmHaplotype.f90
hmmHeader.f90
hmmvariables.f90
main.f90
par_zig_mod.f90
utils.f90

)



if (AlphaHMMEXE)
    add_definitions(-D COMMIT=${VERSION})
    ADD_EXECUTABLE(${AlphaHMMEXE} ${ALPHAHMM_src})
    # add external links here
    TARGET_LINK_LIBRARIES(${AlphaHMMEXE} ${AHLIB})
endif()
