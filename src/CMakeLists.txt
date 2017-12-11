SET(ALPHAHMM_src
AlphaHmmInMod.f90
hmm.f90
hmmInputMod.f90
hmmHaplotype.f90
hmmHeader.f90
hmmvariables.f90
main.f90
ExternalHMMWrappers.f90
HMMInputOutputModule.f90
par_zig_mod.f90
utils.f90
)

if (AlphaHMMEXE)
    add_definitions(-D COMMIT=${VERSION})
    ADD_EXECUTABLE(${AlphaHMMEXE} ${ALPHAHMM_src})
    # add external links here
    TARGET_LINK_LIBRARIES(${AlphaHMMEXE} ${AHLIB})
endif()


if (AHMMLIB)
    ADD_LIBRARY(${AHMMLIB} STATIC  ${ALPHAHMM_src})
    TARGET_LINK_LIBRARIES(${AHMMLIB} ${AHLIB})
endif()
