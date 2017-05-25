MODULE Utils
    use GlobalVariablesHmmMaCH
    IMPLICIT NONE

CONTAINS

    !######################################################################
    SUBROUTINE RemoveInformationGameteSegment(gamete,nStart,nStop)
        ! Remove the allele information of a gamete segment



        INTEGER,INTENT(IN) :: nStart, nStop
        INTEGER(KIND=1),INTENT(OUT),DIMENSION(:) :: gamete

        gamete(nStart:nStop) = ALLELE_MISSING

    END SUBROUTINE RemoveInformationGameteSegment

    !######################################################################
    SUBROUTINE RemoveInformationGamete(gamete)
        ! Remove the allele information of a gamete



        INTEGER(KIND=1),INTENT(OUT),DIMENSION(:) :: gamete

        gamete(:) = ALLELE_MISSING

    END SUBROUTINE RemoveInformationGamete

    !######################################################################
    SUBROUTINE RemoveAlleleInformationIndividual(ToWhom)
        ! Remove the genotype information of an individual

        use AlphaHmmInMod

        INTEGER, INTENT(IN) :: ToWhom
        type(AlphaHmmInput), pointer :: inputParams

        inputParams => defaultInput


        call RemoveAlleleInformationSegment(ToWhom, 1, inputParams%nsnp)

    END SUBROUTINE RemoveAlleleInformationIndividual

    !######################################################################
    SUBROUTINE RemoveAlleleInformationSegment(ToWhom,nStart,nStop)
        ! Remove the genotype information of an individual



        INTEGER, INTENT(IN) :: ToWhom, nStart, nStop

        call RemoveInformationGamete(PhaseHmmMaCH(ToWhom,nStart:nStop,1))
        call RemoveInformationGamete(PhaseHmmMaCH(ToWhom,nStart:nStop,2))

    END SUBROUTINE RemoveAlleleInformationSegment

    !######################################################################
    SUBROUTINE RemoveGenotypeInformationIndividual(ToWhom)
        ! Remove the genotype information of an individual
         use AlphaHmmInMod

        type(AlphaHmmInput), pointer :: inputParams
        INTEGER, INTENT(IN) :: ToWhom

        inputParams => defaultInput

        call RemoveGenotypeInformationIndividualSegment(ToWhom,1,inputParams%nsnp)

    END SUBROUTINE RemoveGenotypeInformationIndividual

    !######################################################################
    SUBROUTINE RemoveGenotypeInformationIndividualSegment(ToWhom,nStart,nStop)
        ! Remove the genotype information of an individual



        INTEGER, INTENT(IN) :: ToWhom,nStart,nStop

        GenosHmmMaCH(ToWhom,nStart:nStop) = MISSING

    END SUBROUTINE RemoveGenotypeInformationIndividualSegment

    !######################################################################
    FUNCTION CountPhasedGametes RESULT( gametesPhased )

        use GlobalVariablesHmmMaCH

        use AlphaHmmInMod

        INTEGER :: gametesPhased, ind
        INTEGER(KIND=1), ALLOCATABLE :: gamete(:)
        type(AlphaHmmInput), pointer :: inputParams

        inputParams => defaultInput

        allocate(gamete(inputParams%nsnp))
        gametesPhased=0
        do ind=1,pedigree%pedigreeSize-pedigree%ndummys
            gamete=imputePhaseHmm(ind,:,1)
            if (float(count(gamete(:)==1 .OR. gamete(:)==0))/inputParams%nsnp >= inputParams%imputedThreshold/100.0) then
                gametesPhased=gametesPhased+1
            endif

            gamete=imputePhaseHmm(ind,:,2)
            if (float(count(gamete(:)==1 .OR. gamete(:)==0))/inputParams%nsnp >= inputParams%imputedThreshold/100.0) then
                gametesPhased=gametesPhased+1
            endif
        enddo
        RETURN
    END FUNCTION CountPhasedGametes

    !######################################################################
    FUNCTION CountPhasedAlleles(gamete) RESULT( allelesPhased )

        use GlobalVariablesHmmMaCH


        INTEGER :: allelesPhased
        INTEGER(KIND=1),INTENT(IN) :: gamete(:)

        allelesPhased = count(gamete(:)==1 .OR. gamete(:)==0)


        RETURN
    END FUNCTION CountPhasedAlleles

    !######################################################################
    FUNCTION CountMissingAlleles(gamete) RESULT( allelesMissing )

        use GlobalVariablesHmmMaCH
        use AlphaHmmInMod

    

        INTEGER :: allelesMissing
        INTEGER(KIND=1),INTENT(IN) :: gamete(:)
        type(AlphaHmmInput), pointer :: inputParams
        inputParams => defaultInput


        allelesMissing = inputParams%nsnp - CountPhasedAlleles(gamete)

        RETURN
    END FUNCTION CountMissingAlleles

    !######################################################################
    FUNCTION CountMissingAllelesByGametes(gamete1,gamete2) RESULT( allelesMissing )
        ! TODO: THIS SUBROUTINE IS NOT WORKING PROPERLY. ALLELESMISSING GETS SHIT

        use GlobalVariablesHmmMaCH
        use AlphaHmmInMod

        INTEGER :: allelesMissing
        INTEGER(KIND=1),INTENT(IN) :: gamete1(:),gamete2(:)
        type(AlphaHmmInput), pointer :: inputParams
        inputParams => defaultInput

        allelesMissing = inputParams%nsnp - CountGenotypedAllelesByGametes(gamete1,gamete2)

        RETURN
    END FUNCTION CountMissingAllelesByGametes

    !######################################################################
    FUNCTION CountGenotypedAllelesByGametes(gamete1, gamete2) RESULT( allelesGenotyped )

        use GlobalVariablesHmmMaCH


        INTEGER(KIND=1),INTENT(IN) :: gamete1(:), gamete2(:)

        ! Local variables
        INTEGER :: allelesGenotyped


        allelesGenotyped = count((gamete1(:)==1 .OR. gamete1(:)==0) .AND. (gamete2(:)==1 .OR. gamete2(:)==0))

        ! do i=1,size(gamete1)
        !     if ((gamete1(i)==1 .OR. gamete1(i)==0) .AND. (gamete2(i)==1 .OR. gamete2(i)==0)) then
        !         allelesGenotyped = allelesGenotyped + 1
        !     endif
        ! enddo

        RETURN
    END FUNCTION CountGenotypedAllelesByGametes

    !######################################################################
    FUNCTION CountGenotypedGenotypesByChromosome(chromosome) RESULT( genotypesMissing )

        use GlobalVariablesHmmMaCH


        INTEGER,INTENT(IN) :: chromosome(:)

        ! Local variables
        INTEGER :: genotypesMissing

        genotypesMissing = count(chromosome(:)==0 .or. chromosome(:)==2 .or. chromosome(:)==1)

        RETURN
    END FUNCTION CountGenotypedGenotypesByChromosome

    !######################################################################
    FUNCTION CountMissingdGenotypesByChromosome(chromosome) RESULT( genotypesMissing )

        use GlobalVariablesHmmMaCH
        use AlphaHmmInMod

        INTEGER,INTENT(IN) :: chromosome(:)

        ! Local variables
        INTEGER :: genotypesMissing
        type(AlphaHmmInput), pointer :: inputParams
        inputParams => defaultInput

        genotypesMissing = inputParams%nsnp - CountGenotypedGenotypesByChromosome(chromosome)

        RETURN
    END FUNCTION CountMissingdGenotypesByChromosome

    !######################################################################

END MODULE
