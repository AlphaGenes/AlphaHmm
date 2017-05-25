module HMMInputOutputModule


contains 

    subroutine WriteProbabilitiesForHMM(outFile, Indexes, nAnims, nSnps)
        use GlobalVariablesHmmMaCH
        use AlphaHmmInMod
        character(len=*), intent(IN) :: outFile
        integer, intent(IN) :: nAnims, nSnps
        integer, intent(IN) :: Indexes(:)
        type(AlphaHmmInput), pointer :: inputParams
        ! character*(20), intent(IN) :: Ids(:)
        ! Local variables
        integer :: i,j, n0, n1, n2
        real :: d
        real, allocatable :: Probs0(:), Probs1(:)

        inputParams => defaultInput
        open (unit=55,file=outFile,status="unknown")

        allocate(Probs0(nSnps))
        allocate(Probs1(nSnps))

        n0=0
        n1=0
        n2=0
        d=0.0

        do i=1,nAnims
            do j=1,nSnps
                n1 = GenosCounts(i,j,1)                           ! Heterozygous
                n2 = GenosCounts(i,j,2)                           ! Homozygous: 2 case
                n0 = (GlobalRoundHmm-inputParams%HmmBurnInRound) - n1 - n2     ! Homozygous: 0 case
                d = n0 + n1 + n2
                Probs0(j)=n0/d
                Probs1(j)=n1/d
            enddo
            write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') pedigree%Pedigree(indexes(i))%originalID,Probs0(:)
            write (55,'(a20,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2,20000f5.2)') pedigree%Pedigree(indexes(i))%originalID,Probs1(:)
            n0=0
            n1=0
            n2=0
            d=0.0
        enddo

        close(55)
        ! write(0,*) "Wrote out file", outFile, "with genotype probabilities"
    end subroutine WriteProbabilitiesForHMM

end module HMMInputOutputModule
