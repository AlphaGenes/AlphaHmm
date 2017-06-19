module ExternalHMMWrappers


    implicit none
    contains

    subroutine AlphaImputeHMMRunner(inputParams, ped, probGenosOut, probPhaseOut, genosCountsOut, fullHOut)
        use AlphaHmmInMod
        use PedigreeModule
        use hmmModule
        use GlobalVariablesHmmMaCH

        type(pedigreeHolder), intent(in), target :: ped
        type(AlphaHmmInput), intent(in)  :: inputParams


        !< all outputs are given in format of genotyped individuals
        real,allocatable,dimension (:,:), intent(out) :: probGenosOut
        real,allocatable,dimension (:,:,:), intent(out) :: probPhaseOut
        integer, allocatable, dimension(:,:,:), intent(out) :: genosCountsOut
        integer(kind=1),allocatable,dimension(:,:,:), intent(out) :: fullHOut !< full phase and imputation information from HMM (nanis, snps, haps[2])

        imputeGenosHMM = ped%getGenotypesAsArray()
        imputePhaseHmm = ped%getPhaseAsArray()
        pedigree => ped

        defaultInput = inputParams

        call HMMController(inputParams%HMMOption) 


        ! imputeGenosOut = ImputeGenosHMM


        genosCountsOut = GenosCounts
        probGenosOut = ProbImputeGenosHmm
        probPhaseOut = ProbImputePhaseHmm
        fullHOut = fullH

    end subroutine AlphaImputeHMMRunner

end module ExternalHMMWrappers