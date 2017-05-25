module ExternalHMMWrappers


    implicit none
    contains

    subroutine AlphaImputeHMMRunner(inputParams, ImputeGenos, ImputePhase, ped, probGenosOut, probPhaseOut, genosCountsOut, fullHOut)
        use AlphaHmmInMod
        use PedigreeModule
        use hmmModule
        use GlobalVariablesHmmMaCH

        type(pedigreeHolder), intent(in), target :: ped
        type(AlphaHmmInput), intent(in)  :: inputParams
        integer(kind=1),allocatable,dimension (:,:),target, intent(in) :: ImputeGenos
        integer(kind=1),allocatable,dimension (:,:,:), target, intent(in) :: ImputePhase

        real,allocatable,dimension (:,:), intent(out) :: probGenosOut
        real,allocatable,dimension (:,:,:), intent(out) :: probPhaseOut
        integer, allocatable, dimension(:,:,:), intent(out) :: genosCountsOut
        integer(kind=1),allocatable,dimension(:,:,:), intent(out) :: fullHOut !< full phase and imputation information from HMM (nanis, snps, haps[2])

        imputeGenosHMM = ImputeGenos
        imputePhaseHmm = ImputePhase
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