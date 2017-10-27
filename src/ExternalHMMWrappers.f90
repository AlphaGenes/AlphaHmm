!-----------------------------------------------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-----------------------------------------------------------------------------------------------------------------------
!
! MODULE:       ExternalHMMWrappers
!
!> @file        ExternalHMMWrappers.f90
!
! DESCRIPTION:
!> @brief       Module holding a wrapper for the HMM in AlphaImpute

!>
!> @details     This file contains a wrapper to translate the data structure from AlphaImpute into that of the HMM
!
!> @author      David Wilson, david.wilson@roslin.ed.ac.uk
!
!> @date        Oct 1, 2017
!
!> @version     0.0.1 (alpha)
!
! REVISION HISTORY:
! 2017.10.01  Dwilson - Initial Version
!
!-----------------------------------------------------------------------------------------------------------------------
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
        imputePhaseHMM = ped%getPhaseAsArray()
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