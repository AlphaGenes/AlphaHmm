!-----------------------------------------------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-----------------------------------------------------------------------------------------------------------------------
!
! MODULE: HmmInputMod
!
!> @file        HmmOutputMod.f90
!
! DESCRIPTION:
!> @brief       Module holding writing subroutines
!>
!> @details     This MODULE contains a class which contains all subroutines to write out a file.
!
!> @author      Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date        May 11, 2017
!
!> @version     0.0.1 (alpha)
!
! REVISION HISTORY:
! 2017.05.11  Rantolin - Initial Version
!
!-----------------------------------------------------------------------------------------------------------------------

module HmmOutputMod
    use iso_fortran_env

    private
    ! public :: WriteOutProbabilities, WriteOutGenotypes, WriteOutDosage

contains

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Write out results
    !
    !> @details    This subroutine encapsulates the different the subroutines
    !>             that write out the final results of the HMM
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date       May 11, 2017
    !
    ! PARAMETERS:
    !> @param[in]
    !---------------------------------------------------------------------------
    subroutine WriteOutput(Output)
    use iso_fortran_env
    implicit none

    integer, intent(in) :: Output

    if (btest(Output,0))
        call WriteOutGenotypes
    endif
    if (btest(Output,1))
        call WriteOutDosage
    endif
    if (btest(Output,2))
        call WriteOutProbabilities
    endif

    end subroutine WriteOutput

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Write out most likely genotypes
    !
    !> @details    Write out most likely genotypes
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date       May 11, 2017
    !
    ! PARAMETERS:
    !> @param[in]
    !---------------------------------------------------------------------------
    subroutine WriteOutGenotypes(nAnisG)
        implicit none

        integer, intent(in) :: nAnisG

        integer :: i
        integer :: GenosFileUnit
        logical :: opened, named
        character(len=300) :: GenosFile="Results/ImputedGenotypes.txt"

        inquire(unit=GenosFileUnit, opened=opened, named=named, name=GenosFile)
        if (.NOT. opened .and. named) then
            open(unit=GenosFileUnit, file=GenosFile, status='unknown')
        else if (.NOT. named) then
            write(0, *) "ERROR - Something went wrong when trying to read the file of pre-phased data"
            open(newunit=GenosFileUnit, file=GenosFile)
        end if

        do i = 1, nAnisG
            ! write (,'(a20,240000i0)') pedigree%pedigree(hmmID)%originalID,ImputeGenos(i,:)
        end do

        close(GenosFileUnit)

    end subroutine WriteOutGenotypes

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Write out genotype dosages
    !
    !> @details    Write out genotype dosages
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date       May 11, 2017
    !
    ! PARAMETERS:
    !> @param[in]
    !---------------------------------------------------------------------------
    subroutine WriteOutDosage(nAnisG)
        implicit none

        integer, intent(in) :: nAnisG

        integer :: i
        integer :: GenosFileUnit
        logical :: opened, named
        character(len=300) :: DosageFile="Results/ImputedDosages.txt"

        inquire(unit=DosageFileUnit, opened=opened, named=named, name=DosageFile)
        if (.NOT. opened .and. named) then
            open(unit=DosageFileUnit, file=DosageFile, status='unknown')
        else if (.NOT. named) then
            write(0, *) "ERROR - Something went wrong when trying to read the file of pre-phased data"
            open(newunit=DosageFileUnit, file=DosageFile)
        end if

        do i = 1, nAnisG
            ! write (,'(a20,240000i0)') pedigree%pedigree(hmmID)%originalID,ImputeGenos(i,:)
        end do

        close(DosageFileUnit)

    end subroutine WriteOutDosage

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Write out genotype probabilities
    !
    !> @details    Write out genotype probabilities
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date       May 11, 2017
    !
    ! PARAMETERS:
    !> @param[in]
    !---------------------------------------------------------------------------
    subroutine WriteOutProbabilities(nAnisG)
        implicit none

        integer, intent(in) :: nAnisG

        integer :: i
        integer :: ProbsFileUnit
        logical :: opened, named
        character(len=300) :: ProbsFile="Results/ImputedDosages.txt"

        inquire(unit=ProbsFileUnit, opened=opened, named=named, name=ProbsFile)
        if (.NOT. opened .and. named) then
            open(unit=ProbsFileUnit, file=ProbsFile, status='unknown')
        else if (.NOT. named) then
            write(0, *) "ERROR - Something went wrong when trying to read the file of pre-phased data"
            open(newunit=ProbsFileUnit, file=ProbsFile)
        end if

        do i = 1, nAnisG
            ! write (,'(a20,240000i0)') pedigree%pedigree(hmmID)%originalID,ImputeGenos(i,:)
            ! write (,'(a20,240000i0)') pedigree%pedigree(hmmID)%originalID,ImputeGenos(i,:)
        end do

        close(ProbsFileUnit)

    end subroutine WriteOutProbabilities

end module HmmOutputMod




