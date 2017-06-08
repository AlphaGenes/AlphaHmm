!-----------------------------------------------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-----------------------------------------------------------------------------------------------------------------------
!
! MODULE: HmmInputMod
!
!> @file        HmmInputMod.f90
!
! DESCRIPTION:
!> @brief       Module holding reading subroutines
!>
!> @details     This MODULE contains a class which contains all subroutines to read in a file.
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

module HmmInputMod
    use iso_fortran_env

    private
    public :: CountInData, ReadInData

contains

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Count the number of individuals
    !
    !> @details    This subroutine encapsulates the two different the subroutines
    !>             that count the number of animals and markers
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date       May 11, 2017
    !
    ! PARAMETERS:
    !> @param[in]
    !---------------------------------------------------------------------------
    function CountInData result(nGenotyped)
        use AlphaHmmInMod
        implicit none

        integer :: nGenotyped

        type(AlphaHmmInput), pointer :: inputParams
        inputParams => defaultInput

        nGenotyped = CountInGenotypeData()
    end function CountInData

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Count the number of individuals
    !
    !> @details    This subroutine counts the number of individuals genotyped.
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date       May 11, 2017
    !
    ! PARAMETERS:
    !> @param[in]
    !---------------------------------------------------------------------------
    function CountInGenotypeData result(nGenotyped)
        use hmmPARAMETERS
        use AlphaHmmInMod

        implicit none

        integer :: nGenotyped       !< Number of genotyped/sequenced individuals

        integer :: k
        character (len=300) :: dumC
        type(AlphaHmmInput), pointer :: inputParams

        inputParams => defaultInput
        rewind(inputParams%genotypeFileUnit)

        do
            read (inputParams%genotypeFileUnit,*,iostat=k) dumC
            nGenotyped = nGenotyped + 1
            if (k/=0) then
                nGenotyped = nGenotyped - 1
                exit
            endif
        enddo

        if (inputParams%hmmoption == RUN_HMM_NGS) then
            if(mod(nGenotyped,2)==0) then
                nGenotyped = nGenotyped / 2
            else
                write(0,*) "Error: The number of lines in the file of reads is not even. Is the file corrupt?"
                write(0,*) "The program will now stop"
                stop
            endif
        endif

        print*, " ", nGenotyped, " individuals in the genotype file"

    end function CountInGenotypeData

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Read files
    !
    !> @details    This subroutine reads in data
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date       May 11, 2017
    !
    ! PARAMETERS:
    !> @param[in]
    !---------------------------------------------------------------------------
    subroutine ReadInData(nGenotyped)
        use hmmPARAMETERS
        use AlphaHmmInMod
        implicit none

        integer, intent(in) :: nGenotyped
        type(AlphaHmmInput), pointer :: inputParams

        inputParams=> defaultInput

        ! Read the genotype file
        rewind(inputParams%genotypeFileUnit)
        if (inputParams%hmmoption /= RUN_HMM_NGS) then
            ! if (.not. inputParams%PlinkFormat) then
                call ReadGenotypes(inputParams%genotypeFileUnit, nGenotyped)
            ! else
            !     call ReadPlink(inputParams%genotypeFileUnit)
            ! end if
        else
            call pedigree%addSequenceFromFile(inputparams%GenotypeFile, inputParams%nsnpRaw)
        endif

        close(inputParams%genotypeFileUnit)
    end subroutine ReadInData

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Read Genotypes
    !
    !> @details    This subroutine reads the genotype file
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date       May 11, 2017
    !
    ! PARAMETERS:
    !> @param[in]
    !---------------------------------------------------------------------------
    subroutine ReadGenotypes(GenoFileUnit, nAnisG)
        use GlobalVariablesHmmMaCH
        use AlphaHmmInMod
        implicit none

        integer, intent(in) :: nAnisG
        integer, intent(inout) :: GenoFileUnit
        integer :: i,j
        integer,allocatable,dimension(:) :: temp
        type(AlphaHmmInput), pointer :: inputParams
        logical :: opened, named
        character(len=300) :: GenoFile
        character*(20) :: dumID

        inputParams => defaultInput

        allocate(temp(inputParams%nSnp))
        allocate(GenosHmmMaCH(nAnisG,inputParams%nSnp))

        ! Initialize genotypes with missing
        GenosHmmMaCH=9

        inquire(unit=GenoFileUnit, opened=opened, named=named, name=GenoFile)
        if (.NOT. opened .and. named) then
            open(unit=GenoFileUnit, file=GenoFile, status='unknown')
        else if (.NOT. named) then
            write(0, *) "ERROR - Something went wrong when trying to read the file of pre-phased data"
            open(newunit=GenoFileUnit, file=inputParams%GenotypeFile)
        end if

        do i=1,nAnisG
            read (GenoFileUnit,*) dumID,Temp(:)
            do j=1,inputParams%nsnp
                if ((Temp(j)<0).or.(Temp(j)>2)) Temp(j)=9
            enddo
            GenosHmmMaCH(i,:)=Temp(:)
        enddo

        close(GenoFileUnit)
        deallocate(temp)
    end subroutine ReadGenotypes

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Read Sequence
    !
    !> @details    This subroutine reads the file that contains sequenced data.
    !>             The format of the sequence data is expected to be the number
    !>             of read of the reference allele and the alternative allele
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date       May 11, 2017
    !
    ! PARAMETERS:
    !> @param[in]
    !---------------------------------------------------------------------------
    subroutine ReadSeq(ReadsFileUnit, nAnisS)
        use GlobalVariablesHmmMaCH
        use hmmPARAMETERS
        use AlphaHmmInMod
        implicit none

        integer, intent(in) :: nAnisS
        integer, intent(inout) :: ReadsFileUnit
        integer :: i,j

        type(AlphaHmmInput), pointer :: inputParams
        logical :: opened, named
        character(len=300) :: ReadsFile
        character*(20) :: dumID

        inputParams => defaultInput

        inquire(unit=ReadsFileUnit, opened=opened, named=named, name=ReadsFile)
        if (.NOT. opened .and. named) then
            open(unit=ReadsFileUnit, file=ReadsFile, status='unknown')
        else if (.NOT. named) then
            write(0, *) "ERROR - Something went wrong when trying to read the file of pre-phased data"
            open(newunit=ReadsFileUnit, file=inputParams%GenotypeFile)
        end if

        ! do i=1,nAnisS
        !     read (ReadsFileUnit,*) dumID, Reads%ReferAllele(i,:)
        !     read (ReadsFileUnit,*) dumID, Reads%AlterAllele(i,:)
        ! end do

        ! do i=1,nAnisS
        !     do j=1,inputParams%nsnp
        !         if (Reads%ReferAllele(i,j)>=MAX_READS_COUNT) Reads%ReferAllele(i,j)=MAX_READS_COUNT-1
        !         if (Reads%AlterAllele(i,j)>=MAX_READS_COUNT) Reads%AlterAllele(i,j)=MAX_READS_COUNT-1
        !     enddo
        ! enddo

        close(ReadsFileUnit)

    end subroutine ReadSeq


end module HmmInputMod