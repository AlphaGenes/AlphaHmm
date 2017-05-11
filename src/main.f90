!-----------------------------------------------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-----------------------------------------------------------------------------------------------------------------------
!
! PROGRAM:      AlphaHMM
!
!> @file        main.f90
!
! DESCRIPTION:
!> @brief       Main program of AlphaHMM
!>
!> @details     This file contains the main program of AlphaHMM
!
!> @author      Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
!
!> @date        May 10, 2017
!
!> @version     0.0.1 (alpha)
!
! REVISION HISTORY:
! 2017.05.10  RAntolin - Initial Version
!
!-----------------------------------------------------------------------------------------------------------------------
Program AlphaHMM
    use hmmModule
    use AlphaHmmInMod
    use InputOutputModule
    use hmmHeader
    implicit none

    inputParams => defaultInput
    if (Command_Argument_Count() > 0) then
        call get_command_argument(1,cmd)
        if (cmd(1:2) .eq. "-v") then
            call PrintVersion
            call exit(0)
        end if
    end if

    if (Command_Argument_Count() > 0) then
        call Get_Command_Argument(1,SpecFile)
    else
        specfile="AlphaImputeSpec.txt"
    end if


    call Titles

    call ReadInParameterFile
    ! call ReadInput

    ! call HMMControler(Parameters)

    ! call WriteOutput

end Program AlphaHMM
