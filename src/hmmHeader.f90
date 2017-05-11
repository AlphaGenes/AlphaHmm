!-----------------------------------------------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-----------------------------------------------------------------------------------------------------------------------
!
! MODULE: hmmHeader
!
!> @file        hmmHeader.f90
!
! DESCRIPTION:
!> @brief       Module for the header of the hmm algorithm
!>
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
#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

module hmmHeader

contains

    subroutine Titles

        call PrintVersion
        print *, ""
        print *, ""
        print *, ""

    end subroutine Titles

    subroutine Header

        print *, ""
        print *, "                              ********************                         "
        print *, "                              *                  *                         "
        print *, "                              *     AlphaHMM     *                         "
        print *, "                              *                  *                         "
        print *, "                              ********************                         "
        print *, "                                                                              "
        print *, "                    Software For Phasing and Imputing Genotypes               "

    end subroutine Header

    subroutine PrintVersion

        call Header
        print *, ""
        print *, "                              Commit:   "//TOSTRING(COMMIT),"                     "
        print *, "                              Compiled: "//__DATE__//", "//__TIME__
        print *, ""

    end subroutine PrintVersion

end module hmmHeader