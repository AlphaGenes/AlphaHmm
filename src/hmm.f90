!-----------------------------------------------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-----------------------------------------------------------------------------------------------------------------------
!
! MODULE: hmmModule
!
!> @file        hmm.f90
!
! DESCRIPTION:
!> @brief       Module the hmm algorithm
!>
!> @details     This MODULE contains all the procedures that implements the HMM described in Li et al. 2010, Appendix
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
MODULE hmmModule

    IMPLICIT NONE

CONTAINS


    !######################################################################
    subroutine HMMController(HMM)
        use ConstantModule
        use GlobalVariablesHmmMaCH
        use hmmParameters

        use Par_Zig_mod
        use AlphaHmmInMod
        use omp_lib
        use random


        implicit none
        integer(kind=1), intent(in) :: HMM
        integer :: i, nprocs, nthreads
        real(4) :: r
        real(8) :: tT,t1,t2
        double precision :: Theta
        integer, allocatable :: seed(:)
        integer :: grainsize, count, secs, seed0
        type(AlphaHmmInput), pointer :: inputParams
        integer, allocatable, dimension(:,:) :: InbredHmmMaCH

        inputParams => defaultInput




        ! Read the phased individuals if HMM Only or Sequence data
#ifdef DEBUG
        write(0,*) 'DEBUG: [MaCHController]'
#endif
        nAnisInbred = 0
        if ( (HMM==RUN_HMM_ONLY .OR. HMM==RUN_HMM_NGS) .AND. inputParams%InbredAnimalsFile/="None") then
            call ReadInbred(inputParams%prePhasedFileUnit, InbredHmmMaCH, nAnisInbred)
        end if

        ! Number of animals in the HMM
        ! TODO nAnisInbred broked
        nIndHmmMaCH = pedigree%nGenotyped + nAnisInbred

        ! ALLOCATE MEMORY
#ifdef DEBUG
        write(0,*) 'DEBUG: [MaCHController] Allocate memory'
#endif

        ! Allocate a matrix to store the diploids of every Animal
        ! Template Diploids Library
        ! NOTE: GenosHmmMaCH can contain either genotype or reads information (if working with sequence data NGS)
        allocate(GenosHmmMaCH(nIndHmmMaCH,inputParams%nsnp))
        allocate(PhaseHmmMaCH(nIndHmmMaCH,inputParams%nsnp,2))

        GenosHmmMaCH = MISSING
        PhaseHmmMaCH = ALLELE_MISSING

        ! Allocate memory to store Animals contributing to the Template
        ! Haplotype Library
        allocate(GlobalHmmMachID(nIndHmmMaCH))

        ! Allocate memory to store Animals Highly Dense Genotyped
        allocate(GlobalHmmHDInd(nIndHmmMaCH))

        ! Allocate memory to store Animals Highly Dense Genotyped
        allocate(GlobalHmmPhasedInd(nIndHmmMaCH,2))
        ! No animal has been HD genotyped YET
        ! WARNING: If this variable only stores 1 and 0, then its type should
        !          logical: GlobalHmmHDInd=.false.
        GlobalHmmHDInd=0

        ! Allocate a matrix to store probabilities of genotypes and
        ! alleles for each animal
        allocate(ProbImputeGenosHmm(nIndHmmMaCH,inputParams%nsnp))
        allocate(ProbImputePhaseHmm(nIndHmmMaCH,inputParams%nsnp,2))
        ! Initialise probabilities to 0
        ProbImputeGenosHmm=0.0
        ProbImputePhaseHmm=0.0

        ! The full Template Haplotype Library, H (Li et al. 2010, Appendix)
        allocate(FullH(nIndHmmMaCH,inputParams%nsnp,2))

        ! Vector of Combination of genotyping error (Li et al. 2010, Appendix)
        ! Epsilon is related with the Penetrance Matrix of the HMM which gives
        ! the emision probabilities for each state/marker/snp.
        allocate(Epsilon(inputParams%nsnp))
        allocate(Penetrance(inputParams%nsnp,0:2,0:2))

        ! Vector of Combination of population recombination (Li et al. 2010, Appendix)
        ! Thetas is related with the transition Matrix of the HMM. Since there
        ! are inputParams%nsnp states, there are inputParams%nsnp transitions between states.
        ! WARNING: Is this correctly implemented throughout the hmm process??
        allocate(Thetas(inputParams%nsnp-1))

        allocate(ErrorUncertainty(inputParams%nsnp))
        allocate(ErrorMatches(inputParams%nsnp))
        allocate(Errormismatches(inputParams%nsnp))

        ! Crossover parameter in order to maximize the investigation of
        ! different mosaic configurations
        ! WARNING: crossovers are related with the transition matrix of the HMM.
        !          If there are inputParams%nsnp states and inputParams%nsnp-1 transitions
        !          between states, why the number of crossovers is inputParams%nsnp??
        allocate(Crossovers(inputParams%nsnp-1))

        if (HMM==RUN_HMM_NGS) then
            allocate(ShotgunErrorMatrix(0:2,0:MAX_READS_COUNT,0:MAX_READS_COUNT))
        endif

        ! Set up number of process to be used
        nprocs = OMP_get_num_procs()
        ! call OMP_set_num_threads(inputParams%useprocs)
        nthreads = OMP_get_num_threads()

        allocate(seed(inputParams%useprocs))

        ! Warm up the random seed generator
        call system_clock(count)
        secs = mod(count,int(1e4))
        do i = 1,secs
            call random_number(r)
        enddo

        ! Set up random process seeds for threads
        ! Feed seed as a function of the milliseconds of the system clock
        call system_clock(count)
        secs = mod(count,int(1e6))
        seed0 = secs*1e5
        inputParams%seed=-abs(seed0)
        do i = 1,inputParams%useprocs
            call random_number(r)
            seed(i) = seed0*r
        enddo

        grainsize = 32
        call par_zigset(inputParams%useprocs, seed, grainsize)

#ifdef DEBUG
        write(0,*) 'DEBUG: [ParseMaCHData] ...'
#endif

        ! Populate genotype and phased data from input
        call ParseMaCHData(HMM, InbredHmmMaCH, pedigree%nGenotyped, nAnisInbred)
        if (nAnisInbred > 0) then
            deallocate(InbredHmmMaCH)
        end if

        ! Initialization of HMM parameters
        Epsilon=EPSILON_ERROR
        Thetas=0.01

        do i=1,inputParams%nsnp
            call SetPenetrance(i, EPSILON_ERROR)
        enddo

        if (HMM==RUN_HMM_NGS) then
            call SetShotgunError(SEQUENCING_ERROR)
        endif

#ifdef DEBUG
        write(0,*) 'DEBUG: [SetUpEquations] ...'
#endif

        ! Set up Reference haplotypes and HMM parameters
        ! TODO: DEBUG THIS CODE
        call SetUpEquations(HMM, pedigree%nGenotyped, nAnisInbred)

        open (unit=6,form='formatted')

        print*, ""
        print*, " Impute genotypes by HMM"
        print*, "    Using", nthreads, "processors of", nprocs

        ! Allocate and set up variables storing allele frequencies
        allocate(GenosCounts(nIndHmmMaCH,inputParams%nsnp,2))
        GenosCounts=0

        tT=0.0

        do GlobalRoundHmm=1,inputParams%nroundshmm
            write(6, 100, advance="no") char(13),"   HMM Round   ",GlobalRoundHmm
            !write(6, 100) char(13),"   HMM Round   ",GlobalRoundHmm
            100 format (a1, a17, i10)
            call flush(6)

            call ResetCrossovers
            ! Parallelise the HMM process in animals
#if DEBUG.EQ.1
            t1 = omp_get_wtime()
            write(0,*) 'DEBUG: Begin paralellisation [MaCHController]'
#endif
            !$OMP PARALLEL DO DEFAULT(shared) &
            !$OMP PRIVATE(i) &
            !$OMP SCHEDULE(static)
            do i=1,nIndHmmMaCH
                ! print *,"i:",i,nIndHmmMaCH
                call MaCHForInd(i, HMM)
                ! print *,"after"
            enddo
            !$OMP END PARALLEL DO

#if DEBUG.EQ.1
            t2 = omp_get_wtime()
            tT = tT + (t2-t1)
            write(0,*) 'DEBUG: End paralellisation. Partial time:', t2-t1
#endif

            Theta = 0.01

            ! Update emission probabilities of the HMM process
            call UpdateThetas

            ! Update transition probabilities of the HMM process
            call UpdateErrorRate(Theta)
        enddo

#ifdef DEBUG
        write(0,*) 'DEBUG: End paralellisation'
#endif

        ! Average genotype probability of the different hmm processes
        ProbImputeGenosHmm=ProbImputeGenosHmm/(inputParams%nroundshmm-inputParams%hmmburninround)
        ProbImputePhaseHmm=ProbImputePhaseHmm/(inputParams%nroundshmm-inputParams%hmmburninround)

    end subroutine HMMController

    !######################################################################
    subroutine ParseMaCHData(HMM, PhasedData, nGenotyped, nInbred)
        use GlobalVariablesHmmMaCH
        use hmmParameters
        use AlphaHmmInMod

        implicit none

        integer(kind=1),intent(in) :: HMM
        integer, intent(in) :: nGenotyped, nInbred
        integer, intent(in) :: PhasedData(:,:)

        integer :: maxHaps
        type(AlphaHmmInput), pointer :: inputParams

        inputParams => defaultInput

        allocate(GlobalInbredInd(nGenotyped+nInbred))
        GlobalInbredInd=.FALSE.

        if (HMM == RUN_HMM_NGS) then
            call ParseMaCHDataNGS(nGenotyped)
        else
            call ParseMaCHDataGenos(nGenotyped)
        endif
        if (nInbred > 0) then
            call ParseMaCHPhased(PhasedData, nGenotyped, nInbred)
        endif

        ! Check if the number of Haplotypes the user has considered in the
        ! Spec file, Sub H (MaCH paper: Li et al. 2010), is reached.
        maxHaps = 2 * (sum(GlobalHmmHDInd(1:nGenotyped))-1) + nInbred
        if (inputParams%nhapinsubh > maxHaps) then
            print*, ""
            print*, "WARNING! Number of individuals highly-covered is too small"
            print*, "         for the number of Haplotypes specified."
            print*, "         Reference haplotypes will be taken from the whole population"
            GlobalHmmHDInd=1
            ! stop
        endif

    end subroutine ParseMaCHData

    !######################################################################
    subroutine ParseMaCHPhased(PhasedData, nGenotyped, nInbred)
        use ISO_Fortran_Env
        use GlobalVariablesHmmMaCH


        implicit none
        integer, intent(in) :: nGenotyped, nInbred
        integer, intent(in) :: PhasedData(:,:)

        integer :: i

        print*, 'DEBUG: [ParseMaCHPhased]'
        do i = 1, nInbred
            GenosHmmMaCH(nGenotyped+i,:) = 2 * PhasedData(i,:)
            PhaseHmmMaCH(nGenotyped+i,:,1) = PhasedData(i,:)
            PhaseHmmMaCH(nGenotyped+i,:,2) = PhasedData(i,:)
        enddo

        GlobalInbredInd(nGenotyped+1 : nGenotyped+nInbred) = .TRUE.
        GlobalHmmHDInd(nGenotyped+1 : nGenotyped+nInbred) = 1


    end subroutine ParseMaCHPhased

    !######################################################################
    subroutine ReadInbred(PhaseFileUnit, PhasedData, nInbred)
        use ISO_Fortran_Env
        use GlobalVariablesHmmMaCH
        use AlphaHmmInMod

        implicit none

        integer, intent(inout) :: PhaseFileUnit
        integer, intent(out), allocatable, dimension(:,:) :: PhasedData
        integer, intent(out) :: nInbred

        integer :: i,k,dumC
        logical :: opened, named
        character(len=300) :: InbredFile
        character(len=20) :: dumID
        type(AlphaHmmInput), pointer :: inputParams
        inputParams => defaultInput

        print *, "DEBUG: [ReadInbred]"
        inquire(unit=PhaseFileUnit, opened=opened, named=named, name=InbredFile)

        if (.NOT. opened .and. named) then
            open(unit=PhaseFileUnit, file=InbredFile, status='unknown')
        else if (.NOT. named) then
            write(0, *) "ERROR - Something went wrong when trying to read the file of pre-phased data"
        end if

        nInbred = 0
        do
            read (PhaseFileUnit,*,iostat=k) dumC
            nInbred=nInbred+1
            if (k/=0) then
                nInbred=nInbred-1
                exit            ! This forces to exit if an error is found
            endif
        enddo
        rewind(PhaseFileUnit)

        allocate(PhasedData(nInbred,inputParams%nsnp))

        do i=1,nInbred
            read (PhaseFileUnit,*) dumID, PhasedData(i,:)
        end do
        close(PhaseFileUnit)

    end subroutine ReadInbred

    !######################################################################
    subroutine getHapList(ListIds, nHaps)
        use GlobalVariablesHmmMaCH
        use AlphaHmmInMod
        implicit none

        ! Dummy Arguments
        character(len=20), allocatable, intent(inout) :: ListIds(:)
        integer, intent(out) :: nHaps

        ! Local Variables
        integer :: i, k
        logical :: opened, exists
        character(len=20) :: dumC
        type(AlphaHmmInput), pointer :: inputParams

        inputParams => defaultInput

        inquire(file=trim(inputParams%HapListFile), opened=opened, exist=exists, number=inputParams%HapListUnit)
        if (exists) then
            if (.not. opened) then
                open(newunit=inputParams%HapListUnit, file=trim(inputParams%HapListFile), status="old")
            end if
        else
            write(0,*) "ERROR: File <", trim(inputParams%HapListFile), "> does not exist"
            stop
        end if

        nHaps = 0
        do
            read(inputParams%HapListUnit, *, iostat=k) dumC
            nHaps=nHaps+1
            if(k/=0) then
                nHaps=nHaps-1
                exit
            endif
        enddo

        rewind(inputParams%HapListUnit)

        allocate(ListIds(nHaps))
        do i=1,nHaps
            read(inputParams%HapListUnit, *) ListIds(i)
        enddo
        close(inputParams%HapListUnit)

    end subroutine getHapList

    !######################################################################
    subroutine ParseMaCHDataNGS(nGenotyped)
        use GlobalVariablesHmmMaCH
        use AlphaHmmInMod

        implicit none
        integer, intent(in) :: nGenotyped
        integer :: i, j, nHaps
        type(AlphaHmmInput), pointer :: inputParams
        character(len=20), allocatable :: HapList(:)

#ifdef DEBUG
        write(0,*) 'DEBUG: [ParseMaCHDataNGS]'
#endif

        inputParams => defaultInput

        if (inputParams%HapList) then
            nHaps=0
            call getHapList(HapList, nHaps)
        end if

        do i=1,nGenotyped
            ! Add animal's diploid to the Diploids Library
            ! Find individuals sequenced with high coverage
            if ( float( count(pedigree%pedigree(pedigree%genotypeMap(i))%ReferAllele(:) + &
                               pedigree%pedigree(pedigree%genotypeMap(i))%AlterAllele(:) /= READ_MISSING) ) / inputParams%nSnp > 0.90 ) then

                GlobalHmmHDInd(i)=1
            endif
            if (inputParams%HapList) then
                do j=1,nHaps
                    if (pedigree%pedigree(pedigree%genotypeMap(i))%originalID == HapList(j)) then
                        GlobalInbredInd(i) = .TRUE.
                        exit
                    endif
                enddo
            endif
        enddo

        ! AlphaImpute does not phase sequence data, thus no individual has been phased.
        nGametesPhased = 0
    end subroutine ParseMaCHDataNGS

    !######################################################################
    subroutine ParseMaCHDataGenos(nGenotyped)
        use GlobalVariablesHmmMaCH

        use Utils
        use AlphaHmmInMod
        use iso_fortran_env

        implicit none

        type(AlphaHmmInput), pointer :: inputParams
        integer, intent(in) :: nGenotyped

        integer :: i,j, NoGenosUnit, nIndvG

        inputParams => defaultInput
! #ifdef DEBUG
        write(0,*) 'DEBUG: [ParseMaCHDataGenos] ...'
! #endif

        NoGenosUnit = 111
        open(unit=NoGenosUnit, file='Miscellaneous/NotGenotypedAnimals.txt', status="replace")

        GlobalHmmPhasedInd=.FALSE.
        i=0
        nIndvG=0
        ! Read both the phased information of AlphaImpute and high-denisty genotypes,
        ! store it in PhaseHmmMaCH and GenosHmmMaCH, and keep track of  which
        ! gamete is phased (GlobalHmmPhasedInd) and the high-denisity
        ! genotyped animal (GlobalHmmHDInd)

        do i=1,nGenotyped
            nIndvG=nIndvG+1

            ! Add animal's diploid to the Diploids Library
            GenosHmmMaCH(i,:)=imputeGenosHMM(i,:)

            ! Take the phased information from AlphaImpute
            PhaseHmmMaCH(i,:,1)=imputePhaseHmm(i,:,1)
            PhaseHmmMaCH(i,:,2)=imputePhaseHmm(i,:,2)

            ! Check if this animal is Highly Dense genotyped
            if ((float(count(GenosHmmMaCH(i,:)==MISSING))/inputParams%nsnp)<0.10) then
                GlobalHmmHDInd(i)=1
            endif

            ! Clean the genotypes and alleles from possible coding errors
            do j=1,inputParams%nsnp
                if ((GenosHmmMaCH(i,j)<0).or.(GenosHmmMaCH(i,j)>2)) GenosHmmMaCH(i,j)=MISSING
                if ((PhaseHmmMaCH(i,j,1)/=0) .AND. (PhaseHmmMaCH(i,j,1)/=1)) PhaseHmmMaCH(i,j,1)=ALLELE_MISSING
                if ((PhaseHmmMaCH(i,j,2)/=0) .AND. (PhaseHmmMaCH(i,j,2)/=1)) PhaseHmmMaCH(i,j,2)=ALLELE_MISSING
            enddo

            ! Check if this individual has its haplotypes phased
            if (float(count(PhaseHmmMaCH(i,:,1)/=ALLELE_MISSING))/inputParams%nsnp >= (imputedThreshold/100.0)) Then
                GlobalHmmPhasedInd(i,1)=.TRUE.
            endif
            if (float(count(PhaseHmmMaCH(i,:,2)/=ALLELE_MISSING))/inputParams%nsnp >= (imputedThreshold/100.0)) Then
                GlobalHmmPhasedInd(i,2)=.TRUE.
            endif

            ! Count the number of phased animals
            if ((GlobalHmmPhasedInd(i,1)==.TRUE.).AND.(GlobalHmmPhasedInd(i,2)==.TRUE.)) Then
                nAnimPhased=nAnimPhased+1
            endif

        end do

        close(NoGenosUnit)

        ! Count the number of phased gametes
        nGametesPhased = 0
        nGametesPhased = CountPhasedGametes()

        ! Check if the number of genotyped animals is correct
        if (nIndvG/=nGenotyped) then
            write (6,*) '   ','WARNING: There are individuals in the genotype file that have'
            write (6,*) '   ','         not been genotyped'
            write (6,*) '   ','         For a list of these individuals look into the file'
            write (6,*) '   ','         Miscellaneous/NotGenotypedAnimals.txt'
            ! stop
        endif

    end subroutine ParseMaCHDataGenos

    !######################################################################
    subroutine MaCHForInd(CurrentInd, HMM)
        ! Create a Template Haplotype Library, H, and create HMM for each
        ! individual

        use GlobalVariablesHmmMaCH
        use hmmParameters
        use hmmHaplotyper
        use Utils
        use random
        use Par_Zig_mod
        use omp_lib
        use AlphaHmmInMod

        implicit none

        type(AlphaHmmInput), pointer :: inputParams
        integer, intent(in) :: CurrentInd
        integer(kind=1), intent(in) :: HMM

        ! Local variables
        integer :: genotypeInt, i, states
        integer :: StartSnp, StopSnp

        inputParams => defaultInput

        ! The number of parameters of the HMM are:
        !   inputParams%nhapinsubh = Number of haplotypes in the template haplotype set, H
        !   inputParams%nhapinsubh*(inputParams%nhapinsubh+1)/2 = Number of states, N = H^2 (Li et al, 2010)
        ! ForwardProbs are the accumulated probabilities, and
        ! ForwardPrbos(:,1) are the prior probabilities
        ! Allocate all possible state sequencies
        states = inputParams%nhapinsubh*(inputParams%nhapinsubh+1)/2
        ! allocate(ForwardProbs(states,inputParams%nsnp))
        allocate(SubH(inputParams%nhapinsubh,inputParams%nsnp))

        call ExtractTemplate(HMM, currentInd, nIndHmmMaCH)

        StartSnp=1
        StopSnp=inputParams%nsnp
        if (HMM==RUN_HMM_ONLY .OR. HMM==RUN_HMM_NGS) then
            if (GlobalInbredInd(CurrentInd)==.TRUE.) then
                allocate(ForwardProbs(inputParams%nhapinsubh,inputParams%nsnp))
                call ForwardAlgorithmForSegmentHaplotype(CurrentInd,1,1,inputParams%nsnp)     ! Paternal haplotype
                call SampleSegmentHaplotypeSource(CurrentInd,1,1,inputParams%nsnp)

                deallocate(ForwardProbs)
                allocate(ForwardProbs(inputParams%nhapinsubh,inputParams%nsnp))
                call ForwardAlgorithmForSegmentHaplotype(CurrentInd,2,1,inputParams%nsnp)     ! Maternal haplotype
                call SampleSegmentHaplotypeSource(CurrentInd,2,1,inputParams%nsnp)
            else
                allocate(ForwardProbs(states,inputParams%nsnp))
                call ForwardAlgorithm(CurrentInd)
                call SampleChromosomes(CurrentInd,StartSnp,StopSnp)
            end if
        else
            if (nGametesPhased/float(2*pedigree%nGenotyped)>inputParams%phasedThreshold/100.0) then
                if (GlobalHmmPhasedInd(CurrentInd,1)/=.TRUE. .AND. GlobalHmmPhasedInd(CurrentInd,2)/=.TRUE.) Then
                    allocate(ForwardProbs(states,inputParams%nsnp))
                    call ForwardAlgorithm(CurrentInd)
                    call SampleChromosomes(CurrentInd,1,inputParams%nsnp)
                else
                    allocate(ForwardProbs(inputParams%nhapinsubh,inputParams%nsnp))
                    call ForwardAlgorithmForSegmentHaplotype(currentInd,1,1,inputParams%nsnp)     ! Paternal haplotype
                    call SampleSegmentHaplotypeSource(CurrentInd,1,1,inputParams%nsnp)
                    deallocate(ForwardProbs)
                    allocate(ForwardProbs(inputParams%nhapinsubh,inputParams%nsnp))
                    call ForwardAlgorithmForSegmentHaplotype(currentInd,2,1,inputParams%nsnp)     ! Paternal haplotype
                    call SampleSegmentHaplotypeSource(CurrentInd,2,1,inputParams%nsnp)
                endif
            else
                allocate(ForwardProbs(states,inputParams%nsnp))
                call ForwardAlgorithm(CurrentInd)
                call SampleChromosomes(CurrentInd,1,inputParams%nsnp)
            endif
        endif

        ! WARNING: The idea of not to use the first inputParams%hmmburninround rounds suggests
        !          the imputation in those rounds aren't accurate enough, which
        !          also suggests that each round a better solution is found.
        !          Better solutions are obtained by improving previous solutions
        !          by means of update HMM parameters (recombinations  rates,
        !          Thetas, and the genotyping errors) as implemented in MaCH
        !          code with functions UpdateThetas, UpdateErrorRate and
        !          TotalCrossovers.
        !
        !          However, each time the MaCHForInd subroutine is called is
        !          independent from the previous call and so, HMM solutions
        !          given by ForwardAlgorithm and SampleChromosomes are also
        !          independent.

#if DEBUG.EQ.1
        write(0,*) 'DEBUG: Calculate genotype counts [MaCHForInd]'
#endif
        if (GlobalRoundHmm>inputParams%hmmburninround) then
            do i=1,inputParams%nsnp
                genotypeInt = FullH(CurrentInd,i,1)+FullH(CurrentInd,i,2)
                if (genotypeInt==2) then
                    GenosCounts(CurrentInd,i,2)=GenosCounts(CurrentInd,i,2)+1
                elseif (genotypeInt==1) then
                    GenosCounts(CurrentInd,i,1)=GenosCounts(CurrentInd,i,1)+1
                endif
            enddo
            ! endif

            ! Cumulative genotype probabilities through hmm processes
            !$omp workshare
            ProbImputeGenosHmm(CurrentInd,:)=ProbImputeGenosHmm(CurrentInd,:)&
                +FullH(CurrentInd,:,1)+FullH(CurrentInd,:,2)
            ProbImputePhaseHmm(CurrentInd,:,1)=ProbImputePhaseHmm(CurrentInd,:,1)&
                +FullH(CurrentInd,:,1)
            ProbImputePhaseHmm(CurrentInd,:,2)=ProbImputePhaseHmm(CurrentInd,:,2)&
                +FullH(CurrentInd,:,2)
            !$omp end workshare
        endif

#if DEBUG.EQ.1
        write(0,*) 'DEBUG: Deallocate Forward variable and Haplotype Template [MaCHForInd]'
#endif
        deallocate(ForwardProbs)
        deallocate(SubH)

    end subroutine MaCHForInd

    !######################################################################
    subroutine SampleChromosomes(CurrentInd,StartSnp,StopSnp)
        use GlobalVariablesHmmMaCH

        use omp_lib
        use Par_Zig_mod
        use AlphaHmmInMod
        use hmmParameters
        implicit none

        integer,intent(in) :: CurrentInd,StartSnp,StopSnp
        type(AlphaHmmInput), pointer :: inputParams

        ! Local variables
        integer :: i,j,k,l,SuperJ,Index,OffOn,State1,State2,TmpJ,TopBot,FirstState,SecondState,Tmp,Thread
        double precision :: Summer,Choice,Sum00,Sum01,Sum10,Sum11
        double precision,dimension(:), allocatable :: Probs
        double precision :: Theta


        inputParams => defaultInput
        allocate(probs(inputParams%nhapinsubh*(inputParams%nhapinsubh+1)/2))
        Summer=0.0
        Index=0
        Probs = ForwardProbs(:,StopSnp)
        Thread = omp_get_thread_num()

        ! Calculate sum over all states. The sum corresponds to all the
        ! forward probabilities, that is, the probability of the sequence of
        ! observed genotypes of animal CurrentInd.
        do i=1,inputParams%nhapinsubh
            do j=1,i
                Index=Index+1
                Summer=Summer+Probs(Index)
            enddo
        enddo

        ! Sample number at random and select state: (State1,State2)
        Choice = par_uni(Thread)*Summer
        Summer=0.0
        Index=0
        OffOn=0

        State1 = 1
        State2 = 1

        Probs = ForwardProbs(:,StopSnp)
        do i=1,inputParams%nhapinsubh
            do j=1,i
                Index=Index+1
                Summer=Summer+Probs(Index)
                if (Summer>Choice) then
                    State1=i
                    State2=j
                    OffOn=1
                    exit
                endif
            enddo
            if (OffOn==1) exit
        enddo

        SuperJ=StopSnp
        do while (SuperJ>StartSnp)
            SuperJ=SuperJ-1
            if (inputParams%HMMOption/=RUN_HMM_NGS) then
                call ImputeAlleles(CurrentInd,SuperJ+1,State1,State2)
            else
                call ImputeAllelesNGS(CurrentInd,SuperJ+1,State1,State2)
            endif
            TmpJ=SuperJ

            ! Cumulative recombination fraction allows us to skip over
            ! uninformative positions: Alleles with missing genotype are skipped
            ! but the recombination information (Thetas(SuperJ) is accumulated
            ! and used in the next location.
            Theta=Thetas(SuperJ)
            do while ((GenosHmmMaCH(CurrentInd,SuperJ)==MISSING).and.SuperJ>1)
                SuperJ=SuperJ-1
                Theta=Theta+Thetas(SuperJ)-Theta*Thetas(SuperJ)
            enddo

            ! When examining the previous location we consider three alternatives:
            !   * states that could be reached when both haplotypes recombine (11),
            !   * states that can be reached when the first (10) or second (01)
            !       haplotype recombines, and
            !   * the states that can be reached without recombination (00).
            Sum00=0.0
            Sum01=0.0
            Sum10=0.0
            Sum11=0.0

            Index=0
            Probs = ForwardProbs(:,SuperJ)
            do k=1,inputParams%nhapinsubh
                do l=1,k
                    Index=Index+1
                    Sum11=Sum11+Probs(Index)
                    if ((State1==k).or.(State1==l))&
                        Sum01=Sum01+Probs(Index)
                    if ((State2==k).or.(State2==l))&
                        Sum10=Sum10+Probs(Index)
                    if (((State1==k).and.(State2==l))&
                        .or.((State1==l).and.(State2==k)))&
                        Sum00=Sum00+Probs(Index)
                enddo
            enddo

            Summer=Sum11*Theta*Theta/(inputParams%nhapinsubh*inputParams%nhapinsubh)&
                +(Sum10+Sum01)*Theta*(1.0-Theta)/inputParams%nhapinsubh&
                +Sum00*(1.0-Theta)*(1.0-Theta)

            ! WARNING: Why is this assignment here?!?! In case it has to exit,
            !          shouldn't it exit before?
            if (SuperJ==1) exit


            !Sample number and decide how many state changes occurred between the
            !two positions
            Choice = par_uni(Thread)*Summer

            !The most likely outcome is that no changes occur ...
            Choice=Choice-(Sum00*(1.0-Theta)*(1.0-Theta))

            if (Choice<=0.0) then
                !Record outcomes for intermediate, uninformative, positions
                TopBot=1
                call FillPath(CurrentInd,SuperJ,TmpJ+1,State1,TopBot)
                TopBot=2
                call FillPath(CurrentInd,SuperJ,TmpJ+1,State2,TopBot)
                cycle
            endif

            !But perhaps the first or second haplotype recombined
            Choice=Choice-(Sum10*Theta*(1.0-Theta)/inputParams%nhapinsubh)

            ! Did the first hap recombine?
            if (Choice<=0.0) then
                !The first haplotype changed ...
                Choice=Choice*inputParams%nhapinsubh/(Theta*(1.0-Theta))
                !Record the original state
                FirstState=State1

                ! Sample number at random and decide haplotype
                do State1=1,inputParams%nhapinsubh
                    if (State1>=State2) then
                        Choice=Choice+ForwardProbs(State1*(State1-1)/2+State2,SuperJ)
                    else
                        Choice=Choice+ForwardProbs(State2*(State2-1)/2+State1,SuperJ)
                    endif
                    if (Choice>=0.0) exit
                enddo

                !Record outcomes for intermediate, uninformative, positions
                TopBot=1
                call SamplePath(CurrentInd,SuperJ,TmpJ+1,State1,FirstState,TopBot)
                TopBot=2
                call FillPath(CurrentInd,SuperJ,TmpJ+1,State2,TopBot)
                cycle
            endif

            Choice=Choice-(Sum01*Theta*(1.0-Theta)/inputParams%nhapinsubh)

            ! Did the second hap recombine?
            if (Choice<=0.0) then
                !The second haplotype changed ...
                Choice=Choice*inputParams%nhapinsubh/(Theta*(1.0-Theta))
                !Save the original state
                SecondState=State2

                ! Sample number at random and decide haplotype
                ! WARNING: State2 variable should be set to 0 before the loop!?!?
                !        do while (State2<inputParams%nhapinsubh)
                !            State2=State2+1
                do State2=1,inputParams%nhapinsubh
                    if (State1>=State2) then
                        Choice=Choice+ForwardProbs(State1*(State1-1)/2+State2,SuperJ)
                    else
                        Choice=Choice+ForwardProbs(State2*(State2-1)/2+State1,SuperJ)
                    endif
                    if (Choice>=0.0) exit
                enddo

                !Record outcomes for intermediate, uninformative, positions
                TopBot=1
                call FillPath(CurrentInd,SuperJ,TmpJ+1,State1,TopBot)
                TopBot=2
                call SamplePath(CurrentInd,SuperJ,TmpJ+1,State2,SecondState,TopBot)
                cycle
            endif

            !Try to select any other state
            Choice=Choice*inputParams%nhapinsubh*inputParams%nhapinsubh/(Theta*Theta)

            !Save the original states
            FirstState=State1
            SecondState=State2

            Summer=0.0
            Index=0
            OffOn=0
            do i=1,inputParams%nhapinsubh
                do j=1,i
                    Index=Index+1
                    Summer=Summer+ForwardProbs(Index,SuperJ)
                    if (Summer>Choice) then
                        State1=i
                        State2=j
                        OffOn=1
                        exit
                    endif
                enddo
                if (OffOn==1) exit
            enddo

            ! Shuffle haplotypes at random
            if (par_uni(Thread)>0.5) then
                Tmp=State1
                State2=State1
                State2=Tmp
            endif

            !Record outcomes for intermediate, uninformative, positions
            TopBot=1
            call SamplePath(CurrentInd,SuperJ,TmpJ+1,State1,FirstState,TopBot)
            TopBot=2
            call SamplePath(CurrentInd,SuperJ,TmpJ+1,State2,SecondState,TopBot)

        enddo

        if (inputParams%HMMOption/=RUN_HMM_NGS) then
            call ImputeAlleles(CurrentInd,StartSnp,State1,State2)
        else
            call ImputeAllelesNGS(CurrentInd,StartSnp,State1,State2)
        endif

        deallocate(probs)
    end subroutine SampleChromosomes

    !######################################################################
    subroutine FillPath(CurrentInd,FromMarker,ToMarker,State,TopBot)
        ! Impute alleles to a haplotype region. The region is from FromMarker
        ! to ToMarker, none of them included.


        implicit none

        integer,intent(in) :: CurrentInd,FromMarker,ToMarker,State,TopBot

        ! Local variable
        integer :: j

        do j=FromMarker+1,ToMarker-1
            call ImputeAllele(CurrentInd,j,State,TopBot)
        enddo

    end subroutine FillPath

    !######################################################################
    subroutine SamplePath(CurrentInd,FromMarker,ToMarker,FromState,ToState,TopBot)
        ! Impute a path between the two end markers, assuming no genotypes
        ! are observed -- the only constraint is that we must start at
        ! fromState and end at toState with at least one intervening recombinant


        use omp_lib
        use Par_Zig_mod
        use AlphaHmmInMod
        use GlobalVariablesHmmMaCH
        implicit none

        integer,intent(in) :: CurrentInd,TopBot,FromMarker,ToMarker,FromState,ToState
        type(AlphaHmmInput), pointer :: inputParams
        ! Local variables
        integer :: i,State,FromMarkerLocal, Thread
        double precision :: Recomb,Theta1
        double precision :: Theta


        inputParams => defaultInput

        Thread = omp_get_thread_num()

        FromMarkerLocal=FromMarker
        Theta=0.0
        ! Calculate overall recombination fraction for the interval,
        ! FromMarker excluded
        do i=FromMarkerLocal,ToMarker-1
            Theta=Thetas(i)+Theta-Theta*Thetas(i)
        enddo


        ! Impute a path between the two end markers
        do while (FromMarkerLocal<ToMarker-1)
            ! Random recombination
            Recomb = par_uni(Thread)*Theta

            ! Recombination fraction of the FromMarkerLocal
            Theta1=Thetas(FromMarkerLocal)

            ! Calculate overall recombination fraction for the interval,
            ! FromMarkerLocal included
            if (Theta < 0.9) then
                !Fast closed formula
                Theta=(Theta-Theta1)/(1.0-Theta1)
            else
                Theta = 0.0
                !More accurate, iterative formula
                do i=FromMarkerLocal+1,ToMarker-1
                    Theta=Thetas(i)+Theta-Theta*Thetas(i)
                enddo
            endif

            if (Recomb>Theta1) then
                ! No recombinant in the first interval =>
                !    => Recombinant in second interval
                FromMarkerLocal=FromMarkerLocal+1
                call ImputeAllele(CurrentInd,FromMarkerLocal,FromState,TopBot)
                cycle
            endif

            ! If there is no recombinant in the second interval, then
            ! there is recombinant in the first...
            !$OMP ATOMIC
            Crossovers(FromMarkerLocal)=Crossovers(FromMarkerLocal)+1

            if (Recomb<Theta1*(1.0-Theta)) then
                ! No recombinant in the second interval
                call FillPath(CurrentInd,FromMarkerLocal,ToMarker,ToState,TopBot);
                return
            else
                ! Recombinants in both intervals, so we must sample
                ! an intervening state
                FromMarkerLocal=FromMarkerLocal+1
                State=int(par_uni(Thread)*inputParams%nhapinsubh)+1

                call ImputeAllele(CurrentInd,FromMarkerLocal,State,TopBot)
            endif

        enddo

        !If we get here, record obligate recombinant between two consecutive markers
        !$OMP ATOMIC
        Crossovers(FromMarkerLocal)=Crossovers(FromMarkerLocal)+1

    end subroutine SamplePath

    !######################################################################
    subroutine ImputeAlleles(CurrentInd,CurrentMarker,State1,State2)
        ! Impute alleles to haplotypes based on the HMM information. Count the
        ! number of uncertainties, matches and mismatches of the imputed
        ! alleles according to the genotype information that the individual
        ! carries.

        use omp_lib
        use Par_Zig_mod
        use GlobalVariablesHmmMaCH
        implicit none

        integer,intent(in) :: CurrentInd,CurrentMarker,State1,State2

        ! Local variables
        integer :: Imputed1,Imputed2,genotypeInt,Differences, Thread

        Thread = omp_get_thread_num()
        ! These will be the observed imputed alleles defined by the state:
        ! Sj=(State1, State2)
        Imputed1=SubH(State1,CurrentMarker)
        Imputed2=SubH(State2,CurrentMarker)

        ! This is the individual observed genotype
        genotypeInt=GenosHmmMaCH(CurrentInd,CurrentMarker)

        if ((genotypeInt/=0).and.(genotypeInt/=2)) then
            FullH(CurrentInd,CurrentMarker,1)=Imputed1
            FullH(CurrentInd,CurrentMarker,2)=Imputed2
        endif

        ! If genotype is missing, skip
        if (genotypeInt==3) return

        ! Difference between the observed genotype and the gentoype implied by
        ! the state S=(State1, State2)
        Differences=abs(genotypeInt - (Imputed1+Imputed2))

        ! If allele is heterozygous, there is uncertainty
        if ((genotypeInt==1).and.(Differences==0)) then
            !$OMP ATOMIC
            ErrorUncertainty(CurrentMarker)=ErrorUncertainty(CurrentMarker)+1

        ! If allele is homozygous or the genotype does not agree the observation

        else
            ! count the number of alleles matching
            !$OMP ATOMIC
            ErrorMatches(CurrentMarker)=ErrorMatches(CurrentMarker)+(2-Differences)
            ! count the number of mismatching alleles
            !$OMP ATOMIC
            Errormismatches(CurrentMarker)=Errormismatches(CurrentMarker)+Differences
        endif

        ! If gentoype is homozygous or missing, the skip
        if (genotypeInt/=1) return

        ! If the observed allele is homozygous but the genotype is heterozygous
        if (Imputed1==Imputed2) then
            ! Impute the allele to the paternal or maternal haplotype at random
            if (par_uni(Thread)>=0.5) then
                FullH(CurrentInd,CurrentMarker,1)=abs(Imputed1-1)
            else
                FullH(CurrentInd,CurrentMarker,2)=abs(Imputed2-1)
            endif
        endif

    end subroutine ImputeAlleles

    !######################################################################
    subroutine ImputeAllele(CurrentInd,CurrentMarker,State,TopBot)
        use GlobalVariablesHmmMaCH
        implicit none

        integer :: CurrentInd,CurrentMarker,State,TopBot

        FullH(CurrentInd,CurrentMarker,TopBot)=SubH(State,CurrentMarker)

    end subroutine ImputeAllele

    !######################################################################
    subroutine ForwardAlgorithm(CurrentInd)
        ! Update the forward variable of the HMM model


        use omp_lib
        use GlobalVariablesHmmMaCH
        use AlphaHmmInMod
        implicit none

        integer, intent(in) :: CurrentInd
        double precision :: Theta

        ! Local variables
        integer :: j,PrecedingMarker
        type(AlphaHmmInput), pointer :: inputParams
        inputParams => defaultInput

        ! Setup the initial state distributions

        call SetUpPrior

        j=1
        ! CurrentInd is the individual being studied and it is necessary to
        ! obtain the genotype of this individual in ConditionaOnData subroutine
        ! For j=1, ConditionaOnData will initialize the variable
        ! ForwardProbs(:,1) with the Prior Probabilities
#if DEBUG.EQ.1
        write(0,*) 'DEBUG: Update forward variable with emission probabilities [ForwardAlgorithm]'
#endif
        call ConditionOnData(CurrentInd,j)

        ! WARNING: This variable, Theta, should be considered as local as is
        !          global through out the HMM code for different purposes.
        !          Look at subroutines Transpose, SampleChromosomes and
        !          SamplePath
        Theta=0.0

        PrecedingMarker=1
#if DEBUG.EQ.1
        write(0,*) 'DEBUG: Calculate Forward variables [ForwardAlgorithm]'
#endif

        do j=2,inputParams%nsnp
            ! Cumulative recombination fraction allows us to skip uninformative positions
            Theta=Theta+Thetas(j-1)-Theta*Thetas(j-1)
            ! Skip over uninformative positions to save time
            if ((GenosHmmMaCH(CurrentInd,j)/=MISSING).or.(j==inputParams%nsnp)) then
                call Transpose(j,PrecedingMarker,Theta)
                call ConditionOnData(CurrentInd,j)
                PrecedingMarker=j
                Theta=0.0
            endif
        enddo

    end subroutine ForwardAlgorithm

    !######################################################################
    subroutine ForwardAlgorithmForSegment(CurrentInd, StartSnp, StopSnp)
        ! Update the forward variable of the HMM model


        use omp_lib
        use GlobalVariablesHmmMaCH
        use AlphaHmmInMod
        implicit none

        integer, intent(in) :: CurrentInd, StartSnp, StopSnp
        double precision :: Theta

        ! Local variables
        integer :: j,PrecedingMarker
        type(AlphaHmmInput), pointer :: inputParams
        inputParams => defaultInput

        ! Setup the initial state distributions
        call SetUpPrior


#if DEBUG.EQ.1
        write(0,*) 'DEBUG: Update forward variable with emission probabilities [ForwardAlgorithm]'
#endif
        ! CurrentInd is the individual being studied and it is necessary to
        ! obtain the genotype of this individual in ConditionaOnData subroutine
        ! For j=StartSnp, ConditionaOnData will initialize the variable
        ! ForwardProbs(:,1) with the Prior Probabilities
        j=StartSnp
        call ConditionOnData(CurrentInd,j)

        ! WARNING: This variable, Theta, should be considered as local as is
        !          global through out the HMM code for different purposes.
        !          Look at subroutines Transpose, SampleChromosomes and
        !          SamplePath
        Theta=0.0

#if DEBUG.EQ.1
        write(0,*) 'DEBUG: Calculate Forward variables [ForwardAlgorithm]'
#endif

        PrecedingMarker=StartSnp
        do j=StartSnp,StopSnp
            ! Cumulative recombination fraction allows us to skip uninformative positions
            Theta=Theta+Thetas(j-1)-Theta*Thetas(j-1)
            ! Skip over uninformative positions to save time
            if ((GenosHmmMaCH(CurrentInd,j)/=MISSING).or.(j==inputParams%nsnp)) then
                call Transpose(j,PrecedingMarker,Theta)
                call ConditionOnData(CurrentInd,j)
                PrecedingMarker=j
                Theta=0.0
            endif
        enddo

    end subroutine ForwardAlgorithmForSegment

    !######################################################################
    subroutine Transpose(CurrentMarker,PrecedingMarker,Theta)
        ! Calculates the probability of get a particular state at CurrentMarker
        ! from any other state at PrecedingMarker using the transition probabilities.
        ! It basically calculates the first term of the forward variable.

        use AlphaHmmInMod
        use GlobalVariablesHmmMaCH
        implicit none

        integer,intent(in) :: CurrentMarker,PrecedingMarker
        double precision,intent(in) :: Theta
        type(AlphaHmmInput), pointer :: inputParams
        !Local variables
        integer :: i,j,Index
        double precision, dimension(:), allocatable :: Marginals
        double precision :: Summer,NoChange,OneChange,TwoChange

        inputParams => defaultInput

        allocate(Marginals((inputParams%nhapinsubh)))
#if DEBUG.EQ.1
        write(0,*) 'DEBUG: [Transpose]'
#endif

        if (Theta==0.0) then
            ForwardProbs(:,CurrentMarker)=ForwardProbs(:,PrecedingMarker)
            return
        endif

        Summer=0.0
        Index=0
        Marginals(:)=0.0

        ! Calculate the sum of all forward variables and the marginal of a
        ! given state.
#if DEBUG.EQ.1
        write(0,*) 'DEBUG: Calculate the acumulate sum of forward variables and the marginals [Transpose]'
#endif
        do i=1,inputParams%nhapinsubh
            do j=1,i-1
                Index=Index+1
                Summer=Summer+ForwardProbs(Index,PrecedingMarker)

                ! Given two states sharing a haplotype, there is only one way to
                ! get one state from the other changing only one haplotype. This
                ! is valid for whatever the shared haplotype. Let's consider the
                ! the sum of all the forward variables corresponding
                ! to the states the state (k1,k2) can be reached from changing
                ! only the first haplotype: Sum{h=1,K}(h,k2). In a similar way,
                ! let's define this other Sum{h=1,k}(k1,h). Since the states
                ! are symmetric and we are working in the upper triangular matrix
                ! both marginal are the same, and each forward variable is summed
                ! twice. Each sum can be called marginal in the first or second
                ! haplotype. Again, since we consider the upper triangular
                ! matrix of the forward variables, we need to define the marginals
                ! in a different way: Sum{h=1,k1}(h,k2)+Sum{h=k1,K}(k1,k). In
                ! this way, the sum of the two marginals will give as a result
                ! that each forward variable is summed twice.
                Marginals(i)=Marginals(i)+ForwardProbs(Index,PrecedingMarker)
                Marginals(j)=Marginals(j)+ForwardProbs(Index,PrecedingMarker)
            enddo
            Index=Index+1
            Summer=Summer+ForwardProbs(Index,PrecedingMarker)

            ! The state (k,k) has to be summed twice
            Marginals(i)=Marginals(i)+(ForwardProbs(Index,PrecedingMarker)*2.0)
        enddo

#if DEBUG.EQ.1
        write(0,*) 'DEBUG: Calculate Probabilities of number of changes [Transpose]'
#endif
        ! NOTE: to the transition probabilities:
        !   (1-Theta) is the probability that a haplotype does not change =>
        !       => Theta is the probability that a haplotype does change =>
        !       => Theta/inputParams%nhapinsubh is the probability that a haplotype change
        !          to a particular haplotype
        !
        ! Given those probabilities, then:
        !   * None hapltoyped have changed
        NoChange=(1.0-Theta)*(1.0-Theta)
        !   * Only one haplotype has changed
        OneChange=(1.0-Theta)*Theta/inputParams%nhapinsubh
        !   * The two haplotypes have changed
        TwoChange=Summer*Theta*Theta/(inputParams%nhapinsubh*inputParams%nhapinsubh)

        !Automatically rescale likelihoods when they get too small
#if DEBUG.EQ.1
        write(0,*) 'DEBUG: Rescale [Transpose]'
#endif
        if (Summer < 1e-15) then
            NoChange=NoChange*1e30
            OneChange=OneChange*1e30
            TwoChange=TwoChange*1e30
        endif

        ! This final loop actually transposes the probabilities for each state,
        ! that is, calculates the final probabilities of getting a particular state.
        Index=0
#if DEBUG.EQ.1
        write(0,*) 'DEBUG: Calculate the final probabilities [Transpose]'
#endif
        do i=1,inputParams%nhapinsubh
            do j=1,i-1
                Index=Index+1
                ! The new forward probability will be the probability of no change,
                ! plus the probability of one change, plus two times the
                ! probability of two changes (prob of (k1,k2) and prob of (k2,k1))
                ForwardProbs(Index,CurrentMarker)&
                    = (ForwardProbs(Index,PrecedingMarker)*NoChange)&
                    + (Marginals(i)*OneChange)&
                    + (Marginals(j)*OneChange)&
                    + (2*TwoChange)
            enddo
            Index=Index+1
            ForwardProbs(Index,CurrentMarker)&
                = (ForwardProbs(Index,PrecedingMarker)*NoChange)&
                + (Marginals(i)*OneChange)&
                + (2*TwoChange)
        enddo

        deallocate(Marginals)

    end subroutine Transpose

    !######################################################################
    subroutine ConditionOnData(CurrentInd,Marker)
        ! Introduce the emission probabilities (the probability of observe a
        ! given genotype) into the forward variable of the HMM.
        ! It basically calculate the second term of the forward variables
        use GlobalVariablesHmmMaCH
        use hmmParameters
        use AlphaHmmInMod
        implicit none

        integer, intent(in) :: CurrentInd, Marker
        type(AlphaHmmInput), pointer :: inputParams
        ! Local variables
        integer :: i, j, Index, genotypeInt, RefAll, AltAll
        double precision :: Factors(0:1), cond_probs(0:2)

        inputParams => defaultInput
        ! We treat missing genotypes as uninformative about the mosaic's
        ! underlying state. If we were to allow for deletions and the like,
        ! that may no longer be true.
        ! NOTE: gentoype can be refer either to genotypes or reads if working with sequence data (NGS)
        ! GlobalInbredInd(CurrentInd)==.TRUE.
        if (defaultInput%HMMOption==RUN_HMM_NGS) then
            RefAll = pedigree%pedigree(pedigree%genotypeMap(currentInd))%ReferAllele(Marker)
            AltAll = pedigree%pedigree(pedigree%genotypeMap(currentInd))%AlterAllele(Marker)

            if (RefAll + AltAll == 0) then
                return
            endif
        else
            genotypeInt = GenosHmmMaCH(CurrentInd,Marker)
            if (genotypeInt==MISSING) then
                return
            endif
        endif

        ! Index keeps track of the states already visited. The total number
        ! of states in this chunk of code is (inputParams%nhapinsubh x (inputParams%nhapinsubh-1)/2)
        Index=0
        if (inputParams%HMMOption==RUN_HMM_NGS) then
            do i=0,2
                cond_probs(i)=Penetrance(Marker,i,0)*shotgunErrorMatrix(0,RefAll,AltAll)&
                             +Penetrance(Marker,i,1)*shotgunErrorMatrix(1,RefAll,AltAll)&
                             +Penetrance(Marker,i,2)*shotgunErrorMatrix(2,RefAll,AltAll)
            enddo
        endif

        do i=1, inputParams%nhapinsubh
            if (inputParams%HMMOption /= RUN_HMM_NGS) then
                ! Probability to observe genotype SubH(i) being the true genotype GenosHmmMaCH in locus Marker
                Factors(0) = Penetrance(Marker,SubH(i,Marker),genotypeInt)
                ! Probability to observe genotype SubH(i)+1 being the true
                ! genotype GenosHmmMaCH in locus Marker
                Factors(1) = Penetrance(Marker,SubH(i,Marker)+1,genotypeInt)
            else
                ! Probability to observe genotype SubH(i) being the true genotype GenosHmmMaCH in locus Marker
                Factors(0) = cond_probs(SubH(i,Marker))
                ! Probability to observe genotype SubH(i)+1 being the true genotype GenosHmmMaCH in locus Marker
                Factors(1) = cond_probs(SubH(i,Marker)+1)
            endif

            do j=1,i
                Index=Index+1
                ForwardProbs(Index,Marker)=&
                    ForwardProbs(Index,Marker)*Factors(SubH(j,Marker))
            enddo
        enddo

    end subroutine ConditionOnData

    !######################################################################
    subroutine CalcPenetrance
        ! Initialize the Penetration matrix of the HMM as the emission
        ! probabilities matrix given in Appendix of Li et al. (2010)

        use omp_lib
        use GlobalVariablesHmmMaCH
        use hmmParameters
        use AlphaHmmInMod
        implicit none

        integer :: j, nprocs
         type(AlphaHmmInput), pointer :: inputParams
        inputParams => defaultInput

        ! Penetrance(j,i,k) = P(P_j|S_j)
        ! i = G_j = {0,1,2}
        ! k = T(S_j) = T(x_j) + T(y_j) = {0,1,2}
        ! allocate(Penetrance(inputParams%nsnp,0:2,0:2))

        nprocs = OMP_get_num_procs()

        !$OMP PARALLEL DO DEFAULT(SHARED)
        do j=1,inputParams%nsnp
            Penetrance(j,0,0)=(1.0-Epsilon(j))**2
            Penetrance(j,0,1)=2.0*(1.0-Epsilon(j))*Epsilon(j)
            Penetrance(j,0,2)=(Epsilon(j)**2)
            Penetrance(j,1,0)=(1.0-Epsilon(j))*Epsilon(j)
            Penetrance(j,1,1)=((1.0-Epsilon(j))**2)+(Epsilon(j)**2)
            Penetrance(j,1,2)=(1.0-Epsilon(j))*Epsilon(j)
            Penetrance(j,2,0)=(Epsilon(j)**2)
            Penetrance(j,2,1)=2.0*(1.0-Epsilon(j))*Epsilon(j)
            Penetrance(j,2,2)=(1.0-Epsilon(j))**2
        enddo
        !$OMP END PARALLEL DO

    end subroutine CalcPenetrance

    !######################################################################
    subroutine SetUpPrior
        ! Set up de initial state distribution, that is, the probability that
        ! the sequence starts with the state Sj:
        !   PIj = P(t1=Sj) = ForwardProb(j,1)


        use AlphaHmmInMod
        use omp_lib
        use GlobalVariablesHmmMaCH
        use hmmParameters

        implicit none
        type(AlphaHmmInput), pointer :: inputParams
        integer :: i,j,state
        double precision :: prior

        inputParams => defaultInput
        prior=1.0/(inputParams%nhapinsubh*inputParams%nhapinsubh)

        ! Initially, every state is equally possible
        !ForwardProbs(:,1)=1.0/(inputParams%nhapinsubh*inputParams%nhapinsubh)

        ! NOTE: In MaCH code, FordwardProbs is the array variable leftMatrices.
        !       In there, every element of the array is another array with
        !       s(s+1)/2 elements, as S(i,j)=S(j,i); where s=inputParams%nhapinsubh.
        !       So the probability of each state is two times 1 divided by
        !       the total number of possible states.

        state=0
        do i=1,inputParams%nhapinsubh
            do j=1,i-1
                state=state+1
                ForwardProbs(state,1)=2.0*prior
            enddo
            state=state+1
            ForwardProbs(state,1)=prior
        enddo

    end subroutine SetUpPrior

    !######################################################################
    subroutine SetUpEquations(HMM, nGenotyped, nInbred)
        use GlobalVariablesHmmMaCH
        use hmmParameters
        implicit none

        integer(kind=1),intent(in) :: HMM
        integer, intent(in) :: nGenotyped, nInbred

        ! Initialize to missing haplotypes of the whole population
        FullH=9

        if (HMM==RUN_HMM_NGS) then
            call SetUpEquationsReads(nGenotyped)
        else if (HMM==RUN_HMM_ONLY) then
            call SetUpEquationsGenotypesDiploid(nGenotyped)
        else
            call SetUpEquationsGenotypesHaploid(nGenotyped)
        endif
        if (nInbred > 0) then
            call SetUpEquationsPhaseHaploid(nGenotyped,nInbred)
        end if

        ErrorUncertainty(:)=0
        ErrorMatches(:)=0
        Errormismatches(:)=0
        Crossovers(:)=0
    end subroutine SetUpEquations

    !######################################################################
    subroutine SetUpEquationsPhaseHaploid(nGenotyped,nInbred)
        use GlobalVariablesHmmMaCH

        use random
        implicit none
        integer, intent(in) :: nGenotyped, nInbred

        integer :: i

        do i = 1, nInbred
            FullH(i+nGenotyped,:,1) = PhaseHmmMaCH(i+nGenotyped,:,1)
            FullH(i+nGenotyped,:,2) = PhaseHmmMaCH(i+nGenotyped,:,2)
        end do
    end subroutine SetUpEquationsPhaseHaploid

    !######################################################################
    subroutine SetUpEquationsGenotypesHaploid(nGenotyped)
        ! Initialize the variables and parameters of the HMM model described in
        ! Li et al. 2010, Appendix

        use GlobalVariablesHmmMaCH

        use random
        use AlphaHmmInMod
        implicit none
        integer, intent(in) :: nGenotyped

        integer :: i,j
        type(AlphaHmmInput), pointer :: inputParams

        inputParams => defaultInput
        ! If the number of phased gametes from AlphaImpute is above a threshold, then
        ! haploytpes produced from AlphaImpute are used in the model (FullH)
        if (nGametesPhased/float(2*pedigree%nGenotyped)>phasedThreshold/100.0) then
            do i=1,nGenotyped
                FullH(i,:,:)=PhaseHmmMaCH(i,:,:)

                ! Missing alleles (and individuals that are not phased) are called at random
                do j=1,inputParams%nsnp
                    if (PhaseHmmMaCH(i,j,1)==ALLELE_MISSING) then
                        if (ran1(inputParams%seed)>=0.5) then
                            FullH(i,j,1)=0
                        else
                            FullH(i,j,1)=1
                        endif
                    endif

                    if (PhaseHmmMaCH(i,j,2)==ALLELE_MISSING) then
                        if (ran1(inputParams%seed)>=0.5) then
                            FullH(i,j,2)=0
                        else
                            FullH(i,j,2)=1
                        endif
                    endif

                enddo
            enddo

        else
            ! If the number of phased gametes from AlphaImpute is below a threshold, then
            ! haplotypes of phased animals from AlphaImpute are used in the model, and
            ! haplotypes of genotyped animals are used otherwise.

            ! Initialise FullH with the Genotype information
            call SetUpEquationsGenotypesDiploid(nGenotyped)

            ! Overwrite haplotypes to use phased data in case phased haplotypes from
            ! AlphaImpute are available

            do i=1,nGenotyped      ! For every Individual in the Genotype file
                if (GlobalHmmPhasedInd(i,1)==.TRUE.) then
                    FullH(i,:,1)=PhaseHmmMaCH(i,:,1)
                    !  If there is missing information in the phased data, called allele at random
                    do j=1,inputParams%nsnp
                        if (PhaseHmmMaCH(i,j,1)==ALLELE_MISSING) then
                            if (ran1(inputParams%seed)>=0.5) then
                                FullH(i,j,1)=0
                            else
                                FullH(i,j,1)=1
                            endif
                        endif
                    enddo
                endif

                if (GlobalHmmPhasedInd(i,2)==.TRUE.) then
                    FullH(i,:,2)=PhaseHmmMaCH(i,:,2)
                    !  If there is missing information in the phased data, called allele at random
                    do j=1,inputParams%nsnp

                        if (PhaseHmmMaCH(i,j,2)==ALLELE_MISSING) then
                            if (ran1(inputParams%seed)>=0.5) then
                                FullH(i,j,2)=0
                            else
                                FullH(i,j,2)=1
                            endif
                        endif
                    enddo
                endif

            enddo
        endif
    end subroutine SetUpEquationsGenotypesHaploid

    !######################################################################
    subroutine SetUpEquationsGenotypesDiploid(nGenotyped)
        ! Initialize the variables and parameters of the HMM model described in
        ! Li et al. 2010, Appendix

        use GlobalVariablesHmmMaCH

        use random
        use AlphaHmmInMod
        implicit none
        integer, intent(in) :: nGenotyped

        integer :: i,j
        type(AlphaHmmInput), pointer :: inputParams

        inputParams => defaultInput
        !Initialise FullH
        do i=1,nGenotyped      ! For every Individual in the Genotype file
            do j=1,inputParams%nsnp      ! For each SNP

                ! Phase homozygose locus
                if (GenosHmmMaCH(i,j)==0) then
                    FullH(i,j,:)=0
                endif

                if (GenosHmmMaCH(i,j)==2) then
                    FullH(i,j,:)=1
                endif

                ! Phase heterozygose case at random
                if (GenosHmmMaCH(i,j)==1) then
                    if (ran1(inputParams%seed)>=0.5) then
                        FullH(i,j,1)=0
                        FullH(i,j,2)=1
                    else
                        FullH(i,j,1)=1
                        FullH(i,j,2)=0
                    endif
                endif

                ! If locus is not genotyped, phase each haplotype at random
                if (GenosHmmMaCH(i,j)==MISSING) then
                    if (ran1(inputParams%seed)>=0.5) then
                        FullH(i,j,1)=0
                    else
                        FullH(i,j,1)=1
                    endif
                    if (ran1(inputParams%seed)>=0.5) then
                        FullH(i,j,2)=0
                    else
                        FullH(i,j,2)=1
                    endif
                endif
            enddo
        enddo
    end subroutine SetUpEquationsGenotypesDiploid

    !###########################################################################
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Get allele frequencies from population
    !
    !> @details    This function calculates allele frequencies from the population
    !>             and returns an array with the frequencies
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date       Dec 11, 2017
    !---------------------------------------------------------------------------
    function GetAlleleFrequenciesFromPopulation() result(frequency)
        use GlobalVariablesHmmMaCH
        use AlphaHmmInMod
        implicit none

        ! Dummy Arguments
        double precision, allocatable, dimension(:) :: frequency        !< Frequencies

        ! Local Variables
        integer :: i, j, readObs, alleles, refAll, altAll
        type(AlphaHmmInput), pointer :: inputParams

        inputParams => defaultInput

        allocate(frequency(inputParams%nSnp))

        ! Read frequencies from file
        do j = 1, inputParams%nSnp
            readObs = 0
            alleles = 0
            do i=1,pedigree%nGenotyped
                refAll = pedigree%pedigree(pedigree%genotypeMap(i))%ReferAllele(j)
                altAll = pedigree%pedigree(pedigree%genotypeMap(i))%AlterAllele(j)
                readObs = readObs + refAll + altAll
                alleles = alleles + altAll
            enddo
            frequency(j) =  dble(alleles) / dble(readObs)
        end do
    end function GetAlleleFrequenciesFromPopulation

    !###########################################################################
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Get allele frequencies
    !
    !> @details    This subroutine reads allele frequencies from a file provided
    !>             by the user or calculate them from the population otherwise
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date       Dec 11, 2017
    !---------------------------------------------------------------------------
    function GetAlleleFrequencies() result(frequency)
        use GlobalVariablesHmmMaCH
        use AlphaHmmInMod
        use HmmInputMod
        implicit none

        ! Dummy Arguments
        double precision, allocatable, dimension(:) :: frequency   !< Frequencies

        ! Local Variables
        logical :: exists
        type(AlphaHmmInput), pointer :: inputParams

        inputParams => defaultInput

        allocate(frequency(inputParams%nSnp))

        inquire(file=trim(inputParams%PriorAllFreqsFile), exist=exists)
        if (exists) then
            frequency = ReadAlleleFrequenciesFromFile()
        else
            frequency = GetAlleleFrequenciesFromPopulation()
        end if

    end function GetAlleleFrequencies

    !######################################################################
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Initialize the variables and parameters of the HMM model
    !
    !> @details    Initialize the variables and parameters of the HMM model
    !>             described in Li et al. 2010 (Appendix) for the NGS data
    !
    !> @author     Roberto Antolin, roberto.antolin@roslin.ed.ac.uk
    !
    !> @date       Dec 11, 2017
    !---------------------------------------------------------------------------
    subroutine SetUpEquationsReads(nGenotyped)
        use GlobalVariablesHmmMaCH
        use random
        use AlphaHmmInMod
        implicit none

        integer, intent(in) :: nGenotyped

        integer :: i, j, RefAll, AltAll
        double precision :: prior_11, prior_12, prior_22
        double precision :: posterior_11, posterior_12, posterior_22
        double precision :: r, summ
        double precision, allocatable, dimension(:) :: frequency
        type(AlphaHmmInput), pointer :: inputParams

        inputParams => defaultInput

        ! Get allele frequencies from file
        frequency = GetAlleleFrequencies()

        !Initialise FullH
        do j=1,inputParams%nsnp      ! For each SNP

            ! Calculate prior probabitlities based on allele frequencies
            prior_11 = (1.0 - frequency(j))**2
            prior_12 = 2.0 * (1.0 - frequency(j)) * frequency(j)
            prior_22 = frequency(j)**2

            do i=1,nGenotyped
                RefAll = pedigree%pedigree(pedigree%genotypeMap(i))%ReferAllele(j)
                AltAll = pedigree%pedigree(pedigree%genotypeMap(i))%AlterAllele(j)

                posterior_11 = prior_11 * ShotgunErrorMatrix(0,RefAll,AltAll)
                posterior_12 = prior_12 * ShotgunErrorMatrix(1,RefAll,AltAll)
                posterior_22 = prior_22 * ShotgunErrorMatrix(2,RefAll,AltAll)

                summ = posterior_11 + posterior_12 + posterior_22
                if (summ==0) then
                    write(0,*) 'There is a problem here!'
                end if

                posterior_11 = posterior_11 / summ
                posterior_12 = posterior_12 / summ
                posterior_22 = posterior_22 / summ

                if (posterior_11 > 0.99) then
                    GenosHmmMaCH(i,j) = 0
                    FullH(i,j,:) = 0
                else if (posterior_12 > 0.99) then
                    GenosHmmMaCH(i,j) = 1
                    if (ran1(inputParams%seed)<0.5) then
                        FullH(i,j,1) = 0
                        FullH(i,j,2) = 1
                    else
                        FullH(i,j,1) = 1
                        FullH(i,j,2) = 0
                    endif
                else if (posterior_22 > 0.99) then
                    GenosHmmMaCH(i,j) = 2
                    FullH(i,j,:) = 1
                else
                    if ((RefAll + AltAll)  == 0) then
                        GenosHmmMaCH(i,j) = MISSING
                    end if
                    r = ran1(inputParams%seed)
                    if (r < posterior_11) then
                        FullH(i,j,:) = 0
                    else if (r < posterior_11 + posterior_12) then
                        if (ran1(inputParams%seed)<0.5) then
                            FullH(i,j,1) = 0
                            FullH(i,j,2) = 1
                        else
                            FullH(i,j,1) = 1
                            FullH(i,j,2) = 0
                        endif
                    else
                        FullH(i,j,:)=1
                    endif
                end if

                if (GlobalInbredInd(i) == .TRUE.) then
                    if (RefAll + AltAll == 0) then
                        GenosHmmMaCH(i,j) = MISSING
                        PhaseHmmMaCH(i,j,:) = ALLELE_MISSING
                    end if

                    if (posterior_11 > posterior_22) then
                        GenosHmmMaCH(i,j) = 0
                        FullH(i,j,:) = 0
                        PhaseHmmMaCH(i,j,:) = 0
                    elseif (posterior_11 < posterior_22) then
                        GenosHmmMaCH(i,j) = 2
                        FullH(i,j,:) = 1
                        PhaseHmmMaCH(i,j,:) = 1
                    else
                        GenosHmmMaCH(i,j) = 1
                        if (ran1(inputParams%seed)<0.5) then
                            PhaseHmmMaCH(i,j,1) = 0
                            PhaseHmmMaCH(i,j,2) = 1
                            FullH(i,j,:) = 0
                            FullH(i,j,:) = 1
                        else
                            PhaseHmmMaCH(i,j,1) = 1
                            PhaseHmmMaCH(i,j,2) = 0
                            FullH(i,j,:) = 1
                            FullH(i,j,:) = 0
                        endif
                    endif
                end if

            enddo
        enddo

    end subroutine SetUpEquationsReads

    !######################################################################
    subroutine UpdateThetas

        use GlobalVariablesHmmMaCH
        use hmmParameters
        use AlphaHmmInMod

        implicit none

        integer :: i, BaseCount=1, BaseIntervals=0
        double precision :: BaseRates, BaseCrossovers=1

        double precision :: Scale
        type(AlphaHmmInput), pointer :: inputParams
        inputParams => defaultInput

        Scale=1.0/(nIndHmmMaCH*2)
        BaseCount=1
        BaseIntervals=0

        BaseCrossovers=1.0

        ! First we estimate a base line rate to be applied to intervals with
        ! 0 or 1 observed "crossovers"
        do i=1,inputParams%nsnp-1
            if (Crossovers(i)<=1) then
                BaseCount=BaseCount+Crossovers(i)
                BaseIntervals=BaseIntervals+1
            endif
        enddo

        if (BaseIntervals==0) then
            BaseRates=BaseCount*Scale
        else
            BaseRates=BaseCount*Scale/BaseIntervals
        endif

        ! Then we update the rate for each interval using either the number
        ! of observed crossovers (if > 1) or the baseline rate
        do i=1,inputParams%nsnp-1
            if (Crossovers(i)>1) then
                Thetas(i)=Crossovers(i)*Scale
            else
                Thetas(i)=BaseRates
            endif
        enddo
        !print*,Thetas(inputParams%nsnp-1)

    end subroutine UpdateThetas

    !######################################################################
    subroutine UpdateErrorRate(rate)
        ! Group markers into those with low error rates, which are estimated
        ! as a group, and those with high error rates, which are estimated
        ! individually

        use GlobalVariablesHmmMaCH
        use AlphaHmmInMod
        implicit none

        double precision,intent(out) :: rate

        ! Local variables
        integer :: i,matches=0,mmatches,uncertain=0
          type(AlphaHmmInput), pointer :: inputParams
        inputParams => defaultInput

        mmatches = 0
        rate=0.0
        do i=1,inputParams%nsnp
            if (Errormismatches(i)<=2) then
                matches=matches+ErrorMatches(i)
                mmatches=mmatches+Errormismatches(i)
                uncertain=uncertain+ErrorUncertainty(i)
            else
                call UpdateError(ErrorMatches(i), Errormismatches(i), ErrorUncertainty(i), rate)
                call SetPenetrance(i,rate)
            endif
        enddo

        call UpdateError(matches, mmatches, uncertain, rate)

        do i=1,inputParams%nsnp
            if (Errormismatches(i)<=2) call SetPenetrance(i, rate)
        enddo

    end subroutine UpdateErrorRate

    !######################################################################
    subroutine SetPenetrance(marker, Err)
        use GlobalVariablesHmmMaCH, only: Epsilon, Penetrance

        implicit none

        integer,intent(in) :: marker
        double precision, intent(in) :: Err

        Epsilon(marker) = Err

        Penetrance(marker,0,0)=(1.0-Err)**2
        Penetrance(marker,0,1)=2.0*(1.0-Err)*Err
        Penetrance(marker,0,2)=Err**2
        Penetrance(marker,1,0)=(1.0-Err)*Err
        Penetrance(marker,1,1)=((1.0-Err)**2)+(Err**2)
        Penetrance(marker,1,2)=(1.0-Err)*Err
        Penetrance(marker,2,0)=Err**2
        Penetrance(marker,2,1)=2.0*(1.0-Err)*Err
        Penetrance(marker,2,2)=(1.0-Err)**2

    end subroutine SetPenetrance


    !######################################################################
    subroutine TotalCrossovers(Total)

        use GlobalVariablesHmmMaCH, only : crossovers
        use AlphaHmmInMod
        implicit none


        integer,intent(out) :: Total
        type(AlphaHmmInput), pointer :: inputParams
         integer :: i
        inputParams => defaultInput
        ! Local variables


        Total=0
        do i=1,inputParams%nsnp-1
            Total=Total+Crossovers(i)
        enddo

    end subroutine TotalCrossovers

    !######################################################################
    subroutine ResetCrossovers
            use GlobalVariablesHmmMaCH, only : crossovers
        implicit none

        Crossovers(:)=0
        call ResetErrors

    end subroutine ResetCrossovers

    !######################################################################
    subroutine UpdateError(matches,mmatches, uncertain, rate)

            use GlobalVariablesHmmMaCH
        implicit none

        integer,intent(in) :: matches,mmatches,uncertain
        double precision, intent(out) :: rate

        ! Local variables
        double precision :: previous=0.0, ratio

        rate=0.0    ! Just in case...
        if (matches+mmatches>0) then
            rate=mmatches/dble(matches+mmatches)
            if (uncertain>0) then
                do while((rate>1e-10).and.(abs(rate-previous) > rate*1e-4))
                    ratio=rate*rate/(rate*rate+(1.0-rate)*(1.0-rate))
                    previous=rate
                    rate=(mmatches+ratio*uncertain*2.0)&
                        /(matches+mmatches+uncertain*2)
                enddo
            endif
        else if (uncertain>0) then
            rate=0.0
        endif

    end subroutine UpdateError

    !######################################################################
    subroutine ResetErrors
        use GlobalVariablesHmmMaCH, only : ErrorUncertainty, Errormismatches,ErrorMatches

        implicit none

        ErrorUncertainty(:)=0
        ErrorMatches(:)=0
        Errormismatches(:)=0

    end subroutine ResetErrors

    !######################################################################
    subroutine GetErrorRate(mean)

        use GlobalVariablesHmmMaCH
        use hmmHaplotyper
        use AlphaHmmInMod
        implicit none
        double precision, intent(out) :: mean

        double precision :: ErrorRate
        integer :: i
          type(AlphaHmmInput), pointer :: inputParams
        inputParams => defaultInput

        mean = 0.0
        do i=1,inputParams%nsnp
            call GetErrorRatebyMarker(i, ErrorRate)
            mean = mean + ErrorRate
        enddo
        mean = mean / inputParams%nsnp

    end subroutine GetErrorRate

    !######################################################################
    subroutine ExtractTemplate(HMM, forWhom, nGenotyped)
        use GlobalVariablesHmmMaCH

        use omp_lib
        use random
        use Par_Zig_mod

        use hmmParameters

        implicit none
        integer(kind=1) ::  HMM
        integer, intent(in) ::forWhom, nGenotyped

        integer :: Shuffle1(nGenotyped), Shuffle2(nGenotyped)
        integer :: thread

        ! Create vectors of random indexes
#if DEBUG.EQ.1
        write(0,*) "DEBUG: Suffle individuals: RandomOrderPar [MaCHForInd]"
#endif
        thread=omp_get_thread_num()
        call RandomOrderPar(Shuffle1,nGenotyped,thread)
        call RandomOrderPar(Shuffle2,nGenotyped,thread)

        ! Extract haps template (SubH) ...
        if (HMM==RUN_HMM_ONLY) then
            call ExtractTemplateHaps(forWhom,Shuffle1,Shuffle2)
        else
            if (nGametesPhased/float(2*pedigree%nGenotyped)>phasedThreshold/100.0) then
                ! If the number of phased gametes with AlphaImpute is above
                ! a threshold, then template is populated with the phased data
                call ExtractTemplateByHaps(forWhom,Shuffle1,Shuffle2)
            else
                ! Otherwise, the template is populated with haplotypes at random
                ! from all the HD animals
                call ExtractTemplateHaps(forWhom,Shuffle1,Shuffle2)
            endif
        endif

        ! ... selecting pairs of haplotypes at random
        !call ExtractTemplateHapsByAnimals(CurrentInd,Shuffle1)
    end subroutine ExtractTemplate

    !######################################################################
    subroutine ExtractTemplateHaps(CurrentInd,Shuffle1,Shuffle2)
        ! Set the Template of Haplotypes used in the HMM model.
        ! It takes the haplotypes from HD animals
        ! It differentiates between paternal (even) and maternal (odd) haps
        use AlphaHmmInMod
        use GlobalVariablesHmmMaCH


        integer, intent(in) :: Shuffle1(nIndHmmMaCH), Shuffle2(nIndHmmMaCH),CurrentInd

        ! Local variables
        integer :: HapCount, ShuffleInd1, ShuffleInd2
        type(AlphaHmmInput), pointer :: inputParams

        inputParams => defaultInput

        HapCount=0
        ShuffleInd1=0
        ShuffleInd2=0

#if DEBUG.EQ.1
        write(0,*) 'DEBUG: Create Haplotypes Template [ExtractTemplateHaps]'
#endif
        ! While the maximum number of haps in the template haplotypes set,
        ! H, is not reached...
        do while (HapCount<inputParams%nhapinsubh)
            if (mod(HapCount,2)==0) then

                ShuffleInd1=ShuffleInd1+1
                if (ShuffleInd1 > nIndHmmMaCH) then
                    write(error_unit,*) "WARNING: not enough haplotypes added as there are too few hd amimals"
                    exit
                endif
                ! Select the paternal haplotype if the individual it belongs
                ! to is genotyped and it is not the current individual
                if ((Shuffle1(ShuffleInd1)/=CurrentInd)&
                    .and.(GlobalHmmHDInd(Shuffle1(ShuffleInd1))==1)) then
                    HapCount=HapCount+1
                    SubH(HapCount,:)=FullH(Shuffle1(ShuffleInd1),:,1)
                endif
            else
                ShuffleInd2=ShuffleInd2+1
                if (ShuffleInd2 > nIndHmmMaCH) then
                    write(error_unit,*) "WARNING: not enough haplotypes added as there are too few hd amimals"
                    exit
                endif
                ! Select the maternal haplotype if the individual it belongs
                ! too is genotyped and it is not the current individual
                if ((Shuffle2(ShuffleInd2)/=CurrentInd)&
                    .and.(GlobalHmmHDInd(Shuffle2(ShuffleInd2))==1)) then
                    HapCount=HapCount+1
                    SubH(HapCount,:)=FullH(Shuffle2(ShuffleInd2),:,2)
                endif
            endif
        enddo

    end subroutine ExtractTemplateHaps

    !######################################################################
    subroutine ExtractTemplateByHaps(CurrentInd,Shuffle1,Shuffle2)
        ! Set the Template of Haplotypes used in the HMM model.
        ! It takes the haplotypes produced by AlphaImpute
        ! It differentiates between paternal (even) and maternal (odd) haps

        use GlobalVariablesHmmMaCH
        use AlphaHmmInMod

        integer, intent(in) :: Shuffle1(nIndHmmMaCH), Shuffle2(nIndHmmMaCH),CurrentInd

        ! Local variables
        integer :: HapCount, ShuffleInd1, ShuffleInd2
        type(AlphaHmmInput), pointer :: inputParams

        inputParams => defaultInput
        HapCount=0
        ShuffleInd1=0
        ShuffleInd2=0

#if DEBUG.EQ.1
        write(0,*) 'DEBUG: Create Haplotypes Template [ExtractTemplateByHaps]'
#endif
        ! While the maximum number of haps in the template haplotypes set,
        ! H, is not reached...
        do while (HapCount<inputParams%nhapinsubh)
            if (mod(HapCount,2)==0) then
                ShuffleInd1=ShuffleInd1+1

                ! Select the paternal haplotype if the individual it belongs
                ! to is genotyped and it is not the current individual
                if (Shuffle1(ShuffleInd1)/=CurrentInd) then

                    if (GlobalHmmPhasedInd(Shuffle1(ShuffleInd1),1)==.TRUE.) then

                        HapCount=HapCount+1
                        SubH(HapCount,:)=FullH(Shuffle1(ShuffleInd1),:,1)
                    endif
                endif
            else
                ShuffleInd2=ShuffleInd2+1

                ! Select the maternal haplotype if the individual it belongs
                ! too is genotyped and it is not the current individual
                if (Shuffle2(ShuffleInd2)/=CurrentInd) then

                    if (GlobalHmmPhasedInd(Shuffle2(ShuffleInd2),2)==.TRUE.) then
                        HapCount=HapCount+1
                        SubH(HapCount,:)=FullH(Shuffle2(ShuffleInd2),:,2)
                    endif
                endif
            endif
        enddo

    end subroutine ExtractTemplateByHaps

    !######################################################################
    subroutine ExtractTemplateHapsByAnimals(CurrentInd,Shuffle)
        ! Extract the Template of Haplotypes used in the HMM model.
        ! It consideres the two haps of a given individual.


        use GlobalVariablesHmmMaCH
        use AlphaHmmInMod


        integer, intent(in) :: Shuffle(nIndHmmMaCH),CurrentInd
        type(AlphaHmmInput), pointer :: inputParams
        ! Local variables
        integer :: HapCount, ShuffleInd

        inputParams => defaultInput
        HapCount=1
        ShuffleInd=0

#if DEBUG.EQ.1
        write(0,*) 'DEBUG: Create Haplotypes Template [ExtractTemplateHapsByAnimals]'
#endif
        ! While the maximum number of haps in the template haplotypes set,
        ! H, is not reached...
        do while (HapCount<inputParams%nhapinsubh)
            ShuffleInd=ShuffleInd+1

            ! Select the paternal and maternal haplotypes from one individual
            ! if this individual is not being phased/imputed
            if ((Shuffle(ShuffleInd)/=CurrentInd)&
                .and.(GlobalHmmHDInd(ShuffleInd)==1)) then
                SubH(HapCount,:)=FullH(Shuffle(ShuffleInd),:,1)
                SubH(HapCount+1,:)=FullH(Shuffle(ShuffleInd),:,2)
                HapCount=HapCount+2
            endif
        enddo

    end subroutine ExtractTemplateHapsByAnimals

    !######################################################################
    subroutine SetShotgunError(ErrorRate)
        use ConstantModule
        use GlobalVariablesHmmMaCH
        use hmmParameters

        implicit none

        double precision, intent(in) :: ErrorRate

        ! Local variables
        double precision :: ProdFactTmp
        double precision :: DFactorialInLog(0:(MAX_READS_COUNT-1)**2)
        integer :: i,k

        do k=0,MAX_READS_COUNT-1
            do i=0,MAX_READS_COUNT-1
                ProdFactTmp=DFactorialInLog(k+i)-(DFactorialInLog(i)+DFactorialInLog(k))
                ShotgunErrorMatrix(0,k,i)=exp(ProdFactTmp+(dfloat(i)*log(ErrorRate))+(dfloat(k)*log(1.0-ErrorRate)))
                ShotgunErrorMatrix(1,k,i)=exp(ProdFactTmp+(dfloat(k+i)*log(0.5)))
                ShotgunErrorMatrix(2,k,i)=exp(ProdFactTmp+(dfloat(i)*log(1.0-ErrorRate))+(dfloat(k)*log(ErrorRate)))
            enddo
        enddo

    end subroutine SetShotgunError

    !######################################################################
    subroutine ImputeAllelesNGS(CurrentInd,CurrentMarker,State1,State2)
        use GlobalVariablesHmmMaCH

        use hmmHaplotyper
        use omp_lib
        use Par_Zig_mod

        implicit none

        integer, intent(in) :: CurrentInd, CurrentMarker, State1, State2

        ! Local variables
        integer :: copied1, copied2, Thread, imputed1, imputed2, Differences, RefAll, AltAll
        double precision :: posterior_11, posterior_12, posterior_22, summ, random, rate

        Thread = omp_get_thread_num()

        copied1 = SubH(State1,CurrentMarker)
        copied2 = SubH(State2,CurrentMarker)

        RefAll = pedigree%pedigree(pedigree%genotypeMap(currentInd))%ReferAllele(CurrentMarker)
        AltAll = pedigree%pedigree(pedigree%genotypeMap(currentInd))%AlterAllele(CurrentMarker)

        posterior_11 = Penetrance(CurrentMarker,copied1+copied2,0)*shotgunErrorMatrix(0,RefAll,AltAll)
        posterior_12 = Penetrance(CurrentMarker,copied1+copied2,1)*shotgunErrorMatrix(1,RefAll,AltAll)
        posterior_22 = Penetrance(CurrentMarker,copied1+copied2,2)*shotgunErrorMatrix(2,RefAll,AltAll)

        summ = posterior_11 + posterior_12 + posterior_22

        posterior_11 = posterior_11 / summ
        posterior_22 = posterior_22 / summ

        random = par_uni(Thread)

        if (random < posterior_11) then
            FullH(CurrentInd,CurrentMarker,1) = 0
            FullH(CurrentInd,CurrentMarker,2) = 0

        elseif (random < posterior_11 + posterior_22) then
            FullH(CurrentInd,CurrentMarker,1) = 1
            FullH(CurrentInd,CurrentMarker,2) = 1

        elseif (copied1 /= copied2) then
            call GetErrorRatebyMarker(CurrentMarker, rate)

            if (par_uni(Thread) < rate*rate / ((rate*rate) + (1-rate)*(1-rate))) then
                copied1 = 1 - copied1
                copied2 = 1 - copied2
            endif
            FullH(CurrentInd,CurrentMarker,1) = copied1
            FullH(CurrentInd,CurrentMarker,2) = copied2

        else
            if (par_uni(Thread)<0.5) then
                FullH(CurrentInd,CurrentMarker,1) = 0
                FullH(CurrentInd,CurrentMarker,2) = 1
            else
                FullH(CurrentInd,CurrentMarker,1) = 1
                FullH(CurrentInd,CurrentMarker,2) = 0
            endif
        endif

        imputed1 = FullH(CurrentInd,CurrentMarker,1)
        imputed2 = FullH(CurrentInd,CurrentMarker,2)

        Differences = abs(copied1 - imputed1) + abs(copied2 - imputed2)
        ! count the number of alleles matching
        !$OMP ATOMIC
        ErrorMatches(CurrentMarker)=ErrorMatches(CurrentMarker)+(2-Differences)
        ! count the number of mismatching alleles
        !$OMP ATOMIC
        Errormismatches(CurrentMarker)=Errormismatches(CurrentMarker)+Differences

    end subroutine ImputeAllelesNGS

    !######################################################################################################################################################

    real(kind=8) function DFactorialInLog(n)

        implicit none
        integer,intent(in) :: n

        ! Local variables
        integer :: i
        real(kind=8) :: Ans

        Ans=0.0
        if (n==0) then
            Ans=1.0
        else
            do i=1,n
                Ans=Ans+log(dfloat(i))
            enddo
        endif

        DFactorialInLog=Ans

    end function DFactorialInLog

    !######################################################################################################################################################
end MODULE hmmModule
