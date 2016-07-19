!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	  UNIVERSIDAD NACIONAL DE COLOMBIA"
!!	  PROF. ANDRES REYES GROUP"
!!	  http://www.qcc.unal.edu.co"
!!	
!!	  UNIVERSIDAD DE GUADALAJARA"
!!	  PROF. ROBERTO FLORES GROUP"
!!	  http://www.cucei.udg.mx/~robertof"
!!	
!!	AUTHORS
!!		E.F. POSADA. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		S.A. GONZALEZ. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		F.S. MONCADA. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		J. ROMERO. UNIVERSIDAD NACIONAL DE COLOMBIA
!!
!!	CONTRIBUTORS
!!		N.F.AGUIRRE. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		GABRIEL MERINO. UNIVERSIDAD DE GUANAJUATO
!!   		J.A. CHARRY UNIVERSIDAD NACIONAL DE COLOMBIA
!!
!!
!!		Todos los derechos reservados, 2011
!!
!!******************************************************************************
                
module ConfigurationInteraction_
  use Exception_
  use Matrix_
  use Vector_
!  use MolecularSystem_
  use Configuration_
  use ReadTransformedIntegrals_
  use MolecularSystem_
  use String_
  use IndexMap_
  implicit none
      
  !>
  !! @brief Configuration Interaction Module, works in spin orbitals
  !!
  !! @author felix
  !!
  !! <b> Creation data : </b> 07-24-12
  !!
  !! <b> History change: </b>
  !!
  !!   - <tt> 07-24-12 </tt>: Felix Moncada ( fsmoncadaa@unal.edu.co )
  !!        -# description.
  !!   - <tt> 07-09-16 </tt>: Jorge Charry ( jacharrym@unal.edu.co )
  !!        -# Add CIS, and Fix CISD.
  !!   - <tt> MM-DD-YYYY </tt>:  authorOfChange ( email@server )
  !!        -# description
  !!
  !<
  type, public :: ConfigurationInteraction
     logical :: isInstanced
     type(matrix) :: hamiltonianMatrix
     type(matrix) :: eigenVectors
     integer :: numberOfConfigurations
     type(vector) :: numberOfOccupiedOrbitals
     type(vector) :: numberOfOrbitals
     type(vector) :: eigenvalues
     type(vector) :: lambda !!Number of particles per orbital, module only works for 1 or 2 particles per orbital
     type(matrix), allocatable :: fourCenterIntegrals(:,:)
     type(matrix), allocatable :: twoCenterIntegrals(:)
     type(matrix), allocatable :: FockMatrix(:)
     type(vector), allocatable :: energyofmolecularorbitals(:)
     type(configuration), allocatable :: configurations(:)
     real(8) :: totalEnergy

     character(20) :: level

  end type ConfigurationInteraction

  type, public :: HartreeFock
        real(8) :: totalEnergy
        real(8) :: puntualInteractionEnergy
        type(matrix) :: coefficientsofcombination 
        type(matrix) :: HcoreMatrix 
  end type HartreeFock
  
  type(ConfigurationInteraction) :: ConfigurationInteraction_instance
  type(HartreeFock) :: HartreeFock_instance

  public :: &
       ConfigurationInteraction_constructor, &
       ConfigurationInteraction_destructor, &
       ConfigurationInteraction_getTotalEnergy, &
       ConfigurationInteraction_run, &
       ConfigurationInteraction_show

  private
contains


  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine ConfigurationInteraction_constructor(level)
    implicit none
    character(*) :: level

    integer :: numberOfSpecies
    integer :: i,j,k,m,n,p,q,c,r,s
    integer :: isLambdaEqual1,lambda
    type(vector) :: order
    type(vector) :: occupiedCode
    type(vector) :: unoccupiedCode
    real(8) :: totalEnergy

    character(50) :: wfnFile
    integer :: wfnUnit
    character(50) :: nameOfSpecie
    integer :: numberOfContractions
    character(50) :: arguments(2)

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

    !! Open file for wavefunction
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

    !! Load results...
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=HartreeFock_instance%totalEnergy, &
         arguments=["TOTALENERGY"])
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=HartreeFock_instance%puntualInteractionEnergy, &
         arguments=["PUNTUALINTERACTIONENERGY"])

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    do i=1, numberOfSpecies
        nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( i ) )
        numberOfContractions = MolecularSystem_getTotalNumberOfContractions( i )

        arguments(2) = nameOfSpecie
        arguments(1) = "HCORE"
        HartreeFock_instance%HcoreMatrix  = &
                  Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
                  columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

        arguments(1) = "COEFFICIENTS"
        HartreeFock_instance%coefficientsofcombination = &
                  Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
                  columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))
    end do

    ConfigurationInteraction_instance%isInstanced=.true.
    ConfigurationInteraction_instance%level=level
    ConfigurationInteraction_instance%numberOfConfigurations=0

    call Vector_constructor (ConfigurationInteraction_instance%numberOfOccupiedOrbitals, numberOfSpecies)
    call Vector_constructor (ConfigurationInteraction_instance%numberOfOrbitals, numberOfSpecies)
    call Vector_constructor (ConfigurationInteraction_instance%lambda, numberOfSpecies)

    do i=1, numberOfSpecies
       !! We are working in spin orbitals not in spatial orbitals!
       ConfigurationInteraction_instance%lambda%values(i) = MolecularSystem_getLambda( i )
       ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)=MolecularSystem_getOcupationNumber( i )* ConfigurationInteraction_instance%lambda%values(i)
       ConfigurationInteraction_instance%numberOfOrbitals%values(i)=MolecularSystem_getTotalNumberOfContractions( i )* ConfigurationInteraction_instance%lambda%values(i)
       !!Uneven occupation number = alpha
       !!Even occupation number = beta     
    end do

    select case ( trim(level) )

    case ( "CIS" )

       !!Count
       
       !!Ground State
       c=1

       do i=1, numberOfSpecies

          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          !!Singles
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
             do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )

                   if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                     c=c+1
                   end if
                end do
             end do
          end if


       end do
       !print *, "number OF confi", c 
       ConfigurationInteraction_instance%numberOfConfigurations = c
       allocate (ConfigurationInteraction_instance%configurations(ConfigurationInteraction_instance%numberOfConfigurations) )

       call Matrix_Constructor(ConfigurationInteraction_instance%hamiltonianMatrix, int(ConfigurationInteraction_instance%numberOfConfigurations,8),int(ConfigurationInteraction_instance%numberOfConfigurations,8),0.0_8)


    case ( "CISD" )

  !!Count
       
       !!Ground State
       c=1

       do i=1, numberOfSpecies

          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          !!Singles
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
             do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )

                   if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                     c=c+1
                   end if
                end do
             end do
          end if


       end do
       !print *, "number OF confi", c 


       do i=1, numberOfSpecies
   !!Doubles of the same specie

          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
             do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                do n=m+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   do p= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                      !print *, "dou", m,n,p
                      !if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                        do q=p+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                   !do q=ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)+1, ConfigurationInteraction_instance%numberOfOrbitals%values(i) 
                      !print *, "dou", m,n,p
                      if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) .and. (p /= q) ) then !! alpha -> alpha, beta -> beta

                             c=c+1
                           !end if
                      end if
                        end do

                   end do
                end do
             end do
          end if

          !!Doubles of different species
          do j=i+1, numberOfSpecies
             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
                do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                   do n=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                      do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                         do q= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                            c=c+1
                         end do
                      end do
                   end do
                end do
             end if
          end do


       end do !! species

!      c = 6
       !print *, "number OF doubles confi", c 

       ConfigurationInteraction_instance%numberOfConfigurations = c
       allocate (ConfigurationInteraction_instance%configurations(ConfigurationInteraction_instance%numberOfConfigurations) )

       call Matrix_Constructor(ConfigurationInteraction_instance%hamiltonianMatrix, int(ConfigurationInteraction_instance%numberOfConfigurations,8),int(ConfigurationInteraction_instance%numberOfConfigurations,8),0.0_8)


    case ( "CIDD" )

       !!Count
       
       !!Ground State
       c=1

       do i=1, numberOfSpecies

          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
   !!Doubles of the same specie
          if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
             do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                do n=m+1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i))
                        do q=p+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i))
                      if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) .and. (p /= q) ) then !! alpha -> alpha, beta -> beta

                             c=c+1
                      end if
                        end do

                   end do
                end do
             end do
          end if


       end do

       ConfigurationInteraction_instance%numberOfConfigurations = c
       allocate (ConfigurationInteraction_instance%configurations(ConfigurationInteraction_instance%numberOfConfigurations) )

       call Matrix_Constructor(ConfigurationInteraction_instance%hamiltonianMatrix, int(ConfigurationInteraction_instance%numberOfConfigurations,8),int(ConfigurationInteraction_instance%numberOfConfigurations,8),0.0_8)

    case ( "FCI" )

  !!Count
       
       !!Ground State
       c=1

       do i=1, numberOfSpecies

          print *, "n particles", MolecularSystem_getNumberOfParticles(i)
          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          !!Singles
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
             do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )

                   if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                     c=c+1
                   end if
                end do
             end do
          end if


       end do
       !print *, "number OF confi", c 


       do i=1, numberOfSpecies
   !!Doubles of the same specie

          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
             do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                do n=m+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   do p= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                      !print *, "dou", m,n,p
                      !if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                        do q=p+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                   !do q=ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)+1, ConfigurationInteraction_instance%numberOfOrbitals%values(i) 
                      !print *, "dou", m,n,p
                      if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) .and. (p /= q) ) then !! alpha -> alpha, beta -> beta

                             c=c+1
                           !end if
                      end if
                        end do

                   end do
                end do
             end do
          end if

          !!Doubles of different species
          do j=i+1, numberOfSpecies
             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
                do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                   do n=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                      do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                         do q= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                            c=c+1
                         end do
                      end do
                   end do
                end do
             end if
          end do

   !!triples of diff specie
          do j=i+1, numberOfSpecies
             do k=j+1, numberOfSpecies

             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. & 
                ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 .and. &
                ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 ) then
                do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                   do n=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                      do r=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) )
                      do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                         do q= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                         do s= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(k))
                            c=c+1
                         end do
                      end do
                         end do
                      end do
                   end do
                end do
             end if
                         end do
          end do


       end do !! species

!      c = 6
       print *, "number OF doubles confi", c 

       ConfigurationInteraction_instance%numberOfConfigurations = c
       allocate (ConfigurationInteraction_instance%configurations(ConfigurationInteraction_instance%numberOfConfigurations) )

       call Matrix_Constructor(ConfigurationInteraction_instance%hamiltonianMatrix, int(ConfigurationInteraction_instance%numberOfConfigurations,8),int(ConfigurationInteraction_instance%numberOfConfigurations,8),0.0_8)



    case ( "FCI-oneSpecie" )

    case default

       call ConfigurationInteraction_exception( ERROR, "Configuration interactor constructor", "Correction level not implemented")

    end select

    close(wfnUnit)

  end subroutine ConfigurationInteraction_constructor


  !>
  !! @brief Destructor por omision
  !!
  !! @param this
  !<
  subroutine ConfigurationInteraction_destructor()
    implicit none
    integer i,j,m,n,p,q,c
    integer numberOfSpecies
    integer :: isLambdaEqual1

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    select case ( trim(ConfigurationInteraction_instance%level) )

    case ( "CIS" )

       !!Destroy configurations
       !!Ground State
       c=1
       call Configuration_destructor(ConfigurationInteraction_instance%configurations(c) )

       !do i=1, numberOfSpecies
          !!Singles
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
             do c=2, ConfigurationInteraction_instance%numberOfConfigurations
               call Configuration_destructor(ConfigurationInteraction_instance%configurations(c) )                
             end do
          end if
       !end do

    deallocate(ConfigurationInteraction_instance%configurations)

    call Matrix_destructor(ConfigurationInteraction_instance%hamiltonianMatrix)
    call Vector_destructor (ConfigurationInteraction_instance%numberOfOccupiedOrbitals)
    call Vector_destructor (ConfigurationInteraction_instance%numberOfOrbitals)
    call Vector_destructor (ConfigurationInteraction_instance%lambda)


    case ( "CISD" )

     !!Destroy configurations
       !!Ground State
       c=1
       call Configuration_destructor(ConfigurationInteraction_instance%configurations(c) )

       !do i=1, numberOfSpecies
          !!Singles
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
             do c=2, ConfigurationInteraction_instance%numberOfConfigurations
               call Configuration_destructor(ConfigurationInteraction_instance%configurations(c) )                
             end do
          end if
       !end do

    deallocate(ConfigurationInteraction_instance%configurations)

    call Matrix_destructor(ConfigurationInteraction_instance%hamiltonianMatrix)
    call Vector_destructor (ConfigurationInteraction_instance%numberOfOccupiedOrbitals)
    call Vector_destructor (ConfigurationInteraction_instance%numberOfOrbitals)
    call Vector_destructor (ConfigurationInteraction_instance%lambda)

    case ( "CIDD" )

       !!Count
     !!Destroy configurations
       !!Ground State
       c=1
       call Configuration_destructor(ConfigurationInteraction_instance%configurations(c) )

       !do i=1, numberOfSpecies
          !!Doubles of the same specie
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
             do c=2, ConfigurationInteraction_instance%numberOfConfigurations
               call Configuration_destructor(ConfigurationInteraction_instance%configurations(c) )                
             end do
          end if
       !end do



       deallocate (ConfigurationInteraction_instance%configurations )
       
       call Matrix_destructor(ConfigurationInteraction_instance%hamiltonianMatrix)
       call Vector_destructor (ConfigurationInteraction_instance%numberOfOccupiedOrbitals)
       call Vector_destructor (ConfigurationInteraction_instance%numberOfOrbitals)
       call Vector_destructor (ConfigurationInteraction_instance%lambda)

    case ( "FCI-oneSpecie" )

    case ( "FCI" )

     !!Destroy configurations
       !!Ground State
       c=1
       call Configuration_destructor(ConfigurationInteraction_instance%configurations(c) )

       !do i=1, numberOfSpecies
          !!Singles
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
             do c=2, ConfigurationInteraction_instance%numberOfConfigurations
               call Configuration_destructor(ConfigurationInteraction_instance%configurations(c) )                
             end do
          end if
       !end do

    deallocate(ConfigurationInteraction_instance%configurations)

    call Matrix_destructor(ConfigurationInteraction_instance%hamiltonianMatrix)
    call Vector_destructor (ConfigurationInteraction_instance%numberOfOccupiedOrbitals)
    call Vector_destructor (ConfigurationInteraction_instance%numberOfOrbitals)
    call Vector_destructor (ConfigurationInteraction_instance%lambda)

    case default

    end select


    ConfigurationInteraction_instance%isInstanced=.false.

  end subroutine ConfigurationInteraction_destructor

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_show()
    implicit none
    type(ConfigurationInteraction) :: this
    integer :: i
    integer :: m
    real(8) :: davidsonCorrection, HFcoefficient, CIcorrection
 
    if ( ConfigurationInteraction_instance%isInstanced ) then

       print *,""
       print *," POST HARTREE-FOCK CALCULATION"
       print *," CONFIGURATION INTERACTION THEORY:"
       print *,"=============================="
       print *,""
       write (6,"(T8,A30, A5)") "LEVEL = ", ConfigurationInteraction_instance%level
       write (6,"(T8,A30, I8)") "NUMBER OF CONFIGURATIONS = ", ConfigurationInteraction_instance%numberOfConfigurations
       do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
        write (6,"(T8,A17,I3,A10, F18.12)") "STATE: ", i, " ENERGY = ", ConfigurationInteraction_instance%eigenvalues%values(i)
       end do
       print *, ""
       CIcorrection = ConfigurationInteraction_instance%eigenvalues%values(1) - &
                HartreeFock_instance%totalEnergy

       write (6,"(T4,A34, F20.12)") "GROUND STATE CORRELATION ENERGY = ", CIcorrection

       if (  ConfigurationInteraction_instance%level == "CISD" ) then
         print *, ""
         write (6,"(T2,A34)") "RENORMALIZED DAVIDSON CORRECTION:"
         print *, ""
         write (6,"(T8,A54)") "E(CISDTQ) \approx E(CISD) + \delta E(Q)               "
         write (6,"(T8,A54)") "\delta E(Q) = (1 - c_0^2) * \delta E(CISD) / c_0^2    "
         print *, ""
         HFcoefficient = ConfigurationInteraction_instance%eigenVectors%values(1,1) 
         davidsonCorrection = ( 1 - HFcoefficient*HFcoefficient) * CIcorrection / (HFcoefficient*HFcoefficient)
  
  
         write (6,"(T8,A19, F20.12)") "HF COEFFICIENT = ", HFcoefficient
         write (6,"(T8,A19, F20.12)") "\delta E(Q) = ", davidsonCorrection
         write (6,"(T8,A19, F20.12)") "E(CISDTQ) ESTIMATE ",  HartreeFock_instance%totalEnergy +&
            CIcorrection + davidsonCorrection
       end if

       ! do i=1, ConfigurationInteraction_instance%numberOfConfigurations
       !    call Configuration_show (ConfigurationInteraction_instance%configurations(i))
       ! end do

    else 

    end if

  end subroutine ConfigurationInteraction_show

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_run()
    implicit none 
    integer :: m
    real(8), allocatable :: eigenValues(:) 

    select case ( trim(ConfigurationInteraction_instance%level) )

    case ( "CIS" )

       print *, ""
       print *, "==============================================="
       print *, "|            BEGIN CIS CALCULATION            |"
       print *, "-----------------------------------------------"
       print *, ""

       call ConfigurationInteraction_getTransformedIntegrals()
       print *, "Building configurations..."

       call ConfigurationInteraction_buildConfigurations()
       print *, "Total number of configurations", ConfigurationInteraction_instance%numberOfConfigurations
       print *, ""
       print *, "Building hamiltonian..."
       call ConfigurationInteraction_buildHamiltonianMatrix()

       call Vector_constructor ( ConfigurationInteraction_instance%eigenvalues, &
                                 ConfigurationInteraction_instance%numberOfConfigurations, 0.0_8)
       print *, "Reference Energy", ConfigurationInteraction_instance%hamiltonianMatrix%values(1,1)

       print *, ""
       print *, "Diagonalizing hamiltonian..."

       select case (CONTROL_instance%CI_DIAGONALIZATION_METHOD) 

       case ("DSYEVX")

         call Matrix_eigen_select (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
              1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
              flags = SYMMETRIC, dm = ConfigurationInteraction_instance%numberOfConfigurations )

       case ("DSYEVR")

         call Matrix_eigen_dsyevr (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
              1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
              flags = SYMMETRIC, dm = ConfigurationInteraction_instance%numberOfConfigurations )

        end select

       print *,""
       print *, "-----------------------------------------------"
       print *, "|              END CIS CALCULATION            |"
       print *, "==============================================="
       print *, ""

    case ( "CISD" )

       print *, ""
       print *, "==============================================="
       print *, "|            BEGIN CISD CALCULATION           |"
       print *, "-----------------------------------------------"
       print *, ""

       call ConfigurationInteraction_getTransformedIntegrals()
       print *, "Building configurations..."

       call ConfigurationInteraction_buildConfigurations()
       print *, "Total number of configurations", ConfigurationInteraction_instance%numberOfConfigurations
       print *, ""
       print *, "Building hamiltonian..."
       call ConfigurationInteraction_buildHamiltonianMatrix()

       call Vector_constructor ( ConfigurationInteraction_instance%eigenvalues, &
                                 ConfigurationInteraction_instance%numberOfConfigurations, 0.0_8)
       print *, "Reference Energy", ConfigurationInteraction_instance%hamiltonianMatrix%values(1,1)

       print *, ""
       print *, "Diagonalizing hamiltonian..."
       print *, "  Using : ", (CONTROL_instance%CI_DIAGONALIZATION_METHOD)

       select case (CONTROL_instance%CI_DIAGONALIZATION_METHOD) 

       case ("DSYEVX")

         call Matrix_constructor (ConfigurationInteraction_instance%eigenVectors, &
              int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
              int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

         call Matrix_eigen_select (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
              1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
              eigenVectors = ConfigurationInteraction_instance%eigenVectors, &
              flags = SYMMETRIC)

!         call Matrix_eigen_select (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
!              1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
!              flags = SYMMETRIC, dm = ConfigurationInteraction_instance%numberOfConfigurations )


       case ("DSYEVR")

         call Matrix_constructor (ConfigurationInteraction_instance%eigenVectors, &
              int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
              int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

         call Matrix_eigen_dsyevr (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
              1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
              eigenVectors = ConfigurationInteraction_instance%eigenVectors, &
              flags = SYMMETRIC)

!        call Matrix_eigen_dsyevr (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
!              1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
!              flags = SYMMETRIC, dm = ConfigurationInteraction_instance%numberOfConfigurations )

        end select

!! Debugging
!!     call Matrix_constructor (ConfigurationInteraction_instance%eigenVectors, &
!!                                 int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
!!                                 int(ConfigurationInteraction_instance%numberOfConfigurations,8), 0.0_8)

!!     call Matrix_eigen (ConfigurationInteraction_instance%hamiltonianMatrix, &
!!                          ConfigurationInteraction_instance%eigenvalues, &
!!                          flags = SYMMETRIC, &
!!                          eigenVectors = ConfigurationInteraction_instance%eigenVectors)

!!      print *, "eigenvalues"
!!      call Vector_show (ConfigurationInteraction_instance%eigenValues)


!!      print *, "hamiltonianMatrix"
!!      call Matrix_show (ConfigurationInteraction_instance%hamiltonianMatrix)
!!      call Vector_show( ConfigurationInteraction_instance%eigenvalues)
       
       print *,""
       print *, "-----------------------------------------------"
       print *, "|              END CISD CALCULATION           |"
       print *, "==============================================="
       print *, ""

    case ( "CIDD" )

       print *, ""
       print *, ""
       print *, "==============================================="
       print *, "|            BEGIN CID CALCULATION            |"
       print *, "|       ONLY EXCITATIONS TO THE SAME ORBITAL  |"
       print *, "-----------------------------------------------"
       print *, ""

       call ConfigurationInteraction_buildConfigurations()

       call ConfigurationInteraction_getTransformedIntegrals()

       call ConfigurationInteraction_buildHamiltonianMatrix()

       call Vector_constructor(ConfigurationInteraction_instance%eigenvalues, ConfigurationInteraction_instance%numberOfConfigurations, 0.0_8)
        

       call Matrix_constructor (ConfigurationInteraction_instance%eigenVectors, &
                                 int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
                                 int(ConfigurationInteraction_instance%numberOfConfigurations,8), 0.0_8)


       call Matrix_eigen (ConfigurationInteraction_instance%hamiltonianMatrix, &
                          ConfigurationInteraction_instance%eigenvalues, &
                          flags = SYMMETRIC, &
                          eigenVectors = ConfigurationInteraction_instance%eigenVectors)

       print *,""
       print *, "-----------------------------------------------"
       print *, "|              END CID CALCULATION            |"
       print *, "==============================================="
       print *, ""

    case ( "FCI-oneSpecie" )

       print *, ""
       print *, ""
       print *, "==============================================="
       print *, "|  Full CI for one specie calculation          |"
       print *, "|  Use fci program to perform the calculation  |"
       print *, "-----------------------------------------------"
       print *, ""
       ! call ConfigurationInteraction_getTransformedIntegrals()
       !call ConfigurationInteraction_printTransformedIntegralsToFile()


    case ( "FCI" )

       print *, ""
       print *, "==============================================="
       print *, "|            BEGIN FCI CALCULATION            |"
       print *, "-----------------------------------------------"
       print *, ""

       call ConfigurationInteraction_getTransformedIntegrals()
       print *, "Building configurations..."

       call ConfigurationInteraction_buildConfigurations()
       print *, "Total number of configurations", ConfigurationInteraction_instance%numberOfConfigurations
       print *, ""
       print *, "Building hamiltonian..."
       call ConfigurationInteraction_buildHamiltonianMatrix()

       call Vector_constructor ( ConfigurationInteraction_instance%eigenvalues, &
                                 ConfigurationInteraction_instance%numberOfConfigurations, 0.0_8)
       print *, "Reference Energy", ConfigurationInteraction_instance%hamiltonianMatrix%values(1,1)

       print *, ""
       print *, "Diagonalizing hamiltonian..."
       print *, "  Using : ", (CONTROL_instance%CI_DIAGONALIZATION_METHOD)

       select case (CONTROL_instance%CI_DIAGONALIZATION_METHOD) 

       case ("DSYEVX")

         call Matrix_constructor (ConfigurationInteraction_instance%eigenVectors, &
              int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
              int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

         call Matrix_eigen_select (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
              1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
              eigenVectors = ConfigurationInteraction_instance%eigenVectors, &
              flags = SYMMETRIC)

!         call Matrix_eigen_select (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
!              1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
!              flags = SYMMETRIC, dm = ConfigurationInteraction_instance%numberOfConfigurations )


       case ("DSYEVR")

         call Matrix_constructor (ConfigurationInteraction_instance%eigenVectors, &
              int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
              int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

         call Matrix_eigen_dsyevr (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
              1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
              eigenVectors = ConfigurationInteraction_instance%eigenVectors, &
              flags = SYMMETRIC)

!        call Matrix_eigen_dsyevr (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
!              1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
!              flags = SYMMETRIC, dm = ConfigurationInteraction_instance%numberOfConfigurations )

        end select

       print *,""
       print *, "-----------------------------------------------"
       print *, "|              END FCI CALCULATION            |"
       print *, "==============================================="
       print *, ""

    case default

       call ConfigurationInteraction_exception( ERROR, "Configuration interactor constructor", "Correction level not implemented")

    end select


  end subroutine ConfigurationInteraction_run


  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_buildConfigurations()
    implicit none

    integer :: numberOfSpecies
    integer :: i,j,k,m,n,p,q,a,b,c,d,cc,r,s
    integer :: isLambdaEqual1
    type(vector) :: order
    type(vector) :: occupiedCode
    type(vector) :: unoccupiedCode
    logical :: sameConfiguration
    real(8) :: CIenergy(2,2)
    integer :: nEquivalentConfigurations, newNumberOfConfigurations
    integer, allocatable :: equivalentConfigurations (:,:), auxArray(:,:), auxvector(:),auxvectorA(:)

    integer :: lambda

    nEquivalentConfigurations = 0
    if (allocated ( equivalentConfigurations )) deallocate ( equivalentConfigurations)
    allocate( equivalentConfigurations(nEquivalentConfigurations,2) )
    equivalentConfigurations = 0

    newNumberOfConfigurations = ConfigurationInteraction_instance%numberOfConfigurations 

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    select case ( trim(ConfigurationInteraction_instance%level) )

    case ( "CIS" )

   !!Build configurations
       !!Ground State
       c=1
       call Vector_constructor (order, numberOfSpecies, 0.0_8)
       call Vector_constructor (occupiedCode, numberOfSpecies, 0.0_8)
       call Vector_constructor (unoccupiedCode, numberOfSpecies, 0.0_8)
       call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations,c) 

       do i=1, numberOfSpecies
          cc = 0
          call Vector_constructor (order, numberOfSpecies, 0.0_8)
          order%values(i)=1
          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital

          !!Singles
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
             do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                call Vector_constructor (occupiedCode, numberOfSpecies, 0.0_8)
                occupiedCode%values(i)=m
                do p=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i))
                   if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                     call Vector_constructor (unoccupiedCode, numberOfSpecies, 0.0_8)
                     unoccupiedCode%values(i)=p
                     c=c+1
                     cc = cc + 1
                     call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                          ConfigurationInteraction_instance%numberOfConfigurations,c) 

                   end if

                end do !! p
             end do !! m
          end if
          write (6, "(T4,A20,A21,A3,I4)") "Singles for species ", trim(MolecularSystem_getNameOfSpecie(i)), " : ", cc
        end do !! species

        !do c = 1,  ConfigurationInteraction_instance%numberOfConfigurations 
        !   call Configuration_show(  ConfigurationInteraction_instance%configurations(c) )
        !end do 

       !! Search for equivalent configurations
       do a=1, ConfigurationInteraction_instance%numberOfConfigurations
          do b=a+1, ConfigurationInteraction_instance%numberOfConfigurations

             call Configuration_checkTwoConfigurations(ConfigurationInteraction_instance%configurations(a), &
                    ConfigurationInteraction_instance%configurations(b), sameConfiguration, numberOfSpecies)

             if ( sameConfiguration .eqv. .True. ) then
               !! append one value....
               nEquivalentConfigurations = nEquivalentConfigurations + 1
               allocate(auxArray(nEquivalentConfigurations,2))
               auxArray = 0
               auxArray(:nEquivalentConfigurations-1,:) = equivalentConfigurations(:nEquivalentConfigurations-1,:)
               deallocate(equivalentConfigurations)
               allocate(equivalentConfigurations(nEquivalentConfigurations,2))
               equivalentConfigurations(:,:) = auxArray(:,:) 
               deallocate(auxArray)

               equivalentConfigurations(nEquivalentConfigurations,1) = a
               equivalentConfigurations(nEquivalentConfigurations,2) = b

               newNumberOfConfigurations = newNumberOfConfigurations -1

             end if
             sameConfiguration = .false.
          end do
       end do


       !! Remove equivalent configurations
        do c = nEquivalentConfigurations, 1, -1

                call Configuration_destructor(ConfigurationInteraction_instance%configurations( &
                     equivalentConfigurations(c,2) ) )                
        end do

        !! Compact the configuration vector... 
        do c = 1,  ConfigurationInteraction_instance%numberOfConfigurations 

          if ( c <= newNumberOfConfigurations ) then

            if ( ConfigurationInteraction_instance%configurations(c)%isInstanced .eqv. .false. ) then
              do d = c +1 , ConfigurationInteraction_instance%numberOfConfigurations 

                if ( ConfigurationInteraction_instance%configurations(d)%isInstanced .eqv. .true. ) then
                  call Configuration_copyConstructor(  ConfigurationInteraction_instance%configurations(d), &
                    ConfigurationInteraction_instance%configurations(c) )

                  call Configuration_destructor(ConfigurationInteraction_instance%configurations(d))                

                  exit
                end if

              end do 
            end if

         end if

        end do 

        ConfigurationInteraction_instance%numberOfConfigurations = newNumberOfConfigurations 

        !do c = 1,  ConfigurationInteraction_instance%numberOfConfigurations 
        !   call Configuration_show(  ConfigurationInteraction_instance%configurations(c) )
        !end do 
         
       call Matrix_destructor(ConfigurationInteraction_instance%hamiltonianMatrix)

       !! Rebuild the hamiltonian matrix without the equivalent configurations
       call Matrix_Constructor(ConfigurationInteraction_instance%hamiltonianMatrix, &
            int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
            int(ConfigurationInteraction_instance%numberOfConfigurations,8),0.0_8 )

    case ( "CISD" )

   !!Build configurations
       !!Ground State
       c=1
       call Vector_constructor (order, numberOfSpecies, 0.0_8)
       call Vector_constructor (occupiedCode, numberOfSpecies, 0.0_8)
       call Vector_constructor (unoccupiedCode, numberOfSpecies, 0.0_8)
       call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations, c) 

       do i=1, numberOfSpecies

          call Vector_constructor (order, numberOfSpecies, 0.0_8)
          order%values(i)=1
          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          cc = 0
          !!Singles
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
             do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                call Vector_constructor (occupiedCode, numberOfSpecies, 0.0_8)
                occupiedCode%values(i)=m
                do p=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1,int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                   if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                     call Vector_constructor (unoccupiedCode, numberOfSpecies, 0.0_8)
                     unoccupiedCode%values(i)=p
                     c=c+1
                     cc=cc+1
                     call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                          ConfigurationInteraction_instance%numberOfConfigurations,c) 

                   end if

                end do !! p
             end do !! m
          end if

          write (6, "(T4,A20,A21,A3,I4)") "Singles for species ", trim(MolecularSystem_getNameOfSpecie(i)), " : ", cc

        end do !! species

       do i=1, numberOfSpecies
          !!Doubles of the same specie
          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          call Vector_constructor (order, numberOfSpecies, 0.0_8)
          order%values(i)=2
          cc = 0
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
             do m=1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                do n=m+1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   call Vector_constructor (occupiedCode, numberOfSpecies, 0.0_8)
                   occupiedCode%values(i)=m*1024+n
                   do p = int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                     do q=p+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                   !do q=ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)+1, ConfigurationInteraction_instance%numberOfOrbitals%values(i) 

                      if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) .and. (p /= q) ) then !! alpha -> alpha, beta -> beta
                            call Vector_constructor (unoccupiedCode, numberOfSpecies, 0.0_8)
                            unoccupiedCode%values(i)=p*1024+q
                            c=c+1
                            cc=cc+1
                            call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), &
                                  occupiedCode, unoccupiedCode, order, &
                                   ConfigurationInteraction_instance%numberOfConfigurations, c) 
                     end if
                       end do
                   end do
                end do
             end do
          end if
          write (6, "(T4,A20,A21,A3,I4)") "Doubles for species ", trim(MolecularSystem_getNameOfSpecie(i)), " : ", cc

          !!Doubles of different species
          do j=i+1, numberOfSpecies

             call Vector_constructor (order, numberOfSpecies, 0.0_8)
             order%values(i)=1
             order%values(j)=1
             cc = 0

             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
                do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   do n=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))

                      call Vector_constructor (occupiedCode, numberOfSpecies, 0.0_8)
                      occupiedCode%values(i)=m
                      occupiedCode%values(j)=n
                      do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i))
                         do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                            call Vector_constructor (unoccupiedCode, numberOfSpecies, 0.0_8)
                            unoccupiedCode%values(i)=p
                            unoccupiedCode%values(j)=q
                            c=c+1
                            cc=cc+1
                            call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations, c) 
                         end do
                      end do
                   end do
                end do
             end if

            write (6, "(T4,A20,A21,A3,I4)") "Doubles for species ", trim(MolecularSystem_getNameOfSpecie(i))//&
                                                   "/"//trim(MolecularSystem_getNameOfSpecie(j)), " : ", cc
          end do

       end do
!!!!

       !! Search for equivalent configurations
!!       do a=1, ConfigurationInteraction_instance%numberOfConfigurations
!!          do b=a+1, ConfigurationInteraction_instance%numberOfConfigurations
!!             call Configuration_checkTwoConfigurations(ConfigurationInteraction_instance%configurations(a), &
!!                    ConfigurationInteraction_instance%configurations(b), sameConfiguration, numberOfSpecies)
!!
!!             if ( sameConfiguration .eqv. .True. ) then
!!
!!               !! append one value....
!!               nEquivalentConfigurations = nEquivalentConfigurations + 1
!!               allocate(auxArray(nEquivalentConfigurations,2))
!!               auxArray = 0
!!               auxArray(:nEquivalentConfigurations-1,:) = equivalentConfigurations(:nEquivalentConfigurations-1,:)
!!               deallocate(equivalentConfigurations)
!!               allocate(equivalentConfigurations(nEquivalentConfigurations,2))
!!               equivalentConfigurations(:,:) = auxArray(:,:) 
!!               deallocate(auxArray)
!!
!!               equivalentConfigurations(nEquivalentConfigurations,1) = a
!!               equivalentConfigurations(nEquivalentConfigurations,2) = b
!!
!!               newNumberOfConfigurations = newNumberOfConfigurations -1
!!
!!             end if
!!             sameConfiguration = .false.
!!          end do
!!       end do
!!
!!
!!       print *, "nEquivalentConfigurations",  nEquivalentConfigurations
!!       allocate(auxvector( ConfigurationInteraction_instance%numberOfConfigurations ))
!!       allocate(auxvectorA( ConfigurationInteraction_instance%numberOfConfigurations ))
!!       auxvector = 0
!!       auxvectorA = 0
!!       !! Remove equivalent configurations
!!        do c = nEquivalentConfigurations, 1, -1
!!                print *, c,  equivalentConfigurations(c,1),  equivalentConfigurations(c,2)
!!                auxvector(equivalentConfigurations(c,2)) = auxvector(equivalentConfigurations(c,2)) + 1
!!!                auxvector(equivalentConfigurations(c,1)) = auxvector(equivalentConfigurations(c,1)) + 1
!!
!!                print *,auxvector(equivalentConfigurations(c,2))
!!                if ( auxvector(equivalentConfigurations(c,2)) < 2 ) then
!!                  call Configuration_destructor(ConfigurationInteraction_instance%configurations( &
!!                     equivalentConfigurations(c,2) ) )                
!!                end if
!!        end do
!!       print *, " newNumberOfConfiguration" ,newNumberOfConfigurations
!!        !! Compact the configuration vector... 
!!        do c = 1,  ConfigurationInteraction_instance%numberOfConfigurations 
!!
!!!          if ( c <= newNumberOfConfigurations ) then
!!
!!            if ( ConfigurationInteraction_instance%configurations(c)%isInstanced .eqv. .false. ) then
!!              do d = c +1 , ConfigurationInteraction_instance%numberOfConfigurations 
!!
!!                if ( ConfigurationInteraction_instance%configurations(d)%isInstanced .eqv. .true. ) then
!!                  call Configuration_copyConstructor(  ConfigurationInteraction_instance%configurations(d), &
!!                    ConfigurationInteraction_instance%configurations(c) )
!!
!!                  call Configuration_destructor(ConfigurationInteraction_instance%configurations(d))                
!!
!!                  exit
!!                end if
!!
!!              end do 
!!            end if
!!
!!!          end if
!!
!!        end do 
!!
!!        newNumberOfConfigurations  = 0
!!
!!        do c = 1,  ConfigurationInteraction_instance%numberOfConfigurations 
!!                if ( ConfigurationInteraction_instance%configurations(c)%isInstanced .eqv. .true. ) then
!!                  newNumberOfConfigurations = newNumberOfConfigurations + 1
!!                end if
!!        end do 
!!
!!        ConfigurationInteraction_instance%numberOfConfigurations = newNumberOfConfigurations 
!!
!!        do c = 1,  ConfigurationInteraction_instance%numberOfConfigurations 
!!           call Configuration_show(  ConfigurationInteraction_instance%configurations(c) )
!!        end do 
!!
!!         
!!       call Matrix_destructor(ConfigurationInteraction_instance%hamiltonianMatrix)
!!
!!       !! Rebuild the hamiltonian matrix without the equivalent configurations
!!       call Matrix_Constructor(ConfigurationInteraction_instance%hamiltonianMatrix, &
!!            int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
!!            int(ConfigurationInteraction_instance%numberOfConfigurations,8),0.0_8 )

!       ConfigurationInteraction_instance%configurations(2)%occupations(1)%values = (/1,0,0,1/)


    case ( "CIDD" ) 

       !!Build configurations
       !!Ground State
       c=1
       call Vector_constructor (order, numberOfSpecies, 0.0_8)
       call Vector_constructor (occupiedCode, numberOfSpecies, 0.0_8)
       call Vector_constructor (unoccupiedCode, numberOfSpecies, 0.0_8)
       call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations, c) 

       do i=1, numberOfSpecies
          !!Doubles of the same specie
          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          call Vector_constructor (order, numberOfSpecies, 0.0_8)
          order%values(i)=2
          cc = 0
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
             do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                do n=m+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   call Vector_constructor (occupiedCode, numberOfSpecies, 0.0_8)
                   occupiedCode%values(i)=m*1024+n
                   do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int( ConfigurationInteraction_instance%numberOfOrbitals%values(i))
                        do q=p+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )

                      if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) .and. (p /= q) ) then !! alpha -> alpha, beta -> beta
                            call Vector_constructor (unoccupiedCode, numberOfSpecies, 0.0_8)
                            unoccupiedCode%values(i)=p*1024+q
                            c=c+1
                            cc=cc+1
                            call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), &
                                  occupiedCode, unoccupiedCode, order, &
                                   ConfigurationInteraction_instance%numberOfConfigurations, c) 
  
                     end if
                       end do
                   end do
                end do
             end do
          end if
          write (6, "(T4,A20,A21,A3,I4)") "Doubles for species ", trim(MolecularSystem_getNameOfSpecie(i)), " : ", cc

       end do

       call Vector_Destructor(order)
       call Vector_Destructor(occupiedCode)
       call Vector_Destructor(unoccupiedCode)

    case ( "FCI" )

   !!Build configurations
       !!Ground State
       c=1
       call Vector_constructor (order, numberOfSpecies, 0.0_8)
       call Vector_constructor (occupiedCode, numberOfSpecies, 0.0_8)
       call Vector_constructor (unoccupiedCode, numberOfSpecies, 0.0_8)
       call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations, c) 

       do i=1, numberOfSpecies

          call Vector_constructor (order, numberOfSpecies, 0.0_8)
          order%values(i)=1
          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          cc = 0
          !!Singles
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
             do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                call Vector_constructor (occupiedCode, numberOfSpecies, 0.0_8)
                occupiedCode%values(i)=m
                do p=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1,int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                   if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                     call Vector_constructor (unoccupiedCode, numberOfSpecies, 0.0_8)
                     unoccupiedCode%values(i)=p
                     c=c+1
                     cc=cc+1
                     call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                          ConfigurationInteraction_instance%numberOfConfigurations,c) 

                   end if

                end do !! p
             end do !! m
          end if

          write (6, "(T4,A20,A21,A3,I4)") "Singles for species ", trim(MolecularSystem_getNameOfSpecie(i)), " : ", cc

        end do !! species

       do i=1, numberOfSpecies
          !!Doubles of the same specie
          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          call Vector_constructor (order, numberOfSpecies, 0.0_8)
          order%values(i)=2
          cc = 0
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
             do m=1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                do n=m+1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   call Vector_constructor (occupiedCode, numberOfSpecies, 0.0_8)
                   occupiedCode%values(i)=m*1024+n
                   do p = int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                     do q=p+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                   !do q=ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)+1, ConfigurationInteraction_instance%numberOfOrbitals%values(i) 

                      if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) .and. (p /= q) ) then !! alpha -> alpha, beta -> beta
                            call Vector_constructor (unoccupiedCode, numberOfSpecies, 0.0_8)
                            unoccupiedCode%values(i)=p*1024+q
                            c=c+1
                            cc=cc+1
                            call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), &
                                  occupiedCode, unoccupiedCode, order, &
                                   ConfigurationInteraction_instance%numberOfConfigurations, c) 
                     end if
                       end do
                   end do
                end do
             end do
          end if
          write (6, "(T4,A20,A21,A3,I4)") "Doubles for species ", trim(MolecularSystem_getNameOfSpecie(i)), " : ", cc

          !!Doubles of different species
          do j=i+1, numberOfSpecies

             call Vector_constructor (order, numberOfSpecies, 0.0_8)
             order%values(i)=1
             order%values(j)=1
             cc = 0

             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
                do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   do n=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))

                      call Vector_constructor (occupiedCode, numberOfSpecies, 0.0_8)
                      occupiedCode%values(i)=m
                      occupiedCode%values(j)=n
                      do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i))
                         do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                            call Vector_constructor (unoccupiedCode, numberOfSpecies, 0.0_8)
                            unoccupiedCode%values(i)=p
                            unoccupiedCode%values(j)=q
                            c=c+1
                            cc=cc+1
                            call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations, c) 
                         end do
                      end do
                   end do
                end do
             end if

            write (6, "(T4,A20,A21,A3,I4)") "Doubles for species ", trim(MolecularSystem_getNameOfSpecie(i))//&
                                                   "/"//trim(MolecularSystem_getNameOfSpecie(j)), " : ", cc
          end do

          !!triples of different species
          do j=i+1, numberOfSpecies

            do k=j+1, numberOfSpecies

             call Vector_constructor (order, numberOfSpecies, 0.0_8)
             order%values(i)=1
             order%values(j)=1
             order%values(k)=1
             cc = 0

             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. &
                 ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 .and. &
                 ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 ) then
                do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   do n=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))
                     do r=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))

                      call Vector_constructor (occupiedCode, numberOfSpecies, 0.0_8)
                      occupiedCode%values(i)=m
                      occupiedCode%values(j)=n
                      occupiedCode%values(k)=r

                      do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i))
                         do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )

                         do s=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(k) )
                            call Vector_constructor (unoccupiedCode, numberOfSpecies, 0.0_8)
                            unoccupiedCode%values(i)=p
                            unoccupiedCode%values(j)=q
                            unoccupiedCode%values(k)=s
                            c=c+1
                            cc=cc+1
                            call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations, c) 

                         end do
                         end do
                      end do
                      end do
                   end do
                end do

             end if

            write (6, "(T4,A20,A21,A3,I4)") "Triples for species ", trim(MolecularSystem_getNameOfSpecie(i))//&
                                                   "/"//trim(MolecularSystem_getNameOfSpecie(j))//& 
                                                   "/"//trim(MolecularSystem_getNameOfSpecie(k)), " : ", cc
            end do
          end do


       end do

    case default

       call ConfigurationInteraction_exception( ERROR, "Configuration interactor constructor", "Correction level not implemented")

    end select

  end subroutine ConfigurationInteraction_buildConfigurations

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_buildHamiltonianMatrix()
    implicit none

!    type(TransformIntegrals) :: repulsionTransformer
    integer :: numberOfSpecies
    integer :: i,j,a,b,pos,neg
    integer :: m,n,l,k,z
    integer :: specieID
    integer :: otherSpecieID
    character(10) :: nameOfSpecie
    character(10) :: nameOfOtherSpecie
    integer :: numberOfOrbitals
    integer :: numberOfSpatialOrbitals
    integer :: numberOfOtherSpecieOrbitals
    integer :: numberOfOtherSpecieSpatialOrbitals
    integer, allocatable :: spin(:)
    integer, allocatable :: spatialOrbital(:)
    integer :: lambda !occupation per orbital
    integer :: lambdaOfOtherSpecie !occupation per orbital
    real(8) :: kappa !positive or negative exchange
    real(8) :: charge
    real(8) :: otherSpecieCharge
    real(8) :: twoParticlesEnergy
    real(8) :: couplingEnergy
    integer(8) :: auxIndex
    type(vector) :: diffAB
    type(vector), allocatable :: differentOrbitals(:)
    
    real(8) :: CIenergy

!!    n = 2
!!    auxIndex = IndexMap_tensorR4ToVector( 1, 1, 1, 1, n )
!!    print *, "1 1 | 1 1", ConfigurationInteraction_instance%fourCenterIntegrals(1,1)%values(auxIndex, 1)


    !a,b configuration iterators
    !i,j specie iterators
    !k,l orbital iterators


    ! print *, "build H"
       do a=1, ConfigurationInteraction_instance%numberOfConfigurations
          do b=a, ConfigurationInteraction_instance%numberOfConfigurations

            !call ConfigurationInteraction_calculateCIenergyB(a,b,CIenergy)
            call ConfigurationInteraction_calculateCIenergy(a,b,CIenergy)
            ConfigurationInteraction_instance%hamiltonianMatrix%values(a,b) = CIenergy
               
          end do
       end do
!     print *, "diag hamiltonian" 
         do a=1, ConfigurationInteraction_instance%numberOfConfigurations
            do b=a, ConfigurationInteraction_instance%numberOfConfigurations
!              if ( a == b ) print *, ConfigurationInteraction_instance%hamiltonianMatrix%values(a,b)
               ConfigurationInteraction_instance%hamiltonianMatrix%values(b,a)=ConfigurationInteraction_instance%hamiltonianMatrix%values(a,b)
            end do
         end do

   deallocate(ConfigurationInteraction_instance%twoCenterIntegrals)
   deallocate(ConfigurationInteraction_instance%fourCenterIntegrals)

  end subroutine ConfigurationInteraction_buildHamiltonianMatrix

  subroutine ConfigurationInteraction_calculateCIenergy(a,b,CIenergy)
    implicit none

!    type(TransformIntegrals) :: repulsionTransformer
    integer :: numberOfSpecies
    integer :: i,j,a,b,pos,neg,ia,ib
    integer :: m,n,l,k,z
    integer :: specieID
    integer :: otherSpecieID
    character(10) :: nameOfSpecie
    character(10) :: nameOfOtherSpecie
    integer :: numberOfOrbitals
    integer :: numberOfSpatialOrbitals
    integer :: numberOfOtherSpecieOrbitals
    integer :: numberOfOtherSpecieSpatialOrbitals
    integer, allocatable :: spin(:)
    integer, allocatable :: spatialOrbital(:)
    integer :: auxSpatialOrbital
    integer :: lambda !occupation per orbital
    integer :: lambdaOfOtherSpecie !occupation per orbital
    integer :: numberOfOccupiedOrbitals
    real(8) :: kappa !positive or negative exchange
    real(8) :: charge
    real(8) :: otherSpecieCharge
    real(8) :: twoParticlesEnergy
    real(8) :: couplingEnergy
    real(8) :: factor, auxfactor
    integer(8) :: auxIndex
    type(vector) :: diffAB
    type(vector), allocatable :: differentOrbitals(:)
    type(vector), allocatable :: occupiedOrbitals(:,:) !! nspecies

    real(8), intent(out) :: CIenergy
    real(8) :: auxCIenergy,prefactor

    CIenergy = 0.0_8

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    if (allocated(differentOrbitals )) deallocate (differentOrbitals)
    allocate(differentOrbitals (numberOfSpecies))

  !  print *, "%%%%%%%Beginng ab", a,b
   
    do ia = 1, ConfigurationInteraction_instance%configurations(a)%nDeterminants 
      do ib = 1, ConfigurationInteraction_instance%configurations(b)%nDeterminants 

        auxCIenergy = 0.0_8
   
        call Vector_constructor (diffAB, numberOfSpecies)
 
        do i=1, numberOfSpecies

           !call vector_show(ConfigurationInteraction_instance%configurations(a)%occupations(i,ia))
           !call vector_show(ConfigurationInteraction_instance%configurations(b)%occupations(i,ib))
           diffAB%values(i)= sum ( abs ( ConfigurationInteraction_instance%configurations(a)%occupations(i,ia)%values- & 
                ConfigurationInteraction_instance%configurations(b)%occupations(i,ib)%values ) )
        end do

        call Configuration_setAtMaximumCoincidence( ConfigurationInteraction_instance%configurations(a),&
                 ConfigurationInteraction_instance%configurations(b), &
                 ia, ib, numberOfSpecies, occupiedOrbitals, factor )

        !print *, "case: ",( sum (diffAB%values) ) 
        !if ( sum (diffAB%values) > 4 ) print *, "Yina case" 
        select case ( int( sum (diffAB%values) ) )

        case (0)
          
              do i=1, numberOfSpecies
                 lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
                 kappa=MolecularSystem_getKappa(i) !exchange sign
                 charge=MolecularSystem_getCharge(i)
                 numberOfOrbitals=ConfigurationInteraction_instance%numberOfOrbitals%values(i)
                 numberOfSpatialOrbitals=numberOfOrbitals/lambda

                 do k=1, numberOfOrbitals
                    if ( ConfigurationInteraction_instance%configurations(a)%occupations(i,ia)%values(k)  > 0.0_8 ) then

                       !!Uneven orbital > alpha (0) spin
                       !!Even orbital > beta(1) spin
                       if (allocated(spin)) deallocate (spin)
                       if (allocated(spatialOrbital)) deallocate (spatialOrbital)
                       allocate(spin(2))
                       allocate(spatialOrbital(2))
                       spin = 0 
                       spatialOrbital = 0 
                       spin(1)= mod(k,lambda)
                       spatialOrbital(1)=int((k+spin(1))/lambda)

                       !One particle terms
                       auxCIenergy= auxCIenergy + &
                            ConfigurationInteraction_instance%twoCenterIntegrals(i)%values(spatialOrbital(1),spatialOrbital(1))

                       !Two particles, same specie
                       twoParticlesEnergy=0

                       do l=k+1, numberOfOrbitals 
                          if ( ConfigurationInteraction_instance%configurations(b)%occupations(i,ib)%values(l)  > 0.0_8 ) then
                             spin(2)= mod(l,lambda)
                             spatialOrbital(2)=int((l+spin(2))/lambda)

                             !Coulomb
                             auxIndex = IndexMap_tensorR4ToVector(spatialOrbital(1),spatialOrbital(1), &
                                  spatialOrbital(2),spatialOrbital(2), numberOfSpatialOrbitals )
                             twoParticlesEnergy=twoParticlesEnergy + &
                                 ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

                             !Exchange, depends on spin
                             if ( spin(1) .eq. spin(2) ) then
                                auxIndex = IndexMap_tensorR4ToVector(spatialOrbital(1),spatialOrbital(2),spatialOrbital(2),spatialOrbital(1), & 
                                     numberOfSpatialOrbitals )

                               !print *, "Exchange 0"
                                TwoParticlesEnergy=TwoParticlesEnergy + &
                                     kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

                             end if
                          end if
                       end do

                       auxCIenergy = auxCIenergy + twoParticlesEnergy*charge*charge

                       ! !Two particles, different species
                       if (numberOfSpecies > 1 ) then
                          do j=i+1, numberOfSpecies
                             lambdaOfOtherSpecie=ConfigurationInteraction_instance%lambda%values(j)
                             numberOfOtherSpecieOrbitals=ConfigurationInteraction_instance%numberOfOrbitals%values(j)
                             numberOfOtherSpecieSpatialOrbitals=numberOfOtherSpecieOrbitals/lambdaOfOtherSpecie
                             otherSpecieCharge=MolecularSystem_getCharge(j)

                             couplingEnergy=0
                             do l=1, numberOfOtherSpecieOrbitals
                                !! (j,ia) this could be a problem latter...
                                if (ConfigurationInteraction_instance%configurations(a)%occupations(j,ia)%values(l) > 0.0_8 ) then
                                   spin(2)= mod(l,lambdaOfOtherSpecie)
                                   spatialOrbital(2)=(l+spin(2))/lambdaOfOtherSpecie

                                   auxIndex = IndexMap_tensorR4ToVector(spatialOrbital(1),spatialOrbital(1),&
                                        spatialOrbital(2),spatialOrbital(2),&
                                        numberOfSpatialOrbitals, numberOfOtherSpecieSpatialOrbitals)

                                   couplingEnergy=couplingEnergy+ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1)

                                end if

                             end do

                             auxCIenergy = auxCIenergy + couplingEnergy*charge*otherSpecieCharge

                          end do

                       end if

                    end if
                 end do
              end do
              !Interaction with point charges
              auxCIenergy= auxCIenergy + HartreeFock_instance%puntualInteractionEnergy

        case (2)
           do i=1, numberOfSpecies
              lambda=ConfigurationInteraction_instance%lambda%values(i)
              kappa=MolecularSystem_getKappa(i)
              charge=MolecularSystem_getCharge(i)
              numberOfOrbitals=ConfigurationInteraction_instance%numberOfOrbitals%values(i)
              numberOfSpatialOrbitals=numberOfOrbitals/lambda

              !Determine different orbitals
              call Vector_constructor (differentOrbitals(i),2)

              differentOrbitals(i)%values=0.0
              z=0
              do k=1, numberOfOrbitals
                 if ( abs(ConfigurationInteraction_instance%configurations(a)%occupations(i,ia)%values(k)-&
                      ConfigurationInteraction_instance%configurations(b)%occupations(i,ib)%values(k)) > 0.0_8 ) then
                    z=z+1
                    differentOrbitals(i)%values(z)=k
                 end if
              end do

              !call vector_show(differentOrbitals(i))
              if ( abs( differentOrbitals(i)%values(2)) > 0.1) then !?
                 if (allocated(spin)) deallocate (spin)
                 if (allocated(spatialOrbital)) deallocate (spatialOrbital)
                 allocate(spin(3))
                 allocate(spatialOrbital(3))
                 spin = 0
                 spatialOrbital = 0

                 spin(1)=mod( int(differentOrbitals(i)%values(1)),lambda)
                 spatialOrbital(1)=int((differentOrbitals(i)%values(1)+spin(1))/lambda)
                 spin(2)=mod( int(differentOrbitals(i)%values(2)),lambda)
                 spatialOrbital(2)=int((differentOrbitals(i)%values(2)+spin(2))/lambda)

                !One particle terms
                 !if (spin(1) .eq. spin(2) ) then
                    auxCIenergy= auxCIenergy +  ConfigurationInteraction_instance%twoCenterIntegrals(i)%values( spatialOrbital(1), spatialOrbital(2) )
!                 print *, "one p", auxCIenergy
                 !end if

                 !Two particles, same specie
                 !Coulomb
                 twoParticlesEnergy=0.0_8
                 !if (spin(1) .eq. spin(2) ) then
                    do l=1, numberOfOrbitals 
                       if ( ConfigurationInteraction_instance%configurations(a)%occupations(i,ia)%values(l) > 0.0_8 .and. &
                            ConfigurationInteraction_instance%configurations(b)%occupations(i,ib)%values(l) > 0.0_8 ) then
                          spin(3)=mod(l,lambda)
                          spatialOrbital(3)=int((l+spin(3))/lambda)
                          auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(1), spatialOrbital(2), spatialOrbital(3), spatialOrbital(3), &
                               numberOfSpatialOrbitals )

                    if (spin(1) .eq. spin(2) ) then
!                  print *, "J",  ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)*charge*charge
                          TwoParticlesEnergy=TwoParticlesEnergy + &
                               ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                       end if

                       end if
                    end do
                 !end if

                 !Exchange
                 do l=1, numberOfOrbitals 
                    spin(3)=mod(l,lambda)
                    spatialOrbital(3)=int((l+spin(3))/lambda)
                    !if (spin(1) .eq. spin(3) .and. spin(2) .eq. spin(3) )
                    if (spin(1) .eq. spin(3) .and. spin(2) .eq. spin(3) ) then
                    !if (spin(1) .eq. spin(3) ) then
                       if ( ConfigurationInteraction_instance%configurations(a)%occupations(i,ia)%values(l) > 0.0_8 .and. &
                            ConfigurationInteraction_instance%configurations(b)%occupations(i,ib)%values(l) > 0.0_8 ) then

                          auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(1), spatialOrbital(3), spatialOrbital(3), spatialOrbital(2), numberOfSpatialOrbitals )

 !                 print *, "K",  ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)*charge*charge
                          TwoParticlesEnergy=TwoParticlesEnergy + &
                               kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                       end if
                    end if
                 end do

                 auxCIenergy=&
                      auxCIenergy + twoParticlesEnergy*charge*charge

                 ! !Two particles, different species

                 if (numberOfSpecies .gt. 1 .and. spin(1) .eq. spin(2) ) then
                    do j=1, numberOfSpecies

                       couplingEnergy=0
                       if (i .ne. j) then
                          lambdaOfOtherSpecie=ConfigurationInteraction_instance%lambda%values(j)
                          numberOfOtherSpecieOrbitals=ConfigurationInteraction_instance%numberOfOrbitals%values(j)
                          numberOfOtherSpecieSpatialOrbitals=numberOfOtherSpecieOrbitals/lambdaOfOtherSpecie
                          otherSpecieCharge=MolecularSystem_getCharge(j)
                          do l=1, numberOfOtherSpecieOrbitals
                             if (ConfigurationInteraction_instance%configurations(a)%occupations(j,ia)%values(l) > 0.0_8 ) then

                                spin(3)=mod(l,lambdaOfOtherSpecie)
                                spatialOrbital(3)=(l+spin(3))/lambdaOfOtherSpecie

                                auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(1), spatialOrbital(2), spatialOrbital(3), spatialOrbital(3), &
                                     numberOfSpatialOrbitals, numberOfOtherSpecieSpatialOrbitals)

                                couplingEnergy=couplingEnergy+ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1) 
                             end if
                          end do
                          auxCIenergy=&
                               auxCIenergy + couplingEnergy*charge*otherSpecieCharge

                       end if
                    end do
                 end if
              end if
              
           end do

        case (4)

           !Determine different orbitals

                 !call Configuration_checkMaximumCoincidence(ConfigurationInteraction_instance%configurations(a),&
                 !       ConfigurationInteraction_instance%configurations(b), &
                 !       factor, ia, ib, numberOfSpecies,occupiedOrbitals )


                 !call Configuration_setAtMaximumCoincidence( ConfigurationInteraction_instance%configurations(a),&
                 !       ConfigurationInteraction_instance%configurations(b), &
                 !       ia, ib, numberOfSpecies, occupiedOrbitals )
          
           do i=1, numberOfSpecies
              numberOfOrbitals=ConfigurationInteraction_instance%numberOfOrbitals%values(i)
              call Vector_constructor (differentOrbitals(i),4)
              lambda=ConfigurationInteraction_instance%lambda%values(i)
              differentOrbitals(i)%values=0.0
              numberOfOccupiedOrbitals=MolecularSystem_getOcupationNumber(i)*lambda

         !     z=0
         !     do k=1, numberOfOrbitals
         !        if ( ConfigurationInteraction_instance%configurations(a)%occupations(i,ia)%values(k)-&
         !             ConfigurationInteraction_instance%configurations(b)%occupations(i,ib)%values(k) > 0.0_8 ) then
         !           z=z+1
         !           differentOrbitals(i)%values(z)=k
         !        end if
         !        if ( ConfigurationInteraction_instance%configurations(a)%occupations(i,ia)%values(k)-&
         !             ConfigurationInteraction_instance%configurations(b)%occupations(i,ib)%values(k) < 0.0_8 ) then
         !           z=z+1
         !           differentOrbitals(i)%values(z)=-k
         !           !Negative with negative, positive with positive, to identify which orbitals belongs to the same configuration 
         !        end if
         !     end do

          ! call vector_show(differentOrbitals(i))

                 !call vector_show (occupiedOrbitals(i,1))
                 !call vector_show (occupiedOrbitals(i,2))
                  z = 1
                  do k = 1, numberOfOccupiedOrbitals
                    if ( z > 4 ) exit
                    !if ( abs(occupiedOrbitals(k,1) - occupiedOrbitals(k,2)) > 0 ) then
                    !  differentOrbitals(i)%values(z) = occupiedOrbitals(k,1)
                    !  differentOrbitals(i)%values(z+1) = -occupiedOrbitals(k,2)
                    !  z = z + 2
                    !end if 

                    if ( abs(occupiedOrbitals(i,1)%values(k) - occupiedOrbitals(i,2)%values(k)) > 0 ) then
                      differentOrbitals(i)%values(z) = occupiedOrbitals(i,1)%values(k)
                      differentOrbitals(i)%values(z+1) = -occupiedOrbitals(i,2)%values(k)
                      z = z + 2
                    end if 

                 end do 

         !  call vector_show(differentOrbitals(i))

           end do
           ! !Two cases: 4 different orbitals of the same species, and 2 and 2 of different species

           do i=1, numberOfSpecies

              lambda=ConfigurationInteraction_instance%lambda%values(i)
              kappa=MolecularSystem_getKappa(i)
              charge=MolecularSystem_getCharge(i)
              numberOfOrbitals=ConfigurationInteraction_instance%numberOfOrbitals%values(i)
              numberOfSpatialOrbitals=numberOfOrbitals/lambda

              if (allocated(spin)) deallocate (spin)
              if (allocated(spatialOrbital)) deallocate (spatialOrbital)
              allocate(spin(4))
              allocate(spatialOrbital(4))
          
              spin=0.0
              spatialOrbital = 0

              if ( abs( differentOrbitals(i)%values(3)) > 0.1) then


                 !Negative with negative, positive with positive, to identify which orbitals belongs to the same configuration
                 z=0
                 pos=0
                 neg=2

                !call vector_show(differentOrbitals(i))

                 do z=1, 4
                    if ( differentOrbitals(i)%values(z) > 0.0_8) then
                       pos=pos+1
                       spin(pos)=mod( int(differentOrbitals(i)%values(z)),lambda)
                       spatialOrbital(pos)=int((differentOrbitals(i)%values(z)+spin(pos))/lambda)
                    end if
                    if ( differentOrbitals(i)%values(z) < 0.0_8) then
                       neg=neg+1
                       spin(neg)=mod( int(-differentOrbitals(i)%values(z)),lambda)
                       spatialOrbital(neg)=int((-differentOrbitals(i)%values(z)+spin(neg))/lambda)
                    end if
                 end do

                  !print *, "spatial1", spatialOrbital(1), spatialOrbital(2), spatialOrbital(3), spatialOrbital(4)
                  !if ( spin(1) /= spin (3) .and. spin(2) /= spin(4) ) then
                  !    auxSpatialOrbital = spatialOrbital(3)
                  !    spatialOrbital(3) = spatialOrbital(4)
                  !    spatialOrbital(4) = auxSpatialOrbital
                  !    spin(3) = spin(1)
                  !    spin(4) = spin(2)
                  !end if


                 !Coulomb
                  !! 12|34
                 if ( spin(1) .eq. spin(3) .and. spin(2) .eq. spin(4) ) then
                    !auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(1), spatialOrbital(4), spatialOrbital(2), spatialOrbital(3), numberOfSpatialOrbitals )


                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(1), spatialOrbital(3), spatialOrbital(2), spatialOrbital(4), numberOfSpatialOrbitals )
!                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(1), spatialOrbital(2), spatialOrbital(3), spatialOrbital(4), numberOfSpatialOrbitals )

                    auxCIenergy = &
                         auxCIenergy + &
                         ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)*charge*charge
                 !print *, "J",  ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)*charge*charge
                 end if

                 !Exchange

                 if ( spin(1) .eq. spin(4) .and. spin(2) .eq. spin(3) ) then

                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(1), spatialOrbital(4), spatialOrbital(2), spatialOrbital(3), numberOfSpatialOrbitals )

                    auxCIenergy = &
                         auxCIenergy + &
                         kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)*charge*charge
                  !print *, "K",  ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)*charge*charge

                 end if

              end if

              !! different species

              do j=i+1, numberOfSpecies
                 if (i .ne. j) then
                    lambdaOfOtherSpecie=ConfigurationInteraction_instance%lambda%values(j)
                    numberOfOtherSpecieOrbitals=ConfigurationInteraction_instance%numberOfOrbitals%values(j)
                    numberOfOtherSpecieSpatialOrbitals=numberOfOtherSpecieOrbitals/lambdaOfOtherSpecie
                    otherSpecieCharge=MolecularSystem_getCharge(j)

                    if ( abs(differentOrbitals(i)%values(2)) .gt. 0.1 .and.  abs(differentOrbitals(j)%values(2)) .gt. 0.1) then

                       !Negative with negative, positive with positive, to identify which orbitals belongs to the same configuration

                       spin(1)=mod( int(abs(differentOrbitals(i)%values(1))),lambda)
                       spatialOrbital(1)=(abs(differentOrbitals(i)%values(1))+spin(1))/lambda

                       spin(3)=mod( int(abs(differentOrbitals(i)%values(2))),lambdaOfOtherSpecie)
                       spatialOrbital(3)=(abs(differentOrbitals(i)%values(2))+spin(3))/lambda

                       spin(2)=mod( int(abs(differentOrbitals(j)%values(1))),lambdaOfOtherSpecie)
                       spatialOrbital(2)=(abs(differentOrbitals(j)%values(1))+spin(2))/lambdaOfOtherSpecie

                       spin(4)=mod( int(abs(differentOrbitals(j)%values(2))),lambdaOfOtherSpecie)
                       spatialOrbital(4)=(abs(differentOrbitals(j)%values(2))+spin(4))/lambdaOfOtherSpecie

                       if ( spin(1) .eq. spin(3) .and. spin(2) .eq. spin(4) ) then
                          auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(1), spatialOrbital(3), &
                               spatialOrbital(2), spatialOrbital(4), numberOfSpatialOrbitals, numberOfOtherSpecieSpatialOrbitals )

                          auxCIenergy = &
                               auxCIenergy + &
                               ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1)*charge*otherSpecieCharge
                       end if


                    end if
                 end if
              end do
           end do

        case default

           auxCIenergy= 0.0_8
        end select

        !call Configuration_checkMaximumCoincidence(ConfigurationInteraction_instance%configurations(a),&
        !        ConfigurationInteraction_instance%configurations(b), &
        !        factor, ia, ib, numberOfSpecies, occupiedOrbitals )

           CIenergy= auxCIenergy * factor + CIenergy

      end do  ! ib
    end do ! ia 

    prefactor =  1.0_8/sqrt(real(ConfigurationInteraction_instance%configurations(a)%nDeterminants) * &!! (1/n!) 
             real(ConfigurationInteraction_instance%configurations(b)%nDeterminants))

    CIenergy = CIenergy * prefactor

  end subroutine ConfigurationInteraction_calculateCIenergy

  subroutine ConfigurationInteraction_calculateCIenergyB(a,b,CIenergy)
    implicit none

!    type(TransformIntegrals) :: repulsionTransformer
    integer :: numberOfSpecies
    integer :: i,j,a,b,pos,neg,ia,ib
    integer :: m,n,l,k,z,za,zb
    integer :: specieID
    integer :: otherSpecieID
    character(10) :: nameOfSpecie
    character(10) :: nameOfOtherSpecie
    integer :: numberOfOrbitals
    integer :: numberOfSpatialOrbitals
    integer :: numberOfOtherSpecieOrbitals
    integer :: numberOfOtherSpecieSpatialOrbitals
    integer, allocatable :: spin(:)
    integer, allocatable :: spatialOrbital(:)
    integer :: lambda !occupation per orbital
    integer :: lambdaOfOtherSpecie !occupation per orbital
    real(8) :: kappa !positive or negative exchange
    real(8) :: charge
    real(8) :: otherSpecieCharge
    real(8) :: twoParticlesEnergy
    real(8) :: couplingEnergy
    real(8) :: factor,energyFactor,fockFactor, integralVal
    integer(8) :: auxIndex
    type(vector) :: diffAB
!   type(vector), allocatable :: differentOrbitals(:)
    type(vector), allocatable :: auxdifferentOrbitals(:)
    type(matrix), allocatable :: differentOrbitals(:)

    real(8), intent(out) :: CIenergy
    real(8) :: auxCIenergy,prefactor

    CIenergy = 0.0_8

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    allocate(differentOrbitals (numberOfSpecies))
    allocate(auxdifferentOrbitals (numberOfSpecies))


    print *, "%%%%%%%New approach ab", a,b
   
    do ia = 1, ConfigurationInteraction_instance%configurations(a)%nDeterminants 
      do ib = 1, ConfigurationInteraction_instance%configurations(b)%nDeterminants 

        auxCIenergy = 0.0_8
   
        call Vector_constructor (diffAB, numberOfSpecies)
 
        do i=1, numberOfSpecies

           !call vector_show(ConfigurationInteraction_instance%configurations(a)%occupations(i,ia))
           !call vector_show(ConfigurationInteraction_instance%configurations(b)%occupations(i,ib))
           diffAB%values(i)= sum ( abs ( ConfigurationInteraction_instance%configurations(a)%occupations(i,ia)%values- & 
                ConfigurationInteraction_instance%configurations(b)%occupations(i,ib)%values ) )
        end do

        print *, "diffAB", ia,ib
        call vector_show(diffAB)

        !Determine different orbitals
        
        do i=1, numberOfSpecies
           numberOfOrbitals=ConfigurationInteraction_instance%numberOfOrbitals%values(i)
           
           call Matrix_constructor (differentOrbitals(i),int(2,8), int(4,8)) !! 1: \langle Phi | , 2: | \Phi \rangle 
           call Vector_constructor (auxdifferentOrbitals(i),8)
           differentOrbitals(i)%values=0.0_8
           auxdifferentOrbitals(i)%values=0.0_8

           !!za=1
           !!zb=2
           !!do k=1, numberOfOrbitals
           !!   if ( ConfigurationInteraction_instance%configurations(a)%occupations(i,ia)%values(k)-&
           !!        ConfigurationInteraction_instance%configurations(1)%occupations(i,1)%values(k) < 0.0_8 ) then
           !!      differentOrbitals(i)%values(1,za)=k
           !!      za=za+2
           !!   end if
           !!   if ( ConfigurationInteraction_instance%configurations(a)%occupations(i,ia)%values(k)-&
           !!        ConfigurationInteraction_instance%configurations(1)%occupations(i,1)%values(k) > 0.0_8 ) then
           !!      differentOrbitals(i)%values(1,zb)=k
           !!      zb=zb+2
           !!   end if

           !!end do

           !!za=1
           !!zb=2
!!         !!  if ( ConfigurationInteraction_instance%configurations(b)%order%values(i) == 2 ) z=0

           !!do k=1, numberOfOrbitals
           !!   if ( ConfigurationInteraction_instance%configurations(b)%occupations(i,ib)%values(k)-&
           !!        ConfigurationInteraction_instance%configurations(1)%occupations(i,1)%values(k) < 0.0_8 ) then
           !!      print *, "za", za
           !!      differentOrbitals(i)%values(2,za)=k
           !!      za=za+2
           !!   end if
           !!   if ( ConfigurationInteraction_instance%configurations(b)%occupations(i,ib)%values(k)-&
           !!        ConfigurationInteraction_instance%configurations(1)%occupations(i,1)%values(k) > 0.0_8 ) then
           !!      print *, "zb", zb
           !!      differentOrbitals(i)%values(2,zb)=k
           !!      zb=zb+2
           !!   end if
           !!end do

          !!print *, "diff orbitals1"
          !!call Matrix_show(differentOrbitals(i))
  
         ! do z = 1, ConfigurationInteraction_instance%configurations(a)%order%values(i)*2
         !     print *, "za", int(z+mod(z,2)/2)
         !     differentOrbitals(i)%values(1,z) = ConfigurationInteraction_instance%configurations(a)%excitations(i,ia, &
         !                                          mod(z,2))%values(z)
         ! end do
          z = 0
          do zb = 1, int(ConfigurationInteraction_instance%configurations(a)%order%values(i))
            do za = 1, 2
               z = z + 1 
               differentOrbitals(i)%values(1,z) = ConfigurationInteraction_instance%configurations(a)%excitations(i,ia, &
                                                   za)%values(zb)
            end do
          end do


          z = 0
          do zb = 1, int(ConfigurationInteraction_instance%configurations(b)%order%values(i))
            do za = 1, 2
               z = z + 1 
               differentOrbitals(i)%values(2,z) = ConfigurationInteraction_instance%configurations(b)%excitations(i,ib, &
                                                   za)%values(zb)
            end do
          end do

          print *, "diff orbitals2"
          call Matrix_show(differentOrbitals(i))
        end do 


        do i=1, numberOfSpecies
           lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
           kappa=MolecularSystem_getKappa(i) !exchange sign
           charge=MolecularSystem_getCharge(i)
           numberOfOrbitals=ConfigurationInteraction_instance%numberOfOrbitals%values(i)
           numberOfSpatialOrbitals=numberOfOrbitals/lambda

           !!Uneven orbital > alpha (0) spin
           !!Even orbital > beta(1) spin
           if (allocated(spin)) deallocate (spin)
           if (allocated(spatialOrbital)) deallocate (spatialOrbital)
           allocate(spin(8))
           allocate(spatialOrbital(8))
           spin = 0 
           spatialOrbital = 0 
           z = 0
           do za = 1, 2
             do zb = 1, 2 !size(differentOrbitals(i)%values, dim=2 )
             !!do zb = 1, size(differentOrbitals(i)%values, dim=2 )
                  if ( differentOrbitals(i)%values(za,zb) /= 0.0_8  ) then
                  z = z + 1
                  !print *, "z", z, differentOrbitals(i)%values(za,zb)

                  auxdifferentOrbitals(i)%values(z) =   differentOrbitals(i)%values(za,zb)
                  spin(z)= mod(int(differentOrbitals(i)%values(za,zb)),lambda)
                  spatialOrbital(z)=int((differentOrbitals(i)%values(za,zb)+spin(z))/lambda)
                  end if
             end do
           end do

           do za = 1, 2
             do zb = 3, 4 !size(differentOrbitals(i)%values, dim=2 )
             !!do zb = 1, size(differentOrbitals(i)%values, dim=2 )
                  if ( differentOrbitals(i)%values(za,zb) /= 0.0_8  ) then
                  z = z + 1
                  !print *, "z", z, differentOrbitals(i)%values(za,zb)

                  auxdifferentOrbitals(i)%values(z) =   differentOrbitals(i)%values(za,zb)
                  spin(z)= mod(int(differentOrbitals(i)%values(za,zb)),lambda)
                  spatialOrbital(z)=int((differentOrbitals(i)%values(za,zb)+spin(z))/lambda)
                  end if
             end do
           end do
            call vector_show(auxdifferentOrbitals(i))
            energyFactor = 0.0_8

          if ( z == 6 ) then !! Single-double

           !! k,c,j,b,i,a
           !! 1,2,5,6,3,4
            fockFactor = 0.0_8

            !!+F_{ia} \delta_{b,c} \delta_{j,k}
            !!+F_{34} \DELTA_{6,2} \DELTA_{5,1}
            if (auxdifferentOrbitals(i)%values(6) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(5) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor + ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(3), spatialOrbital(4))
            !print *, "FFockk01", fockFactor

            !!-F_{ib} \delta_{a,c} \delta_{j,k}
            !!-F_{36} \DELTA_{4,2} \DELTA_{5,1}
            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(5) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor - ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(3), spatialOrbital(6))
           ! print *, "FFockk02", fockFactor

            !!-F_{ja} \delta_{b,c} \delta_{i,k}
            !!-F_{54} \DELTA_{6,2} \DELTA_{3,1}
            if (auxdifferentOrbitals(i)%values(6) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor - ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(5), spatialOrbital(4))
            !print *, "FFockk03", fockFactor

            !!+F_{jb} \delta_{a,c} \delta_{i,k}
            !!+F_{56} \DELTA_{4,2} \DELTA_{3,1}
            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor + ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(5), spatialOrbital(6))
            !print *, "FFockk04", fockFactor

            auxCIenergy = auxCIenergy + fockFactor

!!+1.0 \langle i j || k a \rangle \delta_{b,c} 
!!+1.0 \LANGLE 3 5 || 1 4 \RANGLE \DELTA_{6,2} 
            if (auxdifferentOrbitals(i)%values(6) == auxdifferentOrbitals(i)%values(2) ) then
                  if ( spin(3) == spin(1) .and. spin(5) == spin(4) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(3), spatialOrbital(1), &
                               spatialOrbital(5), spatialOrbital(4), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(3) == spin(4) .and. spin(5) == spin(1) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(3), spatialOrbital(4), &
                                spatialOrbital(5), spatialOrbital(1), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

            !print *, "auxCIenergy + JK01", auxCIenergy

!!-1.0 \langle i j || k b \rangle \delta_{a,c} 
!!-1.0 \LANGLE 3 5 || 1 6 \RANGLE \DELTA_{4,2} 
            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) ) then
                  if ( spin(3) == spin(6) .and. spin(5) == spin(1) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(3), spatialOrbital(6), &
                               spatialOrbital(5), spatialOrbital(1), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(3) == spin(1) .and. spin(5) == spin(6) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(3), spatialOrbital(1), &
                                spatialOrbital(5), spatialOrbital(6), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

            !print *, "auxCIenergy + JK02", auxCIenergy

!!-1.0 \LANGLE 2 3 || 4 6 \RANGLE \DELTA_{5,1} 
!!-1.0 \langle c i || a b \rangle \delta_{j,k} 

            if (auxdifferentOrbitals(i)%values(5) == auxdifferentOrbitals(i)%values(1) ) then
                  if ( spin(2) == spin(4) .and. spin(3) == spin(6) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(4), &
                               spatialOrbital(3), spatialOrbital(6), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy - ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(2) == spin(6) .and. spin(3) == spin(4) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(6), &
                                spatialOrbital(3), spatialOrbital(4), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy - kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

            !print *, "auxCIenergy + JK03", auxCIenergy

!!+1.0 \langle c j || a b \rangle \delta_{i,k} 
!!+1.0 \LANGLE 2 5 || 4 6 \RANGLE \DELTA_{3,1} 
            if (auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) then
                  if ( spin(2) == spin(6) .and. spin(5) == spin(4) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(6), &
                               spatialOrbital(5), spatialOrbital(4), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy - ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  !print *, "auxCIenergy + J04", auxCIenergy
                  !print *, "2 6 | 5 4 ", spatialOrbital(2), spatialOrbital(6), spatialOrbital(5), spatialOrbital(4)
                  end if 
                  if ( spin(2) == spin(4) .and. spin(5) == spin(6) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(4), &
                                spatialOrbital(5), spatialOrbital(6), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy - kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

            !print *, "auxCIenergy + JK04", auxCIenergy



          else 

           !! k,c,l,d,j,b,i,a
           !! 1,2,5,6,7,8,3,4
           !! ad,il,bc,jk
            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(1) ) energyFactor = 1.0_8

            !! bd,il,ac,jk
            if (auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(1) ) energyFactor = -1.0_8

            !! ad,jl,bc,ik
            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) energyFactor = -1.0_8

            !! bd,jl,ac,ik
            if (auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) energyFactor = 1.0_8


!           if ( sum (diffAB%values ) == 0 ) then
            if ( energyFactor /= 0 ) then

                auxCIenergy = auxCIenergy + energyFactor*HartreeFock_instance%totalEnergy/numberOfSpecies

           end if

           ! print *, "auxCIenergy + HF", auxCIenergy

           !! k,c,l,d,j,b,i,a
           !! 1,2,5,6,7,8,3,4
            fockFactor = 0.0_8

            !! +ik ad,jl,bc
            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(2) ) &
                 fockFactor = fockFactor + ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(3), spatialOrbital(1))
            !print *, "FFockk01", fockFactor

            !! -F_{ik}\delta_{b,d} \delta_{j,l} \delta_{a,c}
            !! -F_{31}\DELTA_{8,6} \DELTA_{7,5} \DELTA_{4,2}
            if (auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) ) &
                 fockFactor = fockFactor -ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(3), spatialOrbital(1))
            !print *, "FFockk02", fockFactor

            !!-F_{jk} \delta_{a,d} \delta_{i,l} \delta_{b,c}
            !!-F_{71} \DELTA_{4,6} \DELTA_{3,5} \DELTA_{8,2}

            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(2) ) &
                 fockFactor = fockFactor -ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(7), spatialOrbital(1))

            !print *, "FFockk03", fockFactor
            !!+F_{jk} \delta_{b,d} \delta_{i,l} \delta_{a,c}
            !!+F_{71} \DELTA_{8,6} \DELTA_{3,5} \DELTA_{4,2}

            if (auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) ) &
                 fockFactor = fockFactor + ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(7), spatialOrbital(1))
            !print *, "FFockk04", fockFactor

            !!+F_{jl} \delta_{a,d} \delta_{b,c} \delta_{i,k}
            !!+F_{75} \DELTA_{4,6} \DELTA_{8,2} \DELTA_{3,1}

            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor + ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(7), spatialOrbital(5))

            !print *, "FFockk05", fockFactor
            !!-F_{jl} \delta_{b,d} \delta_{a,c} \delta_{i,k}
            !!-F_{75} \DELTA_{8,6} \DELTA_{4,2} \DELTA_{3,1}

            if (auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor -ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(7), spatialOrbital(5))

            !print *, "FFockk06", fockFactor
            !!-F_{il} \delta_{a,d} \delta_{b,c} \delta_{j,k}
            !!-F_{35} \DELTA_{4,6} \DELTA_{8,2} \DELTA_{7,1}

            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor -ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(3), spatialOrbital(5))

            !print *, "FFockk07", fockFactor
            !!+F_{il} \delta_{b,d} \delta_{a,c} \delta_{j,k}
            !!+F_{35} \DELTA_{8,6} \DELTA_{4,2} \DELTA_{7,1}

            if (auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor + ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(3), spatialOrbital(5))
            !print *, "FFockk08", fockFactor

            !!+F_{da} \delta_{i,l} \delta_{b,c} \delta_{j,k}
            !!+F_{64} \DELTA_{3,5} \DELTA_{8,2} \DELTA_{7,1}

            if (auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor + ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(6), spatialOrbital(4))

            !print *, "FFockk09", fockFactor
            !!-F_{da} \delta_{j,l} \delta_{b,c} \delta_{i,k}
            !!-F_{64} \DELTA_{7,5} \DELTA_{8,2} \DELTA_{3,1}

            if (auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor -ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(6), spatialOrbital(4))
            !print *, "FFockk10", fockFactor

            !!-F_{cb} \delta_{a,d} \delta_{j,l} \delta_{i,k}
            !!-F_{28} \DELTA_{4,6} \DELTA_{7,5} \DELTA_{3,1}

            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor -ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(2), spatialOrbital(8))

            !print *, "FFockk11", fockFactor
            !!+F_{cb} \delta_{a,d} \delta_{i,l} \delta_{j,k}
            !!+F_{28} \DELTA_{4,6} \DELTA_{3,5} \DELTA_{7,1}

            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor + ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(2), spatialOrbital(8))

            !print *, "FFockk12", fockFactor
            !!+F_{db} \delta_{j,l} \delta_{a,c} \delta_{i,k}
            !!+F_{68} \DELTA_{7,5} \DELTA_{4,2} \DELTA_{3,1}

            if (auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor + ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(6), spatialOrbital(8))

            !print *, "FFockk13", fockFactor
            !!-F_{db} \delta_{i,l} \delta_{a,c} \delta_{j,k}
            !!-F_{68} \DELTA_{3,5} \DELTA_{4,2} \DELTA_{7,1}

            if (auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) .and. &
              auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor -ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(6), spatialOrbital(8))

            !print *, "FFockk14", fockFactor
            !!!-F_{ca} \delta_{b,d} \delta_{i,l} \delta_{j,k}
            !!-F_{24} \DELTA_{8,6} \DELTA_{3,5} \DELTA_{7,1}

            if (auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor -ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(2), spatialOrbital(4))
            !print *, "FFockk15", fockFactor

            !!+F_{24} \DELTA_{8,6} \DELTA_{7,5} \DELTA_{3,1}
            !!+F_{ca} \delta_{b,d} \delta_{j,l} \delta_{i,k}

            if (auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(6) .and. &
              auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(5) .and. &
              auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) &
                 fockFactor = fockFactor + ConfigurationInteraction_instance%fockMatrix(i)%values( spatialOrbital(2), spatialOrbital(4))
            !print *, "FFockk16", fockFactor


            if ( fockFactor /= 0 ) then
                auxCIenergy = auxCIenergy + fockFactor

           end if

            !print *, "auxCIenergy + FFockk", auxCIenergy
           !! k,c,l,d,j,b,i,a
           !! 1,2,5,6,7,8,3,4
          !print *, "aux"
          !call vector_show(auxdifferentOrbitals(i))

          if ( product(auxdifferentOrbitals(i)%values(1:4)) /=0 ) then

!!-1.0 \langle c i || a k \rangle \delta_{b,d} \delta_{j,l}
!!-1.0 \LANGLE 2 3 || 4 1 \RANGLE \DELTA_{8,6} \DELTA_{7,5}

            if (auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(6) .and. &
                auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(5) ) then

                  if ( spin(2) == spin(1) .and. spin(3) == spin(4) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(1), &
                               spatialOrbital(3), spatialOrbital(4), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
!print*,auxdifferentOrbitals(i)%values(2), auxdifferentOrbitals(i)%values(1), auxdifferentOrbitals(i)%values(3), auxdifferentOrbitals(i)%values(4)
!            print *, "auxCIenergy + J1", auxCIenergy
                  end if 

                  if ( spin(2) == spin(4) .and. spin(3) == spin(1) ) then
!print*,auxdifferentOrbitals(i)%values(2), auxdifferentOrbitals(i)%values(4), auxdifferentOrbitals(i)%values(3), auxdifferentOrbitals(i)%values(1)
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(4), &
                                spatialOrbital(3), spatialOrbital(1), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if
!            print *, "auxCIenergy + JK1", auxCIenergy
!!+1.0 \langle c i || b k \rangle \delta_{a,d} \delta_{j,l}
!!+1.0 \LANGLE 2 3 || 8 1 \RANGLE \DELTA_{4,6} \DELTA_{7,5}

            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(6) .and. &
                auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(5) ) then
                  if ( spin(2) == spin(8) .and. spin(3) == spin(1) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(8), &
                               spatialOrbital(3), spatialOrbital(1), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(2) == spin(1) .and. spin(3) == spin(8) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(1), &
                                spatialOrbital(3), spatialOrbital(8), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

!            print *, "auxCIenergy + JK2", auxCIenergy
!!+1.0 \langle c j || a k \rangle \delta_{b,d} \delta_{i,l}
!!+1.0 \LANGLE 2 7 || 4 1 \RANGLE \DELTA_{8,6} \DELTA_{3,5}

            if (auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(6) .and. &
                auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(5) ) then
                  if ( spin(2) == spin(4) .and. spin(7) == spin(1) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(4), &
                               spatialOrbital(7), spatialOrbital(1), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(2) == spin(1) .and. spin(7) == spin(4) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(1), &
                                spatialOrbital(7), spatialOrbital(4), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

!            print *, "auxCIenergy + JK3", auxCIenergy
!!-1.0 \langle c j || b k \rangle \delta_{a,d} \delta_{i,l}
!!-1.0 \LANGLE 2 7 || 8 1 \RANGLE \DELTA_{4,6} \DELTA_{3,5}

            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(6) .and. &
                auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(5) ) then
                  if ( spin(2) == spin(1) .and. spin(7) == spin(8) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(1), &
                               spatialOrbital(7), spatialOrbital(8), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(2) == spin(8) .and. spin(7) == spin(1) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(8), &
                                spatialOrbital(7), spatialOrbital(1), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

!            print *, "auxCIenergy + JK4", auxCIenergy
!!-1.0 \langle i j || k l \rangle \delta_{a,d} \delta_{b,c}
!!-1.0 \LANGLE 3 7 || 1 5 \RANGLE \DELTA_{4,6} \DELTA_{8,2}

            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(6) .and. &
                auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(2) ) then
                  if ( spin(3) == spin(5) .and. spin(7) == spin(1) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(3), spatialOrbital(5), &
                               spatialOrbital(7), spatialOrbital(1), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(3) == spin(1) .and. spin(7) == spin(5) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(3), spatialOrbital(1), &
                                spatialOrbital(7), spatialOrbital(5), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

!            print *, "auxCIenergy + JK5", auxCIenergy
!!+1.0 \langle d i || a k \rangle \delta_{j,l} \delta_{b,c}
!!+1.0 \LANGLE 6 3 || 4 1 \RANGLE \DELTA_{7,5} \DELTA_{8,2}

            if (auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(5) .and. &
                auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(2) ) then
                  if ( spin(6) == spin(4) .and. spin(3) == spin(1) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(4), &
                               spatialOrbital(3), spatialOrbital(1), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(6) == spin(1) .and. spin(3) == spin(4) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(1), &
                                spatialOrbital(3), spatialOrbital(4), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

!            print *, "auxCIenergy + JK6", auxCIenergy
!!-1.0 \langle d j || a k \rangle \delta_{i,l} \delta_{b,c}
!!-1.0 \LANGLE 6 7 || 4 1 \RANGLE \DELTA_{3,5} \DELTA_{8,2}

            if (auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(5) .and. &
                auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(2) ) then
                  if ( spin(6) == spin(1) .and. spin(7) == spin(4) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(1), &
                               spatialOrbital(7), spatialOrbital(4), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(6) == spin(4) .and. spin(7) == spin(1) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(4), &
                                spatialOrbital(7), spatialOrbital(1), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

!            print *, "auxCIenergy + JK7", auxCIenergy
!!+1.0 \langle i j || k l \rangle \delta_{b,d} \delta_{a,c}
!!+1.0 \LANGLE 3 7 || 1 5 \RANGLE \DELTA_{8,6} \DELTA_{4,2}
!! **********
            if (auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(6) .and. &
                auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) ) then
                  if ( spin(3) == spin(1) .and. spin(7) == spin(5) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(3), spatialOrbital(1), &
                    spatialOrbital(7), spatialOrbital(5), numberOfSpatialOrbitals )

                    if (auxIndex > 0 ) then
                  
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
!print*,auxdifferentOrbitals(i)%values(3), auxdifferentOrbitals(i)%values(1), auxdifferentOrbitals(i)%values(7), auxdifferentOrbitals(i)%values(5)
!            print *, "auxCIenergy + J8", auxCIenergy
                  end if 
                  if ( spin(3) == spin(5) .and. spin(7) == spin(1) ) then
           auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(3), spatialOrbital(5), &
                               spatialOrbital(7), spatialOrbital(1), numberOfSpatialOrbitals )

                  if (auxIndex > 0 ) then
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  end if 
                  end if 
            end if

!            print *, "auxCIenergy + JK8", auxCIenergy
!!-1.0 \langle d i || b k \rangle \delta_{j,l} \delta_{a,c}
!!-1.0 \LANGLE 6 3 || 8 1 \RANGLE \DELTA_{7,5} \DELTA_{4,2}
!! ****************
            if (auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(5) .and. &
                auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) ) then
                  if ( spin(6) == spin(8) .and. spin(3) == spin(1) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(8), &
                                spatialOrbital(3), spatialOrbital(1), numberOfSpatialOrbitals )
                  if (auxIndex > 0 ) then
                  auxCIenergy = auxCIenergy - ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

!print*,auxdifferentOrbitals(i)%values(6), auxdifferentOrbitals(i)%values(8), auxdifferentOrbitals(i)%values(3), auxdifferentOrbitals(i)%values(1)
!            print *, "auxCIenergy + J9", auxCIenergy
                  end if 
                  if ( spin(6) == spin(1) .and. spin(3) == spin(8) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(1), &
                               spatialOrbital(3), spatialOrbital(8), numberOfSpatialOrbitals )

                  if (auxIndex > 0 ) then
                    auxCIenergy = auxCIenergy - kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  end if 
                  end if 
            end if

!            print *, "auxCIenergy + JK9", auxCIenergy
!!+1.0 \langle d j || b k \rangle \delta_{i,l} \delta_{a,c}
!!+1.0 \LANGLE 6 7 || 8 1 \RANGLE \DELTA_{3,5} \DELTA_{4,2}

            if (auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(5) .and. &
                auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) ) then
                  if ( spin(6) == spin(8) .and. spin(7) == spin(1) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(8), &
                               spatialOrbital(7), spatialOrbital(1), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(6) == spin(1) .and. spin(7) == spin(8) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(1), &
                                spatialOrbital(7), spatialOrbital(8), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

!            print *, "auxCIenergy + JK10", auxCIenergy
!!+1.0 \langle c i || a l \rangle \delta_{b,d} \delta_{j,k}
!!+1.0 \LANGLE 2 3 || 4 5 \RANGLE \DELTA_{8,6} \DELTA_{7,1}

            if (auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(6) .and. &
                auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(1) ) then
                  if ( spin(2) == spin(4) .and. spin(3) == spin(5) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(4), &
                               spatialOrbital(3), spatialOrbital(5), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(2) == spin(5) .and. spin(3) == spin(4) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(5), &
                                spatialOrbital(3), spatialOrbital(4), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

!            print *, "auxCIenergy + JK11", auxCIenergy
!!-1.0 \langle c i || b l \rangle \delta_{a,d} \delta_{j,k}
!!-1.0 \LANGLE 2 3 || 8 5 \RANGLE \DELTA_{4,6} \DELTA_{7,1}

            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(6) .and. &
                auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(1) ) then
                  if ( spin(2) == spin(5) .and. spin(3) == spin(8) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(5), &
                               spatialOrbital(3), spatialOrbital(8), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(2) == spin(8) .and. spin(3) == spin(5) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(8), &
                                spatialOrbital(3), spatialOrbital(5), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

!            print *, "auxCIenergy + JK12", auxCIenergy
!!-1.0 \langle c d || a b \rangle \delta_{i,l} \delta_{j,k}
!!-1.0 \LANGLE 2 6 || 4 8 \RANGLE \DELTA_{3,5} \DELTA_{7,1}

            if (auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(5) .and. &
                auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(1) ) then
                  if ( spin(2) == spin(8) .and. spin(6) == spin(4) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(8), &
                               spatialOrbital(6), spatialOrbital(4), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(2) == spin(4) .and. spin(6) == spin(8) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(4), &
                                spatialOrbital(6), spatialOrbital(8), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

!           print *, "auxCIenergy + JK13", auxCIenergy
!!-1.0 \langle d i || a l \rangle \delta_{b,c} \delta_{j,k}
!!-1.0 \LANGLE 6 3 || 4 5 \RANGLE \DELTA_{8,2} \DELTA_{7,1}

            if (auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(2) .and. &
                auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(1) ) then
                  if ( spin(6) == spin(5) .and. spin(3) == spin(4) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(5), &
                               spatialOrbital(3), spatialOrbital(4), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(6) == spin(4) .and. spin(3) == spin(5) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(4), &
                                spatialOrbital(3), spatialOrbital(5), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

!            print *, "auxCIenergy + JK14", auxCIenergy
!!+1.0 \langle d i || b l \rangle \delta_{a,c} \delta_{j,k}
!!+1.0 \LANGLE 6 3 || 8 5 \RANGLE \DELTA_{4,2} \DELTA_{7,1}

            if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) .and. &
                auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(1) ) then
                  if ( spin(6) == spin(8) .and. spin(3) == spin(5) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(8), &
                               spatialOrbital(3), spatialOrbital(5), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(6) == spin(5) .and. spin(3) == spin(8) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(5), &
                                spatialOrbital(3), spatialOrbital(8), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

            !print *, "auxCIenergy + JK15", auxCIenergy
!!-1.0 \langle c j || a l \rangle \delta_{b,d} \delta_{i,k}
!!-1.0 \LANGLE 2 7 || 4 5 \RANGLE \DELTA_{8,6} \DELTA_{3,1}
!!*************
            if (auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(6) .and. &
                auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) then
                  if ( spin(2) == spin(4) .and. spin(7) == spin(5) ) then
               auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(4), &
                                spatialOrbital(7), spatialOrbital(5), numberOfSpatialOrbitals )
                  if (auxIndex > 0 ) then
                  auxCIenergy = auxCIenergy - ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
!print*,auxdifferentOrbitals(i)%values(2), auxdifferentOrbitals(i)%values(4), auxdifferentOrbitals(i)%values(7), auxdifferentOrbitals(i)%values(5)
            !print *, "auxCIenergy + J16", auxCIenergy
                  end if 
                  if ( spin(2) == spin(5) .and. spin(7) == spin(4) ) then
                                   auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(5), &
                               spatialOrbital(7), spatialOrbital(4), numberOfSpatialOrbitals )

                  if (auxIndex > 0 ) then
                    auxCIenergy = auxCIenergy - kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  end if 
                  end if 
            end if

            !print *, "auxCIenergy + JK16", auxCIenergy
!!+1.0 \langle c j || b l \rangle \delta_{a,d} \delta_{i,k}
!!+1.0 \LANGLE 2 7 || 8 5 \RANGLE \DELTA_{4,6} \DELTA_{3,1}

           if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(6) .and. &
                auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) then
                  if ( spin(2) == spin(8) .and. spin(7) == spin(5) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(8), &
                               spatialOrbital(7), spatialOrbital(5), numberOfSpatialOrbitals )
                  if (auxIndex > 0 ) then
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  end if 
                  if ( spin(2) == spin(7) .and. spin(8) == spin(5) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(7), &
                                spatialOrbital(8), spatialOrbital(5), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

            !print *, "auxCIenergy + JK17", auxCIenergy
!!+1.0 \langle c d || a b \rangle \delta_{j,l} \delta_{i,k}
!!+1.0 \LANGLE 2 6 || 4 8 \RANGLE \DELTA_{7,5} \DELTA_{3,1}

           if (auxdifferentOrbitals(i)%values(7) == auxdifferentOrbitals(i)%values(5) .and. &
                auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) then
                  if ( spin(2) == spin(4) .and. spin(6) == spin(8) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(4), &
                               spatialOrbital(6), spatialOrbital(8), numberOfSpatialOrbitals )
                  if (auxIndex > 0 ) then
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
!print*,auxdifferentOrbitals(i)%values(2), auxdifferentOrbitals(i)%values(4), auxdifferentOrbitals(i)%values(6), auxdifferentOrbitals(i)%values(8)
            !print *, "auxCIenergy + J18", auxCIenergy
                  end if 
                  end if 
                  if ( spin(2) == spin(8) .and. spin(6) == spin(4) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(8), &
                                spatialOrbital(6), spatialOrbital(4), numberOfSpatialOrbitals )
                  if (auxIndex > 0 ) then
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  end if 
            end if

            !print *, "auxCIenergy + JK18", auxCIenergy
!!+1.0 \langle d j || a l \rangle \delta_{b,c} \delta_{i,k}
!!+1.0 \LANGLE 6 7 || 4 5 \RANGLE \DELTA_{8,2} \DELTA_{3,1}

           if (auxdifferentOrbitals(i)%values(8) == auxdifferentOrbitals(i)%values(2) .and. &
                auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) then
                  if ( spin(6) == spin(4) .and. spin(7) == spin(5) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(4), &
                               spatialOrbital(7), spatialOrbital(5), numberOfSpatialOrbitals )
                  auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  if ( spin(6) == spin(5) .and. spin(7) == spin(4) ) then
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(5), &
                                spatialOrbital(7), spatialOrbital(4), numberOfSpatialOrbitals )
                    auxCIenergy = auxCIenergy + kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
            end if

            !print *, "auxCIenergy + JK19", auxCIenergy
!!-1.0 \langle d j || b l \rangle \delta_{a,c} \delta_{i,k}
!!-1.0 \LANGLE 6 7 || 8 5 \RANGLE \DELTA_{4,2} \DELTA_{3,1}

           if (auxdifferentOrbitals(i)%values(4) == auxdifferentOrbitals(i)%values(2) .and. &
                auxdifferentOrbitals(i)%values(3) == auxdifferentOrbitals(i)%values(1) ) then
                  if ( spin(6) == spin(8) .and. spin(7) == spin(5) ) then
                  auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(8), &
                               spatialOrbital(7), spatialOrbital(5), numberOfSpatialOrbitals )
                  if (auxIndex > 0 ) then
                  auxCIenergy = auxCIenergy - ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
!print*,auxdifferentOrbitals(i)%values(6), auxdifferentOrbitals(i)%values(8), auxdifferentOrbitals(i)%values(7), auxdifferentOrbitals(i)%values(5)
            !print *, "auxCIenergy + J20", auxCIenergy
                  end if 
                  end if 
                  if ( spin(6) == spin(5) .and. spin(7) == spin(8) ) then
print*,auxdifferentOrbitals(i)%values(6), auxdifferentOrbitals(i)%values(5), auxdifferentOrbitals(i)%values(7), auxdifferentOrbitals(i)%values(8)
                    auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(6), spatialOrbital(5), &
                                spatialOrbital(7), spatialOrbital(8), numberOfSpatialOrbitals )
                  if (auxIndex > 0 ) then
                    auxCIenergy = auxCIenergy - kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                  end if 
                  end if 
            end if
            !print *, "auxCIenergy + JK20", auxCIenergy

           !!if ( product(auxdifferentOrbitals(i)%values(1:4)) /=0 ) then

           !!         ! call vector_show(differentOrbitals(i))
           !!         ! i -> a = 1->2 , j-> b = 3->4 , [ai|jb] = [ia|jb] -> [12|34], [ab|ji] = [ab|ij] -> [24|13]
           !!         auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(1), spatialOrbital(2), &
           !!                    spatialOrbital(3), spatialOrbital(4), numberOfSpatialOrbitals )

           !!          auxCIenergy = auxCIenergy + &
           !!                    ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)*charge*charge
           !!            print *, "auxCIenergy + J", auxCIenergy

           !!          if ( spin(1) == spin(3) ) then
           !!         auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(2), spatialOrbital(4), &
           !!                    spatialOrbital(1), spatialOrbital(3), numberOfSpatialOrbitals )

           !!          auxCIenergy = auxCIenergy + &
           !!                     kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)*charge*charge

           !!            print *, "auxCIenergy + JK", auxCIenergy
           !!          end if 

           end if
          end if
            !! different species

            do j=i+1, numberOfSpecies
               if (i .ne. j) then
                  lambdaOfOtherSpecie=ConfigurationInteraction_instance%lambda%values(j)
                  numberOfOtherSpecieOrbitals=ConfigurationInteraction_instance%numberOfOrbitals%values(j)
                  numberOfOtherSpecieSpatialOrbitals=numberOfOtherSpecieOrbitals/lambdaOfOtherSpecie
                  otherSpecieCharge=MolecularSystem_getCharge(j)

                  !!! check this!
                  if ( abs(auxdifferentOrbitals(i)%values(2)) .gt. 0.1 .and.  abs(auxdifferentOrbitals(j)%values(4)) .gt. 0.1) then

                     !Negative with negative, positive with positive, to identify which orbitals belongs to the same configuration

                     spin(1)=mod( int(abs(auxdifferentOrbitals(i)%values(1))),lambda)
                     spatialOrbital(1)=(abs(auxdifferentOrbitals(i)%values(1))+spin(1))/lambda

                     spin(3)=mod( int(abs(auxdifferentOrbitals(i)%values(2))),lambdaOfOtherSpecie)
                     spatialOrbital(3)=(abs(auxdifferentOrbitals(i)%values(2))+spin(3))/lambda

                     spin(2)=mod( int(abs(auxdifferentOrbitals(j)%values(3))),lambdaOfOtherSpecie)
                     spatialOrbital(2)=(abs(auxdifferentOrbitals(j)%values(3))+spin(2))/lambdaOfOtherSpecie

                     spin(4)=mod( int(abs(auxdifferentOrbitals(j)%values(4))),lambdaOfOtherSpecie)
                     spatialOrbital(4)=(abs(auxdifferentOrbitals(j)%values(4))+spin(4))/lambdaOfOtherSpecie
                      !! I,A,j,b, <Ij|Ab>

                     if ( spin(1) .eq. spin(3) .and. spin(2) .eq. spin(4) ) then

                        auxIndex = IndexMap_tensorR4ToVector( spatialOrbital(1), spatialOrbital(3), &
                             spatialOrbital(2), spatialOrbital(4), numberOfSpatialOrbitals, numberOfOtherSpecieSpatialOrbitals )

                        auxCIenergy = auxCIenergy + &
                             ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1)*charge*otherSpecieCharge
                     end if


                  end if
               end if
            end do

            !print *, "auxCIenergy + Cou", auxCIenergy

          !print *, "aux energy for specie", i, auxCIenergy
 
        end do !! i

        !Interaction with point charges. This is already inclued in the Fock matrix.
!        auxCIenergy = auxCIenergy + HartreeFock_instance%puntualInteractionEnergy


        !call Configuration_checkMaximumCoincidence(ConfigurationInteraction_instance%configurations(a),&
        !        ConfigurationInteraction_instance%configurations(b), &
        !        factor, ia, ib, numberOfSpecies)

                 CIenergy= auxCIenergy + CIenergy
                 print *, "ia,ib,auxCIenergy",ia,ib,auxCIenergy
                 !CIenergy= auxCIenergy * factor + CIenergy

      end do  ! ib
    end do ! ia 

    prefactor =  1.0_8/sqrt(real(ConfigurationInteraction_instance%configurations(a)%nDeterminants))
    print *, "prefactor1", prefactor
    prefactor =  1.0_8/sqrt(real(ConfigurationInteraction_instance%configurations(b)%nDeterminants))
    print *, "prefactor2", prefactor
    prefactor =  1.0_8/sqrt(real(ConfigurationInteraction_instance%configurations(a)%nDeterminants) * &!! (1/n!) 
             real(ConfigurationInteraction_instance%configurations(b)%nDeterminants))

    CIenergy = CIenergy * prefactor
    print *, "E2",CIenergy

  end subroutine ConfigurationInteraction_calculateCIenergyB



  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_getTransformedIntegrals()
    implicit none

!    type(TransformIntegrals) :: repulsionTransformer
    integer :: numberOfSpecies
    integer :: i,j,m,n,mu,nu
    integer :: specieID
    integer :: otherSpecieID
    character(10) :: nameOfSpecie
    character(10) :: nameOfOtherSpecie
    integer :: ocupationNumber
    integer :: ocupationNumberOfOtherSpecie
    integer :: numberOfContractions
    integer :: numberOfContractionsOfOtherSpecie
    type(Matrix) :: auxMatrix
    type(Matrix) :: molecularCouplingMatrix
    type(Matrix) :: molecularExtPotentialMatrix
    type(Matrix) :: couplingMatrix
    type(Matrix) :: hcoreMatrix
    type(Matrix) :: coefficients

    character(50) :: wfnFile
    character(50) :: arguments(20)
    integer :: wfnUnit

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    allocate(ConfigurationInteraction_instance%twoCenterIntegrals(numberOfSpecies))
    allocate(ConfigurationInteraction_instance%FockMatrix(numberOfSpecies))
    allocate(ConfigurationInteraction_instance%fourCenterIntegrals(numberOfSpecies,numberOfSpecies))
    allocate(ConfigurationInteraction_instance%energyofmolecularorbitals(numberOfSpecies))

!    print *,""
!    print *,"BEGIN INTEGRALS TRANFORMATION:"
!    print *,"========================================"
!    print *,""
!    print *,"--------------------------------------------------"
!    print *,"    Algorithm Four-index integral tranformation"
!    print *,"      Yamamoto, Shigeyoshi; Nagashima, Umpei. "
!    print *,"  Computer Physics Communications, 2005, 166, 58-65"
!    print *,"--------------------------------------------------"
!    print *,""
!
!    call TransformIntegrals_constructor( repulsionTransformer )

    do i=1, numberOfSpecies
        nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( i ) )
        specieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecie )
        ocupationNumber = MolecularSystem_getOcupationNumber( i )
        numberOfContractions = MolecularSystem_getTotalNumberOfContractions( i )

!        write (6,"(T10,A)")"ONE PARTICLE INTEGRALS TRANSFORMATION FOR: "//trim(nameOfSpecie)
        call Matrix_constructor (ConfigurationInteraction_instance%twoCenterIntegrals(i), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )

        call Matrix_constructor (ConfigurationInteraction_instance%FockMatrix(i), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )
        call Matrix_constructor (hcoreMatrix,int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)
        call Matrix_constructor (couplingMatrix,int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)

        call Vector_constructor (ConfigurationInteraction_instance%energyofmolecularorbitals(i), &
                int(numberOfContractions,4), 0.0_8)
        !hcoreMatrix  = HartreeFock_instance%HcoreMatrix 

!        hcoreMatrix%values = WaveFunction_HF_instance( specieID )%IndependentParticleMatrix%values

        !! Open file for wavefunction

        wfnFile = "lowdin.wfn"
        wfnUnit = 20

        open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

        arguments(2) = MolecularSystem_getNameOfSpecie(i)
        arguments(1) = "COEFFICIENTS"

        coefficients = &
          Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
          columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

        arguments(1) = "HCORE"

        hcoreMatrix = &
          Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
          columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

        !arguments(1) = "FOCK"
        !ConfigurationInteraction_instance%FockMatrix(i) = &
        !  Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
        !  columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

        arguments(1) = "ORBITALS"
        call Vector_getFromFile( elementsNum = numberOfContractions, &
          unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
          output =ConfigurationInteraction_instance%energyofmolecularorbitals(i) )     

        do m=1,numberOfContractions
           ConfigurationInteraction_instance%fockMatrix(i)%values(m,m) = &
                ConfigurationInteraction_instance%energyofmolecularorbitals(i)%values(m) 
        end do

        !print *, "fock matrix for species", i
        !call matrix_show ( ConfigurationInteraction_instance%fockMatrix(i) )

        do m=1,numberOfContractions
          do n=m, numberOfContractions
             do mu=1, numberOfContractions
                do nu=1, numberOfContractions
                    ConfigurationInteraction_instance%twoCenterIntegrals(i)%values(m,n) = &
                        ConfigurationInteraction_instance%twoCenterIntegrals(i)%values(m,n) + &
                        coefficients%values(mu,m)* &
                        coefficients%values(nu,n)* &
                        hcoreMatrix%values(mu,nu)
                end do
             end do
          end do
       end do

        !print *, "one particle matrix species", i
        !call matrix_show ( ConfigurationInteraction_instance%twoCenterIntegrals(i) )

!! Not implemented yet
!!       if( WaveFunction_HF_instance( specieID )%isThereExternalPotential ) then
!!          do m=1,numberOfContractions
!!             do n=m, numberOfContractions
!!                do mu=1, numberOfContractions
!!                   do nu=1, numberOfContractions
!!                      ConfigurationInteraction_instance%twoCenterIntegrals(i)%values(m,n) = &
!!                           ConfigurationInteraction_instance%twoCenterIntegrals(i)%values(m,n) + &
!!                           WaveFunction_HF_instance( specieID )%waveFunctionCoefficients%values(mu,m)* &
!!                           WaveFunction_HF_instance( specieID )%waveFunctionCoefficients%values(nu,n) * &
!!                           WaveFunction_HF_instance( specieID )%ExternalPotentialMatrix%values(mu,nu)
!!                   end do
!!                end do
!!             end do
!!          end do
!!       end if

       do m=1,numberOfContractions
          do n=m, numberOfContractions
             ConfigurationInteraction_instance%twoCenterIntegrals(i)%values(n,m)=&
                  ConfigurationInteraction_instance%twoCenterIntegrals(i)%values(m,n)
          end do
       end do

        !print *, "Independent Particle"
        !call Matrix_show ( ConfigurationInteraction_instance%twoCenterIntegrals(i) )

!       write (6,"(T10,A)")"TWO PARTICLES INTEGRALS TRANSFORMATION FOR: "//trim(nameOfSpecie)
!       print *,""

!       call TransformIntegrals_atomicToMolecularOfOneSpecie( repulsionTransformer, MolecularSystem_getEigenvectors(specieID), &
!            ConfigurationInteraction_instance%fourCenterIntegrals(i,i), specieID, trim(nameOfSpecie) )

        call ReadTransformedIntegrals_readOneSpecies( specieID, ConfigurationInteraction_instance%fourCenterIntegrals(i,i)   )
        !print *, "two Particle"
        !call Matrix_show(ConfigurationInteraction_instance%fourCenterIntegrals(i,i))

       if ( numberOfSpecies > 1 ) then
          do j = 1 , numberOfSpecies
             if ( i .ne. j) then
                nameOfOtherSpecie= trim(  MolecularSystem_getNameOfSpecie( j ) )
                otherSpecieID =MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecie )
                ocupationNumberOfOtherSpecie = MolecularSystem_getOcupationNumber( j )
                numberOfContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( j )

!                write (6,"(T10,A)") "INTER-SPECIES INTEGRALS TRANSFORMATION FOR: "//trim(nameOfSpecie)//"/"//trim(nameOfOtherSpecie)
!                print *,""

!                call TransformIntegrals_atomicToMolecularOfTwoSpecies( repulsionTransformer, &
!                     MolecularSystem_getEigenVectors(i), MolecularSystem_getEigenVectors(j), &
!                     ConfigurationInteraction_instance%fourCenterIntegrals(i,j), specieID, nameOfSpecie, otherSpecieID, nameOfOtherSpecie )

                 call ReadTransformedIntegrals_readTwoSpecies( specieID, otherSpecieID, &
                         ConfigurationInteraction_instance%fourCenterIntegrals(i,j) )

             end if
          end do
       end if
    end do
!    print *,"END INTEGRALS TRANFORMATION:"
    close (wfnUnit)
  call Matrix_destructor (hcoreMatrix)
        call Matrix_destructor (couplingMatrix)

  end subroutine ConfigurationInteraction_getTransformedIntegrals


  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
!  subroutine ConfigurationInteraction_printTransformedIntegralsToFile()
!    implicit none
!
!!    type(TransformIntegrals) :: repulsionTransformer
!    integer :: numberOfSpecies
!    integer :: i,j,m,n,mu,nu
!    integer :: a,b,r,s,u, auxIndex
!    integer :: z
!    integer :: stats, recNum
!    character(10) :: nameOfSpecie, auxNameOfSpecie
!    character(10) :: nameOfOtherSpecie
!    integer :: ocupationNumber
!    integer :: ocupationNumberOfOtherSpecie
!    integer :: numberOfContractions
!    integer :: numberOfContractionsOfOtherSpecie
!    type(Matrix) :: auxMatrix
!    type(Matrix) :: molecularCouplingMatrix
!    type(Matrix) :: molecularExtPotentialMatrix
!
!    integer :: spin
!
!    real(8) :: totalCoupEnergy
!    real(8) :: fixedPotEnergy
!    real(8) :: fixedIntEnergy
!    real(8) :: KineticEnergy
!    real(8) :: RepulsionEnergy
!    real(8) :: couplingEnergy


!    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
!
!    print *,""
!    print *,"BEGIN INTEGRALS TRANFORMATION:"
!    print *,"========================================"
!    print *,""
!    print *,"--------------------------------------------------"
!    print *,"    Algorithm Four-index integral tranformation"
!    print *,"      Yamamoto, Shigeyoshi; Nagashima, Umpei. "
!    print *,"  Computer Physics Communications, 2005, 166, 58-65"
!    print *,"--------------------------------------------------"
!    print *,""
!
!    totalCoupEnergy = 0.0_8
!    fixedPotEnergy = 0.0_8
!    fixedIntEnergy = 0.0_8
!    KineticEnergy = 0.0_8
!    RepulsionEnergy = 0.0_8
!    couplingEnergy = 0.0_8
!    spin = 0
!
!    do i=1, numberOfSpecies
!        nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( i ) )
!        numberOfContractions = MolecularSystem_getTotalNumberOfContractions( i )
!        spin = MolecularSystem_getMultiplicity(i) - 1
!
!        if(trim(nameOfSpecie) /= "E-BETA" ) then
!
!           if(trim(nameOfSpecie) /= "U-" ) then 
!
!              open(unit=35, file="FCIDUMP-"//trim(nameOfSpecie)//".com", form="formatted", status="replace")
!
!              write(35,"(A)")"gprint basis"
!              write(35,"(A)")"memory 1000 M"
!              write(35,"(A)")"cartesian"
!              write(35,"(A)")"gthresh twoint=1e-12 prefac=1e-14 energy=1e-10 edens=1e-10 zero=1e-12"
!              write(35,"(A)")"basis={"
!              call ConfigurationInteraction_printBasisSetToFile(35)
!              write(35,"(A)")"}"
!
!              write(35,"(A)")"symmetry nosym"
!              write(35,"(A)")"angstrom"
!              write(35,"(A)")"geometry={"
!              call ConfigurationInteraction_printGeometryToFile(35)
!              write(35,"(A)")"}"
!
!              write(35,"(A)")"import 21500.2 "//trim(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//"jcoup")
!              write(35,"(A)")"import 21510.2 "//trim(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//"icoup")
!              write(35,"(A)")"import 21520.2 "//trim(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//"kin")
!              write(35,"(A)")"import 21530.2 "//trim(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//"coeff")
!
!              if(trim(nameOfSpecie) == "E-ALPHA") then
!
!                 write(35,"(A)")"import 21550.2 "//trim(trim(CONTROL_instance%INPUT_FILE)//"E-BETA"//"."//"coeff")
!
!              end if
!
!              write(35,"(A)")"import 21540.2 "//trim(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//"dens")
!              !write(35,"(A)")"import 21560.2 "//trim(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//"pot")
!
!              write(35,"(A)")"{matrop"
!              write(35,"(A)")"load Jcoup, SQUARE 21500.2"
!              write(35,"(A)")"load Icoup, SQUARE 21510.2"
!              write(35,"(A)")"load K, SQUARE 21520.2"
!              !write(35,"(A)")"load Pot, SQUARE 21560.2"
!              write(35,"(A)")"add H01, K Icoup Jcoup"! Pot"
!              write(35,"(A)")"save H01, 21511.2 H0"
!              write(35,"(A)")"}"
!
!              if(trim(nameOfSpecie) == "E-ALPHA") then
!                 write(35,"(A)")"{matrop"
!                 write(35,"(A)")"load Ca, SQUARE 21530.2"
!                 write(35,"(A)")"load Cb, SQUARE 21550.2"               
!                 write(35,"(A)")"save Ca, 2100.1 ORBITALS alpha"
!                 write(35,"(A)")"save Cb, 2100.1 ORBITALS beta"
!                 write(35,"(A)")"}"
!              else
!                 write(35,"(A)")"{matrop"
!                 write(35,"(A)")"load C, SQUARE 21530.2"
!                 write(35,"(A)")"save C, 2100.1 ORBITALS"
!                 write(35,"(A)")"}"
!              end if
!
!              write(35,"(A)")"{matrop"
!              write(35,"(A)")"load D, SQUARE 21540.2"
!              write(35,"(A)")"save D, 21400.1 DENSITY"
!              write(35,"(A)")"}"
!
!
!              !            write(35,"(A,I3,A,I3,A,I3,A1)")"$FCI NORB=",numberOfContractions, ",NELEC=", MolecularSystem_getNumberOfParticles(i)-spin, ", MS2=", spin,","
!              !
!              !            write(35,"(A)",advance="no") "ORBSYM="
!              !            do z=1, numberOfContractions
!              !                write(35,"(I1,A1)",advance="no") 1,","
!              !            end do
!              !            write(35,"(A)") ""
!              !
!              !            write(35, "(A,I3,A,I9)") "ISYM=",1, ",MEMORY=", 200000000
!              !
!              !            write(35, "(A)") "$"
!              !
!              !            print *, "FOUR CENTER INTEGRALS FOR SPECIE: ", trim(nameOfSpecie)
!              !
!              !            recNum = 0
!              !            do a = 1, numberOfContractions
!              !                n = a
!              !                do b=a, numberOfContractions
!              !                    u = b
!              !                    do r = n, numberOfContractions
!              !                        do s = u, numberOfContractions
!              !
!              !                            auxIndex = IndexMap_tensorR4ToVector( a, b, r, s, numberOfContractions )
!              !                            write(35,"(F20.10,4I3)") ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1), a, b, r, s
!              !
!              !                        end do
!              !                        u=r+1
!              !                    end do
!              !                end do
!              !            end do
!              !
!              !
!              !            print *, "TWO CENTER TRANSFORMED INTEGRALS FOR SPECIE: ", trim(nameOfSpecie)
!              !
!              !            do m=1,numberOfContractions
!              !                do n=1, m
!              !                    write(35,"(F20.10,4I3)") ConfigurationInteraction_instance%twoCenterIntegrals(i)%values(m,n), m, n, 0, 0
!              !                end do
!              !            end do
!
!              !!Calculating the core energy....
!
!
!
!              totalCoupEnergy = MolecularSystem_instance%totalCouplingEnergy
!              fixedPotEnergy = MolecularSystem_instance%puntualInteractionEnergy
!
!              do j = 1, numberOfSpecies
!
!                 auxNameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( j ) )
!
!                 if(trim(auxNameOfSpecie) == "E-ALPHA" .or.  trim(auxNameOfSpecie) == "E-BETA" .or.  trim(auxNameOfSpecie) == "e-") cycle
!
!                 fixedIntEnergy = fixedIntEnergy + MolecularSystem_instance%quantumPuntualInteractionEnergy(j)
!                 KineticEnergy = KineticEnergy + MolecularSystem_instance%kineticEnergy(j)
!                 RepulsionEnergy = RepulsionEnergy + MolecularSystem_instance%repulsionEnergy(j)
!                 couplingEnergy = couplingEnergy + MolecularSystem_instance%couplingEnergy(j)
!
!              end do
!
!              !!COMO SEA QUE SE META LA ENERGIA DE CORE
!              !write(35,"(F20.10,4I3)") (couplingEnergy-totalCoupEnergy+fixedPotEnergy+fixedIntEnergy+KineticEnergy+RepulsionEnergy), 0, 0, 0, 0
!              
!              print*, "COREENERGY ", (couplingEnergy-totalCoupEnergy+fixedPotEnergy+fixedIntEnergy+KineticEnergy+RepulsionEnergy)
!
!              write(35,"(A)")"{hf"
!              write(35,"(A)")"maxit 250"
!              write(35,"(A10,I2,A1,A6,I2,A1,A6,I3)")"wf spin=", spin, ",", "charge=",0, ",", "elec=", MolecularSystem_getNumberOfParticles(i)-spin
!              write(35,"(A)")"start 2100.1"
!              write(35,"(A)")"}"
!
!
!              write(35,"(A)")"{fci"
!              write(35,"(A)")"maxit 250"
!              write(35,"(A)")"dm 21400.1, IGNORE_ERROR"
!              write(35,"(A)")"orbit 2100.1, IGNORE_ERROR"
!              write(35,"(A10,I2,A1,A6,I2,A1,A6,I3)")"wf spin=", spin, ",", "charge=",0, ",", "elec=", MolecularSystem_getNumberOfParticles(i)-spin
!              !            write(35,"(A)")"print, orbital=2 integral = 2"
!              !            write(35,"(A)")"CORE"
!              write(35,"(A)")"}"
!
!              write(35,"(A)")"{matrop"
!              write(35,"(A)")"load D, DEN, 21400.1"
!              !	    write(35,"(A)")"print D"
!              write(35,"(A)")"natorb Norb, D"
!              write(35,"(A)")"save Norb, 21570.2"
!              !	    write(35,"(A)")"print Norb"
!              write(35,"(A)")"}"
!
!              write(35,"(A)")"put molden "//trim(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//"molden")//"; orb, 21570.2"
!
!              close(35)
!
!              print*, ""
!
!              stats = system("molpro "//"FCIDUMP-"//trim(nameOfSpecie)//".com ")
!              stats = system("cat "//"FCIDUMP-"//trim(nameOfSpecie)//".out ")
!
!              print*, ""
!
!              print *,"END"
!              
!           end if
!
!         end if
!
!    end do
    
!  end subroutine ConfigurationInteraction_printTransformedIntegralsToFile

  subroutine ConfigurationInteraction_printGeometryToFile(unit)
    implicit none
    integer :: unit

    integer :: i
    integer :: from, to
    real(8) :: origin(3)
    character(50) :: auxString

    
!    do i = 1, MolecularSystem_getTotalNumberOfParticles()
!       
!       origin = MolecularSystem_getOrigin( iterator = i ) * AMSTRONG
!       auxString = trim( MolecularSystem_getNickName( iterator = i ) )
!       
!       if( String_findSubstring( trim( auxString ), "e-") == 1 ) then
!          if( String_findSubstring( trim( auxString ), "BETA") > 1 ) then
!             cycle
!          end if
!            
!          from =String_findSubstring( trim(auxString), "[")
!          to = String_findSubstring( trim(auxString), "]")
!          auxString = auxString(from+1:to-1)
!          
!       else if( String_findSubstring( trim( auxString ), "_") /= 0 ) then
!          cycle
!       end if
!         
!         
!       write (unit,"(A10,3F20.10)") trim( auxString ), origin(1), origin(2), origin(3)
!       
!    end do

  end subroutine ConfigurationInteraction_printGeometryToFile


  subroutine ConfigurationInteraction_printBasisSetToFile(unit)
    implicit none

    integer :: unit

    integer :: i, j
    character(16) :: auxString


!    do i =1, MolecularSystem_instance%numberOfQuantumSpecies
!       
!       auxString=trim( Map_getKey( MolecularSystem_instance%speciesID, iterator=i ) )
!       
!       if( String_findSubstring( trim(auxString), "e-") == 1 ) then
!          
!          if( String_findSubstring( trim(auxString), "BETA") > 1 ) then
!             
!             cycle
!             
!          end if
!          
!          
!       end if
!       
!       if(trim(auxString)=="U-") cycle
!
!       do j =1, size(MolecularSystem_instance%particlesPtr)
!
!          if (    trim(MolecularSystem_instance%particlesPtr(j)%symbol) == trim( Map_getKey( MolecularSystem_instance%speciesID, iterator=i ) ) &
!               .and. MolecularSystem_instance%particlesPtr(j)%isQuantum ) then
!             
!             call BasisSet_showInMolproForm( MolecularSystem_instance%particlesPtr(j)%basis, trim(MolecularSystem_instance%particlesPtr(j)%nickname), unit=unit )
!             
!          end if
!          
!       end do
!       
!    end do
    
  end subroutine ConfigurationInteraction_printBasisSetToFile


  !**
  ! @ Retorna la energia final com correccion Moller-Plesset de orrden dado
  !**
  function ConfigurationInteraction_getTotalEnergy() result(output)
    implicit none
    real(8) :: output

    output = ConfigurationInteraction_instance%totalEnergy

  end function ConfigurationInteraction_getTotalEnergy


  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine ConfigurationInteraction_exception( typeMessage, description, debugDescription)
    implicit none
    integer :: typeMessage
    character(*) :: description
    character(*) :: debugDescription

    type(Exception) :: ex

    call Exception_constructor( ex , typeMessage )
    call Exception_setDebugDescription( ex, debugDescription )
    call Exception_setDescription( ex, description )
    call Exception_show( ex )
    call Exception_destructor( ex )

  end subroutine ConfigurationInteraction_exception

end module ConfigurationInteraction_
