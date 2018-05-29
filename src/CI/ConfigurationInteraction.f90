!******************************************************************************
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
  use InputCI_
  use omp_lib
  use ArpackInterface_
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
     type(ivector8) :: auxIndexCIMatrix
     type(matrix) :: eigenVectors
     type(matrix) :: initialEigenVectors
     type(vector8) :: initialEigenValues
     integer(8) :: numberOfConfigurations
     type(vector) :: numberOfOccupiedOrbitals
     type(vector) :: numberOfOrbitals
     type(vector) ::  numberOfSpatialOrbitals2 
     type(vector8) :: eigenvalues
     type(vector) :: lambda !!Number of particles per orbital, module only works for 1 or 2 particles per orbital
     type(matrix), allocatable :: fourCenterIntegrals(:,:)
     type(matrix), allocatable :: twoCenterIntegrals(:)
     type(matrix), allocatable :: FockMatrix(:)
     type(imatrix), allocatable :: twoIndexArray(:)
     type(imatrix), allocatable :: fourIndexArray(:)
     type(vector), allocatable :: energyofmolecularorbitals(:)
     type(configuration), allocatable :: configurations(:)
     type (Vector8) :: diagonalHamiltonianMatrix
     real(8) :: totalEnergy
     integer, allocatable :: totalNumberOfContractions(:)

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
       ConfigurationInteraction_diagonalize, &
       ConfigurationInteraction_naturalOrbitals, &
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
    integer :: i,j,k,l,m,n,p,q,cc,r,s
    integer(8) :: c
    integer :: ma,mb,mc,md,me,pa,pb,pc,pd,pe
    integer :: isLambdaEqual1,lambda,otherlambda
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
    call Vector_constructor (ConfigurationInteraction_instance%numberOfSpatialOrbitals2, numberOfSpecies)

    if ( allocated ( ConfigurationInteraction_instance%totalNumberOfContractions ) ) &
    deallocate ( ConfigurationInteraction_instance%totalNumberOfContractions ) 
    allocate ( ConfigurationInteraction_instance%totalNumberOfContractions (numberOfSpecies ) )

    do i=1, numberOfSpecies
       !! We are working in spin orbitals not in spatial orbitals!
       ConfigurationInteraction_instance%lambda%values(i) = MolecularSystem_getLambda( i )
       ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)=MolecularSystem_getOcupationNumber( i )* ConfigurationInteraction_instance%lambda%values(i)
       ConfigurationInteraction_instance%numberOfOrbitals%values(i)=MolecularSystem_getTotalNumberOfContractions( i )* ConfigurationInteraction_instance%lambda%values(i)
       ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(i) = MolecularSystem_getTotalNumberOfContractions( i )
       ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(i) = &
         ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(i) *  ( &
         ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(i) + 1 ) / 2

      
       ConfigurationInteraction_instance%totalNumberOfContractions( i ) = MolecularSystem_getTotalNumberOfContractions( i )

      !! Take the active space from input
      if ( InputCI_Instance(i)%activeOrbitals /= 0 ) then
        ConfigurationInteraction_instance%numberOfOrbitals%values(i) = InputCI_Instance(i)%activeOrbitals * &
                                    ConfigurationInteraction_instance%lambda%values(i)

      end if

       !!Uneven occupation number = alpha
       !!Even occupation number = beta     
    end do

    call Configuration_globalConstructor()

    select case ( trim(level) )

    case ( "CIS" )

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

      ConfigurationInteraction_instance%numberOfConfigurations = c
      allocate (ConfigurationInteraction_instance%configurations(ConfigurationInteraction_instance%numberOfConfigurations) )

      call Matrix_Constructor(ConfigurationInteraction_instance%hamiltonianMatrix, int(ConfigurationInteraction_instance%numberOfConfigurations,8),int(ConfigurationInteraction_instance%numberOfConfigurations,8),0.0_8)


    case ( "CISD" )

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

       do i=1, numberOfSpecies
   !!Doubles of the same specie

          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
             do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                do n=m+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   do p= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                      !if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                        do q=p+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                   !do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                      if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda)  ) then !! alpha -> alpha, beta -> beta

                             c=c+1
                      end if
                        end do

                   end do
                end do
             end do
          end if
       end do !! species

       do i=1, numberOfSpecies
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

       ConfigurationInteraction_instance%numberOfConfigurations = c
       allocate (ConfigurationInteraction_instance%configurations(ConfigurationInteraction_instance%numberOfConfigurations) )

       if ( trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) /= "ARPACK" .and. &
        trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) /= "JADAMILU"  ) then

         call Matrix_Constructor(ConfigurationInteraction_instance%hamiltonianMatrix, &
                int(ConfigurationInteraction_instance%numberOfConfigurations,8), & 
                int(ConfigurationInteraction_instance%numberOfConfigurations,8),0.0_8)

       end if

    case ( "CIDD" )

       !!Count
       
       !!Ground State
       c=1
       do i=1, numberOfSpecies
   !!Doubles of the same specie

          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
             do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                do n=m+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   do p= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                      !if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                        do q=p+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                   !do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                      if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda)  ) then !! alpha -> alpha, beta -> beta

                             c=c+1
                      end if
                        end do

                   end do
                end do
             end do
          end if

       end do !! species

       do i=1, numberOfSpecies
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

       ConfigurationInteraction_instance%numberOfConfigurations = c
       allocate (ConfigurationInteraction_instance%configurations(ConfigurationInteraction_instance%numberOfConfigurations) )

       if ( trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) /= "ARPACK" .and. &
        trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) /= "JADAMILU"  ) then


         call Matrix_Constructor(ConfigurationInteraction_instance%hamiltonianMatrix, &
                int(ConfigurationInteraction_instance%numberOfConfigurations,8), & 
                int(ConfigurationInteraction_instance%numberOfConfigurations,8),0.0_8)

       end if

    case ( "FCI" )

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


       do i=1, numberOfSpecies
   !!Doubles of the same specie

          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
             do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                do n=m+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   do p= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                        do q=p+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                      if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) .and. (p /= q) ) then !! alpha -> alpha, beta -> beta

                             c=c+1
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

       ConfigurationInteraction_instance%numberOfConfigurations = c
       allocate (ConfigurationInteraction_instance%configurations(ConfigurationInteraction_instance%numberOfConfigurations) )


       if ( trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) /= "ARPACK" .and. &
        trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) /= "JADAMILU"  ) then

       call Matrix_Constructor(ConfigurationInteraction_instance%hamiltonianMatrix, int(ConfigurationInteraction_instance%numberOfConfigurations,8),int(ConfigurationInteraction_instance%numberOfConfigurations,8),0.0_8)
        end if

    case ( "CISDT" )

      !!Ground State
      c=1

      do i=1, numberOfSpecies

        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        !!Singles (1)
        if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
          do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
            do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
              int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
              if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                c=c+1
              end if
            end do
          end do
        end if
        
       end do

      !!Doubles of the same specie (2)
       do i=1, numberOfSpecies
         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
           do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
             do n=m+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
               do p= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                 do q=p+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                   if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) .and. (p /= q) ) then !! alpha -> alpha, beta -> beta
                     c=c+1
                   end if
                 end do
               end do
             end do
           end do
         end if
       end do !! species

       !!Doubles of different species (11)
       do i=1, numberOfSpecies
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


       !!Triples of the same specie (3)
       do i=1, numberOfSpecies

         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 3 ) then
           do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
             do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
               do mc=mb+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                 do pa= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                         int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                   do pb=pa+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                     do pc=pb+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                       if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                           mod(mb,lambda) == mod(pb,lambda) .and. &
                           mod(mc,lambda) == mod(pc,lambda) ) then !! alpha -> alpha, beta -> beta
                            c=c+1
                       end if
                     end do
                   end do
                 end do
               end do
             end do
           end do
         end if
       end do !! species


       !!Triples (21) 
       do i=1, numberOfSpecies
         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         do j=1, numberOfSpecies
           if ( j /= i ) then
             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
                .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
               do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                 do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                   do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                     do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                            int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                       do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                         do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, &
                           int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                           if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                mod(mb,lambda) == mod(pb,lambda) .and. &
                                mod(mc,lambda) == mod(pc,lambda) )  then !! alpha -> alpha, beta -> beta
                             c=c+1
                           end if
                         end do
                       end do
                     end do  
                   end do
                 end do
               end do
             end if
           end if
         end do

       end do !! species


      !!Triples (111)
       do i=1, numberOfSpecies
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

       ConfigurationInteraction_instance%numberOfConfigurations = c
       allocate (ConfigurationInteraction_instance%configurations(ConfigurationInteraction_instance%numberOfConfigurations) )

       if ( trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) /= "ARPACK" .and. &
        trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) /= "JADAMILU"  ) then

       call Matrix_Constructor(ConfigurationInteraction_instance%hamiltonianMatrix, int(ConfigurationInteraction_instance%numberOfConfigurations,8),int(ConfigurationInteraction_instance%numberOfConfigurations,8),0.0_8)
        end if

    case ( "CISDTQ" )

      !!Ground State
      c=1

      do i=1, numberOfSpecies

        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        !!Singles (1)
        if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
          do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
            do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
              int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
              if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                c=c+1
              end if
            end do
          end do
        end if
        
       end do

      !!Doubles of the same specie (2)
       do i=1, numberOfSpecies
         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
           do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
             do n=m+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
               do p= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                 do q=p+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                   if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) .and. (p /= q) ) then !! alpha -> alpha, beta -> beta
                     c=c+1
                   end if
                 end do
               end do
             end do
           end do
         end if
       end do !! species

       !!Doubles of different species (11)
       do i=1, numberOfSpecies
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


       !!Triples of the same specie (3)
       do i=1, numberOfSpecies

         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 3 ) then
           do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
             do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
               do mc=mb+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                 do pa= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                         int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                   do pb=pa+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                     do pc=pb+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                       if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                           mod(mb,lambda) == mod(pb,lambda) .and. &
                           mod(mc,lambda) == mod(pc,lambda) ) then !! alpha -> alpha, beta -> beta
                            c=c+1
                       end if
                     end do
                   end do
                 end do
               end do
             end do
           end do
         end if
       end do !! species


       !!Triples (21) 
       do i=1, numberOfSpecies
         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         do j=1, numberOfSpecies
           if ( j /= i ) then
             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
                .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
               do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                 do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                   do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                     do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                            int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                       do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                         do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, &
                           int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                           if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                mod(mb,lambda) == mod(pb,lambda) .and. &
                                mod(mc,lambda) == mod(pc,lambda) )  then !! alpha -> alpha, beta -> beta
                             c=c+1
                           end if
                         end do
                       end do
                     end do  
                   end do
                 end do
               end do
             end if
           end if
         end do

       end do !! species


      !!Triples (111)
       do i=1, numberOfSpecies
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

       !!Quadruples of the same species (4)
       do i=1, numberOfSpecies
         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 4 ) then
           do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
             do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
               do mc=mb+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                 do md=mc+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   do pa= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                           int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                     do pb=pa+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                       do pc=pb+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                         do pd=pc+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                           if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                               mod(mb,lambda) == mod(pb,lambda) .and. &
                               mod(mc,lambda) == mod(pc,lambda) .and. &
                               mod(md,lambda) == mod(pd,lambda) ) then !! alpha -> alpha, beta -> beta
                                c=c+1
                           end if
                         end do
                       end do
                     end do
                   end do
                 end do
               end do
             end do
           end do
         end if

       end do !! species

       !!Quadruples (31)
       do i=1, numberOfSpecies

         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         do j=1, numberOfSpecies
           if ( j /= i ) then
             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 3 &
                .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
                do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                  do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                    do mc=mb+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                      do md=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                        do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                               int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                          do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                            do pc=pb+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                              do pd= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, &
                                int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                                if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                     mod(mb,lambda) == mod(pb,lambda) .and. &
                                     mod(mc,lambda) == mod(pc,lambda) .and. &
                                     mod(md,lambda) == mod(pd,lambda) )  then !! alpha -> alpha, beta -> beta
                                  c=c+1
                                end if
                              end do
                            end do
                          end do
                        end do  
                      end do  
                    end do
                  end do
                end do
             end if
           end if
         end do

       end do !! species

       !!Quadruples (22)
       do i=1, numberOfSpecies
         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         do j=i+1, numberOfSpecies
           otherlambda=ConfigurationInteraction_instance%lambda%values(j) !Particles per orbital
           if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
              .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 2 ) then
              do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                  do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                    do md=mc+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                      do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                             int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                        do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                          do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1,&
                             int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                            do pd=pc+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                              if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                   mod(mb,lambda) == mod(pb,lambda) .and. &
                                   mod(mc,otherlambda) == mod(pc,otherlambda) .and. &
                                   mod(md,otherlambda) == mod(pd,otherlambda) )  then !! alpha -> alpha, beta -> beta
                                c=c+1
                              end if
                            end do
                          end do
                        end do
                      end do  
                    end do  
                  end do
                end do
              end do
           end if
         end do

       end do !! species

       !!quadruples of diff specie (1111)
       do i=1, numberOfSpecies
         if ( numberOfSpecies > 3 ) then
          do j=i+1, numberOfSpecies
            do k=j+1, numberOfSpecies
              do l=k+1, numberOfSpecies

                if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. & 
                   ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 .and. &
                   ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 .and. &
                   ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(l) .ge. 1 ) then
                  do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                    do mb=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                      do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) )
                        do md=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(l) )
                          do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                            do pb= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                              do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(k))
                                do pd= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(l))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(l))
                                  c=c+1
                                end do
                              end do
                            end do
                          end do
                        end do
                      end do
                    end do
                  end do

                end if

              end do
            end do
          end do
        end if

       end do !! species

       !!Quadruples (211)
       do i=1, numberOfSpecies
         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         do j=1, numberOfSpecies
           if ( j /= i ) then
             do k=j+1, numberOfSpecies
               if ( k /= j .and. k /=i ) then
                 if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
                    .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 &
                    .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 ) then
                   do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                     do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                       do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                         do md=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) )
                           do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                                  int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                             do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                               do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1,&
                                  int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                                 do pd= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1,&
                                    int(ConfigurationInteraction_instance%numberOfOrbitals%values(k) )
                                   if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                        mod(mb,lambda) == mod(pb,lambda) .and. &
                                        mod(mc,lambda) == mod(pc,lambda) .and. &
                                        mod(md,lambda) == mod(pd,lambda) )  then !! alpha -> alpha, beta -> beta
                                     c=c+1
                                   end if
                                 end do
                               end do
                             end do
                           end do  
                         end do  
                       end do
                     end do
                   end do
                 end if
               end if
             end do
           end if
         end do

      end do !! species

       ConfigurationInteraction_instance%numberOfConfigurations = c
       allocate (ConfigurationInteraction_instance%configurations(ConfigurationInteraction_instance%numberOfConfigurations) )

       if ( trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) /= "ARPACK" .and. &
        trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) /= "JADAMILU"  ) then

       call Matrix_Constructor(ConfigurationInteraction_instance%hamiltonianMatrix, int(ConfigurationInteraction_instance%numberOfConfigurations,8),int(ConfigurationInteraction_instance%numberOfConfigurations,8),0.0_8)
        end if


    case ( "CISDTQQ" )

      !!Ground State
      c=1

      do i=1, numberOfSpecies

        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        !!Singles (1)
        if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
          do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
            do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
              int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
              if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                c=c+1
              end if
            end do
          end do
        end if
        
       end do

      !!Doubles of the same specie (2)
       do i=1, numberOfSpecies
         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
           do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
             do n=m+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
               do p= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                 do q=p+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                   if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) .and. (p /= q) ) then !! alpha -> alpha, beta -> beta
                     c=c+1
                   end if
                 end do
               end do
             end do
           end do
         end if
       end do !! species

       !!Doubles of different species (11)
       do i=1, numberOfSpecies
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


       !!Triples of the same specie (3)
       do i=1, numberOfSpecies

         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 3 ) then
           do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
             do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
               do mc=mb+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                 do pa= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                         int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                   do pb=pa+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                     do pc=pb+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                       if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                           mod(mb,lambda) == mod(pb,lambda) .and. &
                           mod(mc,lambda) == mod(pc,lambda) ) then !! alpha -> alpha, beta -> beta
                            c=c+1
                       end if
                     end do
                   end do
                 end do
               end do
             end do
           end do
         end if
       end do !! species


       !!Triples (21) 
       do i=1, numberOfSpecies
         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         do j=1, numberOfSpecies
           if ( j /= i ) then
             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
                .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
               do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                 do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                   do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                     do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                            int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                       do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                         do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, &
                           int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                           if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                mod(mb,lambda) == mod(pb,lambda) .and. &
                                mod(mc,lambda) == mod(pc,lambda) )  then !! alpha -> alpha, beta -> beta
                             c=c+1
                           end if
                         end do
                       end do
                     end do  
                   end do
                 end do
               end do
             end if
           end if
         end do

       end do !! species


      !!Triples (111)
       do i=1, numberOfSpecies
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

       !!Quadruples of the same species (4)
       do i=1, numberOfSpecies
         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 4 ) then
           do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
             do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
               do mc=mb+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                 do md=mc+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   do pa= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                           int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                     do pb=pa+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                       do pc=pb+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                         do pd=pc+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                           if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                               mod(mb,lambda) == mod(pb,lambda) .and. &
                               mod(mc,lambda) == mod(pc,lambda) .and. &
                               mod(md,lambda) == mod(pd,lambda) ) then !! alpha -> alpha, beta -> beta
                                c=c+1
                           end if
                         end do
                       end do
                     end do
                   end do
                 end do
               end do
             end do
           end do
         end if

       end do !! species

       !!Quadruples (31)
       do i=1, numberOfSpecies

         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         do j=1, numberOfSpecies
           if ( j /= i ) then
             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 3 &
                .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
                do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                  do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                    do mc=mb+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                      do md=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                        do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                               int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                          do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                            do pc=pb+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                              do pd= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, &
                                int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                                if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                     mod(mb,lambda) == mod(pb,lambda) .and. &
                                     mod(mc,lambda) == mod(pc,lambda) .and. &
                                     mod(md,lambda) == mod(pd,lambda) )  then !! alpha -> alpha, beta -> beta
                                  c=c+1
                                end if
                              end do
                            end do
                          end do
                        end do  
                      end do  
                    end do
                  end do
                end do
             end if
           end if
         end do

       end do !! species

       !!Quadruples (22)
       do i=1, numberOfSpecies
         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         do j=i+1, numberOfSpecies
           otherlambda=ConfigurationInteraction_instance%lambda%values(j) !Particles per orbital
           if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
              .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 2 ) then
              do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                  do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                    do md=mc+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                      do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                             int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                        do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                          do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1,&
                             int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                            do pd=pc+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                              if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                   mod(mb,lambda) == mod(pb,lambda) .and. &
                                   mod(mc,otherlambda) == mod(pc,otherlambda) .and. &
                                   mod(md,otherlambda) == mod(pd,otherlambda) )  then !! alpha -> alpha, beta -> beta
                                c=c+1
                              end if
                            end do
                          end do
                        end do
                      end do  
                    end do  
                  end do
                end do
              end do
           end if
         end do

       end do !! species

       !!quadruples of diff specie (1111)
       do i=1, numberOfSpecies
         if ( numberOfSpecies > 3 ) then
          do j=i+1, numberOfSpecies
            do k=j+1, numberOfSpecies
              do l=k+1, numberOfSpecies

                if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. & 
                   ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 .and. &
                   ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 .and. &
                   ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(l) .ge. 1 ) then
                  do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                    do mb=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                      do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) )
                        do md=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(l) )
                          do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                            do pb= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                              do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(k))
                                do pd= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(l))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(l))
                                  c=c+1
                                end do
                              end do
                            end do
                          end do
                        end do
                      end do
                    end do
                  end do

                end if

              end do
            end do
          end do
        end if

       end do !! species

       !!Quadruples (211)
       do i=1, numberOfSpecies
         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         do j=1, numberOfSpecies
           if ( j /= i ) then
             do k=j+1, numberOfSpecies
               if ( k /= j .and. k /=i ) then
                 if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
                    .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 &
                    .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 ) then
                   do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                     do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                       do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                         do md=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) )
                           do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                                  int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                             do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                               do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1,&
                                  int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                                 do pd= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1,&
                                    int(ConfigurationInteraction_instance%numberOfOrbitals%values(k) )
                                   if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                        mod(mb,lambda) == mod(pb,lambda) .and. &
                                        mod(mc,lambda) == mod(pc,lambda) .and. &
                                        mod(md,lambda) == mod(pd,lambda) )  then !! alpha -> alpha, beta -> beta
                                     c=c+1
                                   end if
                                 end do
                               end do
                             end do
                           end do  
                         end do  
                       end do
                     end do
                   end do
                 end if
               end if
             end do
           end if
         end do

      end do !! species

       !!Quintuplesw (221)
       do i=1, numberOfSpecies
         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         do j=i+1, numberOfSpecies
           otherlambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
           do k=j+1, numberOfSpecies
             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
               .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 2 &
               .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 ) then
               do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                 do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                   do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                     do md=mc+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                       do me=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) )
                         do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                                int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                           do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                             do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1,&
                                int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                               do pd=pc+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                                 do pe= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1,&
                                    int(ConfigurationInteraction_instance%numberOfOrbitals%values(k) )
                                   if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                        mod(mb,lambda) == mod(pb,lambda) .and. &
                                        mod(mc,otherlambda) == mod(pc,otherlambda) .and. &
                                        mod(md,otherlambda) == mod(pd,otherlambda) .and. &
                                        mod(me,lambda) == mod(pe,lambda) )  then !! alpha -> alpha, beta -> beta
                                     c=c+1
                                   end if
                                 end do
                               end do
                             end do
                           end do
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

       ConfigurationInteraction_instance%numberOfConfigurations = c
       allocate (ConfigurationInteraction_instance%configurations(ConfigurationInteraction_instance%numberOfConfigurations) )

       if ( trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) /= "ARPACK" .and. &
        trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) /= "JADAMILU"  ) then

       call Matrix_Constructor(ConfigurationInteraction_instance%hamiltonianMatrix, int(ConfigurationInteraction_instance%numberOfConfigurations,8),int(ConfigurationInteraction_instance%numberOfConfigurations,8),0.0_8)
        end if
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

    !!Destroy configurations
    !!Ground State
    if (allocated(ConfigurationInteraction_instance%configurations)) then
      c=1
      call Configuration_destructor(ConfigurationInteraction_instance%configurations(c) )
  
      do c=2, ConfigurationInteraction_instance%numberOfConfigurations
         call Configuration_destructor(ConfigurationInteraction_instance%configurations(c) )                
      end do
  
      if (allocated(ConfigurationInteraction_instance%configurations)) deallocate(ConfigurationInteraction_instance%configurations)
    end if

    call Matrix_destructor(ConfigurationInteraction_instance%hamiltonianMatrix)
    call Vector_destructor (ConfigurationInteraction_instance%numberOfOccupiedOrbitals)
    call Vector_destructor (ConfigurationInteraction_instance%numberOfOrbitals)
    call Vector_destructor (ConfigurationInteraction_instance%lambda)

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
    real(8) :: davidsonCorrection, HFcoefficient, CIcorrection
    integer numberOfSpecies

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

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
       else 

         print *, ""
         HFcoefficient = ConfigurationInteraction_instance%eigenVectors%values(1,1) 
         write (6,"(T8,A19, F20.12)") "HF COEFFICIENT = ", HFcoefficient

       end if

    else 

    end if

  end subroutine ConfigurationInteraction_show


  !FELIX IS HERE
  subroutine ConfigurationInteraction_naturalOrbitals()
    implicit none
    type(ConfigurationInteraction) :: this
    type(Configuration) :: auxthisA, auxthisB
    integer :: i, j, k, mu, nu
    integer :: factor
    integer :: unit, wfnunit
    integer :: numberOfOrbitals, numberOfOccupiedOrbitals
    integer :: state, specie, orbital, orbitalA, orbitalB
    character(50) :: file, wfnfile, speciesName, auxstring
    character(50) :: arguments(2)
    real(8) :: sumaPrueba
    type(matrix) :: coefficients, densityMatrix
    type(matrix) :: ciOccupationNumbers, ciOccupationMatrix
    integer numberOfSpecies
    type(Vector) :: eigenValues
    type(Matrix) :: eigenVectors, auxMatrix

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    
    if ( ConfigurationInteraction_instance%isInstanced .and. CONTROL_instance%CI_STATES_TO_PRINT .gt. 0 ) then

       print *,""
       print *," FRACTIONAL ORBITAL OCCUPATIONS"
       print *,"=============================="
       print *,"column: state, row: orbital"
       print *,""

       !! Open file - to print natural orbitals
       unit = 29

       file = trim(CONTROL_instance%INPUT_FILE)//"CIOccupations.occ"
       open(unit = unit, file=trim(file), status="new", form="formatted")
       
       ! call Matrix_show (ConfigurationInteraction_instance%eigenVectors)
       
       do specie=1, numberOfSpecies
          
          speciesName = MolecularSystem_getNameOfSpecie(specie)
          numberOfOrbitals = ConfigurationInteraction_instance%numberOfOrbitals%values(specie)
          numberOfOccupiedOrbitals = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(specie)

          wfnFile = "lowdin.wfn"
          wfnUnit = 20
          open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
          arguments(2) = speciesName
          arguments(1) = "COEFFICIENTS"
          coefficients = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfOrbitals,4), &
               columns= int(numberOfOrbitals,4), binary=.true., arguments=arguments(1:2))
          close(wfnUnit)

          do state=1, CONTROL_instance%CI_STATES_TO_PRINT
             !Inicializando la matriz

             call Matrix_constructor ( densityMatrix , &
                  int(numberOfOrbitals,8), &
                  int(numberOfOrbitals,8),  0.0_8 )

             !!Diagonal ground state
             do i=1, ConfigurationInteraction_instance%numberOfConfigurations
                do j=1, numberOfOccupiedOrbitals
                   ! !! Occupied orbitals
                   ! ciOccupationMatrix%values( j, j)= ciOccupationMatrix%values( j, j) -  &
                   !      ConfigurationInteraction_instance%eigenVectors%values(i,1)**2
                   ! !! Unoccupied orbitals
                   orbital = ConfigurationInteraction_instance%configurations(i)%occupations(j,specie) 

                   ! ciOccupationMatrix%values( orbital, orbital)= ciOccupationMatrix%values( orbital, orbital) + &
                   !      ConfigurationInteraction_instance%eigenVectors%values(i,1)**2

                   ! print *, i, i, orbital, orbital, ConfigurationInteraction_instance%eigenVectors%values(i,1)**2

                   do mu = 1 , numberOfOrbitals
                      do nu = 1 , numberOfOrbitals

                         densityMatrix%values(mu,nu) =  &
                              densityMatrix%values(mu,nu) + &
                              ConfigurationInteraction_instance%eigenVectors%values(i,state)**2 *&
                              coefficients%values(mu,orbital)*coefficients%values(nu,orbital)
                      end do
                   end do
                end do
             end do

             !!off-Diagonal ground state
             do i=1, ConfigurationInteraction_instance%numberOfConfigurations-1
                do j=i+1, ConfigurationInteraction_instance%numberOfConfigurations

                   if (Configuration_checkCoincidenceB(ConfigurationInteraction_instance%configurations(i)%occupations,&
                        ConfigurationInteraction_instance%configurations(j)%occupations, numberOfSpecies) .eq. 1) then

                      auxthisA%occupations = ConfigurationInteraction_instance%configurations(i)%occupations
                      auxthisB%occupations = ConfigurationInteraction_instance%configurations(j)%occupations
                      factor = 1
                      call Configuration_setAtMaximumCoincidenceB( auxthisA%occupations,auxthisB%occupations, numberOfSpecies, factor )

                      ! print *, i, j, ConfigurationInteraction_instance%configurations(i)%occupations(:,specie), ConfigurationInteraction_instance%configurations(j)%occupations(:,specie)
                      ! print *, i, j, auxthisA%occupations(:,specie), auxthisB%occupations(:,specie)

                      do k=1, numberOfOccupiedOrbitals

                         if(auxthisA%occupations(k,specie) .ne. auxthisB%occupations(k,specie)) then

                            orbitalA = auxthisA%occupations(k,specie)
                            orbitalB = auxthisB%occupations(k,specie)

                            ! print *, "trololooooo"
                            ! print *, i, j, orbitalA, orbitalB, factor*ConfigurationInteraction_instance%eigenVectors%values(i,1)*ConfigurationInteraction_instance%eigenVectors%values(j,1)

                            ! ciOccupationMatrix%values( orbitalA, orbitalB)= ciOccupationMatrix%values( orbitalA, orbitalB) + &
                            !      factor*ConfigurationInteraction_instance%eigenVectors%values(i,1)*ConfigurationInteraction_instance%eigenVectors%values(j,1)

                            ! ciOccupationMatrix%values( orbitalB, orbitalA)= ciOccupationMatrix%values( orbitalB, orbitalA) + &
                            !      factor*ConfigurationInteraction_instance%eigenVectors%values(i,1)*ConfigurationInteraction_instance%eigenVectors%values(j,1)

                            do mu = 1 , numberOfOrbitals
                               do nu = 1 , numberOfOrbitals

                                  densityMatrix%values(mu,nu) =  &
                                       densityMatrix%values(mu,nu) + &
                                       factor *&
                                       ConfigurationInteraction_instance%eigenVectors%values(i,state) *&
                                       ConfigurationInteraction_instance%eigenVectors%values(j,state) *&
                                       (coefficients%values(mu,orbitalA)*coefficients%values(nu,orbitalB) + coefficients%values(mu,orbitalB)*coefficients%values(nu,orbitalA))
                               end do
                            end do

                         end if
                      end do

                   end if

                end do
             end do

             print *, "density CI", speciesName
             call Matrix_show ( densityMatrix )

             write(auxstring,*) state
             arguments(2) = speciesName
             arguments(1) = "DENSITYMATRIX"//trim(adjustl(auxstring)) 

             call Matrix_writeToFile ( densityMatrix, unit , arguments=arguments(1:2) )

          end do
          
       end do
       close(unit)

    end if

          ! call Vector_constructor(eigenValues, numberOfOrbitals)
          ! call Matrix_constructor(eigenVectors, int(numberOfOrbitals,8), int(numberOfOrbitals,8))
          ! call Matrix_eigen(ciOccupationMatrix, eigenValues, eigenVectors, SYMMETRIC)
          
          ! print *, "Diagonal sum", sum(eigenValues%values)
          ! call Vector_show(eigenValues)

          ! call Matrix_show(eigenVectors)
          ! print *, arguments(1:2)
          ! call Matrix_show ( densityMatrix )

          ! call Matrix_constructor ( ciOccupationNumbers , int(numberOfOrbitals,8) , &
          !      int(CONTROL_instance%CI_STATES_TO_PRINT,8),  0.0_8 )
          
          ! do state=1, CONTROL_instance%CI_STATES_TO_PRINT
          !    sumaPrueba=0
          !    do j=1, numberOfOccupiedOrbitals
          !       ciOccupationNumbers%values(j,state) = 1.0
          !    end do
          
          ! ! !Get occupation numbers from each configuration contribution
             
          !    do i=1, ConfigurationInteraction_instance%numberOfConfigurations
          !       do j=1, numberOfOccupiedOrbitals

          !          !! Occupied orbitals
          !          ciOccupationNumbers%values( j, state)= ciOccupationNumbers%values( j, state) -  &
          !               ConfigurationInteraction_instance%eigenVectors%values(i,state)**2
          !          !! Unoccupied orbitals
          !          orbital = ConfigurationInteraction_instance%configurations(i)%occupations(j,specie) 

          !          ciOccupationNumbers%values( orbital, state)= ciOccupationNumbers%values( orbital, state) + &
          !               ConfigurationInteraction_instance%eigenVectors%values(i,state)**2

          !          ! print *, j, orbital, ConfigurationInteraction_instance%eigenVectors%values(i,state)**2
          !          ! sumaPrueba=sumaPrueba+ConfigurationInteraction_instance%eigenVectors%values(i,state)**2
          !       end do
          !       ! end if

          !    end do

          !    ! print *, "suma", sumaPrueba
          !    !Build a new density matrix (P) in atomic orbitals

          !    call Matrix_constructor ( densityMatrix , &
          !         int(numberOfOrbitals,8), &
          !         int(numberOfOrbitals,8),  0.0_8 )
             
          !    wfnFile = "lowdin.wfn"
          !    wfnUnit = 20

          !    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

          !    arguments(2) = speciesName
          !    arguments(1) = "COEFFICIENTS"

          !    coefficients = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfOrbitals,4), &
          !         columns= int(numberOfOrbitals,4), binary=.true., arguments=arguments(1:2))
             
          !    close(wfnUnit)
             
          !    do mu = 1 , numberOfOrbitals
          !       do nu = 1 , numberOfOrbitals
          !          do k = 1 , numberOfOrbitals

          !             densityMatrix%values(mu,nu) =  &
          !                  densityMatrix%values(mu,nu) + &
          !                  ciOccupationNumbers%values(k, state)**2* &
          !                  coefficients%values(mu,k)*coefficients%values(nu,k)
          !           end do
          !        end do
          !     end do

          !     write(auxstring,*) state
          !     arguments(2) = speciesName
          !     arguments(1) = "DENSITYMATRIX"//trim(adjustl(auxstring)) 
              
          !     call Matrix_writeToFile ( densityMatrix, unit , arguments=arguments(1:2) )

          !     print *, arguments(1:2)
          !     call Matrix_show ( densityMatrix )

          !     call Matrix_destructor(coefficients)          
          !     call Matrix_destructor(densityMatrix)          


          !  end do

          ! !Write occupation numbers to file
          ! write (6,"(T8,A10,A20)") trim(MolecularSystem_getNameOfSpecie(specie)),"OCCUPATIONS:"

          ! call Matrix_show ( ciOccupationNumbers )

          ! arguments(2) = speciesName
          ! arguments(1) = "OCCUPATIONS"

          ! call Matrix_writeToFile ( ciOccupationNumbers, unit , arguments=arguments(1:2) )

          ! call Matrix_destructor(ciOccupationNumbers)          



  end subroutine ConfigurationInteraction_naturalOrbitals



  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_run()
    implicit none 
    integer :: m
    real(8), allocatable :: eigenValues(:) 

!    select case ( trim(ConfigurationInteraction_instance%level) )

       print *, ""
       print *, "==============================================="
       print *, "|         BEGIN ", trim(ConfigurationInteraction_instance%level)," CALCULATION"
       print *, "-----------------------------------------------"
       print *, ""

       print *, "Getting transformed integrals..."
       call ConfigurationInteraction_getTransformedIntegrals()
       print *, "Building configurations..."

       call ConfigurationInteraction_buildConfigurations()
       print *, "Total number of configurations", ConfigurationInteraction_instance%numberOfConfigurations
       print *, ""
       call Vector_constructor8 ( ConfigurationInteraction_instance%eigenvalues, &
                                 ConfigurationInteraction_instance%numberOfConfigurations, 0.0_8)

       select case (trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)))

       case ("ARPACK")

         print *, "Building initial hamiltonian..."
         call ConfigurationInteraction_buildInitialCIMatrix()

         print *, "Building and saving hamiltonian..."
         call ConfigurationInteraction_buildAndSaveCIMatrix()

         !! deallocate transformed integrals
         ! deallocate (ConfigurationInteraction_instance%configurations)
         deallocate(ConfigurationInteraction_instance%twoCenterIntegrals)
         deallocate(ConfigurationInteraction_instance%fourCenterIntegrals)

         call Matrix_constructor (ConfigurationInteraction_instance%eigenVectors, &
              int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
              int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

         print *, ""
         print *, "Diagonalizing hamiltonian..."
         print *, "  Using : ", trim(String_getUppercase((CONTROL_instance%CI_DIAGONALIZATION_METHOD)))
         
         call ConfigurationInteraction_diagonalize(ConfigurationInteraction_instance%numberOfConfigurations, &
              ConfigurationInteraction_instance%numberOfConfigurations, &
              CONTROL_instance%NUMBER_OF_CI_STATES, &
              CONTROL_instance%CI_MAX_NCV, &
              ConfigurationInteraction_instance%eigenvalues, &
              ConfigurationInteraction_instance%eigenVectors )

         case ("JADAMILU")

         print *, "Building initial hamiltonian..."
         call ConfigurationInteraction_buildInitialCIMatrix()

         call Matrix_constructor (ConfigurationInteraction_instance%eigenVectors, &
              int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
              int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

         if ( CONTROL_instance%CI_LOAD_EIGENVECTOR ) then 
           call ConfigurationInteraction_loadEigenVector (ConfigurationInteraction_instance%eigenvalues, &
                  ConfigurationInteraction_instance%eigenVectors) 
         end if 

         print *, ""
         print *, "Diagonalizing hamiltonian..."
         print *, "  Using : ", trim(String_getUppercase((CONTROL_instance%CI_DIAGONALIZATION_METHOD)))
         print *, "============================================================="
         print *, "M. BOLLHÖFER AND Y. NOTAY, JADAMILU:"
         print *, " a software code for computing selected eigenvalues of "
         print *, " large sparse symmetric matrices, "
         print *, "Computer Physics Communications, vol. 177, pp. 951-964, 2007." 
         print *, "============================================================="


         call ConfigurationInteraction_jadamiluInterface(ConfigurationInteraction_instance%numberOfConfigurations, &
              int(CONTROL_instance%NUMBER_OF_CI_STATES,8), &
              ConfigurationInteraction_instance%eigenvalues, &
              ConfigurationInteraction_instance%eigenVectors )

         if ( CONTROL_instance%CI_SAVE_EIGENVECTOR ) then 
           call ConfigurationInteraction_saveEigenVector () 
         end if
       case ("DSYEVX")

         call ConfigurationInteraction_buildHamiltonianMatrix()
         print *, "Reference Energy", ConfigurationInteraction_instance%hamiltonianMatrix%values(1,1)

         call Matrix_constructor (ConfigurationInteraction_instance%eigenVectors, &
              int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
              int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

         !! deallocate transformed integrals
         deallocate(ConfigurationInteraction_instance%twoCenterIntegrals)
         deallocate(ConfigurationInteraction_instance%fourCenterIntegrals)


         print *, ""
         print *, "Diagonalizing hamiltonian..."
         print *, "  Using : ", trim(String_getUppercase((CONTROL_instance%CI_DIAGONALIZATION_METHOD)))


         call Matrix_eigen_select (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
              int(1), int(CONTROL_instance%NUMBER_OF_CI_STATES), &  
              eigenVectors = ConfigurationInteraction_instance%eigenVectors, &
              flags = int(SYMMETRIC,4))

!         call Matrix_eigen_select (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
!              1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
!              flags = SYMMETRIC, dm = ConfigurationInteraction_instance%numberOfConfigurations )


       case ("DSYEVR")

         call ConfigurationInteraction_buildHamiltonianMatrix()
         print *, "Reference Energy", ConfigurationInteraction_instance%hamiltonianMatrix%values(1,1)

         call Matrix_constructor (ConfigurationInteraction_instance%eigenVectors, &
              int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
              int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

         !! deallocate transformed integrals
         deallocate(ConfigurationInteraction_instance%twoCenterIntegrals)
         deallocate(ConfigurationInteraction_instance%fourCenterIntegrals)

         print *, ""
         print *, "Diagonalizing hamiltonian..."
         print *, "  Using : ", trim(String_getUppercase((CONTROL_instance%CI_DIAGONALIZATION_METHOD)))


         call Matrix_eigen_dsyevr (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
              1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
              eigenVectors = ConfigurationInteraction_instance%eigenVectors, &
              flags = SYMMETRIC)

!        call Matrix_eigen_dsyevr (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
!              1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
!              flags = SYMMETRIC, dm = ConfigurationInteraction_instance%numberOfConfigurations )

       case default

         call ConfigurationInteraction_exception( ERROR, "Configuration interactor constructor", "Diagonalization method not implemented")


       end select

       print *,""
       print *, "-----------------------------------------------"
       print *, "|     END ", trim(ConfigurationInteraction_instance%level)," CALCULATION"
       print *, "==============================================="
       print *, ""

         
!    case ( "FCI-oneSpecie" )
!
!       print *, ""
!       print *, ""
!       print *, "==============================================="
!       print *, "|  Full CI for one specie calculation          |"
!       print *, "|  Use fci program to perform the calculation  |"
!       print *, "-----------------------------------------------"
!       print *, ""
!       ! call ConfigurationInteraction_getTransformedIntegrals()
!       !call ConfigurationInteraction_printTransformedIntegralsToFile()
!

  end subroutine ConfigurationInteraction_run

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_buildConfigurations()
    implicit none

    integer :: numberOfSpecies
    integer :: i,ii,j,k,l,m,n,p,q,a,b,d,r,s
    integer(8) :: c, cc
    integer :: ma,mb,mc,md,me,pa,pb,pc,pd,pe
    integer :: isLambdaEqual1
    type(ivector) :: order
    type(vector), allocatable :: occupiedCode(:)
    type(vector), allocatable :: unoccupiedCode(:)
    logical :: sameConfiguration
    integer :: nEquivalentConfigurations, newNumberOfConfigurations
    integer, allocatable :: equivalentConfigurations (:,:), auxArray(:,:), auxvector(:),auxvectorA(:)
    integer :: lambda, otherlambda

    nEquivalentConfigurations = 0
    if (allocated ( equivalentConfigurations )) deallocate ( equivalentConfigurations)
    allocate( equivalentConfigurations(nEquivalentConfigurations,2) )
    equivalentConfigurations = 0

    newNumberOfConfigurations = ConfigurationInteraction_instance%numberOfConfigurations 

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    if ( allocated( occupiedCode ) ) deallocate( occupiedCode )
    allocate (occupiedCode ( numberOfSpecies ) )
    if ( allocated( unoccupiedCode ) ) deallocate( unoccupiedCode )
    allocate (unoccupiedCode ( numberOfSpecies ) )

    select case ( trim(ConfigurationInteraction_instance%level) )

    case ( "CIS" )

   !!Build configurations
       !!Ground State
       c=1
       call Vector_constructorInteger (order, numberOfSpecies, 0 )

       do i=1, numberOfSpecies
         call Vector_constructor (occupiedCode(i), 1, 0.0_8)
         call Vector_constructor (unoccupiedCode(i), 1, 0.0_8)
       end do

       call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations,c) 

       do i=1, numberOfSpecies
          cc = 0
          call Vector_constructorInteger (order, numberOfSpecies, 0 )
          order%values(i)=1
          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital

          !!Singles
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
             do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                !call Vector_constructor (occupiedCode(i), numberOfSpecies, 0.0_8)
                call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                occupiedCode(i)%values(1)=m
                do p=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i))
                   if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                     call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8)
                     unoccupiedCode(i)%values(1)=p
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
        !   call Configuration_sho(w  ConfigurationInteraction_instance%configurations(c) )
        !end do 

!!! Working!! but checktwoConfigurations was commentted
!!       !! Search for equivalent configurations
!!       do a=1, ConfigurationInteraction_instance%numberOfConfigurations
!!          do b=a+1, ConfigurationInteraction_instance%numberOfConfigurations
!!
!!             !call Configuration_checkTwoConfigurations(ConfigurationInteraction_instance%configurations(a), &
!!             !       ConfigurationInteraction_instance%configurations(b), sameConfiguration, numberOfSpecies)
!!
!!             if ( sameConfiguration .eqv. .True. ) then
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
!!       !! Remove equivalent configurations
!!        do c = nEquivalentConfigurations, 1, -1
!!
!!                call Configuration_destructor(ConfigurationInteraction_instance%configurations( &
!!                     equivalentConfigurations(c,2) ) )                
!!        end do
!!
!!        !! Compact the configuration vector... 
!!        do c = 1,  ConfigurationInteraction_instance%numberOfConfigurations 
!!
!!          if ( c <= newNumberOfConfigurations ) then
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
!!         end if
!!
!!        end do 
!!
!!        ConfigurationInteraction_instance%numberOfConfigurations = newNumberOfConfigurations 
        ConfigurationInteraction_instance%numberOfConfigurations = c

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
       call Vector_constructorInteger (order, numberOfSpecies, 0 )

       do i=1, numberOfSpecies
         call Vector_constructor (occupiedCode(i), 1, 0.0_8)
         call Vector_constructor (unoccupiedCode(i), 1, 0.0_8)
       end do

       call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations, c) 

       do i=1, numberOfSpecies

          call Vector_constructorInteger (order, numberOfSpecies, 0 )
          order%values(i)=1
          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          cc = 0
          !!Singles
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
             do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                call Vector_constructor (occupiedCode(i), int( order%values(i),4), 0.0_8)
                occupiedCode(i)%values(1)=m
                do p=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1,int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                   if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                     call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                     unoccupiedCode(i)%values(1)=p
                     c=c+1
                     cc=cc+1
                     call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                          ConfigurationInteraction_instance%numberOfConfigurations,c) 

                   end if

                end do !! p
             end do !! m
          end if

          write (6, "(T4,A20,A21,A3,I8)") "Singles for species ", trim(MolecularSystem_getNameOfSpecie(i)), " : ", cc

        end do !! species

       do i=1, numberOfSpecies
          !!Doubles of the same specie
          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          call Vector_constructorInteger (order, numberOfSpecies, 0 )
          order%values(i)=2
          cc = 0
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
             do m=1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                do n=m+1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8 )
                   !occupiedCode%values(i,1)=m*1024+n
                   occupiedCode(i)%values(1)=m
                   occupiedCode(i)%values(2)=n
                   do p = int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                     do q=p+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                   !do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )

                      if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) ) then !! alpha -> alpha, beta -> beta
                            call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                            !unoccupiedCode%values(i)=p*1024+q
                            unoccupiedCode(i)%values(1)=p
                            unoccupiedCode(i)%values(2)=q
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
          write (6, "(T4,A20,A21,A3,I8)") "Doubles for species ", trim(MolecularSystem_getNameOfSpecie(i)), " : ", cc

          !!Doubles of different species
          do j=i+1, numberOfSpecies

             call Vector_constructorInteger (order, numberOfSpecies, 0 )
             order%values(i)=1
             order%values(j)=1
             cc = 0

             do ii=1, numberOfSpecies
               call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
               call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
             end do


             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
                do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   do n=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))

                      call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                      call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                      occupiedCode(i)%values(1)=m
                      occupiedCode(j)%values(1)=n
                      do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i))
                         do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                            call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                            call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                            unoccupiedCode(i)%values(1)=p
                            unoccupiedCode(j)%values(1)=q
                            c=c+1
                            cc=cc+1
                            call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations, c) 
                         end do
                      end do
                   end do
                end do
             end if

            write (6, "(T4,A20,A21,A3,I8)") "Doubles for species ", trim(MolecularSystem_getNameOfSpecie(i))//&
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
       call Vector_constructorInteger (order, numberOfSpecies, 0 )

       do i=1, numberOfSpecies
         call Vector_constructor (occupiedCode(i), 1, 0.0_8)
         call Vector_constructor (unoccupiedCode(i), 1, 0.0_8)
       end do

       call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations, c) 

       do i=1, numberOfSpecies
          !!Doubles of the same specie
          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          call Vector_constructorInteger (order, numberOfSpecies, 0 )
          order%values(i)=2
          cc = 0
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
             do m=1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                do n=m+1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8 )
                   !occupiedCode%values(i,1)=m*1024+n
                   occupiedCode(i)%values(1)=m
                   occupiedCode(i)%values(2)=n
                   do p = int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                     do q=p+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                   !do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )

                      if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) ) then !! alpha -> alpha, beta -> beta
                            call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                            !unoccupiedCode%values(i)=p*1024+q
                            unoccupiedCode(i)%values(1)=p
                            unoccupiedCode(i)%values(2)=q
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
          write (6, "(T4,A20,A21,A3,I8)") "Doubles for species ", trim(MolecularSystem_getNameOfSpecie(i)), " : ", cc

          !!Doubles of different species
          do j=i+1, numberOfSpecies

             call Vector_constructorInteger (order, numberOfSpecies, 0 )
             order%values(i)=1
             order%values(j)=1
             cc = 0

             do ii=1, numberOfSpecies
               call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
               call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
             end do


             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
                do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   do n=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))

                      call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                      call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                      occupiedCode(i)%values(1)=m
                      occupiedCode(j)%values(1)=n
                      do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i))
                         do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                            call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                            call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                            unoccupiedCode(i)%values(1)=p
                            unoccupiedCode(j)%values(1)=q
                            c=c+1
                            cc=cc+1
                            call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations, c) 
                         end do
                      end do
                   end do
                end do
             end if

            write (6, "(T4,A20,A21,A3,I8)") "Doubles for species ", trim(MolecularSystem_getNameOfSpecie(i))//&
                                                   "/"//trim(MolecularSystem_getNameOfSpecie(j)), " : ", cc
          end do
       end do

    case ( "FCI" )

   !!Build configurations
       !!Ground State
       c=1
       call Vector_constructorInteger (order, numberOfSpecies, 0 )

       do i=1, numberOfSpecies
         call Vector_constructor (occupiedCode(i), 1, 0.0_8)
         call Vector_constructor (unoccupiedCode(i), 1, 0.0_8)
       end do

       call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations, c) 

       do i=1, numberOfSpecies

          call Vector_constructorInteger (order, numberOfSpecies, 0 )
          order%values(i)=1
          lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
          cc = 0
          !!Singles
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
             do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                call Vector_constructor (occupiedCode(i), int( order%values(i),4), 0.0_8)
                occupiedCode(i)%values(1)=m
                do p=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1,int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                   if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                     call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                     unoccupiedCode(i)%values(1)=p
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
          call Vector_constructorInteger (order, numberOfSpecies, 0 )
          order%values(i)=2
          cc = 0
          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
             do m=1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                do n=m+1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8 )
                   !occupiedCode%values(i,1)=m*1024+n
                   occupiedCode(i)%values(1)=m
                   occupiedCode(i)%values(2)=n
                   do p = int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                     do q=p+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                   !do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )

                      if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) ) then !! alpha -> alpha, beta -> beta
                            call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                            !unoccupiedCode%values(i)=p*1024+q
                            unoccupiedCode(i)%values(1)=p
                            unoccupiedCode(i)%values(2)=q
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

             call Vector_constructorInteger (order, numberOfSpecies, 0 )
             order%values(i)=1
             order%values(j)=1
             cc = 0

             do ii=1, numberOfSpecies
               call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
               call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
             end do


             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
                do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   do n=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))

                      call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                      call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                      occupiedCode(i)%values(1)=m
                      occupiedCode(j)%values(1)=n
                      do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i))
                         do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                            call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                            call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                            unoccupiedCode(i)%values(1)=p
                            unoccupiedCode(j)%values(1)=q
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

             call Vector_constructorInteger (order, numberOfSpecies, 0)
             order%values(i)=1
             order%values(j)=1
             order%values(k)=1
             cc = 0

             do ii=1, numberOfSpecies
               call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
               call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
             end do


             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. &
                 ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 .and. &
                 ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 ) then
                do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                   do n=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))
                     do r=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))

                      call Vector_constructor (occupiedCode(i), order%values(i), 0.0_8)
                      call Vector_constructor (occupiedCode(j), order%values(j), 0.0_8)
                      call Vector_constructor (occupiedCode(k), order%values(k), 0.0_8)
                      occupiedCode(i)%values(1)=m
                      occupiedCode(j)%values(1)=n
                      occupiedCode(k)%values(1)=r

                      do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i))
                         do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )

                         do s=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(k) )
                            call Vector_constructor (unoccupiedCode(i), order%values(i), 0.0_8)
                            call Vector_constructor (unoccupiedCode(j), order%values(j), 0.0_8)
                            call Vector_constructor (unoccupiedCode(k), order%values(k), 0.0_8)
                            unoccupiedCode(i)%values(1)=p
                            unoccupiedCode(j)%values(1)=q
                            unoccupiedCode(k)%values(1)=s
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

            write (6, "(T4,A20,A21,A3,I6)") "Triples for species ", trim(MolecularSystem_getNameOfSpecie(i))//&
                                                   "/"//trim(MolecularSystem_getNameOfSpecie(j))//& 
                                                   "/"//trim(MolecularSystem_getNameOfSpecie(k)), " : ", cc
            end do
          end do


       end do

    case ("CISDT")
    !!Build configurations
    !!Ground State
      c=1
      call Vector_constructorInteger (order, numberOfSpecies, 0 )

      do i=1, numberOfSpecies
        call Vector_constructor (occupiedCode(i), 1, 0.0_8)
        call Vector_constructor (unoccupiedCode(i), 1, 0.0_8)
      end do

      call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations, c) 

      do i=1, numberOfSpecies

        call Vector_constructorInteger (order, numberOfSpecies, 0 )
        order%values(i)=1
        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        cc = 0
        !!Singles
        if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
          do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
            call Vector_constructor (occupiedCode(i), int( order%values(i),4), 0.0_8)
            occupiedCode(i)%values(1)=m
            do p=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1,int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
              if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                unoccupiedCode(i)%values(1)=p
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

      !!Doubles of the same specie
      do i=1, numberOfSpecies
        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        call Vector_constructorInteger (order, numberOfSpecies, 0 )
        order%values(i)=2
        cc = 0
        if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
           do m=1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
              do n=m+1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                 call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8 )
                 !occupiedCode%values(i,1)=m*1024+n
                 occupiedCode(i)%values(1)=m
                 occupiedCode(i)%values(2)=n
                 do p = int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                   do q=p+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                 !do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )

                    if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) ) then !! alpha -> alpha, beta -> beta
                          call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                          !unoccupiedCode%values(i)=p*1024+q
                          unoccupiedCode(i)%values(1)=p
                          unoccupiedCode(i)%values(2)=q
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

      end do !! species

      !!Doubles of different species
      do i=1, numberOfSpecies

        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital

        do j=i+1, numberOfSpecies

          call Vector_constructorInteger (order, numberOfSpecies, 0 )

          order%values(i)=1
          order%values(j)=1
          cc = 0

             do ii=1, numberOfSpecies
               call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
               call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
             end do


          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
            do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
              do n=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))

                call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                occupiedCode(i)%values(1)=m
                occupiedCode(j)%values(1)=n
                do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i))
                  do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                    call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                    call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                    unoccupiedCode(i)%values(1)=p
                    unoccupiedCode(j)%values(1)=q
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
      end do !! species

      !!Triples of the same specie (3)
      do i=1, numberOfSpecies

        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        call Vector_constructorInteger (order, numberOfSpecies, 0 )
        order%values(i)=3

        if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 3 ) then
          do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
            do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
              do mc=mb+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))

                call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                occupiedCode(i)%values(1)=ma
                occupiedCode(i)%values(2)=mb
                occupiedCode(i)%values(3)=mc

                do pa= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                        int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                  do pb=pa+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                    do pc=pb+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                      if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                          mod(mb,lambda) == mod(pb,lambda) .and. &
                          mod(mc,lambda) == mod(pc,lambda) ) then !! alpha -> alpha, beta -> beta

                          call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                          unoccupiedCode(i)%values(1)=pa
                          unoccupiedCode(i)%values(2)=pb
                          unoccupiedCode(i)%values(3)=pc
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
            end do
          end do
        end if
      end do !! species

       !!Triples (21)  problem here!!
      do i=1, numberOfSpecies

        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        call Vector_constructorInteger (order, numberOfSpecies, 0 )

        do j=1, numberOfSpecies
          if ( j /= i ) then

            order%values = 0
            order%values(i)=2
            order%values(j)=1

            do ii=1, numberOfSpecies
              call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
              call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
            end do


            if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
               .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
              do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                  do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )

                    call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                    occupiedCode(i)%values(1)=ma
                    occupiedCode(i)%values(2)=mb
                    call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                    occupiedCode(j)%values(1)=mc

                    do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                           int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                      do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                        do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, &
                          int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                          if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                               mod(mb,lambda) == mod(pb,lambda) .and. &
                               mod(mc,lambda) == mod(pc,lambda) )  then !! alpha -> alpha, beta -> beta

                            call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                            unoccupiedCode(i)%values(1)=pa
                            unoccupiedCode(i)%values(2)=pb
                            call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                            unoccupiedCode(j)%values(1)=pc
                            c=c+1
                            cc=cc+1
                            call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                              ConfigurationInteraction_instance%numberOfConfigurations,c) 

                          end if
                        end do
                      end do
                    end do  
                  end do
                end do
              end do
            end if
          end if
        end do

      end do !! species

      !!Triples (111)
       do i=1, numberOfSpecies

         call Vector_constructorInteger (order, numberOfSpecies, 0 )

         do j=i+1, numberOfSpecies
           order%values(j)=1
           do k=j+1, numberOfSpecies

             order%values=0
             order%values(i)=1
             order%values(j)=1
             order%values(k)=1
             do ii=1, numberOfSpecies
               call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
               call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
             end do

             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. & 
               ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 .and. &
               ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 ) then
               do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                 do n=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                   do r=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) )

                     call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                     occupiedCode(i)%values(1)=m
                     call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                     occupiedCode(j)%values(1)=n
                     call Vector_constructor (occupiedCode(k), int(order%values(k),4), 0.0_8)
                     occupiedCode(k)%values(1)=r

                     do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                       do q= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                         do s= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(k))
                            call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                            unoccupiedCode(i)%values(1)=p
                            call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                            unoccupiedCode(j)%values(1)=q
                            call Vector_constructor (unoccupiedCode(k), int(order%values(k),4), 0.0_8 )
                            unoccupiedCode(k)%values(1)=s
                            c=c+1
                            cc=cc+1
                            call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                              ConfigurationInteraction_instance%numberOfConfigurations,c) 


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

    case ("CISDTQ")
    !!Build configurations
    !!Ground State
      c=1
      call Vector_constructorInteger (order, numberOfSpecies, 0 )

      do i=1, numberOfSpecies
        call Vector_constructor (occupiedCode(i), 1, 0.0_8)
        call Vector_constructor (unoccupiedCode(i), 1, 0.0_8)
      end do

      call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations, c) 

      do i=1, numberOfSpecies

        call Vector_constructorInteger (order, numberOfSpecies, 0 )
        order%values(i)=1
        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        cc = 0
        !!Singles
        if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
          do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
            call Vector_constructor (occupiedCode(i), int( order%values(i),4), 0.0_8)
            occupiedCode(i)%values(1)=m
            do p=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1,int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
              if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                unoccupiedCode(i)%values(1)=p
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

      !!Doubles of the same specie
      do i=1, numberOfSpecies
        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        call Vector_constructorInteger (order, numberOfSpecies, 0 )
        order%values(i)=2
        cc = 0
        if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
           do m=1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
              do n=m+1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                 call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8 )
                 !occupiedCode%values(i,1)=m*1024+n
                 occupiedCode(i)%values(1)=m
                 occupiedCode(i)%values(2)=n
                 do p = int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                   do q=p+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                 !do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )

                    if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) ) then !! alpha -> alpha, beta -> beta
                          call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                          !unoccupiedCode%values(i)=p*1024+q
                          unoccupiedCode(i)%values(1)=p
                          unoccupiedCode(i)%values(2)=q
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

      end do !! species

      !!Doubles of different species
      do i=1, numberOfSpecies

        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital

        do j=i+1, numberOfSpecies

          call Vector_constructorInteger (order, numberOfSpecies, 0 )

          order%values(i)=1
          order%values(j)=1
          cc = 0

             do ii=1, numberOfSpecies
               call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
               call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
             end do


          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
            do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
              do n=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))

                call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                occupiedCode(i)%values(1)=m
                occupiedCode(j)%values(1)=n
                do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i))
                  do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                    call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                    call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                    unoccupiedCode(i)%values(1)=p
                    unoccupiedCode(j)%values(1)=q
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
      end do !! species

      !!Triples of the same specie (3)
      do i=1, numberOfSpecies

        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        call Vector_constructorInteger (order, numberOfSpecies, 0 )
        order%values(i)=3

        if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 3 ) then
          do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
            do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
              do mc=mb+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))

                call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                occupiedCode(i)%values(1)=ma
                occupiedCode(i)%values(2)=mb
                occupiedCode(i)%values(3)=mc

                do pa= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                        int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                  do pb=pa+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                    do pc=pb+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                      if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                          mod(mb,lambda) == mod(pb,lambda) .and. &
                          mod(mc,lambda) == mod(pc,lambda) ) then !! alpha -> alpha, beta -> beta

                          call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                          unoccupiedCode(i)%values(1)=pa
                          unoccupiedCode(i)%values(2)=pb
                          unoccupiedCode(i)%values(3)=pc
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
            end do
          end do
        end if
      end do !! species

       !!Triples (21)  problem here!!
      do i=1, numberOfSpecies

        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        call Vector_constructorInteger (order, numberOfSpecies, 0 )

        do j=1, numberOfSpecies
          if ( j /= i ) then

            order%values = 0
            order%values(i)=2
            order%values(j)=1

            do ii=1, numberOfSpecies
              call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
              call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
            end do


            if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
               .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
              do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                  do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )

                    call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                    occupiedCode(i)%values(1)=ma
                    occupiedCode(i)%values(2)=mb
                    call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                    occupiedCode(j)%values(1)=mc

                    do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                           int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                      do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                        do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, &
                          int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                          if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                               mod(mb,lambda) == mod(pb,lambda) .and. &
                               mod(mc,lambda) == mod(pc,lambda) )  then !! alpha -> alpha, beta -> beta

                            call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                            unoccupiedCode(i)%values(1)=pa
                            unoccupiedCode(i)%values(2)=pb
                            call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                            unoccupiedCode(j)%values(1)=pc
                            c=c+1
                            cc=cc+1
                            call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                              ConfigurationInteraction_instance%numberOfConfigurations,c) 

                          end if
                        end do
                      end do
                    end do  
                  end do
                end do
              end do
            end if
          end if
        end do

      end do !! species

      !!Triples (111)
       do i=1, numberOfSpecies

         call Vector_constructorInteger (order, numberOfSpecies, 0 )

         do j=i+1, numberOfSpecies
           order%values(j)=1
           do k=j+1, numberOfSpecies

             order%values=0
             order%values(i)=1
             order%values(j)=1
             order%values(k)=1
             do ii=1, numberOfSpecies
               call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
               call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
             end do

             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. & 
               ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 .and. &
               ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 ) then
               do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                 do n=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                   do r=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) )

                     call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                     occupiedCode(i)%values(1)=m
                     call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                     occupiedCode(j)%values(1)=n
                     call Vector_constructor (occupiedCode(k), int(order%values(k),4), 0.0_8)
                     occupiedCode(k)%values(1)=r

                     do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                       do q= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                         do s= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(k))
                            call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                            unoccupiedCode(i)%values(1)=p
                            call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                            unoccupiedCode(j)%values(1)=q
                            call Vector_constructor (unoccupiedCode(k), int(order%values(k),4), 0.0_8 )
                            unoccupiedCode(k)%values(1)=s
                            c=c+1
                            cc=cc+1
                            call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                              ConfigurationInteraction_instance%numberOfConfigurations,c) 


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

       !!Quadruples of the same species (4)
       do i=1, numberOfSpecies

        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        call Vector_constructorInteger (order, numberOfSpecies, 0 )
        order%values(i)=4

         if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 4 ) then
           do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
             do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
               do mc=mb+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                 do md=mc+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))

                     call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                     occupiedCode(i)%values(1)=ma
                     occupiedCode(i)%values(2)=mb
                     occupiedCode(i)%values(3)=mc
                     occupiedCode(i)%values(4)=md


                   do pa= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                           int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                     do pb=pa+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                       do pc=pb+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                         do pd=pc+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                           if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                               mod(mb,lambda) == mod(pb,lambda) .and. &
                               mod(mc,lambda) == mod(pc,lambda) .and. &
                               mod(md,lambda) == mod(pd,lambda) ) then !! alpha -> alpha, beta -> beta
                             call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                             unoccupiedCode(i)%values(1)=pa
                             unoccupiedCode(i)%values(2)=pb
                             unoccupiedCode(i)%values(3)=pc
                             unoccupiedCode(i)%values(4)=pd
                             cc=cc+1
                             c=c+1
                             call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                              ConfigurationInteraction_instance%numberOfConfigurations,c) 

                           end if
                         end do
                       end do
                     end do
                   end do
                 end do
               end do
             end do
           end do
         end if

       end do !! species

       !!Quadruples (31)
       do i=1, numberOfSpecies

         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         call Vector_constructorInteger (order, numberOfSpecies, 0 )

         do j=1, numberOfSpecies
           if ( j /= i ) then

             order%values=0
             order%values(i)=3
             order%values(j)=1
             do ii=1, numberOfSpecies
               call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
               call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
             end do


             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 3 &
                .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
                do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                  do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                    do mc=mb+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                      do md=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )

                        call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                        occupiedCode(i)%values(1)=ma
                        occupiedCode(i)%values(2)=mb
                        occupiedCode(i)%values(3)=mc
                        call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                        occupiedCode(j)%values(1)=md

                        do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                               int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                          do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                            do pc=pb+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                              do pd= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, &
                                int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                                if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                     mod(mb,lambda) == mod(pb,lambda) .and. &
                                     mod(mc,lambda) == mod(pc,lambda) .and. &
                                     mod(md,lambda) == mod(pd,lambda) )  then !! alpha -> alpha, beta -> beta
                                  call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                                  unoccupiedCode(i)%values(1)=pa
                                  unoccupiedCode(i)%values(2)=pb
                                  unoccupiedCode(i)%values(3)=pc
                                  call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                                  unoccupiedCode(j)%values(1)=pd
                                  cc=cc+1
                                  c=c+1
                                  call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                                    ConfigurationInteraction_instance%numberOfConfigurations,c) 
                                end if
                              end do
                            end do
                          end do
                        end do  
                      end do  
                    end do
                  end do
                end do
             end if
           end if
         end do

       end do !! species

       !!Quadruples (22)
       do i=1, numberOfSpecies

         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         call Vector_constructorInteger (order, numberOfSpecies, 0 )

         do j=i+1, numberOfSpecies

           otherlambda=ConfigurationInteraction_instance%lambda%values(j) !Particles per orbital
           !call Vector_constructorInteger (order, numberOfSpecies, 0 )
           order%values=0
           order%values(i)=2
           order%values(j)=2
           do ii=1, numberOfSpecies
             call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
             call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
           end do


           if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
              .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 2 ) then
              do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                  do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                    do md=mc+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )

                      call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                      occupiedCode(i)%values(1)=ma
                      occupiedCode(i)%values(2)=mb
                      call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                      occupiedCode(j)%values(1)=mc
                      occupiedCode(j)%values(2)=md

                      do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                             int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                        do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                          do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1,&
                             int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                            do pd=pc+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                              if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                   mod(mb,lambda) == mod(pb,lambda) .and. &
                                   mod(mc,otherlambda) == mod(pc,otherlambda) .and. &
                                   mod(md,otherlambda) == mod(pd,otherlambda) )  then !! alpha -> alpha, beta -> beta

                                  call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                                  unoccupiedCode(i)%values(1)=pa
                                  unoccupiedCode(i)%values(2)=pb
                                  call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                                  unoccupiedCode(j)%values(1)=pc
                                  unoccupiedCode(j)%values(2)=pd
                                  cc=cc+1
                                  c=c+1
                                  call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                                   ConfigurationInteraction_instance%numberOfConfigurations,c) 

                              end if
                            end do
                          end do
                        end do
                      end do  
                    end do  
                  end do
                end do
              end do
           end if
         end do

       end do !! species

       !!quadruples of diff specie (1111)
       do i=1, numberOfSpecies

         call Vector_constructorInteger (order, numberOfSpecies, 0 )

         if ( numberOfSpecies > 3 ) then
          do j=i+1, numberOfSpecies
            do k=j+1, numberOfSpecies
              do l=k+1, numberOfSpecies

                order%values=0
                order%values(i)=1
                order%values(j)=1
                order%values(k)=1
                order%values(l)=1
                do ii=1, numberOfSpecies
                  call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
                  call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
                end do



                if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. & 
                   ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 .and. &
                   ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 .and. &
                   ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(l) .ge. 1 ) then
                  do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                    do mb=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                      do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) )
                        do md=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(l) )

                          call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                          occupiedCode(i)%values(1)=ma
                          call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                          occupiedCode(j)%values(1)=mb
                          call Vector_constructor (occupiedCode(k), int(order%values(k),4), 0.0_8)
                          occupiedCode(k)%values(1)=mc
                          call Vector_constructor (occupiedCode(l), int(order%values(l),4), 0.0_8)
                          occupiedCode(l)%values(1)=md


                          do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                            do pb= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                              do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(k))
                                do pd= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(l))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(l))
                                  call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                                  unoccupiedCode(i)%values(1)=pa
                                  call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                                  unoccupiedCode(j)%values(1)=pb
                                  call Vector_constructor (unoccupiedCode(k), int(order%values(k),4), 0.0_8 )
                                  unoccupiedCode(k)%values(1)=pc
                                  call Vector_constructor (unoccupiedCode(l), int(order%values(l),4), 0.0_8 )
                                  unoccupiedCode(l)%values(1)=pd
                                  cc=cc+1
                                  c=c+1
                                  call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                                    ConfigurationInteraction_instance%numberOfConfigurations,c) 

                                end do
                              end do
                            end do
                          end do
                        end do
                      end do
                    end do
                  end do

                end if

              end do
            end do
          end do
        end if

       end do !! species

       !!Quadruples (211)
       do i=1, numberOfSpecies

         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         call Vector_constructorInteger (order, numberOfSpecies, 0 )

         do j=1, numberOfSpecies
           if ( j /= i ) then
             do k=j+1, numberOfSpecies
               if ( k /= j .and. k /=i ) then
                 order%values=0
                 order%values(i)=2
                 order%values(j)=1
                 order%values(k)=1
                 do ii=1, numberOfSpecies
                   call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
                   call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
                 end do

                 if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
                    .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 &
                    .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 ) then
                   do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                     do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                       do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                         do md=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) )

                            call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                            occupiedCode(i)%values(1)=ma
                            occupiedCode(i)%values(2)=mb
                            call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                            occupiedCode(j)%values(1)=mc
                            call Vector_constructor (occupiedCode(k), int(order%values(k),4), 0.0_8)
                            occupiedCode(k)%values(1)=md

                           do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                                  int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                             do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                               do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1,&
                                  int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                                 do pd= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1,&
                                    int(ConfigurationInteraction_instance%numberOfOrbitals%values(k) )
                                   if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                        mod(mb,lambda) == mod(pb,lambda) .and. &
                                        mod(mc,lambda) == mod(pc,lambda) .and. &
                                        mod(md,lambda) == mod(pd,lambda) )  then !! alpha -> alpha, beta -> beta
                                     call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                                     unoccupiedCode(i)%values(1)=pa
                                     unoccupiedCode(i)%values(2)=pb
                                     call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                                     unoccupiedCode(j)%values(1)=pc
                                     call Vector_constructor (unoccupiedCode(k), int(order%values(k),4), 0.0_8 )
                                     unoccupiedCode(k)%values(1)=pd
                                     cc=cc+1
                                     c=c+1
                                     call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                                       ConfigurationInteraction_instance%numberOfConfigurations,c) 
                                   end if
                                 end do
                               end do
                             end do
                           end do  
                         end do  
                       end do
                     end do
                   end do
                 end if
               end if
             end do
           end if
         end do

      end do !! species

    case ("CISDTQQ")
    !!Build configurations
    !!Ground State
      c=1
      call Vector_constructorInteger (order, numberOfSpecies, 0 )

      do i=1, numberOfSpecies
        call Vector_constructor (occupiedCode(i), 1, 0.0_8)
        call Vector_constructor (unoccupiedCode(i), 1, 0.0_8)
      end do

      call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, ConfigurationInteraction_instance%numberOfConfigurations, c) 

      do i=1, numberOfSpecies

        call Vector_constructorInteger (order, numberOfSpecies, 0 )
        order%values(i)=1
        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        cc = 0
        !!Singles
        if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 ) then
          do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
            call Vector_constructor (occupiedCode(i), int( order%values(i),4), 0.0_8)
            occupiedCode(i)%values(1)=m
            do p=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1,int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
              if ( mod(m,lambda) == mod(p,lambda) ) then !! alpha -> alpha, beta -> beta
                call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                unoccupiedCode(i)%values(1)=p
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

      !!Doubles of the same specie
      do i=1, numberOfSpecies
        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        call Vector_constructorInteger (order, numberOfSpecies, 0 )
        order%values(i)=2
        cc = 0
        if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 ) then
           do m=1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
              do n=m+1,int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                 call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8 )
                 !occupiedCode%values(i,1)=m*1024+n
                 occupiedCode(i)%values(1)=m
                 occupiedCode(i)%values(2)=n
                 do p = int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                   do q=p+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                 !do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )

                    if ( mod(m,lambda) == mod(p,lambda) .and.  mod(n,lambda) == mod(q,lambda) ) then !! alpha -> alpha, beta -> beta
                          call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                          !unoccupiedCode%values(i)=p*1024+q
                          unoccupiedCode(i)%values(1)=p
                          unoccupiedCode(i)%values(2)=q
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

      end do !! species

      !!Doubles of different species
      do i=1, numberOfSpecies

        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital

        do j=i+1, numberOfSpecies

          call Vector_constructorInteger (order, numberOfSpecies, 0 )

          order%values(i)=1
          order%values(j)=1
          cc = 0

             do ii=1, numberOfSpecies
               call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
               call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
             end do


          if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
            do m=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
              do n=1, int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))

                call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                occupiedCode(i)%values(1)=m
                occupiedCode(j)%values(1)=n
                do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i))
                  do q=int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                    call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                    call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                    unoccupiedCode(i)%values(1)=p
                    unoccupiedCode(j)%values(1)=q
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
      end do !! species

      !!Triples of the same specie (3)
      do i=1, numberOfSpecies

        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        call Vector_constructorInteger (order, numberOfSpecies, 0 )
        order%values(i)=3

        if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 3 ) then
          do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
            do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
              do mc=mb+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))

                call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                occupiedCode(i)%values(1)=ma
                occupiedCode(i)%values(2)=mb
                occupiedCode(i)%values(3)=mc

                do pa= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                        int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                  do pb=pa+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                    do pc=pb+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                      if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                          mod(mb,lambda) == mod(pb,lambda) .and. &
                          mod(mc,lambda) == mod(pc,lambda) ) then !! alpha -> alpha, beta -> beta

                          call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                          unoccupiedCode(i)%values(1)=pa
                          unoccupiedCode(i)%values(2)=pb
                          unoccupiedCode(i)%values(3)=pc
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
            end do
          end do
        end if
      end do !! species

       !!Triples (21)  problem here!!
      do i=1, numberOfSpecies

        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        call Vector_constructorInteger (order, numberOfSpecies, 0 )

        do j=1, numberOfSpecies
          if ( j /= i ) then

            order%values = 0
            order%values(i)=2
            order%values(j)=1

            do ii=1, numberOfSpecies
              call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
              call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
            end do


            if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
               .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
              do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                  do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )

                    call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                    occupiedCode(i)%values(1)=ma
                    occupiedCode(i)%values(2)=mb
                    call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                    occupiedCode(j)%values(1)=mc

                    do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                           int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                      do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                        do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, &
                          int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                          if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                               mod(mb,lambda) == mod(pb,lambda) .and. &
                               mod(mc,lambda) == mod(pc,lambda) )  then !! alpha -> alpha, beta -> beta

                            call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                            unoccupiedCode(i)%values(1)=pa
                            unoccupiedCode(i)%values(2)=pb
                            call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                            unoccupiedCode(j)%values(1)=pc
                            c=c+1
                            cc=cc+1
                            call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                              ConfigurationInteraction_instance%numberOfConfigurations,c) 

                          end if
                        end do
                      end do
                    end do  
                  end do
                end do
              end do
            end if
          end if
        end do

      end do !! species

      !!Triples (111)
       do i=1, numberOfSpecies

         call Vector_constructorInteger (order, numberOfSpecies, 0 )

         do j=i+1, numberOfSpecies
           order%values(j)=1
           do k=j+1, numberOfSpecies

             order%values=0
             order%values(i)=1
             order%values(j)=1
             order%values(k)=1
             do ii=1, numberOfSpecies
               call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
               call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
             end do

             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. & 
               ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 .and. &
               ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 ) then
               do m=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                 do n=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                   do r=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) )

                     call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                     occupiedCode(i)%values(1)=m
                     call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                     occupiedCode(j)%values(1)=n
                     call Vector_constructor (occupiedCode(k), int(order%values(k),4), 0.0_8)
                     occupiedCode(k)%values(1)=r

                     do p= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                       do q= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                         do s= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(k))
                            call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                            unoccupiedCode(i)%values(1)=p
                            call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                            unoccupiedCode(j)%values(1)=q
                            call Vector_constructor (unoccupiedCode(k), int(order%values(k),4), 0.0_8 )
                            unoccupiedCode(k)%values(1)=s
                            c=c+1
                            cc=cc+1
                            call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                              ConfigurationInteraction_instance%numberOfConfigurations,c) 


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

       !!Quadruples of the same species (4)
       do i=1, numberOfSpecies

        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        call Vector_constructorInteger (order, numberOfSpecies, 0 )
        order%values(i)=4

         if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 4 ) then
           do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
             do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
               do mc=mb+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
                 do md=mc+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))

                     call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                     occupiedCode(i)%values(1)=ma
                     occupiedCode(i)%values(2)=mb
                     occupiedCode(i)%values(3)=mc
                     occupiedCode(i)%values(4)=md


                   do pa= int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                           int( ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                     do pb=pa+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                       do pc=pb+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                         do pd=pc+1,int( ConfigurationInteraction_instance%numberOfOrbitals%values(i)) 
                           if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                               mod(mb,lambda) == mod(pb,lambda) .and. &
                               mod(mc,lambda) == mod(pc,lambda) .and. &
                               mod(md,lambda) == mod(pd,lambda) ) then !! alpha -> alpha, beta -> beta
                             call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                             unoccupiedCode(i)%values(1)=pa
                             unoccupiedCode(i)%values(2)=pb
                             unoccupiedCode(i)%values(3)=pc
                             unoccupiedCode(i)%values(4)=pd
                             cc=cc+1
                             c=c+1
                             call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                              ConfigurationInteraction_instance%numberOfConfigurations,c) 

                           end if
                         end do
                       end do
                     end do
                   end do
                 end do
               end do
             end do
           end do
         end if

       end do !! species

       !!Quadruples (31)
       do i=1, numberOfSpecies

         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         call Vector_constructorInteger (order, numberOfSpecies, 0 )

         do j=1, numberOfSpecies
           if ( j /= i ) then

             order%values=0
             order%values(i)=3
             order%values(j)=1
             do ii=1, numberOfSpecies
               call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
               call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
             end do


             if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 3 &
                .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 ) then
                do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                  do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                    do mc=mb+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                      do md=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )

                        call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                        occupiedCode(i)%values(1)=ma
                        occupiedCode(i)%values(2)=mb
                        occupiedCode(i)%values(3)=mc
                        call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                        occupiedCode(j)%values(1)=md

                        do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                               int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                          do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                            do pc=pb+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                              do pd= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, &
                                int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                                if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                     mod(mb,lambda) == mod(pb,lambda) .and. &
                                     mod(mc,lambda) == mod(pc,lambda) .and. &
                                     mod(md,lambda) == mod(pd,lambda) )  then !! alpha -> alpha, beta -> beta
                                  call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                                  unoccupiedCode(i)%values(1)=pa
                                  unoccupiedCode(i)%values(2)=pb
                                  unoccupiedCode(i)%values(3)=pc
                                  call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                                  unoccupiedCode(j)%values(1)=pd
                                  cc=cc+1
                                  c=c+1
                                  call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                                    ConfigurationInteraction_instance%numberOfConfigurations,c) 
                                end if
                              end do
                            end do
                          end do
                        end do  
                      end do  
                    end do
                  end do
                end do
             end if
           end if
         end do

       end do !! species

       !!Quadruples (22)
       do i=1, numberOfSpecies

         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         call Vector_constructorInteger (order, numberOfSpecies, 0 )

         do j=i+1, numberOfSpecies

           otherlambda=ConfigurationInteraction_instance%lambda%values(j) !Particles per orbital
           !call Vector_constructorInteger (order, numberOfSpecies, 0 )
           order%values=0
           order%values(i)=2
           order%values(j)=2
           do ii=1, numberOfSpecies
             call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
             call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
           end do


           if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
              .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 2 ) then
              do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                  do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                    do md=mc+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )

                      call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                      occupiedCode(i)%values(1)=ma
                      occupiedCode(i)%values(2)=mb
                      call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                      occupiedCode(j)%values(1)=mc
                      occupiedCode(j)%values(2)=md

                      do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                             int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                        do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                          do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1,&
                             int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                            do pd=pc+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                              if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                   mod(mb,lambda) == mod(pb,lambda) .and. &
                                   mod(mc,otherlambda) == mod(pc,otherlambda) .and. &
                                   mod(md,otherlambda) == mod(pd,otherlambda) )  then !! alpha -> alpha, beta -> beta

                                  call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                                  unoccupiedCode(i)%values(1)=pa
                                  unoccupiedCode(i)%values(2)=pb
                                  call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                                  unoccupiedCode(j)%values(1)=pc
                                  unoccupiedCode(j)%values(2)=pd
                                  cc=cc+1
                                  c=c+1
                                  call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                                   ConfigurationInteraction_instance%numberOfConfigurations,c) 

                              end if
                            end do
                          end do
                        end do
                      end do  
                    end do  
                  end do
                end do
              end do
           end if
         end do

       end do !! species

       !!quadruples of diff specie (1111)
       do i=1, numberOfSpecies

         call Vector_constructorInteger (order, numberOfSpecies, 0 )

         if ( numberOfSpecies > 3 ) then
          do j=i+1, numberOfSpecies
            do k=j+1, numberOfSpecies
              do l=k+1, numberOfSpecies

                order%values=0
                order%values(i)=1
                order%values(j)=1
                order%values(k)=1
                order%values(l)=1
                do ii=1, numberOfSpecies
                  call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
                  call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
                end do



                if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 1 .and. & 
                   ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 .and. &
                   ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 .and. &
                   ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(l) .ge. 1 ) then
                  do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                    do mb=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                      do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) )
                        do md=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(l) )

                          call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                          occupiedCode(i)%values(1)=ma
                          call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                          occupiedCode(j)%values(1)=mb
                          call Vector_constructor (occupiedCode(k), int(order%values(k),4), 0.0_8)
                          occupiedCode(k)%values(1)=mc
                          call Vector_constructor (occupiedCode(l), int(order%values(l),4), 0.0_8)
                          occupiedCode(l)%values(1)=md


                          do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                            do pb= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j))
                              do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(k))
                                do pd= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(l))+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(l))
                                  call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                                  unoccupiedCode(i)%values(1)=pa
                                  call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                                  unoccupiedCode(j)%values(1)=pb
                                  call Vector_constructor (unoccupiedCode(k), int(order%values(k),4), 0.0_8 )
                                  unoccupiedCode(k)%values(1)=pc
                                  call Vector_constructor (unoccupiedCode(l), int(order%values(l),4), 0.0_8 )
                                  unoccupiedCode(l)%values(1)=pd
                                  cc=cc+1
                                  c=c+1
                                  call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                                    ConfigurationInteraction_instance%numberOfConfigurations,c) 

                                end do
                              end do
                            end do
                          end do
                        end do
                      end do
                    end do
                  end do

                end if

              end do
            end do
          end do
        end if

       end do !! species

       !!Quadruples (211)
       do i=1, numberOfSpecies

         lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
         call Vector_constructorInteger (order, numberOfSpecies, 0 )

         do j=1, numberOfSpecies
           if ( j /= i ) then
             do k=j+1, numberOfSpecies
               if ( k /= j .and. k /=i ) then
                 order%values=0
                 order%values(i)=2
                 order%values(j)=1
                 order%values(k)=1
                 do ii=1, numberOfSpecies
                   call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
                   call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
                 end do

                 if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
                    .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 1 &
                    .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 ) then
                   do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                     do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                       do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                         do md=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) )

                            call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                            occupiedCode(i)%values(1)=ma
                            occupiedCode(i)%values(2)=mb
                            call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                            occupiedCode(j)%values(1)=mc
                            call Vector_constructor (occupiedCode(k), int(order%values(k),4), 0.0_8)
                            occupiedCode(k)%values(1)=md

                           do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                                  int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                             do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                               do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1,&
                                  int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                                 do pd= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1,&
                                    int(ConfigurationInteraction_instance%numberOfOrbitals%values(k) )
                                   if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                        mod(mb,lambda) == mod(pb,lambda) .and. &
                                        mod(mc,lambda) == mod(pc,lambda) .and. &
                                        mod(md,lambda) == mod(pd,lambda) )  then !! alpha -> alpha, beta -> beta
                                     call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                                     unoccupiedCode(i)%values(1)=pa
                                     unoccupiedCode(i)%values(2)=pb
                                     call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                                     unoccupiedCode(j)%values(1)=pc
                                     call Vector_constructor (unoccupiedCode(k), int(order%values(k),4), 0.0_8 )
                                     unoccupiedCode(k)%values(1)=pd
                                     cc=cc+1
                                     c=c+1
                                     call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                                       ConfigurationInteraction_instance%numberOfConfigurations,c) 
                                   end if
                                 end do
                               end do
                             end do
                           end do  
                         end do  
                       end do
                     end do
                   end do
                 end if
               end if
             end do
           end if
         end do

      end do !! species

      !!Quintuples (221)
      do i=1, numberOfSpecies

        lambda=ConfigurationInteraction_instance%lambda%values(i) !Particles per orbital
        call Vector_constructorInteger (order, numberOfSpecies, 0 )

        do j=i+1, numberOfSpecies

          otherlambda=ConfigurationInteraction_instance%lambda%values(j) !Particles per orbital

          do k=j+1, numberOfSpecies

            order%values=0
            order%values(i)=2
            order%values(j)=2
            order%values(k)=1
            do ii=1, numberOfSpecies
              call Vector_constructor (occupiedCode(ii), 1, 0.0_8)
              call Vector_constructor (unoccupiedCode(ii), 1, 0.0_8)
            end do

            if (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) .ge. 2 &
              .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) .ge. 2 &
              .and. ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) .ge. 1 ) then
              do ma=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                do mb=ma+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
                  do mc=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                    do md=mc+1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j) )
                      do me=1, int( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k) )

                         call Vector_constructor (occupiedCode(i), int(order%values(i),4), 0.0_8)
                         occupiedCode(i)%values(1)=ma
                         occupiedCode(i)%values(2)=mb
                         call Vector_constructor (occupiedCode(j), int(order%values(j),4), 0.0_8)
                         occupiedCode(j)%values(1)=mc
                         occupiedCode(j)%values(2)=md
                         call Vector_constructor (occupiedCode(k), int(order%values(k),4), 0.0_8)
                         occupiedCode(k)%values(1)=me


                        do pa= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))+1, &
                               int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                          do pb=pa+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
                            do pc= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j))+1,&
                               int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                              do pd=pc+1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(j) )
                                do pe= int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(k))+1,&
                                   int(ConfigurationInteraction_instance%numberOfOrbitals%values(k) )
                                  if ( mod(ma,lambda) == mod(pa,lambda) .and. & 
                                       mod(mb,lambda) == mod(pb,lambda) .and. &
                                       mod(mc,otherlambda) == mod(pc,otherlambda) .and. &
                                       mod(md,otherlambda) == mod(pd,otherlambda) .and. &
                                       mod(me,lambda) == mod(pe,lambda) )  then !! alpha -> alpha, beta -> beta
                                     call Vector_constructor (unoccupiedCode(i), int(order%values(i),4), 0.0_8 )
                                     unoccupiedCode(i)%values(1)=pa
                                     unoccupiedCode(i)%values(2)=pb
                                     call Vector_constructor (unoccupiedCode(j), int(order%values(j),4), 0.0_8 )
                                     unoccupiedCode(j)%values(1)=pc
                                     unoccupiedCode(j)%values(2)=pd
                                     call Vector_constructor (unoccupiedCode(k), int(order%values(k),4), 0.0_8 )
                                     unoccupiedCode(k)%values(1)=pe
                                     cc=cc+1
                                     c=c+1
                                     call Configuration_constructor(ConfigurationInteraction_instance%configurations(c), occupiedCode, unoccupiedCode, order, &
                                      ConfigurationInteraction_instance%numberOfConfigurations,c) 
                                  end if
                                end do
                              end do
                            end do
                          end do
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

    case default

       call ConfigurationInteraction_exception( ERROR, "Configuration interactor constructor", "Correction level not implemented")

    end select

    ! do c = 1,  ConfigurationInteraction_instance%numberOfConfigurations 
    !    call Configuration_show(  ConfigurationInteraction_instance%configurations(c) )
    ! end do


  end subroutine ConfigurationInteraction_buildConfigurations

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_buildHamiltonianMatrix()
    implicit none

    type(Configuration) :: auxConfigurationA, auxConfigurationB
    !integer(2), allocatable :: auxConfiguration(:,:), auxConfigurationA(:,:)
    integer(8) :: a,b,c
    integer :: size1, size2
    real(8) :: timeA, timeB
    real(8) :: CIenergyb, CIenergy
    integer(2) coupingCoefficient 
    integer(2), allocatable :: auxMatrix( :,:)

    !timeA = omp_get_wtime()
    !a,b configuration iterators

    size1 = size(ConfigurationInteraction_instance%configurations(1)%occupations, dim=1)
    size2 = size(ConfigurationInteraction_instance%configurations(1)%occupations, dim=2) 
    ! allocate(auxMatrix(size1,size2))
    !auxMatrix = 0

    !do a=1, ConfigurationInteraction_instance%numberOfConfigurations
    !  do b=a, ConfigurationInteraction_instance%numberOfConfigurations
!
    !    coupingCoefficient = ConfigurationInteraction_calculateCoupling( & 
    !      auxMatrix, &
    !      auxMatrix )
    !  end do
    !end do

    !timeB = omp_get_wtime()
    !write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for Building Couping coefficients : ", timeB - timeA ," (s)"

    !call Configuration_copyConstructor ( ConfigurationInteraction_instance%configurations(1), auxConfigurationA )
    !call Configuration_copyConstructor ( ConfigurationInteraction_instance%configurations(1), auxConfigurationB )
    !auxConfigurationA%occupations = 0
    !auxConfigurationB%occupations = 0

    !allocate ( auxConfiguration( size(ConfigurationInteraction_instance%configurations(1)%occupations, dim=1), &
    !   size(ConfigurationInteraction_instance%configurations(1)%occupations, dim=2)  ) )
    !size1 = size(ConfigurationInteraction_instance%configurations(1)%occupations, dim=1)
    !size2 = size(ConfigurationInteraction_instance%configurations(1)%occupations, dim=2) 

    !auxConfiguration = 0
    !auxConfigurationA = 0


    !allocate ( auxMatrix ( ConfigurationInteraction_instance%numberOfConfigurations, &
    !      ConfigurationInteraction_instance%numberOfConfigurations ))
    !auxMatrix=0

    call omp_set_num_threads(omp_get_max_threads())

    timeA = omp_get_wtime()

    do a=1, ConfigurationInteraction_instance%numberOfConfigurations
!$omp parallel & 
!$omp& private(b,CIenergy),&
!$omp& shared(ConfigurationInteraction_instance, HartreeFock_instance, a, size1, size2)
!$omp do 
      do b=a, ConfigurationInteraction_instance%numberOfConfigurations

          CIenergyb = ConfigurationInteraction_calculateCIenergyB( & 
          ConfigurationInteraction_instance%configurations(a)%occupations, &
          ConfigurationInteraction_instance%configurations(b)%occupations, size1, size2 )

          ConfigurationInteraction_instance%hamiltonianMatrix%values(a,b) = CIenergyB

      end do
!$omp end do nowait
!$omp end parallel
    end do

    !! symmetrize

  timeB = omp_get_wtime()
    
    do a=1, ConfigurationInteraction_instance%numberOfConfigurations
      do b=a, ConfigurationInteraction_instance%numberOfConfigurations
         ConfigurationInteraction_instance%hamiltonianMatrix%values(b,a)=ConfigurationInteraction_instance%hamiltonianMatrix%values(a,b)
      end do
    end do

  
  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for Building CI matrix : ", timeB - timeA ," (s)"

  end subroutine ConfigurationInteraction_buildHamiltonianMatrix

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_buildInitialCIMatrix()
    implicit none

    type(Configuration) :: auxConfigurationA, auxConfigurationB
    type (Vector8) :: diagonalHamiltonianMatrix
!    type (Vector) :: initialEigenValues
    type (Matrix) :: initialHamiltonianMatrix
    integer :: a,b,c,aa,bb
    real(8) :: timeA, timeB
    real(8) :: CIenergy
    integer :: initialCIMatrixSize 
    integer :: nproc

    initialCIMatrixSize = CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX 

    timeA = omp_get_wtime()
    !a,b configuration iterators
    call Configuration_copyConstructor ( ConfigurationInteraction_instance%configurations(1), auxConfigurationA )
    call Configuration_copyConstructor ( ConfigurationInteraction_instance%configurations(1), auxConfigurationB )

    call Vector_constructorInteger8 ( ConfigurationInteraction_instance%auxIndexCIMatrix, &
                              ConfigurationInteraction_instance%numberOfConfigurations, 0_8 ) 

    call Vector_constructor8 ( diagonalHamiltonianMatrix, &
                              ConfigurationInteraction_instance%numberOfConfigurations, 0.0_8 ) 

    print *, "  OMP Number of threads: " , omp_get_max_threads()
    nproc = omp_get_max_threads()

    call omp_set_num_threads(omp_get_max_threads())
    call omp_set_num_threads(nproc)

!$omp parallel & 
!$omp& private(a,b,CIenergy,auxConfigurationA,auxConfigurationB),&
!$omp& shared(ConfigurationInteraction_instance, HartreeFock_instance,diagonalHamiltonianMatrix)
!$omp do 
    do a=1, ConfigurationInteraction_instance%numberOfConfigurations
      auxConfigurationA%occupations = ConfigurationInteraction_instance%configurations(a)%occupations
      auxConfigurationB%occupations = ConfigurationInteraction_instance%configurations(a)%occupations
      CIenergy = ConfigurationInteraction_calculateCIenergyC(&
                      auxConfigurationA, auxConfigurationB )

      diagonalHamiltonianMatrix%values(a) = CIenergy
    end do
!$omp end do nowait
!$omp end parallel
    timeB = omp_get_wtime()
    write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time for Building Initial CI matrix : ", timeB - timeA ," (s)"


   !! save the unsorted diagonal Matrix )
   if ( trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) == "JADAMILU" .or. &
      trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) == "ARPACK"  ) then

     call Vector_copyConstructor8 (  ConfigurationInteraction_instance%diagonalHamiltonianMatrix, diagonalHamiltonianMatrix)

   end if

   !! To get only the lowest 300 values.
   call Vector_reverseSortElements8(diagonalHamiltonianMatrix, ConfigurationInteraction_instance%auxIndexCIMatrix, int(initialCIMatrixSize,8))

   call Matrix_constructor ( initialHamiltonianMatrix, int(initialCIMatrixSize,8) , &
                               int(initialCIMatrixSize,8) , 0.0_8 ) 


    do a=1, initialCIMatrixSize 
      initialHamiltonianMatrix%values(a,a) = diagonalHamiltonianMatrix%values(a)
      aa = ConfigurationInteraction_instance%auxIndexCIMatrix%values(a)
      do b=a+1, initialCIMatrixSize 
        bb = ConfigurationInteraction_instance%auxIndexCIMatrix%values(b)
        auxConfigurationA%occupations = ConfigurationInteraction_instance%configurations(aa)%occupations
        auxConfigurationB%occupations = ConfigurationInteraction_instance%configurations(bb)%occupations
        initialHamiltonianMatrix%values(a,b) = ConfigurationInteraction_calculateCIenergyC( &
                      auxConfigurationA, auxConfigurationB )

      end do
    end do

    !! diagonalize the initial matrix
    call Vector_constructor8 ( ConfigurationInteraction_instance%initialEigenValues, int(CONTROL_instance%NUMBER_OF_CI_STATES,8),  0.0_8)

    call Matrix_constructor (ConfigurationInteraction_instance%initialEigenVectors, &
           int(initialCIMatrixSize,8), &
           int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

    call Matrix_eigen_select ( initialHamiltonianMatrix, ConfigurationInteraction_instance%initialEigenValues, &
           1, int(CONTROL_instance%NUMBER_OF_CI_STATES,4), &  
           eigenVectors = ConfigurationInteraction_instance%initialEigenVectors, &
           flags = int(SYMMETRIC,4))

    !! cleaning
    call Vector_destructor8 ( diagonalHamiltonianMatrix )
!    call Vector_destructor ( initialEigenValues )
    call Matrix_destructor ( initialHamiltonianMatrix )

    timeB = omp_get_wtime()
    write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time for Building Initial CI matrix : ", timeB - timeA ," (s)"

  end subroutine ConfigurationInteraction_buildInitialCIMatrix

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_buildAndSaveCIMatrix()
    implicit none
    type(Configuration) :: auxConfigurationA, auxConfigurationB
    integer :: a,b,c,d,n, nproc,cc
    real(8) :: timeA, timeB
    character(50) :: CIFile
    integer :: CIUnit
    real(8) :: CIenergy
    integer, allocatable :: indexArray(:),auxIndexArray(:)
    real(8), allocatable :: energyArray(:),auxEnergyArray(:)
    integer :: starting, ending, step, maxConfigurations
    character(50) :: fileNumberA, fileNumberB
    integer, allocatable :: cmax(:)
    integer :: maxStackSize, i, ia, ib, ssize, ci,cj

    maxStackSize = CONTROL_instance%CI_STACK_SIZE 


    call Configuration_copyConstructor ( ConfigurationInteraction_instance%configurations(1), auxConfigurationA )
    call Configuration_copyConstructor ( ConfigurationInteraction_instance%configurations(1), auxConfigurationB )


    timeA = omp_get_wtime()

    CIFile = "lowdin.ci"
    CIUnit = 4

#ifdef intel
    open(unit=CIUnit, file=trim(CIFile), action = "write", form="unformatted", BUFFERED="YES")
#else
    open(unit=CIUnit, file=trim(CIFile), action = "write", form="unformatted")
#endif
    
    print *, "  OMP Number of threads: " , omp_get_max_threads()
    nproc = omp_get_max_threads()

    call omp_set_num_threads(omp_get_max_threads())
    call omp_set_num_threads(nproc)

    if (allocated(cmax)) deallocate(cmax)
    allocate(cmax(nproc))
    cmax = 0

    do a=1, ConfigurationInteraction_instance%numberOfConfigurations

      maxConfigurations = ConfigurationInteraction_instance%numberOfConfigurations - a + 1 
      step = ceiling ( real( maxConfigurations )/real(nproc))
      if (allocated(indexArray )) deallocate(indexArray)
      allocate (indexArray(maxConfigurations))
      indexArray = 0
      if (allocated(energyArray )) deallocate(energyArray)
      allocate (energyArray(maxConfigurations))
      energyArray = 0
      cmax = 0
!$omp parallel & 
!$omp& private(n,ending,starting,b,CIenergy),&
!$omp& private(auxConfigurationA ),&
!$omp& private(auxConfigurationB ),&
!$omp& firstprivate(nproc,maxConfigurations),&
!$omp& shared(cmax,indexArray,energyArray, HartreeFock_instance),&
!$omp& shared(ConfigurationInteraction_instance)
!$omp do 
      do n = 1, nproc

        ending = n * step  + a - 1
        starting = ending - step + 1

        if( ending > ConfigurationInteraction_instance%numberOfConfigurations ) then
         ending =ConfigurationInteraction_instance%numberOfConfigurations 
        end if  

        do b= starting, ending

          auxConfigurationA%occupations = ConfigurationInteraction_instance%configurations(a)%occupations
          auxConfigurationB%occupations = ConfigurationInteraction_instance%configurations(b)%occupations

          CIenergy = ConfigurationInteraction_calculateCIenergyC( & 
                      auxConfigurationA, auxConfigurationB )
          if ( abs(CIenergy) > 1E-9 ) then
            cmax(n) = cmax(n) +1   
            indexArray(b-a+1) = b
            energyArray(b-a+1) = CIenergy
          end if
        end do
      end do
!$omp end do nowait
!$omp end parallel
       
      c = sum(cmax)

      write(CIUnit) c
      write(CIUnit) a

      allocate (auxEnergyArray(c))
      allocate (auxIndexArray(c))

      cj = 0
      do ci = 1, maxConfigurations
        if ( indexArray(ci) > 0 ) then
          cj = cj + 1
          auxIndexArray(cj) = indexArray(ci)
          auxEnergyArray(cj) = energyArray(ci)
        end if
      end do

      do i = 1, ceiling(real(c) / real(maxStackSize) )
        ib = maxStackSize * i  
        ia = ib - maxStackSize + 1
        if ( ib > c ) ib = c
        write(CIUnit) auxIndexArray(ia:ib)
      end do
      deallocate(auxIndexArray)

      do i = 1, ceiling(real(c) / real(maxStackSize) )
        ib = maxStackSize * i  
        ia = ib - maxStackSize + 1
        if ( ib >  c ) ib = c
        write(CIUnit) auxEnergyArray(ia:ib)
      end do
      deallocate (auxEnergyArray)

!      write(CIUnit) pack( indexArray, indexArray /= 0 )
!      write(CIUnit) pack( energyArray, energyArray /= 0.0_8)

    end do

    write(CIUnit) -1

    close(CIUnit)

    deallocate(indexArray)
    deallocate(energyArray)

    timeB = omp_get_wtime()
  
    write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time for Building CI matrix : ", timeB - timeA ," (s)"

  end subroutine ConfigurationInteraction_buildAndSaveCIMatrix


  function ConfigurationInteraction_calculateCoupling(auxthisA, auxthisB) result (numberOfDiffOrbitals)
    implicit none
    !type(Configuration) :: auxthisA, auxthisB
    integer(2), intent(in) :: auxthisA(:,:), auxthisB(:,:)
    integer(2) :: numberOfDiffOrbitals
    integer :: numberOfSpecies

    numberOfSpecies = MolecularSystem_instance%numberOfQuantumSpecies
    numberOfDiffOrbitals = Configuration_checkCoincidenceB( auxthisA, auxthisB, numberOfSpecies )
    numberOfDiffOrbitals=1

  end function ConfigurationInteraction_calculateCoupling

  function ConfigurationInteraction_calculateCIenergyB( thisA, auxthisB, m, n) result (auxCIenergy)
    implicit none
    !type(Configuration), intent(in) :: thisA
    !type(Configuration), intent(in) :: auxthisB
    integer(2), intent(in) :: thisA(:,:), auxthisB(:,:)
    integer, intent(in) :: m,n
    integer :: i,j,s
    integer :: l,k,z,kk,ll
    integer :: numberOfSpecies
    integer :: numberOfOccupiedOrbitals
    real(8) :: kappa !positive or negative exchange
    real(8) :: twoParticlesEnergy
    real(8) :: couplingEnergy
    integer :: factor
    integer :: numberOfDiffOrbitals
    integer :: auxnumberOfOtherSpecieSpatialOrbitals
    integer :: auxIndex1, auxIndex2, auxIndex
    integer :: diffOrb1, diffOrb2, diffOrb3, diffOrb4, otherdiffOrb1,otherdiffOrb3
    real(8) :: auxCIenergy
    integer :: auxOcc
    integer(2) :: score, auxscore
    integer(2) :: diagonal
    logical(1) :: swap
    integer(2) :: auxthisA(m,n)

    numberOfSpecies = MolecularSystem_instance%numberOfQuantumSpecies

    !allocate(differentOrbitals (numberOfSpecies))

    !do ia = 1, thisA%nDeterminants 
    !  do ib = 1, thisB%nDeterminants 

        auxCIenergy = 0.0_8
   
        numberOfDiffOrbitals = Configuration_checkCoincidenceB( thisA, auxthisB, numberOfSpecies )

        factor = 1
        if  (  numberOfDiffOrbitals <= 2  ) then
          auxthisA = thisA
        else
          return
        end if

        if  (  numberOfDiffOrbitals == 1 .or. numberOfDiffOrbitals == 2  ) then
          !call Configuration_setAtMaximumCoincidenceC( auxthisA,auxthisB%occupations, m, n, numberOfSpecies, factor )
          !factor = Configuration_setAtMaximumCoincidenceB( numberOfSpecies)
          do s = 1, numberOfSpecies
              numberOfOccupiedOrbitals = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(s) 
        
              score = 0
              auxscore = 0
              diagonal = 0
      
              do i = 1, numberOfOccupiedOrbitals 
                  do j = 1, numberOfOccupiedOrbitals
                     if ( auxthisA(i,s) == auxthisB(j,s) ) then
                       auxscore = 1
                     else 
                       auxscore = 0
                     end if 
                     if ( i == j ) diagonal = diagonal + auxscore
                     score = score + auxscore
                  end do 
               end do 
        
              !do while ( (score) > diagonal )
              do k = 1, numberOfOccupiedOrbitals 

                  if ((score) > diagonal ) then
                  swap = .false. 
                  do i = 1, numberOfOccupiedOrbitals 
                    do j = 1, numberOfOccupiedOrbitals
                      if ( i /= j ) then
                        if ( auxthisA(i,s) == auxthisB(j,s) ) then
      
                          auxOcc = auxthisA(i,s)
                          auxthisA(i,s) = auxthisA(j,s)
                          auxthisA(j,s) = auxOcc
                          swap = .true.  
                        end if
                      end if
      
                      if ( swap .eqv. .true. ) exit
                    end do 
                    if ( swap .eqv. .true. ) exit
                  end do 
      
                diagonal = 0
                score = 0
                do i = 1, numberOfOccupiedOrbitals 
                  do j = 1, numberOfOccupiedOrbitals
      
                     if ( auxthisA(i,s) == auxthisB(j,s) ) then
                        auxscore = 1
                     else  
                       auxscore = 0
                     end if 
                     if ( i == j ) diagonal = diagonal + auxscore
                     score = score + auxscore
                  end do 
                end do 
                  factor = -1 * factor 
                else
                  exit
                end if
              end do! while
      
           end do
        end if

        !numberOfDiffOrbitals = 3
        select case (  numberOfDiffOrbitals )

        case (0)

          
              do i=1, numberOfSpecies

                 kappa = MolecularSystem_instance%species(i)%kappa
                 do kk=1, int( MolecularSystem_instance%species(i)%ocupationNumber) !! 1 is from a and 2 from b

                    k = auxthisA(kk,i)

                       !One particle terms
                       auxCIenergy= auxCIenergy + &
                            ConfigurationInteraction_instance%twoCenterIntegrals(i)%values( k, k )

                       !Two particles, same specie
                       twoParticlesEnergy=0

                       !auxIndex1 = IndexMap_tensorR2ToVectorC( k, k, numberOfSpatialOrbitals )
                       auxIndex1= ConfigurationInteraction_instance%twoIndexArray(i)%values(k,k)
                       do ll=kk+1, int( MolecularSystem_instance%species(i)%ocupationNumber )!! 1 is from a and 2 from b
                             l = auxthisA(ll,i)
                             !auxIndex2 = IndexMap_tensorR2ToVectorC( l, l, numberOfSpatialOrbitals )
                             auxIndex2= ConfigurationInteraction_instance%twoIndexArray(i)%values(l,l)
                             auxIndex = ConfigurationInteraction_instance%fourIndexArray(i)%values(auxIndex1,auxIndex2) 
                             !auxIndex = IndexMap_tensorR2ToVectorC( auxIndex1, auxIndex2, &
                             !          auxnumberOfSpatialOrbitals  )

                             !Coulomb
                             twoParticlesEnergy=twoParticlesEnergy + &
                                 ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

                             !Exchange, depends on spin
                             !if ( spin(1) .eq. spin(2) ) then

                             auxIndex = ConfigurationInteraction_instance%fourIndexArray(i)%values( &
                                           ConfigurationInteraction_instance%twoIndexArray(i)%values(k,l), &
                                           ConfigurationInteraction_instance%twoIndexArray(i)%values(l,k) )
                                            !k, l, l, k, numberOfSpatialOrbitals )

                                TwoParticlesEnergy=TwoParticlesEnergy + &
                                     kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

                             !end if
                          !end if
                       end do

                       auxCIenergy = auxCIenergy + twoParticlesEnergy

                       ! !Two particles, different species
                       if (numberOfSpecies > 1 ) then
                          do j=i+1, numberOfSpecies

                             !numberOfOtherSpecieSpatialOrbitals= ConfigurationInteraction_instance%totalNumberOfContractions ( j )
                             auxnumberOfOtherSpecieSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(j) 

                             couplingEnergy=0

                             do ll=1, int( MolecularSystem_instance%species(j)%ocupationNumber )!! 1 is from a and 2 from b
                                   l = auxthisA(ll,j)
                                   auxIndex2= ConfigurationInteraction_instance%twoIndexArray(j)%values(l,l)
                                   !auxIndex2 = IndexMap_tensorR2ToVectorC( l, & 
                                   !             l, numberOfOtherSpecieSpatialOrbitals )
                                   auxIndex = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) + auxIndex2

                                   couplingEnergy=couplingEnergy+ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1)

                             end do

                             auxCIenergy = auxCIenergy + couplingEnergy

                          end do

                       end if

                    !end if
                 end do
              end do

             !Interaction with point charges
              auxCIenergy= auxCIenergy + HartreeFock_instance%puntualInteractionEnergy

        case (1)

           !if ( allocated (differentOrbitals) ) deallocate (differentOrbitals)
           !allocate (differentOrbitals (numberOfSpecies,2 ) )
           !differentOrbitals= 0

           do i=1, numberOfSpecies

              diffOrb1 = 0
              diffOrb2 = 0

              kappa = MolecularSystem_instance%species(i)%kappa
              !numberOfSpatialOrbitals = ConfigurationInteraction_instance%totalNumberOfContractions ( i )
              !auxnumberOfSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(i) 

              !Determine different orbitals
              !call Vector_constructorInteger (differentOrbitals(i), 2)

              !differentOrbitals(i)%values = 0

              do kk=1, int( MolecularSystem_instance%species(i)%ocupationNumber )!! 1 is from a and 2 from b
                if ( abs (auxthisA(kk,i) - auxthisB(kk,i) ) > 0 ) then
                  diffOrb1= auxthisA(kk,i) 
                  diffOrb2= auxthisB(kk,i) 

                end if
              end do

              if (  diffOrb2 > 0 ) then 

                !One particle terms
                auxCIenergy= auxCIenergy +  ConfigurationInteraction_instance%twoCenterIntegrals(i)%values( &
                                  !differentOrbitals(i)%values(1), differentOrbitals(i)%values(2) )
                                  !differentOrbitals(i,1), differentOrbitals(i,2) )
                                  diffOrb1, diffOrb2 )

                 twoParticlesEnergy=0.0_8
                 !if (spin(1) .eq. spin(2) ) then

                    !auxIndex1 = IndexMap_tensorR2ToVectorC(differentOrbitals(i)%values(1), differentOrbitals(i)%values(2), &
                    !                              numberOfSpatialOrbitals )
                    auxIndex1= ConfigurationInteraction_instance%twoIndexArray(i)%values( & 
                                !differentOrbitals(i)%values(1), differentOrbitals(i)%values(2))
                                diffOrb1, diffOrb2)

                    do ll=1, int( MolecularSystem_instance%species(i)%ocupationNumber) !! 1 is from a and 2 from b
                      if ( abs (auxthisA(ll,i) - auxthisB(ll,i) ) == 0 ) then
                          l = auxthisA(ll,i) !! or b

                          !auxIndex2 = IndexMap_tensorR2ToVectorC( l, l, numberOfSpatialOrbitals )
                          auxIndex2 = ConfigurationInteraction_instance%twoIndexArray(i)%values( l,l) 

                          auxIndex = ConfigurationInteraction_instance%fourIndexArray(i)%values( auxIndex1, auxIndex2 )
                          !auxIndex = IndexMap_tensorR2ToVectorC( auxIndex1, auxIndex2, &
                          !             auxnumberOfSpatialOrbitals )

                          TwoParticlesEnergy=TwoParticlesEnergy + &
                               ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

                       end if
                    end do
                 !end if

                 !Exchange
                 do ll=1, int( MolecularSystem_instance%species(i)%ocupationNumber) !! 1 is from a and 2 from b
                    if ( abs (auxthisA(ll,i) - auxthisB(ll,i) ) == 0 ) then
                          l = auxthisA(ll,i) !! or b

                          auxIndex = ConfigurationInteraction_instance%fourIndexArray(i)%values( &
                                       ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                                         diffOrb1,l), &
                                       ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                                         l,diffOrb2) ) 

                          TwoParticlesEnergy=TwoParticlesEnergy + &
                               kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                       !end if
                    end if
                 end do

                 auxCIenergy = auxCIenergy + twoParticlesEnergy

                 ! !Two particles, different species

                 if (numberOfSpecies .gt. 1 ) then !.and. spin(1) .eq. spin(2) ) then
                    do j=1, numberOfSpecies

                       couplingEnergy=0
                       if (i .ne. j) then

                          !numberOfOtherSpecieSpatialOrbitals = ConfigurationInteraction_instance%totalNumberOfContractions ( j )
                          auxnumberOfOtherSpecieSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(j) 

                          do ll=1, int( MolecularSystem_instance%species(j)%ocupationNumber )!! 1 is from a and 2 from b
                                l = auxthisA(ll,j)

                               ! auxIndex2 = IndexMap_tensorR2ToVectorC( l, l, numberOfOtherSpecieSpatialOrbitals )
                                auxIndex2 = ConfigurationInteraction_instance%twoIndexArray(j)%values( l,l) 
                                auxIndex = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) + auxIndex2

                                couplingEnergy=couplingEnergy+ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1) 
                          end do
                          auxCIenergy = auxCIenergy + couplingEnergy

                       end if
                    end do
                 end if
              end if
              
           end do

        case (2)

           do i=1, numberOfSpecies

              numberOfOccupiedOrbitals = MolecularSystem_instance%species(i)%ocupationNumber!*lambda

              diffOrb1 = 0
              diffOrb2 = 0
              diffOrb3 = 0
              diffOrb4 = 0

              z = 1
              do k = 1, numberOfOccupiedOrbitals
                if ( z > 2 ) exit
                if ( abs(auxthisA(k,i) - auxthisB(k,i)) > 0 ) then
                  if ( z == 1 ) then
                    diffOrb1 = auxthisA(k,i)
                    diffOrb3 = auxthisB(k,i)
                  else if ( z == 2 ) then
                    diffOrb2 = auxthisA(k,i)
                    diffOrb4 = auxthisB(k,i)
                  end if 
                  z = z + 1
                end if 
              end do 

           ! !Two cases: 4 different orbitals of the same species, and 2 and 2 of different species
           !do i=1, numberOfSpecies

              kappa = MolecularSystem_instance%species(i)%kappa
              !numberOfSpatialOrbitals = ConfigurationInteraction_instance%totalNumberOfContractions ( i )
              !auxnumberOfSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(i) 

              if (  diffOrb2 > 0 ) then

                 !Coulomb
                  !! 12|34
                 !if ( spin(1) .eq. spin(3) .and. spin(2) .eq. spin(4) ) then

                    !auxIndex = IndexMap_tensorR4ToVectorC( differentOrbitals(i)%values(1), differentOrbitals(i)%values(3), &
                    !            differentOrbitals(i)%values(2), differentOrbitals(i)%values(4), numberOfSpatialOrbitals )

                     auxIndex = ConfigurationInteraction_instance%fourIndexArray(i)%values( &
                                 ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                                    diffOrb1,diffOrb3),&
                                 ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                                    diffOrb2,diffOrb4) )
 

                    auxCIenergy = auxCIenergy + &
                         ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                 !end if

                 !Exchange
                 !if ( spin(1) .eq. spin(4) .and. spin(2) .eq. spin(3) ) then
                    !auxIndex = IndexMap_tensorR4ToVectorC(  differentOrbitals(i)%values(1),  differentOrbitals(i)%values(4), &
                    !             differentOrbitals(i)%values(2),  differentOrbitals(i)%values(3), numberOfSpatialOrbitals )

                    auxIndex = ConfigurationInteraction_instance%fourIndexArray(i)%values( &
                                 ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                                    diffOrb1,diffOrb4),&
                                 ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                                    diffOrb2,diffOrb3) )
                    auxCIenergy = auxCIenergy + &
                         kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                 !end if

              end if

              !! different species

              do j=i+1, numberOfSpecies
                    !numberOfOtherSpecieSpatialOrbitals = ConfigurationInteraction_instance%totalNumberOfContractions ( j )
                    auxnumberOfOtherSpecieSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(j) 
                    numberOfOccupiedOrbitals = MolecularSystem_instance%species(j)%ocupationNumber!* &
                                                !ConfigurationInteraction_instance%lambda%values(j)
                    otherdiffOrb1 = 0
                    otherdiffOrb3 = 0

                    z = 1
                    do k = 1, numberOfOccupiedOrbitals
                      if ( z > 1 ) exit
                      if ( abs(auxthisA(k,j) - auxthisB(k,j)) > 0 ) then
                        otherdiffOrb1 = auxthisA(k,j)
                        otherdiffOrb3 = auxthisB(k,j)
                        z = z + 1
                      end if 
                    end do 


                    !if ( differentOrbitals(i)%values(3) .gt. 0 .and.  differentOrbitals(j)%values(3) .gt. 0 ) then
                    !if ( differentOrbitals(i,3) .gt. 0 .and.  differentOrbitals(j,3) .gt. 0 ) then
                    if ( diffOrb3 .gt. 0 .and. otherdiffOrb3 .gt. 0 ) then

                       !if ( spin(1) .eq. spin(3) .and. spin(2) .eq. spin(4) ) then
                          !auxIndex = IndexMap_tensorR4ToVectorC( differentOrbitals(i)%values(1), differentOrbitals(i)%values(3), &
                          !             differentOrbitals(j)%values(1), differentOrbitals(j)%values(3), &
                          !             numberOfSpatialOrbitals, numberOfOtherSpecieSpatialOrbitals )
                          auxIndex1 = ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                                          diffOrb1,diffOrb3 )
                          auxIndex2 = ConfigurationInteraction_instance%twoIndexArray(j)%values(&
                                          otherdiffOrb1,otherdiffOrb3 )
                          auxIndex = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) + auxIndex2

                          auxCIenergy = auxCIenergy + &
                               ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1)
                       !end if

                    end if
              end do
           end do

    case default

      auxCIenergy= 0.0_8

    end select
    auxCIenergy= auxCIenergy * factor

  end function ConfigurationInteraction_calculateCIenergyB



  function ConfigurationInteraction_calculateCIenergyC(auxthisA, auxthisB) result (CIenergy)
    implicit none
    type(Configuration) :: auxthisA, auxthisB
    integer :: i,j,a,b,ia,ib
    integer :: l,k,z,kk,ll
    integer :: numberOfSpecies
    integer :: numberOfSpatialOrbitals
    integer :: auxnumberOfSpatialOrbitals
    integer :: numberOfOtherSpecieSpatialOrbitals
    integer :: auxnumberOfOtherSpecieSpatialOrbitals
    integer :: lambda !occupation per orbital
    integer :: numberOfOccupiedOrbitals
    real(8) :: kappa !positive or negative exchange
    real(8) :: twoParticlesEnergy
    real(8) :: couplingEnergy
    integer :: factor
    integer :: numberOfDiffOrbitals
    integer :: auxIndex1, auxIndex2, auxIndex
    !type(Ivector), allocatable :: differentOrbitals(:)
!    integer, allocatable :: differentOrbitals(:,:)
    integer :: diffOrb1, diffOrb2, diffOrb3, diffOrb4, otherdiffOrb1,otherdiffOrb3
    !type(Ivector), allocatable :: occupiedOrbitals(:,:) !! nspecies
    !integer, allocatable :: occupiedOrbitalsA(:,:) !! nspecies
   ! integer, allocatable :: occupiedOrbitalsB(:,:) !! nspecies
    real(8) :: CIenergy
    real(8) :: auxCIenergy

    CIenergy = 0.0_8

    numberOfSpecies = MolecularSystem_instance%numberOfQuantumSpecies

    !allocate(differentOrbitals (numberOfSpecies))

    !do ia = 1, thisA%nDeterminants 
    !  do ib = 1, thisB%nDeterminants 

        auxCIenergy = 0.0_8
   
        numberOfDiffOrbitals = Configuration_checkCoincidenceB( auxthisA%occupations, auxthisB%occupations, numberOfSpecies )

        factor = 1
        if  (  numberOfDiffOrbitals == 1 .or. numberOfDiffOrbitals == 2  ) then
          call Configuration_setAtMaximumCoincidenceB( auxthisA%occupations,auxthisB%occupations, numberOfSpecies, factor )
          !call Configuration_setAtMaximumCoincidenceB( numberOfSpecies, factor )
          !factor = Configuration_setAtMaximumCoincidenceB( numberOfSpecies)
        !print *, "c", factor
        else if ( numberOfDiffOrbitals > 2 ) then
          return
        end if

        !numberOfDiffOrbitals = 3
        select case (  numberOfDiffOrbitals )

        case (0)
          
              do i=1, numberOfSpecies

                 kappa = MolecularSystem_instance%species(i)%kappa
                 numberOfSpatialOrbitals = ConfigurationInteraction_instance%totalNumberOfContractions ( i )
                 auxnumberOfSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(i) 

                 do kk=1, int( MolecularSystem_instance%species(i)%ocupationNumber) !! 1 is from a and 2 from b

                    k = auxthisA%occupations(kk,i)

                       !One particle terms
                       auxCIenergy= auxCIenergy + &
                            ConfigurationInteraction_instance%twoCenterIntegrals(i)%values( k, k )

                       !Two particles, same specie
                       twoParticlesEnergy=0

                       !auxIndex1 = IndexMap_tensorR2ToVectorC( k, k, numberOfSpatialOrbitals )
                       auxIndex1= ConfigurationInteraction_instance%twoIndexArray(i)%values(k,k)
                       do ll=kk+1, int( MolecularSystem_instance%species(i)%ocupationNumber )!! 1 is from a and 2 from b
                             l = auxthisA%occupations(ll,i)
                             !auxIndex2 = IndexMap_tensorR2ToVectorC( l, l, numberOfSpatialOrbitals )
                             auxIndex2= ConfigurationInteraction_instance%twoIndexArray(i)%values(l,l)
                             auxIndex = IndexMap_tensorR2ToVectorC( auxIndex1, auxIndex2, &
                                       auxnumberOfSpatialOrbitals  )

                             !Coulomb
                             twoParticlesEnergy=twoParticlesEnergy + &
                                 ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

                             !Exchange, depends on spin
                             !if ( spin(1) .eq. spin(2) ) then

                                auxIndex = IndexMap_tensorR2ToVectorC(&
                                           ConfigurationInteraction_instance%twoIndexArray(i)%values(k,l), &
                                           ConfigurationInteraction_instance%twoIndexArray(i)%values(l,k), &
                                           auxnumberOfSpatialOrbitals )
                                            !k, l, l, k, numberOfSpatialOrbitals )

                                TwoParticlesEnergy=TwoParticlesEnergy + &
                                     kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

                             !end if
                          !end if
                       end do

                       auxCIenergy = auxCIenergy + twoParticlesEnergy

                       ! !Two particles, different species
                       if (numberOfSpecies > 1 ) then
                          do j=i+1, numberOfSpecies

                             numberOfOtherSpecieSpatialOrbitals= ConfigurationInteraction_instance%totalNumberOfContractions ( j )
                             auxnumberOfOtherSpecieSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(j) 

                             couplingEnergy=0

                             do ll=1, int( MolecularSystem_instance%species(j)%ocupationNumber )!! 1 is from a and 2 from b
                                   l = auxthisA%occupations(ll,j)
                                   auxIndex2= ConfigurationInteraction_instance%twoIndexArray(j)%values(l,l)
                                   !auxIndex2 = IndexMap_tensorR2ToVectorC( l, & 
                                   !             l, numberOfOtherSpecieSpatialOrbitals )
                                   auxIndex = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) + auxIndex2

                                   couplingEnergy=couplingEnergy+ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1)

                             end do

                             auxCIenergy = auxCIenergy + couplingEnergy

                          end do

                       end if

                    !end if
                 end do
              end do

             !Interaction with point charges
              auxCIenergy= auxCIenergy + HartreeFock_instance%puntualInteractionEnergy

        case (1)

           !if ( allocated (differentOrbitals) ) deallocate (differentOrbitals)
           !allocate (differentOrbitals (numberOfSpecies,2 ) )
           !differentOrbitals= 0

           do i=1, numberOfSpecies

              diffOrb1 = 0
              diffOrb2 = 0

              kappa = MolecularSystem_instance%species(i)%kappa
              numberOfSpatialOrbitals = ConfigurationInteraction_instance%totalNumberOfContractions ( i )
              auxnumberOfSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(i) 

              !Determine different orbitals
              !call Vector_constructorInteger (differentOrbitals(i), 2)

              !differentOrbitals(i)%values = 0

              do kk=1, int( MolecularSystem_instance%species(i)%ocupationNumber )!! 1 is from a and 2 from b
                if ( abs (auxthisA%occupations(kk,i) - &
                     auxthisB%occupations(kk,i) ) > 0 ) then
                  !differentOrbitals(i)%values(1)= occupiedOrbitals(i,1)%values(kk) 
                  !differentOrbitals(i)%values(2)= occupiedOrbitals(i,2)%values(kk) 
                  !differentOrbitals(i,1)= auxthisA%occupations(kk,i) 
                  !differentOrbitals(i,2)= auxthisB%occupations(kk,i) 
                  diffOrb1= auxthisA%occupations(kk,i) 
                  diffOrb2= auxthisB%occupations(kk,i) 

                end if
              end do

              !call vector_show(differentOrbitals(i))
              !if (  differentOrbitals(i)%values(2) > 0 ) then !?
              !if (  differentOrbitals(i,2) > 0 ) then !?
              if (  diffOrb2 > 0 ) then !?

                !One particle terms
                auxCIenergy= auxCIenergy +  ConfigurationInteraction_instance%twoCenterIntegrals(i)%values( &
                                  !differentOrbitals(i)%values(1), differentOrbitals(i)%values(2) )
                                  !differentOrbitals(i,1), differentOrbitals(i,2) )
                                  diffOrb1, diffOrb2 )

                 twoParticlesEnergy=0.0_8
                 !if (spin(1) .eq. spin(2) ) then

                    !auxIndex1 = IndexMap_tensorR2ToVectorC(differentOrbitals(i)%values(1), differentOrbitals(i)%values(2), &
                    !                              numberOfSpatialOrbitals )
                    auxIndex1= ConfigurationInteraction_instance%twoIndexArray(i)%values( & 
                                !differentOrbitals(i)%values(1), differentOrbitals(i)%values(2))
                                diffOrb1, diffOrb2)

                    do ll=1, int( MolecularSystem_instance%species(i)%ocupationNumber) !! 1 is from a and 2 from b
                      if ( abs (auxthisA%occupations(ll,i) - &
                            auxthisB%occupations(ll,i) ) == 0 ) then
                          l = auxthisA%occupations(ll,i) !! or b

                          !auxIndex2 = IndexMap_tensorR2ToVectorC( l, l, numberOfSpatialOrbitals )
                          auxIndex2 = ConfigurationInteraction_instance%twoIndexArray(i)%values( l,l) 

                          auxIndex = IndexMap_tensorR2ToVectorC( auxIndex1, auxIndex2, &
                                       auxnumberOfSpatialOrbitals )

                          TwoParticlesEnergy=TwoParticlesEnergy + &
                               ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

                       end if
                    end do
                 !end if

                 !Exchange
                 do ll=1, int( MolecularSystem_instance%species(i)%ocupationNumber) !! 1 is from a and 2 from b
                    if ( abs (auxthisA%occupations(ll,i) -&
                         auxthisB%occupations(ll,i) ) == 0 ) then
                       l = auxthisA%occupations(ll,i) !! or b
                       !if (spin(1) .eq. spin(3) .and. spin(2) .eq. spin(3) ) then
                          !auxIndex = IndexMap_tensorR4ToVectorC( differentOrbitals(i)%values(1), l,  &
                          !             l, differentOrbitals(i)%values(2), & 
                          !             numberOfSpatialOrbitals )
                          auxIndex = IndexMap_tensorR2ToVectorC(&
                                       ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                                         !differentOrbitals(i)%values(1),l), &
                                         !differentOrbitals(i,1),l), &
                                         diffOrb1,l), &
                                       ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                                         !l,differentOrbitals(i)%values(2)), &
                                         !l,differentOrbitals(i,2)), &
                                         l,diffOrb2), &
                                       auxnumberOfSpatialOrbitals  )

                          TwoParticlesEnergy=TwoParticlesEnergy + &
                               kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                       !end if
                    end if
                 end do

                 auxCIenergy = auxCIenergy + twoParticlesEnergy

                 ! !Two particles, different species

                 if (numberOfSpecies .gt. 1 ) then !.and. spin(1) .eq. spin(2) ) then
                    do j=1, numberOfSpecies

                       couplingEnergy=0
                       if (i .ne. j) then

                          numberOfOtherSpecieSpatialOrbitals = ConfigurationInteraction_instance%totalNumberOfContractions ( j )
                          auxnumberOfOtherSpecieSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(j) 

                          do ll=1, int( MolecularSystem_instance%species(j)%ocupationNumber )!! 1 is from a and 2 from b
                                l = auxthisA%occupations(ll,j)

                               ! auxIndex2 = IndexMap_tensorR2ToVectorC( l, l, numberOfOtherSpecieSpatialOrbitals )
                                auxIndex2 = ConfigurationInteraction_instance%twoIndexArray(j)%values( l,l) 
                                auxIndex = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) + auxIndex2

                                couplingEnergy=couplingEnergy+ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1) 
                          end do
                          auxCIenergy = auxCIenergy + couplingEnergy

                       end if
                    end do
                 end if
              end if
              
           end do

        case (2)

           !if ( allocated (differentOrbitals) ) deallocate (differentOrbitals)
           !allocate (differentOrbitals (numberOfSpecies,4 ) )
           !differentOrbitals = 0


           do i=1, numberOfSpecies
              !call Vector_constructorInteger (differentOrbitals(i),4)
              !lambda=ConfigurationInteraction_instance%lambda%values(i)
              !differentOrbitals(i)%values = 0
              numberOfOccupiedOrbitals = MolecularSystem_instance%species(i)%ocupationNumber!*lambda

              diffOrb1 = 0
              diffOrb2 = 0
              diffOrb3 = 0
              diffOrb4 = 0

              z = 1
              do k = 1, numberOfOccupiedOrbitals
                if ( z > 2 ) exit
                if ( abs(auxthisA%occupations(k,i) - auxthisB%occupations(k,i)) > 0 ) then
                  if ( z == 1 ) then
                    diffOrb1 = auxthisA%occupations(k,i)
                    diffOrb3 = auxthisB%occupations(k,i)
                  else if ( z == 2 ) then
                    diffOrb2 = auxthisA%occupations(k,i)
                    diffOrb4 = auxthisB%occupations(k,i)
                  end if 
                  z = z + 1
                end if 
              end do 

              !print *, "c", diffOrb1, diffOrb2, diffOrb3, diffOrb4
           !end do
           ! !Two cases: 4 different orbitals of the same species, and 2 and 2 of different species

           !do i=1, numberOfSpecies

              kappa = MolecularSystem_instance%species(i)%kappa
              numberOfSpatialOrbitals = ConfigurationInteraction_instance%totalNumberOfContractions ( i )
              auxnumberOfSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(i) 

              if (  diffOrb2 > 0 ) then
              !if (  differentOrbitals(i,2) > 0 ) then
              !if (  differentOrbitals(i)%values(2) > 0 ) then

                 !Coulomb
                  !! 12|34
                 !if ( spin(1) .eq. spin(3) .and. spin(2) .eq. spin(4) ) then

                    !auxIndex = IndexMap_tensorR4ToVectorC( differentOrbitals(i)%values(1), differentOrbitals(i)%values(3), &
                    !            differentOrbitals(i)%values(2), differentOrbitals(i)%values(4), numberOfSpatialOrbitals )

                    auxIndex = IndexMap_tensorR2ToVectorC(&
                                 ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                                    !differentOrbitals(i)%values(1),differentOrbitals(i)%values(3)),&
                                    !differentOrbitals(i,1),differentOrbitals(i,3)),&
                                    diffOrb1,diffOrb3),&
                                 ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                                    !differentOrbitals(i)%values(2),differentOrbitals(i)%values(4)),&
                                    !differentOrbitals(i,2),differentOrbitals(i,4)),&
                                    diffOrb2,diffOrb4),&
                                 auxnumberOfSpatialOrbitals )


                    auxCIenergy = auxCIenergy + &
                         ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                 !end if

                 !Exchange
                 !if ( spin(1) .eq. spin(4) .and. spin(2) .eq. spin(3) ) then
                    !auxIndex = IndexMap_tensorR4ToVectorC(  differentOrbitals(i)%values(1),  differentOrbitals(i)%values(4), &
                    !             differentOrbitals(i)%values(2),  differentOrbitals(i)%values(3), numberOfSpatialOrbitals )
                    auxIndex = IndexMap_tensorR2ToVectorC(&
                                 ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                                    diffOrb1,diffOrb4),&
                                    !differentOrbitals(i,1),differentOrbitals(i,4)),&
                                    !differentOrbitals(i)%values(1),differentOrbitals(i)%values(4)),&
                                 ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                                    diffOrb2,diffOrb3),&
                                    !differentOrbitals(i,2),differentOrbitals(i,3)),&
                                    !differentOrbitals(i)%values(2),differentOrbitals(i)%values(3)),&
                                 auxnumberOfSpatialOrbitals  )

                    auxCIenergy = auxCIenergy + &
                         kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
                 !end if

              end if

              !! different species

              do j=i+1, numberOfSpecies
                    numberOfOtherSpecieSpatialOrbitals = ConfigurationInteraction_instance%totalNumberOfContractions ( j )
                    auxnumberOfOtherSpecieSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(j) 
                    numberOfOccupiedOrbitals = MolecularSystem_instance%species(j)%ocupationNumber!* &
                                                !ConfigurationInteraction_instance%lambda%values(j)
                    otherdiffOrb1 = 0
                    otherdiffOrb3 = 0

                    z = 1
                    do k = 1, numberOfOccupiedOrbitals
                      if ( z > 1 ) exit
                      if ( abs(auxthisA%occupations(k,j) - auxthisB%occupations(k,j)) > 0 ) then
                        otherdiffOrb1 = auxthisA%occupations(k,j)
                        otherdiffOrb3 = auxthisB%occupations(k,j)
                        z = z + 1
                      end if 
                    end do 


                    !if ( differentOrbitals(i)%values(3) .gt. 0 .and.  differentOrbitals(j)%values(3) .gt. 0 ) then
                    !if ( differentOrbitals(i,3) .gt. 0 .and.  differentOrbitals(j,3) .gt. 0 ) then
                    if ( diffOrb3 .gt. 0 .and. otherdiffOrb3 .gt. 0 ) then

                       !if ( spin(1) .eq. spin(3) .and. spin(2) .eq. spin(4) ) then
                          !auxIndex = IndexMap_tensorR4ToVectorC( differentOrbitals(i)%values(1), differentOrbitals(i)%values(3), &
                          !             differentOrbitals(j)%values(1), differentOrbitals(j)%values(3), &
                          !             numberOfSpatialOrbitals, numberOfOtherSpecieSpatialOrbitals )
                          auxIndex1 = ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                                          diffOrb1,diffOrb3 )
                                          !differentOrbitals(i,1),differentOrbitals(i,3) )
                                          !differentOrbitals(i)%values(1),differentOrbitals(i)%values(3) )
                          auxIndex2 = ConfigurationInteraction_instance%twoIndexArray(j)%values(&
                                          otherdiffOrb1,otherdiffOrb3 )
                                          !differentOrbitals(j,1),differentOrbitals(j,3) )
                                          !differentOrbitals(j)%values(1),differentOrbitals(j)%values(3) )
                          auxIndex = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) + auxIndex2

                          auxCIenergy = auxCIenergy + &
                               ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1)
                       !end if

                    end if
              end do
           end do

        case default

           auxCIenergy= 0.0_8

        end select
        CIenergy= auxCIenergy * factor + CIenergy
        !print *, "factor", factor
    !  end do  ! ib
    !end do ! ia 
!    thisA%occupations2 = thisA%occupations
!    thisB%occupations2 = thisB%occupations
    !deallocate (differentOrbitals)

  end function ConfigurationInteraction_calculateCIenergyC



  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_getTransformedIntegrals()
    implicit none

!    type(TransformIntegrals) :: repulsionTransformer
    integer :: numberOfSpecies
    integer :: i,j,m,n,mu,nu,a,b,c
    integer :: specieID
    integer :: otherSpecieID
    character(10) :: nameOfSpecie
    character(10) :: nameOfOtherSpecie
    integer :: ocupationNumber
    integer :: ocupationNumberOfOtherSpecie
    integer :: numberOfContractions
    integer :: numberOfContractionsOfOtherSpecie
!    type(Matrix) :: auxMatrix
!    type(Matrix) :: molecularCouplingMatrix
!    type(Matrix) :: molecularExtPotentialMatrix
!    type(Matrix) :: couplingMatrix
    type(Matrix) :: hcoreMatrix
    type(Matrix) :: coefficients
    real(8) :: charge
    real(8) :: otherSpecieCharge

    integer :: ssize1, ssize2
    type(Matrix) :: externalPotential

    character(50) :: wfnFile
    character(50) :: arguments(20)
    integer :: wfnUnit

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    allocate(ConfigurationInteraction_instance%twoCenterIntegrals(numberOfSpecies))
    !allocate(ConfigurationInteraction_instance%FockMatrix(numberOfSpecies))
    allocate(ConfigurationInteraction_instance%fourCenterIntegrals(numberOfSpecies,numberOfSpecies))
    !allocate(ConfigurationInteraction_instance%energyofmolecularorbitals(numberOfSpecies))

    allocate(ConfigurationInteraction_instance%twoIndexArray(numberOfSpecies))
    allocate(ConfigurationInteraction_instance%fourIndexArray(numberOfSpecies))

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
        charge=MolecularSystem_getCharge(i)

!        write (6,"(T10,A)")"ONE PARTICLE INTEGRALS TRANSFORMATION FOR: "//trim(nameOfSpecie)
        call Matrix_constructor (ConfigurationInteraction_instance%twoCenterIntegrals(i), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )

        !call Matrix_constructor (ConfigurationInteraction_instance%FockMatrix(i), int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )
        call Matrix_constructor (hcoreMatrix,int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)

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

        !arguments(1) = "ORBITALS"
        !call Vector_getFromFile( elementsNum = numberOfContractions, &
        !  unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
        !  output =ConfigurationInteraction_instance%energyofmolecularorbitals(i) )     

        !do m=1,numberOfContractions
        !   ConfigurationInteraction_instance%fockMatrix(i)%values(m,m) = &
        !        ConfigurationInteraction_instance%energyofmolecularorbitals(i)%values(m) 
        !end do

        if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
          arguments(1) = "EXTERNAL_POTENTIAL"

          externalPotential = &
            Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
            columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

          hcoreMatrix%values = hcoreMatrix%values + externalPotential%values
        end if
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

        call Matrix_constructorInteger(ConfigurationInteraction_instance%twoIndexArray(i), &
                          int( numberOfContractions,8), int( numberOfContractions,8) , 0 )

       c = 0
       do a=1,numberOfContractions
         do b=a, numberOfContractions
           c = c + 1
           ConfigurationInteraction_instance%twoIndexArray(i)%values(a,b) = c !IndexMap_tensorR2ToVectorC( a, b, numberOfContractions )
           ConfigurationInteraction_instance%twoIndexArray(i)%values(b,a) = ConfigurationInteraction_instance%twoIndexArray(i)%values(a,b)
         end do 
       end do


          ssize1 = MolecularSystem_getTotalNumberOfContractions( i )
          ssize1 = ( ssize1 * ( ssize1 + 1 ) ) / 2

          call Matrix_constructorInteger(ConfigurationInteraction_instance%fourIndexArray(i), &
                          int( ssize1,8), int( ssize1,8) , 0 )
          c = 0
          do a=1, ssize1
            do b=a, ssize1
              c = c + 1
              !print *, a,b, c
              ConfigurationInteraction_instance%fourIndexArray(i)%values(a,b) = c! IndexMap_tensorR2ToVectorC( a, b, numberOfContractions )
              ConfigurationInteraction_instance%fourIndexArray(i)%values(b,a) = &
                ConfigurationInteraction_instance%fourIndexArray(i)%values(a,b)
             end do 
           end do


        call ReadTransformedIntegrals_readOneSpecies( specieID, ConfigurationInteraction_instance%fourCenterIntegrals(i,i)   )
        ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values = &
           ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values * charge * charge

       if ( numberOfSpecies > 1 ) then
          do j = 1 , numberOfSpecies
             if ( i .ne. j) then
                nameOfOtherSpecie= trim(  MolecularSystem_getNameOfSpecie( j ) )
                otherSpecieID =MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecie )
                ocupationNumberOfOtherSpecie = MolecularSystem_getOcupationNumber( j )
                numberOfContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( j )
                otherSpecieCharge=MolecularSystem_getCharge(j)

                 call ReadTransformedIntegrals_readTwoSpecies( specieID, otherSpecieID, &
                         ConfigurationInteraction_instance%fourCenterIntegrals(i,j) )
                 ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values = &
                   ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values * charge * otherSpeciecharge


             end if
          end do
       end if
    end do
    close (wfnUnit)
    call Matrix_destructor (hcoreMatrix)
!   call Matrix_destructor (couplingMatrix)

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

!  subroutine ConfigurationInteraction_printGeometryToFile(unit)
!    implicit none
!    integer :: unit
!
!    integer :: i
!    integer :: from, to
!    real(8) :: origin(3)
!    character(50) :: auxString
!
!    
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

!  end subroutine ConfigurationInteraction_printGeometryToFile


!  subroutine ConfigurationInteraction_printBasisSetToFile(unit)
!    implicit none
!
!    integer :: unit
!
!    integer :: i, j
!    character(16) :: auxString
!
!
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
    
!  end subroutine ConfigurationInteraction_printBasisSetToFile


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

  subroutine ConfigurationInteraction_saveEigenVector () 
    implicit none
    character(50) :: nameFile
    integer :: unitFile
    integer(8) :: i, ia
    integer :: ib, nonzero
    integer, allocatable :: auxIndexArray(:)
    real(8), allocatable :: auxArray(:)
    integer :: maxStackSize

    maxStackSize = CONTROL_instance%CI_STACK_SIZE 
    nameFile = "lowdin.civec"
    unitFile = 20

    nonzero = 0
    do i = 1, ConfigurationInteraction_instance%numberOfConfigurations
      if ( abs(ConfigurationInteraction_instance%eigenVectors%values(i,1) ) >= 1E-8 ) nonzero = nonzero + 1
    end do 

    print *, "nonzero", nonzero

    allocate(auxArray(nonzero))
    allocate(auxIndexArray(nonzero))

    ia = 0
    do i = 1, ConfigurationInteraction_instance%numberOfConfigurations
      if ( abs(ConfigurationInteraction_instance%eigenVectors%values(i,1) ) >= 1E-8 ) then 
        ia = ia + 1
        auxIndexArray(ia) = i 
        auxArray(ia) = ConfigurationInteraction_instance%eigenVectors%values(i,1) 
      end if
    end do 

    open(unit=unitFile, file=trim(nameFile), status="replace", form="unformatted")

    write(unitFile) ConfigurationInteraction_instance%eigenValues%values(1)
    write(unitFile) nonzero

    do i = 1, ceiling(real(nonzero) / real(maxStackSize) )
      ib = maxStackSize * i  
      ia = ib - maxStackSize + 1
      if ( ib > nonzero ) ib = nonzero
      write(unitFile) auxIndexArray(ia:ib)
    end do
    deallocate(auxIndexArray)

    do i = 1, ceiling(real(nonzero) / real(maxStackSize) )
      ib = maxStackSize * i  
      ia = ib - maxStackSize + 1
      if ( ib > nonzero ) ib = nonzero
      write(unitFile) auxArray(ia:ib)
    end do
    deallocate(auxArray)

    close(unitFile)

  end subroutine ConfigurationInteraction_saveEigenVector

  subroutine ConfigurationInteraction_loadEigenVector (eigenValues,eigenVectors) 
    implicit none
    type(Vector8) :: eigenValues
    type(Matrix) :: eigenVectors
    character(50) :: nameFile
    integer :: unitFile
    integer :: i, ia, ib, nonzero
    real(8) :: eigenValue
    integer, allocatable :: auxIndexArray(:)
    real(8), allocatable :: auxArray(:)
    integer :: maxStackSize

    maxStackSize = CONTROL_instance%CI_STACK_SIZE 
 

    nameFile = "lowdin.civec"
    unitFile = 20


    open(unit=unitFile, file=trim(nameFile), status="old", action="read", form="unformatted")

    readvectors : do
      read (unitFile) eigenValue
      read (unitFile) nonzero
      print *, "eigenValue", eigenValue
      print *, "nonzero", nonzero

      allocate (auxIndexArray(nonzero))
      auxIndexArray = 0

      do i = 1, ceiling(real(nonZero) / real(maxStackSize) )
        ib = maxStackSize * i  
        ia = ib - maxStackSize + 1
        if ( ib >  nonZero ) ib = nonZero
       read (unitFile) auxIndexArray(ia:ib)
      end do

      allocate (auxArray(nonzero))
      auxArray = 0

      do i = 1, ceiling(real(nonZero) / real(maxStackSize) )
        ib = maxStackSize * i  
        ia = ib - maxStackSize + 1
        if ( ib >  nonZero ) ib = nonZero
       read (unitFile) auxArray(ia:ib)
      end do
      exit readvectors
    end do readvectors

    eigenValues%values(1) = eigenValue
    do i = 1, nonzero
      eigenVectors%values(auxIndexArray(i),1) = auxArray(i)
    end do

    deallocate (auxIndexArray )
    deallocate (auxArray )


    close(unitFile)

  end subroutine ConfigurationInteraction_loadEigenVector


  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_diagonalize(maxn,ldv, maxnev, maxncv, eigenValues, eigenVectors)
    implicit none

  !*******************************************************************************
  !
  !! SSSIMP is a simple program to call ARPACK for a symmetric eigenproblem.
  !
  !    This code shows how to use ARPACK to find a few eigenvalues
  !    LAMBDA and corresponding eigenvectors X for the standard
  !    eigenvalue problem:
  !
  !      A * X = LAMBDA * X
  !
  !    where A is an N by N real symmetric matrix.
  !
  !    The only things that must be supplied in order to use this
  !    routine on your problem is:
  !
  !    * to change the array dimensions appropriately, 
  !    * to specify WHICH eigenvalues you want to compute
  !    * to supply a matrix-vector product
  !      w <- A * v
  !      in place of the call to AV( ) below.
  !
  !  Example by:
  !
  !    Richard Lehoucq, Danny Sorensen, Chao Yang,
  !    Department of Computational and Applied Mathematics,
  !    Rice University,
  !    Houston, Texas.
  !
  !  Storage:
  ! 
  !    The maximum dimensions for all arrays are set here to accommodate 
  !    a problem size of N <= MAXN
  !
  !    NEV is the number of eigenvalues requested.
  !    See specifications for ARPACK usage below.
  !
  !    NCV is the largest number of basis vectors that will be used in 
  !    the Implicitly Restarted Arnoldi Process.  Work per major iteration is
  !    proportional to N*NCV*NCV.
  !
  !    You must set: 
  ! 
  !    MAXN:   Maximum dimension of the A allowed. 
  !    MAXNEV: Maximum NEV allowed. 
  !    MAXNCV: Maximum NCV allowed. 

    integer(8) :: maxn 
    integer :: maxnev 
    integer :: maxncv 
    integer(8) :: ldv 
    integer :: iter
  
!    intrinsic abs
    character(1) bmat  
    character ( len = 2 ) which
    integer ido,ierr,info,iparam(11),ipntr(11),ishfts,j,lworkl,maxitr,mode1,n,nconv,ncv,nev,nx
    logical rvec
    external saxpy
    real(8) sigma
    real(8) tol
    real(8), parameter :: zero = 1E-08

!    real(8) :: v(ldv,maxncv) 
!    real(8) :: ax(maxn)
!    real(8), external :: snrm2

    !! arrays
!    real(8) d(maxnev)
!    real(8) :: resid(maxn)
!    logical select(maxncv)
!    real(8) :: z(ldv,maxnev) 
!    real(8) workl(maxncv*(maxncv+8))
!    real(8) :: workd(3*maxn)
    real(8), allocatable :: v(:,:) 
    real(8), allocatable :: d(:)
    real(8), allocatable :: resid(:),residi(:)
    real(8), allocatable :: z(:,:) 
    real(8), allocatable :: workl(:)
    real(8), allocatable :: workd(:)
    logical, allocatable :: select(:)

    type(Vector8), intent(inout) :: eigenValues
    type(Matrix), intent(inout) :: eigenVectors
    integer :: ii, jj, ia

    if (allocated(v) ) deallocate(v)
    allocate (v(ldv,maxncv))
    if (allocated(d) ) deallocate(d)
    allocate (d(maxnev))
    if (allocated(resid) ) deallocate(resid)
    allocate (resid(maxn))
    if (allocated(residi) ) deallocate(residi)
    allocate (residi(maxn))

    if (allocated(z) ) deallocate(z)
    allocate (z(ldv,maxnev))
    if (allocated(workl) ) deallocate(workl)
    allocate (workl(maxncv*(maxncv+8)))
    if (allocated(workd) ) deallocate(workd)
    allocate (workd(3*maxn))
    if (allocated(select) ) deallocate(select)
    allocate (select(maxncv))

    v = 0.0_8
    d = 0.0_8
    resid = 0.0_8
    residi = 0.0_8
    z = 0.0_8
    workl = 0.0_8
    workd = 0.0_8

  !
  !  The following include statement and assignments control trace output 
  !  from the internal actions of ARPACK.  See debug.doc in the
  !  DOCUMENTS directory for usage.  
  !
  !  Initially, the most useful information will be a breakdown of
  !  time spent in the various stages of computation given by setting 
  !  msaupd = 1.
  !
  
  !  ndigit = -3
  !  logfil = 6
  !  msgets = 0
  !  msaitr = 0
  !  msapps = 0
  !  msaupd = 1
  !  msaup2 = 0
  !  mseigt = 0
  !  mseupd = 0
  !
  !  Set dimensions for this problem.
  !
    nx = ConfigurationInteraction_instance%numberOfConfigurations
    !n = nx * nx  
    !print *, "nnn", n

  !
  !  Specifications for ARPACK usage are set below:
  !
  !  1) NEV = 4 asks for 4 eigenvalues to be computed.                            !
  !  2) NCV = 20 sets the length of the Arnoldi factorization.
  !
  !  3) This is a standard problem(indicated by bmat  = 'I')
  !
  !  4) Ask for the NEV eigenvalues of largest magnitude
  !     (indicated by which = 'LM')
  !
  !  See documentation in SSAUPD for the other options SM, LA, SA, LI, SI.
  !
  !  NEV and NCV must satisfy the following conditions:
  !
  !    NEV <= MAXNEV
  !    NEV + 1 <= NCV <= MAXNCV
  !
    nev = CONTROL_instance%NUMBER_OF_CI_STATES 
    ncv = maxncv !! 
    bmat = 'I'
    which = 'SA'

!
!  WHICH   Character*2.  (INPUT)
!          Specify which of the Ritz values of OP to compute.
!
!          'LA' - compute the NEV largest (algebraic) eigenvalues.
!          'SA' - compute the NEV smallest (algebraic) eigenvalues.
!          'LM' - compute the NEV largest (in magnitude) eigenvalues.
!          'SM' - compute the NEV smallest (in magnitude) eigenvalues. 
!          'BE' - compute NEV eigenvalues, half from each end of the
!                 spectrum.  When NEV is odd, compute one more from the
!                 high end than from the low end.
!           (see remark 1 below)

  
!    if ( maxn < n ) then
!      write ( *, '(a)' ) ' '
!      write ( *, '(a)' ) 'SSSIMP - Fatal error!'
!      write ( *, '(a)' ) '  N is greater than MAXN '
!      !stop
!    else if ( maxnev < nev ) then
!      write ( *, '(a)' ) ' '
!      write ( *, '(a)' ) 'SSSIMP - Fatal error!'
!      write ( *, '(a)' ) '  NEV is greater than MAXNEV '
!      stop
!    else if ( maxncv < ncv ) then
!      write ( *, '(a)' ) ' '
!      write ( *, '(a)' ) 'SSSIMP - Fatal error!'
!      write ( *, '(a)' ) '  NCV is greater than MAXNCV '
!      stop
!    end if
  !
  !  TOL determines the stopping criterion.  Expect
  !    abs(lambdaC - lambdaT) < TOL*abs(lambdaC)
  !  computed   true
  !  If TOL <= 0, then TOL <- macheps (machine precision) is used.
  !
  !  IDO is the REVERSE COMMUNICATION parameter
  !  used to specify actions to be taken on return
  !  from SSAUPD. (See usage below.)
  !  It MUST initially be set to 0 before the first
  !  call to SSAUPD.
  !
  !  INFO on entry specifies starting vector information
  !  and on return indicates error codes
  !  Initially, setting INFO=0 indicates that a
  !  random starting vector is requested to 
  !  start the ARNOLDI iteration.  Setting INFO to
  !  a nonzero value on the initial call is used 
  !  if you want to specify your own starting 
  !  vector. (This vector must be placed in RESID.)
  !
  !  The work array WORKL is used in SSAUPD as workspace.  Its dimension
  !  LWORKL is set as illustrated below. 
  !
    lworkl = ncv * ( ncv + 8 )
    !tol = zero
    TOL = CONTROL_instance%CI_CONVERGENCE !1.0d-4 !    tolerance for the eigenvector residual
    ido = 0
  !
  !  Specification of Algorithm Mode:
  !
  !  This program uses the exact shift strategy
  !  (indicated by setting PARAM(1) = 1).
  !
  !  IPARAM(3) specifies the maximum number of Arnoldi iterations allowed.  
  !
  !  Mode 1 of SSAUPD is used (IPARAM(7) = 1). 
  !
  !  All these options can be changed by the user.  For details see the
  !  documentation in SSAUPD.
  !
    ishfts = 1
    maxitr = 300
    mode1 = 1
  
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode1
  !
  !  MAIN LOOP (Reverse communication loop)
  !
  !  Repeatedly call SSAUPD and take actions indicated by parameter 
  !  IDO until convergence is indicated or MAXITR is exceeded.
  !

    !call omp_set_num_threads(omp_get_max_threads())

    !! starting vector
    resid = 0
    info = 1

    do ii = 1, CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX 
      resid(ConfigurationInteraction_instance%auxIndexCIMatrix%values(ii)) = ConfigurationInteraction_instance%initialEigenVectors%values(ii,1)
      residi(ConfigurationInteraction_instance%auxIndexCIMatrix%values(ii)) = ConfigurationInteraction_instance%initialEigenVectors%values(ii,1)
    end do

    print *, ConfigurationInteraction_instance%initialEigenValues%values(1)
    iter = 0
    do
  
      call dsaupd ( ido, bmat, nx, which, nev, tol, resid, &
        maxncv, v, int(ldv,4), iparam, ipntr, workd, workl, &
        lworkl, info )
      iter = iter + 1
      if ( ido /= -1 .and. ido /= 1 ) then
        exit
      end if
  !
  !  Perform matrix vector multiplication
  !
  !    y <--- OP*x
  !
  !  The user supplies a matrix-vector multiplication routine that takes
  !  workd(ipntr(1)) as the input, and return the result to workd(ipntr(2)).
  !
      call av ( nx, workd(ipntr(1)), workd(ipntr(2)) )
      !call matvec ( nx, residi,workd(ipntr(1)), workd(ipntr(2)), iter )
  
     end do
  !
  !  Either we have convergence or there is an error.
  !
    if ( info < 0 ) then
  !
  !  Error message. Check the documentation in SSAUPD.
  !
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SSSIMP - Fatal error!'
      write ( *, '(a,i6)' ) '  Error with DSAUPD, INFO = ', info
      write ( *, '(a)' ) '  Check documentation in SSAUPD.'
  
    else
  !
  !  No fatal errors occurred.
  !  Post-Process using SSEUPD.
  !
  !  Computed eigenvalues may be extracted.
  !
  !  Eigenvectors may be also computed now if
  !  desired.  (indicated by rvec = .true.)
  !
  !  The routine SSEUPD now called to do this
  !  post processing (Other modes may require
  !  more complicated post processing than mode1.)
  !
      rvec = .true.
  
      call dseupd ( rvec, 'A', select, d, z,int( ldv,4), sigma, &
        bmat, nx, which, nev, tol, resid, ncv, v, int(ldv,4), &
        iparam, ipntr, workd, workl, lworkl, ierr )
  !
  !  Eigenvalues are returned in the first column of the two dimensional 
  !  array D and the corresponding eigenvectors are returned in the first 
  !  NCONV (=IPARAM(5)) columns of the two dimensional array V if requested.
  !
  !  Otherwise, an orthogonal basis for the invariant subspace corresponding 
  !  to the eigenvalues in D is returned in V.
  !
   !! Saving the eigenvalues !! Saving the eigenvectors

      do ii = 1, maxnev
        eigenValues%values(ii) = d(ii)
      end do

      ia = 0
      do ii = 1, nx
        do jj = 1, CONTROL_instance%NUMBER_OF_CI_STATES
          eigenVectors%values(ii,jj) = z(ii,jj)
          if ( abs(z(ii,jj)) > 1E-6) ia = ia + 1 
        end do
      end do 

      if ( ierr /= 0 ) then
  
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SSSIMP - Fatal error!'
        write ( *, '(a,i6)' ) '  Error with SSEUPD, IERR = ', ierr
        write ( *, '(a)' ) '  Check the documentation of SSEUPD.'
  !
  !  Compute the residual norm
  !
  !    ||  A*x - lambda*x ||
  ! 
  !  for the NCONV accurately computed eigenvalues and 
  !  eigenvectors.  (iparam(5) indicates how many are 
  !  accurate to the requested tolerance)
  !
      else
  
        nconv =  iparam(5)
!  
!        do j = 1, nconv
!          call av ( nx, v(1,j), ax )
!          call saxpy ( nx, -d(j,1), v(1,j), 1, ax, 1 )
!          d(j,2) = snrm2 ( nx, ax, 1)
!          d(j,2) = d(j,2) / abs ( d(j,1) )
!        end do
!  !
!  !  Display computed residuals.
!  !
!        call smout ( 6, nconv, 2, d, maxncv, -6, &
!          'Ritz values and relative residuals' )
  
      end if
  !
  !  Print additional convergence information.
  !
      if ( info == 1) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Maximum number of iterations reached.'
      else if ( info == 3) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  No shifts could be applied during implicit' &
          // ' Arnoldi update, try increasing NCV.'
      end if
  
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '====== '
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Size of the matrix is ', nx
      write ( *, '(a,i6)' ) '  The number of Ritz values requested is ', nev
      write ( *, '(a,i6)' ) &
        '  The number of Arnoldi vectors generated (NCV) is ', ncv
      write ( *, '(a)' ) '  What portion of the spectrum: ' // which
      write ( *, '(a,i6)' ) &
        '  The number of converged Ritz values is ', nconv
      write ( *, '(a,i8)' ) &
        '  The number of Implicit Arnoldi update iterations taken is ', iparam(3)
      write ( *, '(a,i6)' ) '  The number of OP*x is ', iparam(9)
     write ( *, '(a,g14.6)' ) '  The convergence criterion is ', tol
 
    end if
  
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Normal end of execution.'
  
    write ( *, '(a)' ) ' '
    !call timestamp ( )
  
  end subroutine ConfigurationInteraction_diagonalize

  subroutine av ( nx, v, w)
  
  !*******************************************************************************
  !! AV computes w <- A * V where A is a discretized Laplacian.
  !  Parameters:
  !    Input, integer NX, the length of the vectors.
  !    Input, real V(NX), the vector to be operated on by A.
  !    Output, real W(NX), the result of A*V.
  !
    implicit none
  
    integer nx
    real(8) v(nx)
    real(8) w(nx)
    character(50) :: CIFile
    integer :: CIUnit
    integer, allocatable :: jj(:)
    real(8), allocatable :: CIEnergy(:)
    integer :: nonzero,ii, kk
    integer :: maxStackSize, i, ia, ib

    CIFile = "lowdin.ci"
    CIUnit = 20
    nonzero = 0
    maxStackSize = CONTROL_instance%CI_STACK_SIZE 

    w = 0
#ifdef intel
    open(unit=CIUnit, file=trim(CIFile), action = "read", form="unformatted", BUFFERED="YES")
#else
    open(unit=CIUnit, file=trim(CIFile), action = "read", form="unformatted")
#endif

    readmatrix : do
      read (CIUnit) nonzero
      if (nonzero > 0 ) then

        read (CIUnit) ii

        if ( allocated(jj)) deallocate (jj)
        allocate (jj(nonzero))
        jj = 0

        if ( allocated(CIEnergy)) deallocate (CIEnergy)
        allocate (CIEnergy(nonzero))
        CIEnergy = 0

        do i = 1, ceiling(real(nonZero) / real(maxStackSize) )

          ib = maxStackSize * i  
          ia = ib - maxStackSize + 1
          if ( ib >  nonZero ) ib = nonZero
          read (CIUnit) jj(ia:ib)
    
        end do

        do i = 1, ceiling(real(nonZero) / real(maxStackSize) )

          ib = maxStackSize * i  
          ia = ib - maxStackSize + 1
          if ( ib >  nonZero ) ib = nonZero
          read (CIUnit) CIEnergy(ia:ib)
    
        end do

        w(ii) = w(ii) + CIEnergy(1)*v(jj(1)) !! disk
        do kk = 2, nonzero
          !w(ii) = w(ii) + ConfigurationInteraction_calculateCIenergy(ii,jj(kk))*v(jj(kk))  !! direct
          w(ii) = w(ii) + CIEnergy(kk)*v(jj(kk)) !! disk
          w(jj(kk)) = w(jj(kk)) + CIEnergy(kk)*v(ii) !! disk
        end do

      else if ( nonzero == -1 ) then
        exit readmatrix
      end if
    end do readmatrix

!! memory
!    do i = 1, nx
!        w(:) = w(:) + ConfigurationInteraction_instance%hamiltonianMatrix%values(:,i)*v(i)
!    end do 

     close(CIUnit)

    return
  end subroutine av


  subroutine ConfigurationInteraction_jadamiluInterface(n,  maxeig, eigenValues, eigenVectors)
    implicit none
    integer(8) :: maxnev
    integer(8) :: a
    real(8) :: CIenergy
    integer(8) :: nproc
    type(Vector8), intent(inout) :: eigenValues
    type(Matrix), intent(inout) :: eigenVectors
    type (Vector8) :: diagonalHamiltonianMatrix

!   N: size of the problem
!   MAXEIG: max. number of wanteg eig (NEIG<=MAXEIG)
!   MAXSP: max. value of MADSPACE
    integer(8) :: n, maxeig, MAXSP
    integer(8) :: LX
    real(8), allocatable :: EIGS(:), RES(:), X(:)!, D(:)
!   arguments to pass to the routines
    integer(8) :: NEIG, MADSPACE, ISEARCH, NINIT
    integer(8) :: ICNTL(5)
    integer(8) :: ITER, IPRINT, INFO
    real(8) :: SIGMA, TOL, GAP, MEM, DROPTOL, SHIFT
    integer(8) :: NDX1, NDX2, NDX3
    integer(8) :: IJOB!   some local variables
    integer(8) :: I,J,K
    integer(4) :: iiter

    maxsp = 20
    LX = N*(3*MAXSP+MAXEIG+1)+4*MAXSP*MAXSP

    if ( allocated ( eigs ) ) deallocate ( eigs )
    allocate ( eigs ( maxeig ) )
    eigs = 0.0_8
    if ( allocated ( res ) ) deallocate ( res )
    allocate ( res ( maxeig ) )
    res = 0.0_8
    if ( allocated ( x ) ) deallocate ( x )
    allocate ( x ( lx ) )
    x = 0.0_8

!    set input variables
!    the matrix is already in the required format

     IPRINT = -6 !     standard report on standard output
     ISEARCH = 1 !    we want the smallest eigenvalues
     NEIG = maxeig !    number of wanted eigenvalues
     !NINIT = 0 !    no initial approximate eigenvectors
     NINIT = 1 !    initial approximate eigenvectors
     MADSPACE = maxsp !    desired size of the search space
     ITER = 100 !    maximum number of iteration steps
     TOL = CONTROL_instance%CI_CONVERGENCE !1.0d-4 !    tolerance for the eigenvector residual

      NDX1 = 0
      NDX2 = 0
      MEM = 0

!    additional parameters set to default
     ICNTL(1)=0
     ICNTL(2)=0
     ICNTL(3)=0
     ICNTL(4)=0
     ICNTL(5)=1

    DROPTOL = 1E-4

     IJOB=0

     ! set initial eigenpairs
     if ( CONTROL_instance%CI_LOAD_EIGENVECTOR ) then 
       do i = 1, CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX 
         X(ConfigurationInteraction_instance%auxIndexCIMatrix%values(i)) = eigenVectors%values(i,1)
       end do

       do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
         EIGS(i) = eigenValues%values(i)
        print *,  eigenValues%values(i), EIGS(I)

       end do
     else
       do i = 1, CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX 
         X(ConfigurationInteraction_instance%auxIndexCIMatrix%values(i)) = ConfigurationInteraction_instance%initialEigenVectors%values(i,1)
       end do

       do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
         EIGS(i) = ConfigurationInteraction_instance%initialEigenValues%values(i)
       end do
     end if

     SIGMA = EIGS(1)
     gap= 0
     SHIFT = EIGS(1)
      print *, "Eigenvalue(1)", eigs(1), "Eigenvector(1)", x(1)
      iiter = 0
10   CALL DPJDREVCOM( N, ConfigurationInteraction_instance%diagonalHamiltonianMatrix%values ,-1_8,-1_8,EIGS, RES, X, LX, NEIG, &
                      SIGMA, ISEARCH, NINIT, MADSPACE, ITER, TOL, &
                      SHIFT, DROPTOL, MEM, ICNTL, &
                      IJOB, NDX1, NDX2, IPRINT, INFO, GAP)
!    your private matrix-vector multiplication

      iiter = iiter +1
      IF (IJOB.EQ.1) THEN
!       X(NDX1) input,  X(NDX2) output
         !call matvec(N,X(1),X(NDX1),X(NDX2),iiter)
         call matvec(N,X(1),X(NDX1),X(NDX2),iiter)
         GOTO 10
      END IF

!    release internal memory and discard preconditioner
      CALL PJDCLEANUP

      !! saving the eigenvalues
      eigenValues%values = EIGS

      !! saving the eigenvectors
      k = 0
      do j = 1, maxeig
         do i = 1, N
          k = k + 1
          eigenVectors%values(i,j) = X(k)
        end do
      end do
      

  end subroutine ConfigurationInteraction_jadamiluInterface

  subroutine matvec ( nx, y, v, w, iter)
  
  !*******************************************************************************
  !! AV computes w <- A * V where A is a discretized Laplacian.
  !  Parameters:
  !    Input, integer NX, the length of the vectors.
  !    Input, real V(NX), the vector to be operated on by A.
  !    Output, real W(NX), the result of A*V.
  !
    implicit none
  
    integer(8) nx
    real(8) y(nx)
    real(8) v(nx)
    real(8) w(nx)
    real(8) :: CIEnergy
    integer(8) :: nonzero
    integer(8) :: i, j, ia, ib, ii, jj
    integer(4) :: nproc
    real(8) :: wi
    real(8) :: timeA, timeB
    type(Configuration) :: auxConfigurationI, auxConfigurationJ
    real(8) :: tol
    integer(4) :: iter, size1, size2
    integer(8), allocatable :: indexArray(:)
    tol = 1E-8

    call Configuration_copyConstructor ( ConfigurationInteraction_instance%configurations(1), auxConfigurationI )
    call Configuration_copyConstructor ( ConfigurationInteraction_instance%configurations(1), auxConfigurationJ )

    w = 0
!! memory
!    do i = 1, nx
!        w(:) = w(:) + ConfigurationInteraction_instance%hamiltonianMatrix%values(:,i)*v(i)
!    end do 

    nproc = omp_get_max_threads()

    call omp_set_num_threads(omp_get_max_threads())
    call omp_set_num_threads(nproc)

    CIenergy = 0.0_8
    nonzero = 0

      do i = 1 , nx
        if ( abs(y(i)+v(i) ) > tol) nonzero = nonzero + 1
      end do
  
      allocate(indexArray(nonzero))
      indexArray = 0
  
      ia = 0
      do i = 1 , nx
        if ( abs(y(i)+v(i) ) > tol) then
          ia = ia + 1
          indexArray(ia) = i
        end if
      end do

    ib = 0

    size1 = size(ConfigurationInteraction_instance%configurations(1)%occupations, dim=1)
    size2 = size(ConfigurationInteraction_instance%configurations(1)%occupations, dim=2) 


    timeA = omp_get_wtime()
    if ( iter == -1 ) then
    do i = 1, nx
  
      w(i) = w(i) + ConfigurationInteraction_instance%diagonalHamiltonianMatrix%values(i)*v(i)  !! direct
      wi = 0
!$omp parallel &
!$omp& private(j,CIEnergy),&
!$omp& private(auxConfigurationI ),&
!$omp& private(auxConfigurationJ ),&
!$omp& shared(i,ConfigurationInteraction_instance, HartreeFock_instance,v,nx,w,y,size1,size2) reduction (+:wi)
!$omp do 
        do j = i+1 , nx
          if ( abs(v(i) ) > tol .or. abs(v(j)) > tol ) then
            !auxConfigurationI%occupations = ConfigurationInteraction_instance%configurations(i)%occupations
            !auxConfigurationJ%occupations = ConfigurationInteraction_instance%configurations(j)%occupations
            !CIenergy = ConfigurationInteraction_calculateCIenergyC( & 
            !        auxConfigurationI, auxConfigurationJ )
            CIenergy = ConfigurationInteraction_calculateCIenergyB( & 
                    ConfigurationInteraction_instance%configurations(i)%occupations, &
                    ConfigurationInteraction_instance%configurations(j)%occupations, &
                    size1, size2  )

            w(j) = w(j) + CIEnergy*v(i)  !! direct
            wi = wi + CIEnergy*v(j)  !! direct
            ib = ib + 1
          end if
        end do 
!$omp end do nowait
!$omp end parallel
        w(i) = w(i) + wi
      end do 
else
    do i = 1, nonzero
      ii = indexArray(i)
!$omp parallel &
!$omp& private(j,jj,CIEnergy),&
!$omp& shared(i,ConfigurationInteraction_instance, HartreeFock_instance,v,nx,w,size1,size2) 
!$omp do 
         do j = 1 , nx
          jj = j

            CIenergy = ConfigurationInteraction_calculateCIenergyB( & 
                    ConfigurationInteraction_instance%configurations(ii)%occupations, &
                    ConfigurationInteraction_instance%configurations(jj)%occupations, &
                    size1, size2  )

            w(jj) = w(jj) + CIEnergy*v(ii)  !! direct
        end do 
!$omp end do nowait
!$omp end parallel
end do 
end if

      deallocate(indexArray)

      !print *, "ia ib ",ia,ib
      timeB = omp_get_wtime()
      write(*,"(A,I2,A,E10.3,A2,I8)") "  ", iter, "  ", timeB - timeA ,"  ", ia
      !w = w2

    return

  end subroutine matvec

end module ConfigurationInteraction_
