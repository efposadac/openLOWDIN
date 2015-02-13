!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief This module contains the common parameters for ALL programs (from APMO)
!! @author E. F. Posada, 2013
!! @version 1.0
module CONTROL_
  use Exception_
  use omp_lib
  implicit none

  type, public :: CONTROL
     !!***************************************************************************
     !! Parameter to control Integrals library
     !!
     real(8) :: TV
     real(8) :: INTEGRAL_THRESHOLD
     integer :: INTEGRAL_STACK_SIZE
     character(20) :: INTEGRAL_DESTINY
     character(20) :: INTEGRAL_SCHEME

     !!***************************************************************************
     !! Parameter to control SCF program
     !!
     real(8) :: SCF_NONELECTRONIC_ENERGY_TOLERANCE
     real(8) :: SCF_ELECTRONIC_ENERGY_TOLERANCE
     real(8) :: NONELECTRONIC_DENSITY_MATRIX_TOLERANCE
     real(8) :: ELECTRONIC_DENSITY_MATRIX_TOLERANCE
     real(8) :: TOTAL_ENERGY_TOLERANCE
     real(8) :: STRONG_ENERGY_TOLERANCE	!< Permite controlar la convergencia exaustiva de todas las especies
     real(8) :: DENSITY_FACTOR_THRESHOLD !< define cuando recalcula un elemeto Gij de acuedo con el valor de Pij
     real(8) :: DIIS_SWITCH_THRESHOLD
     real(8) :: DIIS_SWITCH_THRESHOLD_BKP
     real(8) :: ELECTRONIC_LEVEL_SHIFTING
     real(8) :: NONELECTRONIC_LEVEL_SHIFTING
     real(8) :: WAVE_FUNCTION_SCALE
     integer :: SCF_NONELECTRONIC_MAX_ITERATIONS
     integer :: SCF_ELECTRONIC_MAX_ITERATIONS
     integer :: SCF_MAX_ITERATIONS !< Limita el numero global de itereaciones para cualquier especie
     integer :: SCF_GLOBAL_MAXIMUM_ITERATIONS
     integer :: LISTS_SIZE
     integer :: CONVERGENCE_METHOD
     integer :: DIIS_DIMENSIONALITY
     integer :: ITERATION_SCHEME
     character(10) :: SCF_ELECTRONIC_TYPE_GUESS
     character(10) :: SCF_NONELECTRONIC_TYPE_GUESS
     character(10) ::  SCF_CONVERGENCE_CRITERIUM !< Define si el criterio de convergencia esta dado por cambios de  energia o densidad
     logical :: DIIS_ERROR_IN_DAMPING
     logical :: ACTIVATE_LEVEL_SHIFTING
     logical :: EXCHANGE_ORBITALS_IN_SCF
     logical :: DEBUG_SCFS

     !!***************************************************************************
     !! Hartree-Fock options
     !!
     character(20) :: FROZEN_PARTICLE(5)
     logical :: FREEZE_NON_ELECTRONIC_ORBITALS
     logical :: HARTREE_PRODUCT_GUESS
     logical :: READ_COEFFICIENTS
     logical :: NO_SCF
     logical :: FINITE_MASS_CORRECTION
     logical :: REMOVE_TRANSLATIONAL_CONTAMINATION
     logical :: BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE
     logical :: BUILD_MIXED_DENSITY_MATRIX
     logical :: ONLY_ELECTRONIC_EFFECT
     logical :: ELECTRONIC_WAVEFUNCTION_ANALYSIS
     logical :: IS_OPEN_SHELL

     !!***************************************************************************
     !! Parameter to control geometry optimization
     !!
     real(8) :: NUMERICAL_DERIVATIVE_DELTA
     real(8) :: MINIMIZATION_INITIAL_STEP_SIZE
     real(8) :: MINIMIZATION_LINE_TOLERANCE
     real(8) :: MINIMIZATION_TOLERANCE_GRADIENT
     integer :: MINIMIZATION_MAX_ITERATION
     character(10) :: MINIMIZATION_METHOD
     character(10) :: MINIMIZATION_LIBRARY
     character(50) :: COORDINATES
     character(10) :: ENERGY_CALCULATOR
     logical :: ANALYTIC_GRADIENT
     logical :: MINIMIZATION_WITH_SINGLE_POINT !< Realiza proceso de minimizacion empleando parametros de un calculo de punto sencillo
     logical :: USE_SYMMETRY_IN_MATRICES
     logical :: RESTART_OPTIMIZATION
     logical :: OPTIMIZE_WITH_CP_CORRECTION
     logical :: CP_CORRECTION
     logical :: TDHF
     logical :: OPTIMIZE
     logical :: FIRST_STEP
     logical :: OPTIMIZE_WITH_MP
     logical :: PROJECT_HESSIANE

     !!***************************************************************************
     !! Parameter of atomic conectivity
     !!
     real(8) :: BOND_DISTANCE_FACTOR
     real(8) :: BOND_ANGLE_THRESHOLD
     real(8) :: DIHEDRAL_ANGLE_THRESHOLD

     !!***************************************************************************
     !! Parameter to control MPn theory
     !!
     integer :: MOLLER_PLESSET_CORRECTION
     integer :: MP_FROZEN_CORE_BOUNDARY
     logical :: MP_ONLY_ELECTRONIC_CORRECTION

     !!***************************************************************************
     !! Parameter to control cosmo  
     !!
     logical :: COSMO
	 	 real(8) :: COSMO_SOLVENT_DIALECTRIC
	   integer :: COSMO_MIN_BEM
	   integer :: COSMO_MAX_BEM
	   real(8) :: COSMO_RSOLV

     !!***************************************************************************

     !! Parameter to control the propagator theory module
     !!
     logical :: PT_ONLY_ONE_SPECIE_CORRECTION
     logical :: PT_SELF_ENERGY_SCAN
     logical :: PT_TRANSITION_OPERATOR
     logical :: PT_JUST_ONE_ORBITAL
     real(8) :: PT_SELF_ENERGY_SPACING  !< Control variable for the spacing in Self-Energy scanning
     real(8) :: PT_SELF_ENERGY_RANGE  !< Control variable for range of Self-Energy scannning
     integer :: PT_ORDER !< order of propagator correction
     integer :: PT_MAX_ITERATIONS
     integer :: PT_ITERATION_METHOD_2_LIMIT
     integer :: PT_ITERATION_SCHEME
     integer :: PT_MAX_NUMBER_POLES_SEARCHED


     !!***************************************************************************
     !! Control print level and units
     !!
     integer :: FORMAT_NUMBER_OF_COLUMNS
     integer :: UNIT_FOR_OUTPUT_FILE
     integer :: UNIT_FOR_MP2_INTEGRALS_FILE
     integer :: UNIT_FOR_MOLECULAR_ORBITALS_FILE
     integer :: PRINT_LEVEL
     character(50) :: UNITS
     real(8) :: DOUBLE_ZERO_THRESHOLD

     !!***************************************************************************
     !! CISD - FCI
     !!
     character(20) :: CONFIGURATION_INTERACTION_LEVEL

     !!*****************************************************
     !! Parameter to general control
     !!
     character(50) :: METHOD
     logical :: TRANSFORM_TO_CENTER_OF_MASS
     logical :: ARE_THERE_DUMMY_ATOMS
     logical :: IS_THERE_EXTERNAL_POTENTIAL
     logical :: IS_THERE_INTERPARTICLE_POTENTIAL
     logical :: IS_THERE_OUTPUT
     logical :: IS_THERE_FROZEN_PARTICLE
     integer :: DIMENSIONALITY

     !!*****************************************************
     !! Density Functional Theory Options
     !!
     character(50) :: ELECTRON_CORRELATION_FUNCTIONAL
     character(50) :: ELECTRON_EXCHANGE_FUNCTIONAL
     character(50) :: ELECTRON_NUCLEAR_CORRELATION_FUNCTIONAL
     integer :: POLARIZATION_ORDER
     integer :: NUMBER_OF_BLOCKS_IN_AUXILIARY_FUNCTIONS
     logical :: FUKUI_FUNCTIONS
     logical :: AUXILIARY_DENSITY
     logical :: STORE_THREE_CENTER_ELECTRON_INTEGRALS
     logical :: CALL_DFT

     !!*****************************************************
     !! External Potential Options
     !!
     logical :: NUMERICAL_INTEGRATION_FOR_EXTERNAL_POTENTIAL
     logical :: NUMERICAL_INTEGRATION_FOR_OVERLAP
     real(8) :: MAX_INTERVAL_IN_NUMERICAL_INTEGRATION
     real(8) :: RELATIVE_ERROR_IN_NUMERICAL_INTEGRATION
     integer :: INITIAL_NUMBER_OF_EVALUATIONS
     integer :: INCREASE_NUMBER_OF_EVALUATIONS
     integer :: MINIMUM_NUMBER_OF_EVALUATIONS
     integer :: MAXIMUM_NUMBER_OF_EVALUATIONS
     real(8) :: STEP_IN_NUMERICAL_INTEGRATION
     real(8) :: COEFFICIENT_FOR_GAUSSIAN_EXTERNAL_POTENTIAL
     real(8) :: EXPONENT_FOR_GAUSSIAN_EXTERNAL_POTENTIAL
     real(8) :: ORIGIN_OF_GAUSSIAN_EXTERNAL_POTENTIAL
     character(50) :: NUMERICAL_INTEGRATION_METHOD

     !!*****************************************************
     !! Graphs Options
     !!
     integer :: NUMBER_OF_POINTS_PER_DIMENSION

     !!*****************************************************
     !! Cubes Options
     !!
     integer :: CUBE_POINTS_DENSITY
     real(8) :: VOLUME_DENSITY_THRESHOLD

     !!***************************************************** 
     !! Molecular Mechanics Options                                                        
     character(50) :: FORCE_FIELD
     logical :: ELECTROSTATIC_MM
     logical :: CHARGES_MM
     logical :: PRINT_MM

     !!***************************************************** 
     !! Output Options                                                        
     logical :: MOLDEN_FILE

     !!*****************************************************
     !! Properties Options
     logical :: CALCULATE_INTERPARTICLE_DISTANCES
     logical :: CALCULATE_DENSITY_VOLUME

     !!*****************************************************
     !! Miscelaneous Options
     !!
     real(8) :: MO_FRACTION_OCCUPATION
     integer :: IONIZE_MO
     character(50) :: IONIZE_SPECIE(10)
     character(50) :: EXCITE_SPECIE
     integer :: NUMBER_OF_CORES

     !!***************************************************************************
     !! Environment variables
     !!
     character(255) :: HOME_DIRECTORY="NONE"
     character(255) :: DATA_DIRECTORY="NONE"
     character(255) :: EXTERNAL_COMMAND="NONE"
     character(30) :: EXTERNAL_SOFTWARE_NAME="NONE"
     character(255) :: UFF_PARAMETERS_DATABASE="NONE"
     character(255) :: ATOMIC_ELEMENTS_DATABASE="NONE"
     character(255) :: BASIS_SET_DATABASE="NONE"
     character(255) :: POTENTIALS_DATABASE="NONE"
     character(255) :: ELEMENTAL_PARTICLES_DATABASE="NONE"
     character(100) :: INPUT_FILE=""

  end type CONTROL

  !< Namelist definition
  
  !!***************************************************************************
  !! Parameter to control Integrals library
  !!  
  real(8) :: LowdinParameters_tv
  real(8) :: LowdinParameters_integralThreshold
  integer :: LowdinParameters_integralStackSize
  character(20) :: LowdinParameters_integralDestiny
  character(20) :: LowdinParameters_integralScheme

  !!***************************************************************************
  !! Parameter to control SCF program
  !!
  real(8) :: LowdinParameters_scfNonelectronicEnergyTolerance
  real(8) :: LowdinParameters_scfElectronicEnergyTolerance
  real(8) :: LowdinParameters_nonelectronicDensityMatrixTolerance
  real(8) :: LowdinParameters_electronicDensityMatrixTolerance
  real(8) :: LowdinParameters_totalEnergyTolerance
  real(8) :: LowdinParameters_strongEnergyTolerance
  real(8) :: LowdinParameters_densityFactorThreshold
  real(8) :: LowdinParameters_diisSwitchThreshold
  real(8) :: LowdinParameters_diisSwitchThreshold_bkp
  real(8) :: LowdinParameters_electronicLevelShifting
  real(8) :: LowdinParameters_nonelectronicLevelShifting
  real(8) :: LowdinParameters_waveFunctionScale
  integer :: LowdinParameters_scfNonelectronicMaxIterations
  integer :: LowdinParameters_scfElectronicMaxIterations
  integer :: LowdinParameters_scfMaxIterations
  integer :: LowdinParameters_scfGlobalMaxIterations
  integer :: LowdinParameters_listSize
  integer :: LowdinParameters_convergenceMethod
  integer :: LowdinParameters_diisDimensionality
  integer :: LowdinParameters_iterationScheme
  character(10) :: LowdinParameters_scfElectronicTypeGuess
  character(10) :: LowdinParameters_scfNonelectronicTypeGuess
  character(10) :: LowdinParameters_scfConvergenceCriterium
  logical :: LowdinParameters_diisErrorInDamping
  logical :: LowdinParameters_activateLevelShifting
  logical :: LowdinParameters_exchangeOrbitalsInSCF
  logical :: LowdinParameters_debugScfs

  !!*****************************************************
  !! Hartree-Fock Options
  !!
  character(20) :: LowdinParameters_frozen(5)
  logical :: LowdinParameters_freezeNonElectronicOrbitals
  logical :: LowdinParameters_hartreeProductGuess
  logical :: LowdinParameters_readCoefficients
  logical :: LowdinParameters_noSCF
  logical :: LowdinParameters_finiteMassCorrection
  logical :: LowdinParameters_removeTranslationalContamination
  logical :: LowdinParameters_buildTwoParticlesMatrixForOneParticle
  logical :: LowdinParameters_buildMixedDensityMatrix
  logical :: LowdinParameters_onlyElectronicEffect
  logical :: LowdinParameters_electronicWaveFunctionAnalysis
  logical :: LowdinParameters_isOpenShell

  !!***************************************************************************
  !! Parameter to control geometry optimization
  !!
  real(8) :: LowdinParameters_numericalDerivativeDelta
  real(8) :: LowdinParameters_minimizationInitialStepSize
  real(8) :: LowdinParameters_minimizationLineTolerance
  real(8) :: LowdinParameters_minimizationToleranceGradient
  integer :: LowdinParameters_minimizationMaxIteration
  character(10) :: LowdinParameters_minimizationMethod
  character(10) :: LowdinParameters_minimizationLibrary
  character(50) :: LowdinParameters_coordinates
  character(20) :: LowdinParameters_energyCalculator
  logical :: LowdinParameters_analyticGradient
  logical :: LowdinParameters_minimizationWithSinglePoint
  logical :: LowdinParameters_useSymmetryInMatrices
  logical :: LowdinParameters_restartOptimization
  logical :: LowdinParameters_optimizeWithCpCorrection
  logical :: LowdinParameters_cpCorrection
  logical :: LowdinParameters_TDHF
  logical :: LowdinParameters_optimize
  logical :: LowdinParameters_optimizeGeometryWithMP
  logical :: LowdinParameters_projectHessiane

  !!***************************************************************************
  !! Parameter of atomic conectivity
  !!
  real(8) :: LowdinParameters_bondDistanceFactor
  real(8) :: LowdinParameters_bondAngleThreshold
  real(8) :: LowdinParameters_dihedralAngleThreshold

  !!***************************************************************************
  !! Parameter to control MPn theory
  !!
  integer :: LowdinParameters_mpCorrection
  integer :: LowdinParameters_mpFrozenCoreBoundary
  logical :: LowdinParameters_mpOnlyElectronicCorrection

  !!***************************************************************************
  !! Parameter to control cosmo theory
  !!
  logical :: LowdinParameters_cosmo
  real(8) :: LowdinParameters_cosmo_solvent_dialectric
  integer :: LowdinParameters_cosmo_min_bem
  integer :: LowdinParameters_cosmo_max_bem
  real(8) :: LowdinParameters_cosmo_rsolv
  
  !!***************************************************************************
  !! Parameter to control the propagator theory module
  !!
  logical :: LowdinParameters_ptOnlyOneSpecieCorrection
  logical :: LowdinParameters_selfEnergyScan
  logical :: LowdinParameters_ptTransitionOperator
  logical :: LowdinParameters_ptJustOneOrbital
  real(8) :: LowdinParameters_selfEnergySpacing
  real(8) :: LowdinParameters_selfEnergyRange
  integer :: LowdinParameters_ptOrder
  integer :: LowdinParameters_ptMaxIterations
  integer :: LowdinParameters_ptIterationMethod2Limit
  integer :: LowdinParameters_ptIterationScheme
  integer :: LowdinParameters_ptMaxNumberOfPolesSearched


  !!***************************************************************************
  !! Control print level and units
  !!
  integer :: LowdinParameters_formatNumberOfColumns
  integer :: LowdinParameters_unitForOutputFile
  integer :: LowdinParameters_unitForMolecularOrbitalsFile
  integer :: LowdinParameters_unitForMP2IntegralsFile
  integer :: LowdinParameters_printLevel
  character(50) :: LowdinParameters_units    
  real(8) :: LowdinParameters_doubleZeroThreshold

  !!***************************************************************************
  !! CISD - FCI
  !!
  character(20) :: LowdinParameters_configurationInteractionLevel

  !!*****************************************************
  !! Parameter to general control
  !!
  character(50) :: LowdinParameters_method
  logical :: LowdinParameters_transformToCenterOfMass
  logical :: LowdinParameters_areThereDummyAtoms
  logical :: LowdinParameters_isThereExternalPotential
  logical :: LowdinParameters_isThereInterparticlePotential
  logical :: LowdinParameters_isThereOutput
  logical :: LowdinParameters_isThereFrozenParticle
  integer :: LowdinParameters_dimensionality

  !!*****************************************************
  !! Density Functional Theory Options
  !!
  character(50) :: LowdinParameters_electronCorrelationFunctional
  character(50) :: LowdinParameters_electronExchangeFunctional
  character(50) :: LowdinParameters_electronNuclearCorrelationFunctional
  integer :: LowdinParameters_polarizationOrder
  integer :: LowdinParameters_numberOfBlocksInAuxiliaryFunctions
  logical :: LowdinParameters_fukuiFunctions
  logical :: LowdinParameters_auxiliaryDensity
  logical :: LowdinParameters_storeThreeCenterElectronIntegrals
  logical :: LowdinParameters_callDft

  !!*****************************************************
  !! External Potential Options
  !!
  logical :: LowdinParameters_numericalIntegrationForExternalPotential
  logical :: LowdinParameters_numericalIntegrationForOverlap  
  real(8) :: LowdinParameters_maxIntervalInNumericalIntegration
  real(8) :: LowdinParameters_relativeErrorInNumericalIntegration
  integer :: LowdinParameters_initialNumberOfEvaluations
  integer :: LowdinParameters_increaseNumberOfEvaluations
  integer :: LowdinParameters_minimumNumberOfEvaluations
  integer :: LowdinParameters_maximumNumberOfEvaluations
  real(8) :: LowdinParameters_stepInNumericalIntegration
  real(8) :: LowdinParameters_coefficientForGaussianExternalPotential
  real(8) :: LowdinParameters_exponentForGaussianExternalPotential
  real(8) :: LowdinParameters_originOfGaussianExternalPotential
  character(50) :: LowdinParameters_numericalIntegrationMethod

  !!*****************************************************
  !! Graphs Options
  !!
  integer :: LowdinParameters_numberOfPointsPerDimension

  !!*****************************************************
  !! Cubes Options
  !!
  integer :: LowdinParameters_cubePointsDensity
  real(8) :: LowdinParameters_volumeDensityThreshold
  
  !!***************************************************** 
  !! Molecular Mechanics Options                                                        
  character(50) :: LowdinParameters_forceField
  logical :: LowdinParameters_electrostaticMM
  logical :: LowdinParameters_chargesMM
  logical :: LowdinParameters_printMM

  !!*****************************************************                                                           
  !! Output Options                                                                                                 
  logical :: LowdinParameters_moldenFile

  !!*****************************************************
  !! Properties Options
  logical :: LowdinParameters_calculateInterparticleDistances
  logical :: LowdinParameters_calculateDensityVolume

  !!*****************************************************
  !! Miscelaneous Options
  !!
  real(8) :: LowdinParameters_MOFractionOccupation
  integer :: LowdinParameters_ionizeMO
  character(50) :: LowdinParameters_ionizeSpecie(10)
  character(50) :: LowdinParameters_exciteSpecie
  integer :: LowdinParameters_numberOfCores

  !!***************************************************************************
  !! Environment variables
  !!
  character(255) :: LowdinParameters_homeDirectory
  character(255) :: LowdinParameters_dataDirectory
  character(255) :: LowdinParameters_externalCommand
  character(30) :: LowdinParameters_externalSoftwareName
  character(255) :: LowdinParameters_uffParametersDataBase
  character(255) :: LowdinParameters_atomicElementsDataBase
  character(255) :: LowdinParameters_basisSetDataBase
  character(255) :: LowdinParameters_potentialsDataBase
  character(255) :: LowdinParameters_elementalParticlesDataBase
  character(100) :: LowdinParameters_inputFile

  
  NAMELIST /LowdinParameters/ &
       !!***************************************************************************
       !! Parameter to control Integrals library
       !!  
       LowdinParameters_tv,&
       LowdinParameters_integralThreshold,&
       LowdinParameters_integralStackSize,&
       LowdinParameters_integralDestiny,&
       LowdinParameters_integralScheme,&

       !!***************************************************************************
       !! Parameter to control SCF program
       !!
       LowdinParameters_scfNonelectronicEnergyTolerance,&
       LowdinParameters_scfElectronicEnergyTolerance,&
       LowdinParameters_nonelectronicDensityMatrixTolerance,&
       LowdinParameters_electronicDensityMatrixTolerance,&
       LowdinParameters_totalEnergyTolerance,&
       LowdinParameters_strongEnergyTolerance,&
       LowdinParameters_densityFactorThreshold,&
       LowdinParameters_diisSwitchThreshold,&
       LowdinParameters_diisSwitchThreshold_bkp,&
       LowdinParameters_electronicLevelShifting,&
       LowdinParameters_nonelectronicLevelShifting,&
       LowdinParameters_waveFunctionScale,&
       LowdinParameters_scfNonelectronicMaxIterations,&
       LowdinParameters_scfElectronicMaxIterations,&
       LowdinParameters_scfMaxIterations,&
       LowdinParameters_scfGlobalMaxIterations,&
       LowdinParameters_listSize,&
       LowdinParameters_convergenceMethod,&
       LowdinParameters_diisDimensionality,&
       LowdinParameters_iterationScheme,&
       LowdinParameters_scfElectronicTypeGuess,&
       LowdinParameters_scfNonelectronicTypeGuess,&
       LowdinParameters_scfConvergenceCriterium,&
       LowdinParameters_diisErrorInDamping,&
       LowdinParameters_activateLevelShifting,&
       LowdinParameters_exchangeOrbitalsInSCF,&
       LowdinParameters_debugScfs,&
       
       !!*****************************************************
       !! Hartree-Fock Options
       !!
       LowdinParameters_frozen,&
       LowdinParameters_freezeNonElectronicOrbitals,&
       LowdinParameters_hartreeProductGuess,&
       LowdinParameters_readCoefficients,&
       LowdinParameters_noSCF,&
       LowdinParameters_finiteMassCorrection,&
       LowdinParameters_removeTranslationalContamination,&
       LowdinParameters_buildTwoParticlesMatrixForOneParticle,&
       LowdinParameters_buildMixedDensityMatrix,&
       LowdinParameters_onlyElectronicEffect,&
       LowdinParameters_electronicWaveFunctionAnalysis,&
       LowdinParameters_isOpenShell, &

       !!***************************************************************************
       !! Parameter to control geometry optimization
       !!
       LowdinParameters_numericalDerivativeDelta,&
       LowdinParameters_minimizationInitialStepSize,&
       LowdinParameters_minimizationLineTolerance,&
       LowdinParameters_minimizationToleranceGradient,&
       LowdinParameters_minimizationMaxIteration,&
       LowdinParameters_minimizationMethod,&
       LowdinParameters_minimizationLibrary,&
       LowdinParameters_coordinates,&
       LowdinParameters_energyCalculator,&
       LowdinParameters_analyticGradient,&
       LowdinParameters_minimizationWithSinglePoint,&
       LowdinParameters_useSymmetryInMatrices,&
       LowdinParameters_restartOptimization,&
       LowdinParameters_optimizeWithCpCorrection,&
       LowdinParameters_cpCorrection,&
       LowdinParameters_TDHF,&
       LowdinParameters_optimize,&
       LowdinParameters_optimizeGeometryWithMP,&
       LowdinParameters_projectHessiane,&
       
       !!***************************************************************************
       !! Parameter of atomic conectivity
       !!
       LowdinParameters_bondDistanceFactor,&
       LowdinParameters_bondAngleThreshold,&
       LowdinParameters_dihedralAngleThreshold,&
       
       !!***************************************************************************
       !! Parameter to control MPn theory
       !!
       LowdinParameters_mpCorrection,&
       LowdinParameters_mpFrozenCoreBoundary,&
       LowdinParameters_mpOnlyElectronicCorrection,&
       
       !!***************************************************************************
       !! Parameter to control cosmo theory
       !!
       LowdinParameters_cosmo,& 
       LowdinParameters_cosmo_solvent_dialectric,& 
       LowdinParameters_cosmo_min_bem,&  
       LowdinParameters_cosmo_max_bem,&  
       LowdinParameters_cosmo_rsolv,& 
       
       !!***************************************************************************
       !! Parameter to control the propagator theory module
       !!
       LowdinParameters_ptOnlyOneSpecieCorrection,&
       LowdinParameters_selfEnergyScan,&
       LowdinParameters_ptTransitionOperator,&
       LowdinParameters_ptJustOneOrbital,&
       LowdinParameters_selfEnergySpacing,&
       LowdinParameters_selfEnergyRange,&
       LowdinParameters_ptOrder,&
       LowdinParameters_ptMaxIterations,&
       LowdinParameters_ptIterationMethod2Limit,&
       LowdinParameters_ptIterationScheme,&
       LowdinParameters_ptMaxNumberOfPolesSearched,&
       
       !!***************************************************************************
       !! Control print level and units
       !!
       LowdinParameters_formatNumberOfColumns,&
       LowdinParameters_unitForOutputFile,&
       LowdinParameters_unitForMolecularOrbitalsFile,&
       LowdinParameters_unitForMP2IntegralsFile,&
       LowdinParameters_printLevel,&
       LowdinParameters_units    ,&
       LowdinParameters_doubleZeroThreshold,&
       
       !!***************************************************************************
       !! CISD - FCI
       !!
       LowdinParameters_configurationInteractionLevel,&
       
       !!*****************************************************
       !! Parameter to general control
       !!
       LowdinParameters_method,&
       LowdinParameters_transformToCenterOfMass,&
       LowdinParameters_areThereDummyAtoms,&
       LowdinParameters_isThereExternalPotential,&
       LowdinParameters_isThereInterparticlePotential,&
       LowdinParameters_isThereOutput,&
       LowdinParameters_isThereFrozenParticle,&
       LowdinParameters_dimensionality,&
       
       !!*****************************************************
       !! Density Functional Theory Options
       !!
       LowdinParameters_electronCorrelationFunctional,&
       LowdinParameters_electronExchangeFunctional,&
       LowdinParameters_electronNuclearCorrelationFunctional,&
       LowdinParameters_polarizationOrder,&
       LowdinParameters_numberOfBlocksInAuxiliaryFunctions,&
       LowdinParameters_fukuiFunctions,&
       LowdinParameters_auxiliaryDensity,&
       LowdinParameters_storeThreeCenterElectronIntegrals,&
       LowdinParameters_callDft,&
       
       !!*****************************************************
       !! External Potential Options
       !!
       LowdinParameters_numericalIntegrationForExternalPotential,&
       LowdinParameters_numericalIntegrationForOverlap  ,&
       LowdinParameters_maxIntervalInNumericalIntegration,&
       LowdinParameters_relativeErrorInNumericalIntegration,&
       LowdinParameters_initialNumberOfEvaluations,&
       LowdinParameters_increaseNumberOfEvaluations,&
       LowdinParameters_minimumNumberOfEvaluations,&
       LowdinParameters_maximumNumberOfEvaluations,&
       LowdinParameters_stepInNumericalIntegration,&
       LowdinParameters_coefficientForGaussianExternalPotential,&
       LowdinParameters_exponentForGaussianExternalPotential,&
       LowdinParameters_originOfGaussianExternalPotential,&
       LowdinParameters_numericalIntegrationMethod,&
       
       !!*****************************************************
       !! Graphs Options
       !!
       LowdinParameters_numberOfPointsPerDimension,&
       
       !!*****************************************************
       !! Cubes Options
       !!
       LowdinParameters_cubePointsDensity,&
       LowdinParameters_volumeDensityThreshold,&

       !!***************************************************** 
       !! Molecular Mechanics Options                                                        
       LowdinParameters_forceField,&
       LowdinParameters_electrostaticMM,&
       LowdinParameters_chargesMM,&
       LowdinParameters_printMM,&

       !!*****************************************************                                                      
       !! Output Options                                                                                            
       LowdinParameters_moldenFile,&
       
       !!*****************************************************
       !! Properties Options
       LowdinParameters_calculateInterparticleDistances,&
       LowdinParameters_calculateDensityVolume,&
       
       !!*****************************************************
       !! Miscelaneous Options
       !!
       LowdinParameters_MOFractionOccupation,&
       LowdinParameters_ionizeMO,&
       LowdinParameters_ionizeSpecie,&
       LowdinParameters_exciteSpecie,&
       LowdinParameters_numberOfCores,&

       !!***************************************************************************
       !! Variables de ambiente al sistema de archivos del programa
       !!
       LowdinParameters_inputFile, &
       LowdinParameters_homeDirectory,&
       LowdinParameters_dataDirectory,&
       LowdinParameters_externalCommand,&
       LowdinParameters_externalSoftwareName,&
       LowdinParameters_uffParametersDataBase,&
       LowdinParameters_atomicElementsDataBase,&
       LowdinParameters_basisSetDataBase,&
       LowdinParameters_potentialsDataBase,&
       LowdinParameters_elementalParticlesDataBase


  public :: &
       CONTROL_start, &
       CONTROL_load, &
       CONTROL_save, &
       CONTROL_copy, &
       CONTROL_show

  private :: &
       CONTROL_getHomeDirectory, &
       CONTROL_getDataDirectory, &
       CONTROL_getExternalCommand, &
       CONTROL_getExternalSoftwareName, &
       CONTROL_exception
  
  !> Singleton
  type(CONTROL), save, public :: CONTROL_instance
  
contains

  
  !>
  !! @brief Set defaults
  subroutine CONTROL_start()
    implicit none

    !! Set defaults for namelist
    
    !!***************************************************************************
    !! Parameter to control Integrals library
    !!  
    LowdinParameters_tv = 1.0E-6
    LowdinParameters_integralThreshold = 1.0E-10
    LowdinParameters_integralStackSize = 30000
    LowdinParameters_integralDestiny = "MEMORY" !! "MEMORY" or "DISK"
    LowdinParameters_integralScheme = "LIBINT" !! LIBINT or RYS

    !!***************************************************************************
    !! Parameter to control SCF program
    !!
    LowdinParameters_scfNonelectronicEnergyTolerance = 1.0E-5
    LowdinParameters_scfElectronicEnergyTolerance =  1.0E-6
    LowdinParameters_nonelectronicDensityMatrixTolerance =  5.0E-4
    LowdinParameters_electronicDensityMatrixTolerance = 1.0E-6
    LowdinParameters_totalEnergyTolerance = 1.0E-7
    LowdinParameters_strongEnergyTolerance = 1.0E-9
    LowdinParameters_densityFactorThreshold = 1.0E-8
    LowdinParameters_diisSwitchThreshold = 0.5
    LowdinParameters_diisSwitchThreshold_bkp = 0.5 
    LowdinParameters_electronicLevelShifting = 0.0
    LowdinParameters_nonelectronicLevelShifting = 0.0
    LowdinParameters_waveFunctionScale = 1000.0
    LowdinParameters_scfNonelectronicMaxIterations = 10000
    LowdinParameters_scfElectronicMaxIterations = 20000
    LowdinParameters_scfMaxIterations = 20000
    LowdinParameters_scfGlobalMaxIterations = 5000 
    LowdinParameters_listSize = -20
    LowdinParameters_convergenceMethod = 1 !!(0) NONE, (1) DAMPING, (2) DIIS, (3) LEVEL SHIFTING (4) DAMPING/DIIS
    LowdinParameters_diisDimensionality = 10
    LowdinParameters_iterationScheme = 3 !!(0) NONELECRONIC FULLY / e- (1) ELECTRONIC FULLY (2) CONVERGED INDIVIDIALLY (3) SCHEMESIMULTANEOUS
    LowdinParameters_scfElectronicTypeGuess = "HCORE"
    LowdinParameters_scfNonelectronicTypeGuess = "HCORE"
    LowdinParameters_scfConvergenceCriterium = "ENERGY"
    LowdinParameters_diisErrorInDamping = .false.
    LowdinParameters_activateLevelShifting = .false.
    LowdinParameters_exchangeOrbitalsInSCF = .false.
    LowdinParameters_debugScfs = .false.

    !!*****************************************************
    !! Hartree-Fock Options
    !!
    LowdinParameters_frozen = "NONE"
    LowdinParameters_freezeNonElectronicOrbitals = .false.
    LowdinParameters_hartreeProductGuess = .false.
    LowdinParameters_readCoefficients = .false.
    LowdinParameters_noSCF = .false.
    LowdinParameters_finiteMassCorrection = .false.
    LowdinParameters_removeTranslationalContamination = .false.
    LowdinParameters_buildTwoParticlesMatrixForOneParticle = .false.
    LowdinParameters_buildMixedDensityMatrix = .false.
    LowdinParameters_onlyElectronicEffect = .false.
    LowdinParameters_electronicWaveFunctionAnalysis = .false.
    LowdinParameters_isOpenShell = .false.

    !!***************************************************************************
    !! Parameter to control geometry optimization
    !!
    LowdinParameters_numericalDerivativeDelta = 1.0E-3
    LowdinParameters_minimizationInitialStepSize = 0.1_8
    LowdinParameters_minimizationLineTolerance = 0.1_8
    LowdinParameters_minimizationToleranceGradient = 1.0E-5
    LowdinParameters_minimizationMaxIteration = 100
    LowdinParameters_minimizationMethod = "TR"
    LowdinParameters_minimizationLibrary = "GENERIC"
    LowdinParameters_coordinates = "CARTESIAN"
    LowdinParameters_energyCalculator = "INTERNAL"
    LowdinParameters_analyticGradient = .true.
    LowdinParameters_minimizationWithSinglePoint = .true.
    LowdinParameters_useSymmetryInMatrices = .false.
    LowdinParameters_restartOptimization = .false.
    LowdinParameters_optimizeWithCpCorrection = .false.
    LowdinParameters_cpCorrection = .false.
    LowdinParameters_TDHF = .false.
    LowdinParameters_optimize = .false.
    LowdinParameters_optimizeGeometryWithMP = .false.
    LowdinParameters_projectHessiane = .true.

    !!***************************************************************************
    !! Parameter of atomic conectivity
    !!
    LowdinParameters_bondDistanceFactor = 1.3_8
    LowdinParameters_bondAngleThreshold = 170.0_8
    LowdinParameters_dihedralAngleThreshold = 170.0_8

    !!***************************************************************************
    !! Parameter to control MPn theory
    !!
    LowdinParameters_mpCorrection = 1
    LowdinParameters_mpFrozenCoreBoundary = 1
    LowdinParameters_mpOnlyElectronicCorrection = .false.

    !!***************************************************************************
    !! Parameter to control cosmo theory
    !!
    LowdinParameters_cosmo = .false.
    LowdinParameters_cosmo_solvent_dialectric = 78.4d+00
    LowdinParameters_cosmo_min_bem = 2
    LowdinParameters_cosmo_max_bem = 3
    LowdinParameters_cosmo_rsolv =0.5d+00
    
    !!***************************************************************************
    !! Parameter to control the propagator theory module
    !!
    LowdinParameters_ptOnlyOneSpecieCorrection = .false.
    LowdinParameters_selfEnergyScan = .false.
    LowdinParameters_ptTransitionOperator = .false.
    LowdinParameters_ptJustOneOrbital = .false.
    LowdinParameters_selfEnergySpacing = 0.5_8
    LowdinParameters_selfEnergyRange = 5.0_8
    LowdinParameters_ptOrder = 1
    LowdinParameters_ptMaxIterations = 50
    LowdinParameters_ptIterationMethod2Limit = 1
    LowdinParameters_ptIterationScheme = 1
    LowdinParameters_ptMaxNumberOfPolesSearched = 10



    !!***************************************************************************
    !! Control print level and units
    !!
    LowdinParameters_formatNumberOfColumns = 5
    LowdinParameters_unitForOutputFile = 6
    LowdinParameters_unitForMolecularOrbitalsFile = 8 
    LowdinParameters_unitForMP2IntegralsFile = 7
    LowdinParameters_printLevel =  1 !! (1) normal output, (5) method (6) metod and WF (7) method, WF and GLOBAL(8) method, WF, GLOBAL, SCF
    LowdinParameters_units = "ANGS"
    LowdinParameters_doubleZeroThreshold = 1.0E-12

    !!***************************************************************************
    !! CISD - FCI
    !!
    LowdinParameters_configurationInteractionLevel = "NONE"

    !!*****************************************************
    !! Parameter to general control
    !!
    LowdinParameters_method = "NONE"
    LowdinParameters_transformToCenterOfMass = .false.
    LowdinParameters_areThereDummyAtoms = .false.
    LowdinParameters_isThereExternalPotential = .false.
    LowdinParameters_isThereInterparticlePotential = .false.
    LowdinParameters_isThereOutput = .false.
    LowdinParameters_isThereFrozenParticle = .false. 
    LowdinParameters_dimensionality = 3

    !!*****************************************************
    !! Density Functional Theory Options
    !!
    LowdinParameters_electronCorrelationFunctional = "NONE"
    LowdinParameters_electronExchangeFunctional = "NONE"
    LowdinParameters_electronNuclearCorrelationFunctional = "NONE"
    LowdinParameters_polarizationOrder = 1
    LowdinParameters_numberOfBlocksInAuxiliaryFunctions = 3
    LowdinParameters_fukuiFunctions = .false.
    LowdinParameters_auxiliaryDensity = .false.
    LowdinParameters_storeThreeCenterElectronIntegrals = .true.
    LowdinParameters_callDft = .false.

    !!*****************************************************
    !! External Potential Options
    !!
    LowdinParameters_numericalIntegrationForExternalPotential = .false. 
    LowdinParameters_numericalIntegrationForOverlap   = .false.
    LowdinParameters_maxIntervalInNumericalIntegration = 10.0_8
    LowdinParameters_relativeErrorInNumericalIntegration = 1E-10
    LowdinParameters_initialNumberOfEvaluations = 5000
    LowdinParameters_increaseNumberOfEvaluations = 5000
    LowdinParameters_minimumNumberOfEvaluations = 10000
    LowdinParameters_maximumNumberOfEvaluations = 20000
    LowdinParameters_stepInNumericalIntegration = 0.1
    LowdinParameters_coefficientForGaussianExternalPotential = 0.0_8
    LowdinParameters_exponentForGaussianExternalPotential = 0.0_8
    LowdinParameters_originOfGaussianExternalPotential = 0.0_8
    LowdinParameters_numericalIntegrationMethod = "NONE"

    !!*****************************************************
    !! Graphs Options
    !!
    LowdinParameters_numberOfPointsPerDimension = 50

    !!*****************************************************
    !! Cubes Options
    !!
    LowdinParameters_cubePointsDensity = 125
    LowdinParameters_volumeDensityThreshold = 1E-3

    !!***************************************************** 
    !! Molecular Mechanics Options                                                        
    LowdinParameters_forceField = "UFF"
    LowdinParameters_electrostaticMM = .false.
    LowdinParameters_chargesMM = .false.
    LowdinParameters_printMM = .false.

    !!*****************************************************                                                       
    !! Output Options          
    LowdinParameters_moldenFile = .false.

    !!*****************************************************
    !! Properties Options
    LowdinParameters_calculateInterparticleDistances = .false.
    LowdinParameters_calculateDensityVolume = .false.

    !!*****************************************************
    !! Miscelaneous Options
    !!
    LowdinParameters_MOFractionOccupation = 1.0_8
    LowdinParameters_ionizeMO = 0
    LowdinParameters_ionizeSpecie = "NONE"
    LowdinParameters_exciteSpecie = "NONE"
    LowdinParameters_numberOfCores = 1
    LowdinParameters_exciteSpecie = "NONE"     
    !$OMP PARALLEL
    LowdinParameters_numberOfCores = OMP_get_thread_num() + 1
    !$OMP END PARALLEL 

    !!***************************************************************************
    !! Variables de ambiente al sistema de archivos del programa
    !!
    LowdinParameters_homeDirectory = CONTROL_getHomeDirectory()
    LowdinParameters_dataDirectory = CONTROL_getDataDirectory()
    LowdinParameters_externalCommand = CONTROL_getExternalCommand()
    LowdinParameters_externalSoftwareName = CONTROL_getExternalSoftwareName()
    LowdinParameters_uffParametersDataBase = "/dataBases/uffParameters.lib"
    LowdinParameters_atomicElementsDataBase = "/dataBases/atomicElements.lib"
    LowdinParameters_basisSetDataBase = "/basis/"
    LowdinParameters_potentialsDataBase = "/potentials/"
    LowdinParameters_elementalParticlesDataBase = "/dataBases/elementalParticles.lib"
    LowdinParameters_inputFile = CONTROL_instance%INPUT_FILE

    !! Set defaults for CONTROL Object
    
    !!***************************************************************************
    !!***************************************************************************
    !!***************************************************************************
    !!***************************************************************************

    !!***************************************************************************    
    !! Parameter to control Integrals library                       
    !!
    CONTROL_instance%TV = 1.0E-6
    CONTROL_instance%INTEGRAL_THRESHOLD = 1.0E-10
    CONTROL_instance%INTEGRAL_STACK_SIZE = 30000
    CONTROL_instance%INTEGRAL_DESTINY = "MEMORY" !! "MEMORY" or "DISK"
    CONTROL_instance%INTEGRAL_SCHEME = "LIBINT" !! LIBINT or Rys

    !!***************************************************************************
    !! Parameter to control SCF program
    !!
    CONTROL_instance%SCF_NONELECTRONIC_ENERGY_TOLERANCE = 1.0E-5
    CONTROL_instance%SCF_ELECTRONIC_ENERGY_TOLERANCE =  1.0E-6
    CONTROL_instance%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE =  5.0E-4
    CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE = 1.0E-6
    CONTROL_instance%TOTAL_ENERGY_TOLERANCE = 1.0E-7
    CONTROL_instance%STRONG_ENERGY_TOLERANCE = 1.0E-9
    CONTROL_instance%DENSITY_FACTOR_THRESHOLD = 1.0E-8
    CONTROL_instance%DIIS_SWITCH_THRESHOLD = 0.5
    CONTROL_instance%DIIS_SWITCH_THRESHOLD_BKP = 0.5 
    CONTROL_instance%ELECTRONIC_LEVEL_SHIFTING = 0.0
    CONTROL_instance%NONELECTRONIC_LEVEL_SHIFTING = 0.0
    CONTROL_instance%WAVE_FUNCTION_SCALE = 1000.0
    CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS = 10000
    CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS = 20000
    CONTROL_instance%SCF_MAX_ITERATIONS = 20000
    CONTROL_instance%SCF_GLOBAL_MAXIMUM_ITERATIONS = 5000 
    CONTROL_instance%LISTS_SIZE = -20
    CONTROL_instance%CONVERGENCE_METHOD = 1 !!(0) NONE, (1) DAMPING, (2) DIIS, (3) LEVEL SHIFTING (4) DAMPING/DIIS
    CONTROL_instance%DIIS_DIMENSIONALITY = 10
    CONTROL_instance%ITERATION_SCHEME = 3 !!(0) NONELECRONIC FULLY / e- (1) ELECTRONIC FULLY (2) CONVERGED INDIVIDIALLY (3) SCHEMESIMULTANEOUS
    CONTROL_instance%SCF_ELECTRONIC_TYPE_GUESS = "HCORE"
    CONTROL_instance%SCF_NONELECTRONIC_TYPE_GUESS = "HCORE"
    CONTROL_instance%SCF_CONVERGENCE_CRITERIUM = "ENERGY"
    CONTROL_instance%DIIS_ERROR_IN_DAMPING = .false.
    CONTROL_instance%ACTIVATE_LEVEL_SHIFTING = .false.
    CONTROL_instance%EXCHANGE_ORBITALS_IN_SCF = .false.
    CONTROL_instance%DEBUG_SCFS = .false.

    !!***************************************************************************                                              
    !! Hartree-Fock options                                                                                                    
    !!                                                                                                                         
    CONTROL_instance%FROZEN_PARTICLE = "NONE"
    CONTROL_instance%FREEZE_NON_ELECTRONIC_ORBITALS = .false.
    CONTROL_instance%HARTREE_PRODUCT_GUESS = .false.
    CONTROL_instance%READ_COEFFICIENTS = .false.
    CONTROL_instance%NO_SCF = .false.
    CONTROL_instance%FINITE_MASS_CORRECTION = .false.
    CONTROL_instance%REMOVE_TRANSLATIONAL_CONTAMINATION = .false.
    CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE = .false.
    CONTROL_instance%BUILD_MIXED_DENSITY_MATRIX = .false.
    CONTROL_instance%ONLY_ELECTRONIC_EFFECT = .false.
    CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS = .false.
    CONTROL_instance%IS_OPEN_SHELL = .false.

    !!***************************************************************************                                              
    !! Parameter to control geometry optimization                                                                              
    !!                                                                                                                         
    CONTROL_instance%NUMERICAL_DERIVATIVE_DELTA = 1.0E-3
    CONTROL_instance%MINIMIZATION_INITIAL_STEP_SIZE = 0.1_8
    CONTROL_instance%MINIMIZATION_LINE_TOLERANCE = 0.1_8
    CONTROL_instance%MINIMIZATION_TOLERANCE_GRADIENT = 1.0E-5
    CONTROL_instance%MINIMIZATION_MAX_ITERATION = 100
    CONTROL_instance%MINIMIZATION_METHOD = "TR"
    CONTROL_instance%MINIMIZATION_LIBRARY = "GENERIC"
    CONTROL_instance%COORDINATES = "CARTESIAN"
    CONTROL_instance%ENERGY_CALCULATOR = "INTERNAL"
    CONTROL_instance%ANALYTIC_GRADIENT = .true.
    CONTROL_instance%MINIMIZATION_WITH_SINGLE_POINT = .true.
    CONTROL_instance%USE_SYMMETRY_IN_MATRICES = .false.
    CONTROL_instance%RESTART_OPTIMIZATION = .false.
    CONTROL_instance%OPTIMIZE_WITH_CP_CORRECTION = .false.
    CONTROL_instance%CP_CORRECTION = .false.
    CONTROL_instance%TDHF = .false.
    CONTROL_instance%OPTIMIZE = .false.
    CONTROL_instance%FIRST_STEP = .true.
    CONTROL_instance%OPTIMIZE_WITH_MP = .false.
    CONTROL_instance%PROJECT_HESSIANE = .true.

    !!***************************************************************************                                              
    !! Parameter of atomic conectivity                                                                                         
    !!                                                                                                                         
    CONTROL_instance%BOND_DISTANCE_FACTOR = 1.3_8
    CONTROL_instance%BOND_ANGLE_THRESHOLD = 180.0_8
    CONTROL_instance%DIHEDRAL_ANGLE_THRESHOLD = 180.0_8

    !!***************************************************************************                                              
    !! Parameter to control MPn theory                                                                                         
    !!                                                                                                                         
    CONTROL_instance%MOLLER_PLESSET_CORRECTION = 1
    CONTROL_instance%MP_FROZEN_CORE_BOUNDARY = 1
    CONTROL_instance%MP_ONLY_ELECTRONIC_CORRECTION = .false.

    !!***************************************************************************                                              
    !! Parameter to control cosmo method                                                                                         
    !!                                                                                                                         
    CONTROL_instance%COSMO = .false.
    CONTROL_instance%COSMO_SOLVENT_DIALECTRIC= 78.4d+00 
    CONTROL_instance%COSMO_MIN_BEM= 2
    CONTROL_instance%COSMO_MAX_BEM= 3
    CONTROL_instance%COSMO_RSOLV= 0.5d+00

    !!***************************************************************************                                              
    !! Parameter to control the propagator theory module                                                                       
    !!                                                                                                                         
    CONTROL_instance%PT_ONLY_ONE_SPECIE_CORRECTION = .false.
    CONTROL_instance%PT_SELF_ENERGY_SCAN = .false.
    CONTROL_instance%PT_TRANSITION_OPERATOR = .false.
    CONTROL_instance%PT_JUST_ONE_ORBITAL = .false.
    CONTROL_instance%PT_SELF_ENERGY_SPACING = 0.5_8
    CONTROL_instance%PT_SELF_ENERGY_RANGE = 5.0_8
    CONTROL_instance%PT_ORDER = 1
    CONTROL_instance%PT_MAX_ITERATIONS = 50
    CONTROL_instance%PT_ITERATION_METHOD_2_LIMIT = 1
    CONTROL_instance%PT_ITERATION_SCHEME = 1
    CONTROL_instance%PT_MAX_NUMBER_POLES_SEARCHED = 10


    !!***************************************************************************                                              
    !! Control print level and units                                                                                           
    !!                                                                                                                         
    CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS = 5
    CONTROL_instance%UNIT_FOR_OUTPUT_FILE = 6
    CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE = 8 
    CONTROL_instance%UNIT_FOR_MOLECULAR_ORBITALS_FILE = 7
    CONTROL_instance%PRINT_LEVEL =  1 !! (1) normal output, (5) method (6) metod and WF (7) method, WF and GLOBAL(8) method, WF, GLOBAL, SCF
    CONTROL_instance%UNITS = "ANGS"
    CONTROL_instance%DOUBLE_ZERO_THRESHOLD = 1.0E-12

    !!***************************************************************************                                              
    !! CISD - FCI                                                                                                              
    !!                                                                                                                         
    CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL = "NONE"

    !!*****************************************************                                                                    
    !! Parameter to general control                                                                                            
    !!                                                                                                                         
    CONTROL_instance%METHOD = "NONE"
    CONTROL_instance%TRANSFORM_TO_CENTER_OF_MASS = .false.
    CONTROL_instance%ARE_THERE_DUMMY_ATOMS = .false.
    CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL = .false.
    CONTROL_instance%IS_THERE_INTERPARTICLE_POTENTIAL = .false.
    CONTROL_instance%IS_THERE_OUTPUT = .false.
    CONTROL_instance%IS_THERE_FROZEN_PARTICLE = .false. 
    CONTROL_instance%DIMENSIONALITY = 3

    !!*****************************************************                                                                    
    !! Density Functional Theory Options                                                                                       
    !!                                                                                                                         
    CONTROL_instance%ELECTRON_CORRELATION_FUNCTIONAL = "NONE"
    CONTROL_instance%ELECTRON_EXCHANGE_FUNCTIONAL = "NONE"
    CONTROL_instance%ELECTRON_NUCLEAR_CORRELATION_FUNCTIONAL = "NONE"
    CONTROL_instance%POLARIZATION_ORDER = 1
    CONTROL_instance%NUMBER_OF_BLOCKS_IN_AUXILIARY_FUNCTIONS = 3
    CONTROL_instance%FUKUI_FUNCTIONS = .false.
    CONTROL_instance%AUXILIARY_DENSITY = .false.
    CONTROL_instance%STORE_THREE_CENTER_ELECTRON_INTEGRALS = .true.
    CONTROL_instance%CALL_DFT = .false.

    !!*****************************************************                                                                    
    !! External Potential Options                                                                                              
    !!                                                                                                                         
    CONTROL_instance%NUMERICAL_INTEGRATION_FOR_EXTERNAL_POTENTIAL = .false. 
    CONTROL_instance%NUMERICAL_INTEGRATION_FOR_OVERLAP = .false.
    CONTROL_instance%MAX_INTERVAL_IN_NUMERICAL_INTEGRATION = 10.0_8
    CONTROL_instance%RELATIVE_ERROR_IN_NUMERICAL_INTEGRATION = 1E-10
    CONTROL_instance%INITIAL_NUMBER_OF_EVALUATIONS = 5000
    CONTROL_instance%INCREASE_NUMBER_OF_EVALUATIONS = 5000
    CONTROL_instance%MINIMUM_NUMBER_OF_EVALUATIONS = 10000
    CONTROL_instance%MAXIMUM_NUMBER_OF_EVALUATIONS = 20000
    CONTROL_instance%STEP_IN_NUMERICAL_INTEGRATION = 0.1
    CONTROL_instance%COEFFICIENT_FOR_GAUSSIAN_EXTERNAL_POTENTIAL = 0.0_8
    CONTROL_instance%EXPONENT_FOR_GAUSSIAN_EXTERNAL_POTENTIAL = 0.0_8
    CONTROL_instance%ORIGIN_OF_GAUSSIAN_EXTERNAL_POTENTIAL = 0.0_8
    CONTROL_instance%NUMERICAL_INTEGRATION_METHOD = "NONE"

    !!*****************************************************                                                                    
    !! Graphs Options                                                                                                          
    !!                                                                                                                         
    CONTROL_instance%NUMBER_OF_POINTS_PER_DIMENSION = 50

    !!*****************************************************                                                                    
    !! Cubes Options                                                                                                           
    !!                                                                                                                         
    CONTROL_instance%CUBE_POINTS_DENSITY = 125
    CONTROL_instance%VOLUME_DENSITY_THRESHOLD = 1E-3

    !!***************************************************** 
    !! Molecular Mechanics Options                                                        
    CONTROL_instance%FORCE_FIELD = "UFF"
    CONTROL_instance%ELECTROSTATIC_MM = .false.
    CONTROL_instance%CHARGES_MM = .false. 
    CONTROL_instance%PRINT_MM = .false.

    !!*****************************************************  
    !! Output Options     
    CONTROL_instance%MOLDEN_FILE = .false.

    !!*****************************************************                                                                    
    !! Properties Options                                                                                                      
    CONTROL_instance%CALCULATE_INTERPARTICLE_DISTANCES = .false.
    CONTROL_instance%CALCULATE_DENSITY_VOLUME = .false.

    !!*****************************************************                                                                    
    !! Miscelaneous Options                                                                                                    
    !!                                                                                                                         
    CONTROL_instance%MO_FRACTION_OCCUPATION = 1.0_8
    CONTROL_instance%IONIZE_MO = 0
    CONTROL_instance%IONIZE_SPECIE = "NONE"
    CONTROL_instance%EXCITE_SPECIE = "NONE"                                                            
    CONTROL_instance%NUMBER_OF_CORES = 1
    !$OMP PARALLEL
    CONTROL_instance%NUMBER_OF_CORES = OMP_get_thread_num() + 1
    !$OMP END PARALLEL 

    !!***************************************************************************                                              
    !! Environment variables                                                                                                   
    !!                                                                                                                         
    CONTROL_instance%HOME_DIRECTORY = CONTROL_getHomeDirectory()
    CONTROL_instance%DATA_DIRECTORY = CONTROL_getDataDirectory()
    CONTROL_instance%EXTERNAL_COMMAND = CONTROL_getExternalCommand()
    CONTROL_instance%EXTERNAL_SOFTWARE_NAME = CONTROL_getExternalSoftwareName()
    CONTROL_instance%UFF_PARAMETERS_DATABASE = "/dataBases/uffParameters.lib"
    CONTROL_instance%ATOMIC_ELEMENTS_DATABASE = "/dataBases/atomicElements.lib"
    CONTROL_instance%BASIS_SET_DATABASE = "/basis/"
    CONTROL_instance%POTENTIALS_DATABASE = "/potentials/"
    CONTROL_instance%ELEMENTAL_PARTICLES_DATABASE = "/dataBases/elementalParticles.lib"
    CONTROL_instance%INPUT_FILE = CONTROL_instance%INPUT_FILE
    
  end subroutine CONTROL_start

  !<
  !! @brief load the LowdinParameter name list and set all values for parameters object
  subroutine CONTROL_load(unit)
    implicit none

    integer, optional :: unit
    
    integer :: uunit
    integer:: stat
    
    uunit = 4
    if(present(unit)) uunit = unit
    
    !! Reload file
    rewind(uunit)

    !! Reads name-list
    read(uunit,NML=LowdinParameters, iostat=stat)

    !! Check the process
    if(stat > 0 ) then

       call CONTROL_exception( ERROR, "Class object CONTROL in the load function", &
            "check the CONTROL block in your input file")
    end if
    
    !! Fixing load...
    if ( LowdinParameters_diisSwitchThreshold > LowdinParameters_doubleZeroThreshold  ) then
       LowdinParameters_diisSwitchThreshold_bkp = LowdinParameters_diisSwitchThreshold
    end if

    if ( LowdinParameters_frozen(1) /= "NONE" ) then
       LowdinParameters_isThereFrozenParticle = .true.
    end if

    if ( LowdinParameters_convergenceMethod >= 2) then

       LowdinParameters_electronicDensityMatrixTolerance = LowdinParameters_electronicDensityMatrixTolerance / 2.5_8
       LowdinParameters_scfNonelectronicEnergyTolerance = LowdinParameters_scfNonelectronicEnergyTolerance / 2.5_8
       LowdinParameters_scfElectronicEnergyTolerance = LowdinParameters_scfElectronicEnergyTolerance / 2.5_8
       LowdinParameters_nonelectronicDensityMatrixTolerance = LowdinParameters_nonelectronicDensityMatrixTolerance / 2.5_8
       LowdinParameters_totalEnergyTolerance = LowdinParameters_totalEnergyTolerance / 2.5_8

    end if
    
    !!***************************************************************************      
    !! Parameter to control Integrals library                                          
    !!                                                                                 
    CONTROL_instance%TV = LowdinParameters_tv
    CONTROL_instance%INTEGRAL_THRESHOLD = LowdinParameters_integralThreshold
    CONTROL_instance%INTEGRAL_STACK_SIZE = LowdinParameters_integralStackSize
    CONTROL_instance%INTEGRAL_DESTINY = LowdinParameters_integralDestiny
    CONTROL_instance%INTEGRAL_SCHEME =  LowdinParameters_integralScheme
    
    !!***************************************************************************      
    !! Parameter to control SCF program                                                
    !!                                                                                 
    CONTROL_instance%SCF_NONELECTRONIC_ENERGY_TOLERANCE = LowdinParameters_scfNonelectronicEnergyTolerance
    CONTROL_instance%SCF_ELECTRONIC_ENERGY_TOLERANCE = LowdinParameters_scfElectronicEnergyTolerance
    CONTROL_instance%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE = LowdinParameters_nonelectronicDensityMatrixTolerance
    CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE = LowdinParameters_electronicDensityMatrixTolerance
    CONTROL_instance%TOTAL_ENERGY_TOLERANCE = LowdinParameters_totalEnergyTolerance
    CONTROL_instance%STRONG_ENERGY_TOLERANCE = LowdinParameters_strongEnergyTolerance
    CONTROL_instance%DENSITY_FACTOR_THRESHOLD = LowdinParameters_densityFactorThreshold
    CONTROL_instance%DIIS_SWITCH_THRESHOLD = LowdinParameters_diisSwitchThreshold
    CONTROL_instance%DIIS_SWITCH_THRESHOLD_BKP = LowdinParameters_diisSwitchThreshold_bkp
    CONTROL_instance%ELECTRONIC_LEVEL_SHIFTING = LowdinParameters_electronicLevelShifting
    CONTROL_instance%NONELECTRONIC_LEVEL_SHIFTING = LowdinParameters_nonelectronicLevelShifting
    CONTROL_instance%WAVE_FUNCTION_SCALE = LowdinParameters_waveFunctionScale
    CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS = LowdinParameters_scfNonelectronicMaxIterations
    CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS = LowdinParameters_scfElectronicMaxIterations
    CONTROL_instance%SCF_MAX_ITERATIONS = LowdinParameters_scfMaxIterations
    CONTROL_instance%SCF_GLOBAL_MAXIMUM_ITERATIONS = LowdinParameters_scfGlobalMaxIterations
    CONTROL_instance%LISTS_SIZE = LowdinParameters_listSize
    CONTROL_instance%CONVERGENCE_METHOD = LowdinParameters_convergenceMethod
    CONTROL_instance%DIIS_DIMENSIONALITY = LowdinParameters_diisDimensionality
    CONTROL_instance%ITERATION_SCHEME = LowdinParameters_iterationScheme
    CONTROL_instance%SCF_ELECTRONIC_TYPE_GUESS = LowdinParameters_scfElectronicTypeGuess
    CONTROL_instance%SCF_NONELECTRONIC_TYPE_GUESS = LowdinParameters_scfNonelectronicTypeGuess
    CONTROL_instance%SCF_CONVERGENCE_CRITERIUM = LowdinParameters_scfConvergenceCriterium
    CONTROL_instance%DIIS_ERROR_IN_DAMPING = LowdinParameters_diisErrorInDamping
    CONTROL_instance%ACTIVATE_LEVEL_SHIFTING = LowdinParameters_activateLevelShifting
    CONTROL_instance%EXCHANGE_ORBITALS_IN_SCF = LowdinParameters_exchangeOrbitalsInSCF
    CONTROL_instance%DEBUG_SCFS = LowdinParameters_debugScfs
                                                                                                                                                                                          
    !!*****************************************************                            
    !! Hartree-Fock Options                                                            
    !!                                                                                 
    CONTROL_instance%FROZEN_PARTICLE = LowdinParameters_frozen
    CONTROL_instance%FREEZE_NON_ELECTRONIC_ORBITALS = LowdinParameters_freezeNonElectronicOrbitals
    CONTROL_instance%HARTREE_PRODUCT_GUESS = LowdinParameters_hartreeProductGuess
    CONTROL_instance%READ_COEFFICIENTS = LowdinParameters_readCoefficients
    CONTROL_instance%NO_SCF = LowdinParameters_noSCF
    CONTROL_instance%FINITE_MASS_CORRECTION = LowdinParameters_finiteMassCorrection
    CONTROL_instance%REMOVE_TRANSLATIONAL_CONTAMINATION = LowdinParameters_removeTranslationalContamination
    CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE = LowdinParameters_buildTwoParticlesMatrixForOneParticle
    CONTROL_instance%BUILD_MIXED_DENSITY_MATRIX = LowdinParameters_buildMixedDensityMatrix
    CONTROL_instance%ONLY_ELECTRONIC_EFFECT = LowdinParameters_onlyElectronicEffect
    CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS = LowdinParameters_electronicWaveFunctionAnalysis
    CONTROL_instance%IS_OPEN_SHELL = LowdinParameters_isOpenShell
                                                                                                                                                                                          
    !!***************************************************************************      
    !! Parameter to control geometry optimization                                      
    !!                                                                                 
    CONTROL_instance%NUMERICAL_DERIVATIVE_DELTA = LowdinParameters_numericalDerivativeDelta
    CONTROL_instance%MINIMIZATION_INITIAL_STEP_SIZE = LowdinParameters_minimizationInitialStepSize
    CONTROL_instance%MINIMIZATION_LINE_TOLERANCE = LowdinParameters_minimizationLineTolerance
    CONTROL_instance%MINIMIZATION_TOLERANCE_GRADIENT = LowdinParameters_minimizationToleranceGradient
    CONTROL_instance%MINIMIZATION_MAX_ITERATION = LowdinParameters_minimizationMaxIteration
    CONTROL_instance%MINIMIZATION_METHOD = LowdinParameters_minimizationMethod
    CONTROL_instance%MINIMIZATION_LIBRARY = LowdinParameters_minimizationLibrary
    CONTROL_instance%COORDINATES = LowdinParameters_coordinates
    CONTROL_instance%ENERGY_CALCULATOR = LowdinParameters_energyCalculator
    CONTROL_instance%ANALYTIC_GRADIENT = LowdinParameters_analyticGradient
    CONTROL_instance%MINIMIZATION_WITH_SINGLE_POINT = LowdinParameters_minimizationWithSinglePoint
    CONTROL_instance%USE_SYMMETRY_IN_MATRICES = LowdinParameters_useSymmetryInMatrices
    CONTROL_instance%RESTART_OPTIMIZATION = LowdinParameters_restartOptimization
    CONTROL_instance%OPTIMIZE_WITH_CP_CORRECTION = LowdinParameters_optimizeWithCpCorrection
    CONTROL_instance%CP_CORRECTION = LowdinParameters_cpCorrection
    CONTROL_instance%TDHF = LowdinParameters_TDHF
    CONTROL_instance%OPTIMIZE = LowdinParameters_optimize
    CONTROL_instance%OPTIMIZE_WITH_MP = LowdinParameters_optimizeGeometryWithMP
    CONTROL_instance%PROJECT_HESSIANE = LowdinParameters_projectHessiane
                                                                                                                                                                                          
    !!***************************************************************************      
    !! Parameter of atomic conectivity                                                 
    !!                                                                                 
    CONTROL_instance%BOND_DISTANCE_FACTOR = LowdinParameters_bondDistanceFactor
    CONTROL_instance%BOND_ANGLE_THRESHOLD = LowdinParameters_bondAngleThreshold
    CONTROL_instance%DIHEDRAL_ANGLE_THRESHOLD = LowdinParameters_dihedralAngleThreshold
                                                                                                                                                                                          
    !!***************************************************************************      
    !! Parameter to control MPn theory                                                 
    !!                                                                                 
    CONTROL_instance%MOLLER_PLESSET_CORRECTION = LowdinParameters_mpCorrection
    CONTROL_instance%MP_FROZEN_CORE_BOUNDARY = LowdinParameters_mpFrozenCoreBoundary
    CONTROL_instance%MP_ONLY_ELECTRONIC_CORRECTION = LowdinParameters_mpOnlyElectronicCorrection
                                                                                                                                                                                          
    !!***************************************************************************      
    !! Parameter to control cosmo method                                               
    !!                                                                                 
    CONTROL_instance%COSMO = LowdinParameters_cosmo
  	CONTROL_instance%COSMO_SOLVENT_DIALECTRIC= LowdinParameters_cosmo_solvent_dialectric
  	CONTROL_instance%COSMO_MIN_BEM=LowdinParameters_cosmo_min_bem
  	CONTROL_instance%COSMO_MAX_BEM=LowdinParameters_cosmo_max_bem
  	CONTROL_instance%COSMO_RSOLV=LowdinParameters_cosmo_rsolv

    !!***************************************************************************      
    !! Parameter to control the propagator theory module                               
    !!                                                                                 
    CONTROL_instance%PT_ONLY_ONE_SPECIE_CORRECTION = LowdinParameters_ptOnlyOneSpecieCorrection
    CONTROL_instance%PT_SELF_ENERGY_SCAN = LowdinParameters_selfEnergyScan
    CONTROL_instance%PT_TRANSITION_OPERATOR = LowdinParameters_ptTransitionOperator
    CONTROL_instance%PT_JUST_ONE_ORBITAL = LowdinParameters_ptJustOneOrbital
    CONTROL_instance%PT_SELF_ENERGY_SPACING = LowdinParameters_selfEnergySpacing
    CONTROL_instance%PT_SELF_ENERGY_RANGE = LowdinParameters_selfEnergyRange
    CONTROL_instance%PT_ORDER = LowdinParameters_ptOrder
    CONTROL_instance%PT_MAX_ITERATIONS = LowdinParameters_ptMaxIterations
    CONTROL_instance%PT_ITERATION_METHOD_2_LIMIT = LowdinParameters_ptIterationMethod2Limit
    CONTROL_instance%PT_ITERATION_SCHEME = LowdinParameters_ptIterationScheme
    CONTROL_instance%PT_MAX_NUMBER_POLES_SEARCHED = LowdinParameters_ptMaxNumberOfPolesSearched
                                                                                                                                                                                          
    !!***************************************************************************      
    !! Control print level and units                                                   
    !!                                                                                 
    CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS = LowdinParameters_formatNumberOfColumns
    CONTROL_instance%UNIT_FOR_OUTPUT_FILE = LowdinParameters_unitForOutputFile
    CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE = LowdinParameters_unitForMolecularOrbitalsFile
    CONTROL_instance%UNIT_FOR_MOLECULAR_ORBITALS_FILE = LowdinParameters_unitForMP2IntegralsFile
    CONTROL_instance%PRINT_LEVEL = LowdinParameters_printLevel
    CONTROL_instance%UNITS = LowdinParameters_units
    CONTROL_instance%DOUBLE_ZERO_THRESHOLD = LowdinParameters_doubleZeroThreshold
                                                                                                                                                                                          
    !!***************************************************************************      
    !! CISD - FCI                                                                      
    !!                                                                                 
    CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL = LowdinParameters_configurationInteractionLevel
                                                                                                                                                                                          
    !!*****************************************************                            
    !! Parameter to general control                                                    
    !!                                                                                 
    CONTROL_instance%METHOD = LowdinParameters_method
    CONTROL_instance%TRANSFORM_TO_CENTER_OF_MASS = LowdinParameters_transformToCenterOfMass
    CONTROL_instance%ARE_THERE_DUMMY_ATOMS = LowdinParameters_areThereDummyAtoms
    CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL = LowdinParameters_isThereExternalPotential
    CONTROL_instance%IS_THERE_INTERPARTICLE_POTENTIAL = LowdinParameters_isThereInterparticlePotential
    CONTROL_instance%IS_THERE_OUTPUT = LowdinParameters_isThereOutput
    CONTROL_instance%IS_THERE_FROZEN_PARTICLE = LowdinParameters_isThereFrozenParticle
    CONTROL_instance%DIMENSIONALITY = LowdinParameters_dimensionality
                                                                                                                                                                                          
    !!*****************************************************                            
    !! Density Functional Theory Options                                               
    !!                                                                                 
    CONTROL_instance%ELECTRON_CORRELATION_FUNCTIONAL = LowdinParameters_electronCorrelationFunctional
    CONTROL_instance%ELECTRON_EXCHANGE_FUNCTIONAL = LowdinParameters_electronExchangeFunctional
    CONTROL_instance%ELECTRON_NUCLEAR_CORRELATION_FUNCTIONAL = LowdinParameters_electronNuclearCorrelationFunctional
    CONTROL_instance%POLARIZATION_ORDER = LowdinParameters_polarizationOrder
    CONTROL_instance%NUMBER_OF_BLOCKS_IN_AUXILIARY_FUNCTIONS = LowdinParameters_numberOfBlocksInAuxiliaryFunctions
    CONTROL_instance%FUKUI_FUNCTIONS = LowdinParameters_fukuiFunctions
    CONTROL_instance%AUXILIARY_DENSITY = LowdinParameters_auxiliaryDensity
    CONTROL_instance%STORE_THREE_CENTER_ELECTRON_INTEGRALS = LowdinParameters_storeThreeCenterElectronIntegrals
    CONTROL_instance%CALL_DFT = LowdinParameters_callDft
                                                                                                                                                                                          
    !!*****************************************************                            
    !! External Potential Options                                                      
    !!                                                                                 
    CONTROL_instance%NUMERICAL_INTEGRATION_FOR_EXTERNAL_POTENTIAL = LowdinParameters_numericalIntegrationForExternalPotential
    CONTROL_instance%NUMERICAL_INTEGRATION_FOR_OVERLAP = LowdinParameters_numericalIntegrationForOverlap  
    CONTROL_instance%MAX_INTERVAL_IN_NUMERICAL_INTEGRATION = LowdinParameters_maxIntervalInNumericalIntegration
    CONTROL_instance%RELATIVE_ERROR_IN_NUMERICAL_INTEGRATION = LowdinParameters_relativeErrorInNumericalIntegration
    CONTROL_instance%INITIAL_NUMBER_OF_EVALUATIONS = LowdinParameters_initialNumberOfEvaluations
    CONTROL_instance%INCREASE_NUMBER_OF_EVALUATIONS = LowdinParameters_increaseNumberOfEvaluations
    CONTROL_instance%MINIMUM_NUMBER_OF_EVALUATIONS = LowdinParameters_minimumNumberOfEvaluations
    CONTROL_instance%MAXIMUM_NUMBER_OF_EVALUATIONS = LowdinParameters_maximumNumberOfEvaluations
    CONTROL_instance%STEP_IN_NUMERICAL_INTEGRATION = LowdinParameters_stepInNumericalIntegration
    CONTROL_instance%COEFFICIENT_FOR_GAUSSIAN_EXTERNAL_POTENTIAL = LowdinParameters_coefficientForGaussianExternalPotential
    CONTROL_instance%EXPONENT_FOR_GAUSSIAN_EXTERNAL_POTENTIAL = LowdinParameters_exponentForGaussianExternalPotential
    CONTROL_instance%ORIGIN_OF_GAUSSIAN_EXTERNAL_POTENTIAL = LowdinParameters_originOfGaussianExternalPotential
    CONTROL_instance%NUMERICAL_INTEGRATION_METHOD = LowdinParameters_numericalIntegrationMethod
                                                                                                                                                                                          
    !!*****************************************************                            
    !! Graphs Options                                                                  
    !!                                                                                 
    CONTROL_instance%NUMBER_OF_POINTS_PER_DIMENSION = LowdinParameters_numberOfPointsPerDimension
                                                                                                                                                                                          
    !!*****************************************************                            
    !! Cubes Options                                                                   
    !!                                                                                 
    CONTROL_instance%CUBE_POINTS_DENSITY = LowdinParameters_cubePointsDensity
    CONTROL_instance%VOLUME_DENSITY_THRESHOLD = LowdinParameters_volumeDensityThreshold
                                                                                                                                
    !!***************************************************** 
    !! Molecular Mechanics Options                                                        
    CONTROL_instance%FORCE_FIELD = LowdinParameters_forceField
    CONTROL_instance%ELECTROSTATIC_MM = LowdinParameters_electrostaticMM
    CONTROL_instance%CHARGES_MM = LowdinParameters_chargesMM
    CONTROL_instance%PRINT_MM = LowdinParameters_printMM
    !!*****************************************************   
    !! Output Options                                                               
    CONTROL_instance%MOLDEN_FILE = LowdinParameters_moldenFile
                                                          
    !!*****************************************************                            
    !! Properties Options                                                              
    CONTROL_instance%CALCULATE_INTERPARTICLE_DISTANCES = LowdinParameters_calculateInterparticleDistances
    CONTROL_instance%CALCULATE_DENSITY_VOLUME = LowdinParameters_calculateDensityVolume
                                                                                                                                                                                          
    !!*****************************************************                            
    !! Miscelaneous Options                                                            
    !!                                                                                 
    CONTROL_instance%MO_FRACTION_OCCUPATION = LowdinParameters_MOFractionOccupation
    CONTROL_instance%IONIZE_MO = LowdinParameters_ionizeMO
    CONTROL_instance%IONIZE_SPECIE = LowdinParameters_ionizeSpecie
    CONTROL_instance%EXCITE_SPECIE = LowdinParameters_exciteSpecie
    CONTROL_instance%NUMBER_OF_CORES = LowdinParameters_numberOfCores

    !!***************************************************************************      
    !! Variables de ambiente al sistema de archivos del programa                       
    !!                                                                                 
    CONTROL_instance%HOME_DIRECTORY = LowdinParameters_homeDirectory
    CONTROL_instance%DATA_DIRECTORY = LowdinParameters_dataDirectory
    CONTROL_instance%EXTERNAL_COMMAND = LowdinParameters_externalCommand
    CONTROL_instance%EXTERNAL_SOFTWARE_NAME = LowdinParameters_externalSoftwareName
    CONTROL_instance%UFF_PARAMETERS_DATABASE = LowdinParameters_uffParametersDataBase
    CONTROL_instance%ATOMIC_ELEMENTS_DATABASE = LowdinParameters_atomicElementsDataBase
    CONTROL_instance%BASIS_SET_DATABASE = LowdinParameters_basisSetDataBase
    CONTROL_instance%POTENTIALS_DATABASE = LowdinParameters_potentialsDataBase
    CONTROL_instance%ELEMENTAL_PARTICLES_DATABASE = LowdinParameters_elementalParticlesDataBase
    CONTROL_instance%INPUT_FILE = LowdinParameters_inputFile
                                                                                              

  end subroutine CONTROL_load

  !> 
  !! @brief Save all options in file
  subroutine CONTROL_save( unit )
    implicit none
    
    integer :: unit
    
    !! Saving de control parameters on the name list.
    
    !!***************************************************************************      
    !! Parameter to control Integrals library                                          
    !!                                                                                 
    LowdinParameters_tv = CONTROL_instance%TV
    LowdinParameters_integralThreshold = CONTROL_instance%INTEGRAL_THRESHOLD
    LowdinParameters_integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE
    LowdinParameters_integralDestiny = CONTROL_instance%INTEGRAL_DESTINY
    LowdinParameters_integralScheme = CONTROL_instance%INTEGRAL_SCHEME
    
    !!***************************************************************************      
    !! Parameter to control SCF program                                                
    !!                                                                                 
    LowdinParameters_scfNonelectronicEnergyTolerance = CONTROL_instance%SCF_NONELECTRONIC_ENERGY_TOLERANCE
    LowdinParameters_scfElectronicEnergyTolerance = CONTROL_instance%SCF_ELECTRONIC_ENERGY_TOLERANCE
    LowdinParameters_nonelectronicDensityMatrixTolerance = CONTROL_instance%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE
    LowdinParameters_electronicDensityMatrixTolerance = CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE
    LowdinParameters_totalEnergyTolerance = CONTROL_instance%TOTAL_ENERGY_TOLERANCE
    LowdinParameters_strongEnergyTolerance = CONTROL_instance%STRONG_ENERGY_TOLERANCE
    LowdinParameters_densityFactorThreshold = CONTROL_instance%DENSITY_FACTOR_THRESHOLD
    LowdinParameters_diisSwitchThreshold = CONTROL_instance%DIIS_SWITCH_THRESHOLD
    LowdinParameters_diisSwitchThreshold_bkp = CONTROL_instance%DIIS_SWITCH_THRESHOLD_BKP
    LowdinParameters_electronicLevelShifting = CONTROL_instance%ELECTRONIC_LEVEL_SHIFTING
    LowdinParameters_nonelectronicLevelShifting = CONTROL_instance%NONELECTRONIC_LEVEL_SHIFTING
    LowdinParameters_waveFunctionScale = CONTROL_instance%WAVE_FUNCTION_SCALE
    LowdinParameters_scfNonelectronicMaxIterations = CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS
    LowdinParameters_scfElectronicMaxIterations = CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS
    LowdinParameters_scfMaxIterations = CONTROL_instance%SCF_MAX_ITERATIONS
    LowdinParameters_scfGlobalMaxIterations = CONTROL_instance%SCF_GLOBAL_MAXIMUM_ITERATIONS
    LowdinParameters_listSize = CONTROL_instance%LISTS_SIZE
    LowdinParameters_convergenceMethod = CONTROL_instance%CONVERGENCE_METHOD
    LowdinParameters_diisDimensionality = CONTROL_instance%DIIS_DIMENSIONALITY
    LowdinParameters_iterationScheme = CONTROL_instance%ITERATION_SCHEME
    LowdinParameters_scfElectronicTypeGuess = CONTROL_instance%SCF_ELECTRONIC_TYPE_GUESS
    LowdinParameters_scfNonelectronicTypeGuess = CONTROL_instance%SCF_NONELECTRONIC_TYPE_GUESS
    LowdinParameters_scfConvergenceCriterium = CONTROL_instance%SCF_CONVERGENCE_CRITERIUM
    LowdinParameters_diisErrorInDamping = CONTROL_instance%DIIS_ERROR_IN_DAMPING
    LowdinParameters_activateLevelShifting = CONTROL_instance%ACTIVATE_LEVEL_SHIFTING
    LowdinParameters_exchangeOrbitalsInSCF = CONTROL_instance%EXCHANGE_ORBITALS_IN_SCF
    LowdinParameters_debugScfs = CONTROL_instance%DEBUG_SCFS
                                                                                                                                                                                          
    !!*****************************************************                            
    !! Hartree-Fock Options                                                            
    !!                                                                                 
    LowdinParameters_frozen = CONTROL_instance%FROZEN_PARTICLE
    LowdinParameters_freezeNonElectronicOrbitals = CONTROL_instance%FREEZE_NON_ELECTRONIC_ORBITALS
    LowdinParameters_hartreeProductGuess = CONTROL_instance%HARTREE_PRODUCT_GUESS
    LowdinParameters_readCoefficients = CONTROL_instance%READ_COEFFICIENTS
    LowdinParameters_noSCF = CONTROL_instance%NO_SCF
    LowdinParameters_finiteMassCorrection = CONTROL_instance%FINITE_MASS_CORRECTION
    LowdinParameters_removeTranslationalContamination = CONTROL_instance%REMOVE_TRANSLATIONAL_CONTAMINATION
    LowdinParameters_buildTwoParticlesMatrixForOneParticle = CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE
    LowdinParameters_buildMixedDensityMatrix = CONTROL_instance%BUILD_MIXED_DENSITY_MATRIX
    LowdinParameters_onlyElectronicEffect = CONTROL_instance%ONLY_ELECTRONIC_EFFECT
    LowdinParameters_electronicWaveFunctionAnalysis = CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS
    LowdinParameters_isOpenShell = CONTROL_instance%IS_OPEN_SHELL
                                                                                                                                                                                          
    !!***************************************************************************      
    !! Parameter to control geometry optimization                                      
    !!                                                                                 
    LowdinParameters_numericalDerivativeDelta = CONTROL_instance%NUMERICAL_DERIVATIVE_DELTA
    LowdinParameters_minimizationInitialStepSize = CONTROL_instance%MINIMIZATION_INITIAL_STEP_SIZE
    LowdinParameters_minimizationLineTolerance = CONTROL_instance%MINIMIZATION_LINE_TOLERANCE
    LowdinParameters_minimizationToleranceGradient = CONTROL_instance%MINIMIZATION_TOLERANCE_GRADIENT
    LowdinParameters_minimizationMaxIteration = CONTROL_instance%MINIMIZATION_MAX_ITERATION
    LowdinParameters_minimizationMethod = CONTROL_instance%MINIMIZATION_METHOD
    LowdinParameters_minimizationLibrary = CONTROL_instance%MINIMIZATION_LIBRARY
    LowdinParameters_coordinates = CONTROL_instance%COORDINATES
    LowdinParameters_energyCalculator = CONTROL_instance%ENERGY_CALCULATOR
    LowdinParameters_analyticGradient = CONTROL_instance%ANALYTIC_GRADIENT
    LowdinParameters_minimizationWithSinglePoint = CONTROL_instance%MINIMIZATION_WITH_SINGLE_POINT
    LowdinParameters_useSymmetryInMatrices = CONTROL_instance%USE_SYMMETRY_IN_MATRICES
    LowdinParameters_restartOptimization = CONTROL_instance%RESTART_OPTIMIZATION
    LowdinParameters_optimizeWithCpCorrection = CONTROL_instance%OPTIMIZE_WITH_CP_CORRECTION
    LowdinParameters_cpCorrection = CONTROL_instance%CP_CORRECTION
    LowdinParameters_TDHF = CONTROL_instance%TDHF
    LowdinParameters_optimize = CONTROL_instance%OPTIMIZE
    LowdinParameters_optimizeGeometryWithMP = CONTROL_instance%OPTIMIZE_WITH_MP
    LowdinParameters_projectHessiane = CONTROL_instance%PROJECT_HESSIANE
                                                                                                                                                                                          
    !!***************************************************************************      
    !! Parameter of atomic conectivity                                                 
    !!                                                                                 
    LowdinParameters_bondDistanceFactor = CONTROL_instance%BOND_DISTANCE_FACTOR
    LowdinParameters_bondAngleThreshold = CONTROL_instance%BOND_ANGLE_THRESHOLD
    LowdinParameters_dihedralAngleThreshold = CONTROL_instance%DIHEDRAL_ANGLE_THRESHOLD
                                                                                                                                                                                          
    !!***************************************************************************      
    !! Parameter to control MPn theory                                                 
    !!                                                                                 
    LowdinParameters_mpCorrection = CONTROL_instance%MOLLER_PLESSET_CORRECTION
    LowdinParameters_mpFrozenCoreBoundary = CONTROL_instance%MP_FROZEN_CORE_BOUNDARY
    LowdinParameters_mpOnlyElectronicCorrection = CONTROL_instance%MP_ONLY_ELECTRONIC_CORRECTION
                                                                                                                                                                                          
    !!***************************************************************************      
    !! Parameter to control cosmo method                                                 
    !!                                                                                 
    LowdinParameters_cosmo = CONTROL_instance%COSMO
		LowdinParameters_cosmo_solvent_dialectric = CONTROL_instance%COSMO_SOLVENT_DIALECTRIC 
		LowdinParameters_cosmo_min_bem = CONTROL_instance%COSMO_MIN_BEM
  	LowdinParameters_cosmo_max_bem = CONTROL_instance%COSMO_MAX_BEM
 		LowdinParameters_cosmo_rsolv = CONTROL_instance%COSMO_RSOLV

                                                                                                                                                                                          
    !!***************************************************************************      
    !! Parameter to control the propagator theory module                               
    !!                                                                                 
    LowdinParameters_ptOnlyOneSpecieCorrection = CONTROL_instance%PT_ONLY_ONE_SPECIE_CORRECTION
    LowdinParameters_selfEnergyScan = CONTROL_instance%PT_SELF_ENERGY_SCAN
    LowdinParameters_ptTransitionOperator = CONTROL_instance%PT_TRANSITION_OPERATOR
    LowdinParameters_ptJustOneOrbital = CONTROL_instance%PT_JUST_ONE_ORBITAL
    LowdinParameters_selfEnergySpacing = CONTROL_instance%PT_SELF_ENERGY_SPACING
    LowdinParameters_selfEnergyRange = CONTROL_instance%PT_SELF_ENERGY_RANGE
    LowdinParameters_ptOrder = CONTROL_instance%PT_ORDER
    LowdinParameters_ptMaxIterations = CONTROL_instance%PT_MAX_ITERATIONS
    LowdinParameters_ptIterationMethod2Limit = CONTROL_instance%PT_ITERATION_METHOD_2_LIMIT
    LowdinParameters_ptIterationScheme = CONTROL_instance%PT_ITERATION_SCHEME
    LowdinParameters_ptMaxNumberOfPolesSearched = CONTROL_instance%PT_MAX_NUMBER_POLES_SEARCHED
                                                                                                                                                                                          
                                                                                  
                                                                                       
    !!***************************************************************************      
    !! Control print level and units                                                   
    !!                                                                                 
    LowdinParameters_formatNumberOfColumns = CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS
    LowdinParameters_unitForOutputFile = CONTROL_instance%UNIT_FOR_OUTPUT_FILE
    LowdinParameters_unitForMolecularOrbitalsFile = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE
    LowdinParameters_unitForMP2IntegralsFile = CONTROL_instance%UNIT_FOR_MOLECULAR_ORBITALS_FILE
    LowdinParameters_printLevel = CONTROL_instance%PRINT_LEVEL
    LowdinParameters_units = CONTROL_instance%UNITS
    LowdinParameters_doubleZeroThreshold = CONTROL_instance%DOUBLE_ZERO_THRESHOLD
                                                                                                                                                                                          
    !!***************************************************************************      
    !! CISD - FCI                                                                      
    !!                                                                                 
    LowdinParameters_configurationInteractionLevel = CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL
                                                                                                                                                                                          
    !!*****************************************************                            
    !! Parameter to general control                                                    
    !!                                                                                 
    LowdinParameters_method = CONTROL_instance%METHOD
    LowdinParameters_transformToCenterOfMass = CONTROL_instance%TRANSFORM_TO_CENTER_OF_MASS
    LowdinParameters_areThereDummyAtoms = CONTROL_instance%ARE_THERE_DUMMY_ATOMS
    LowdinParameters_isThereExternalPotential = CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL
    LowdinParameters_isThereInterparticlePotential = CONTROL_instance%IS_THERE_INTERPARTICLE_POTENTIAL
    LowdinParameters_isThereOutput = CONTROL_instance%IS_THERE_OUTPUT
    LowdinParameters_isThereFrozenParticle = CONTROL_instance%IS_THERE_FROZEN_PARTICLE
    LowdinParameters_dimensionality = CONTROL_instance%DIMENSIONALITY
                                                                                                                                                                                          
    !!*****************************************************                            
    !! Density Functional Theory Options                                               
    !!                                                                                 
    LowdinParameters_electronCorrelationFunctional = CONTROL_instance%ELECTRON_CORRELATION_FUNCTIONAL
    LowdinParameters_electronExchangeFunctional = CONTROL_instance%ELECTRON_EXCHANGE_FUNCTIONAL
    LowdinParameters_electronNuclearCorrelationFunctional = CONTROL_instance%ELECTRON_NUCLEAR_CORRELATION_FUNCTIONAL
    LowdinParameters_polarizationOrder = CONTROL_instance%POLARIZATION_ORDER
    LowdinParameters_numberOfBlocksInAuxiliaryFunctions = CONTROL_instance%NUMBER_OF_BLOCKS_IN_AUXILIARY_FUNCTIONS
    LowdinParameters_fukuiFunctions = CONTROL_instance%FUKUI_FUNCTIONS
    LowdinParameters_auxiliaryDensity = CONTROL_instance%AUXILIARY_DENSITY
    LowdinParameters_storeThreeCenterElectronIntegrals = CONTROL_instance%STORE_THREE_CENTER_ELECTRON_INTEGRALS
    LowdinParameters_callDft = CONTROL_instance%CALL_DFT
                                                                                                                                                                                          
    !!*****************************************************                            
    !! External Potential Options                                                      
    !!                                                                                 
    LowdinParameters_numericalIntegrationForExternalPotential = CONTROL_instance%NUMERICAL_INTEGRATION_FOR_EXTERNAL_POTENTIAL
    LowdinParameters_numericalIntegrationForOverlap   = CONTROL_instance%NUMERICAL_INTEGRATION_FOR_OVERLAP
    LowdinParameters_maxIntervalInNumericalIntegration = CONTROL_instance%MAX_INTERVAL_IN_NUMERICAL_INTEGRATION
    LowdinParameters_relativeErrorInNumericalIntegration = CONTROL_instance%RELATIVE_ERROR_IN_NUMERICAL_INTEGRATION
    LowdinParameters_initialNumberOfEvaluations = CONTROL_instance%INITIAL_NUMBER_OF_EVALUATIONS
    LowdinParameters_increaseNumberOfEvaluations = CONTROL_instance%INCREASE_NUMBER_OF_EVALUATIONS
    LowdinParameters_minimumNumberOfEvaluations = CONTROL_instance%MINIMUM_NUMBER_OF_EVALUATIONS
    LowdinParameters_maximumNumberOfEvaluations = CONTROL_instance%MAXIMUM_NUMBER_OF_EVALUATIONS
    LowdinParameters_stepInNumericalIntegration = CONTROL_instance%STEP_IN_NUMERICAL_INTEGRATION
    LowdinParameters_coefficientForGaussianExternalPotential = CONTROL_instance%COEFFICIENT_FOR_GAUSSIAN_EXTERNAL_POTENTIAL
    LowdinParameters_exponentForGaussianExternalPotential = CONTROL_instance%EXPONENT_FOR_GAUSSIAN_EXTERNAL_POTENTIAL
    LowdinParameters_originOfGaussianExternalPotential = CONTROL_instance%ORIGIN_OF_GAUSSIAN_EXTERNAL_POTENTIAL
    LowdinParameters_numericalIntegrationMethod = CONTROL_instance%NUMERICAL_INTEGRATION_METHOD
                                                                                                                                                                                          
    !!*****************************************************                            
    !! Graphs Options                                                                  
    !!                                                                                 
    LowdinParameters_numberOfPointsPerDimension = CONTROL_instance%NUMBER_OF_POINTS_PER_DIMENSION
                                                                                                                                                                                          
    !!*****************************************************                            
    !! Cubes Options                                                                   
    !!                                                                                 
    LowdinParameters_cubePointsDensity = CONTROL_instance%CUBE_POINTS_DENSITY
    LowdinParameters_volumeDensityThreshold = CONTROL_instance%VOLUME_DENSITY_THRESHOLD

    !!***************************************************** 
    !! Molecular Mechanics Options                                                        
     LowdinParameters_forceField = CONTROL_instance%FORCE_FIELD
     LowdinParameters_electrostaticMM = CONTROL_instance%ELECTROSTATIC_MM
     LowdinParameters_chargesMM = CONTROL_instance%CHARGES_MM
     LowdinParameters_printMM = CONTROL_instance%PRINT_MM
    !!*****************************************************      
    !! Output Options                               
    LowdinParameters_moldenFile = CONTROL_instance%MOLDEN_FILE                                                                                                                                                                                          
    !!*****************************************************                            
    !! Properties Options                                                              
    LowdinParameters_calculateInterparticleDistances = CONTROL_instance%CALCULATE_INTERPARTICLE_DISTANCES
    LowdinParameters_calculateDensityVolume = CONTROL_instance%CALCULATE_DENSITY_VOLUME
                                                                                                                                                                                          
    !!*****************************************************                            
    !! Miscelaneous Options                                                            
    !!                                                                                 
    LowdinParameters_MOFractionOccupation = CONTROL_instance%MO_FRACTION_OCCUPATION
    LowdinParameters_ionizeMO = CONTROL_instance%IONIZE_MO
    LowdinParameters_ionizeSpecie = CONTROL_instance%IONIZE_SPECIE
    LowdinParameters_exciteSpecie = CONTROL_instance%EXCITE_SPECIE
    LowdinParameters_numberOfCores = CONTROL_instance%NUMBER_OF_CORES

    !!***************************************************************************      
    !! Variables de ambiente al sistema de archivos del programa                       
    !!                                                                                 
    LowdinParameters_homeDirectory = CONTROL_instance%HOME_DIRECTORY
    LowdinParameters_dataDirectory = CONTROL_instance%DATA_DIRECTORY
    LowdinParameters_externalCommand = CONTROL_instance%EXTERNAL_COMMAND
    LowdinParameters_externalSoftwareName = CONTROL_instance%EXTERNAL_SOFTWARE_NAME
    LowdinParameters_uffParametersDataBase = CONTROL_instance%UFF_PARAMETERS_DATABASE
    LowdinParameters_atomicElementsDataBase = CONTROL_instance%ATOMIC_ELEMENTS_DATABASE
    LowdinParameters_basisSetDataBase = CONTROL_instance%BASIS_SET_DATABASE
    LowdinParameters_potentialsDataBase = CONTROL_instance%POTENTIALS_DATABASE
    LowdinParameters_elementalParticlesDataBase = CONTROL_instance%ELEMENTAL_PARTICLES_DATABASE
    LowdinParameters_inputFile = CONTROL_instance%INPUT_FILE
                                                                                                        
    !! Write the name list in the specified unit.
    write(unit, NML=LowdinParameters)
    
  end subroutine CONTROL_save

  !>
  !! @brief Copy the original parameters in other object
  !! @param this
  subroutine CONTROL_copy(this, otherThis)
    implicit none
    type(CONTROL) :: this
    type(CONTROL) :: otherThis
    integer :: i

    !!*****************************************************
    !! Variables para control de integrales
    !!
    otherThis%TV = this%TV 
    otherThis%INTEGRAL_THRESHOLD = this%INTEGRAL_THRESHOLD 
    otherThis%INTEGRAL_DESTINY = this%INTEGRAL_DESTINY 
    otherThis%INTEGRAL_SCHEME = this%INTEGRAL_SCHEME
    otherThis%INTEGRAL_STACK_SIZE = this%INTEGRAL_STACK_SIZE 
    !!***************************************************************************
    !! Parametros para control de proceso de minizacion de energia mediante
    !! metodo SCF
    !!
    otherThis%SCF_NONELECTRONIC_ENERGY_TOLERANCE = this%SCF_NONELECTRONIC_ENERGY_TOLERANCE 
    otherThis%SCF_ELECTRONIC_ENERGY_TOLERANCE = this%SCF_ELECTRONIC_ENERGY_TOLERANCE 
    otherThis%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE = this%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE 
    otherThis%ELECTRONIC_DENSITY_MATRIX_TOLERANCE = this%ELECTRONIC_DENSITY_MATRIX_TOLERANCE 
    otherThis%TOTAL_ENERGY_TOLERANCE = this%TOTAL_ENERGY_TOLERANCE 
    otherThis%STRONG_ENERGY_TOLERANCE = this%STRONG_ENERGY_TOLERANCE 
    otherThis%DENSITY_FACTOR_THRESHOLD = this%DENSITY_FACTOR_THRESHOLD 
    otherThis%DIIS_SWITCH_THRESHOLD = this%DIIS_SWITCH_THRESHOLD 
    otherThis%DIIS_SWITCH_THRESHOLD_BKP = this%DIIS_SWITCH_THRESHOLD_BKP 
    otherThis%ELECTRONIC_LEVEL_SHIFTING = this%ELECTRONIC_LEVEL_SHIFTING 
    otherThis%NONELECTRONIC_LEVEL_SHIFTING = this%NONELECTRONIC_LEVEL_SHIFTING 
    otherThis%WAVE_FUNCTION_SCALE= this%WAVE_FUNCTION_SCALE
    otherThis%SCF_NONELECTRONIC_MAX_ITERATIONS = this%SCF_NONELECTRONIC_MAX_ITERATIONS 
    otherThis%SCF_ELECTRONIC_MAX_ITERATIONS = this%SCF_ELECTRONIC_MAX_ITERATIONS 
    otherThis%SCF_MAX_ITERATIONS= this%SCF_MAX_ITERATIONS
    otherThis%SCF_GLOBAL_MAXIMUM_ITERATIONS = this%SCF_GLOBAL_MAXIMUM_ITERATIONS 
    otherThis%LISTS_SIZE = this%LISTS_SIZE 
    otherThis%CONVERGENCE_METHOD = this%CONVERGENCE_METHOD 
    otherThis%DIIS_DIMENSIONALITY = this%DIIS_DIMENSIONALITY 
    otherThis%ITERATION_SCHEME = this%ITERATION_SCHEME 
    otherThis%SCF_ELECTRONIC_TYPE_GUESS = this%SCF_ELECTRONIC_TYPE_GUESS 
    otherThis%SCF_NONELECTRONIC_TYPE_GUESS = this%SCF_NONELECTRONIC_TYPE_GUESS 
    otherThis%SCF_CONVERGENCE_CRITERIUM = this%SCF_CONVERGENCE_CRITERIUM 
    otherThis%DIIS_ERROR_IN_DAMPING = this%DIIS_ERROR_IN_DAMPING 
    otherThis%ACTIVATE_LEVEL_SHIFTING = this%ACTIVATE_LEVEL_SHIFTING 
    otherThis%EXCHANGE_ORBITALS_IN_SCF = this%EXCHANGE_ORBITALS_IN_SCF 
    otherThis%DEBUG_SCFS = this%DEBUG_SCFS 
    !!***************************************************************************
    !! Parametros para control Hartree-Fock
    !!
    otherThis%FROZEN_PARTICLE = this%FROZEN_PARTICLE 
    otherThis%FREEZE_NON_ELECTRONIC_ORBITALS = this%FREEZE_NON_ELECTRONIC_ORBITALS 
    otherThis%HARTREE_PRODUCT_GUESS = this%HARTREE_PRODUCT_GUESS 
    otherThis%READ_COEFFICIENTS = this%READ_COEFFICIENTS 
    otherThis%NO_SCF = this%NO_SCF 
    otherThis%FINITE_MASS_CORRECTION = this%FINITE_MASS_CORRECTION 
    otherThis%REMOVE_TRANSLATIONAL_CONTAMINATION = this%REMOVE_TRANSLATIONAL_CONTAMINATION 
    otherThis%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE = this%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE 
    otherThis%BUILD_MIXED_DENSITY_MATRIX = this%BUILD_MIXED_DENSITY_MATRIX 
    otherThis%ONLY_ELECTRONIC_EFFECT = this%ONLY_ELECTRONIC_EFFECT 
    otherThis%ELECTRONIC_WAVEFUNCTION_ANALYSIS = this%ELECTRONIC_WAVEFUNCTION_ANALYSIS 
    otherThis%IS_OPEN_SHELL = this%IS_OPEN_SHELL    
    !!***************************************************************************
    !! Parametros para control de proceso de minimizacion multidimensional
    !!
    otherThis%NUMERICAL_DERIVATIVE_DELTA = this%NUMERICAL_DERIVATIVE_DELTA 
    otherThis%MINIMIZATION_INITIAL_STEP_SIZE = this%MINIMIZATION_INITIAL_STEP_SIZE 
    otherThis%MINIMIZATION_LINE_TOLERANCE = this%MINIMIZATION_LINE_TOLERANCE 
    otherThis%MINIMIZATION_TOLERANCE_GRADIENT = this%MINIMIZATION_TOLERANCE_GRADIENT 
    otherThis%MINIMIZATION_MAX_ITERATION = this%MINIMIZATION_MAX_ITERATION 
    otherThis%MINIMIZATION_LIBRARY = this%MINIMIZATION_LIBRARY 
    otherThis%COORDINATES = this%COORDINATES 
    otherThis%ENERGY_CALCULATOR = this%ENERGY_CALCULATOR 
    otherThis%MINIMIZATION_METHOD = this%MINIMIZATION_METHOD 
    otherThis%ANALYTIC_GRADIENT = this%ANALYTIC_GRADIENT 
    otherThis%MINIMIZATION_WITH_SINGLE_POINT = this%MINIMIZATION_WITH_SINGLE_POINT 
    otherThis%USE_SYMMETRY_IN_MATRICES = this%USE_SYMMETRY_IN_MATRICES 
    otherThis%RESTART_OPTIMIZATION = this%RESTART_OPTIMIZATION 
    otherThis%OPTIMIZE_WITH_CP_CORRECTION = this%OPTIMIZE_WITH_CP_CORRECTION 
    otherThis%CP_CORRECTION = this%CP_CORRECTION 
    otherThis%TDHF = this%TDHF 
    otherThis%OPTIMIZE = this%OPTIMIZE 
    otherThis%OPTIMIZE_WITH_MP = this%OPTIMIZE_WITH_MP 
    otherThis%PROJECT_HESSIANE = this%PROJECT_HESSIANE 
    !!***************************************************************************
    !! Criterios de conectividad atomica
    !!
    otherThis%BOND_DISTANCE_FACTOR = this%BOND_DISTANCE_FACTOR 
    otherThis%BOND_ANGLE_THRESHOLD = this%BOND_ANGLE_THRESHOLD 
    otherThis%DIHEDRAL_ANGLE_THRESHOLD  = this%DIHEDRAL_ANGLE_THRESHOLD  
    !!*****************************************************
    !! Control de parametros de teoria de pertubaciones
    !!
    otherThis%MOLLER_PLESSET_CORRECTION = this%MOLLER_PLESSET_CORRECTION 
    otherThis%MP_FROZEN_CORE_BOUNDARY = this%MP_FROZEN_CORE_BOUNDARY 
    otherThis%MP_ONLY_ELECTRONIC_CORRECTION = this%MP_ONLY_ELECTRONIC_CORRECTION 
    !!*****************************************************
    !! Control de parametros de metodo cosmo
    !!
    otherThis%COSMO  = this%COSMO                   
    otherThis%COSMO_SOLVENT_DIALECTRIC = this%COSMO_SOLVENT_DIALECTRIC
    otherThis%COSMO_MIN_BEM = this%COSMO_MIN_BEM           
    otherThis%COSMO_MAX_BEM = this%COSMO_MAX_BEM           
    otherThis%COSMO_RSOLV = this%COSMO_RSOLV             

    !!*****************************************************
    ! Control this for Propagator Theory= !! Control this for Propagator Theory
    !!
    otherThis%PT_ONLY_ONE_SPECIE_CORRECTION = this%PT_ONLY_ONE_SPECIE_CORRECTION 
    otherThis%PT_SELF_ENERGY_SCAN = this%PT_SELF_ENERGY_SCAN 
    otherThis%PT_TRANSITION_OPERATOR = this%PT_TRANSITION_OPERATOR 
    otherThis%PT_JUST_ONE_ORBITAL = this%PT_JUST_ONE_ORBITAL 
    otherThis%PT_SELF_ENERGY_SPACING = this%PT_SELF_ENERGY_SPACING 
    otherThis%PT_SELF_ENERGY_RANGE = this%PT_SELF_ENERGY_RANGE 
    otherThis%PT_MAX_ITERATIONS= this%PT_MAX_ITERATIONS
    otherThis%PT_ITERATION_METHOD_2_LIMIT = this%PT_ITERATION_METHOD_2_LIMIT 
    otherThis%PT_MAX_NUMBER_POLES_SEARCHED = this%PT_MAX_NUMBER_POLES_SEARCHED 
    otherThis%PT_ITERATION_SCHEME = this%PT_ITERATION_SCHEME 
    otherThis%PT_ORDER = this%PT_ORDER 

    !!*****************************************************
    !! Control parametros de formato
    !!
    otherThis%FORMAT_NUMBER_OF_COLUMNS = this%FORMAT_NUMBER_OF_COLUMNS 
    otherThis%UNITS = this%UNITS 
    otherThis%UNIT_FOR_OUTPUT_FILE = this%UNIT_FOR_OUTPUT_FILE 
    otherThis%UNIT_FOR_MP2_INTEGRALS_FILE = this%UNIT_FOR_MP2_INTEGRALS_FILE 
    otherThis%UNIT_FOR_MOLECULAR_ORBITALS_FILE = this%UNIT_FOR_MOLECULAR_ORBITALS_FILE 
    otherThis%PRINT_LEVEL = this%PRINT_LEVEL 
    !!***************************************************************************
    !! CISD - FCI
    !!
    otherThis%CONFIGURATION_INTERACTION_LEVEL = this%CONFIGURATION_INTERACTION_LEVEL 
    !!*****************************************************
    !! Parametros de control general
    !!
    otherThis%DIMENSIONALITY = this%DIMENSIONALITY 
    otherThis%DOUBLE_ZERO_THRESHOLD  = this%DOUBLE_ZERO_THRESHOLD  
    otherThis%METHOD = this%METHOD 
    otherThis%TRANSFORM_TO_CENTER_OF_MASS = this%TRANSFORM_TO_CENTER_OF_MASS 
    otherThis%ARE_THERE_DUMMY_ATOMS = this%ARE_THERE_DUMMY_ATOMS 
    otherThis%IS_THERE_EXTERNAL_POTENTIAL = this%IS_THERE_EXTERNAL_POTENTIAL 
    otherThis%IS_THERE_FROZEN_PARTICLE = this%IS_THERE_FROZEN_PARTICLE 
    !!*****************************************************
    !! Density Functional Theory Options
    !!
    otherThis%AUXILIARY_DENSITY = this%AUXILIARY_DENSITY 
    otherThis%CALL_DFT = this%CALL_DFT 
    otherThis%ELECTRON_CORRELATION_FUNCTIONAL = this%ELECTRON_CORRELATION_FUNCTIONAL 
    otherThis%ELECTRON_EXCHANGE_FUNCTIONAL = this%ELECTRON_EXCHANGE_FUNCTIONAL 
    otherThis%ELECTRON_NUCLEAR_CORRELATION_FUNCTIONAL = this%ELECTRON_NUCLEAR_CORRELATION_FUNCTIONAL 
    otherThis%STORE_THREE_CENTER_ELECTRON_INTEGRALS = this%STORE_THREE_CENTER_ELECTRON_INTEGRALS 
    otherThis%POLARIZATION_ORDER = this%POLARIZATION_ORDER 
    otherThis%FUKUI_FUNCTIONS = this%FUKUI_FUNCTIONS 
    otherThis%NUMBER_OF_BLOCKS_IN_AUXILIARY_FUNCTIONS = this%NUMBER_OF_BLOCKS_IN_AUXILIARY_FUNCTIONS 
    !!*****************************************************
    !! External Potential Options
    !!
    otherThis%NUMERICAL_INTEGRATION_FOR_EXTERNAL_POTENTIAL = this%NUMERICAL_INTEGRATION_FOR_EXTERNAL_POTENTIAL 
    otherThis%NUMERICAL_INTEGRATION_FOR_OVERLAP = this%NUMERICAL_INTEGRATION_FOR_OVERLAP 
    otherThis%MAX_INTERVAL_IN_NUMERICAL_INTEGRATION = this%MAX_INTERVAL_IN_NUMERICAL_INTEGRATION 
    otherThis%STEP_IN_NUMERICAL_INTEGRATION = this%STEP_IN_NUMERICAL_INTEGRATION 
    otherThis%COEFFICIENT_FOR_GAUSSIAN_EXTERNAL_POTENTIAL = this%COEFFICIENT_FOR_GAUSSIAN_EXTERNAL_POTENTIAL 
    otherThis%EXPONENT_FOR_GAUSSIAN_EXTERNAL_POTENTIAL = this%EXPONENT_FOR_GAUSSIAN_EXTERNAL_POTENTIAL 
    otherThis%ORIGIN_OF_GAUSSIAN_EXTERNAL_POTENTIAL = this%ORIGIN_OF_GAUSSIAN_EXTERNAL_POTENTIAL 
    otherThis%NUMERICAL_INTEGRATION_METHOD = this%NUMERICAL_INTEGRATION_METHOD 
    otherThis%RELATIVE_ERROR_IN_NUMERICAL_INTEGRATION = this%RELATIVE_ERROR_IN_NUMERICAL_INTEGRATION 
    otherThis%INITIAL_NUMBER_OF_EVALUATIONS = this%INITIAL_NUMBER_OF_EVALUATIONS 
    otherThis%INCREASE_NUMBER_OF_EVALUATIONS = this%INCREASE_NUMBER_OF_EVALUATIONS 
    otherThis%MINIMUM_NUMBER_OF_EVALUATIONS = this%MINIMUM_NUMBER_OF_EVALUATIONS 
    otherThis%MAXIMUM_NUMBER_OF_EVALUATIONS = this%MAXIMUM_NUMBER_OF_EVALUATIONS 
    !!*****************************************************
    !! Graphs Options
    !!
    otherThis%NUMBER_OF_POINTS_PER_DIMENSION = this%NUMBER_OF_POINTS_PER_DIMENSION 
    !!*****************************************************
    !! Cubes Options
    !!
    otherThis%CUBE_POINTS_DENSITY = this%CUBE_POINTS_DENSITY 
    otherThis%VOLUME_DENSITY_THRESHOLD = this%VOLUME_DENSITY_THRESHOLD 
    !!***************************************************** 
    !! Molecular Mechanics Options                                                        
    otherThis%FORCE_FIELD = this%FORCE_FIELD
    otherThis%ELECTROSTATIC_MM = this%ELECTROSTATIC_MM
    otherThis%CHARGES_MM = this%CHARGES_MM
    otherThis%PRINT_MM = this%PRINT_MM
    !!***************************************************** 
    !! Output Options   
    otherThis%MOLDEN_FILE = this%MOLDEN_FILE
    !!*****************************************************
    !! Properties Options
    otherThis%CALCULATE_INTERPARTICLE_DISTANCES  = this%CALCULATE_INTERPARTICLE_DISTANCES  
    otherThis%CALCULATE_DENSITY_VOLUME  = this%CALCULATE_DENSITY_VOLUME  
    !!*****************************************************
    !! Miscelaneous Options
    !!
    otherThis%MO_FRACTION_OCCUPATION = this%MO_FRACTION_OCCUPATION 
    otherThis%IONIZE_MO = this%IONIZE_MO 
    otherThis%IONIZE_SPECIE = this%IONIZE_SPECIE 
    otherThis%EXCITE_SPECIE = this%EXCITE_SPECIE 
    otherThis%NUMBER_OF_CORES = this%NUMBER_OF_CORES

    !!***************************************************************************
    !! Variables de ambiente al sistema de archivos del programa
    !!
    otherThis%INPUT_FILE = this%INPUT_FILE 
    otherThis%HOME_DIRECTORY = this%HOME_DIRECTORY 
    otherThis%DATA_DIRECTORY = this%DATA_DIRECTORY 
    otherThis%EXTERNAL_COMMAND = this%EXTERNAL_COMMAND 
    otherThis%EXTERNAL_SOFTWARE_NAME = this%EXTERNAL_SOFTWARE_NAME 
    otherThis%UFF_PARAMETERS_DATABASE = this%UFF_PARAMETERS_DATABASE 
    otherThis%ATOMIC_ELEMENTS_DATABASE = this%ATOMIC_ELEMENTS_DATABASE 
    otherThis%BASIS_SET_DATABASE = this%BASIS_SET_DATABASE 
    otherThis%POTENTIALS_DATABASE = this%POTENTIALS_DATABASE 
    otherThis%ELEMENTAL_PARTICLES_DATABASE = this%ELEMENTAL_PARTICLES_DATABASE 

  end subroutine CONTROL_copy

  !>
  !! @brief Shows running parameters
  !! @author S. A. Gonzalez
  subroutine CONTROL_show()
    implicit none

    integer:: i
    
    print *,""
    print *,"LOWDIN IS RUNNING WITH NEXT PARAMETERS: "
    print *,"----------------------------------------"
    print *,""
    
    write (*,"(T10,A)") "METHOD TYPE:  "//trim(CONTROL_instance%METHOD)
    write (*,"(T10,A,I5)") "NUMBER OF CORES: ",CONTROL_instance%NUMBER_OF_CORES
    
    
    if(CONTROL_instance%METHOD=="RKS" .or. CONTROL_instance%METHOD=="UKS" .or. CONTROL_instance%METHOD=="ROKS" ) then

       if(CONTROL_instance%AUXILIARY_DENSITY) write (*,"(T10,A)") "USING AUXILIARY DENSITY"
       
       write (*,"(T10,A)") "ELECTRON CORRELATION FUNCTIONAL: "//trim(CONTROL_instance%ELECTRON_CORRELATION_FUNCTIONAL)
       write (*,"(T10,A)") "ELECTRON EXCHANGE FUNCTIONAL: "//trim(CONTROL_instance%ELECTRON_EXCHANGE_FUNCTIONAL)
       write (*,"(T10,A)") "ELECTRON-NUCLEAR CORRELATION FUNCTIONAL: "//trim(CONTROL_instance%ELECTRON_NUCLEAR_CORRELATION_FUNCTIONAL)

       if(CONTROL_instance%STORE_THREE_CENTER_ELECTRON_INTEGRALS) then

          write (*,"(T10,A)") "STORING THREE CENTER ELECTRON INTEGRALS IN DISK"

       else

          write (*,"(T10,A)") "CALCULATING THREE CENTER ELECTRON INTEGRALS ON THE FLY"

       end if

    end if

    if(CONTROL_instance%METHOD=="MM") then
       write (*,"(T10,A)") "FORCE FIELD: "//trim(CONTROL_instance%FORCE_FIELD)
       if(CONTROL_instance%ELECTROSTATIC_MM) then
          write (*,"(T10,A)") "CALCULATE ELECTROSTATIC ENERGY: YES"
          write (*,"(T10,A)") "PARTIAL CHARGES METHOD: EQeq(Extended Charge Equilibration)"
       else
          write (*,"(T10,A)") "CALCULATE ELECTROSTATIC ENERGY: NO"
       end if
    end if

    if(CONTROL_instance%MOLLER_PLESSET_CORRECTION>=2) then

       write (*,"(T10,A,I5)") "MOLLER PLESSET CORRECTION:  ",CONTROL_instance%MOLLER_PLESSET_CORRECTION

    end if

    if(CONTROL_instance%COSMO) then

       write (*,"(T10,A)") "COSMO:  T "

    end if

    if(CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE" ) then

       write (*,"(T10,A,A)") "CONFIGURATION INTERACTION LEVEL:  ", CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL

    end if


    if(CONTROL_instance%PT_ORDER>=2) then

       write (*,"(T10,A,I5)") "PROPAGATOR THEORY ORDER:  ",CONTROL_instance%PT_ORDER

    end if

    if((CONTROL_instance%IONIZE_SPECIE(1)) /= "NONE") then 

       write (*,"(T10,A,I5)") "MOLECULAR ORBITAL TO BE IONIZED: ", CONTROL_instance%IONIZE_MO
       do i = 1, size(CONTROL_instance%IONIZE_SPECIE)
         if ( CONTROL_instance%IONIZE_SPECIE(i) /= "NONE" ) then
         write (*,"(T10,A,A)") "FOR SPECIE0: ", (CONTROL_instance%IONIZE_SPECIE(i))
         write (*,"(T10,A,ES15.5)") "IONIZED MOLECULAR ORBITAL OCCUPATION: ",CONTROL_instance%MO_FRACTION_OCCUPATION
         end if
      end do 
    end if

    if(CONTROL_instance%POLARIZATION_ORDER > 1) then 

       write (*,"(T10,A,I5)") "POLARIZATION ORDER TO BE CALCULATED: ", CONTROL_instance%POLARIZATION_ORDER

    end if

    if(CONTROL_instance%FUKUI_FUNCTIONS) then 

       write (*,"(T10,A,I5)") "CALCULATING FUKUI FUNCTIONS"

    end if

    if(CONTROL_instance%EXCITE_SPECIE /= "NONE") then 

       write (*,"(T10,A,T10)") "CALCULATING", trim(CONTROL_instance%EXCITE_SPECIE) ,"IN THE FIRST EXCITED STATE"

    end if

    if(CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE) then 

       write (*,"(T10,A)") "BUILDING TWO PARTICLES MATRIX FOR ONE PARTICLE"

    end if

    if (CONTROL_instance%OPTIMIZE) then

       write (*,"(T10,A)") "GEOMETRY OPTIMIZATION:  T"
       write (*,"(T10,A)") "OPTIMIZATION METHOD: "//trim(CONTROL_instance%MINIMIZATION_METHOD)
       write (*,"(T10,A,ES15.5)") "GRADIENT THRESHOLD FOR THE MINIMIZATION: ",CONTROL_instance%MINIMIZATION_TOLERANCE_GRADIENT

       if(CONTROL_instance%MOLLER_PLESSET_CORRECTION>=2) then

          CONTROL_instance%ANALYTIC_GRADIENT=.false.

       end if

       if( CONTROL_instance%ANALYTIC_GRADIENT) then

          write (*,"(T10,A)") "TYPE OF GRADIENT : ANALYTIC"

       else

          write (*,"(T10,A)") "TYPE OF GRADIENT : NUMERIC"
          write (*,"(T10,A,E15.5)") "NUMERICAL STEP: ",CONTROL_instance%NUMERICAL_DERIVATIVE_DELTA

       end if

    else if (CONTROL_instance%TDHF) then

       write (*,"(T10,A)") "TIME DEPENDENT HATREE-FOCK:  T"
       write (*,"(T10,A)") "EVOLUTION METHOD: "//trim(CONTROL_instance%MINIMIZATION_METHOD)
       write (*,"(T10,A,ES15.5)") "GRADIENT THRESHOLD FOR THE MINIMIZATION: ",CONTROL_instance%MINIMIZATION_TOLERANCE_GRADIENT

       if(CONTROL_instance%MOLLER_PLESSET_CORRECTION>=2) then

          CONTROL_instance%ANALYTIC_GRADIENT=.false.

       end if

       if( CONTROL_instance%ANALYTIC_GRADIENT) then

          write (*,"(T10,A)") "TYPE OF GRADIENT : ANALYTIC"

       else

          write (*,"(T10,A)") "TYPE OF GRADIENT : NUMERIC"
          write (*,"(T10,A,E15.5)") "NUMERICAL STEP: ",CONTROL_instance%NUMERICAL_DERIVATIVE_DELTA

       end if

       write (*,"(T10,A,A)") "SINGLE POINT CALCULATION"

    end if

    if(CONTROL_instance%METHOD/="MM") then
       write (*,"(T10,A,E15.5)") "NONELECTRONIC ENERGY TOLERANCE IN SCFs: ",CONTROL_instance%SCF_NONELECTRONIC_ENERGY_TOLERANCE
       write (*,"(T10,A,E15.5)") "NONELECTRONIC DENSITY MATRIX TOLERANCE IN SCFs: ",CONTROL_instance%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE
       write (*,"(T10,A,E15.5)") "ELECTRONIC ENERGY TOLERANCE IN SCFs: ",CONTROL_instance%SCF_ELECTRONIC_ENERGY_TOLERANCE
       write (*,"(T10,A,E15.5)") "ELECTRONIC DENSITY MATRIX TOLERANCE IN SCFs: ",CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE
       write (*,"(T10,A,E15.5)") "TOTAL ENERGY TOLERANCE IN SCFs: ",CONTROL_instance%TOTAL_ENERGY_TOLERANCE
       write (*,"(T10,A,I5)") "SCF MAX. ITERATIONS - NONELECTRONICS : ",CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS
       write (*,"(T10,A,I5)") "SCF MAX. ITERATIONS - ELECTRONICS : ",CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS
       write (*,"(T10,A,I5)") "SCF MAX. ITERATIONS - INTERSPECIES : ",CONTROL_instance%SCF_GLOBAL_MAXIMUM_ITERATIONS
       write (*,"(T10,A)") "CRITERIUM OF CONVERGENCE: "//trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM)
       write (*,"(T10,A)") "NONELECTRONIC DENSITY GUESS: "//trim(CONTROL_instance%SCF_NONELECTRONIC_TYPE_GUESS)
       write (*,"(T10,A)") "ELECTRONIC DENSITY GUESS: "//trim(CONTROL_instance%SCF_ELECTRONIC_TYPE_GUESS)
    end if

    if (CONTROL_instance%NO_SCF) write (*,"(T10,A)") "NO SCF WILL BE PERFORMED"

    if (CONTROL_instance%FREEZE_NON_ELECTRONIC_ORBITALS) write (*,"(T10,A)") "Electrons will be frozen during SCF calculation"

    if ( CONTROL_instance%HARTREE_PRODUCT_GUESS) then

       write (*,"(T10,A)") "HARTREE PRODUCT GUESS: T"

    end if

    if(CONTROL_instance%METHOD/="MM") then
       write (*,"(T10,A,I5)") "SCHEME OF ITERATION: ",CONTROL_instance%ITERATION_SCHEME
       write (*,"(T10,A)") "INTEGRAL DESTINY: "//trim(CONTROL_instance%INTEGRAL_DESTINY)
       write (*,"(T10,A,I5)") "STACK SIZE FOR ERIS : ", CONTROL_instance%INTEGRAL_STACK_SIZE

       select case(CONTROL_instance%CONVERGENCE_METHOD)

       case(1)

          write(*,"(T10,A)") "METHOD OF SCF CONVERGENCE: DAMPING"

       case(2)

          write(*,"(T10,A)") "METHOD OF SCF CONVERGENCE: DIIS"
          write(*,"(T10,A,E15.5)") "DIIS SWITCH THRESHOLD", CONTROL_instance%DIIS_SWITCH_THRESHOLD
          write(*,"(T10,A,I5)") "DIIS DIMENSIONALITY: ", CONTROL_instance%DIIS_DIMENSIONALITY

       case(4)

          write(*,"(T10,A)") "METHOD OF SCF CONVERGENCE: DAMPING/DIIS"
          write(*,"(T10,A,E15.5)") "DIIS SWITCH THRESHOLD", CONTROL_instance%DIIS_SWITCH_THRESHOLD
          write(*,"(T10,A,I5)") "DIIS DIMENSIONALITY: ", CONTROL_instance%DIIS_DIMENSIONALITY

       end select
    end if

  end subroutine CONTROL_show

  !>
  !!  @brief returs root directory of LOWDIN
  function CONTROL_getHomeDirectory() result ( directory )
    implicit none
    character(255):: directory

    call getenv( "LOWDIN_HOME", directory )

    directory = trim( directory )

  end function CONTROL_getHomeDirectory

  !>
  !!  @brief returns library path
  function CONTROL_getDataDirectory() result ( directory )
    implicit none

    character(255):: directory

    call getenv( "LOWDIN_DATA", directory )

    directory = trim( directory )

  end function CONTROL_getDataDirectory

  !>
  !!  @brief returns some environment variables
  function CONTROL_getExternalCommand() result ( output )
    implicit none
    character(255):: output

    call getenv( "EXTERNAL_COMMAND", output )

    output = trim( output )

  end function CONTROL_getExternalCommand

  !>
  !!  @brief retorna la ruta al directorio raiz de LOWDIN
  function CONTROL_getExternalSoftwareName() result ( output )
    implicit none
    character(255):: output

    call getenv( "EXTERNAL_SOFTWARE_NAME", output )

    output = trim( output )

  end function CONTROL_getExternalSoftwareName

  !>
  !! @brief  handle exceptions
  subroutine CONTROL_exception( typeMessage, description, debugDescription)
    implicit none
    integer:: typeMessage
    character(*):: description
    character(*):: debugDescription

    type(Exception):: ex

    call Exception_constructor( ex , typeMessage )
    call Exception_setDebugDescription( ex, debugDescription )
    call Exception_setDescription( ex, description )
    call Exception_show( ex )
    call Exception_destructor( ex )

  end subroutine CONTROL_exception
  
end module CONTROL_
