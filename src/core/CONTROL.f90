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
  use Units_
  use Exception_
  use omp_lib
  implicit none

  type, public :: CONTROL

     !!***************************************************************************
     !! Dummy variables, just for debugging. 
     !!
     real(8) :: DUMMY_REAL(10)
     integer :: DUMMY_INTEGER(10)
     logical :: DUMMY_LOGICAL(10)
     character(50) :: DUMMY_CHARACTER(10)

     !!***************************************************************************
     !! Parameter to control Integrals library
     !!
     real(8) :: TV
     real(8) :: INTEGRAL_THRESHOLD
     integer :: INTEGRAL_STACK_SIZE
     character(20) :: INTEGRAL_STORAGE
     character(20) :: INTEGRAL_SCHEME
     logical :: SCHWARZ_INEQUALITY
     real(8) :: HARMONIC_CONSTANT

     !!***************************************************************************
     !! Parameter to control SCF program
     !!
     real(8) :: NONELECTRONIC_ENERGY_TOLERANCE
     real(8) :: ELECTRONIC_ENERGY_TOLERANCE
     real(8) :: NONELECTRONIC_DENSITY_MATRIX_TOLERANCE
     real(8) :: ELECTRONIC_DENSITY_MATRIX_TOLERANCE
     real(8) :: TOTAL_ENERGY_TOLERANCE
     real(8) :: TOTAL_DENSITY_MATRIX_TOLERANCE
     real(8) :: DENSITY_FACTOR_THRESHOLD !< define cuando recalcula un elemeto Gij de acuedo con el valor de Pij
     real(8) :: DIIS_SWITCH_THRESHOLD
     real(8) :: DIIS_SWITCH_THRESHOLD_BKP
     real(8) :: ELECTRONIC_LEVEL_SHIFTING
     real(8) :: NONELECTRONIC_LEVEL_SHIFTING
     real(8) :: EXCHANGE_ORBITAL_THRESHOLD
     real(8) :: WAVE_FUNCTION_SCALE
     integer :: SCF_NONELECTRONIC_MAX_ITERATIONS
     integer :: SCF_ELECTRONIC_MAX_ITERATIONS
     integer :: SCF_GLOBAL_MAX_ITERATIONS
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
     logical :: FORCE_CLOSED_SHELL
     logical :: DEBUG_SCFS
     character(10) ::  SCF_GHOST_SPECIES

     !!***************************************************************************
     !! Hartree-Fock options
     !!
     character(20) :: FROZEN_PARTICLE(5)
     logical :: FREEZE_NON_ELECTRONIC_ORBITALS
     logical :: FREEZE_ELECTRONIC_ORBITALS
     logical :: HARTREE_PRODUCT_GUESS
     logical :: READ_COEFFICIENTS
     logical :: READ_FCHK
     logical :: WRITE_COEFFICIENTS_IN_BINARY
     logical :: READ_EIGENVALUES
     logical :: READ_EIGENVALUES_IN_BINARY
     logical :: WRITE_EIGENVALUES_IN_BINARY
     logical :: NO_SCF
     logical :: FINITE_MASS_CORRECTION
     logical :: REMOVE_TRANSLATIONAL_CONTAMINATION
     logical :: BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE
     logical :: BUILD_MIXED_DENSITY_MATRIX
     logical :: ONLY_ELECTRONIC_EFFECT
     logical :: ELECTRONIC_WAVEFUNCTION_ANALYSIS
     logical :: IS_OPEN_SHELL
     logical :: GET_GRADIENTS
     logical :: HF_PRINT_EIGENVALUES
     character(20) :: HF_PRINT_EIGENVECTORS
     real(8) :: OVERLAP_EIGEN_THRESHOLD
     real(8) :: ELECTRIC_FIELD(3)
     integer :: MULTIPOLE_ORDER

     !!***************************************************************************
     !! Parameter to control geometry optimization
     !!
     real(8) :: NUMERICAL_DERIVATIVE_DELTA
     real(8) :: MINIMIZATION_INITIAL_STEP_SIZE
     real(8) :: MINIMIZATION_LINE_TOLERANCE
     real(8) :: MINIMIZATION_TOLERANCE_GRADIENT
     integer :: MINIMIZATION_MAX_ITERATION
     integer :: MINIMIZATION_METHOD
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
     logical :: LAST_STEP
     logical :: OPTIMIZE_WITH_MP
     logical :: PROJECT_HESSIANE

     !!***************************************************************************
     !! Parameter of atomic conectivity
     !!
     real(8) :: BOND_DISTANCE_FACTOR
     real(8) :: BOND_ANGLE_THRESHOLD
     real(8) :: DIHEDRAL_ANGLE_THRESHOLD

     !!***************************************************************************
     !! Parameter to control MBP theory
     !!
     integer :: MOLLER_PLESSET_CORRECTION
     integer :: MP_FROZEN_CORE_BOUNDARY
     logical :: MP_ONLY_ELECTRONIC_CORRECTION

     integer :: EPSTEIN_NESBET_CORRECTION

     !!***************************************************************************
     !! Parameter to control cosmo  
     !!
     logical :: COSMO
     real(8) :: COSMO_SOLVENT_DIELECTRIC
     real(8) :: COSMO_SCALING

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
     real(8) :: PT_FACTOR_SS
     real(8) :: PT_FACTOR_OS
     character(10) :: PT_P3_METHOD(7)


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
     !! CI
     !!
     character(20) :: CONFIGURATION_INTERACTION_LEVEL
     integer :: NUMBER_OF_CI_STATES
     character(20) :: CI_DIAGONALIZATION_METHOD
     character(20) :: CI_PRINT_EIGENVECTORS_FORMAT
     character(20) :: CI_DIAGONAL_DRESSED_SHIFT
     real(8) :: CI_PRINT_THRESHOLD
     integer :: CI_STATES_TO_PRINT
     integer :: CI_ACTIVE_SPACE
     integer :: CI_MAX_NCV
     integer :: CI_SIZE_OF_GUESS_MATRIX
     integer :: CI_STACK_SIZE
     real(8) :: CI_CONVERGENCE
     real(8) :: CI_MATVEC_TOLERANCE
     logical :: CI_SAVE_EIGENVECTOR
     logical :: CI_LOAD_EIGENVECTOR
     logical :: CI_JACOBI
     logical :: CI_BUILD_FULL_MATRIX
     integer :: CI_MADSPACE
     logical :: CI_NATURAL_ORBITALS
     integer :: CI_SCI_CORE_SPACE
     integer :: CI_SCI_TARGET_SPACE

     !!***************************************************************************
     !! Non-orthogonal CI
     !!
     logical :: NONORTHOGONAL_CONFIGURATION_INTERACTION
     integer :: TRANSLATION_SCAN_GRID(3)
     integer :: ROTATIONAL_SCAN_GRID
     integer :: NESTED_ROTATIONAL_GRIDS
     real(8) :: ROTATION_AROUND_Z_MAX_ANGLE
     real(8) :: ROTATION_AROUND_Z_STEP
     real(8) :: TRANSLATION_STEP
     real(8) :: NESTED_GRIDS_DISPLACEMENT
     real(8) :: CONFIGURATION_ENERGY_THRESHOLD
     real(8) :: CONFIGURATION_OVERLAP_THRESHOLD
     real(8) :: CONFIGURATION_MAX_DISPLACEMENT(3)
     real(8) :: CONFIGURATION_MIN_DISPLACEMENT(3)
     real(8) :: CONFIGURATION_MAX_NP_DISTANCE
     real(8) :: CONFIGURATION_MIN_PP_DISTANCE
     real(8) :: CONFIGURATION_MAX_PP_DISTANCE
     real(8) :: CONFIGURATION_EQUIVALENCE_DISTANCE
     real(8) :: EMPIRICAL_OVERLAP_PARAMETER_A
     real(8) :: EMPIRICAL_OVERLAP_PARAMETER_B
     real(8) :: EMPIRICAL_OVERLAP_PARAMETER_E0
     real(8) :: EMPIRICAL_OVERLAP_PARAMETER_SC
     logical :: CONFIGURATION_USE_SYMMETRY
     logical :: READ_NOCI_GEOMETRIES
     logical :: EMPIRICAL_OVERLAP_CORRECTION
     logical :: ONLY_FIRST_NOCI_ELEMENTS
     logical :: NOCI_KINETIC_APPROXIMATION
     logical :: COMPUTE_ROCI_FORMULA
     logical :: REMOVE_QDO_IN_CI

     !!***************************************************************************
     !! CCSD Parameters
     !!
     character(20) :: COUPLED_CLUSTER_LEVEL

     !!*****************************************************
     !! Parameter to general control
     !!
     character(50) :: METHOD
     logical :: TRANSFORM_TO_CENTER_OF_MASS
     logical :: ARE_THERE_DUMMY_ATOMS
     logical :: ARE_THERE_QDO_POTENTIALS
     logical :: SET_QDO_ENERGY_ZERO
     logical :: IS_THERE_EXTERNAL_POTENTIAL
     logical :: IS_THERE_INTERPARTICLE_POTENTIAL
     logical :: IS_THERE_OUTPUT
     logical :: IS_THERE_FROZEN_PARTICLE
     integer :: DIMENSIONALITY

     !!*****************************************************
     !! Density Functional Theory Options
     !!
     character(50) :: GRID_STORAGE
     character(50) :: ELECTRON_CORRELATION_FUNCTIONAL
     character(50) :: ELECTRON_EXCHANGE_FUNCTIONAL
     character(50) :: ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL
     character(50) :: NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL
     character(50) :: POSITRON_ELECTRON_CORRELATION_FUNCTIONAL
     character(50) :: BETA_FUNCTION
     integer :: GRID_RADIAL_POINTS
     integer :: GRID_ANGULAR_POINTS
     integer :: GRID_NUMBER_OF_SHELLS
     integer :: FINAL_GRID_RADIAL_POINTS
     integer :: FINAL_GRID_ANGULAR_POINTS
     integer :: FINAL_GRID_NUMBER_OF_SHELLS
     integer :: POLARIZATION_ORDER
     integer :: NUMBER_OF_BLOCKS_IN_AUXILIARY_FUNCTIONS
     logical :: FUKUI_FUNCTIONS
     logical :: AUXILIARY_DENSITY
     logical :: STORE_THREE_CENTER_ELECTRON_INTEGRALS
     logical :: CALL_LIBXC
     real(8) :: NUCLEAR_ELECTRON_DENSITY_THRESHOLD
     real(8) :: ELECTRON_DENSITY_THRESHOLD
     real(8) :: GRID_WEIGHT_THRESHOLD
     real(8) :: BETA_PARAMETER_A
     real(8) :: BETA_PARAMETER_B
     real(8) :: BETA_PARAMETER_C

     !!*****************************************************
     !! Subsystem embedding Options
     !!
     logical :: SUBSYSTEM_EMBEDDING
     logical :: LOCALIZE_ORBITALS
     real(8) :: SUBSYSTEM_LEVEL_SHIFTING
     real(8) :: SUBSYSTEM_ORBITAL_THRESHOLD
     real(8) :: SUBSYSTEM_BASIS_THRESHOLD
     character(50) :: ERKALE_LOCALIZATION_METHOD

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
     character(50) :: MOLDEN_FILE_FORMAT

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
     logical :: AMBER_FILE

     !!*****************************************************
     !! Properties Options
     logical :: CALCULATE_INTERPARTICLE_DISTANCES
     logical :: CALCULATE_DENSITY_VOLUME

     !!*****************************************************
     !! Miscelaneous Options
     !!
     real(8) :: MO_FRACTION_OCCUPATION(10)
     integer :: IONIZE_MO(10)
     character(50) :: IONIZE_SPECIES(10)
     character(50) :: EXCITE_SPECIES
     integer :: NUMBER_OF_CORES

     !!*****************************************************
     !! Integrals transformation options
     !!
     character(10) :: INTEGRALS_TRANSFORMATION_METHOD
     integer :: IT_BUFFERSIZE

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
  !! Dummy variables, just for debugging. 
  !!
  real(8) :: LowdinParameters_dummyReal(10)
  integer :: LowdinParameters_dummyInteger(10)
  logical :: LowdinParameters_dummyLogical(10)
  character(50) :: LowdinParameters_dummyCharacter(10)

  !!***************************************************************************
  !! Parameter to control Integrals library
  !!  
  real(8) :: LowdinParameters_tv
  real(8) :: LowdinParameters_integralThreshold
  integer :: LowdinParameters_integralStackSize
  character(20) :: LowdinParameters_integralStorage
  character(20) :: LowdinParameters_integralScheme
  logical :: LowdinParameters_schwarzInequality
  real(8) :: LowdinParameters_harmonicConstant

  !!***************************************************************************
  !! Parameter to control SCF program
  !!
  real(8) :: LowdinParameters_nonElectronicEnergyTolerance
  real(8) :: LowdinParameters_electronicEnergyTolerance
  real(8) :: LowdinParameters_nonelectronicDensityMatrixTolerance
  real(8) :: LowdinParameters_electronicDensityMatrixTolerance
  real(8) :: LowdinParameters_totalEnergyTolerance
  real(8) :: LowdinParameters_totalDensityMatrixTolerance
  real(8) :: LowdinParameters_densityFactorThreshold
  real(8) :: LowdinParameters_diisSwitchThreshold
  real(8) :: LowdinParameters_diisSwitchThreshold_bkp
  real(8) :: LowdinParameters_electronicLevelShifting
  real(8) :: LowdinParameters_nonelectronicLevelShifting
  real(8) :: LowdinParameters_exchangeOrbitalThreshold
  real(8) :: LowdinParameters_waveFunctionScale
  integer :: LowdinParameters_scfNonelectronicMaxIterations
  integer :: LowdinParameters_scfElectronicMaxIterations
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
  logical :: LowdinParameters_forceClosedShell
  logical :: LowdinParameters_debugScfs
  character(10) :: LowdinParameters_scfGhostSpecies

  !!*****************************************************
  !! Hartree-Fock Options
  !!
  character(20) :: LowdinParameters_frozen(5)
  logical :: LowdinParameters_freezeNonElectronicOrbitals
  logical :: LowdinParameters_freezeElectronicOrbitals
  logical :: LowdinParameters_hartreeProductGuess
  logical :: LowdinParameters_readCoefficients
  logical :: LowdinParameters_readFchk
  logical :: LowdinParameters_writeCoefficientsInBinary
  logical :: LowdinParameters_readCoefficientsInBinary
  logical :: LowdinParameters_readEigenvalues
  logical :: LowdinParameters_readEigenvaluesInBinary
  logical :: LowdinParameters_writeEigenvaluesInBinary
  logical :: LowdinParameters_noSCF
  logical :: LowdinParameters_finiteMassCorrection
  logical :: LowdinParameters_removeTranslationalContamination
  logical :: LowdinParameters_buildTwoParticlesMatrixForOneParticle
  logical :: LowdinParameters_buildMixedDensityMatrix
  logical :: LowdinParameters_onlyElectronicEffect
  logical :: LowdinParameters_electronicWaveFunctionAnalysis
  logical :: LowdinParameters_isOpenShell
  logical :: LowdinParameters_getGradients
  logical :: LowdinParameters_HFprintEigenvalues
  character(20) :: LowdinParameters_HFprintEigenvectors
  real(8) :: LowdinParameters_overlapEigenThreshold
  real(8) :: LowdinParameters_electricField(3)
  integer :: LowdinParameters_multipoleOrder

  !!***************************************************************************
  !! Parameter to control geometry optimization
  !!
  real(8) :: LowdinParameters_numericalDerivativeDelta
  real(8) :: LowdinParameters_minimizationInitialStepSize
  real(8) :: LowdinParameters_minimizationLineTolerance
  real(8) :: LowdinParameters_minimizationToleranceGradient
  integer :: LowdinParameters_minimizationMaxIteration
  integer :: LowdinParameters_minimizationMethod
  character(10) :: LowdinParameters_minimizationLibrary
  character(50) :: LowdinParameters_coordinates
  character(20) :: LowdinParameters_energyCalculator
  logical :: LowdinParameters_analyticGradient
  logical :: LowdinParameters_minimizationWithSinglePoint
  logical :: LowdinParameters_useSymmetryInMatrices
  logical :: LowdinParameters_restartOptimization
  logical :: LowdinParameters_firstStep
  logical :: LowdinParameters_lastStep
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
  !! Parameter to control MBPn theory
  !!
  integer :: LowdinParameters_mpCorrection
  integer :: LowdinParameters_mpFrozenCoreBoundary
  logical :: LowdinParameters_mpOnlyElectronicCorrection
  integer :: LowdinParameters_epsteinNesbetCorrection 

  !!***************************************************************************
  !! Parameter to control cosmo theory
  !!
  logical :: LowdinParameters_cosmo
  real(8) :: LowdinParameters_cosmo_solvent_dielectric
  real(8) :: LowdinParameters_cosmo_scaling

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
  real(8) :: LowdinParameters_ptFactorSS 
  real(8) :: LowdinParameters_ptFactorOS 
  character(10) :: LowdinParameters_ptP3Method(7)


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
  integer :: LowdinParameters_numberOfCIStates
  character(20) :: LowdinParameters_CIdiagonalizationMethod
  character(20) :: LowdinParameters_CIPrintEigenVectorsFormat
  character(20) :: LowdinParameters_CIdiagonalDressedShift 
  real(8) :: LowdinParameters_CIPrintThreshold
  integer :: LowdinParameters_CIactiveSpace
  integer :: LowdinParameters_CIstatesToPrint
  integer :: LowdinParameters_CImaxNCV
  integer :: LowdinParameters_CIsizeOfGuessMatrix
  integer :: LowdinParameters_CIstackSize
  real(8) :: LowdinParameters_CIConvergence
  real(8) :: LowdinParameters_CImatvecTolerance
  logical :: LowdinParameters_CISaveEigenVector
  logical :: LowdinParameters_CILoadEigenVector
  logical :: LowdinParameters_CIJacobi
  logical :: LowdinParameters_CIBuildFullMatrix
  integer :: LowdinParameters_CIMadSpace
  logical :: LowdinParameters_CINaturalOrbitals
  integer :: LowdinParameters_CISCICoreSpace 
  integer :: LowdinParameters_CISCITargetSpace 

  !!***************************************************************************
  !! Non-orthogonal CI
  !!
  logical :: LowdinParameters_nonOrthogonalConfigurationInteraction
  integer :: LowdinParameters_translationScanGrid(3)
  integer :: LowdinParameters_rotationalScanGrid
  integer :: LowdinParameters_nestedRotationalGrids
  real(8) :: LowdinParameters_rotationAroundZMaxAngle
  real(8) :: LowdinParameters_rotationAroundZStep
  real(8) :: LowdinParameters_translationStep
  real(8) :: LowdinParameters_nestedGridsDisplacement
  real(8) :: LowdinParameters_configurationEnergyThreshold
  real(8) :: LowdinParameters_configurationOverlapThreshold
  real(8) :: LowdinParameters_configurationMaxDisplacement(3)
  real(8) :: LowdinParameters_configurationMinDisplacement(3)
  real(8) :: LowdinParameters_configurationMaxNPDistance
  real(8) :: LowdinParameters_configurationMinPPDistance
  real(8) :: LowdinParameters_configurationMaxPPDistance
  real(8) :: LowdinParameters_configurationEquivalenceDistance
  real(8) :: LowdinParameters_empiricalOverlapParameterA
  real(8) :: LowdinParameters_empiricalOverlapParameterB
  real(8) :: LowdinParameters_empiricalOverlapParameterE0
  real(8) :: LowdinParameters_empiricalOverlapParameterSc
  logical :: LowdinParameters_configurationUseSymmetry
  logical :: LowdinParameters_readNOCIGeometries
  logical :: LowdinParameters_empiricalOverlapCorrection
  logical :: LowdinParameters_onlyFirstNOCIelements
  logical :: LowdinParameters_NOCIKineticApproximation
  logical :: LowdinParameters_computeROCIformula
  logical :: LowdinParameters_removeQDOinCI

  !!***************************************************************************
  !! CCSD
  !! 
  character(20) :: LowdinParameters_coupledClusterLevel

  !!*****************************************************
  !! Parameter to general control
  !!
  character(50) :: LowdinParameters_method
  logical :: LowdinParameters_transformToCenterOfMass
  logical :: LowdinParameters_areThereDummyAtoms
  logical :: LowdinParameters_areThereQDOPotentials
  logical :: LowdinParameters_setQDOEnergyZero
  logical :: LowdinParameters_isThereExternalPotential
  logical :: LowdinParameters_isThereInterparticlePotential
  logical :: LowdinParameters_isThereOutput
  logical :: LowdinParameters_isThereFrozenParticle
  integer :: LowdinParameters_dimensionality

  !!*****************************************************
  !! Density Functional Theory Options
  !!
  character(50) :: LowdinParameters_gridStorage
  character(50) :: LowdinParameters_electronCorrelationFunctional
  character(50) :: LowdinParameters_electronExchangeFunctional
  character(50) :: LowdinParameters_electronExchangeCorrelationFunctional
  character(50) :: LowdinParameters_nuclearElectronCorrelationFunctional
  character(50) :: LowdinParameters_positronElectronCorrelationFunctional
  character(50) :: LowdinParameters_betaFunction
  integer :: LowdinParameters_gridRadialPoints
  integer :: LowdinParameters_gridAngularPoints
  integer :: LowdinParameters_gridNumberOfShells
  integer :: LowdinParameters_finalGridRadialPoints
  integer :: LowdinParameters_finalGridAngularPoints
  integer :: LowdinParameters_finalGridNumberOfShells
  integer :: LowdinParameters_polarizationOrder
  integer :: LowdinParameters_numberOfBlocksInAuxiliaryFunctions
  logical :: LowdinParameters_fukuiFunctions
  logical :: LowdinParameters_auxiliaryDensity
  logical :: LowdinParameters_storeThreeCenterElectronIntegrals
  logical :: LowdinParameters_callLibxc
  real(8) :: LowdinParameters_nuclearElectronDensityThreshold
  real(8) :: LowdinParameters_electronDensityThreshold
  real(8) :: LowdinParameters_gridWeightThreshold
  real(8) :: LowdinParameters_betaParameterA
  real(8) :: LowdinParameters_betaParameterB
  real(8) :: LowdinParameters_betaParameterC

  !!*****************************************************
  !! Subsystem embedding Options
  !!
  logical :: LowdinParameters_subsystemEmbedding
  logical :: LowdinParameters_localizeOrbitals
  real(8) :: LowdinParameters_subsystemLevelShifting
  real(8) :: LowdinParameters_subsystemOrbitalThreshold
  real(8) :: LowdinParameters_subsystemBasisThreshold
  character(50) :: LowdinParameters_erkaleLocalizationMethod

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
  character(50) :: LowdinParameters_moldenFileFormat

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
  logical :: LowdinParameters_amberFile

  !!*****************************************************
  !! Properties Options
  logical :: LowdinParameters_calculateInterparticleDistances
  logical :: LowdinParameters_calculateDensityVolume

  !!*****************************************************
  !! Miscelaneous Options
  !!
  real(8) :: LowdinParameters_MOFractionOccupation(10)
  integer :: LowdinParameters_ionizeMO(10)
  character(50) :: LowdinParameters_ionizeSpecies(10)
  character(50) :: LowdinParameters_exciteSpecies

  !!*****************************************************
  !! Integrals transformation options
  !!
  character(10) :: LowdinParameters_integralsTransformationMethod
  integer :: LowdinParameters_ITBuffersize

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
                                !! Dummy variables, just for debugging. 
                                !!
       LowdinParameters_dummyReal,&
       LowdinParameters_dummyInteger,&
       LowdinParameters_dummyLogical,&
       LowdinParameters_dummyCharacter,&
       
       
                                !!***************************************************************************
                                !! Parameter to control Integrals library
                                !!  
       LowdinParameters_tv,&
       LowdinParameters_integralThreshold,&
       LowdinParameters_integralStackSize,&
       LowdinParameters_integralStorage,&
       LowdinParameters_integralScheme,&
       LowdinParameters_schwarzInequality, &
       LowdinParameters_harmonicConstant, &
       
                                !!***************************************************************************
                                !! Parameter to control SCF program
                                !!
       LowdinParameters_nonElectronicEnergyTolerance,&
       LowdinParameters_electronicEnergyTolerance,&
       LowdinParameters_nonelectronicDensityMatrixTolerance,&
       LowdinParameters_electronicDensityMatrixTolerance,&
       LowdinParameters_totalEnergyTolerance,&
       LowdinParameters_totalDensityMatrixTolerance,&
       LowdinParameters_densityFactorThreshold,&
       LowdinParameters_diisSwitchThreshold,&
       LowdinParameters_diisSwitchThreshold_bkp,&
       LowdinParameters_electronicLevelShifting,&
       LowdinParameters_nonelectronicLevelShifting,&
       LowdinParameters_exchangeOrbitalThreshold,&
       LowdinParameters_waveFunctionScale,&
       LowdinParameters_scfNonelectronicMaxIterations,&
       LowdinParameters_scfElectronicMaxIterations,&
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
       LowdinParameters_forceClosedShell,&
       LowdinParameters_debugScfs,&
       LowdinParameters_scfGhostSpecies, &
       
                                !!*****************************************************
                                !! Hartree-Fock Options
                                !!
       LowdinParameters_frozen,&
       LowdinParameters_freezeNonElectronicOrbitals,&
       LowdinParameters_freezeElectronicOrbitals,&
       LowdinParameters_hartreeProductGuess,&
       LowdinParameters_noSCF,&
       LowdinParameters_readCoefficients,&
       LowdinParameters_readFchk,&
       LowdinParameters_readCoefficientsInBinary, &
       LowdinParameters_writeCoefficientsInBinary, &
       LowdinParameters_readEigenvalues,&
       LowdinParameters_readEigenvaluesInBinary, &
       LowdinParameters_writeEigenvaluesInBinary, &
       LowdinParameters_finiteMassCorrection,&
       LowdinParameters_removeTranslationalContamination,&
       LowdinParameters_buildTwoParticlesMatrixForOneParticle,&
       LowdinParameters_buildMixedDensityMatrix,&
       LowdinParameters_onlyElectronicEffect,&
       LowdinParameters_electronicWaveFunctionAnalysis,&
       LowdinParameters_isOpenShell, &
       LowdinParameters_getGradients, &
       LowdinParameters_HFprintEigenvalues, &
       LowdinParameters_HFprintEigenvectors, &
       LowdinParameters_overlapEigenThreshold, &
       LowdinParameters_electricField, &
       
       
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
       LowdinParameters_firstStep,&
       LowdinParameters_lastStep,&
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
                                !! Parameter to control MBPn theory
                                !!
       LowdinParameters_mpCorrection,&
       LowdinParameters_mpFrozenCoreBoundary,&
       LowdinParameters_mpOnlyElectronicCorrection,&
       LowdinParameters_epsteinNesbetCorrection, &
       
                                !!***************************************************************************
                                !! Parameter to control cosmo theory
                                !!
       LowdinParameters_cosmo,& 
       LowdinParameters_cosmo_solvent_dielectric,& 
       LowdinParameters_cosmo_scaling,& 
       
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
       LowdinParameters_ptFactorSS, &
       LowdinParameters_ptFactorOS, &
       LowdinParameters_ptP3Method, &
       
       
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
       LowdinParameters_numberOfCIStates, &
       LowdinParameters_CIdiagonalizationMethod, &
       LowdinParameters_CIdiagonalDressedShift, &
       LowdinParameters_CIactiveSpace, &
       LowdinParameters_CIstatesToPrint, &
       LowdinParameters_CImaxNCV, &
       LowdinParameters_CIsizeOfGuessMatrix, &
       LowdinParameters_CIstackSize, &
       LowdinParameters_CIConvergence, &
       LowdinParameters_CImatvecTolerance, &
       LowdinParameters_CISaveEigenVector, &
       LowdinParameters_CILoadEigenVector, &
       LowdinParameters_CIJacobi, &
       LowdinParameters_CIBuildFullMatrix, &
       LowdinParameters_CIMadSpace, &
       LowdinParameters_CINaturalOrbitals, &
       LowdinParameters_CIPrintEigenVectorsFormat, &
       LowdinParameters_CIPrintThreshold, &
       LowdinParameters_CISCICoreSpace, &
       LowdinParameters_CISCITargetSpace, &


       
                                !!***************************************************************************
                                !! Non-orthogonal CI
                                !!
       LowdinParameters_nonOrthogonalConfigurationInteraction,&
       LowdinParameters_translationScanGrid,&
       LowdinParameters_rotationalScanGrid,&
       LowdinParameters_rotationAroundZMaxAngle,&
       LowdinParameters_rotationAroundZStep,&
       LowdinParameters_nestedRotationalGrids,&
       LowdinParameters_translationStep,&
       LowdinParameters_nestedGridsDisplacement,&
       LowdinParameters_configurationEnergyThreshold,&
       LowdinParameters_configurationOverlapThreshold,&
       LowdinParameters_configurationMinDisplacement,&       
       LowdinParameters_configurationMaxDisplacement,&       
       LowdinParameters_configurationMaxNPDistance,&
       LowdinParameters_configurationMinPPDistance,&
       LowdinParameters_configurationMaxPPDistance,&
       LowdinParameters_configurationEquivalenceDistance,&
       LowdinParameters_empiricalOverlapParameterA,&
       LowdinParameters_empiricalOverlapParameterB,&
       LowdinParameters_empiricalOverlapParameterE0,&
       LowdinParameters_empiricalOverlapParameterSc,&
       LowdinParameters_configurationUseSymmetry,&
       LowdinParameters_readNOCIGeometries,&
       LowdinParameters_empiricalOverlapCorrection,&
       LowdinParameters_onlyFirstNOCIelements,&
       LowdinParameters_NOCIKineticApproximation,&
       LowdinParameters_computeROCIformula,&
       LowdinParameters_removeQDOinCI,&
       !!***************************************************************************
                                !! CCSD 
                                !!
       LowdinParameters_coupledClusterLevel,&
       
                                !!*****************************************************
                                !! Parameter to general control
                                !!
       LowdinParameters_method,&
       LowdinParameters_transformToCenterOfMass,&
       LowdinParameters_areThereDummyAtoms,&
       LowdinParameters_areThereQDOPotentials,&
       LowdinParameters_setQDOEnergyZero, &
       LowdinParameters_isThereExternalPotential,&
       LowdinParameters_isThereInterparticlePotential,&
       LowdinParameters_isThereOutput,&
       LowdinParameters_isThereFrozenParticle,&
       LowdinParameters_dimensionality,&
       
                                !!*****************************************************
                                !! Density Functional Theory Options
                                !!
       LowdinParameters_gridStorage,&
       LowdinParameters_electronCorrelationFunctional,&
       LowdinParameters_electronExchangeFunctional,&
       LowdinParameters_electronExchangeCorrelationFunctional,&
       LowdinParameters_nuclearElectronCorrelationFunctional,&
       LowdinParameters_positronElectronCorrelationFunctional,&
       LowdinParameters_betaFunction,&
       LowdinParameters_gridRadialPoints,&
       LowdinParameters_gridAngularPoints,&
       LowdinParameters_gridNumberOfShells,&
       LowdinParameters_finalGridRadialPoints,&
       LowdinParameters_finalGridAngularPoints,&
       LowdinParameters_finalGridNumberOfShells,&
       LowdinParameters_polarizationOrder,&
       LowdinParameters_numberOfBlocksInAuxiliaryFunctions,&
       LowdinParameters_fukuiFunctions,&
       LowdinParameters_auxiliaryDensity,&
       LowdinParameters_storeThreeCenterElectronIntegrals,&
       LowdinParameters_callLibxc,&
       LowdinParameters_nuclearElectronDensityThreshold,&
       LowdinParameters_electronDensityThreshold,&
       LowdinParameters_gridWeightThreshold,&
       LowdinParameters_betaParameterA,&
       LowdinParameters_betaParameterB,&
       LowdinParameters_betaParameterC,&
       
                                !!*****************************************************
                                !! Subsystem embedding Options
                                !!
       LowdinParameters_subsystemEmbedding,&
       LowdinParameters_localizeOrbitals,&
       LowdinParameters_subsystemLevelShifting,&
       LowdinParameters_subsystemOrbitalThreshold,&
       LowdinParameters_subsystemBasisThreshold,&
       LowdinParameters_erkaleLocalizationMethod,&

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
       LowdinParameters_moldenFileFormat, &
       
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
       LowdinParameters_amberFile,&       
       
                                !!*****************************************************
                                !! Properties Options
       LowdinParameters_calculateInterparticleDistances,&
       LowdinParameters_calculateDensityVolume,&
       
                                !!*****************************************************
                                !! Miscelaneous Options
                                !!
       LowdinParameters_MOFractionOccupation,&
       LowdinParameters_ionizeMO,&
       LowdinParameters_ionizeSpecies,&
       LowdinParameters_exciteSpecies,&
       
                                !!*****************************************************
                                !! Integrals transformation options
                                !!
       LowdinParameters_integralsTransformationMethod, &
       LowdinParameters_ITBuffersize, &
       
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

    integer:: nthreads, proc
    !! Set defaults for namelist

    !!***************************************************************************
    !! Dummy variables, just for debugging. 
    !!
    LowdinParameters_dummyReal(:) = 0.0_8
    LowdinParameters_dummyInteger(:) = 0
    LowdinParameters_dummyLogical(:) = .false.
    LowdinParameters_dummyCharacter(:) = ""

    !!***************************************************************************
    !! Parameter to control Integrals library
    !!  
    LowdinParameters_tv = 1.0E-6
    LowdinParameters_integralThreshold = 1.0E-10
    LowdinParameters_integralStackSize = 30000
    LowdinParameters_integralStorage = "DISK" !! "MEMORY" or "DISK" or "DIRECT"
    LowdinParameters_integralScheme = "LIBINT" !! LIBINT or RYS
    LowdinParameters_schwarzInequality = .false.
    LowdinParameters_harmonicConstant = 0.0_8

    !!***************************************************************************
    !! Parameter to control SCF program
    !!
    LowdinParameters_nonElectronicEnergyTolerance = 1.0E-8
    LowdinParameters_electronicEnergyTolerance =  1.0E-8
    LowdinParameters_nonelectronicDensityMatrixTolerance =  1.0E-6
    LowdinParameters_electronicDensityMatrixTolerance = 1.0E-6
    LowdinParameters_totalEnergyTolerance = 1.0E-8
    LowdinParameters_totalDensityMatrixTolerance = 1.0E-6
    LowdinParameters_densityFactorThreshold = 1.0E-8
    LowdinParameters_diisSwitchThreshold = 0.5
    LowdinParameters_diisSwitchThreshold_bkp = 0.5 
    LowdinParameters_electronicLevelShifting = 0.0
    LowdinParameters_nonelectronicLevelShifting = 0.0
    LowdinParameters_exchangeOrbitalThreshold = 0.8
    LowdinParameters_waveFunctionScale = 1000.0
    LowdinParameters_scfNonelectronicMaxIterations = 50
    LowdinParameters_scfElectronicMaxIterations = 50
    LowdinParameters_scfGlobalMaxIterations = 200
    LowdinParameters_listSize = -20
    LowdinParameters_convergenceMethod = 1 !!(0) NONE, (1) DAMPING, (2) DIIS, (3) LEVEL SHIFTING (4) DAMPING/DIIS
    LowdinParameters_diisDimensionality = 10
    LowdinParameters_iterationScheme = 3 !!(0) NONELECRONIC FULLY / e- (1) ELECTRONIC FULLY (2) CONVERGED INDIVIDIALLY (3) SCHEMESIMULTANEOUS
    LowdinParameters_scfElectronicTypeGuess = "HCORE"
    LowdinParameters_scfNonelectronicTypeGuess = "HCORE"
    LowdinParameters_scfConvergenceCriterium = "ENERGY" !ENERGY, DENSITY, BOTH
    LowdinParameters_diisErrorInDamping = .false.
    LowdinParameters_activateLevelShifting = .false.
    LowdinParameters_exchangeOrbitalsInSCF = .false.
    LowdinParameters_forceClosedShell = .false.
    LowdinParameters_debugScfs = .false.
    LowdinParameters_scfGhostSpecies = "NONE"

    !!*****************************************************
    !! Hartree-Fock Options
    !!
    LowdinParameters_frozen = "NONE"
    LowdinParameters_freezeNonElectronicOrbitals = .false.
    LowdinParameters_freezeElectronicOrbitals = .false.
    LowdinParameters_hartreeProductGuess = .false.
    LowdinParameters_readCoefficients = .true.
    LowdinParameters_readFchk = .false.
    LowdinParameters_writeCoefficientsInBinary = .true.
    LowdinParameters_readEigenvalues = .false.
    LowdinParameters_readEigenvaluesInBinary = .true.
    LowdinParameters_writeEigenvaluesInBinary = .true.
    LowdinParameters_noSCF = .false.
    LowdinParameters_finiteMassCorrection = .false.
    LowdinParameters_removeTranslationalContamination = .false.
    LowdinParameters_buildTwoParticlesMatrixForOneParticle = .false.
    LowdinParameters_buildMixedDensityMatrix = .false.
    LowdinParameters_onlyElectronicEffect = .false.
    LowdinParameters_electronicWaveFunctionAnalysis = .false.
    LowdinParameters_isOpenShell = .false.
    LowdinParameters_getGradients = .false.
    LowdinParameters_HFprintEigenvalues = .true.
    LowdinParameters_HFprintEigenvectors = "OCCUPIED"
    LowdinParameters_overlapEigenThreshold = 1.0E-8_8
    LowdinParameters_electricField(:) = 0.0_8
    LowdinParameters_multipoleOrder = 0

    !!***************************************************************************
    !! Parameter to control geometry optimization
    !!
    LowdinParameters_numericalDerivativeDelta = 1.0E-3
    LowdinParameters_minimizationInitialStepSize = 0.5_8
    LowdinParameters_minimizationLineTolerance = 0.001_8
    LowdinParameters_minimizationToleranceGradient = 0.00001_8
    LowdinParameters_minimizationMaxIteration = 200
    LowdinParameters_minimizationMethod = 4
    LowdinParameters_minimizationLibrary = "GENERIC"
    LowdinParameters_coordinates = "CARTESIAN"
    LowdinParameters_energyCalculator = "INTERNAL"
    LowdinParameters_analyticGradient = .true.
    LowdinParameters_minimizationWithSinglePoint = .true.
    LowdinParameters_useSymmetryInMatrices = .false.
    LowdinParameters_restartOptimization = .false.
    LowdinParameters_firstStep = .true.
    LowdinParameters_lastStep = .true.
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
    !! Parameter to control MBPn theory
    !!
    LowdinParameters_mpCorrection = 1
    LowdinParameters_mpFrozenCoreBoundary = 0
    LowdinParameters_mpOnlyElectronicCorrection = .false.
    LowdinParameters_epsteinNesbetCorrection = 1

    !!***************************************************************************
    !! Parameter to control cosmo theory
    !!
    LowdinParameters_cosmo = .false.
    LowdinParameters_cosmo_solvent_dielectric = 78.3553d+00
    LowdinParameters_cosmo_scaling =0.0d+00

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

    LowdinParameters_ptFactorSS = 0 
    LowdinParameters_ptFactorOS = 0
    LowdinParameters_ptP3Method = "NONE"
    LowdinParameters_ptP3Method(1) = "ALL"

    !!***************************************************************************
    !! Control print level and units
    !!
    LowdinParameters_formatNumberOfColumns = 5
    LowdinParameters_unitForOutputFile = 6
    LowdinParameters_unitForMolecularOrbitalsFile = 8 
    LowdinParameters_unitForMP2IntegralsFile = 7
    LowdinParameters_printLevel =  1 !! (0) no output (1) normal output, (5) method (6) metod and WF (7) method, WF and GLOBAL(8) method, WF, GLOBAL, SCF
    LowdinParameters_units = "ANGS"
    LowdinParameters_doubleZeroThreshold = 1.0E-12

    !!***************************************************************************
    !! CISD - FCI
    !!
    LowdinParameters_configurationInteractionLevel = "NONE"
    LowdinParameters_numberOfCIStates = 1
    LowdinParameters_CIdiagonalizationMethod = "DSYEVR"
    LowdinParameters_CIdiagonalDressedShift = "NONE"
    LowdinParameters_CIactiveSpace = 0 !! Full
    LowdinParameters_CIstatesToPrint = 1
    LowdinParameters_CImaxNCV = 30
    LowdinParameters_CIsizeOfGuessMatrix = 300
    LowdinParameters_CIstackSize = 5000
    LowdinParameters_CIConvergence = 1E-4
    LowdinParameters_CImatvecTolerance = 1E-10
    LowdinParameters_CISaveEigenVector = .false.
    LowdinParameters_CILoadEigenVector = .false.
    LowdinParameters_CIJacobi = .false.
    LowdinParameters_CIBuildFullMatrix = .false. 
    LowdinParameters_CIMadSpace = 5
    LowdinParameters_CINaturalOrbitals = .false.
    LowdinParameters_CIPrintEigenVectorsFormat = "OCCUPIED"
    LowdinParameters_CIPrintThreshold = 1E-1

    !!***************************************************************************
    !! Non-orthogonal CI
    !!
    LowdinParameters_nonOrthogonalConfigurationInteraction=.false.
    LowdinParameters_translationScanGrid(:)=0
    LowdinParameters_rotationalScanGrid=0
    LowdinParameters_rotationAroundZMaxAngle=360.0_8
    LowdinParameters_rotationAroundZStep=0
    LowdinParameters_nestedRotationalGrids=1
    LowdinParameters_translationStep=0.0
    LowdinParameters_nestedGridsDisplacement=0.0
    LowdinParameters_configurationEnergyThreshold=1.0
    LowdinParameters_configurationOverlapThreshold=1.0E-8
    LowdinParameters_configurationMaxDisplacement(:)=0.0
    LowdinParameters_configurationMinDisplacement(:)=0.0
    LowdinParameters_configurationMaxNPDistance=1.0E8
    LowdinParameters_configurationMinPPDistance=0.0
    LowdinParameters_configurationMaxPPDistance=1.0E8
    LowdinParameters_configurationEquivalenceDistance=1.0E-8
    LowdinParameters_empiricalOverlapParameterA=0.0604
    LowdinParameters_empiricalOverlapParameterB=0.492
    LowdinParameters_empiricalOverlapParameterE0=0.0
    LowdinParameters_empiricalOverlapParameterSc=0.0
    LowdinParameters_configurationUseSymmetry=.false.
    LowdinParameters_readNOCIgeometries=.false.
    LowdinParameters_empiricalOverlapCorrection=.false.
    LowdinParameters_onlyFirstNOCIelements=.false.
    LowdinParameters_NOCIKineticApproximation=.false.
    LowdinParameters_computeROCIformula=.false.
    LowdinParameters_removeQDOinCI=.false.
    !!***************************************************************************
    !! CCSD
    !!
    LowdinParameters_coupledClusterLevel = "NONE"

    !!*****************************************************
    !! Parameter to general control
    !!
    LowdinParameters_method = "NONE"
    LowdinParameters_transformToCenterOfMass = .false.
    LowdinParameters_areThereDummyAtoms = .false.
    LowdinParameters_areThereQDOPotentials = .false.
    LowdinParameters_setQDOEnergyZero = .false.
    LowdinParameters_isThereExternalPotential = .false.
    LowdinParameters_isThereInterparticlePotential = .false.
    LowdinParameters_isThereOutput = .false.
    LowdinParameters_isThereFrozenParticle = .false. 
    LowdinParameters_dimensionality = 3

    !!*****************************************************
    !! Density Functional Theory Options
    !!
    LowdinParameters_gridStorage="DISK"
    LowdinParameters_electronCorrelationFunctional = "NONE"
    LowdinParameters_electronExchangeFunctional = "NONE"
    LowdinParameters_electronExchangeCorrelationFunctional = "NONE"
    LowdinParameters_nuclearElectronCorrelationFunctional = "NONE"
    LowdinParameters_positronElectronCorrelationFunctional = "NONE"
    LowdinParameters_betaFunction = "NONE"
    LowdinParameters_gridRadialPoints=35
    LowdinParameters_gridAngularPoints=110
    LowdinParameters_gridNumberOfShells=5
    LowdinParameters_finalGridRadialPoints=50
    LowdinParameters_finalGridAngularPoints=302
    LowdinParameters_finalGridNumberOfShells=5
    LowdinParameters_polarizationOrder = 1
    LowdinParameters_numberOfBlocksInAuxiliaryFunctions = 3
    LowdinParameters_fukuiFunctions = .false.
    LowdinParameters_auxiliaryDensity = .false.
    LowdinParameters_storeThreeCenterElectronIntegrals = .true.
    LowdinParameters_callLibxc = .true.
    LowdinParameters_nuclearElectronDensityThreshold = 1E-10
    LowdinParameters_electronDensityThreshold = 1E-10
    LowdinParameters_gridWeightThreshold = 1E-10
    LowdinParameters_betaParameterA = 0.0
    LowdinParameters_betaParameterB = 0.0
    LowdinParameters_betaParameterC = 0.0

    !!*****************************************************
    !! Subsystem embedding Options
    !!
    LowdinParameters_subsystemEmbedding = .false.
    LowdinParameters_localizeOrbitals = .false.
    LowdinParameters_subsystemLevelShifting = 1.0E6
    LowdinParameters_subsystemOrbitalThreshold = 0.1
    LowdinParameters_subsystemBasisThreshold = 0.0001
    LowdinParameters_erkaleLocalizationMethod = "MU"

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
    LowdinParameters_moldenFileFormat = "MIXED" 

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
    LowdinParameters_amberFile = .false.    

    !!*****************************************************
    !! Properties Options
    LowdinParameters_calculateInterparticleDistances = .false.
    LowdinParameters_calculateDensityVolume = .false.

    !!*****************************************************
    !! Miscelaneous Options
    !!
    LowdinParameters_MOFractionOccupation = 1.0_8
    LowdinParameters_ionizeMO = 0
    LowdinParameters_ionizeSpecies = "NONE"
    LowdinParameters_exciteSpecies = "NONE"

    !!*****************************************************
    !! Integrals transformation options
    !!
    LowdinParameters_integralsTransformationMethod = "C"
    LowdinParameters_ITBuffersize = 1024

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
    !! Dummy variables, just for debugging. 
    !!
    CONTROL_instance%DUMMY_REAL(:) = 0 
    CONTROL_instance%DUMMY_INTEGER(:) = 0
    CONTROL_instance%DUMMY_LOGICAL(:) = .false.
    CONTROL_instance%DUMMY_CHARACTER(:) = ""

    !!***************************************************************************    
    !! Parameter to control Integrals library                       
    !!
    CONTROL_instance%TV = 1.0E-6
    CONTROL_instance%INTEGRAL_THRESHOLD = 1.0E-10
    CONTROL_instance%INTEGRAL_STACK_SIZE = 30000
    CONTROL_instance%INTEGRAL_STORAGE = "DISK" !! "DISK" or "DIRECT"
    CONTROL_instance%INTEGRAL_SCHEME = "LIBINT" !! LIBINT or Rys
    CONTROL_instance%SCHWARZ_INEQUALITY = .false.
    CONTROL_instance%HARMONIC_CONSTANT = 0.0_8

    !!***************************************************************************
    !! Parameter to control SCF program
    !!
    CONTROL_instance%NONELECTRONIC_ENERGY_TOLERANCE = 1.0E-8
    CONTROL_instance%ELECTRONIC_ENERGY_TOLERANCE =  1.0E-8
    CONTROL_instance%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE =  1.0E-6
    CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE = 1.0E-6
    CONTROL_instance%TOTAL_ENERGY_TOLERANCE = 1.0E-8
    CONTROL_instance%TOTAL_DENSITY_MATRIX_TOLERANCE = 1.0E-6
    CONTROL_instance%DENSITY_FACTOR_THRESHOLD = 1.0E-8
    CONTROL_instance%DIIS_SWITCH_THRESHOLD = 0.5
    CONTROL_instance%DIIS_SWITCH_THRESHOLD_BKP = 0.5 
    CONTROL_instance%ELECTRONIC_LEVEL_SHIFTING = 0.0
    CONTROL_instance%NONELECTRONIC_LEVEL_SHIFTING = 0.0
    CONTROL_instance%EXCHANGE_ORBITAL_THRESHOLD = 0.8
    CONTROL_instance%WAVE_FUNCTION_SCALE = 1000.0
    CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS = 50
    CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS = 50
    CONTROL_instance%SCF_GLOBAL_MAX_ITERATIONS = 200 
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
    CONTROL_instance%FORCE_CLOSED_SHELL = .true.
    CONTROL_instance%DEBUG_SCFS = .false.
    CONTROL_instance%SCF_GHOST_SPECIES = "NONE"
    ! CONTROL_instance%DEBUG_SCFS = .true.

    !!***************************************************************************                                              
    !! Hartree-Fock options                                                                                                    
    !!                                                                                                                         
    CONTROL_instance%FROZEN_PARTICLE = "NONE"
    CONTROL_instance%FREEZE_NON_ELECTRONIC_ORBITALS = .false.
    CONTROL_instance%FREEZE_ELECTRONIC_ORBITALS = .false.
    CONTROL_instance%HARTREE_PRODUCT_GUESS = .false.
    CONTROL_instance%READ_COEFFICIENTS = .true.
    CONTROL_instance%READ_FCHK=.false.
    CONTROL_instance%WRITE_COEFFICIENTS_IN_BINARY = .true.
    CONTROL_instance%NO_SCF = .false.
    CONTROL_instance%FINITE_MASS_CORRECTION = .false.
    CONTROL_instance%REMOVE_TRANSLATIONAL_CONTAMINATION = .false.
    CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE = .false.
    CONTROL_instance%BUILD_MIXED_DENSITY_MATRIX = .false.
    CONTROL_instance%ONLY_ELECTRONIC_EFFECT = .false.
    CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS = .false.
    CONTROL_instance%IS_OPEN_SHELL = .false.
    CONTROL_instance%GET_GRADIENTS = .false.
    CONTROL_instance%HF_PRINT_EIGENVALUES = .true.
    CONTROL_instance%HF_PRINT_EIGENVECTORS = "OCCUPIED"
    CONTROL_instance%OVERLAP_EIGEN_THRESHOLD = 1.0E-8_8
    CONTROL_instance%ELECTRIC_FIELD(:) = 0.0_8
    CONTROL_instance%MULTIPOLE_ORDER = 0
    !!***************************************************************************                                              
    !! Parameter to control geometry optimization                                                                              
    !!                                                                                                                         
    CONTROL_instance%NUMERICAL_DERIVATIVE_DELTA = 1.0E-3
    CONTROL_instance%MINIMIZATION_INITIAL_STEP_SIZE = 0.1_8
    CONTROL_instance%MINIMIZATION_LINE_TOLERANCE = 0.1_8
    CONTROL_instance%MINIMIZATION_TOLERANCE_GRADIENT = 1.0E-5
    CONTROL_instance%MINIMIZATION_MAX_ITERATION = 100
    CONTROL_instance%MINIMIZATION_METHOD = 4
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
    CONTROL_instance%LAST_STEP = .true.
    CONTROL_instance%OPTIMIZE_WITH_MP = .false.
    CONTROL_instance%PROJECT_HESSIANE = .true.

    !!***************************************************************************                                              
    !! Parameter of atomic conectivity                                                                                         
    !!                                                                                                                         
    CONTROL_instance%BOND_DISTANCE_FACTOR = 1.3_8
    CONTROL_instance%BOND_ANGLE_THRESHOLD = 180.0_8
    CONTROL_instance%DIHEDRAL_ANGLE_THRESHOLD = 180.0_8

    !!***************************************************************************                                              
    !! Parameter to control MBPn theory                                                                                         
    !!                                                                                                                         
    CONTROL_instance%MOLLER_PLESSET_CORRECTION = 1
    CONTROL_instance%MP_FROZEN_CORE_BOUNDARY = 0
    CONTROL_instance%MP_ONLY_ELECTRONIC_CORRECTION = .false.

    CONTROL_instance%EPSTEIN_NESBET_CORRECTION = 1

    !!***************************************************************************                                              
    !! Parameter to control cosmo method                                                                                         
    !!                                                                                                                         
    CONTROL_instance%COSMO = .false.
    CONTROL_instance%COSMO_SOLVENT_DIELECTRIC= 78.3553d+00 
    CONTROL_instance%COSMO_SCALING= 0.0d+00

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
    CONTROL_instance%PT_FACTOR_SS = 0
    CONTROL_instance%PT_FACTOR_OS = 0 
    CONTROL_instance%PT_P3_METHOD = "NONE"
    CONTROL_instance%PT_P3_METHOD(1) = "ALL"

    !!***************************************************************************                                              
    !! Control print level and units                                                                                           
    !!                                                                                                                         
    CONTROL_instance%FORMAT_NUMBER_OF_COLUMNS = 5
    CONTROL_instance%UNIT_FOR_OUTPUT_FILE = 6
    CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE = 8 
    CONTROL_instance%UNIT_FOR_MOLECULAR_ORBITALS_FILE = 7
    CONTROL_instance%PRINT_LEVEL =  1 !! (0) no output, (1) normal output, (5) method (6) metod and WF (7) method, WF and GLOBAL(8) method, WF, GLOBAL, SCF
    CONTROL_instance%UNITS = "ANGS"
    CONTROL_instance%DOUBLE_ZERO_THRESHOLD = 1.0E-12

    !!***************************************************************************                                              
    !! CISD - FCI                                                                                                              
    !!                                                                                                                         
    CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL = "NONE"
    CONTROL_instance%NUMBER_OF_CI_STATES= 1
    CONTROL_instance%CI_DIAGONALIZATION_METHOD = "DSYEVR"
    CONTROL_instance%CI_DIAGONAL_DRESSED_SHIFT = "NONE"
    CONTROL_instance%CI_ACTIVE_SPACE = 0 !! Full
    CONTROL_instance%CI_STATES_TO_PRINT = 1
    CONTROL_instance%CI_MAX_NCV = 30 
    CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX = 300
    CONTROL_instance%CI_STACK_SIZE = 5000
    CONTROL_instance%CI_CONVERGENCE = 1E-4
    CONTROL_instance%CI_MATVEC_TOLERANCE = 1E-10
    CONTROL_instance%CI_SAVE_EIGENVECTOR = .FALSE.
    CONTROL_instance%CI_LOAD_EIGENVECTOR = .FALSE.
    CONTROL_instance%CI_JACOBI = .False.
    CONTROL_instance%CI_BUILD_FULL_MATRIX = .FALSE. 
    CONTROL_instance%CI_MADSPACE = 5
    CONTROL_instance%CI_NATURAL_ORBITALS=.FALSE.
    CONTROL_instance%CI_PRINT_EIGENVECTORS_FORMAT = "OCCUPIED"
    CONTROL_instance%CI_PRINT_THRESHOLD = 1E-1
    CONTROL_instance%CI_SCI_CORE_SPACE = 100
    CONTROL_instance%CI_SCI_TARGET_SPACE = 10000

    !!***************************************************************************
    !! Non-orthogonal CI
    !!
    CONTROL_instance%NONORTHOGONAL_CONFIGURATION_INTERACTION=.FALSE.
    CONTROL_instance%TRANSLATION_SCAN_GRID(:)=0
    CONTROL_instance%ROTATIONAL_SCAN_GRID=0
    CONTROL_instance%ROTATION_AROUND_Z_MAX_ANGLE=360
    CONTROL_instance%ROTATION_AROUND_Z_STEP=0
    CONTROL_instance%NESTED_ROTATIONAL_GRIDS=1
    CONTROL_instance%TRANSLATION_STEP=0.0
    CONTROL_instance%NESTED_GRIDS_DISPLACEMENT=0.0
    CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD=1.0
    CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD=1.0E-8
    CONTROL_instance%CONFIGURATION_MAX_DISPLACEMENT(:)=0.0
    CONTROL_instance%CONFIGURATION_MIN_DISPLACEMENT(:)=0.0
    CONTROL_instance%CONFIGURATION_MAX_NP_DISTANCE=1.0E8
    CONTROL_instance%CONFIGURATION_MIN_PP_DISTANCE=0.0
    CONTROL_instance%CONFIGURATION_MAX_PP_DISTANCE=1.0E8
    CONTROL_instance%CONFIGURATION_EQUIVALENCE_DISTANCE=1.0E-8
    CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_A=0.0604
    CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_B=0.492
    CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_E0=0.0
    CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_Sc=0.0
    CONTROL_instance%CONFIGURATION_USE_SYMMETRY=.false.
    CONTROL_instance%READ_NOCI_GEOMETRIES=.false.
    CONTROL_instance%EMPIRICAL_OVERLAP_CORRECTION=.false.
    CONTROL_instance%ONLY_FIRST_NOCI_ELEMENTS=.false.
    CONTROL_instance%NOCI_KINETIC_APPROXIMATION=.false.
    CONTROL_instance%COMPUTE_ROCI_FORMULA=.false.
    CONTROL_instance%REMOVE_QDO_IN_CI=.false.
    !!***************************************************************************                                              
    !! CCSD                                                                                                              
    !!                                                                                                                         
    CONTROL_instance%COUPLED_CLUSTER_LEVEL = "NONE"

    !!*****************************************************                                                                    
    !! Parameter to general control                                                                                            
    !!                                                                                                                         
    CONTROL_instance%METHOD = "NONE"
    CONTROL_instance%TRANSFORM_TO_CENTER_OF_MASS = .false.
    CONTROL_instance%ARE_THERE_DUMMY_ATOMS = .false.
    CONTROL_instance%ARE_THERE_QDO_POTENTIALS = .false.
    CONTROL_instance%SET_QDO_ENERGY_ZERO = .false.
    CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL = .false.
    CONTROL_instance%IS_THERE_INTERPARTICLE_POTENTIAL = .false.
    CONTROL_instance%IS_THERE_OUTPUT = .false.
    CONTROL_instance%IS_THERE_FROZEN_PARTICLE = .false. 
    CONTROL_instance%DIMENSIONALITY = 3

    !!*****************************************************                                                                    
    !! Density Functional Theory Options                                                                                       
    !!                                                                                                                         
    CONTROL_instance%GRID_STORAGE="DISK"
    CONTROL_instance%ELECTRON_CORRELATION_FUNCTIONAL = "NONE"
    CONTROL_instance%ELECTRON_EXCHANGE_FUNCTIONAL = "NONE"
    CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL = "NONE"
    CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL = "NONE"
    CONTROL_instance%POSITRON_ELECTRON_CORRELATION_FUNCTIONAL = "NONE"
    CONTROL_instance%BETA_FUNCTION = "NONE"
    CONTROL_instance%GRID_RADIAL_POINTS= 35
    CONTROL_instance%GRID_ANGULAR_POINTS= 110
    CONTROL_instance%GRID_NUMBER_OF_SHELLS= 5
    CONTROL_instance%FINAL_GRID_RADIAL_POINTS= 50
    CONTROL_instance%FINAL_GRID_ANGULAR_POINTS= 302
    CONTROL_instance%FINAL_GRID_NUMBER_OF_SHELLS= 5
    CONTROL_instance%POLARIZATION_ORDER = 1
    CONTROL_instance%NUMBER_OF_BLOCKS_IN_AUXILIARY_FUNCTIONS = 3
    CONTROL_instance%FUKUI_FUNCTIONS = .false.
    CONTROL_instance%AUXILIARY_DENSITY = .false.
    CONTROL_instance%STORE_THREE_CENTER_ELECTRON_INTEGRALS = .true.
    CONTROL_instance%CALL_LIBXC = .true.
    CONTROL_instance%NUCLEAR_ELECTRON_DENSITY_THRESHOLD = 1E-10
    CONTROL_instance%ELECTRON_DENSITY_THRESHOLD = 1E-10
    CONTROL_instance%GRID_WEIGHT_THRESHOLD = 1E-10
    CONTROL_instance%BETA_PARAMETER_A=0.0
    CONTROL_instance%BETA_PARAMETER_B=0.0
    CONTROL_instance%BETA_PARAMETER_C=0.0

    !!*****************************************************
    !! Subsystem embedding Options
    !!
    CONTROL_instance%SUBSYSTEM_EMBEDDING = .false.
    CONTROL_instance%LOCALIZE_ORBITALS = .false.
    CONTROL_instance%SUBSYSTEM_LEVEL_SHIFTING = 1.0E6
    CONTROL_instance%SUBSYSTEM_ORBITAL_THRESHOLD = 0.1
    CONTROL_instance%SUBSYSTEM_BASIS_THRESHOLD = 0.0001
    CONTROL_instance%ERKALE_LOCALIZATION_METHOD = "MU"

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
    LowdinParameters_moldenFileFormat = "MIXED" 

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
    CONTROL_instance%AMBER_FILE = .false.    

    !!*****************************************************                                                                    
    !! Properties Options                                                                                                      
    CONTROL_instance%CALCULATE_INTERPARTICLE_DISTANCES = .false.
    CONTROL_instance%CALCULATE_DENSITY_VOLUME = .false.

    !!*****************************************************                                                                    
    !! Miscelaneous Options                                                                                                    
    !!                                                                                                                         
    CONTROL_instance%MO_FRACTION_OCCUPATION = 1.0_8
    CONTROL_instance%IONIZE_MO = 0
    CONTROL_instance%IONIZE_SPECIES = "NONE"
    CONTROL_instance%EXCITE_SPECIES = "NONE"                                                            
    !$OMP PARALLEL private(nthreads, proc)
    proc = OMP_GET_THREAD_NUM()
    if(proc == 0) then
      nthreads = OMP_GET_NUM_THREADS()
      CONTROL_instance%NUMBER_OF_CORES=nthreads
    end if
    !$OMP END PARALLEL
    
    
    !!*****************************************************
    !! Integrals transformation options
    !!
    CONTROL_instance%INTEGRALS_TRANSFORMATION_METHOD = "C"
    CONTROL_instance%IT_BUFFERSIZE = 8192

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
    integer :: stat
    character(1000) :: line
    
    uunit = 4
    if(present(unit)) uunit = unit

    !! Reload file
    rewind(uunit)

    !! Reads name-list
    read(uunit,NML=LowdinParameters, iostat=stat)

    !! Check the process
    if(stat > 0 ) then

       write (*,'(A)') 'Error reading LowdinParameters'
       ! write (*, '(a, i0)') 'iostat was:', stat
       backspace(uunit)
       read(uunit,fmt='(A)') line
       write(*,'(A)') 'Invalid line : '//trim(line)
         
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

    !!***************************************************************************
    !! Dummy variables, just for debugging. 
    !!
    CONTROL_instance%DUMMY_REAL(:) = LowdinParameters_dummyReal(:)
    CONTROL_instance%DUMMY_INTEGER(:) = LowdinParameters_dummyInteger(:)
    CONTROL_instance%DUMMY_LOGICAL(:) = LowdinParameters_dummyLogical(:)
    CONTROL_instance%DUMMY_CHARACTER(:) = LowdinParameters_dummyCharacter(:)


    !!***************************************************************************      
    !! Parameter to control Integrals library                                          
    !!                                                                                 
    CONTROL_instance%TV = LowdinParameters_tv
    CONTROL_instance%INTEGRAL_THRESHOLD = LowdinParameters_integralThreshold
    CONTROL_instance%INTEGRAL_STACK_SIZE = LowdinParameters_integralStackSize
    CONTROL_instance%INTEGRAL_STORAGE = LowdinParameters_integralStorage
    CONTROL_instance%INTEGRAL_SCHEME =  LowdinParameters_integralScheme
    CONTROL_instance%SCHWARZ_INEQUALITY = LowdinParameters_schwarzInequality
    CONTROL_instance%HARMONIC_CONSTANT = LowdinParameters_harmonicConstant


    !!***************************************************************************      
    !! Parameter to control SCF program                                                
    !!                                                                                 
    CONTROL_instance%NONELECTRONIC_ENERGY_TOLERANCE = LowdinParameters_nonElectronicEnergyTolerance
    CONTROL_instance%ELECTRONIC_ENERGY_TOLERANCE = LowdinParameters_electronicEnergyTolerance
    CONTROL_instance%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE = LowdinParameters_nonelectronicDensityMatrixTolerance
    CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE = LowdinParameters_electronicDensityMatrixTolerance
    CONTROL_instance%TOTAL_ENERGY_TOLERANCE = LowdinParameters_totalEnergyTolerance
    CONTROL_instance%TOTAL_DENSITY_MATRIX_TOLERANCE = LowdinParameters_totalDensityMatrixTolerance
    CONTROL_instance%DENSITY_FACTOR_THRESHOLD = LowdinParameters_densityFactorThreshold
    CONTROL_instance%DIIS_SWITCH_THRESHOLD = LowdinParameters_diisSwitchThreshold
    CONTROL_instance%DIIS_SWITCH_THRESHOLD_BKP = LowdinParameters_diisSwitchThreshold_bkp
    CONTROL_instance%ELECTRONIC_LEVEL_SHIFTING = LowdinParameters_electronicLevelShifting
    CONTROL_instance%NONELECTRONIC_LEVEL_SHIFTING = LowdinParameters_nonelectronicLevelShifting
    CONTROL_instance%EXCHANGE_ORBITAL_THRESHOLD = LowdinParameters_exchangeOrbitalThreshold
    CONTROL_instance%WAVE_FUNCTION_SCALE = LowdinParameters_waveFunctionScale
    CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS = LowdinParameters_scfNonelectronicMaxIterations
    CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS = LowdinParameters_scfElectronicMaxIterations
    CONTROL_instance%SCF_GLOBAL_MAX_ITERATIONS = LowdinParameters_scfGlobalMaxIterations
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
    CONTROL_instance%FORCE_CLOSED_SHELL = LowdinParameters_forceClosedShell
    CONTROL_instance%DEBUG_SCFS = LowdinParameters_debugScfs
    CONTROL_instance%SCF_GHOST_SPECIES = LowdinParameters_scfGhostSpecies

    !!*****************************************************                            
    !! Hartree-Fock Options                                                            
    !!                                                                                 
    CONTROL_instance%FROZEN_PARTICLE = LowdinParameters_frozen
    CONTROL_instance%FREEZE_NON_ELECTRONIC_ORBITALS = LowdinParameters_freezeNonElectronicOrbitals
    CONTROL_instance%FREEZE_ELECTRONIC_ORBITALS = LowdinParameters_freezeElectronicOrbitals
    CONTROL_instance%HARTREE_PRODUCT_GUESS = LowdinParameters_hartreeProductGuess
    CONTROL_instance%READ_COEFFICIENTS = LowdinParameters_readCoefficients
    CONTROL_instance%READ_FCHK = LowdinParameters_readFchk
    CONTROL_instance%WRITE_COEFFICIENTS_IN_BINARY = LowdinParameters_writeCoefficientsInBinary
    CONTROL_instance%READ_EIGENVALUES = LowdinParameters_readEigenvalues
    CONTROL_instance%READ_EIGENVALUES_IN_BINARY =  LowdinParameters_readEigenvaluesInBinary
    CONTROL_instance%WRITE_EIGENVALUES_IN_BINARY = LowdinParameters_writeEigenvaluesInBinary
    CONTROL_instance%NO_SCF = LowdinParameters_noSCF
    CONTROL_instance%FINITE_MASS_CORRECTION = LowdinParameters_finiteMassCorrection
    CONTROL_instance%REMOVE_TRANSLATIONAL_CONTAMINATION = LowdinParameters_removeTranslationalContamination
    CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE = LowdinParameters_buildTwoParticlesMatrixForOneParticle
    CONTROL_instance%BUILD_MIXED_DENSITY_MATRIX = LowdinParameters_buildMixedDensityMatrix
    CONTROL_instance%ONLY_ELECTRONIC_EFFECT = LowdinParameters_onlyElectronicEffect
    CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS = LowdinParameters_electronicWaveFunctionAnalysis
    CONTROL_instance%IS_OPEN_SHELL = LowdinParameters_isOpenShell
    CONTROL_instance%GET_GRADIENTS = LowdinParameters_getGradients
    CONTROL_instance%HF_PRINT_EIGENVALUES = LowdinParameters_HFprintEigenvalues
    CONTROL_instance%HF_PRINT_EIGENVECTORS = LowdinParameters_HFprintEigenvectors
    CONTROL_instance%OVERLAP_EIGEN_THRESHOLD = LowdinParameters_overlapEigenThreshold 

    CONTROL_instance%ELECTRIC_FIELD = LowdinParameters_electricField
    CONTROL_instance%MULTIPOLE_ORDER = LowdinParameters_multipoleOrder
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
    CONTROL_instance%FIRST_STEP = LowdinParameters_firstStep
    CONTROL_instance%LAST_STEP = LowdinParameters_lastStep
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
    !! Parameter to control MBPn theory                                                 
    !!                                                                                 
    CONTROL_instance%MOLLER_PLESSET_CORRECTION = LowdinParameters_mpCorrection
    CONTROL_instance%MP_FROZEN_CORE_BOUNDARY = LowdinParameters_mpFrozenCoreBoundary
    CONTROL_instance%MP_ONLY_ELECTRONIC_CORRECTION = LowdinParameters_mpOnlyElectronicCorrection
    CONTROL_instance%EPSTEIN_NESBET_CORRECTION = LowdinParameters_epsteinNesbetCorrection

    !!***************************************************************************      
    !! Parameter to control cosmo method                                               
    !!                                                                                 
    CONTROL_instance%COSMO = LowdinParameters_cosmo
    CONTROL_instance%COSMO_SOLVENT_DIELECTRIC= LowdinParameters_cosmo_solvent_dielectric
    CONTROL_instance%COSMO_SCALING=LowdinParameters_cosmo_SCALING

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
    CONTROL_instance%PT_FACTOR_SS = LowdinParameters_ptFactorSS
    CONTROL_instance%PT_FACTOR_OS = LowdinParameters_ptFactorOS
    CONTROL_instance%PT_P3_METHOD = LowdinParameters_ptP3Method


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
    CONTROL_instance%NUMBER_OF_CI_STATES       = LowdinParameters_numberOfCIStates
    CONTROL_instance%CI_DIAGONALIZATION_METHOD = LowdinParameters_CIdiagonalizationMethod
    CONTROL_instance%CI_DIAGONAL_DRESSED_SHIFT = LowdinParameters_CIdiagonalDressedShift
    CONTROL_instance%CI_ACTIVE_SPACE = LowdinParameters_CIactiveSpace  
    CONTROL_instance%CI_STATES_TO_PRINT = LowdinParameters_CIstatesToPrint
    if(CONTROL_instance%CI_STATES_TO_PRINT .gt. CONTROL_instance%NUMBER_OF_CI_STATES) &
         CONTROL_instance%NUMBER_OF_CI_STATES=CONTROL_instance%CI_STATES_TO_PRINT
    CONTROL_instance%CI_MAX_NCV = LowdinParameters_CImaxNCV
    CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX = LowdinParameters_CIsizeOfGuessMatrix
    CONTROL_instance%CI_STACK_SIZE = LowdinParameters_CIstackSize
    CONTROL_instance%CI_CONVERGENCE = LowdinParameters_CIConvergence
    CONTROL_instance%CI_MATVEC_TOLERANCE = LowdinParameters_CIMatvecTolerance
    CONTROL_instance%CI_SAVE_EIGENVECTOR = LowdinParameters_CISaveEigenVector
    CONTROL_instance%CI_LOAD_EIGENVECTOR = LowdinParameters_CILoadEigenVector
    CONTROL_instance%CI_JACOBI = LowdinParameters_CIJacobi
    CONTROL_instance%CI_BUILD_FULL_MATRIX = LowdinParameters_CIBuildFullMatrix 
    CONTROL_instance%CI_MADSPACE = LowdinParameters_CIMadSpace
    CONTROL_instance%CI_NATURAL_ORBITALS= LowdinParameters_CINaturalOrbitals
    CONTROL_instance%CI_PRINT_EIGENVECTORS_FORMAT = LowdinParameters_CIPrintEigenVectorsFormat 
    CONTROL_instance%CI_PRINT_THRESHOLD = LowdinParameters_CIPrintThreshold 
    CONTROL_instance%CI_SCI_CORE_SPACE = LowdinParameters_CISCICoreSpace
    CONTROL_instance%CI_SCI_TARGET_SPACE = LowdinParameters_CISCITargetSpace



    !!***************************************************************************
    !! Non-orthogonal CI
    !!
    CONTROL_instance%NONORTHOGONAL_CONFIGURATION_INTERACTION=LowdinParameters_nonOrthogonalConfigurationInteraction
    CONTROL_instance%TRANSLATION_SCAN_GRID=LowdinParameters_translationScanGrid
    CONTROL_instance%ROTATIONAL_SCAN_GRID=LowdinParameters_rotationalScanGrid
    CONTROL_instance%ROTATION_AROUND_Z_MAX_ANGLE=LowdinParameters_rotationAroundZMaxAngle
    CONTROL_instance%ROTATION_AROUND_Z_STEP=LowdinParameters_rotationAroundZStep
    CONTROL_instance%NESTED_ROTATIONAL_GRIDS=LowdinParameters_nestedRotationalGrids
    CONTROL_instance%TRANSLATION_STEP=LowdinParameters_translationStep
    CONTROL_instance%NESTED_GRIDS_DISPLACEMENT=LowdinParameters_nestedGridsDisplacement
    CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD=LowdinParameters_configurationEnergyThreshold
    CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD=LowdinParameters_configurationOverlapThreshold
    CONTROL_instance%CONFIGURATION_MIN_DISPLACEMENT=LowdinParameters_configurationMinDisplacement
    if(sum(LowdinParameters_configurationMaxDisplacement).ne.0.0) then
       CONTROL_instance%CONFIGURATION_MAX_DISPLACEMENT=LowdinParameters_configurationMaxDisplacement
    else
       CONTROL_instance%CONFIGURATION_MAX_DISPLACEMENT=1.001*CONTROL_instance%TRANSLATION_STEP*CONTROL_instance%TRANSLATION_SCAN_GRID/2
    end if
    CONTROL_instance%CONFIGURATION_MAX_NP_DISTANCE=LowdinParameters_configurationMaxNPDistance
    CONTROL_instance%CONFIGURATION_MIN_PP_DISTANCE=LowdinParameters_configurationMinPPDistance
    CONTROL_instance%CONFIGURATION_MAX_PP_DISTANCE=LowdinParameters_configurationMaxPPDistance
    CONTROL_instance%CONFIGURATION_EQUIVALENCE_DISTANCE=LowdinParameters_configurationEquivalenceDistance
    CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_A=LowdinParameters_empiricalOverlapParameterA
    CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_B=LowdinParameters_empiricalOverlapParameterB
    CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_E0=LowdinParameters_empiricalOverlapParameterE0
    CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_SC=LowdinParameters_empiricalOverlapParameterSc
    CONTROL_instance%CONFIGURATION_USE_SYMMETRY=LowdinParameters_configurationUseSymmetry
    CONTROL_instance%READ_NOCI_GEOMETRIES=LowdinParameters_readNOCIGeometries
    CONTROL_instance%EMPIRICAL_OVERLAP_CORRECTION=LowdinParameters_empiricalOverlapCorrection
    CONTROL_instance%ONLY_FIRST_NOCI_ELEMENTS=LowdinParameters_onlyFirstNOCIelements
    CONTROL_instance%NOCI_KINETIC_APPROXIMATION=LowdinParameters_NOCIKineticApproximation
    CONTROL_instance%COMPUTE_ROCI_FORMULA=LowdinParameters_computeROCIformula
    CONTROL_instance%REMOVE_QDO_IN_CI=LowdinParameters_removeQDOinCI

    !!***************************************************************************      
    !! CCSD                                                                       
    !!                                                                                 
    CONTROL_instance%COUPLED_CLUSTER_LEVEL = LowdinParameters_coupledClusterLevel

    !!*****************************************************                            
    !! Parameter to general control                                                    
    !!                                                                                 
    CONTROL_instance%METHOD = LowdinParameters_method
    CONTROL_instance%TRANSFORM_TO_CENTER_OF_MASS = LowdinParameters_transformToCenterOfMass
    CONTROL_instance%ARE_THERE_DUMMY_ATOMS = LowdinParameters_areThereDummyAtoms
    CONTROL_instance%ARE_THERE_QDO_POTENTIALS = LowdinParameters_areThereQDOPotentials
    CONTROL_instance%SET_QDO_ENERGY_ZERO = LowdinParameters_setQDOEnergyZero
    CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL = LowdinParameters_isThereExternalPotential
    CONTROL_instance%IS_THERE_INTERPARTICLE_POTENTIAL = LowdinParameters_isThereInterparticlePotential
    CONTROL_instance%IS_THERE_OUTPUT = LowdinParameters_isThereOutput
    CONTROL_instance%IS_THERE_FROZEN_PARTICLE = LowdinParameters_isThereFrozenParticle
    CONTROL_instance%DIMENSIONALITY = LowdinParameters_dimensionality

    !!*****************************************************                            
    !! Density Functional Theory Options                                               
    !!                                                                                 
    CONTROL_instance%GRID_STORAGE=LowdinParameters_gridStorage
    CONTROL_instance%ELECTRON_CORRELATION_FUNCTIONAL = LowdinParameters_electronCorrelationFunctional
    CONTROL_instance%ELECTRON_EXCHANGE_FUNCTIONAL = LowdinParameters_electronExchangeFunctional
    CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL = LowdinParameters_electronExchangeCorrelationFunctional
    CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL = LowdinParameters_nuclearElectronCorrelationFunctional
    CONTROL_instance%POSITRON_ELECTRON_CORRELATION_FUNCTIONAL = LowdinParameters_positronElectronCorrelationFunctional
    CONTROL_instance%BETA_FUNCTION = LowdinParameters_betaFunction
    CONTROL_instance%GRID_RADIAL_POINTS= LowdinParameters_gridRadialPoints
    CONTROL_instance%GRID_ANGULAR_POINTS= LowdinParameters_gridAngularPoints
    CONTROL_instance%GRID_NUMBER_OF_SHELLS= LowdinParameters_gridNumberOfShells
    if(LowdinParameters_finalGridRadialPoints*LowdinParameters_finalGridAngularPoints .gt. LowdinParameters_gridRadialPoints*LowdinParameters_gridAngularPoints) then
       CONTROL_instance%FINAL_GRID_RADIAL_POINTS= LowdinParameters_finalGridRadialPoints
       CONTROL_instance%FINAL_GRID_ANGULAR_POINTS= LowdinParameters_finalGridAngularPoints
       CONTROL_instance%FINAL_GRID_NUMBER_OF_SHELLS= LowdinParameters_finalGridNumberOfShells
    else
       CONTROL_instance%FINAL_GRID_RADIAL_POINTS= LowdinParameters_gridRadialPoints
       CONTROL_instance%FINAL_GRID_ANGULAR_POINTS= LowdinParameters_gridAngularPoints
       CONTROL_instance%FINAL_GRID_NUMBER_OF_SHELLS= LowdinParameters_gridNumberOfShells
    end if

    CONTROL_instance%POLARIZATION_ORDER = LowdinParameters_polarizationOrder
    CONTROL_instance%NUMBER_OF_BLOCKS_IN_AUXILIARY_FUNCTIONS = LowdinParameters_numberOfBlocksInAuxiliaryFunctions
    CONTROL_instance%FUKUI_FUNCTIONS = LowdinParameters_fukuiFunctions
    CONTROL_instance%AUXILIARY_DENSITY = LowdinParameters_auxiliaryDensity
    CONTROL_instance%STORE_THREE_CENTER_ELECTRON_INTEGRALS = LowdinParameters_storeThreeCenterElectronIntegrals
    CONTROL_instance%CALL_LIBXC = LowdinParameters_callLibxc
    CONTROL_instance%NUCLEAR_ELECTRON_DENSITY_THRESHOLD = LowdinParameters_nuclearElectronDensityThreshold
    CONTROL_instance%ELECTRON_DENSITY_THRESHOLD = LowdinParameters_electronDensityThreshold
    CONTROL_instance%GRID_WEIGHT_THRESHOLD = LowdinParameters_gridWeightThreshold
    CONTROL_instance%BETA_PARAMETER_A = LowdinParameters_betaParameterA
    CONTROL_instance%BETA_PARAMETER_B = LowdinParameters_betaParameterB
    CONTROL_instance%BETA_PARAMETER_C = LowdinParameters_betaParameterC

    !!*****************************************************
    !! Subsystem embedding Options
    !!
    CONTROL_instance%SUBSYSTEM_EMBEDDING = LowdinParameters_subsystemEmbedding
    CONTROL_instance%LOCALIZE_ORBITALS = (LowdinParameters_localizeOrbitals .or. LowdinParameters_subsystemEmbedding)
    CONTROL_instance%SUBSYSTEM_LEVEL_SHIFTING = LowdinParameters_subsystemLevelShifting
    CONTROL_instance%SUBSYSTEM_ORBITAL_THRESHOLD = LowdinParameters_subsystemOrbitalThreshold
    CONTROL_instance%SUBSYSTEM_BASIS_THRESHOLD = LowdinParameters_subsystemBasisThreshold
    CONTROL_instance%ERKALE_LOCALIZATION_METHOD = LowdinParameters_erkaleLocalizationMethod

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
    CONTROL_instance%MOLDEN_FILE_FORMAT = LowdinParameters_moldenFileFormat


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
    CONTROL_instance%AMBER_FILE = LowdinParameters_amberFile    

    !!*****************************************************                            
    !! Properties Options                                                              
    CONTROL_instance%CALCULATE_INTERPARTICLE_DISTANCES = LowdinParameters_calculateInterparticleDistances
    CONTROL_instance%CALCULATE_DENSITY_VOLUME = LowdinParameters_calculateDensityVolume

    !!*****************************************************                            
    !! Miscelaneous Options                                                            
    !!                                                                                 
    CONTROL_instance%MO_FRACTION_OCCUPATION = LowdinParameters_MOFractionOccupation
    CONTROL_instance%IONIZE_MO = LowdinParameters_ionizeMO
    CONTROL_instance%IONIZE_SPECIES = LowdinParameters_ionizeSpecies
    CONTROL_instance%EXCITE_SPECIES = LowdinParameters_exciteSpecies

    !!*****************************************************
    !! Integrals transformation options
    !!
    CONTROL_instance%INTEGRALS_TRANSFORMATION_METHOD = LowdinParameters_integralsTransformationMethod
    CONTROL_instance%IT_BUFFERSIZE = LowdinParameters_ITBuffersize

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
  subroutine CONTROL_save( unit, lastStep, firstStep )
    implicit none

    integer :: unit
    logical, optional :: lastStep 
    logical, optional :: firstStep

    !! Saving de control parameters on the name list.

    !!***************************************************************************
    !! Dummy variables, just for debugging. 
    !!
    LowdinParameters_dummyReal(:) = CONTROL_instance%DUMMY_REAL(:) 
    LowdinParameters_dummyInteger(:) = CONTROL_instance%DUMMY_INTEGER(:)
    LowdinParameters_dummyLogical(:) = CONTROL_instance%DUMMY_LOGICAL(:)  
    LowdinParameters_dummyCharacter(:) = CONTROL_instance%DUMMY_CHARACTER(:)

    !!***************************************************************************      
    !! Parameter to control Integrals library                                          
    !!                                                                                 
    LowdinParameters_tv = CONTROL_instance%TV
    LowdinParameters_integralThreshold = CONTROL_instance%INTEGRAL_THRESHOLD
    LowdinParameters_integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE
    LowdinParameters_integralStorage = CONTROL_instance%INTEGRAL_STORAGE
    LowdinParameters_integralScheme = CONTROL_instance%INTEGRAL_SCHEME
    LowdinParameters_schwarzInequality = CONTROL_instance%SCHWARZ_INEQUALITY
    LowdinParameters_harmonicConstant = CONTROL_instance%HARMONIC_CONSTANT 

    !!***************************************************************************      
    !! Parameter to control SCF program                                                
    !!                                                                                 
    LowdinParameters_nonElectronicEnergyTolerance = CONTROL_instance%NONELECTRONIC_ENERGY_TOLERANCE
    LowdinParameters_electronicEnergyTolerance = CONTROL_instance%ELECTRONIC_ENERGY_TOLERANCE
    LowdinParameters_nonelectronicDensityMatrixTolerance = CONTROL_instance%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE
    LowdinParameters_electronicDensityMatrixTolerance = CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE
    LowdinParameters_totalEnergyTolerance = CONTROL_instance%TOTAL_ENERGY_TOLERANCE
    LowdinParameters_totalDensityMatrixTolerance = CONTROL_instance%TOTAL_DENSITY_MATRIX_TOLERANCE
    LowdinParameters_densityFactorThreshold = CONTROL_instance%DENSITY_FACTOR_THRESHOLD
    LowdinParameters_diisSwitchThreshold = CONTROL_instance%DIIS_SWITCH_THRESHOLD
    LowdinParameters_diisSwitchThreshold_bkp = CONTROL_instance%DIIS_SWITCH_THRESHOLD_BKP
    LowdinParameters_electronicLevelShifting = CONTROL_instance%ELECTRONIC_LEVEL_SHIFTING
    LowdinParameters_nonelectronicLevelShifting = CONTROL_instance%NONELECTRONIC_LEVEL_SHIFTING
    LowdinParameters_exchangeOrbitalThreshold = CONTROL_instance%EXCHANGE_ORBITAL_THRESHOLD
    LowdinParameters_waveFunctionScale = CONTROL_instance%WAVE_FUNCTION_SCALE
    LowdinParameters_scfNonelectronicMaxIterations = CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS
    LowdinParameters_scfElectronicMaxIterations = CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS
    LowdinParameters_scfGlobalMaxIterations = CONTROL_instance%SCF_GLOBAL_MAX_ITERATIONS
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
    LowdinParameters_forceClosedShell = CONTROL_instance%FORCE_CLOSED_SHELL 
    LowdinParameters_debugScfs = CONTROL_instance%DEBUG_SCFS

    LowdinParameters_scfGhostSpecies = CONTROL_instance%SCF_GHOST_SPECIES

    !!*****************************************************                            
    !! Hartree-Fock Options                                                            
    !!                                                                                 
    LowdinParameters_frozen = CONTROL_instance%FROZEN_PARTICLE
    LowdinParameters_freezeNonElectronicOrbitals = CONTROL_instance%FREEZE_NON_ELECTRONIC_ORBITALS
    LowdinParameters_freezeElectronicOrbitals = CONTROL_instance%FREEZE_ELECTRONIC_ORBITALS
    LowdinParameters_hartreeProductGuess = CONTROL_instance%HARTREE_PRODUCT_GUESS
    LowdinParameters_readCoefficients = CONTROL_instance%READ_COEFFICIENTS
    LowdinParameters_readFchk = CONTROL_instance%READ_FCHK
    LowdinParameters_writeCoefficientsInBinary = CONTROL_instance%WRITE_COEFFICIENTS_IN_BINARY
    LowdinParameters_readEigenvalues = CONTROL_instance%READ_EIGENVALUES
    LowdinParameters_readEigenvaluesInBinary = CONTROL_instance%READ_EIGENVALUES_IN_BINARY
    LowdinParameters_writeEigenvaluesInBinary = CONTROL_instance%WRITE_EIGENVALUES_IN_BINARY
    LowdinParameters_noSCF = CONTROL_instance%NO_SCF
    LowdinParameters_finiteMassCorrection = CONTROL_instance%FINITE_MASS_CORRECTION
    LowdinParameters_removeTranslationalContamination = CONTROL_instance%REMOVE_TRANSLATIONAL_CONTAMINATION
    LowdinParameters_buildTwoParticlesMatrixForOneParticle = CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE
    LowdinParameters_buildMixedDensityMatrix = CONTROL_instance%BUILD_MIXED_DENSITY_MATRIX
    LowdinParameters_onlyElectronicEffect = CONTROL_instance%ONLY_ELECTRONIC_EFFECT
    LowdinParameters_electronicWaveFunctionAnalysis = CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS
    LowdinParameters_isOpenShell = CONTROL_instance%IS_OPEN_SHELL
    LowdinParameters_getGradients = CONTROL_instance%GET_GRADIENTS
    LowdinParameters_HFprintEigenvalues = CONTROL_instance%HF_PRINT_EIGENVALUES 
    LowdinParameters_HFprintEigenvectors = CONTROL_instance%HF_PRINT_EIGENVECTORS
    LowdinParameters_overlapEigenThreshold = CONTROL_instance%OVERLAP_EIGEN_THRESHOLD 
    LowdinParameters_electricField = CONTROL_instance%ELECTRIC_FIELD 
    LowdinParameters_multipoleOrder = CONTROL_instance%MULTIPOLE_ORDER 
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
    if(present(firstStep)) then
       LowdinParameters_firstStep = firstStep
    else
       LowdinParameters_firstStep = CONTROL_instance%FIRST_STEP
    end if
    if(present(lastStep)) then
       LowdinParameters_lastStep = lastStep
    else
       LowdinParameters_lastStep = CONTROL_instance%LAST_STEP
    end if
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
    !! Parameter to control MBPn theory                                                 
    !!                                                                                 
    LowdinParameters_mpCorrection = CONTROL_instance%MOLLER_PLESSET_CORRECTION
    LowdinParameters_mpFrozenCoreBoundary = CONTROL_instance%MP_FROZEN_CORE_BOUNDARY
    LowdinParameters_mpOnlyElectronicCorrection = CONTROL_instance%MP_ONLY_ELECTRONIC_CORRECTION
    LowdinParameters_epsteinNesbetCorrection = CONTROL_instance%EPSTEIN_NESBET_CORRECTION 
    !!***************************************************************************      
    !! Parameter to control cosmo method                                                 
    !!                                                                                 
    LowdinParameters_cosmo = CONTROL_instance%COSMO
    LowdinParameters_cosmo_solvent_dielectric = CONTROL_instance%COSMO_SOLVENT_DIELECTRIC 
    LowdinParameters_cosmo_scaling = CONTROL_instance%COSMO_SCALING

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
    LowdinParameters_ptFactorSS = CONTROL_instance%PT_FACTOR_SS 
    LowdinParameters_ptFactorOS = CONTROL_instance%PT_FACTOR_OS 
    LowdinParameters_ptP3Method =CONTROL_instance%PT_P3_METHOD 


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
    LowdinParameters_numberOfCIStates        = CONTROL_instance%NUMBER_OF_CI_STATES
    LowdinParameters_CIdiagonalizationMethod = CONTROL_instance%CI_DIAGONALIZATION_METHOD
    LowdinParameters_CIdiagonalDressedShift = CONTROL_instance%CI_DIAGONAL_DRESSED_SHIFT

    LowdinParameters_CIactiveSpace = CONTROL_instance%CI_ACTIVE_SPACE 
    LowdinParameters_CIstatesToPrint = CONTROL_instance%CI_STATES_TO_PRINT
    LowdinParameters_CImaxNCV = CONTROL_instance%CI_MAX_NCV 
    LowdinParameters_CIsizeOfGuessMatrix = CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX  
    LowdinParameters_CIstackSize = CONTROL_instance%CI_STACK_SIZE 
    LowdinParameters_CIJacobi = CONTROL_instance%CI_JACOBI
    LowdinParameters_CIBuildFullMatrix = CONTROL_instance%CI_BUILD_FULL_MATRIX 
    LowdinParameters_CIMadSpace = CONTROL_instance%CI_MADSPACE
    LowdinParameters_CINaturalOrbitals = CONTROL_instance%CI_NATURAL_ORBITALS
    LowdinParameters_CIPrintEigenVectorsFormat = CONTROL_instance%CI_PRINT_EIGENVECTORS_FORMAT 
    LowdinParameters_CIPrintThreshold = CONTROL_instance%CI_PRINT_THRESHOLD 
    LowdinParameters_CISCICoreSpace = CONTROL_instance%CI_SCI_CORE_SPACE 
    LowdinParameters_CISCITargetSpace = CONTROL_instance%CI_SCI_TARGET_SPACE 

    !!***************************************************************************
    !! Non-orthogonal CI
    !!
    LowdinParameters_nonOrthogonalConfigurationInteraction=CONTROL_instance%NONORTHOGONAL_CONFIGURATION_INTERACTION
    LowdinParameters_translationScanGrid=CONTROL_instance%TRANSLATION_SCAN_GRID
    LowdinParameters_rotationalScanGrid=CONTROL_instance%ROTATIONAL_SCAN_GRID
    LowdinParameters_rotationAroundZMaxAngle=CONTROL_instance%ROTATION_AROUND_Z_MAX_ANGLE
    LowdinParameters_rotationAroundZStep=CONTROL_instance%ROTATION_AROUND_Z_STEP
    LowdinParameters_nestedRotationalGrids=CONTROL_instance%NESTED_ROTATIONAL_GRIDS
    LowdinParameters_translationStep=CONTROL_instance%TRANSLATION_STEP
    LowdinParameters_nestedGridsDisplacement=CONTROL_instance%NESTED_GRIDS_DISPLACEMENT
    LowdinParameters_configurationEnergyThreshold=CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD
    LowdinParameters_configurationOverlapThreshold=CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD
    LowdinParameters_configurationMaxDisplacement=CONTROL_instance%CONFIGURATION_MAX_DISPLACEMENT
    LowdinParameters_configurationMinDisplacement=CONTROL_instance%CONFIGURATION_MIN_DISPLACEMENT
    LowdinParameters_configurationMaxNPDistance=CONTROL_instance%CONFIGURATION_MAX_NP_DISTANCE
    LowdinParameters_configurationMinPPDistance=CONTROL_instance%CONFIGURATION_MIN_PP_DISTANCE
    LowdinParameters_configurationMaxPPDistance=CONTROL_instance%CONFIGURATION_MAX_PP_DISTANCE
    LowdinParameters_configurationEquivalenceDistance=CONTROL_instance%CONFIGURATION_EQUIVALENCE_DISTANCE
    LowdinParameters_empiricalOverlapParameterA=CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_A
    LowdinParameters_empiricalOverlapParameterB=CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_B
    LowdinParameters_empiricalOverlapParameterE0=CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_E0
    LowdinParameters_empiricalOverlapParameterSc=CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_SC
    LowdinParameters_configurationUseSymmetry=CONTROL_instance%CONFIGURATION_USE_SYMMETRY
    LowdinParameters_readNOCIGeometries=CONTROL_instance%READ_NOCI_GEOMETRIES
    LowdinParameters_empiricalOverlapCorrection=CONTROL_instance%EMPIRICAL_OVERLAP_CORRECTION
    LowdinParameters_onlyFirstNOCIelements=CONTROL_instance%ONLY_FIRST_NOCI_ELEMENTS
    LowdinParameters_NOCIKineticApproximation=CONTROL_instance%NOCI_KINETIC_APPROXIMATION
    LowdinParameters_computeROCIformula=CONTROL_instance%COMPUTE_ROCI_FORMULA
    LowdinParameters_removeQDOinCI=CONTROL_instance%REMOVE_QDO_IN_CI

    !!***************************************************************************      
    !! CCSD                                                                      
    !!                                                                                 
    LowdinParameters_coupledClusterLevel = CONTROL_instance%COUPLED_CLUSTER_LEVEL

    !!*****************************************************                            
    !! Parameter to general control                                                    
    !!                                                                                 
    LowdinParameters_method = CONTROL_instance%METHOD
    LowdinParameters_transformToCenterOfMass = CONTROL_instance%TRANSFORM_TO_CENTER_OF_MASS
    LowdinParameters_areThereDummyAtoms = CONTROL_instance%ARE_THERE_DUMMY_ATOMS
    LowdinParameters_areThereQDOPotentials = CONTROL_instance%ARE_THERE_QDO_POTENTIALS
    LowdinParameters_setQDOEnergyZero = CONTROL_instance%SET_QDO_ENERGY_ZERO
    LowdinParameters_isThereExternalPotential = CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL
    LowdinParameters_isThereInterparticlePotential = CONTROL_instance%IS_THERE_INTERPARTICLE_POTENTIAL
    LowdinParameters_isThereOutput = CONTROL_instance%IS_THERE_OUTPUT
    LowdinParameters_isThereFrozenParticle = CONTROL_instance%IS_THERE_FROZEN_PARTICLE
    LowdinParameters_dimensionality = CONTROL_instance%DIMENSIONALITY

    !!*****************************************************                            
    !! Density Functional Theory Options                                               
    !!                                                                                 
    LowdinParameters_gridStorage=CONTROL_instance%GRID_STORAGE
    LowdinParameters_electronCorrelationFunctional = CONTROL_instance%ELECTRON_CORRELATION_FUNCTIONAL
    LowdinParameters_electronExchangeFunctional = CONTROL_instance%ELECTRON_EXCHANGE_FUNCTIONAL
    LowdinParameters_electronExchangeCorrelationFunctional = CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL
    LowdinParameters_nuclearElectronCorrelationFunctional = CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL
    LowdinParameters_positronElectronCorrelationFunctional = CONTROL_instance%POSITRON_ELECTRON_CORRELATION_FUNCTIONAL
    LowdinParameters_betaFunction = CONTROL_instance%BETA_FUNCTION
    LowdinParameters_gridRadialPoints = CONTROL_instance%GRID_RADIAL_POINTS
    LowdinParameters_gridAngularPoints = CONTROL_instance%GRID_ANGULAR_POINTS
    LowdinParameters_gridNumberOfShells = CONTROL_instance%GRID_NUMBER_OF_SHELLS
    LowdinParameters_finalGridRadialPoints = CONTROL_instance%FINAL_GRID_RADIAL_POINTS
    LowdinParameters_finalGridAngularPoints = CONTROL_instance%FINAL_GRID_ANGULAR_POINTS
    LowdinParameters_finalGridNumberOfShells = CONTROL_instance%FINAL_GRID_NUMBER_OF_SHELLS
    LowdinParameters_polarizationOrder = CONTROL_instance%POLARIZATION_ORDER
    LowdinParameters_numberOfBlocksInAuxiliaryFunctions = CONTROL_instance%NUMBER_OF_BLOCKS_IN_AUXILIARY_FUNCTIONS
    LowdinParameters_fukuiFunctions = CONTROL_instance%FUKUI_FUNCTIONS
    LowdinParameters_auxiliaryDensity = CONTROL_instance%AUXILIARY_DENSITY
    LowdinParameters_storeThreeCenterElectronIntegrals = CONTROL_instance%STORE_THREE_CENTER_ELECTRON_INTEGRALS
    LowdinParameters_callLibxc = CONTROL_instance%CALL_LIBXC
    LowdinParameters_nuclearElectronDensityThreshold=CONTROL_instance%NUCLEAR_ELECTRON_DENSITY_THRESHOLD
    LowdinParameters_electronDensityThreshold=CONTROL_instance%ELECTRON_DENSITY_THRESHOLD
    LowdinParameters_gridWeightThreshold=CONTROL_instance%GRID_WEIGHT_THRESHOLD
    LowdinParameters_betaParameterA=CONTROL_instance%BETA_PARAMETER_A 
    LowdinParameters_betaParameterB=CONTROL_instance%BETA_PARAMETER_B
    LowdinParameters_betaParameterC=CONTROL_instance%BETA_PARAMETER_C

    !!*****************************************************
    !! Subsystem embedding Options
    !!
    LowdinParameters_subsystemEmbedding=CONTROL_instance%SUBSYSTEM_EMBEDDING
    LowdinParameters_localizeOrbitals=CONTROL_instance%LOCALIZE_ORBITALS
    LowdinParameters_subsystemLevelShifting=CONTROL_instance%SUBSYSTEM_LEVEL_SHIFTING
    LowdinParameters_subsystemOrbitalThreshold=CONTROL_instance%SUBSYSTEM_ORBITAL_THRESHOLD
    LowdinParameters_subsystemBasisThreshold=CONTROL_instance%SUBSYSTEM_BASIS_THRESHOLD
    LowdinParameters_erkaleLocalizationMethod=CONTROL_instance%ERKALE_LOCALIZATION_METHOD

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
    LowdinParameters_moldenFileFormat = CONTROL_instance%MOLDEN_FILE_FORMAT 
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
    LowdinParameters_amberFile = CONTROL_instance%AMBER_FILE     

    !!*****************************************************                            
    !! Properties Options                                                              
    LowdinParameters_calculateInterparticleDistances = CONTROL_instance%CALCULATE_INTERPARTICLE_DISTANCES
    LowdinParameters_calculateDensityVolume = CONTROL_instance%CALCULATE_DENSITY_VOLUME

    !!*****************************************************                            
    !! Miscelaneous Options                                                            
    !!                                                                                 
    LowdinParameters_MOFractionOccupation = CONTROL_instance%MO_FRACTION_OCCUPATION
    LowdinParameters_ionizeMO = CONTROL_instance%IONIZE_MO
    LowdinParameters_ionizeSpecies = CONTROL_instance%IONIZE_SPECIES
    LowdinParameters_exciteSpecies = CONTROL_instance%EXCITE_SPECIES

    !!*****************************************************
    !! Integrals transformation options
    !!
    LowdinParameters_integralsTransformationMethod = CONTROL_instance%INTEGRALS_TRANSFORMATION_METHOD 
    LowdinParameters_ITBuffersize = CONTROL_instance%IT_BUFFERSIZE 

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

    !!***************************************************************************
    !! Dummy variables, just for debugging. 
    !!
    otherThis%DUMMY_REAL(:) = this%DUMMY_REAL(:)
    otherThis%DUMMY_INTEGER(:) = this%DUMMY_INTEGER(:)
    otherThis%DUMMY_LOGICAL(:) = this%DUMMY_LOGICAL(:)
    otherThis%DUMMY_CHARACTER(:) = this%DUMMY_CHARACTER(:)

    !!*****************************************************
    !! Variables para control de integrales
    !!
    otherThis%TV = this%TV 
    otherThis%INTEGRAL_THRESHOLD = this%INTEGRAL_THRESHOLD 
    otherThis%INTEGRAL_STORAGE = this%INTEGRAL_STORAGE 
    otherThis%INTEGRAL_SCHEME = this%INTEGRAL_SCHEME
    otherThis%INTEGRAL_STACK_SIZE = this%INTEGRAL_STACK_SIZE 
    otherThis%SCHWARZ_INEQUALITY = this%SCHWARZ_INEQUALITY
    otherThis%HARMONIC_CONSTANT = this%HARMONIC_CONSTANT 

    !!***************************************************************************
    !! Parametros para control de proceso de minizacion de energia mediante
    !! metodo SCF
    !!
    otherThis%NONELECTRONIC_ENERGY_TOLERANCE = this%NONELECTRONIC_ENERGY_TOLERANCE 
    otherThis%ELECTRONIC_ENERGY_TOLERANCE = this%ELECTRONIC_ENERGY_TOLERANCE 
    otherThis%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE = this%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE 
    otherThis%ELECTRONIC_DENSITY_MATRIX_TOLERANCE = this%ELECTRONIC_DENSITY_MATRIX_TOLERANCE 
    otherThis%TOTAL_ENERGY_TOLERANCE = this%TOTAL_ENERGY_TOLERANCE 
    otherThis%TOTAL_DENSITY_MATRIX_TOLERANCE = this%TOTAL_DENSITY_MATRIX_TOLERANCE 
    otherThis%DENSITY_FACTOR_THRESHOLD = this%DENSITY_FACTOR_THRESHOLD 
    otherThis%DIIS_SWITCH_THRESHOLD = this%DIIS_SWITCH_THRESHOLD 
    otherThis%DIIS_SWITCH_THRESHOLD_BKP = this%DIIS_SWITCH_THRESHOLD_BKP 
    otherThis%ELECTRONIC_LEVEL_SHIFTING = this%ELECTRONIC_LEVEL_SHIFTING 
    otherThis%NONELECTRONIC_LEVEL_SHIFTING = this%NONELECTRONIC_LEVEL_SHIFTING 
    otherthis%EXCHANGE_ORBITAL_THRESHOLD = this%EXCHANGE_ORBITAL_THRESHOLD
    otherThis%WAVE_FUNCTION_SCALE= this%WAVE_FUNCTION_SCALE
    otherThis%SCF_NONELECTRONIC_MAX_ITERATIONS = this%SCF_NONELECTRONIC_MAX_ITERATIONS 
    otherThis%SCF_ELECTRONIC_MAX_ITERATIONS = this%SCF_ELECTRONIC_MAX_ITERATIONS 
    otherThis%SCF_GLOBAL_MAX_ITERATIONS = this%SCF_GLOBAL_MAX_ITERATIONS 
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
    otherThis%FORCE_CLOSED_SHELL = this%FORCE_CLOSED_SHELL
    otherThis%DEBUG_SCFS = this%DEBUG_SCFS 
    otherThis%SCF_GHOST_SPECIES = this%SCF_GHOST_SPECIES 
    !!***************************************************************************
    !! Parametros para control Hartree-Fock
    !!
    otherThis%FROZEN_PARTICLE = this%FROZEN_PARTICLE 
    otherThis%FREEZE_NON_ELECTRONIC_ORBITALS = this%FREEZE_NON_ELECTRONIC_ORBITALS 
    otherThis%FREEZE_ELECTRONIC_ORBITALS = this%FREEZE_ELECTRONIC_ORBITALS 
    otherThis%HARTREE_PRODUCT_GUESS = this%HARTREE_PRODUCT_GUESS 
    otherThis%READ_COEFFICIENTS = this%READ_COEFFICIENTS 
    otherThis%READ_FCHK = this%READ_FCHK
    otherThis%WRITE_COEFFICIENTS_IN_BINARY = this%WRITE_COEFFICIENTS_IN_BINARY
    otherThis%NO_SCF = this%NO_SCF 
    otherThis%FINITE_MASS_CORRECTION = this%FINITE_MASS_CORRECTION 
    otherThis%REMOVE_TRANSLATIONAL_CONTAMINATION = this%REMOVE_TRANSLATIONAL_CONTAMINATION 
    otherThis%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE = this%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE 
    otherThis%BUILD_MIXED_DENSITY_MATRIX = this%BUILD_MIXED_DENSITY_MATRIX 
    otherThis%ONLY_ELECTRONIC_EFFECT = this%ONLY_ELECTRONIC_EFFECT 
    otherThis%ELECTRONIC_WAVEFUNCTION_ANALYSIS = this%ELECTRONIC_WAVEFUNCTION_ANALYSIS 
    otherThis%IS_OPEN_SHELL = this%IS_OPEN_SHELL    
    otherThis%GET_GRADIENTS = this%GET_GRADIENTS    
    otherThis%HF_PRINT_EIGENVALUES = this%HF_PRINT_EIGENVALUES 
    otherThis%HF_PRINT_EIGENVECTORS = this%HF_PRINT_EIGENVECTORS
    otherThis%OVERLAP_EIGEN_THRESHOLD = this%OVERLAP_EIGEN_THRESHOLD 
    otherThis%ELECTRIC_FIELD = this%ELECTRIC_FIELD 
    otherThis%MULTIPOLE_ORDER = this%MULTIPOLE_ORDER 
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
    otherThis%LAST_STEP = this%LAST_STEP
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
    otherThis%EPSTEIN_NESBET_CORRECTION = this%EPSTEIN_NESBET_CORRECTION
    !!*****************************************************
    !! Control de parametros de metodo cosmo
    !!
    otherThis%COSMO  = this%COSMO                   
    otherThis%COSMO_SOLVENT_DIELECTRIC = this%COSMO_SOLVENT_DIELECTRIC
    otherThis%COSMO_SCALING = this%COSMO_SCALING

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
    otherThis%PT_FACTOR_SS = this%PT_FACTOR_SS
    otherThis%PT_FACTOR_OS = this%PT_FACTOR_OS
    otherThis%PT_P3_METHOD = this%PT_P3_METHOD 

    !!*****************************************************
    !! Control parametros de formato
    !!
    otherThis%FORMAT_NUMBER_OF_COLUMNS = this%FORMAT_NUMBER_OF_COLUMNS 
    otherThis%PRINT_LEVEL = this%PRINT_LEVEL
    otherThis%UNITS = this%UNITS 
    otherThis%UNIT_FOR_OUTPUT_FILE = this%UNIT_FOR_OUTPUT_FILE 
    otherThis%UNIT_FOR_MP2_INTEGRALS_FILE = this%UNIT_FOR_MP2_INTEGRALS_FILE 
    otherThis%UNIT_FOR_MOLECULAR_ORBITALS_FILE = this%UNIT_FOR_MOLECULAR_ORBITALS_FILE 
    !!***************************************************************************
    !! CISD - FCI
    !!
    otherThis%CONFIGURATION_INTERACTION_LEVEL = this%CONFIGURATION_INTERACTION_LEVEL 
    otherThis%NUMBER_OF_CI_STATES       = this%NUMBER_OF_CI_STATES
    otherThis%CI_DIAGONALIZATION_METHOD = this%CI_DIAGONALIZATION_METHOD
    otherThis%CI_DIAGONAL_DRESSED_SHIFT = this%CI_DIAGONAL_DRESSED_SHIFT
    otherThis%CI_ACTIVE_SPACE =  this%CI_ACTIVE_SPACE 
    otherThis%CI_STATES_TO_PRINT =  this%CI_STATES_TO_PRINT
    otherThis%CI_MAX_NCV = this%CI_MAX_NCV
    otherThis%CI_SIZE_OF_GUESS_MATRIX = this%CI_SIZE_OF_GUESS_MATRIX
    otherThis%CI_STACK_SIZE = this%CI_STACK_SIZE 
    otherThis%CI_CONVERGENCE = this%CI_CONVERGENCE
    otherThis%CI_MATVEC_TOLERANCE = this%CI_MATVEC_TOLERANCE 
    otherThis%CI_SAVE_EIGENVECTOR = this%CI_SAVE_EIGENVECTOR
    otherThis%CI_LOAD_EIGENVECTOR = this%CI_LOAD_EIGENVECTOR
    otherThis%CI_JACOBI = this%CI_JACOBI
    otherThis%CI_BUILD_FULL_MATRIX = this%CI_BUILD_FULL_MATRIX
    otherThis%CI_MADSPACE = this%CI_MADSPACE
    otherThis%CI_NATURAL_ORBITALS = this%CI_NATURAL_ORBITALS
    otherThis%CI_PRINT_EIGENVECTORS_FORMAT = this%CI_PRINT_EIGENVECTORS_FORMAT 
    otherThis%CI_PRINT_THRESHOLD = this%CI_PRINT_THRESHOLD 
    otherThis%CI_SCI_CORE_SPACE = this%CI_SCI_CORE_SPACE 
    otherThis%CI_SCI_TARGET_SPACE = this%CI_SCI_TARGET_SPACE 

    !!***************************************************************************
    !! Non-orthogonal CI
    !!
    otherThis%NONORTHOGONAL_CONFIGURATION_INTERACTION = this%NONORTHOGONAL_CONFIGURATION_INTERACTION
    otherThis%TRANSLATION_SCAN_GRID = this%TRANSLATION_SCAN_GRID
    otherThis%ROTATIONAL_SCAN_GRID = this%ROTATIONAL_SCAN_GRID
    otherThis%ROTATION_AROUND_Z_MAX_ANGLE=this%ROTATION_AROUND_Z_MAX_ANGLE
    otherThis%ROTATION_AROUND_Z_STEP=this%ROTATION_AROUND_Z_STEP
    otherThis%NESTED_ROTATIONAL_GRIDS = this%NESTED_ROTATIONAL_GRIDS
    otherThis%TRANSLATION_STEP = this%TRANSLATION_STEP
    otherThis%NESTED_GRIDS_DISPLACEMENT = this%NESTED_GRIDS_DISPLACEMENT
    otherThis%CONFIGURATION_ENERGY_THRESHOLD=this%CONFIGURATION_ENERGY_THRESHOLD
    otherThis%CONFIGURATION_OVERLAP_THRESHOLD=this%CONFIGURATION_OVERLAP_THRESHOLD
    otherThis%CONFIGURATION_MAX_DISPLACEMENT=this%CONFIGURATION_MAX_DISPLACEMENT
    otherThis%CONFIGURATION_MIN_DISPLACEMENT=this%CONFIGURATION_MIN_DISPLACEMENT
    otherThis%CONFIGURATION_MAX_NP_DISTANCE=this%CONFIGURATION_MAX_NP_DISTANCE
    otherThis%CONFIGURATION_MIN_PP_DISTANCE=this%CONFIGURATION_MIN_PP_DISTANCE
    otherThis%CONFIGURATION_MAX_PP_DISTANCE=this%CONFIGURATION_MAX_PP_DISTANCE
    otherThis%CONFIGURATION_EQUIVALENCE_DISTANCE=this%CONFIGURATION_EQUIVALENCE_DISTANCE
    otherThis%CONFIGURATION_USE_SYMMETRY=this%CONFIGURATION_USE_SYMMETRY
    otherThis%READ_NOCI_GEOMETRIES=this%READ_NOCI_GEOMETRIES
    otherThis%ONLY_FIRST_NOCI_ELEMENTS=this%ONLY_FIRST_NOCI_ELEMENTS
    otherThis%NOCI_KINETIC_APPROXIMATION=this%NOCI_KINETIC_APPROXIMATION
    otherThis%COMPUTE_ROCI_FORMULA=this%COMPUTE_ROCI_FORMULA
    otherThis%REMOVE_QDO_IN_CI=this%REMOVE_QDO_IN_CI

    !!***************************************************************************
    !! CCSD
    !!
    otherThis%COUPLED_CLUSTER_LEVEL = this%COUPLED_CLUSTER_LEVEL 
    !!*****************************************************
    !! Parametros de control general
    !!
    otherThis%DIMENSIONALITY = this%DIMENSIONALITY 
    otherThis%DOUBLE_ZERO_THRESHOLD  = this%DOUBLE_ZERO_THRESHOLD  
    otherThis%METHOD = this%METHOD 
    otherThis%TRANSFORM_TO_CENTER_OF_MASS = this%TRANSFORM_TO_CENTER_OF_MASS 
    otherThis%ARE_THERE_DUMMY_ATOMS = this%ARE_THERE_DUMMY_ATOMS 
    otherThis%ARE_THERE_QDO_POTENTIALS = this%ARE_THERE_QDO_POTENTIALS 
    otherThis%SET_QDO_ENERGY_ZERO = this%SET_QDO_ENERGY_ZERO
    otherThis%IS_THERE_EXTERNAL_POTENTIAL = this%IS_THERE_EXTERNAL_POTENTIAL 
    otherThis%IS_THERE_INTERPARTICLE_POTENTIAL = this%IS_THERE_INTERPARTICLE_POTENTIAL
    otherThis%IS_THERE_OUTPUT = this%IS_THERE_OUTPUT
    otherThis%IS_THERE_FROZEN_PARTICLE = this%IS_THERE_FROZEN_PARTICLE 
    !!*****************************************************
    !! Density Functional Theory Options
    !!
    otherThis%AUXILIARY_DENSITY = this%AUXILIARY_DENSITY 
    otherThis%CALL_LIBXC = this%CALL_LIBXC
    otherThis%GRID_STORAGE=this%GRID_STORAGE
    otherThis%ELECTRON_CORRELATION_FUNCTIONAL = this%ELECTRON_CORRELATION_FUNCTIONAL 
    otherThis%ELECTRON_EXCHANGE_FUNCTIONAL = this%ELECTRON_EXCHANGE_FUNCTIONAL 
    otherThis%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL = this%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL 
    otherThis%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL = this%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL 
    otherThis%POSITRON_ELECTRON_CORRELATION_FUNCTIONAL = this%POSITRON_ELECTRON_CORRELATION_FUNCTIONAL 
    otherThis%BETA_FUNCTION = this%BETA_FUNCTION
    otherThis%GRID_RADIAL_POINTS=this%GRID_RADIAL_POINTS
    otherThis%GRID_ANGULAR_POINTS=this%GRID_ANGULAR_POINTS
    otherThis%GRID_NUMBER_OF_SHELLS=this%GRID_NUMBER_OF_SHELLS
    otherThis%FINAL_GRID_RADIAL_POINTS=this%FINAL_GRID_RADIAL_POINTS
    otherThis%FINAL_GRID_ANGULAR_POINTS=this%FINAL_GRID_ANGULAR_POINTS
    otherThis%FINAL_GRID_NUMBER_OF_SHELLS=this%FINAL_GRID_NUMBER_OF_SHELLS
    otherThis%STORE_THREE_CENTER_ELECTRON_INTEGRALS = this%STORE_THREE_CENTER_ELECTRON_INTEGRALS 
    otherThis%POLARIZATION_ORDER = this%POLARIZATION_ORDER 
    otherThis%FUKUI_FUNCTIONS = this%FUKUI_FUNCTIONS 
    otherThis%NUMBER_OF_BLOCKS_IN_AUXILIARY_FUNCTIONS = this%NUMBER_OF_BLOCKS_IN_AUXILIARY_FUNCTIONS 
    otherThis%NUCLEAR_ELECTRON_DENSITY_THRESHOLD= this%NUCLEAR_ELECTRON_DENSITY_THRESHOLD
    otherThis%ELECTRON_DENSITY_THRESHOLD= this%ELECTRON_DENSITY_THRESHOLD
    otherThis%GRID_WEIGHT_THRESHOLD= this%GRID_WEIGHT_THRESHOLD
    !!*****************************************************
    !! Subsystem embedding Options
    !!
    otherThis%SUBSYSTEM_EMBEDDING = this%SUBSYSTEM_EMBEDDING
    otherThis%LOCALIZE_ORBITALS = this%LOCALIZE_ORBITALS
    otherThis%SUBSYSTEM_LEVEL_SHIFTING = this%SUBSYSTEM_LEVEL_SHIFTING
    otherThis%SUBSYSTEM_ORBITAL_THRESHOLD = this%SUBSYSTEM_ORBITAL_THRESHOLD
    otherThis%SUBSYSTEM_BASIS_THRESHOLD = this%SUBSYSTEM_BASIS_THRESHOLD
    otherThis%ERKALE_LOCALIZATION_METHOD = this%ERKALE_LOCALIZATION_METHOD

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
    otherThis%MOLDEN_FILE_FORMAT = this%MOLDEN_FILE_FORMAT 
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
    otherThis%AMBER_FILE = this%AMBER_FILE    
    !!*****************************************************
    !! Properties Options
    otherThis%CALCULATE_INTERPARTICLE_DISTANCES  = this%CALCULATE_INTERPARTICLE_DISTANCES  
    otherThis%CALCULATE_DENSITY_VOLUME  = this%CALCULATE_DENSITY_VOLUME  
    !!*****************************************************
    !! Miscelaneous Options
    !!
    otherThis%MO_FRACTION_OCCUPATION = this%MO_FRACTION_OCCUPATION 
    otherThis%IONIZE_MO = this%IONIZE_MO 
    otherThis%IONIZE_SPECIES = this%IONIZE_SPECIES
    otherThis%EXCITE_SPECIES = this%EXCITE_SPECIES

    !!*****************************************************
    !! Integrals transformation options
    !!
    otherThis%INTEGRALS_TRANSFORMATION_METHOD = this%INTEGRALS_TRANSFORMATION_METHOD 
    otherThis%IT_BUFFERSIZE = this%IT_BUFFERSIZE 

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

    integer:: i,j !, nthreads, proc

    print *,""
    print *,"LOWDIN IS RUNNING WITH NEXT PARAMETERS: "
    print *,"----------------------------------------"
    print *,""

    write (*,"(T10,A)") "METHOD TYPE:  "//trim(CONTROL_instance%METHOD)

    write (*,"(T10,A,I5)") "NUMBER OF CORES: ", CONTROL_instance%NUMBER_OF_CORES

    if(CONTROL_instance%METHOD=="RKS" .or. CONTROL_instance%METHOD=="UKS" .or. CONTROL_instance%METHOD=="ROKS" ) then

       if(CONTROL_instance%AUXILIARY_DENSITY) write (*,"(T10,A)") "USING AUXILIARY DENSITY"

       if(CONTROL_instance%GRID_STORAGE .eq. "DISK") then
          write (*,"(T10,A)") "STORING DENSITY GRIDS IN DISK"
       else
          write (*,"(T10,A)") "STORING DENSITY GRIDS IN MEMORY"          
       end if
       
       if(CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL .ne. "NONE") then
          write (*,"(T10,A)") "ELECTRON EXCHANGE CORRELATION FUNCTIONAL: "//trim(CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL)
       else
          write (*,"(T10,A)") "ELECTRON CORRELATION FUNCTIONAL: "//trim(CONTROL_instance%ELECTRON_CORRELATION_FUNCTIONAL)
          write (*,"(T10,A)") "ELECTRON EXCHANGE FUNCTIONAL: "//trim(CONTROL_instance%ELECTRON_EXCHANGE_FUNCTIONAL)
       end if
       
       write (*,"(T10,A)") "ELECTRON-NUCLEAR CORRELATION FUNCTIONAL: "//trim(CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL)
       write (*,"(T10,A)") "ELECTRON-POSITRON CORRELATION FUNCTIONAL: "//trim(CONTROL_instance%POSITRON_ELECTRON_CORRELATION_FUNCTIONAL)
       write (*,"(T10,A,I5,A,I5)") "SCF ATOMIC RADIALxANGULAR GRID SIZE:",CONTROL_instance%GRID_RADIAL_POINTS,"x",CONTROL_instance%GRID_ANGULAR_POINTS
       if( CONTROL_instance%FINAL_GRID_ANGULAR_POINTS*CONTROL_instance%FINAL_GRID_RADIAL_POINTS  .gt. &
            CONTROL_instance%GRID_ANGULAR_POINTS*CONTROL_instance%GRID_RADIAL_POINTS) then
          write (*,"(T10,A,I5,A,I5)") "FINAL ATOMIC RADIALxANGULAR GRID SIZE:",CONTROL_instance%FINAL_GRID_RADIAL_POINTS,"x",CONTROL_instance%FINAL_GRID_ANGULAR_POINTS
       end if
       

       ! if(CONTROL_instance%STORE_THREE_CENTER_ELECTRON_INTEGRALS) then

       !    write (*,"(T10,A)") "STORING THREE CENTER ELECTRON INTEGRALS IN DISK"

       ! else

       !    write (*,"(T10,A)") "CALCULATING THREE CENTER ELECTRON INTEGRALS ON THE FLY"

       ! end if

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

    if(CONTROL_instance%EPSTEIN_NESBET_CORRECTION>=2) then

       write (*,"(T10,A,I5)") "EPSTEIN NESBET CORRECTION:  ",CONTROL_instance%EPSTEIN_NESBET_CORRECTION

    end if



    if(CONTROL_instance%COSMO) then

       write (*,"(T10,A)") "COSMO:  T "

    end if

    
    if(CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE" ) then

       write (*,"(T10,A,A)") "CONFIGURATION INTERACTION LEVEL:  ", CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL
       ! CONTROL_instance%ELECTRONIC_ENERGY_TOLERANCE = 1E-08
       ! CONTROL_instance%NONELECTRONIC_ENERGY_TOLERANCE = 1E-08

    end if

    !!***************************************************************************
    !! Non-orthogonal CI
    !!
    if(CONTROL_instance%NONORTHOGONAL_CONFIGURATION_INTERACTION ) then
       print *, ""
       write (*,"(T10,A)") "PERFORMING A NONORTHOGONAL CONFIGURATION INTERACTION CALCULATION"
       write (*,"(T10,A)") "THAT MIXES HF CALCULATIONS WITH DIFFERENT BASIS SET CENTERS "

       if(CONTROL_instance%UNITS .eq. "ANGS") then
          CONTROL_instance%TRANSLATION_STEP = CONTROL_instance%TRANSLATION_STEP / ANGSTROM
          CONTROL_instance%NESTED_GRIDS_DISPLACEMENT = CONTROL_instance%NESTED_GRIDS_DISPLACEMENT / ANGSTROM
          CONTROL_instance%CONFIGURATION_MAX_DISPLACEMENT = CONTROL_instance%CONFIGURATION_MAX_DISPLACEMENT / ANGSTROM
          CONTROL_instance%CONFIGURATION_MIN_DISPLACEMENT = CONTROL_instance%CONFIGURATION_MIN_DISPLACEMENT / ANGSTROM
          CONTROL_instance%CONFIGURATION_EQUIVALENCE_DISTANCE = CONTROL_instance%CONFIGURATION_EQUIVALENCE_DISTANCE / ANGSTROM
          CONTROL_instance%CONFIGURATION_MAX_NP_DISTANCE = CONTROL_instance%CONFIGURATION_MAX_NP_DISTANCE / ANGSTROM
          CONTROL_instance%CONFIGURATION_MIN_PP_DISTANCE = CONTROL_instance%CONFIGURATION_MIN_PP_DISTANCE / ANGSTROM
          CONTROL_instance%CONFIGURATION_MAX_PP_DISTANCE = CONTROL_instance%CONFIGURATION_MAX_PP_DISTANCE / ANGSTROM
       end if
       
       if(sum(CONTROL_instance%TRANSLATION_SCAN_GRID) .gt. 0 ) then
          write (*,"(T10,A,I6,A)") "THE BASIS FUNCTIONS AT EACH TRANSLATION CENTER WILL BE DISPLACED FORMING ", CONTROL_instance%TRANSLATION_SCAN_GRID(1)*CONTROL_instance%TRANSLATION_SCAN_GRID(2)*CONTROL_instance%TRANSLATION_SCAN_GRID(3) ," BODY-CENTERED-CUBIC CELLS"
          write (*,"(T10,A,I3,I3,I3,A)") "IN A RECTANGULAR ARRAY OF", CONTROL_instance%TRANSLATION_SCAN_GRID(1:3), " CELLS ALONG THE X,Y,Z AXIS"
          write (*,"(T10,A,F6.3,A10)") "WITH A SEPARATION BETWEEN POINTS OF", CONTROL_instance%TRANSLATION_STEP, " BOHRS"
       end if

       if(CONTROL_instance%ROTATIONAL_SCAN_GRID .gt. 0 ) then
          if(CONTROL_instance%NESTED_ROTATIONAL_GRIDS .gt. 1 ) then
             write (*,"(T10,I3,A,I6,A)") CONTROL_instance%NESTED_ROTATIONAL_GRIDS, " LEBEDEV GRIDS OF", CONTROL_instance%ROTATIONAL_SCAN_GRID, " BASIS FUNCTIONS WILL BE PLACED AROUND EACH ROTATIONAL CENTER"
             write (*,"(T10,A,F6.3,A10)") "WITH A RADIAL SEPARATION  OF", CONTROL_instance%NESTED_GRIDS_DISPLACEMENT,  " BOHRS"
          else
             write (*,"(T10,A,I6,A)") "A LEBEDEV GRID OF", CONTROL_instance%ROTATIONAL_SCAN_GRID, " BASIS FUNCTIONS WILL BE PLACED AROUND EACH ROTATIONAL CENTER"
          end if
       end if

       if(CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD .gt. 0.0) &
            write (*,"(T10,A,ES15.5,A)") "GEOMETRIES WITH ENERGY HIGHER THAN ", CONTROL_instance%CONFIGURATION_ENERGY_THRESHOLD," A.U. COMPARED TO THE HF REFERENCE WILL NOT BE INCLUDED IN THE CONFIGURATION SPACE"  
       
       if(CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD .gt. 0.0) &
            write (*,"(T10,A,ES15.5)") "SKIPPING HAMILTONIAN MATRIX ELEMENTS FROM CONFIGURATIONS WITH OVERLAP LOWER THAN ",  CONTROL_instance%CONFIGURATION_OVERLAP_THRESHOLD

       if(sum(CONTROL_instance%CONFIGURATION_MAX_DISPLACEMENT) .gt. 0.0) &
            write (*,"(T10,A,3ES15.5,A10,A)") "SKIPPING GEOMETRIES OUTSIDE AN ELLIPSOID OF SEMI-AXES ",  CONTROL_instance%CONFIGURATION_MAX_DISPLACEMENT,  " BOHRS", " FROM THE REFERENCE CENTER"

       if(sum(CONTROL_instance%CONFIGURATION_MIN_DISPLACEMENT) .gt. 0.0) &
            write (*,"(T10,A,3ES15.5,A10,A)") "SKIPPING GEOMETRIES INSIDE AN ELLIPSOID OF SEMI-AXES ",  CONTROL_instance%CONFIGURATION_MIN_DISPLACEMENT,  " BOHRS", " FROM THE REFERENCE CENTER"

       if(CONTROL_instance%CONFIGURATION_MAX_NP_DISTANCE .lt. 1.0E8) &
            write (*,"(T10,A,ES15.5,A10,A)") "SKIPPING GEOMETRIES WITH SEPARATION BETWEEN POSITIVE AND NEGATIVE BASIS SET CENTERS HIGHER THAN ",  CONTROL_instance%CONFIGURATION_MAX_NP_DISTANCE, " BOHRS"

       if(CONTROL_instance%CONFIGURATION_MIN_PP_DISTANCE .gt. 0.0) &
            write (*,"(T10,A,ES15.5,A10,A)") "SKIPPING GEOMETRIES WITH SEPARATION BETWEEN SAME CHARGE BASIS SET CENTERS LOWER THAN ",  CONTROL_instance%CONFIGURATION_MIN_PP_DISTANCE, " BOHRS"

       if(CONTROL_instance%CONFIGURATION_MAX_PP_DISTANCE .lt. 1.0E8) &
            write (*,"(T10,A,ES15.5,A10,A)") "SKIPPING GEOMETRIES WITH SEPARATION BETWEEN SAME CHARGE BASIS SET CENTERS HIGHER THAN ",  CONTROL_instance%CONFIGURATION_MAX_PP_DISTANCE, " BOHRS"
       
       if(CONTROL_instance%CONFIGURATION_EQUIVALENCE_DISTANCE .gt. 0.0) &
            write (*,"(T10,A,ES15.5,A10,A)") "GEOMETRIES WITH DISTANCE MATRIX THAT DIFFER IN LESS THAN ",  CONTROL_instance%CONFIGURATION_EQUIVALENCE_DISTANCE, " BOHRS", " WILL BE REGARDED AS EQUIVALENT"
       
       if(CONTROL_instance%CONFIGURATION_USE_SYMMETRY) &
            write (*,"(T10,A)") "CONFIGURATION PAIRS WILL BE CLASSIFIED ACCORDING TO THEIR MIXED GEOMETRY DISTANCE MATRIX AND ONLY UNIQUE ELEMENTS WILL BE COMPUTED"

       if(CONTROL_instance%READ_NOCI_GEOMETRIES) &
            write (*,"(T10,A,A,A)") "GEOMETRIES FOR THE NOCI EXPANSION WILL BE READ FROM ",trim(CONTROL_instance%INPUT_FILE)//"NOCI.coords" ," FILE"

       if(CONTROL_instance%EMPIRICAL_OVERLAP_CORRECTION) then
          write (*,"(T10,A,F8.5,A,F8.5)") &
               "SCALING NOCI OVERLAP AND HAMILTONIAN ELEMENTS ACCORDING TO S_I,II'=a(S_I,II)^b with a=",&
               CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_A,&
               " and b=",&
               CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_B
       end if
       
       if(CONTROL_instance%ONLY_FIRST_NOCI_ELEMENTS) &
            write (*,"(T10,A)") "COMPUTING NOCI ELEMENTS ONLY WITH RESPECT TO THE FIRST GEOMETRY - YOU HAVE TO SOLVE THE CI EQUATION MANUALLY!"

       if(CONTROL_instance%NOCI_KINETIC_APPROXIMATION) &
            write (*,"(T10,A)") "IN THIS NOCI CALCULATION, ONLY OVERLAP AND KINETIC ENERGY WILL BE COMPUTED, OTHER ENERGY CONTRIBUTIONS WILL BE APPROXIMATED AS EAB=SAB/2(EA+EB)" 
       
       if(CONTROL_instance%COMPUTE_ROCI_FORMULA) then
          CONTROL_instance%ONLY_FIRST_NOCI_ELEMENTS=.true.
          if(CONTROL_instance%ROTATION_AROUND_Z_MAX_ANGLE .gt. 180 ) CONTROL_instance%ROTATION_AROUND_Z_MAX_ANGLE=180
          write (*,"(T10,A)") "COMPUTING ROTATIONAL ENERGIES FROM THE FIRST GEOMETRY NOCI ELEMENTS"
          if(CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_E0 .gt. 0.0 .or. CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_B .gt. 0.0) then
             write (*,"(T10,A,F8.5,A,F8.5)") &
                  "EMPLOYING EMPIRICAL SCALE FACTORS E0=",&
                  CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_E0,&
                  " AND Sc=",&
                  CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_SC
          end if
          print *, ""
       end if

       if(CONTROL_instance%REMOVE_QDO_IN_CI) &
            write (*,"(T10,A)") "COMPUTING SCF WITH QDO POTENTIALS AND REMOVING THEM IN THE POST-SCF CALCULATION"

       if(CONTROL_instance%ROTATION_AROUND_Z_STEP .gt. 0 ) then
          ! if(CONTROL_instance%NESTED_ROTATIONAL_GRIDS .gt. 1 ) then
          !    write (*,"(T10,I3,A,I6,A)") CONTROL_instance%NESTED_ROTATIONAL_GRIDS, "  GRIDS OF", CONTROL_instance%ROTATIONAL_SCAN_GRID_AROUND_Z, " BASIS FUNCTIONS WILL BE PLACED AROUND EACH ROTATIONAL CENTER"
          !    write (*,"(T10,A,F6.3,A10)") "WITH A RADIAL SEPARATION  OF", CONTROL_instance%NESTED_GRIDS_DISPLACEMENT,  " BOHRS"
          ! else
          write (*,"(T10,A,F8.2,A,I6,A)") "THE MOLECULAR SYSTEM WILL BE ROTATED AROUND THE Z AXIS IN STEPS OF", CONTROL_instance%ROTATION_AROUND_Z_STEP, " DEGREES UP TO ",  CONTROL_instance%ROTATION_AROUND_Z_MAX_ANGLE, " DEGREES"
          ! end if
       end if
       

       print *, ""

    end if

    if(CONTROL_instance%COUPLED_CLUSTER_LEVEL /= "NONE" ) then

       write (*,"(T10,A,A)") "COUPLED CLUSTER LEVEL:  ", CONTROL_instance%COUPLED_CLUSTER_LEVEL

    end if

    if(CONTROL_instance%PT_ORDER>=2) then

       write (*,"(T10,A,I5)") "PROPAGATOR THEORY ORDER:  ",CONTROL_instance%PT_ORDER

    end if

    if((CONTROL_instance%IONIZE_SPECIES(1)) /= "NONE") then 
       ! print *, "size ionizepsecie", size(CONTROL_instance%IONIZE_SPECIE)
       ! print *, "ionizepsecie", CONTROL_instance%IONIZE_SPECIE
       write (*,"(T10,A)") "MOLECULAR ORBITALS TO BE IONIZED "
       do i = 1, size(CONTROL_instance%IONIZE_SPECIES)
          if ( trim(CONTROL_instance%IONIZE_SPECIES(i)) /= "NONE" ) then
             write (*,"(T15,A,A)") "FOR SPECIES: ", trim(CONTROL_instance%IONIZE_SPECIES(i))
             do j = 1, size(CONTROL_instance%IONIZE_MO)
                if(CONTROL_instance%IONIZE_MO(j) .gt. 0) &
                     write (*,"(T20,A,I4,A,ES15.5)")  "ORBITAL", CONTROL_instance%IONIZE_MO(j) ," OCCUPATION:",CONTROL_instance%MO_FRACTION_OCCUPATION(j)
             end do
          end if
       end do
    end if

    if(CONTROL_instance%POLARIZATION_ORDER > 1) then 

       write (*,"(T10,A,I5)") "POLARIZATION ORDER TO BE CALCULATED: ", CONTROL_instance%POLARIZATION_ORDER

    end if

    if(CONTROL_instance%FUKUI_FUNCTIONS) then 

       write (*,"(T10,A,I5)") "CALCULATING FUKUI FUNCTIONS"

    end if

    if(CONTROL_instance%EXCITE_SPECIES /= "NONE") then 

       write (*,"(T10,A,T10)") "CALCULATING", trim(CONTROL_instance%EXCITE_SPECIES) ,"IN THE FIRST EXCITED STATE"

    end if

    if(CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE) then 

       write (*,"(T10,A)") "BUILDING TWO PARTICLES MATRIX FOR ONE PARTICLE"

    end if

    if (CONTROL_instance%OPTIMIZE) then

       write (*,"(T10,A)") "GEOMETRY OPTIMIZATION:  T"
       write (*,"(T10,A,I5)") "OPTIMIZATION METHOD: ",CONTROL_instance%MINIMIZATION_METHOD
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
       write (*,"(T10,A,I5)") "EVOLUTION METHOD: ",CONTROL_instance%MINIMIZATION_METHOD
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
       write (*,"(T10,A)") "NONELECTRONIC DENSITY GUESS: "//trim(CONTROL_instance%SCF_NONELECTRONIC_TYPE_GUESS)
       write (*,"(T10,A)") "ELECTRONIC DENSITY GUESS: "//trim(CONTROL_instance%SCF_ELECTRONIC_TYPE_GUESS)

       write (*,"(T10,A)") "CRITERIUM OF CONVERGENCE: "//trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM)
       
       select case(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM)
       case("ENERGY")
          write (*,"(T10,A,E15.5)") "NONELECTRONIC ENERGY TOLERANCE IN SCFs: ",CONTROL_instance%NONELECTRONIC_ENERGY_TOLERANCE
          write (*,"(T10,A,E15.5)") "ELECTRONIC ENERGY TOLERANCE IN SCFs: ",CONTROL_instance%ELECTRONIC_ENERGY_TOLERANCE
          write (*,"(T10,A,E15.5)") "TOTAL ENERGY TOLERANCE IN SCFs: ",CONTROL_instance%TOTAL_ENERGY_TOLERANCE
          CONTROL_instance%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE=1.0
          CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE=1.0
          CONTROL_instance%TOTAL_DENSITY_MATRIX_TOLERANCE=1.0
       case("DENSITY")
          write (*,"(T10,A,E15.5)") "NONELECTRONIC DENSITY MATRIX TOLERANCE IN SCFs: ",CONTROL_instance%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE
          write (*,"(T10,A,E15.5)") "ELECTRONIC DENSITY MATRIX TOLERANCE IN SCFs: ",CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE
          write (*,"(T10,A,E15.5)") "TOTAL DENSITY TOLERANCE IN SCFs: ",CONTROL_instance%TOTAL_DENSITY_MATRIX_TOLERANCE
          CONTROL_instance%NONELECTRONIC_ENERGY_TOLERANCE=1.0
          CONTROL_instance%ELECTRONIC_ENERGY_TOLERANCE=1.0
          CONTROL_instance%TOTAL_ENERGY_TOLERANCE=1.0
       case ("BOTH")
          write (*,"(T10,A,E15.5)") "NONELECTRONIC ENERGY TOLERANCE IN SCFs: ",CONTROL_instance%NONELECTRONIC_ENERGY_TOLERANCE
          write (*,"(T10,A,E15.5)") "NONELECTRONIC DENSITY MATRIX TOLERANCE IN SCFs: ",CONTROL_instance%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE
          write (*,"(T10,A,E15.5)") "ELECTRONIC ENERGY TOLERANCE IN SCFs: ",CONTROL_instance%ELECTRONIC_ENERGY_TOLERANCE
          write (*,"(T10,A,E15.5)") "ELECTRONIC DENSITY MATRIX TOLERANCE IN SCFs: ",CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE
          write (*,"(T10,A,E15.5)") "TOTAL ENERGY TOLERANCE IN SCFs: ",CONTROL_instance%TOTAL_ENERGY_TOLERANCE
          write (*,"(T10,A,E15.5)") "TOTAL DENSITY TOLERANCE IN SCFs: ",CONTROL_instance%TOTAL_DENSITY_MATRIX_TOLERANCE
       case default
          call CONTROL_exception( ERROR, "unknown convergence criterium chosen", "at core program, CONTROL module")
       end select

       select case(CONTROL_instance%ITERATION_SCHEME)          
       case(0)
          write (*,"(T10,A)") "SCHEME OF ITERATION: NONELECTRONIC FULLY PER GLOBAL ITERATION"
          CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS=1
       case(1)
          write (*,"(T10,A)") "SCHEME OF ITERATION: ELECTRONIC FULLY PER GLOBAL ITERATION"
          CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS=1
       case(2)
          write (*,"(T10,A)") "SCHEME OF ITERATION: EACH SPECIES FULLY PER GLOBAL ITERATION"
       case(3)
          write (*,"(T10,A)") "SCHEME OF ITERATION: GLOBAL ITERATIONS"
          CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS=1
          CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS=1
       end select

       write (*,"(T10,A,I5)") "SCF MAX. SUBITERATIONS - NONELECTRONS : ",CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS
       write (*,"(T10,A,I5)") "SCF MAX. SUBITERATIONS - ELECTRONS : ",CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS
       write (*,"(T10,A,I5)") "SCF MAX. ITERATIONS - INTERSPECIES : ",CONTROL_instance%SCF_GLOBAL_MAX_ITERATIONS

       if (CONTROL_instance%NO_SCF) write (*,"(T10,A)") "COEFFICIENTS WILL BE READ AND NO SCF WILL BE PERFORMED"

       if (CONTROL_instance%FREEZE_NON_ELECTRONIC_ORBITALS) write (*,"(T10,A)") "Electrons will be frozen during SCF calculation"

       if ( CONTROL_instance%HARTREE_PRODUCT_GUESS) write (*,"(T10,A)") "HARTREE PRODUCT GUESS: T"



       write (*,"(T10,A)") "INTEGRAL STORAGE: "//trim(CONTROL_instance%INTEGRAL_STORAGE)
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

    if ( CONTROL_instance%ELECTRONIC_LEVEL_SHIFTING .ne. 0.0_8 .or. CONTROL_instance%NONELECTRONIC_LEVEL_SHIFTING .ne. 0.0_8 ) CONTROL_instance%ACTIVATE_LEVEL_SHIFTING=.true.
       
    if ( CONTROL_instance%ACTIVATE_LEVEL_SHIFTING .eqv. .true. ) then

       if ( CONTROL_instance%ELECTRONIC_LEVEL_SHIFTING .ne. 0.0_8 ) &
            write(*,"(T10,A,F10.6)") "SHIFTING ELECTRONIC VIRTUAL ORBITALS IN SCF BY:", CONTROL_instance%ELECTRONIC_LEVEL_SHIFTING

       if ( CONTROL_instance%NONELECTRONIC_LEVEL_SHIFTING .ne. 0.0_8 ) &
            write(*,"(T10,A,F10.6)") "SHIFTING NON-ELECTRONIC VIRTUAL ORBITALS IN SCF BY:", CONTROL_instance%NONELECTRONIC_LEVEL_SHIFTING

    end if

    if ( CONTROL_instance%EXCHANGE_ORBITALS_IN_SCF .eqv. .true. ) then
       write(*,"(T10,A)") "OCCUPIED ORBITALS IN SCF WILL BE REORDERED TO MAXIMIZE THE OVERLAP WITH RESPECT TO THE INITIAL GUESS"
       write(*,"(T10,A,F8.3)") "OVERLAP THRESHOLD FOR ORBITAL EXCHANGE:", CONTROL_instance%EXCHANGE_ORBITAL_THRESHOLD
    end if

    
    if(CONTROL_instance%LOCALIZE_ORBITALS) then
       select case (trim(CONTROL_instance%ERKALE_LOCALIZATION_METHOD))
       case("MU")
          write (*,"(T10,A)") "OCCUPIED ORBITALS WILL BE LOCALIZED WITH THE PIPEK-MEZEY SCHEME WITH MULLIKEN CHARGES"
       case("IAO")
          write (*,"(T10,A)") "OCCUPIED ORBITALS WILL BE LOCALIZED WITH THE PIPEK-MEZEY SCHEME WITH IBO CHARGES"
       case default
          write (*,"(T10,A)") "OCCUPIED ORBITALS WILL BE LOCALIZED WITH THE ", trim(CONTROL_instance%ERKALE_LOCALIZATION_METHOD) ," SCHEME"
       end select
    end if
       
         
    
    if(CONTROL_instance%SUBSYSTEM_EMBEDDING) then
       print *, "  "
       write (*,"(T10,A,A)") "TWO SCF CALCULATIONS WILL BE PERFORMED, FIRST FOR THE COMPLETE SYSTEM WITH ", CONTROL_instance%METHOD
       write (*,"(T10,A)") "THEN WITH HF FOR FRAGMENT ONE IN THE INPUT (SUBSYSTEM A)"
       write (*,"(T10,A)") "EMBEDDED IN THE POTENTIAL GENERATED BY THE OTHER FRAGMENTS (SUBSYSTEM B)."
       write (*,"(T10,A)") "POST-SCF CORRECTIONS ONLY WILL BE PERFORMED IN THE SECOND CALCULATION"
       print *, "  "

       write (*,"(T10,A,E8.1,A)") "SUBSYSTEM B IS BUILT FROM ORBITALS WITH POPULATION LOWER THAN ", CONTROL_instance%SUBSYSTEM_ORBITAL_THRESHOLD ," OVER FRAGMENT ONE ATOMS"
       if(CONTROL_instance%SUBSYSTEM_BASIS_THRESHOLD .gt. 0.0) &
            write (*,"(T10,A,E8.1,A)") "SHELLS WITH MULLIKEN POPULATION LOWER THAN ", CONTROL_instance%SUBSYSTEM_BASIS_THRESHOLD ," WILL BE REMOVED FROM SUBSYSTEM A BASIS SET"
            
    end if

    if(CONTROL_instance%SET_QDO_ENERGY_ZERO) then
       write (*, "(A)") "Setting the energy zero to the kinetic and potential energy of the free QDOs (3/2*omega)"
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
