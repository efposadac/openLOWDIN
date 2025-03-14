!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	PROF. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief This module is an interface to the libxc library for the evaluation of different electronic exchange correlation functionals. 
!! @author F. Moncada, 2017

!! Copyright (C) 2016 Micael Oliveira
!! All rights reserved.
!!
!! This file is dual-licensed under a GPL and a BSD license
!!
!! MPL License:
!!
!! This Source Code Form is subject to the terms of the Mozilla Public
!! License, v. 2.0. If a copy of the MPL was not distributed with this
!! file, You can obtain one at mozilla.org/MPL/2.0/.
!!
!! BSD License:
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!! 1. Redistributions of source code must retain the above copyright
!! notice, this list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above
!! copyright notice, this list of conditions and the following
!! disclaimer in the documentation and/or other materials provided
!! with the distribution.
!!
!! 3. Neither the name of the copyright holder nor the names of its
!! contributors may be used to endorse or promote products derived
!! from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
!! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
!! HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
!! STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
!! OF THE POSSIBILITY OF SUCH DAMAGE.

module LibxcInterface_

  use, intrinsic :: iso_c_binding
  implicit none

  private
  public :: &
                                ! version
       xc_f03_version, &
       xc_f03_version_string, &
                                ! func_info
       xc_f03_func_info_t, &
       xc_f03_func_info_get_number, &
       xc_f03_func_info_get_kind, &
       xc_f03_func_info_get_name, &
       xc_f03_func_info_get_family, &
       xc_f03_func_info_get_references, &
       xc_f03_func_info_get_flags, &
       xc_f03_func_info_get_n_ext_params, &
       xc_f03_func_info_get_ext_params_name, &
       xc_f03_func_info_get_ext_params_description, &
       xc_f03_func_info_get_ext_params_default_value, &
                                ! func_reference
       xc_f03_func_reference_t, &
       xc_f03_func_reference_get_ref, &
       xc_f03_func_reference_get_doi, &
       xc_f03_func_reference_get_bibtex, &
                                ! func
       xc_f03_func_t, &
       xc_f03_func_init, &
       xc_f03_func_end, &
       xc_f03_func_get_info, &
       xc_f03_functional_get_name, &
       xc_f03_functional_get_number, &
       xc_f03_family_from_id, &
       xc_f03_number_of_functionals, &
       xc_f03_maximum_name_length, &
       xc_f03_available_functional_numbers, &
       xc_f03_available_functional_names, &
       xc_f03_func_set_dens_threshold, &
       xc_f03_func_set_ext_params, &
                                ! lda
       xc_f03_lda, &
       xc_f03_lda_exc, &
       xc_f03_lda_exc_vxc, &
       xc_f03_lda_exc_vxc_fxc, &
       xc_f03_lda_exc_vxc_fxc_kxc, &
       xc_f03_lda_vxc, &
       xc_f03_lda_vxc_fxc, &
       xc_f03_lda_vxc_fxc_kxc, &
       xc_f03_lda_fxc, &
       xc_f03_lda_kxc, &
                                ! gga
       xc_f03_gga, &
       xc_f03_gga_exc, &
       xc_f03_gga_exc_vxc, &
       xc_f03_gga_exc_vxc_fxc, &
       xc_f03_gga_exc_vxc_fxc_kxc, &
       xc_f03_gga_vxc, &
       xc_f03_gga_vxc_fxc, &
       xc_f03_gga_vxc_fxc_kxc, &
       xc_f03_gga_fxc, &
       xc_f03_gga_kxc, &
       xc_f03_gga_ak13_get_asymptotic, &
       xc_f03_hyb_type, &
       xc_f03_hyb_exx_coef, &
       xc_f03_hyb_cam_coef, &
       xc_f03_nlc_coef, &
                                ! mgga
       xc_f03_mgga, &
       xc_f03_mgga_exc, &
       xc_f03_mgga_exc_vxc, &
       xc_f03_mgga_exc_vxc_fxc, &
       xc_f03_mgga_exc_vxc_fxc_kxc, &
       xc_f03_mgga_vxc, &
       xc_f03_mgga_vxc_fxc, &
       xc_f03_mgga_vxc_fxc_kxc, &
       xc_f03_mgga_fxc, &
       xc_f03_mgga_kxc

  integer(c_int), parameter, public :: &
       XC_UNPOLARIZED          =   1,     &  ! Spin unpolarized
       XC_POLARIZED            =   2         ! Spin polarized

  integer(c_int), parameter, public :: &
       XC_NON_RELATIVISTIC     =   0,     &  ! Functional includes or not relativistic
       XC_RELATIVISTIC         =   1         ! corrections. Only available in some functionals.

  ! Kinds
  integer(c_int), parameter, public :: &
       XC_EXCHANGE             =   0,     &
       XC_CORRELATION          =   1,     &
       XC_EXCHANGE_CORRELATION =   2,     &
       XC_KINETIC              =   3

  ! Different types of hybrid functionals.
  integer(c_int), parameter, public :: &
       XC_HYB_SEMILOCAL      =  0, &  !/* Standard semi-local functional (not a hybrid) */
       XC_HYB_HYBRID         =  1, &  !/* Standard hybrid functional */
       XC_HYB_CAM            =  2, &  !/* Coulomb attenuated hybrid */
       XC_HYB_CAMY           =  3, &  !/* Coulomb attenuated hybrid with a Yukawa screening */
       XC_HYB_CAMG           =  4, &  !/* Coulomb attenuated hybrid with a Gaussian screening */
       XC_HYB_DOUBLE_HYBRID  =  5, &  !/* Double hybrid */
       XC_HYB_MIXTURE       = 32768  !/* More complicated mixture (have to check individual terms) */
  
  ! Families of xc functionals
  integer(c_int), parameter, public :: &
       XC_FAMILY_UNKNOWN       =  -1,     &
       XC_FAMILY_NONE          =   0,     &
       XC_FAMILY_LDA           =   1,     &
       XC_FAMILY_GGA           =   2,     &
       XC_FAMILY_MGGA          =   4,     &
       XC_FAMILY_LCA           =   8,     &
       XC_FAMILY_OEP           =  16

  integer(c_int), parameter, public :: &
       XC_FLAGS_HAVE_EXC        =     1,   &
       XC_FLAGS_HAVE_VXC        =     2,   &
       XC_FLAGS_HAVE_FXC        =     4,   &
       XC_FLAGS_HAVE_KXC        =     8,   &
       XC_FLAGS_HAVE_LXC        =    16,   &
       XC_FLAGS_HAVE_ALL        =    31,   & ! The most common case
       XC_FLAGS_1D              =    32,   &
       XC_FLAGS_2D              =    64,   &
       XC_FLAGS_3D              =   128,   &
       XC_FLAGS_VV10            =  1024,   &
       XC_FLAGS_STABLE          =  8192,   &
       XC_FLAGS_DEVELOPMENT     = 16384,   &
       XC_FLAGS_NEEDS_LAPLACIAN = 32768

  integer(c_int), parameter, public :: &
       XC_TAU_EXPLICIT         =     0,   &
       XC_TAU_EXPANSION        =     1

  integer(c_int), parameter, public :: &
       XC_MAX_REFERENCES       =     5

  ! List of functionals
  !#include "libxc_inc.f90"
  integer(c_int), parameter, public :: XC_LDA_X                       =   1  ! Exchange
  integer(c_int), parameter, public :: XC_LDA_C_WIGNER                =   2  ! Wigner parametrization
  integer(c_int), parameter, public :: XC_LDA_C_RPA                   =   3  ! Random Phase Approximation
  integer(c_int), parameter, public :: XC_LDA_C_HL                    =   4  ! Hedin & Lundqvist
  integer(c_int), parameter, public :: XC_LDA_C_GL                    =   5  !  Gunnarson & Lundqvist
  integer(c_int), parameter, public :: XC_LDA_C_XALPHA                =   6  !  Slater Xalpha
  integer(c_int), parameter, public :: XC_LDA_C_VWN                   =   7  ! Vosko, Wilk, & Nusair (5)
  integer(c_int), parameter, public :: XC_LDA_C_VWN_RPA               =   8  ! Vosko, Wilk, & Nusair (RPA)
  integer(c_int), parameter, public :: XC_LDA_C_PZ                    =   9  ! Perdew & Zunger
  integer(c_int), parameter, public :: XC_LDA_C_PZ_MOD                =  10  !  Perdew & Zunger (Modified)
  integer(c_int), parameter, public :: XC_LDA_C_OB_PZ                 =  11  !  Ortiz & Ballone (PZ)
  integer(c_int), parameter, public :: XC_LDA_C_PW                    =  12  ! Perdew & Wang
  integer(c_int), parameter, public :: XC_LDA_C_PW_MOD                =  13  !  Perdew & Wang (Modified)
  integer(c_int), parameter, public :: XC_LDA_C_OB_PW                 =  14  !  Ortiz & Ballone (PW)
  integer(c_int), parameter, public :: XC_LDA_C_2D_AMGB               =  15  ! Attaccalite et al
  integer(c_int), parameter, public :: XC_LDA_C_2D_PRM                =  16  ! Pittalis, Rasanen & Marques correlation in 2D
  integer(c_int), parameter, public :: XC_LDA_C_VBH                   =  17  !  von Barth & Hedin
  integer(c_int), parameter, public :: XC_LDA_C_1D_CSC                =  18  ! Casula, Sorella, and Senatore 1D correlation
  integer(c_int), parameter, public :: XC_LDA_X_2D                    =  19  ! Exchange in 2D
  integer(c_int), parameter, public :: XC_LDA_XC_TETER93              =  20  ! Teter 93 parametrization
  integer(c_int), parameter, public :: XC_LDA_X_1D_SOFT               =  21  ! Exchange in 1D for a soft-Coulomb interaction
  integer(c_int), parameter, public :: XC_LDA_C_ML1                   =  22  ! Modified LSD (version 1) of Proynov and Salahub
  integer(c_int), parameter, public :: XC_LDA_C_ML2                   =  23  !  Modified LSD (version 2) of Proynov and Salahub
  integer(c_int), parameter, public :: XC_LDA_C_GOMBAS                =  24  ! Gombas parametrization
  integer(c_int), parameter, public :: XC_LDA_C_PW_RPA                =  25  !  Perdew & Wang fit of the RPA
  integer(c_int), parameter, public :: XC_LDA_C_1D_LOOS               =  26  ! P-F Loos correlation LDA
  integer(c_int), parameter, public :: XC_LDA_C_RC04                  =  27  ! Ragot-Cortona
  integer(c_int), parameter, public :: XC_LDA_C_VWN_1                 =  28  ! Vosko, Wilk, & Nusair (1)
  integer(c_int), parameter, public :: XC_LDA_C_VWN_2                 =  29  ! Vosko, Wilk, & Nusair (2)
  integer(c_int), parameter, public :: XC_LDA_C_VWN_3                 =  30  ! Vosko, Wilk, & Nusair (3)
  integer(c_int), parameter, public :: XC_LDA_C_VWN_4                 =  31  ! Vosko, Wilk, & Nusair (4)
  integer(c_int), parameter, public :: XC_LDA_XC_ZLP                  =  43  ! Zhao, Levy & Parr, Eq. (20)
  integer(c_int), parameter, public :: XC_LDA_K_TF                    =  50  ! Thomas-Fermi kinetic energy functional
  integer(c_int), parameter, public :: XC_LDA_K_LP                    =  51  !  Lee and Parr Gaussian ansatz
  integer(c_int), parameter, public :: XC_HYB_LDA_XC_LDA0             = 177  !  LDA0: hybrid LDA exchange
  integer(c_int), parameter, public :: XC_HYB_LDA_XC_CAM_LDA0         = 178  ! CAM version of LDA0
  integer(c_int), parameter, public :: XC_LDA_XC_KSDT                 = 259  ! Karasiev et al. parametrization
  integer(c_int), parameter, public :: XC_LDA_C_CHACHIYO              = 287  ! Chachiyo simple 2 parameter correlation
  integer(c_int), parameter, public :: XC_LDA_C_LP96                  = 289  ! Liu-Parr correlation
  integer(c_int), parameter, public :: XC_LDA_C_CHACHIYO_MOD          = 307  ! Chachiyo simple 2 parameter correlation with modified scaling
  integer(c_int), parameter, public :: XC_LDA_C_KARASIEV_MOD          = 308  !  Karasiev reparameterization of Chachiyo with modified scaling
  integer(c_int), parameter, public :: XC_LDA_X_REL                   = 532  ! Relativistic exchange
  integer(c_int), parameter, public :: XC_LDA_XC_1D_EHWLRG_1          = 536  ! LDA constructed from slab-like systems of 1 electron
  integer(c_int), parameter, public :: XC_LDA_XC_1D_EHWLRG_2          = 537  !  LDA constructed from slab-like systems of 2 electrons
  integer(c_int), parameter, public :: XC_LDA_XC_1D_EHWLRG_3          = 538  !  LDA constructed from slab-like systems of 3 electrons
  integer(c_int), parameter, public :: XC_LDA_X_ERF                   = 546  ! Attenuated exchange LDA (erf)
  integer(c_int), parameter, public :: XC_LDA_XC_LP_A                 = 547  !  Lee-Parr reparametrization B
  integer(c_int), parameter, public :: XC_LDA_XC_LP_B                 = 548  !  Lee-Parr reparametrization B
  integer(c_int), parameter, public :: XC_LDA_X_RAE                   = 549  !  Rae self-energy corrected exchange
  integer(c_int), parameter, public :: XC_LDA_K_ZLP                   = 550  ! kinetic energy version of ZLP
  integer(c_int), parameter, public :: XC_LDA_C_MCWEENY               = 551  !  McWeeny 76
  integer(c_int), parameter, public :: XC_LDA_C_BR78                  = 552  !  Brual & Rothstein 78
  integer(c_int), parameter, public :: XC_LDA_C_PK09                  = 554  ! Proynov and Kong 2009
  integer(c_int), parameter, public :: XC_LDA_C_OW_LYP                = 573  !  Wigner with corresponding LYP parameters
  integer(c_int), parameter, public :: XC_LDA_C_OW                    = 574  !  Optimized Wigner
  integer(c_int), parameter, public :: XC_LDA_XC_GDSMFB               = 577  !  Groth et al. parametrization
  integer(c_int), parameter, public :: XC_LDA_C_GK72                  = 578  ! Gordon and Kim 1972
  integer(c_int), parameter, public :: XC_LDA_C_KARASIEV              = 579  !  Karasiev reparameterization of Chachiyo
  integer(c_int), parameter, public :: XC_LDA_K_LP96                  = 580  !  Liu-Parr kinetic
  integer(c_int), parameter, public :: XC_LDA_XC_BN05                 = 588  ! Baer and Neuhauser, gamma=1
  integer(c_int), parameter, public :: XC_LDA_C_PMGB06                = 590  ! Long-range LDA correlation functional
  integer(c_int), parameter, public :: XC_LDA_XC_TIH                  = 599  ! Neural network LDA from Tozer et al
  integer(c_int), parameter, public :: XC_LDA_X_1D_EXPONENTIAL        = 600  ! Exchange in 1D for an exponentially screened interaction
  integer(c_int), parameter, public :: XC_LDA_C_UPW92                 = 683  !  Ruggeri, Rios, and Alavi unrestricted fit
  integer(c_int), parameter, public :: XC_LDA_C_RPW92                 = 684  !  Ruggeri, Rios, and Alavi restricted fit
  integer(c_int), parameter, public :: XC_LDA_X_SLOC                  = 692  ! simple local model for Slater potential
  integer(c_int), parameter, public :: XC_GGA_X_GAM                   =  32  !  GAM functional from Minnesota
  integer(c_int), parameter, public :: XC_GGA_C_GAM                   =  33  !  GAM functional from Minnesota
  integer(c_int), parameter, public :: XC_GGA_X_HCTH_A                =  34  ! HCTH-A
  integer(c_int), parameter, public :: XC_GGA_X_EV93                  =  35  ! Engel and Vosko
  integer(c_int), parameter, public :: XC_GGA_X_BCGP                  =  38  !  Burke, Cancio, Gould, and Pittalis
  integer(c_int), parameter, public :: XC_GGA_C_ACGGA                 =  39  ! acGGA, asymptotically corrected GGA
  integer(c_int), parameter, public :: XC_GGA_X_LAMBDA_OC2_N          =  40  !  lambda_OC2(N) version of PBE
  integer(c_int), parameter, public :: XC_GGA_X_B86_R                 =  41  !  Revised Becke 86 Xalpha,beta,gamma (with mod. grad. correction)
  integer(c_int), parameter, public :: XC_GGA_X_LAMBDA_CH_N           =  44  !  lambda_CH(N) version of PBE
  integer(c_int), parameter, public :: XC_GGA_X_LAMBDA_LO_N           =  45  !  lambda_LO(N) version of PBE
  integer(c_int), parameter, public :: XC_GGA_X_HJS_B88_V2            =  46  ! HJS screened exchange corrected B88 version
  integer(c_int), parameter, public :: XC_GGA_C_Q2D                   =  47  ! Chiodo et al
  integer(c_int), parameter, public :: XC_GGA_X_Q2D                   =  48  ! Chiodo et al
  integer(c_int), parameter, public :: XC_GGA_X_PBE_MOL               =  49  !  Del Campo, Gazquez, Trickey and Vela (PBE-like)
  integer(c_int), parameter, public :: XC_GGA_K_TFVW                  =  52  ! Thomas-Fermi plus von Weiszaecker correction
  integer(c_int), parameter, public :: XC_GGA_K_REVAPBEINT            =  53  !  interpolated version of REVAPBE
  integer(c_int), parameter, public :: XC_GGA_K_APBEINT               =  54  ! interpolated version of APBE
  integer(c_int), parameter, public :: XC_GGA_K_REVAPBE               =  55  !  revised APBE
  integer(c_int), parameter, public :: XC_GGA_X_AK13                  =  56  ! Armiento & Kuemmel 2013
  integer(c_int), parameter, public :: XC_GGA_K_MEYER                 =  57  ! Meyer,  Wang, and Young
  integer(c_int), parameter, public :: XC_GGA_X_LV_RPW86              =  58  ! Berland and Hyldgaard
  integer(c_int), parameter, public :: XC_GGA_X_PBE_TCA               =  59  !  PBE revised by Tognetti et al
  integer(c_int), parameter, public :: XC_GGA_X_PBEINT                =  60  ! PBE for hybrid interfaces
  integer(c_int), parameter, public :: XC_GGA_C_ZPBEINT               =  61  ! spin-dependent gradient correction to PBEint
  integer(c_int), parameter, public :: XC_GGA_C_PBEINT                =  62  !  PBE for hybrid interfaces
  integer(c_int), parameter, public :: XC_GGA_C_ZPBESOL               =  63  !  spin-dependent gradient correction to PBEsol
  integer(c_int), parameter, public :: XC_GGA_XC_OPBE_D               =  65  !  oPBE_D functional of Goerigk and Grimme
  integer(c_int), parameter, public :: XC_GGA_XC_OPWLYP_D             =  66  !  oPWLYP-D functional of Goerigk and Grimme
  integer(c_int), parameter, public :: XC_GGA_XC_OBLYP_D              =  67  ! oBLYP-D functional of Goerigk and Grimme
  integer(c_int), parameter, public :: XC_GGA_X_VMT84_GE              =  68  !  VMT{8,4} with constraint satisfaction with mu = mu_GE
  integer(c_int), parameter, public :: XC_GGA_X_VMT84_PBE             =  69  ! VMT{8,4} with constraint satisfaction with mu = mu_PBE
  integer(c_int), parameter, public :: XC_GGA_X_VMT_GE                =  70  !  Vela, Medel, and Trickey with mu = mu_GE
  integer(c_int), parameter, public :: XC_GGA_X_VMT_PBE               =  71  ! Vela, Medel, and Trickey with mu = mu_PBE
  integer(c_int), parameter, public :: XC_GGA_C_N12_SX                =  79  !  N12-SX functional from Minnesota
  integer(c_int), parameter, public :: XC_GGA_C_N12                   =  80  ! N12 functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_GGA_X_N12_SX            =  81  !  N12-SX functional from Minnesota
  integer(c_int), parameter, public :: XC_GGA_X_N12                   =  82  ! N12 functional from Minnesota
  integer(c_int), parameter, public :: XC_GGA_C_REGTPSS               =  83  ! Regularized TPSS correlation (ex-VPBE)
  integer(c_int), parameter, public :: XC_GGA_C_OP_XALPHA             =  84  ! one-parameter progressive functional (XALPHA version)
  integer(c_int), parameter, public :: XC_GGA_C_OP_G96                =  85  ! one-parameter progressive functional (G96 version)
  integer(c_int), parameter, public :: XC_GGA_C_OP_PBE                =  86  ! one-parameter progressive functional (PBE version)
  integer(c_int), parameter, public :: XC_GGA_C_OP_B88                =  87  ! one-parameter progressive functional (B88 version)
  integer(c_int), parameter, public :: XC_GGA_C_FT97                  =  88  ! Filatov & Thiel correlation
  integer(c_int), parameter, public :: XC_GGA_C_SPBE                  =  89  !  PBE correlation to be used with the SSB exchange
  integer(c_int), parameter, public :: XC_GGA_X_SSB_SW                =  90  ! Swart, Sola and Bickelhaupt correction to PBE
  integer(c_int), parameter, public :: XC_GGA_X_SSB                   =  91  !  Swart, Sola and Bickelhaupt
  integer(c_int), parameter, public :: XC_GGA_X_SSB_D                 =  92  !  Swart, Sola and Bickelhaupt dispersion
  integer(c_int), parameter, public :: XC_GGA_XC_HCTH_407P            =  93  !  HCTH/407+
  integer(c_int), parameter, public :: XC_GGA_XC_HCTH_P76             =  94  !  HCTH p=7/6
  integer(c_int), parameter, public :: XC_GGA_XC_HCTH_P14             =  95  !  HCTH p=1/4
  integer(c_int), parameter, public :: XC_GGA_XC_B97_GGA1             =  96  !  Becke 97 GGA-1
  integer(c_int), parameter, public :: XC_GGA_C_HCTH_A                =  97  ! HCTH-A
  integer(c_int), parameter, public :: XC_GGA_X_BPCCAC                =  98  ! BPCCAC (GRAC for the energy)
  integer(c_int), parameter, public :: XC_GGA_C_REVTCA                =  99  ! Tognetti, Cortona, Adamo (revised)
  integer(c_int), parameter, public :: XC_GGA_C_TCA                   = 100  ! Tognetti, Cortona, Adamo
  integer(c_int), parameter, public :: XC_GGA_X_PBE                   = 101  ! Perdew, Burke & Ernzerhof exchange
  integer(c_int), parameter, public :: XC_GGA_X_PBE_R                 = 102  !  Perdew, Burke & Ernzerhof exchange (revised)
  integer(c_int), parameter, public :: XC_GGA_X_B86                   = 103  ! Becke 86 Xalpha,beta,gamma
  integer(c_int), parameter, public :: XC_GGA_X_HERMAN                = 104  ! Herman et al original GGA
  integer(c_int), parameter, public :: XC_GGA_X_B86_MGC               = 105  !  Becke 86 Xalpha,beta,gamma (with mod. grad. correction)
  integer(c_int), parameter, public :: XC_GGA_X_B88                   = 106  ! Becke 88
  integer(c_int), parameter, public :: XC_GGA_X_G96                   = 107  ! Gill 96
  integer(c_int), parameter, public :: XC_GGA_X_PW86                  = 108  ! Perdew & Wang 86
  integer(c_int), parameter, public :: XC_GGA_X_PW91                  = 109  ! Perdew & Wang 91
  integer(c_int), parameter, public :: XC_GGA_X_OPTX                  = 110  ! Handy & Cohen OPTX 01
  integer(c_int), parameter, public :: XC_GGA_X_DK87_R1               = 111  ! dePristo & Kress 87 (version R1)
  integer(c_int), parameter, public :: XC_GGA_X_DK87_R2               = 112  !  dePristo & Kress 87 (version R2)
  integer(c_int), parameter, public :: XC_GGA_X_LG93                  = 113  ! Lacks & Gordon 93
  integer(c_int), parameter, public :: XC_GGA_X_FT97_A                = 114  ! Filatov & Thiel 97 (version A)
  integer(c_int), parameter, public :: XC_GGA_X_FT97_B                = 115  !  Filatov & Thiel 97 (version B)
  integer(c_int), parameter, public :: XC_GGA_X_PBE_SOL               = 116  !  Perdew, Burke & Ernzerhof exchange (solids)
  integer(c_int), parameter, public :: XC_GGA_X_RPBE                  = 117  ! Hammer, Hansen & Norskov (PBE-like)
  integer(c_int), parameter, public :: XC_GGA_X_WC                    = 118  ! Wu & Cohen
  integer(c_int), parameter, public :: XC_GGA_X_MPW91                 = 119  !  Modified form of PW91 by Adamo & Barone
  integer(c_int), parameter, public :: XC_GGA_X_AM05                  = 120  ! Armiento & Mattsson 05 exchange
  integer(c_int), parameter, public :: XC_GGA_X_PBEA                  = 121  ! Madsen (PBE-like)
  integer(c_int), parameter, public :: XC_GGA_X_MPBE                  = 122  ! Adamo & Barone modification to PBE
  integer(c_int), parameter, public :: XC_GGA_X_XPBE                  = 123  !  xPBE reparametrization by Xu & Goddard
  integer(c_int), parameter, public :: XC_GGA_X_2D_B86_MGC            = 124  ! Becke 86 MGC for 2D systems
  integer(c_int), parameter, public :: XC_GGA_X_BAYESIAN              = 125  ! Bayesian best fit for the enhancement factor
  integer(c_int), parameter, public :: XC_GGA_X_PBE_JSJR              = 126  !  JSJR reparametrization by Pedroza, Silva & Capelle
  integer(c_int), parameter, public :: XC_GGA_X_2D_B88                = 127  ! Becke 88 in 2D
  integer(c_int), parameter, public :: XC_GGA_X_2D_B86                = 128  ! Becke 86 Xalpha, beta, gamma
  integer(c_int), parameter, public :: XC_GGA_X_2D_PBE                = 129  ! Perdew, Burke & Ernzerhof exchange in 2D
  integer(c_int), parameter, public :: XC_GGA_C_PBE                   = 130  ! Perdew, Burke & Ernzerhof correlation
  integer(c_int), parameter, public :: XC_GGA_C_LYP                   = 131  ! Lee, Yang & Parr
  integer(c_int), parameter, public :: XC_GGA_C_P86                   = 132  ! Perdew 86
  integer(c_int), parameter, public :: XC_GGA_C_PBE_SOL               = 133  !  Perdew, Burke & Ernzerhof correlation SOL
  integer(c_int), parameter, public :: XC_GGA_C_PW91                  = 134  ! Perdew & Wang 91
  integer(c_int), parameter, public :: XC_GGA_C_AM05                  = 135  ! Armiento & Mattsson 05 correlation
  integer(c_int), parameter, public :: XC_GGA_C_XPBE                  = 136  !  xPBE reparametrization by Xu & Goddard
  integer(c_int), parameter, public :: XC_GGA_C_LM                    = 137  ! Langreth and Mehl correlation
  integer(c_int), parameter, public :: XC_GGA_C_PBE_JRGX              = 138  !  JRGX reparametrization by Pedroza, Silva & Capelle
  integer(c_int), parameter, public :: XC_GGA_X_OPTB88_VDW            = 139  !  Becke 88 reoptimized to be used with vdW functional of Dion et al
  integer(c_int), parameter, public :: XC_GGA_X_PBEK1_VDW             = 140  !  PBE reparametrization for vdW
  integer(c_int), parameter, public :: XC_GGA_X_OPTPBE_VDW            = 141  !  PBE reparametrization for vdW
  integer(c_int), parameter, public :: XC_GGA_X_RGE2                  = 142  ! Regularized PBE
  integer(c_int), parameter, public :: XC_GGA_C_RGE2                  = 143  !  Regularized PBE
  integer(c_int), parameter, public :: XC_GGA_X_RPW86                 = 144  !  refitted Perdew & Wang 86
  integer(c_int), parameter, public :: XC_GGA_X_KT1                   = 145  ! Exchange part of Keal and Tozer version 1
  integer(c_int), parameter, public :: XC_GGA_XC_KT2                  = 146  !  Keal and Tozer version 2
  integer(c_int), parameter, public :: XC_GGA_C_WL                    = 147  ! Wilson & Levy
  integer(c_int), parameter, public :: XC_GGA_C_WI                    = 148  !  Wilson & Ivanov
  integer(c_int), parameter, public :: XC_GGA_X_MB88                  = 149  !  Modified Becke 88 for proton transfer
  integer(c_int), parameter, public :: XC_GGA_X_SOGGA                 = 150  !  Second-order generalized gradient approximation
  integer(c_int), parameter, public :: XC_GGA_X_SOGGA11               = 151  ! Second-order generalized gradient approximation 2011
  integer(c_int), parameter, public :: XC_GGA_C_SOGGA11               = 152  ! SOGGA11 correlation
  integer(c_int), parameter, public :: XC_GGA_C_WI0                   = 153  ! Wilson & Ivanov initial version
  integer(c_int), parameter, public :: XC_GGA_XC_TH1                  = 154  !  Tozer and Handy v. 1
  integer(c_int), parameter, public :: XC_GGA_XC_TH2                  = 155  ! Tozer and Handy v. 2
  integer(c_int), parameter, public :: XC_GGA_XC_TH3                  = 156  ! Tozer and Handy v. 3
  integer(c_int), parameter, public :: XC_GGA_XC_TH4                  = 157  !  Tozer and Handy v. 4
  integer(c_int), parameter, public :: XC_GGA_X_C09X                  = 158  ! C09x to be used with the VdW of Rutgers-Chalmers
  integer(c_int), parameter, public :: XC_GGA_C_SOGGA11_X             = 159  !  SOGGA11-X correlation
  integer(c_int), parameter, public :: XC_GGA_X_LB                    = 160  ! van Leeuwen & Baerends
  integer(c_int), parameter, public :: XC_GGA_XC_HCTH_93              = 161  !  HCTH functional fitted to  93 molecules
  integer(c_int), parameter, public :: XC_GGA_XC_HCTH_120             = 162  !  HCTH functional fitted to 120 molecules
  integer(c_int), parameter, public :: XC_GGA_XC_HCTH_147             = 163  !  HCTH functional fitted to 147 molecules
  integer(c_int), parameter, public :: XC_GGA_XC_HCTH_407             = 164  !  HCTH functional fitted to 407 molecules
  integer(c_int), parameter, public :: XC_GGA_XC_EDF1                 = 165  ! Empirical functionals from Adamson, Gill, and Pople
  integer(c_int), parameter, public :: XC_GGA_XC_XLYP                 = 166  ! XLYP functional
  integer(c_int), parameter, public :: XC_GGA_XC_KT1                  = 167  !  Keal and Tozer version 1
  integer(c_int), parameter, public :: XC_GGA_X_LSPBE                 = 168  ! PW91-like exchange with simple analytical form
  integer(c_int), parameter, public :: XC_GGA_X_LSRPBE                = 169  ! PW91-like modification of RPBE
  integer(c_int), parameter, public :: XC_GGA_XC_B97_D                = 170  !  Grimme functional to be used with C6 vdW term
  integer(c_int), parameter, public :: XC_GGA_X_OPTB86B_VDW           = 171  !  Becke 86 reoptimized for use with vdW functional of Dion et al
  integer(c_int), parameter, public :: XC_GGA_XC_PBE1W                = 173  !  Functionals fitted for water
  integer(c_int), parameter, public :: XC_GGA_XC_MPWLYP1W             = 174  !  Functionals fitted for water
  integer(c_int), parameter, public :: XC_GGA_XC_PBELYP1W             = 175  !  Functionals fitted for water
  integer(c_int), parameter, public :: XC_GGA_C_ACGGAP                = 176  ! Asymptotically corrected GGA +
  integer(c_int), parameter, public :: XC_GGA_X_B88_6311G             = 179  !  Becke 88 reoptimized with 6-311G** basis set
  integer(c_int), parameter, public :: XC_GGA_X_NCAP                  = 180  ! Nearly correct asymptotic potential
  integer(c_int), parameter, public :: XC_GGA_XC_NCAP                 = 181  !  Nearly correct asymptotic potential + P86 correlation
  integer(c_int), parameter, public :: XC_GGA_X_LBM                   = 182  !  van Leeuwen & Baerends modified
  integer(c_int), parameter, public :: XC_GGA_X_OL2                   = 183  ! Exchange form based on Ou-Yang and Levy v.2
  integer(c_int), parameter, public :: XC_GGA_X_APBE                  = 184  !  mu fixed from the semiclassical neutral atom
  integer(c_int), parameter, public :: XC_GGA_K_APBE                  = 185  ! mu fixed from the semiclassical neutral atom
  integer(c_int), parameter, public :: XC_GGA_C_APBE                  = 186  !  mu fixed from the semiclassical neutral atom
  integer(c_int), parameter, public :: XC_GGA_K_TW1                   = 187  !  Tran and Wesolowski set 1 (Table II)
  integer(c_int), parameter, public :: XC_GGA_K_TW2                   = 188  !  Tran and Wesolowski set 2 (Table II)
  integer(c_int), parameter, public :: XC_GGA_K_TW3                   = 189  !  Tran and Wesolowski set 3 (Table II)
  integer(c_int), parameter, public :: XC_GGA_K_TW4                   = 190  !  Tran and Wesolowski set 4 (Table II)
  integer(c_int), parameter, public :: XC_GGA_X_HTBS                  = 191  ! Haas, Tran, Blaha, and Schwarz
  integer(c_int), parameter, public :: XC_GGA_X_AIRY                  = 192  ! Constantin et al based on the Airy gas
  integer(c_int), parameter, public :: XC_GGA_X_LAG                   = 193  ! Local Airy Gas
  integer(c_int), parameter, public :: XC_GGA_XC_MOHLYP               = 194  !  Functional for organometallic chemistry
  integer(c_int), parameter, public :: XC_GGA_XC_MOHLYP2              = 195  !  Functional for barrier heights
  integer(c_int), parameter, public :: XC_GGA_XC_TH_FL                = 196  ! Tozer and Handy v. FL
  integer(c_int), parameter, public :: XC_GGA_XC_TH_FC                = 197  !  Tozer and Handy v. FC
  integer(c_int), parameter, public :: XC_GGA_XC_TH_FCFO              = 198  !  Tozer and Handy v. FCFO
  integer(c_int), parameter, public :: XC_GGA_XC_TH_FCO               = 199  !  Tozer and Handy v. FCO
  integer(c_int), parameter, public :: XC_GGA_C_OPTC                  = 200  ! Optimized correlation functional of Cohen and Handy
  integer(c_int), parameter, public :: XC_GGA_X_ECMV92                = 215  !  Engel, Chevary, Macdonald, and Vosko
  integer(c_int), parameter, public :: XC_GGA_C_PBE_VWN               = 216  ! Perdew, Burke & Ernzerhof correlation based on VWN LDA
  integer(c_int), parameter, public :: XC_GGA_C_P86_FT                = 217  !  Perdew 86 with a more accurate value for ftilde
  integer(c_int), parameter, public :: XC_GGA_C_PBELOC                = 246  ! Semilocal dynamical correlation
  integer(c_int), parameter, public :: XC_GGA_XC_VV10                 = 255  ! Vydrov and Van Voorhis
  integer(c_int), parameter, public :: XC_GGA_C_PBEFE                 = 258  !  PBE for formation energies
  integer(c_int), parameter, public :: XC_GGA_C_OP_PW91               = 262  ! one-parameter progressive functional (PW91 version)
  integer(c_int), parameter, public :: XC_GGA_X_PBEFE                 = 265  !  PBE for formation energies
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B97_1P           = 266  !  version of B97 by Cohen and Handy
  integer(c_int), parameter, public :: XC_GGA_X_CAP                   = 270  ! Correct Asymptotic Potential
  integer(c_int), parameter, public :: XC_GGA_X_EB88                  = 271  !  Non-empirical (excogitated) B88 functional of Becke and Elliott
  integer(c_int), parameter, public :: XC_GGA_C_PBE_MOL               = 272  !  Del Campo, Gazquez, Trickey and Vela (PBE-like)
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_PBE_MOL0         = 273  !  PBEmol0
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_PBE_SOL0         = 274  !  PBEsol0
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_PBEB0            = 275  !  PBEbeta0
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_PBE_MOLB0        = 276  !  PBEmolbeta0
  integer(c_int), parameter, public :: XC_GGA_K_ABSP3                 = 277  !  gamma-TFvW form by Acharya et al [g = 1 - 1.513/N^0.35]
  integer(c_int), parameter, public :: XC_GGA_K_ABSP4                 = 278  !  gamma-TFvW form by Acharya et al [g = l = 1/(1 + 1.332/N^(1/3))]
  integer(c_int), parameter, public :: XC_GGA_C_BMK                   = 280  !  Boese-Martin for kinetics
  integer(c_int), parameter, public :: XC_GGA_C_TAU_HCTH              = 281  !  correlation part of tau-hcth
  integer(c_int), parameter, public :: XC_GGA_C_HYB_TAU_HCTH          = 283  !  correlation part of hyb_tau-hcth
  integer(c_int), parameter, public :: XC_GGA_X_BEEFVDW               = 285  ! BEEF-vdW exchange
  integer(c_int), parameter, public :: XC_GGA_XC_BEEFVDW              = 286  !  BEEF-vdW exchange-correlation
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_PBE50            = 290  !  PBE0 with 50% exx
  integer(c_int), parameter, public :: XC_GGA_X_PBETRANS              = 291  ! Gradient-based interpolation between PBE and revPBE
  integer(c_int), parameter, public :: XC_GGA_X_CHACHIYO              = 298  ! Chachiyo exchange
  integer(c_int), parameter, public :: XC_GGA_C_CHACHIYO              = 309  ! Chachiyo simple GGA correlation
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LC_BLYP          = 400  ! Long-range corrected BLYP
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B3PW91           = 401  ! The original (ACM) hybrid of Becke
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B3LYP            = 402  !  The (in)famous B3LYP
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B3P86            = 403  !  Perdew 86 hybrid similar to B3PW91
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_O3LYP            = 404  ! hybrid using the optx functional
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_MPW1K            = 405  !  mixture of mPW91 and PW91 optimized for kinetics
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_PBEH             = 406  ! aka PBE0 or PBE1PBE
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B97              = 407  ! Becke 97
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B97_1            = 408  !  Becke 97-1
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_APF              = 409  !  APF hybrid density functional
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B97_2            = 410  !  Becke 97-2
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_X3LYP            = 411  !  hybrid by Xu and Goddard
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B1WC             = 412  ! Becke 1-parameter mixture of WC and PBE
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B97_K            = 413  !  Boese-Martin for Kinetics
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B97_3            = 414  !  Becke 97-3
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_MPW3PW           = 415  !  mixture with the mPW functional
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B1LYP            = 416  !  Becke 1-parameter mixture of B88 and LYP
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B1PW91           = 417  !  Becke 1-parameter mixture of B88 and PW91
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_MPW1PW           = 418  !  Becke 1-parameter mixture of mPW91 and PW91
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_MPW3LYP          = 419  !  mixture of mPW and LYP
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_SB98_1A          = 420  !  Schmider-Becke 98 parameterization 1a
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_SB98_1B          = 421  !  Schmider-Becke 98 parameterization 1b
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_SB98_1C          = 422  !  Schmider-Becke 98 parameterization 1c
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_SB98_2A          = 423  !  Schmider-Becke 98 parameterization 2a
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_SB98_2B          = 424  !  Schmider-Becke 98 parameterization 2b
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_SB98_2C          = 425  !  Schmider-Becke 98 parameterization 2c
  integer(c_int), parameter, public :: XC_HYB_GGA_X_SOGGA11_X         = 426  !  Hybrid based on SOGGA11 form
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HSE03            = 427  ! the 2003 version of the screened hybrid HSE
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HSE06            = 428  !  the 2006 version of the screened hybrid HSE
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HJS_PBE          = 429  !  HJS hybrid screened exchange PBE version
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HJS_PBE_SOL      = 430  !  HJS hybrid screened exchange PBE_SOL version
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HJS_B88          = 431  !  HJS hybrid screened exchange B88 version
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HJS_B97X         = 432  !  HJS hybrid screened exchange B97x version
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_CAM_B3LYP        = 433  ! CAM version of B3LYP
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_TUNED_CAM_B3LYP  = 434  !  CAM version of B3LYP tuned for excitations
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_BHANDH           = 435  !  Becke half-and-half
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_BHANDHLYP        = 436  !  Becke half-and-half with B88 exchange
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_MB3LYP_RC04      = 437  !  B3LYP with RC04 LDA
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_MPWLYP1M         = 453  !  MPW with 1 par. for metals/LYP
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_REVB3LYP         = 454  !  Revised B3LYP
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_CAMY_BLYP        = 455  ! BLYP with yukawa screening
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_PBE0_13          = 456  !  PBE0-1/3
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B3LYPS           = 459  !  B3LYP* functional
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_QTP17            = 460  !  global hybrid for vertical ionization potentials
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B3LYP_MCM1       = 461  !  B3LYP reoptimized in 6-31+G(2df,p) for enthalpies of formation
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B3LYP_MCM2       = 462  !  B3LYP reoptimized in 6-31+G(2df,p) for enthalpies of formation
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_WB97             = 463  ! Chai and Head-Gordon
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_WB97X            = 464  !  Chai and Head-Gordon
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LRC_WPBEH        = 465  !  Long-range corrected functional by Rorhdanz et al
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_WB97X_V          = 466  !  Mardirossian and Head-Gordon
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LCY_PBE          = 467  ! PBE with yukawa screening
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LCY_BLYP         = 468  ! BLYP with yukawa screening
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LC_VV10          = 469  !  Vydrov and Van Voorhis
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_CAMY_B3LYP       = 470  ! B3LYP with Yukawa screening
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_WB97X_D          = 471  !  Chai and Head-Gordon
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HPBEINT          = 472  !  hPBEint
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LRC_WPBE         = 473  !  Long-range corrected functional by Rorhdanz et al
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B3LYP5           = 475  !  B3LYP with VWN functional 5 instead of RPA
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_EDF2             = 476  ! Empirical functional from Lin, George and Gill
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_CAP0             = 477  !  Correct Asymptotic Potential hybrid
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LC_WPBE          = 478  !  Long-range corrected functional by Vydrov and Scuseria
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HSE12            = 479  !  HSE12 by Moussa, Schultz and Chelikowsky
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HSE12S           = 480  !  Short-range HSE12 by Moussa, Schultz, and Chelikowsky
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HSE_SOL          = 481  !  HSEsol functional by Schimka, Harl, and Kresse
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_CAM_QTP_01       = 482  !  CAM-QTP-01
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_MPW1LYP          = 483  !  Becke 1-parameter mixture of mPW91 and LYP
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_MPW1PBE          = 484  !  Becke 1-parameter mixture of mPW91 and PBE
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_KMLYP            = 485  !  Kang-Musgrave hybrid
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LC_WPBE_WHS      = 486  !  Long-range corrected functional by Weintraub, Henderson and Scuseria
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LC_WPBEH_WHS     = 487  !  Long-range corrected functional by Weintraub, Henderson and Scuseria
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LC_WPBE08_WHS    = 488  !  Long-range corrected functional by Weintraub, Henderson and Scuseria
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LC_WPBESOL_WHS   = 489  !  Long-range corrected functional by Weintraub, Henderson and Scuseria
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_CAM_QTP_00       = 490  ! CAM-QTP-00
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_CAM_QTP_02       = 491  !  CAM-QTP-02
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_LC_QTP           = 492  !  LC-QTP
  integer(c_int), parameter, public :: XC_GGA_X_S12G                  = 495  ! Swart 2012 GGA exchange
  integer(c_int), parameter, public :: XC_HYB_GGA_X_S12H              = 496  !  Swart 2012 GGA hybrid exchange
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_BLYP35           = 499  !  Becke 1-parameter mixture for mixed-valence systems
  integer(c_int), parameter, public :: XC_GGA_K_VW                    = 500  !  von Weiszaecker functional
  integer(c_int), parameter, public :: XC_GGA_K_GE2                   = 501  !  Second-order gradient expansion (l = 1/9)
  integer(c_int), parameter, public :: XC_GGA_K_GOLDEN                = 502  !  TF-lambda-vW form by Golden (l = 13/45)
  integer(c_int), parameter, public :: XC_GGA_K_YT65                  = 503  !  TF-lambda-vW form by Yonei and Tomishima (l = 1/5)
  integer(c_int), parameter, public :: XC_GGA_K_BALTIN                = 504  !  TF-lambda-vW form by Baltin (l = 5/9)
  integer(c_int), parameter, public :: XC_GGA_K_LIEB                  = 505  !  TF-lambda-vW form by Lieb (l = 0.185909191)
  integer(c_int), parameter, public :: XC_GGA_K_ABSP1                 = 506  !  gamma-TFvW form by Acharya et al [g = 1 - 1.412/N^(1/3)]
  integer(c_int), parameter, public :: XC_GGA_K_ABSP2                 = 507  !  gamma-TFvW form by Acharya et al [g = 1 - 1.332/N^(1/3)]
  integer(c_int), parameter, public :: XC_GGA_K_GR                    = 508  !  gamma-TFvW form by Gazquez and Robles
  integer(c_int), parameter, public :: XC_GGA_K_LUDENA                = 509  !  gamma-TFvW form by Ludena
  integer(c_int), parameter, public :: XC_GGA_K_GP85                  = 510  !  gamma-TFvW form by Ghosh and Parr
  integer(c_int), parameter, public :: XC_GGA_K_PEARSON               = 511  ! Pearson
  integer(c_int), parameter, public :: XC_GGA_K_OL1                   = 512  ! Ou-Yang and Levy v.1
  integer(c_int), parameter, public :: XC_GGA_K_OL2                   = 513  ! Ou-Yang and Levy v.2
  integer(c_int), parameter, public :: XC_GGA_K_FR_B88                = 514  !  Fuentealba & Reyes (B88 version)
  integer(c_int), parameter, public :: XC_GGA_K_FR_PW86               = 515  ! Fuentealba & Reyes (PW86 version)
  integer(c_int), parameter, public :: XC_GGA_K_DK                    = 516  ! DePristo and Kress
  integer(c_int), parameter, public :: XC_GGA_K_PERDEW                = 517  !  Perdew
  integer(c_int), parameter, public :: XC_GGA_K_VSK                   = 518  !  Vitos, Skriver, and Kollar
  integer(c_int), parameter, public :: XC_GGA_K_VJKS                  = 519  !  Vitos, Johansson, Kollar, and Skriver
  integer(c_int), parameter, public :: XC_GGA_K_ERNZERHOF             = 520  !  Ernzerhof
  integer(c_int), parameter, public :: XC_GGA_K_LC94                  = 521  ! Lembarki & Chermette
  integer(c_int), parameter, public :: XC_GGA_K_LLP                   = 522  ! Lee, Lee & Parr
  integer(c_int), parameter, public :: XC_GGA_K_THAKKAR               = 523  ! Thakkar 1992
  integer(c_int), parameter, public :: XC_GGA_X_WPBEH                 = 524  ! short-range version of the PBE
  integer(c_int), parameter, public :: XC_GGA_X_HJS_PBE               = 525  ! HJS screened exchange PBE version
  integer(c_int), parameter, public :: XC_GGA_X_HJS_PBE_SOL           = 526  !  HJS screened exchange PBE_SOL version
  integer(c_int), parameter, public :: XC_GGA_X_HJS_B88               = 527  !  HJS screened exchange B88 version
  integer(c_int), parameter, public :: XC_GGA_X_HJS_B97X              = 528  !  HJS screened exchange B97x version
  integer(c_int), parameter, public :: XC_GGA_X_ITYH                  = 529  ! short-range recipe B88 functionals - erf
  integer(c_int), parameter, public :: XC_GGA_X_SFAT                  = 530  ! short-range recipe for PBE functional
  integer(c_int), parameter, public :: XC_GGA_X_SG4                   = 533  ! Semiclassical GGA at fourth order
  integer(c_int), parameter, public :: XC_GGA_C_SG4                   = 534  ! Semiclassical GGA at fourth order
  integer(c_int), parameter, public :: XC_GGA_X_GG99                  = 535  ! Gilbert and Gill 1999
  integer(c_int), parameter, public :: XC_GGA_X_PBEPOW                = 539  ! PBE power
  integer(c_int), parameter, public :: XC_GGA_X_KGG99                 = 544  !  Gilbert and Gill 1999 (mixed)
  integer(c_int), parameter, public :: XC_GGA_XC_HLE16                = 545  !  high local exchange 2016
  integer(c_int), parameter, public :: XC_GGA_C_SCAN_E0               = 553  ! GGA component of SCAN
  integer(c_int), parameter, public :: XC_GGA_C_GAPC                  = 555  ! GapC
  integer(c_int), parameter, public :: XC_GGA_C_GAPLOC                = 556  ! Gaploc
  integer(c_int), parameter, public :: XC_GGA_C_ZVPBEINT              = 557  ! another spin-dependent correction to PBEint
  integer(c_int), parameter, public :: XC_GGA_C_ZVPBESOL              = 558  !  another spin-dependent correction to PBEsol
  integer(c_int), parameter, public :: XC_GGA_C_TM_LYP                = 559  !  Takkar and McCarthy reparametrization
  integer(c_int), parameter, public :: XC_GGA_C_TM_PBE                = 560  !  Thakkar and McCarthy reparametrization
  integer(c_int), parameter, public :: XC_GGA_C_W94                   = 561  ! Wilson 94 (Eq. 25)
  integer(c_int), parameter, public :: XC_GGA_C_CS1                   = 565  ! A dynamical correlation functional
  integer(c_int), parameter, public :: XC_GGA_X_B88M                  = 570  !  Becke 88 reoptimized to be used with mgga_c_tau1
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B5050LYP         = 572  !  Like B3LYP but more exact exchange
  integer(c_int), parameter, public :: XC_GGA_XC_KT3                  = 587  !  Keal and Tozer version 3
  integer(c_int), parameter, public :: XC_GGA_XC_LB07                 = 589  !  Livshits and Baer, empirical functional
  integer(c_int), parameter, public :: XC_GGA_K_GDS08                 = 591  ! Combined analytical theory with Monte Carlo sampling
  integer(c_int), parameter, public :: XC_GGA_K_GHDS10                = 592  !  As GDS08 but for an electron gas with spin
  integer(c_int), parameter, public :: XC_GGA_K_GHDS10R               = 593  !  Reparametrized GHDS10
  integer(c_int), parameter, public :: XC_GGA_K_TKVLN                 = 594  !  Trickey, Karasiev, and Vela
  integer(c_int), parameter, public :: XC_GGA_K_PBE3                  = 595  ! Three parameter PBE-like expansion
  integer(c_int), parameter, public :: XC_GGA_K_PBE4                  = 596  !  Four  parameter PBE-like expansion
  integer(c_int), parameter, public :: XC_GGA_K_EXP4                  = 597  ! Intermediate form between PBE3 and PBE4
  integer(c_int), parameter, public :: XC_GGA_X_SFAT_PBE              = 601  ! short-range recipe for PBE functional
  integer(c_int), parameter, public :: XC_GGA_X_FD_LB94               = 604  ! Functional derivative recovered from the stray LB94 potential
  integer(c_int), parameter, public :: XC_GGA_X_FD_REVLB94            = 605  !  Revised FD_LB94
  integer(c_int), parameter, public :: XC_GGA_C_ZVPBELOC              = 606  ! PBEloc variation with enhanced compatibility with exact exchange
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_APBE0            = 607  !  Hybrid based on APBE
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_HAPBE            = 608  !  Hybrid based in APBE and zvPBEloc
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_RCAM_B3LYP       = 610  !  Similar to CAM-B3LYP, but trying to reduce the many-electron self-interaction
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_WC04             = 611  !  hybrid fitted to carbon NMR shifts
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_WP04             = 612  !  hybrid fitted to proton NMR shifts
  integer(c_int), parameter, public :: XC_GGA_K_LKT                   = 613  ! Luo-Karasiev-Trickey kinetic GGA
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_CAM_PBEH         = 681  !  CAM version of PBEH
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_CAMY_PBEH        = 682  !  PBEH with Yukawa screening
  integer(c_int), parameter, public :: XC_HYB_GGA_X_LCGAU             = 708  ! Long-range Gaussian
  integer(c_int), parameter, public :: XC_HYB_GGA_X_LCGAU_CORE        = 709  !  Long-range Gaussian fitted to core excitations
  integer(c_int), parameter, public :: XC_HYB_GGA_X_LC2GAU            = 710  !  Long-range Gaussian 2
  integer(c_int), parameter, public :: XC_GGA_C_MGGAC                 = 712  !  beta fitted to LC20 to be used with MGGAC
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_B2PLYP           = 713  ! Double hybrid of Grimme
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_SRC1_BLYP        = 714  ! Hybrid with two range separations (form 1)
  integer(c_int), parameter, public :: XC_HYB_GGA_XC_SRC2_BLYP        = 715  !  Hybrid with two range separations (form 2)
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_DLDF             =  36  ! Dispersionless Density Functional
  integer(c_int), parameter, public :: XC_MGGA_C_DLDF                 =  37  !  Dispersionless Density Functional
  integer(c_int), parameter, public :: XC_MGGA_XC_ZLP                 =  42  ! Zhao, Levy & Parr, Eq. (21)
  integer(c_int), parameter, public :: XC_MGGA_XC_OTPSS_D             =  64  ! oTPSS_D functional of Goerigk and Grimme
  integer(c_int), parameter, public :: XC_MGGA_C_CS                   =  72  ! Colle and Salvetti
  integer(c_int), parameter, public :: XC_MGGA_C_MN12_SX              =  73  !  MN12-SX correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_MN12_L               =  74  !  MN12-L correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_M11_L                =  75  !  M11-L correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_M11                  =  76  !  M11 correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_M08_SO               =  77  !  M08-SO correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_M08_HX               =  78  ! M08-HX correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_REVM11               = 172  !  Revised M11 correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_X_LTA                  = 201  ! Local tau approximation of Ernzerhof & Scuseria
  integer(c_int), parameter, public :: XC_MGGA_X_TPSS                 = 202  ! Tao, Perdew, Staroverov & Scuseria exchange
  integer(c_int), parameter, public :: XC_MGGA_X_M06_L                = 203  ! M06-L exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_X_GVT4                 = 204  ! GVT4 from Van Voorhis and Scuseria
  integer(c_int), parameter, public :: XC_MGGA_X_TAU_HCTH             = 205  ! tau-HCTH from Boese and Handy
  integer(c_int), parameter, public :: XC_MGGA_X_BR89                 = 206  ! Becke-Roussel 89, gamma = 0.8
  integer(c_int), parameter, public :: XC_MGGA_X_BJ06                 = 207  ! Becke & Johnson correction to Becke-Roussel 89
  integer(c_int), parameter, public :: XC_MGGA_X_TB09                 = 208  !  Tran & Blaha correction to Becke & Johnson
  integer(c_int), parameter, public :: XC_MGGA_X_RPP09                = 209  !  Rasanen, Pittalis, and Proetto correction to Becke & Johnson
  integer(c_int), parameter, public :: XC_MGGA_X_2D_PRHG07            = 210  ! Pittalis, Rasanen, Helbig, Gross Exchange Functional
  integer(c_int), parameter, public :: XC_MGGA_X_2D_PRHG07_PRP10      = 211  ! PRGH07 with PRP10 correction
  integer(c_int), parameter, public :: XC_MGGA_X_REVTPSS              = 212  !  revised Tao, Perdew, Staroverov & Scuseria exchange
  integer(c_int), parameter, public :: XC_MGGA_X_PKZB                 = 213  ! Perdew, Kurth, Zupan, and Blaha
  integer(c_int), parameter, public :: XC_MGGA_X_BR89_1               = 214  !  Becke-Roussel 89, gamma = 1.0
  integer(c_int), parameter, public :: XC_MGGA_X_MS0                  = 221  ! MS exchange of Sun, Xiao, and Ruzsinszky
  integer(c_int), parameter, public :: XC_MGGA_X_MS1                  = 222  !  MS1 exchange of Sun, et al
  integer(c_int), parameter, public :: XC_MGGA_X_MS2                  = 223  !  MS2 exchange of Sun, et al
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_MS2H             = 224  !  MS2 hybrid exchange of Sun, et al
  integer(c_int), parameter, public :: XC_MGGA_X_TH                   = 225  ! Tsuneda and Hirao
  integer(c_int), parameter, public :: XC_MGGA_X_M11_L                = 226  ! M11-L exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_X_MN12_L               = 227  ! MN12-L exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_X_MS2_REV              = 228  !  MS2 exchange of Sun, et al with a revised value for c
  integer(c_int), parameter, public :: XC_MGGA_XC_CC06                = 229  ! Cancio and Chou 2006
  integer(c_int), parameter, public :: XC_MGGA_X_MK00                 = 230  !  Exchange for accurate virtual orbital energies
  integer(c_int), parameter, public :: XC_MGGA_C_TPSS                 = 231  ! Tao, Perdew, Staroverov & Scuseria correlation
  integer(c_int), parameter, public :: XC_MGGA_C_VSXC                 = 232  ! VSxc from Van Voorhis and Scuseria (correlation part)
  integer(c_int), parameter, public :: XC_MGGA_C_M06_L                = 233  ! M06-L correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_M06_HF               = 234  !  M06-HF correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_M06                  = 235  !  M06 correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_M06_2X               = 236  !  M06-2X correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_M05                  = 237  ! M05 correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_M05_2X               = 238  !  M05-2X correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_PKZB                 = 239  ! Perdew, Kurth, Zupan, and Blaha
  integer(c_int), parameter, public :: XC_MGGA_C_BC95                 = 240  ! Becke correlation 95
  integer(c_int), parameter, public :: XC_MGGA_C_REVTPSS              = 241  ! revised TPSS correlation
  integer(c_int), parameter, public :: XC_MGGA_XC_TPSSLYP1W           = 242  !  Functionals fitted for water
  integer(c_int), parameter, public :: XC_MGGA_X_MK00B                = 243  !  Exchange for accurate virtual orbital energies (v. B)
  integer(c_int), parameter, public :: XC_MGGA_X_BLOC                 = 244  !  functional with balanced localization
  integer(c_int), parameter, public :: XC_MGGA_X_MODTPSS              = 245  !  Modified Tao, Perdew, Staroverov & Scuseria exchange
  integer(c_int), parameter, public :: XC_MGGA_C_TPSSLOC              = 247  ! Semilocal dynamical correlation
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_MN12_SX          = 248  !  MN12-SX hybrid exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_X_MBEEF                = 249  ! mBEEF exchange
  integer(c_int), parameter, public :: XC_MGGA_X_MBEEFVDW             = 250  ! mBEEF-vdW exchange
  integer(c_int), parameter, public :: XC_MGGA_C_TM                   = 251  !  Tao and Mo 2016 correlation
  integer(c_int), parameter, public :: XC_MGGA_XC_B97M_V              = 254  ! Mardirossian and Head-Gordon
  integer(c_int), parameter, public :: XC_MGGA_X_MVS                  = 257  ! MVS exchange of Sun, Perdew, and Ruzsinszky
  integer(c_int), parameter, public :: XC_MGGA_X_MN15_L               = 260  !  MN15-L exhange functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_MN15_L               = 261  !  MN15-L correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_X_SCAN                 = 263  ! SCAN exchange of Sun, Ruzsinszky, and Perdew
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_SCAN0            = 264  !  SCAN hybrid exchange
  integer(c_int), parameter, public :: XC_MGGA_C_SCAN                 = 267  ! SCAN correlation
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_MN15             = 268  !  MN15 hybrid exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_MN15                 = 269  !  MN15 correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_BMK              = 279  !  Boese-Martin for kinetics
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_TAU_HCTH         = 282  !  Hybrid version of tau-HCTH
  integer(c_int), parameter, public :: XC_MGGA_X_B00                  = 284  !  Becke 2000
  integer(c_int), parameter, public :: XC_MGGA_XC_HLE17               = 288  ! high local exchange 2017
  integer(c_int), parameter, public :: XC_MGGA_C_SCAN_RVV10           = 292  !  SCAN correlation + rVV10 correlation
  integer(c_int), parameter, public :: XC_MGGA_X_REVM06_L             = 293  !  revised M06-L exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_REVM06_L             = 294  !  Revised M06-L correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_M08_HX           = 295  ! M08-HX exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_M08_SO           = 296  !  M08-SO exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_M11              = 297  ! M11 hybrid exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_X_RTPSS                = 299  ! Revised TPSS exchange by Garza, Bell and Head-Gordon
  integer(c_int), parameter, public :: XC_MGGA_X_MS2B                 = 300  ! MS2beta exchange by Furness and Sun
  integer(c_int), parameter, public :: XC_MGGA_X_MS2BS                = 301  !  MS2beta* exchange by Furness and Sun
  integer(c_int), parameter, public :: XC_MGGA_X_MVSB                 = 302  ! MVSBeta exchange of Furness and Sun
  integer(c_int), parameter, public :: XC_MGGA_X_MVSBS                = 303  !  MVSBeta* exchange of Furness and Sun
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_REVM11           = 304  !  revM11 hybrid exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_REVM06           = 305  !  revised M06 hybrid exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_MGGA_C_REVM06               = 306  !  Revised M06 correlation functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_M05              = 438  ! M05 hybrid exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_M05_2X           = 439  !  M05-2X hybrid exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_B88B95          = 440  ! Mixture of B88 with BC95 (B1B95)
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_B86B95          = 441  !  Mixture of B86 with BC95
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_PW86B95         = 442  !  Mixture of PW86 with BC95
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_BB1K            = 443  !  Mixture of B88 with BC95 from Zhao and Truhlar
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_M06_HF           = 444  !  M06-HF hybrid exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_MPW1B95         = 445  !  Mixture of mPW91 with BC95 from Zhao and Truhlar
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_MPWB1K          = 446  !  Mixture of mPW91 with BC95 for kinetics
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_X1B95           = 447  !  Mixture of X with BC95
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_XB1K            = 448  !  Mixture of X with BC95 for kinetics
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_M06              = 449  !  M06 hybrid exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_M06_2X           = 450  !  M06-2X hybrid exchange functional from Minnesota
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_PW6B95          = 451  !  Mixture of PW91 with BC95 from Zhao and Truhlar
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_PWB6K           = 452  !  Mixture of PW91 with BC95 from Zhao and Truhlar for kinetics
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_TPSSH           = 457  ! TPSS hybrid
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_REVTPSSH        = 458  !  revTPSS hybrid
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_MVSH             = 474  ! MVSh hybrid
  integer(c_int), parameter, public :: XC_MGGA_X_RSCAN                = 493  ! Regularized SCAN exchange
  integer(c_int), parameter, public :: XC_MGGA_C_RSCAN                = 494  ! Regularized SCAN correlation
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_WB97M_V         = 531  ! Mardirossian and Head-Gordon
  integer(c_int), parameter, public :: XC_MGGA_X_TM                   = 540  ! Tao and Mo 2016 exchange
  integer(c_int), parameter, public :: XC_MGGA_X_VT84                 = 541  ! meta-GGA version of VT{8,4} GGA
  integer(c_int), parameter, public :: XC_MGGA_X_SA_TPSS              = 542  ! TPSS with correct surface asymptotics
  integer(c_int), parameter, public :: XC_MGGA_K_PC07                 = 543  ! Perdew and Constantin 2007
  integer(c_int), parameter, public :: XC_MGGA_C_KCIS                 = 562  ! Krieger, Chen, Iafrate, and Savin
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_B0KCIS          = 563  !  Hybrid based on KCIS
  integer(c_int), parameter, public :: XC_MGGA_XC_LP90                = 564  ! Lee & Parr, Eq. (56)
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_MPW1KCIS        = 566  ! Modified Perdew-Wang + KCIS hybrid
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_MPWKCIS1K       = 567  !  Modified Perdew-Wang + KCIS hybrid with more exact exchange
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_PBE1KCIS        = 568  !  Perdew-Burke-Ernzerhof + KCIS hybrid
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_TPSS1KCIS       = 569  !  TPSS hybrid with KCIS correlation
  integer(c_int), parameter, public :: XC_MGGA_C_B88                  = 571  ! Meta-GGA correlation by Becke
  integer(c_int), parameter, public :: XC_MGGA_X_GX                   = 575  ! GX functional of Loos
  integer(c_int), parameter, public :: XC_MGGA_X_PBE_GX               = 576  ! PBE-GX functional of Loos
  integer(c_int), parameter, public :: XC_MGGA_X_REVSCAN              = 581  !  revised SCAN
  integer(c_int), parameter, public :: XC_MGGA_C_REVSCAN              = 582  ! revised SCAN correlation
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_REVSCAN0         = 583  !  revised SCAN hybrid exchange
  integer(c_int), parameter, public :: XC_MGGA_C_SCAN_VV10            = 584  !  SCAN correlation +  VV10 correlation
  integer(c_int), parameter, public :: XC_MGGA_C_REVSCAN_VV10         = 585  !  revised SCAN correlation
  integer(c_int), parameter, public :: XC_MGGA_X_BR89_EXPLICIT        = 586  ! Becke-Roussel 89 with an explicit inversion of x(y), gamma = 0.8
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_B98             = 598  ! Becke 98
  integer(c_int), parameter, public :: XC_MGGA_X_BR89_EXPLICIT_1      = 602  !  Becke-Roussel 89 with an explicit inversion of x(y), gamma = 1.0
  integer(c_int), parameter, public :: XC_MGGA_X_REGTPSS              = 603  ! Regularized TPSS
  integer(c_int), parameter, public :: XC_MGGA_X_2D_JS17              = 609  ! JS17 meta-GGA for 2D
  integer(c_int), parameter, public :: XC_MGGA_X_TLDA                 = 685  !  LDA-type exchange with tau-dependent potential
  integer(c_int), parameter, public :: XC_MGGA_X_EDMGGA               = 686  ! Tao 2001
  integer(c_int), parameter, public :: XC_MGGA_X_GDME_NV              = 687  ! Generalized density-matrix with a=1/2
  integer(c_int), parameter, public :: XC_MGGA_X_RLDA                 = 688  ! Reparametrized local-density approximation
  integer(c_int), parameter, public :: XC_MGGA_X_GDME_0               = 689  !  Generalized density-matrix with a=0
  integer(c_int), parameter, public :: XC_MGGA_X_GDME_KOS             = 690  !  Generalized density-matrix with a=0.00638
  integer(c_int), parameter, public :: XC_MGGA_X_GDME_VT              = 691  !   Varied-terms (VT) mGGA of Koehl, Odom, and Scuseria
  integer(c_int), parameter, public :: XC_MGGA_X_REVTM                = 693  ! revised Tao and Mo 2016 exchange
  integer(c_int), parameter, public :: XC_MGGA_C_REVTM                = 694  !  revised Tao and Mo 2016 correlation
  integer(c_int), parameter, public :: XC_HYB_MGGA_XC_EDMGGAH         = 695  !  Tao 2001 hybrid
  integer(c_int), parameter, public :: XC_MGGA_X_MBRXC_BG             = 696  ! Modified Becke-Roussel for band gaps - cuspless hole
  integer(c_int), parameter, public :: XC_MGGA_X_MBRXH_BG             = 697  ! Modified Becke-Roussel for band gaps - hydrogen hole
  integer(c_int), parameter, public :: XC_MGGA_X_SCANL                = 700  ! Deorbitalized SCAN exchange
  integer(c_int), parameter, public :: XC_MGGA_X_REVSCANL             = 701  !  Deorbitalized revSCAN exchange
  integer(c_int), parameter, public :: XC_MGGA_C_SCANL                = 702  ! SCAN correlation
  integer(c_int), parameter, public :: XC_MGGA_C_SCANL_RVV10          = 703  !  SCAN correlation + rVV10 correlation
  integer(c_int), parameter, public :: XC_MGGA_C_SCANL_VV10           = 704  !  SCAN correlation +  VV10 correlation
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_JS18             = 705  ! a screened version of TM
  integer(c_int), parameter, public :: XC_HYB_MGGA_X_PJS18            = 706  ! a screened version of TM
  integer(c_int), parameter, public :: XC_MGGA_X_TASK                 = 707  ! TASK exchange of Aschebrock and Kuemmel
  integer(c_int), parameter, public :: XC_MGGA_X_MGGAC                = 711  ! MGGAC of Patras et al
  integer(c_int), parameter, public :: XC_MGGA_X_MBR                  = 716  ! modified Becke-Roussel

  ! These are old names kept for compatibility
  integer(c_int), parameter, public :: &
       XC_LDA_X_1D             =  21,     &
       XC_GGA_X_BGCP           =  38,     &
       XC_GGA_C_BGCP           =  39,     &
       XC_GGA_C_BCGP           =  39,     &
       XC_GGA_C_VPBE           =  83,     &
       XC_GGA_XC_LB            = 160,     &
       XC_MGGA_C_CC06          = 229,     &
       XC_GGA_K_ABSR1          = 506,     &
       XC_GGA_K_ABSR2          = 507,     &
       XC_LDA_C_LP_A           = 547,     &
       XC_LDA_C_LP_B           = 548,     &
       XC_MGGA_C_LP90          = 564

  !----------------------------------------------------------------
  interface
     subroutine xc_version(major, minor, micro) bind(c)
       import
       integer(c_int), intent(out) :: major, minor, micro
     end subroutine xc_version

     type(c_ptr) function xc_version_string() bind(c)
       import
     end function xc_version_string
  end interface


  !----------------------------------------------------------------
  type :: xc_f03_func_info_t
     private
     type(c_ptr) :: ptr = C_NULL_PTR
  end type xc_f03_func_info_t

  interface
     integer(c_int) function xc_func_info_get_number(info) bind(c)
       import
       type(c_ptr), value :: info
     end function xc_func_info_get_number

     integer(c_int) function xc_func_info_get_kind(info) bind(c)
       import
       type(c_ptr), value :: info
     end function xc_func_info_get_kind

     type(c_ptr) function xc_func_info_get_name(info) bind(c)
       import
       type(c_ptr), value :: info
     end function xc_func_info_get_name

     integer(c_int) function xc_func_info_get_family(info) bind(c)
       import
       type(c_ptr), value :: info
     end function xc_func_info_get_family

     integer(c_int) function xc_func_info_get_flags(info) bind(c)
       import
       type(c_ptr), value :: info
     end function xc_func_info_get_flags

     type(c_ptr) function xc_func_info_get_references(info, number) bind(c)
       import
       type(c_ptr),    value :: info
       integer(c_int), value :: number
     end function xc_func_info_get_references

     integer(c_int) function xc_func_info_get_n_ext_params(info) bind(c)
       import
       type(c_ptr), value :: info
     end function xc_func_info_get_n_ext_params

     type(c_ptr) function xc_func_info_get_ext_params_name(info, number) bind(c)
       import
       type(c_ptr),    value :: info
       integer(c_int), value :: number
     end function xc_func_info_get_ext_params_name

     type(c_ptr) function xc_func_info_get_ext_params_description(info, number) bind(c)
       import
       type(c_ptr),    value :: info
       integer(c_int), value :: number
     end function xc_func_info_get_ext_params_description

     real(c_double) function xc_func_info_get_ext_params_default_value(info, number) bind(c)
       import
       type(c_ptr),    value :: info
       integer(c_int), value :: number
     end function xc_func_info_get_ext_params_default_value

  end interface

  !----------------------------------------------------------------
  type :: xc_f03_func_reference_t
     private
     type(c_ptr) :: ptr = C_NULL_PTR
  end type xc_f03_func_reference_t

  interface
     type(c_ptr) function xc_func_reference_get_ref(reference) bind(c)
       import
       type(c_ptr), value :: reference
     end function xc_func_reference_get_ref

     type(c_ptr) function xc_func_reference_get_doi(reference) bind(c)
       import
       type(c_ptr), value :: reference
     end function xc_func_reference_get_doi

     type(c_ptr) function xc_func_reference_get_bibtex(reference) bind(c)
       import
       type(c_ptr), value :: reference
     end function xc_func_reference_get_bibtex
  end interface

  !----------------------------------------------------------------
  type :: xc_f03_func_t
     private
     type(c_ptr) :: ptr = C_NULL_PTR
  end type xc_f03_func_t

  interface
     type(c_ptr) function xc_func_alloc() bind(c)
       import
     end function xc_func_alloc

     integer(c_int) function xc_func_init(p, functional, nspin) bind(c)
       import
       type(c_ptr),    value :: p
       integer(c_int), value :: functional, nspin
     end function xc_func_init

     subroutine xc_func_end(p) bind(c)
       import
       type(c_ptr), value :: p
     end subroutine xc_func_end

     subroutine xc_func_free(p) bind(c)
       import
       type(c_ptr), value :: p
     end subroutine xc_func_free

     type(c_ptr) function xc_func_get_info(p) bind(c)
       import
       type(c_ptr), value :: p
     end function xc_func_get_info

     type(c_ptr) function xc_functional_get_name(number) bind(c)
       import
       integer(c_int), value :: number
     end function xc_functional_get_name

     integer(c_int) function xc_functional_get_number(func_string) bind(c)
       import
       character(kind=c_char), intent(in) :: func_string(*)
     end function xc_functional_get_number

     integer(c_int) function xc_family_from_id(id, family, number) bind(c)
       import
       integer(c_int), value :: id
       type(c_ptr),    value :: family, number
     end function xc_family_from_id

     integer(c_int) function xc_f03_number_of_functionals() bind(c, name="xc_number_of_functionals")
       import
     end function xc_f03_number_of_functionals

     integer(c_int) function xc_f03_maximum_name_length() bind(c, name="xc_maximum_name_length")
       import
     end function xc_f03_maximum_name_length

     subroutine xc_f03_available_functional_numbers(list) bind(c, name="xc_available_functional_numbers")
       import
       integer(c_int), intent(out) :: list(*)
     end subroutine xc_f03_available_functional_numbers

     subroutine xc_available_functional_names(list) bind(c)
       import
       type(c_ptr) :: list(*)
     end subroutine xc_available_functional_names

     subroutine xc_func_set_dens_threshold(p, dens_threshold) bind(c)
       import
       type(c_ptr), value :: p
       real(c_double), value :: dens_threshold
     end subroutine xc_func_set_dens_threshold

     subroutine xc_func_set_ext_params(p, ext_params) bind(c)
       import
       type(c_ptr), value      :: p
       real(c_double), intent(in) :: ext_params(*)
     end subroutine xc_func_set_ext_params
  end interface

  ! LDAs
  !----------------------------------------------------------------
  interface
     subroutine xc_lda(p, np, rho, zk, vrho, v2rho2, v3rho3, v4rho4) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*)
       real(c_double),    intent(out) :: zk(*), vrho(*), v2rho2(*), v3rho3(*), v4rho4(*)
     end subroutine xc_lda

     subroutine xc_lda_exc(p, np, rho, zk) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*)
       real(c_double),    intent(out) :: zk(*)
     end subroutine xc_lda_exc

     subroutine xc_lda_exc_vxc(p, np, rho, zk, vrho) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*)
       real(c_double),    intent(out) :: zk(*), vrho(*)
     end subroutine xc_lda_exc_vxc

     subroutine xc_lda_exc_vxc_fxc(p, np, rho, zk, vrho, v2rho2) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*)
       real(c_double),    intent(out) :: zk(*), vrho(*), v2rho2(*)
     end subroutine xc_lda_exc_vxc_fxc

     subroutine xc_lda_exc_vxc_fxc_kxc(p, np, rho, zk, vrho, v2rho2, v3rho3) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*)
       real(c_double),    intent(out) :: zk(*), vrho(*), v2rho2(*), v3rho3(*)
     end subroutine xc_lda_exc_vxc_fxc_kxc

     subroutine xc_lda_vxc(p, np, rho, vrho) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*)
       real(c_double),    intent(out) :: vrho(*)
     end subroutine xc_lda_vxc

     subroutine xc_lda_vxc_fxc(p, np, rho, vrho, v2rho2) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*)
       real(c_double),    intent(out) :: vrho(*), v2rho2(*)
     end subroutine xc_lda_vxc_fxc

     subroutine xc_lda_vxc_fxc_kxc(p, np, rho, vrho, v2rho2, v3rho3) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*)
       real(c_double),    intent(out) :: vrho(*), v2rho2(*), v3rho3(*)
     end subroutine xc_lda_vxc_fxc_kxc

     subroutine xc_lda_fxc(p, np, rho, v2rho2) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*)
       real(c_double),    intent(out) :: v2rho2(*)
     end subroutine xc_lda_fxc

     subroutine xc_lda_kxc(p, np, rho, v3rho3) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*)
       real(c_double),    intent(out) :: v3rho3(*)
     end subroutine xc_lda_kxc

     subroutine xc_lda_lxc(p, np, rho, v4rho4) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*)
       real(c_double),    intent(out) :: v4rho4(*)
     end subroutine xc_lda_lxc
  end interface


  ! GGAs
  !----------------------------------------------------------------
  interface
     subroutine xc_gga(p, np, rho, sigma, zk, vrho, vsigma,        &
          v2rho2, v2rhosigma, v2sigma2,                            &
          v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,              &
          v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4 &
          ) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*)
       real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*)
       real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
       real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
       real(c_double),    intent(out) :: v4rho4(*), v4rho3sigma(*), v4rho2sigma2(*), v4rhosigma3(*), v4sigma4(*)
     end subroutine xc_gga

     subroutine xc_gga_exc(p, np, rho, sigma, zk) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*)
       real(c_double),    intent(out) :: zk(*)
     end subroutine xc_gga_exc

     subroutine xc_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*)
       real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*)
     end subroutine xc_gga_exc_vxc

     subroutine xc_gga_exc_vxc_fxc(p, np, rho, sigma, zk, vrho, vsigma,        &
          v2rho2, v2rhosigma, v2sigma2) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*)
       real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*)
       real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
     end subroutine xc_gga_exc_vxc_fxc

     subroutine xc_gga_exc_vxc_fxc_kxc(p, np, rho, sigma, zk, vrho, vsigma,        &
          v2rho2, v2rhosigma, v2sigma2,                            &
          v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*)
       real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*)
       real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
       real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
     end subroutine xc_gga_exc_vxc_fxc_kxc

     subroutine xc_gga_vxc(p, np, rho, sigma, vrho, vsigma) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*)
       real(c_double),    intent(out) :: vrho(*), vsigma(*)
     end subroutine xc_gga_vxc

     subroutine xc_gga_vxc_fxc(p, np, rho, sigma, vrho, vsigma,    &
          v2rho2, v2rhosigma, v2sigma2) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*)
       real(c_double),    intent(out) :: vrho(*), vsigma(*)
       real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
     end subroutine xc_gga_vxc_fxc

     subroutine xc_gga_vxc_fxc_kxc(p, np, rho, sigma, vrho, vsigma, &
          v2rho2, v2rhosigma, v2sigma2,                             &
          v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*)
       real(c_double),    intent(out) :: vrho(*), vsigma(*)
       real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
       real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
     end subroutine xc_gga_vxc_fxc_kxc

     subroutine xc_gga_fxc(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*)
       real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
     end subroutine xc_gga_fxc

     subroutine xc_gga_kxc(p, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*)
       real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
     end subroutine xc_gga_kxc

     subroutine xc_gga_lxc(p, np, rho, sigma, v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*)
       real(c_double),    intent(out) :: v4rho4(*), v4rho3sigma(*), v4rho2sigma2(*), v4rhosigma3(*), v4sigma4(*)
     end subroutine xc_gga_lxc

  end interface


  interface
     real(c_double) function xc_gga_ak13_get_asymptotic(homo) bind(c)
       import
       real(c_double), value :: homo
     end function xc_gga_ak13_get_asymptotic
  end interface


  interface
     integer(c_int) function xc_hyb_type(p) bind(c)
       import
       type(c_ptr), value :: p
     end function xc_hyb_type

     real(c_double) function xc_hyb_exx_coef(p) bind(c)
       import
       type(c_ptr), value :: p
     end function xc_hyb_exx_coef

     subroutine xc_hyb_cam_coef(p, omega, alpha, beta) bind(c)
       import
       type(c_ptr), value       :: p
       real(c_double), intent(out) :: omega, alpha, beta
     end subroutine xc_hyb_cam_coef

     subroutine xc_nlc_coef(p, nlc_b, nlc_c) bind(c)
       import
       type(c_ptr), value       :: p
       real(c_double), intent(out) :: nlc_b, nlc_c
     end subroutine xc_nlc_coef
  end interface


  ! the meta-GGAs
  !----------------------------------------------------------------
  interface
     subroutine xc_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,     &
          v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
          v2lapl2, v2lapltau, v2tau2,                                                    &
          v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
          v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
          v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
          v3lapltau2, v3tau3,                                                            &
          v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho2sigma2, v4rho2sigmalapl,     &
          v4rho2sigmatau, v4rho2lapl2, v4rho2lapltau, v4rho2tau2, v4rhosigma3,           &
          v4rhosigma2lapl, v4rhosigma2tau, v4rhosigmalapl2, v4rhosigmalapltau,           &
          v4rhosigmatau2, v4rholapl3, v4rholapl2tau, v4rholapltau2, v4rhotau3, v4sigma4, &
          v4sigma3lapl, v4sigma3tau, v4sigma2lapl2, v4sigma2lapltau, v4sigma2tau2,       &
          v4sigmalapl3, v4sigmalapl2tau, v4sigmalapltau2, v4sigmatau3, v4lapl4,          &
          v4lapl3tau, v4lapl2tau2, v4lapltau3, v4tau4 &
          ) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
       real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
       real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
            v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
       real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
            v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),                  &
            v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),         &
            v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),       &
            v3lapltau2(*), v3tau3(*)
       real(c_double),    intent(out) :: &
            v4rho4(*), v4rho3sigma(*), v4rho3lapl(*), v4rho3tau(*), v4rho2sigma2(*), v4rho2sigmalapl(*),     &
            v4rho2sigmatau(*), v4rho2lapl2(*), v4rho2lapltau(*), v4rho2tau2(*), v4rhosigma3(*),              &
            v4rhosigma2lapl(*), v4rhosigma2tau(*), v4rhosigmalapl2(*), v4rhosigmalapltau(*),                 &
            v4rhosigmatau2(*), v4rholapl3(*), v4rholapl2tau(*), v4rholapltau2(*), v4rhotau3(*), v4sigma4(*), &
            v4sigma3lapl(*), v4sigma3tau(*), v4sigma2lapl2(*), v4sigma2lapltau(*), v4sigma2tau2(*),          &
            v4sigmalapl3(*), v4sigmalapl2tau(*), v4sigmalapltau2(*), v4sigmatau3(*), v4lapl4(*),             &
            v4lapl3tau(*), v4lapl2tau2(*), v4lapltau3(*), v4tau4(*)
     end subroutine xc_mgga

     subroutine xc_mgga_exc(p, np, rho, sigma, lapl, tau, zk) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
       real(c_double),    intent(out) :: zk(*)
     end subroutine xc_mgga_exc

     subroutine xc_mgga_exc_vxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
       real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
     end subroutine xc_mgga_exc_vxc

     subroutine xc_mgga_exc_vxc_fxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau, &
          v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
          v2lapl2, v2lapltau, v2tau2) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
       real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
       real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
            v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
     end subroutine xc_mgga_exc_vxc_fxc

     subroutine xc_mgga_exc_vxc_fxc_kxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,     &
          v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
          v2lapl2, v2lapltau, v2tau2,                                                    &
          v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
          v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
          v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
          v3lapltau2, v3tau3) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
       real(c_double),    intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
       real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
            v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
       real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
            v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),                  &
            v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),         &
            v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),       &
            v3lapltau2(*), v3tau3(*)
     end subroutine xc_mgga_exc_vxc_fxc_kxc

     subroutine xc_mgga_vxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
       real(c_double),    intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)
     end subroutine xc_mgga_vxc

     subroutine xc_mgga_vxc_fxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau, &
          v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
          v2lapl2, v2lapltau, v2tau2) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
       real(c_double),    intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)
       real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
            v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
     end subroutine xc_mgga_vxc_fxc

     subroutine xc_mgga_vxc_fxc_kxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau,     &
          v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
          v2lapl2, v2lapltau, v2tau2,                                                    &
          v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
          v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
          v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
          v3lapltau2, v3tau3) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
       real(c_double),    intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)
       real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
            v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
       real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
            v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),                  &
            v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),         &
            v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),       &
            v3lapltau2(*), v3tau3(*)
     end subroutine xc_mgga_vxc_fxc_kxc

     subroutine xc_mgga_fxc(p, np, rho, sigma, lapl, tau, &
          v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
          v2lapl2, v2lapltau, v2tau2) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
       real(c_double),    intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
            v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
     end subroutine xc_mgga_fxc

     subroutine xc_mgga_kxc(p, np, rho, sigma, lapl, tau, &
          v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
          v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
          v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
          v3lapltau2, v3tau3) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
       real(c_double),    intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
            v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),               &
            v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),      &
            v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),    &
            v3lapltau2(*), v3tau3(*)
     end subroutine xc_mgga_kxc

     subroutine xc_mgga_lxc(p, np, rho, sigma, lapl, tau, &
          v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho2sigma2, v4rho2sigmalapl,     &
          v4rho2sigmatau, v4rho2lapl2, v4rho2lapltau, v4rho2tau2, v4rhosigma3,           &
          v4rhosigma2lapl, v4rhosigma2tau, v4rhosigmalapl2, v4rhosigmalapltau,           &
          v4rhosigmatau2, v4rholapl3, v4rholapl2tau, v4rholapltau2, v4rhotau3, v4sigma4, &
          v4sigma3lapl, v4sigma3tau, v4sigma2lapl2, v4sigma2lapltau, v4sigma2tau2,       &
          v4sigmalapl3, v4sigmalapl2tau, v4sigmalapltau2, v4sigmatau3, v4lapl4,          &
          v4lapl3tau, v4lapl2tau2, v4lapltau3, v4tau4 &
          ) bind(c)
       import
       type(c_ptr),       value       :: p
       integer(c_size_t), value       :: np
       real(c_double),    intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
       real(c_double),    intent(out) :: &
            v4rho4(*), v4rho3sigma(*), v4rho3lapl(*), v4rho3tau(*), v4rho2sigma2(*), v4rho2sigmalapl(*),     &
            v4rho2sigmatau(*), v4rho2lapl2(*), v4rho2lapltau(*), v4rho2tau2(*), v4rhosigma3(*),              &
            v4rhosigma2lapl(*), v4rhosigma2tau(*), v4rhosigmalapl2(*), v4rhosigmalapltau(*),                 &
            v4rhosigmatau2(*), v4rholapl3(*), v4rholapl2tau(*), v4rholapltau2(*), v4rhotau3(*), v4sigma4(*), &
            v4sigma3lapl(*), v4sigma3tau(*), v4sigma2lapl2(*), v4sigma2lapltau(*), v4sigma2tau2(*),          &
            v4sigmalapl3(*), v4sigmalapl2tau(*), v4sigmalapltau2(*), v4sigmatau3(*), v4lapl4(*),             &
            v4lapl3tau(*), v4lapl2tau2(*), v4lapltau3(*), v4tau4(*)
     end subroutine xc_mgga_lxc
  end interface

contains

  !----------------------------------------------------------------
  subroutine xc_f03_version(major, minor, micro)
    integer(c_int), intent(out) :: major, minor, micro

    call xc_version(major, minor, micro)

  end subroutine xc_f03_version

  subroutine xc_f03_version_string(version)
    character(len=*), intent(out) :: version

    type(c_ptr) :: c_version

    c_version = xc_version_string()
    call c_to_f_string_ptr(c_version, version)

  end subroutine xc_f03_version_string


  !----------------------------------------------------------------
  integer(c_int) function xc_f03_func_info_get_number(info) result(number)
    type(xc_f03_func_info_t), intent(in) :: info

    number = xc_func_info_get_number(info%ptr)

  end function xc_f03_func_info_get_number

  integer(c_int) function xc_f03_func_info_get_kind(info) result(kind)
    type(xc_f03_func_info_t), intent(in) :: info

    kind = xc_func_info_get_kind(info%ptr)

  end function xc_f03_func_info_get_kind

  character(len=128) function xc_f03_func_info_get_name(info) result(name)
    type(xc_f03_func_info_t), intent(in) :: info

    call c_to_f_string_ptr(xc_func_info_get_name(info%ptr), name)

  end function xc_f03_func_info_get_name

  integer(c_int) function xc_f03_func_info_get_family(info) result(family)
    type(xc_f03_func_info_t), intent(in) :: info

    family = xc_func_info_get_family(info%ptr)

  end function xc_f03_func_info_get_family

  integer(c_int) function xc_f03_func_info_get_flags(info) result(flags)
    type(xc_f03_func_info_t), intent(in) :: info

    flags = xc_func_info_get_flags(info%ptr)

  end function xc_f03_func_info_get_flags

  type(xc_f03_func_reference_t) function xc_f03_func_info_get_references(info, number) result(reference)
    type(xc_f03_func_info_t), intent(in)    :: info
    integer(c_int),           intent(inout) :: number ! number of the reference. Must be 0 in the first call

    type(c_ptr) :: next_ref

    reference%ptr = xc_func_info_get_references(info%ptr, number)
    next_ref = xc_func_info_get_references(info%ptr, number + 1)
    if (c_associated(next_ref)) then
       number = number + 1
    else
       number = -1
    end if

  end function xc_f03_func_info_get_references

  integer(c_int) function xc_f03_func_info_get_n_ext_params(info) result(n_ext_params)
    type(xc_f03_func_info_t), intent(in) :: info

    n_ext_params = xc_func_info_get_n_ext_params(info%ptr)

  end function xc_f03_func_info_get_n_ext_params

  character(len=128) function xc_f03_func_info_get_ext_params_name(info, number) result(name)
    type(xc_f03_func_info_t), intent(in) :: info
    integer(c_int),           intent(in) :: number

    call c_to_f_string_ptr(xc_func_info_get_ext_params_name(info%ptr, number), name)

  end function xc_f03_func_info_get_ext_params_name

  character(len=128) function xc_f03_func_info_get_ext_params_description(info, number) result(description)
    type(xc_f03_func_info_t), intent(in) :: info
    integer(c_int),           intent(in) :: number

    call c_to_f_string_ptr(xc_func_info_get_ext_params_description(info%ptr, number), description)

  end function xc_f03_func_info_get_ext_params_description

  real(c_double) function xc_f03_func_info_get_ext_params_default_value(info, number) result(val)
    type(xc_f03_func_info_t), intent(in) :: info
    integer(c_int),           intent(in) :: number

    val = xc_func_info_get_ext_params_default_value(info%ptr, number)

  end function xc_f03_func_info_get_ext_params_default_value

  !----------------------------------------------------------------
  character(len=120) function xc_f03_func_reference_get_ref(reference) result(ref)
    type(xc_f03_func_reference_t), intent(in) :: reference

    call c_to_f_string_ptr(xc_func_reference_get_ref(reference%ptr), ref)

  end function xc_f03_func_reference_get_ref

  character(len=120) function xc_f03_func_reference_get_doi(reference) result(doi)
    type(xc_f03_func_reference_t), intent(in) :: reference

    call c_to_f_string_ptr(xc_func_reference_get_doi(reference%ptr), doi)

  end function xc_f03_func_reference_get_doi

  character(len=120) function xc_f03_func_reference_get_bibtex(reference) result(bibtex)
    type(xc_f03_func_reference_t), intent(in) :: reference

    call c_to_f_string_ptr(xc_func_reference_get_bibtex(reference%ptr), bibtex)

  end function xc_f03_func_reference_get_bibtex


  !----------------------------------------------------------------
  subroutine xc_f03_func_init(p, functional, nspin, err)
    type(xc_f03_func_t),      intent(inout) :: p
    integer(c_int),           intent(in)    :: functional
    integer(c_int),           intent(in)    :: nspin
    integer(c_int), optional, intent(out)   :: err

    integer(c_int) :: ierr

    p%ptr = xc_func_alloc()
    ierr = xc_func_init(p%ptr, functional, nspin)

    if(present(err)) err = ierr
  end subroutine xc_f03_func_init

  subroutine xc_f03_func_end(p)
    type(xc_f03_func_t), intent(inout) :: p

    call xc_func_end(p%ptr)
    call xc_func_free(p%ptr)

  end subroutine xc_f03_func_end

  type(xc_f03_func_info_t) function xc_f03_func_get_info(p) result(info)
    type(xc_f03_func_t), intent(in) :: p

    info%ptr = xc_func_get_info(p%ptr)

  end function xc_f03_func_get_info

  character(len=128) function xc_f03_functional_get_name(number) result(name)
    integer(c_int), intent(in) :: number

    call c_to_f_string_ptr(xc_functional_get_name(number), name)

  end function xc_f03_functional_get_name

  integer(c_int) function xc_f03_functional_get_number(func_string) result(number)
    character(len=*), intent(in) :: func_string

    number = xc_functional_get_number(f_to_c_string(func_string))

  end function xc_f03_functional_get_number

  integer(c_int) function xc_f03_family_from_id(id, family, number)
    integer(c_int), intent(in)                    :: id
    integer(c_int), intent(out), optional, target :: family, number

    type(c_ptr) c_family, c_number
    integer(c_int), pointer :: f_family, f_number

    if (present(family)) then
       f_family => family
       call c_f_pointer(c_family, f_family)
    else
       c_family = C_NULL_PTR
    end if
    if (present(number)) then
       f_number => number
       call c_f_pointer(c_number, f_number)
    else
       c_number = C_NULL_PTR
    end if

    xc_f03_family_from_id = xc_family_from_id(id, c_family, c_number)

  end function xc_f03_family_from_id

  subroutine xc_f03_available_functional_names(list)
    character(len=*), intent(out) :: list(*)

    integer(c_int) :: n, i, maxlen
    character(kind=c_char), allocatable, target :: names(:,:)
    type(c_ptr), allocatable :: c_list(:)

    n = xc_f03_number_of_functionals()
    maxlen = xc_f03_maximum_name_length()

    allocate(names(maxlen, n))
    allocate(c_list(n))
    do i = 1, n
       c_list(i) = c_loc(names(1,i))
    end do

    call xc_available_functional_names(c_list)

    do i = 1, n
       call c_to_f_string_ptr(c_list(i), list(i))
    end do

    deallocate(c_list)
    deallocate(names)

  end subroutine xc_f03_available_functional_names


  subroutine xc_f03_func_set_dens_threshold(p, dens_threshold)
    type(xc_f03_func_t), intent(in) :: p
    real(c_double),       intent(in) :: dens_threshold

    call xc_func_set_dens_threshold(p%ptr, dens_threshold)

  end subroutine xc_f03_func_set_dens_threshold

  subroutine xc_f03_func_set_ext_params(p, ext_params)
    type(xc_f03_func_t), intent(in) :: p
    real(c_double),       intent(in) :: ext_params(*)

    call xc_func_set_ext_params(p%ptr, ext_params)

  end subroutine xc_f03_func_set_ext_params

  ! LDAs
  !----------------------------------------------------------------
  subroutine xc_f03_lda(p, np, rho, zk, vrho, v2rho2, v3rho3, v4rho4)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), v2rho2(*), v3rho3(*), v4rho4(*)

    call xc_lda(p%ptr, np, rho, zk, vrho, v2rho2, v3rho3, v4rho4)

  end subroutine xc_f03_lda

  subroutine xc_f03_lda_exc(p, np, rho, zk)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: zk(*)

    call xc_lda_exc(p%ptr, np, rho, zk)

  end subroutine xc_f03_lda_exc

  subroutine xc_f03_lda_exc_vxc(p, np, rho, zk, vrho)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: zk(*), vrho(*)

    call xc_lda_exc_vxc(p%ptr, np, rho, zk, vrho)

  end subroutine xc_f03_lda_exc_vxc

  subroutine xc_f03_lda_exc_vxc_fxc(p, np, rho, zk, vrho, v2rho2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), v2rho2(*)

    call xc_lda_exc_vxc_fxc(p%ptr, np, rho, zk, vrho, v2rho2)

  end subroutine xc_f03_lda_exc_vxc_fxc

  subroutine xc_f03_lda_exc_vxc_fxc_kxc(p, np, rho, zk, vrho, v2rho2, v3rho3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), v2rho2(*), v3rho3(*)

    call xc_lda_exc_vxc_fxc_kxc(p%ptr, np, rho, zk, vrho, v2rho2, v3rho3)

  end subroutine xc_f03_lda_exc_vxc_fxc_kxc

  subroutine xc_f03_lda_vxc(p, np, rho, vrho)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: vrho(*)

    call xc_lda_vxc(p%ptr, np, rho, vrho)

  end subroutine xc_f03_lda_vxc

  subroutine xc_f03_lda_vxc_fxc(p, np, rho, vrho, v2rho2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: vrho(*), v2rho2(*)

    call xc_lda_vxc_fxc(p%ptr, np, rho, vrho, v2rho2)

  end subroutine xc_f03_lda_vxc_fxc

  subroutine xc_f03_lda_vxc_fxc_kxc(p, np, rho, vrho, v2rho2, v3rho3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: vrho(*), v2rho2(*), v3rho3(*)

    call xc_lda_vxc_fxc_kxc(p%ptr, np, rho, vrho, v2rho2, v3rho3)

  end subroutine xc_f03_lda_vxc_fxc_kxc

  subroutine xc_f03_lda_fxc(p, np, rho, v2rho2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: v2rho2(*)

    call xc_lda_fxc(p%ptr, np, rho, v2rho2)

  end subroutine xc_f03_lda_fxc

  subroutine xc_f03_lda_kxc(p, np, rho, v3rho3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: v3rho3(*)

    call xc_lda_kxc(p%ptr, np, rho, v3rho3)

  end subroutine xc_f03_lda_kxc

  subroutine xc_f03_lda_lxc(p, np, rho, v4rho4)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*)
    real(c_double),      intent(out) :: v4rho4(*)

    call xc_lda_lxc(p%ptr, np, rho, v4rho4)

  end subroutine xc_f03_lda_lxc

  ! GGAs
  !----------------------------------------------------------------
  subroutine xc_f03_gga(p, np, rho, sigma, zk, vrho, vsigma,    &
       v2rho2, v2rhosigma, v2sigma2,                            &
       v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,              &
       v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4 &
       )
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)
    real(c_double),      intent(out) :: v4rho4(*), v4rho3sigma(*), v4rho2sigma2(*), v4rhosigma3(*), v4sigma4(*)

    call xc_gga(p%ptr, np, rho, sigma, zk, vrho, vsigma,          &
         v2rho2, v2rhosigma, v2sigma2,                            &
         v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3,              &
         v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4 &
         )

  end subroutine xc_f03_gga

  subroutine xc_f03_gga_exc(p, np, rho, sigma, zk)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: zk(*)

    call xc_gga_exc(p%ptr, np, rho, sigma, zk)

  end subroutine xc_f03_gga_exc

  subroutine xc_f03_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*)

    call xc_gga_exc_vxc(p%ptr, np, rho, sigma, zk, vrho, vsigma)

  end subroutine xc_f03_gga_exc_vxc

  subroutine xc_f03_gga_exc_vxc_fxc(p, np, rho, sigma, zk, vrho, vsigma,    &
       v2rho2, v2rhosigma, v2sigma2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)

    call xc_gga_exc_vxc_fxc(p%ptr, np, rho, sigma, zk, vrho, vsigma,          &
         v2rho2, v2rhosigma, v2sigma2)

  end subroutine xc_f03_gga_exc_vxc_fxc

  subroutine xc_f03_gga_exc_vxc_fxc_kxc(p, np, rho, sigma, zk, vrho, vsigma,    &
       v2rho2, v2rhosigma, v2sigma2,                            &
       v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)

    call xc_gga_exc_vxc_fxc_kxc(p%ptr, np, rho, sigma, zk, vrho, vsigma,          &
         v2rho2, v2rhosigma, v2sigma2,                            &
         v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)

  end subroutine xc_f03_gga_exc_vxc_fxc_kxc

  subroutine xc_f03_gga_vxc(p, np, rho, sigma, vrho, vsigma)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: vrho(*), vsigma(*)

    call xc_gga_vxc(p%ptr, np, rho, sigma, vrho, vsigma)

  end subroutine xc_f03_gga_vxc

  subroutine xc_f03_gga_vxc_fxc(p, np, rho, sigma, vrho, vsigma,    &
       v2rho2, v2rhosigma, v2sigma2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: vrho(*), vsigma(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)

    call xc_gga_vxc_fxc(p%ptr, np, rho, sigma, vrho, vsigma, &
         v2rho2, v2rhosigma, v2sigma2)

  end subroutine xc_f03_gga_vxc_fxc

  subroutine xc_f03_gga_vxc_fxc_kxc(p, np, rho, sigma, vrho, vsigma,    &
       v2rho2, v2rhosigma, v2sigma2,                            &
       v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: vrho(*), vsigma(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)

    call xc_gga_vxc_fxc_kxc(p%ptr, np, rho, sigma, vrho, vsigma,  &
         v2rho2, v2rhosigma, v2sigma2,                            &
         v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)

  end subroutine xc_f03_gga_vxc_fxc_kxc

  subroutine xc_f03_gga_fxc(p, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2sigma2(*)

    call xc_gga_fxc(p%ptr, np, rho, sigma, v2rho2, v2rhosigma, v2sigma2)

  end subroutine xc_f03_gga_fxc

  subroutine xc_f03_gga_kxc(p, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rhosigma2(*), v3sigma3(*)

    call xc_gga_kxc(p%ptr, np, rho, sigma, v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3)

  end subroutine xc_f03_gga_kxc

  subroutine xc_f03_gga_lxc(p, np, rho, sigma, v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*)
    real(c_double),      intent(out) :: v4rho4(*), v4rho3sigma(*), v4rho2sigma2(*), v4rhosigma3(*), v4sigma4(*)

    call xc_gga_lxc(p%ptr, np, rho, sigma, v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4)
  end subroutine xc_f03_gga_lxc

  real(c_double) function xc_f03_gga_ak13_get_asymptotic(homo) result(asymptotic)
    real(c_double), intent(in) :: homo

    asymptotic = xc_gga_ak13_get_asymptotic(homo)

  end function xc_f03_gga_ak13_get_asymptotic

  integer(c_int) function xc_f03_hyb_type(p) result(type)
    type(xc_f03_func_t), intent(in) :: p

    type = xc_hyb_type(p%ptr)

  end function xc_f03_hyb_type

  real(c_double) function xc_f03_hyb_exx_coef(p) result(coef)
    type(xc_f03_func_t), intent(in) :: p

    coef = xc_hyb_exx_coef(p%ptr)

  end function xc_f03_hyb_exx_coef

  subroutine xc_f03_hyb_cam_coef(p, omega, alpha, beta)
    type(xc_f03_func_t), intent(in)  :: p
    real(c_double),       intent(out) :: omega, alpha, beta

    call xc_hyb_cam_coef(p%ptr, omega, alpha, beta)

  end subroutine xc_f03_hyb_cam_coef

  subroutine xc_f03_nlc_coef(p, nlc_b, nlc_c)
    type(xc_f03_func_t), intent(in)  :: p
    real(c_double),       intent(out) :: nlc_b, nlc_c

    call xc_nlc_coef(p%ptr, nlc_b, nlc_c)

  end subroutine xc_f03_nlc_coef


  ! the meta-GGAs
  !----------------------------------------------------------------
  subroutine xc_f03_mgga(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,   &
       v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
       v2lapl2, v2lapltau, v2tau2,                                                    &
       v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
       v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
       v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
       v3lapltau2, v3tau3,                                                            &
       v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho2sigma2, v4rho2sigmalapl,     &
       v4rho2sigmatau, v4rho2lapl2, v4rho2lapltau, v4rho2tau2, v4rhosigma3,           &
       v4rhosigma2lapl, v4rhosigma2tau, v4rhosigmalapl2, v4rhosigmalapltau,           &
       v4rhosigmatau2, v4rholapl3, v4rholapl2tau, v4rholapltau2, v4rhotau3, v4sigma4, &
       v4sigma3lapl, v4sigma3tau, v4sigma2lapl2, v4sigma2lapltau, v4sigma2tau2,       &
       v4sigmalapl3, v4sigmalapl2tau, v4sigmalapltau2, v4sigmatau3, v4lapl4,          &
       v4lapl3tau, v4lapl2tau2, v4lapltau3, v4tau4 &
       )
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
         v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
         v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),               &
         v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),      &
         v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),    &
         v3lapltau2(*), v3tau3(*)
    real(c_double),    intent(out) :: &
         v4rho4(*), v4rho3sigma(*), v4rho3lapl(*), v4rho3tau(*), v4rho2sigma2(*), v4rho2sigmalapl(*),     &
         v4rho2sigmatau(*), v4rho2lapl2(*), v4rho2lapltau(*), v4rho2tau2(*), v4rhosigma3(*),              &
         v4rhosigma2lapl(*), v4rhosigma2tau(*), v4rhosigmalapl2(*), v4rhosigmalapltau(*),                 &
         v4rhosigmatau2(*), v4rholapl3(*), v4rholapl2tau(*), v4rholapltau2(*), v4rhotau3(*), v4sigma4(*), &
         v4sigma3lapl(*), v4sigma3tau(*), v4sigma2lapl2(*), v4sigma2lapltau(*), v4sigma2tau2(*),          &
         v4sigmalapl3(*), v4sigmalapl2tau(*), v4sigmalapltau2(*), v4sigmatau3(*), v4lapl4(*),             &
         v4lapl3tau(*), v4lapl2tau2(*), v4lapltau3(*), v4tau4(*)

    call xc_mgga(p%ptr, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,       &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2,                                                    &
         v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
         v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
         v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
         v3lapltau2, v3tau3,                                                            &
         v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho2sigma2, v4rho2sigmalapl,     &
         v4rho2sigmatau, v4rho2lapl2, v4rho2lapltau, v4rho2tau2, v4rhosigma3,           &
         v4rhosigma2lapl, v4rhosigma2tau, v4rhosigmalapl2, v4rhosigmalapltau,           &
         v4rhosigmatau2, v4rholapl3, v4rholapl2tau, v4rholapltau2, v4rhotau3, v4sigma4, &
         v4sigma3lapl, v4sigma3tau, v4sigma2lapl2, v4sigma2lapltau, v4sigma2tau2,       &
         v4sigmalapl3, v4sigmalapl2tau, v4sigmalapltau2, v4sigmatau3, v4lapl4,          &
         v4lapl3tau, v4lapl2tau2, v4lapltau3, v4tau4 &
         )
  end subroutine xc_f03_mgga

  subroutine xc_f03_mgga_exc(p, np, rho, sigma, lapl, tau, zk)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: zk(*)

    call xc_mgga_exc(p%ptr, np, rho, sigma, lapl, tau, zk)

  end subroutine xc_f03_mgga_exc

  subroutine xc_f03_mgga_exc_vxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)

    call xc_mgga_exc_vxc(p%ptr, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau)

  end subroutine xc_f03_mgga_exc_vxc

  subroutine xc_f03_mgga_exc_vxc_fxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,   &
       v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
       v2lapl2, v2lapltau, v2tau2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
         v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)

    call xc_mgga_exc_vxc_fxc(p%ptr, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,       &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2)
  end subroutine xc_f03_mgga_exc_vxc_fxc

  subroutine xc_f03_mgga_exc_vxc_fxc_kxc(p, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,   &
       v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
       v2lapl2, v2lapltau, v2tau2,                                                    &
       v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
       v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
       v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
       v3lapltau2, v3tau3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: zk(*), vrho(*), vsigma(*), vlapl(*), vtau(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
         v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
         v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),               &
         v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),      &
         v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),    &
         v3lapltau2(*), v3tau3(*)

    call xc_mgga_exc_vxc_fxc_kxc(p%ptr, np, rho, sigma, lapl, tau, zk, vrho, vsigma, vlapl, vtau,       &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2,                                                    &
         v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
         v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
         v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
         v3lapltau2, v3tau3)
  end subroutine xc_f03_mgga_exc_vxc_fxc_kxc

  subroutine xc_f03_mgga_vxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)

    call xc_mgga_vxc(p%ptr, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)

  end subroutine xc_f03_mgga_vxc

  subroutine xc_f03_mgga_vxc_fxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau,   &
       v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
       v2lapl2, v2lapltau, v2tau2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
         v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)

    call xc_mgga_vxc_fxc(p%ptr, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau,       &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2)
  end subroutine xc_f03_mgga_vxc_fxc

  subroutine xc_f03_mgga_vxc_fxc_kxc(p, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau,   &
       v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
       v2lapl2, v2lapltau, v2tau2,                                                    &
       v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
       v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
       v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
       v3lapltau2, v3tau3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: vrho(*), vsigma(*), vlapl(*), vtau(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
         v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
         v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),               &
         v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),      &
         v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),    &
         v3lapltau2(*), v3tau3(*)

    call xc_mgga_vxc_fxc_kxc(p%ptr, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau,       &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
         v2lapl2, v2lapltau, v2tau2,                                                    &
         v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
         v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
         v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
         v3lapltau2, v3tau3)
  end subroutine xc_f03_mgga_vxc_fxc_kxc

  subroutine xc_f03_mgga_fxc(p, np, rho, sigma, lapl, tau, &
       v2rho2, v2rhosigma, v2rholapl, v2rhotau, v2sigma2, v2sigmalapl, v2sigmatau,    &
       v2lapl2, v2lapltau, v2tau2)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: v2rho2(*), v2rhosigma(*), v2rholapl(*), v2rhotau(*),    &
         v2sigma2(*), v2sigmalapl(*), v2sigmatau(*), v2lapl2(*), v2lapltau(*), v2tau2(*)

    call xc_mgga_fxc(p%ptr, np, rho, sigma, lapl, tau,   &
         v2rho2, v2rhosigma, v2rholapl, v2rhotau,           &
         v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2)

  end subroutine xc_f03_mgga_fxc

  subroutine xc_f03_mgga_kxc(p, np, rho, sigma, lapl, tau, &
       v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
       v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
       v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
       v3lapltau2, v3tau3)
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: v3rho3(*), v3rho2sigma(*), v3rho2lapl(*), v3rho2tau(*), &
         v3rhosigma2(*), v3rhosigmalapl(*), v3rhosigmatau(*), v3rholapl2(*),               &
         v3rholapltau(*), v3rhotau2(*), v3sigma3(*), v3sigma2lapl(*), v3sigma2tau(*),      &
         v3sigmalapl2(*), v3sigmalapltau(*), v3sigmatau2(*), v3lapl3(*), v3lapl2tau(*),    &
         v3lapltau2(*), v3tau3(*)

    call xc_mgga_kxc(p%ptr, np, rho, sigma, lapl, tau,  &
         v3rho3, v3rho2sigma, v3rho2lapl, v3rho2tau, v3rhosigma2, v3rhosigmalapl,       &
         v3rhosigmatau, v3rholapl2, v3rholapltau, v3rhotau2, v3sigma3, v3sigma2lapl,    &
         v3sigma2tau, v3sigmalapl2, v3sigmalapltau, v3sigmatau2, v3lapl3, v3lapl2tau,   &
         v3lapltau2, v3tau3)
  end subroutine xc_f03_mgga_kxc

  subroutine xc_f03_mgga_lxc(p, np, rho, sigma, lapl, tau, &
       v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho2sigma2, v4rho2sigmalapl,     &
       v4rho2sigmatau, v4rho2lapl2, v4rho2lapltau, v4rho2tau2, v4rhosigma3,           &
       v4rhosigma2lapl, v4rhosigma2tau, v4rhosigmalapl2, v4rhosigmalapltau,           &
       v4rhosigmatau2, v4rholapl3, v4rholapl2tau, v4rholapltau2, v4rhotau3, v4sigma4, &
       v4sigma3lapl, v4sigma3tau, v4sigma2lapl2, v4sigma2lapltau, v4sigma2tau2,       &
       v4sigmalapl3, v4sigmalapl2tau, v4sigmalapltau2, v4sigmatau3, v4lapl4,          &
       v4lapl3tau, v4lapl2tau2, v4lapltau3, v4tau4 &
       )
    type(xc_f03_func_t), intent(in)  :: p
    integer(c_size_t),   intent(in)  :: np
    real(c_double),      intent(in)  :: rho(*), sigma(*), lapl(*), tau(*)
    real(c_double),      intent(out) :: &
         v4rho4(*), v4rho3sigma(*), v4rho3lapl(*), v4rho3tau(*), v4rho2sigma2(*), v4rho2sigmalapl(*),     &
         v4rho2sigmatau(*), v4rho2lapl2(*), v4rho2lapltau(*), v4rho2tau2(*), v4rhosigma3(*),              &
         v4rhosigma2lapl(*), v4rhosigma2tau(*), v4rhosigmalapl2(*), v4rhosigmalapltau(*),                 &
         v4rhosigmatau2(*), v4rholapl3(*), v4rholapl2tau(*), v4rholapltau2(*), v4rhotau3(*), v4sigma4(*), &
         v4sigma3lapl(*), v4sigma3tau(*), v4sigma2lapl2(*), v4sigma2lapltau(*), v4sigma2tau2(*),          &
         v4sigmalapl3(*), v4sigmalapl2tau(*), v4sigmalapltau2(*), v4sigmatau3(*), v4lapl4(*),             &
         v4lapl3tau(*), v4lapl2tau2(*), v4lapltau3(*), v4tau4(*)

    call xc_mgga_lxc(p%ptr, np, rho, sigma, lapl, tau,  &
         v4rho4, v4rho3sigma, v4rho3lapl, v4rho3tau, v4rho2sigma2, v4rho2sigmalapl,     &
         v4rho2sigmatau, v4rho2lapl2, v4rho2lapltau, v4rho2tau2, v4rhosigma3,           &
         v4rhosigma2lapl, v4rhosigma2tau, v4rhosigmalapl2, v4rhosigmalapltau,           &
         v4rhosigmatau2, v4rholapl3, v4rholapl2tau, v4rholapltau2, v4rhotau3, v4sigma4, &
         v4sigma3lapl, v4sigma3tau, v4sigma2lapl2, v4sigma2lapltau, v4sigma2tau2,       &
         v4sigmalapl3, v4sigmalapl2tau, v4sigmalapltau2, v4sigmatau3, v4lapl4,          &
         v4lapl3tau, v4lapl2tau2, v4lapltau3, v4tau4 &
         )
  end subroutine xc_f03_mgga_lxc


  ! Helper functions to convert between C and Fortran strings
  ! Based on the routines by Joseph M. Krahn
  function f_to_c_string(f_string) result(c_string)
    character(len=*), intent(in) :: f_string
    character(kind=c_char,len=1) :: c_string(len_trim(f_string)+1)

    integer :: i, strlen

    strlen = len_trim(f_string)

    forall (i=1:strlen)
       c_string(i) = f_string(i:i)
    end forall
    c_string(strlen+1) = C_NULL_CHAR

  end function f_to_c_string

  subroutine c_to_f_string(c_string, f_string)
    character(kind=c_char,len=1), intent(in)  :: c_string(*)
    character(len=*),             intent(out) :: f_string

    integer :: i

    i = 1
    do while(c_string(i) /= C_NULL_CHAR .and. i <= len(f_string))
       f_string(i:i) = c_string(i)
       i = i + 1
    end do
    if (i < len(f_string)) f_string(i:) = ' '

  end subroutine c_to_f_string

  subroutine c_to_f_string_ptr(c_string, f_string)
    type(c_ptr),      intent(in)  :: c_string
    character(len=*), intent(out) :: f_string

    character(len=1, kind=c_char), pointer :: p_chars(:)
    integer :: i

    if (.not. c_associated(c_string)) then
       f_string = ' '
    else
       call c_f_pointer(c_string, p_chars, [huge(0)])
       i = 1
       do while(p_chars(i) /= C_NULL_CHAR .and. i <= len(f_string))
          f_string(i:i) = p_chars(i)
          i = i + 1
       end do
       if (i < len(f_string)) f_string(i:) = ' '
    end if

  end subroutine c_to_f_string_ptr

end module LibxcInterface_

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

