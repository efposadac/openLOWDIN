
&LowdinParameters
/

&InputTasks
	InputTasks_method = "RHF"
	InputTasks_effectiveCorePotentials = T
/

&InputParticle
	InputParticle_name = "e-[Na]"
	InputParticle_basisSetName = "LANL2DZ_ECP"
	InputParticle_origin = 0.000000000000E+00 0.000000000000E+00 -6.657500000000E-02 
/

&InputParticle
	InputParticle_name = "e-[H]"
	InputParticle_basisSetName = "STO-3G"
	InputParticle_addParticles = 1
	InputParticle_origin = 0.000000000000E+00 7.541750000000E-01 5.283810000000E-01 
/

&InputParticle
	InputParticle_name = "Na"
	InputParticle_basisSetName = "ECP"
	InputParticle_origin = 0.000000000000E+00 0.000000000000E+00 -6.657500000000E-02 
/

&InputParticle
	InputParticle_name = "H"
	InputParticle_basisSetName = "dirac"
	InputParticle_origin = 0.000000000000E+00 7.541750000000E-01 5.283810000000E-01 
/

&InputSystem
	InputSystem_numberOfParticles = 4
	InputSystem_numberOfExternalPots =  0
	InputSystem_numberOfInterPots =  0
	InputSystem_numberOfOutputs =  0
	InputSystem_description = 'NaH con seudopotenciales'

/
