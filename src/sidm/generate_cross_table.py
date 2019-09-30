#script to generate DM self-interaction cross sections for different models:
#writes out file(s) with two columns: rel. velocity (km/s) and cross section per mass (cm^2/g)
#this file is then use to set up internally the cross section tables

import numpy as np

#this should not be changed
CROSS_VBINS                = 1000000 
CROSS_VMAX                 = 1e5
CROSS_VMIN                 = 1e-2
CrossTable                 = np.zeros(CROSS_VBINS)
VelTable                   = np.zeros(CROSS_VBINS)
Dvlog                      = np.log(CROSS_VMAX / CROSS_VMIN) / CROSS_VBINS

#constants
light_speed = 299972.0                           # [km/s]
conv_factor = (1.97327e-17)**2 / 1.782662e-21    # (TeV^-2 ) / TeV   --->   cm^2 / gr



#select model here (note that model specific parameters must be set below)
MODEL=3#  2:INSIDM  3:ETHOS-4    


###################################################
#MODEL 0:                                         #
#constan cross section                            #
###################################################
if (MODEL==0):
        CrossSectionPerMass_in_cgs = 10.0/2.0
        for i in range(0, CROSS_VBINS):
                rel_vel    = np.exp(Dvlog * (i + 0.5) + np.log(CROSS_VMIN))
		VelTable[i]   = rel_vel
                CrossTable[i] = CrossSectionPerMass_in_cgs #cm^2/g 

        data = np.vstack((VelTable, CrossTable)).T
        np.savetxt("sidm_cross_reaction_0.txt", data,delimiter="\t")


########################################################################
#MODEL 1:                                                              #
#plasma-physics inspired formula (Loeb and Wiener) effective parameters#
########################################################################
if (MODEL==1):
	PeakSigma                  = 1.0
	CrossSectionPerMass_in_cgs = 1.0
	for i in range(0, CROSS_VBINS):
		rel_vel    = np.exp(Dvlog * (i + 0.5) + np.log(CROSS_VMIN))
		beta       = np.pi * (PeakSigma/rel_vel) * (PeakSigma/rel_vel)

		sigma_norm = 0.0

		if (beta < 0.1):
			sigma_norm = (4.0 * np.pi / 22.7) * beta * beta * np.log(1.0 + 1.0 / beta)
		if (beta < 1e3) & (beta > 0.1):
			sigma_norm = (8.0 * np.pi / 22.7) * beta * beta / (1.0 + 1.5 * beta**1.65)
		if (beta > 1e3):
			sigma_norm = (np.pi / 22.7) * (np.log(beta) + 1.0 - 0.5 / np.log(beta))**2.0

		VelTable[i]   = rel_vel                                  #km/s
		CrossTable[i] = CrossSectionPerMass_in_cgs * sigma_norm  #cm^2/g 

	data = np.vstack((VelTable, CrossTable)).T
	np.savetxt("sidm_cross_reaction_0.txt", data,delimiter="\t")


###############################################################################
#MODEL 2:                                                                     #
#plasma-physics inspired formula (Loeb and Wiener) particle physics parameters#
###############################################################################
if (MODEL==2):
	alphac = 0.1          #1
	mphi   = 30.0/1000.0  #mass*c^2 in GeV
        mchi   = 10.0         #mass*c^2 in GeV
        for i in range(0, CROSS_VBINS):
                rel_vel    = np.exp(Dvlog * (i + 0.5) + np.log(CROSS_VMIN))
                beta       = 2 * alphac * mphi / (mchi * rel_vel/light_speed ** 2)
		#note: 1GeV^-2=0.389 * 10**(-27) cm^2, 1GeV=1.8 * 10**(-24) see natural units http://www.phys.ufl.edu/~korytov/phz4390/note_01_NaturalUnits_SMsummary.pdf 
		sigma_max  = (22.7 / (mphi**2) * 0.389 * 10**(-27)) / (mchi * 1.8 * 10**(-24))      #cm^2/g 
                sigma_norm = 0.0

                if (beta < 0.1):
                        sigma_norm = (4.0 * np.pi / 22.7) * beta * beta * np.log(1.0 + 1.0 / beta)
                if (beta < 1e3) & (beta > 0.1):
                        sigma_norm = (8.0 * np.pi / 22.7) * beta * beta / (1.0 + 1.5 * beta**1.65)
                if (beta > 1e3):
                        sigma_norm = (np.pi / 22.7) * (np.log(beta) + 1.0 - 0.5 / np.log(beta))**2.0

                VelTable[i]   = rel_vel                                 #km/s
                CrossTable[i] = sigma_max * sigma_norm  #cm^2/g 

        data = np.vstack((VelTable, CrossTable)).T
        np.savetxt("sidm_cross_reaction_0.txt", data,delimiter="\t")




#################################################################################################################################################
#MODEL 3:                                                                                                                                       #
#effective cross section formulae SIDM_EFF=1 (based on Torsten's fits for attractive and repulsive potentials in the classical regime: 22/10/14)#
#################################################################################################################################################
if (MODEL==3):
	SIDM_EFF_alpha_chi   =  0.5
	SIDM_EFF_m_chi       =  3.7 
	SIDM_EFF_m_phi       =  5e-6 

	for i in range(0, CROSS_VBINS):
		rel_vel = np.exp(Dvlog * (i + 0.5) + np.log(CROSS_VMIN))    # km/s
		rel_vel_nat = rel_vel / light_speed			    # dimensioness

		beta = 2.0 * SIDM_EFF_alpha_chi * SIDM_EFF_m_phi / SIDM_EFF_m_chi / pow(rel_vel_nat, 2.0);    # dimensionless

		if(beta < 0.01):
			sigma_norm = 2.0 * beta**2 * np.log(1. + (beta**-2.0))
		if((beta < 1.0e4) & (beta >= 0.01)):
        		sigma_norm = 8.0 *(beta**1.8) / (1. + 5 * (beta**0.9) + 0.85 *(beta**1.6))
		if(beta >= 1.0e4):
        		sigma_norm = (np.log(2.0 * beta) - np.log(np.log(2.0 * beta)))**2.0

		sigma_norm = np.pi * (SIDM_EFF_m_phi**(-2.0) / SIDM_EFF_m_chi )* sigma_norm      # TeV^-3

		if (beta < 0.01):
        		sigma_norm_att = 2.0 * (beta**2.0) * np.log(1. + (beta**(-2.0)))
		if((beta < 1.0e2) & (beta >= 0.01)):
        		sigma_norm_att = 7.0 * ((beta**1.8) + 280.0 * (beta / 10.0)**10.3) / (1. + 1.4 * beta + 0.006 * (beta**4.0) + 160.0 * (beta / 10.0)**10.0)
		if(beta >= 1.0e2):
			sigma_norm_att = 0.81 * (1.0 + np.log(beta) - 0.5 / np.log(beta))**2.0

		sigma_norm_att = np.pi * (SIDM_EFF_m_phi**(-2.0)) / SIDM_EFF_m_chi * sigma_norm_att      # TeV^-3

		VelTable[i]   = rel_vel                                           #km/s
		CrossTable[i] = conv_factor * 0.5 * (sigma_norm + sigma_norm_att) # cm^2/gr 

	data = np.vstack((VelTable, CrossTable)).T
	np.savetxt("sidm_cross_reaction_0.txt", data,delimiter="\t")



#################################################################################################################################################
#MODEL 4:                                                                                                                                       #
#effective cross section formulae SIDM_EFF=2 (based on Torsten's fits for attractive and repulsive potentials in the classical regime: 22/10/14)#
#################################################################################################################################################
if (MODEL==4):
        SIDM_EFF_alpha_chi   =  0.5
        SIDM_EFF_m_chi       =  3.7
        SIDM_EFF_m_phi       =  5e-6

        for i in range(0, CROSS_VBINS):
                rel_vel = np.exp(Dvlog * (i + 0.5) + np.log(CROSS_VMIN))    # km/s
                rel_vel_nat = rel_vel / light_speed                         # dimensioness

                beta = 2.0 * SIDM_EFF_alpha_chi * SIDM_EFF_m_phi / SIDM_EFF_m_chi / pow(rel_vel_nat, 2.0);    # dimensionless

                if(beta < 0.01):
                        sigma_norm = 2.0 * beta**2 * np.log(1. + (beta**-2.0))
                if((beta < 1.0e4) & (beta >= 0.01)):
                        sigma_norm = 8.0 *(beta**1.8) / (1. + 5 * (beta**0.9) + 0.85 *(beta**1.6))
                if(beta >= 1.0e4):
                        sigma_norm = (np.log(2.0 * beta) - np.log(np.log(2.0 * beta)))**2.0

                sigma_norm = np.pi * (SIDM_EFF_m_phi**(-2.0) / SIDM_EFF_m_chi )* sigma_norm      # TeV^-3

		sigma_norm_att = sigma_norm 

                VelTable[i]   = rel_vel                                           #km/s
                CrossTable[i] = conv_factor * 0.5 * (sigma_norm + sigma_norm_att) # cm^2/gr 

        data = np.vstack((VelTable, CrossTable)).T
        np.savetxt("sidm_cross_reaction_0.txt", data,delimiter="\t")




##############################
#MODEL 5:                    #
#formula from Francis-Yan ADM#
##############################
if (MODEL==5):
	light_speed = 299972.0  # [km/s]
	m_DA = 1.03	        # dark atom mass [TeV]
	alpha_DA = 0.009      	# fine-structure constant
	B_DA = 1e-9           	# binding energy [TeV]

	lambd=1e12


        for i in range(0, CROSS_VBINS):
                rel_vel    = np.exp(Dvlog * (i + 0.5) + np.log(CROSS_VMIN))  #km/s
		rel_vel_nat = rel_vel / light_speed                          #dimensionless

		yy = m_DA * (rel_vel_nat**2.0) / (8. * B_DA);
		bb = 400. / (np.log(1. + (lambd * m_DA / 2. * (rel_vel_nat**2.0) / alpha_DA)**2.0))
	
		sigma_norm = 50. * np.pi * (alpha_DA**2.0) / ((B_DA**2.0) * m_DA) / (1. + 50. * yy + bb * (yy**2.0));   # TeV^-3

                VelTable[i]   = rel_vel                         #km/s
                CrossTable[i] = conv_factor * 0.5 * sigma_norm  # cm^2/gr 

        data = np.vstack((VelTable, CrossTable)).T
        np.savetxt("sidm_cross_reaction_0.txt", data,delimiter="\t")



