#!/usr/bin/python
'''
    This script takes in the weight enrichment for the fuel's uranium and the
    smeared density percentage then proceeds to find the nuclide densities for
    U-235, U-238, Zr-90, Zr-91, Zr-92, Zr-94, Zr-96 and Na-23. The fuel assumed 
    is a metallic U-10Zr, an alloy of uranium and zirconium that is 10 percent
    Zr by weight with pure sodium bond.
'''
# Begin calculation ############################################################

# Poll for uranium weight enrichment from user 
w_enrichment = float(input('Enter fuel\'s % weight enrichment: '))
w_enrichment = w_enrichment/100;

# Atomic masses (amu)
Na = 6.022*pow(10,23) # Avagadro's number
m_u235 = 235.0439299
m_u238 = 238.0507883
m_zr = 91.224
# Mass fraction of alloy that is U and Zr
w_u = 0.9;
w_zr = 0.1;
# Find densities of U-235 and U-238 metal by first finding the number density of 
# uranium metal
N_u = 19.1*6.022*pow(10,23)/(238.02891)
p_u235 = m_u235*N_u/Na
p_u238 = m_u238*N_u/Na
# Density of zirconium
p_zr = 6.51

# Find atomic percent of uranium that is U-235
a_enrichment = (w_enrichment*m_u238)/(((1-w_enrichment)*m_u235)+ \
               (w_enrichment*m_u238))

# Find atomic percent of alloy that is U and Zr
x_u = (w_u/(a_enrichment*m_u235+(1-a_enrichment)*m_u238))/(w_zr/m_zr+w_u/ \
      (a_enrichment*m_u235+(1-a_enrichment)*m_u238))
x_zr = (w_zr/m_zr)/(w_zr/m_zr+w_u/(a_enrichment*m_u235+(1-a_enrichment)*m_u238))

# Find atomic percent of alloy that is U-235 and U_238
x_u235 = a_enrichment*x_u
x_u238 = (1-a_enrichment)*x_u

# Find mass-density of alloy, molar mass of alloy, and then number density of 
# alloy
p_alloy = x_zr*p_zr+x_u235*p_u235+x_u238*p_u238
m_alloy = x_zr*m_zr+x_u235*m_u235+x_u238*m_u238
n_alloy = p_alloy*Na/m_alloy

# Find number densities of fuel nuclides in barn-cm
n_u235 = n_alloy*x_u235*pow(10,-24)
n_u238 = n_alloy*x_u238*pow(10,-24)
n_zr = n_alloy*x_zr*pow(10,-24)
# Hardcoded natural abundances
n_zr90 = n_zr*0.5145
n_zr91 = n_zr*0.1122
n_zr92 = n_zr*0.1715
n_zr94 = n_zr*0.1738
n_zr96 = n_zr*0.0280

# Output calculated values #####################################################

# Open new file
import os.path
from os import path
out_file_name = 'run_densities.txt'
i = 1
while path.exists(out_file_name):
  out_file_name = 'run_densities_'+str(i)+'.txt'
  i = i+1
out_file = open(out_file_name,'w')

# Print results to output file, write file explicity for error checking
out_file.write('Run Characteristics\n')
out_file.write('Enrichment Weight Percent = '+str(w_enrichment)+' \n')
out_file.write('Enrichment Weight Percent Split (Outer/Inner) = '+str(1)+' \n')
out_file.write('Smeared Density Percent = '+str(75)+' \n')
out_file.write('U-10Zr alloy density = '+str(round(p_alloy,5))+' g/cm^3 \n')
out_file.write('MCC3ID     Density [#/barn-cm]\n')
out_file.write('U235_7     '+str('{:1.5E}'.format(n_u235))+'\n')
out_file.write('U238_7     '+str('{:1.5E}'.format(n_u238))+'\n')
out_file.write('ZR90_7     '+str('{:1.5E}'.format(n_zr90))+'\n')
out_file.write('ZR91_7     '+str('{:1.5E}'.format(n_zr91))+'\n')
out_file.write('ZR92_7     '+str('{:1.5E}'.format(n_zr92))+'\n')
out_file.write('ZR94_7     '+str('{:1.5E}'.format(n_zr94))+'\n')
out_file.write('ZR96_7     '+str('{:1.5E}'.format(n_zr96))+'\n')
out_file.close()

print('Nuclide Densities Printed To '+out_file_name)