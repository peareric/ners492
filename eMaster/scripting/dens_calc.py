# This script is the Python equivalent of Prof. Yang's nuclide density calculating spreadsheet
import math

# This helper function calculates nuclide densities at both 20C and 600C
def n_density_solver(wpct, mass, wfrac, density, alpha, delta_t, A_num, string):
    n_20c = A_num*density*wfrac/100*(wpct/100)/mass
    n_600c = n_20c/pow((1+alpha*delta_t),3)
    return [string, n_20c, n_600c]


# Inputs
# Enter fuel enrichment in weight percent
enr = float(input('Enter fuel enrichment in weight percent: '))
sd_roomtemperature = float(input('Enter fuel smeared density at room temperature in percent: '))
sd_hot = float(input('Enter fuel smeared density at operating temperature in percent: '))


A_num = 6.0221367*10**(23) # Avagadro's number

##################################################################################
####                     Part 1: Material nuclide densities                   ####
##################################################################################

#########################################################
####                    Fuel Data                    ####
#########################################################


# Fuel constants for U-10Zr
p_fuel = 16.245 # fuel density (g/cc)
alpha_fuel = 1.32*10**(-5) # Thermal expansion coefficient (1/C)
u_wfrac = 90 # Fraction of fuel that is uranium by weight
zr_wfrac = 10


# Isotopic mass data for U-10Zr
u235_mass = 235.043940317917
u238_mass = 238.050770432513
zr90_mass = 89.9047247588784
zr91_mass = 90.9056229550252
zr92_mass = 91.9050081537980
zr94_mass = 93.9063002136336
zr96_mass = 95.9082983389104


# Weight percents for U-10Zr
u235_wpct = enr
u238_wpct = 100 - u235_wpct
zr90_wpct = 50.70613
zr91_wpct = 11.18088
zr92_wpct = 17.27810
zr94_wpct = 17.89110
zr96_wpct = 2.94379


#########################################################
####                     HT9 Data                    ####
#########################################################


p_clad = 7.76
alpha_clad = 1.39*10**(-5)
fe_wfrac = 84.45
ni_wfrac = 0.6
cr_wfrac = 11.95
mo_wfrac = 1
mn_wfrac = 0.6
c_wfrac  = 0.2
si_wfrac = 0.38
w_wfrac  = 0.52
v_wfrac  = 0.3


# Isotopic masses in HT9
fe54_mass  = 53.939365048016
fe56_mass  = 55.934504251864
fe57_mass  = 56.935099848536
fe58_mass  = 57.933678115376
ni58_mass  = 57.935695445208
ni60_mass  = 59.930834649056
ni61_mass  = 60.931430245728
ni62_mass  = 61.927991182736
ni64_mass  = 63.928173711164
cr50_mass  = 49.946060645572
cr52_mass  = 51.940191184504
cr53_mass  = 52.940786781176
cr54_mass  = 53.939365048016
mo92_mass  = 91.9068237506468
mo94_mass  = 93.9050898157344
mo95_mass  = 94.9058871453896
mo96_mass  = 95.9046671452128
mo97_mass  = 96.905968807326
mo98_mass  = 97.9053540060988
mo100_mass = 99.907251264884
mn55_mass  = 54.9380441813476
c_mass     = 12.001095170568
si28_mass  = 27.977338775092
si29_mass  = 28.976925706848
si30_mass  = 29.973486643856
w180_mass  = 179.946829679316
w182_mass  = 181.95306419724
w183_mass  = 182.95164246408
w184_mass  = 183.95022073092
w186_mass  = 185.95746391376
v50_mass   = 49.9471628
v51_mass   = 50.9439637


# Weight percents for HT9
fe54_wpct  = 5.64557
fe56_wpct  = 91.90150
fe57_wpct  = 2.16037
fe58_wpct  = 0.29255
ni58_wpct  = 67.19793
ni60_wpct  = 26.77577
ni61_wpct  = 1.18347
ni62_wpct  = 3.83426
ni64_wpct  = 1.00859
cr50_wpct  = 4.17371
cr52_wpct  = 83.69924
cr53_wpct  = 9.67366
cr54_wpct  = 2.45340
mo92_wpct  = 14.21744
mo94_wpct  = 9.05463
mo95_wpct  = 15.74984
mo96_wpct  = 16.67538
mo97_wpct  = 9.64703
mo98_wpct  = 24.62655
mo100_wpct = 10.02913
mn55_wpct  = 100.0
c_wpct     = 100.0
si28_wpct  = 91.87350
si29_wpct  = 4.81816
si30_wpct  = 3.30834
w180_wpct  = 0.11746
w182_wpct  = 26.22547
w183_wpct  = 14.24453
w184_wpct  = 30.65973
w186_wpct  = 28.75281
v50_wpct   = 0.24512
v51_wpct   = 99.75488


#########################################################
####                    Sodium Data                  ####
#########################################################


p_sodium = 0.847856


# Isotopic mass for sodium
na23_mass = 22.989490765472


#########################################################
####                 Temperature Data                ####
#########################################################


t_fuel = 600
t_clad = 550
delta_t_fuel = t_fuel - 20
delta_t_clad = t_clad - 20
t_in = 360
t_out = 510
t_avg = (t_in+t_out)/2


# Solve fuel number densities
u235_n = n_density_solver(u235_wpct, u235_mass, u_wfrac, p_fuel, alpha_fuel, delta_t_fuel, A_num,  'U235_7     ')
u238_n = n_density_solver(u238_wpct, u238_mass, u_wfrac, p_fuel, alpha_fuel, delta_t_fuel, A_num,  'U238_7     ')
zr90_n = n_density_solver(zr90_wpct, zr90_mass, zr_wfrac, p_fuel, alpha_fuel, delta_t_fuel, A_num, 'ZR90_7     ')
zr91_n = n_density_solver(zr91_wpct, zr91_mass, zr_wfrac, p_fuel, alpha_fuel, delta_t_fuel, A_num, 'ZR91_7     ')
zr92_n = n_density_solver(zr92_wpct, zr92_mass, zr_wfrac, p_fuel, alpha_fuel, delta_t_fuel, A_num, 'ZR92_7     ')
zr94_n = n_density_solver(zr94_wpct, zr94_mass, zr_wfrac, p_fuel, alpha_fuel, delta_t_fuel, A_num, 'ZR94_7     ')
zr96_n = n_density_solver(zr96_wpct, zr96_mass, zr_wfrac, p_fuel, alpha_fuel, delta_t_fuel, A_num, 'ZR96_7     ')


# Solve cladding number densities
fe54_n = n_density_solver(fe54_wpct, fe54_mass, fe_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'FE54_7     ')
fe56_n = n_density_solver(fe56_wpct, fe56_mass, fe_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'FE56_7     ')
fe57_n = n_density_solver(fe57_wpct, fe57_mass, fe_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'FE57_7     ')
fe58_n = n_density_solver(fe58_wpct, fe58_mass, fe_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'FE58_7     ')
ni58_n = n_density_solver(ni58_wpct, ni58_mass, ni_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'NI58_7     ')
ni60_n = n_density_solver(ni60_wpct, ni60_mass, ni_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'NI60_7     ')
ni61_n = n_density_solver(ni61_wpct, ni61_mass, ni_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'NI61_7     ')
ni62_n = n_density_solver(ni62_wpct, ni62_mass, ni_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'NI62_7     ')
ni64_n = n_density_solver(ni64_wpct, ni64_mass, ni_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'NI64_7     ')
cr50_n = n_density_solver(cr50_wpct, cr50_mass, cr_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'CR50_7     ')
cr52_n = n_density_solver(cr52_wpct, cr52_mass, cr_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'CR52_7     ')
cr53_n = n_density_solver(cr53_wpct, cr53_mass, cr_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'CR53_7     ')
cr54_n = n_density_solver(cr54_wpct, cr54_mass, cr_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'CR54_7     ')
mo92_n = n_density_solver(mo92_wpct, mo92_mass, mo_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'MO92_7     ')
mo94_n = n_density_solver(mo94_wpct, mo94_mass, mo_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'MO94_7     ')
mo95_n = n_density_solver(mo95_wpct, mo95_mass, mo_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'MO95_7     ')
mo96_n = n_density_solver(mo96_wpct, mo96_mass, mo_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'MO96_7     ')
mo97_n = n_density_solver(mo97_wpct, mo97_mass, mo_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'MO97_7     ')
mo98_n = n_density_solver(mo98_wpct, mo98_mass, mo_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'MO98_7     ')
mo100_n = n_density_solver(mo100_wpct, mo100_mass, mo_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'MO1007     ')
mn55_n = n_density_solver(mn55_wpct, mn55_mass, mn_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'MN55_7     ')
c_n = n_density_solver(c_wpct, c_mass, c_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'C')
si28_n = n_density_solver(si28_wpct, si28_mass, si_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'Si-28')
si29_n = n_density_solver(si29_wpct, si29_mass, si_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'Si-29')
si30_n = n_density_solver(si30_wpct, si30_mass, si_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'Si-30')
w180_n = n_density_solver(w180_wpct, w180_mass, w_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'W-180')
w182_n = n_density_solver(w182_wpct, w182_mass, w_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'W-182')
w183_n = n_density_solver(w183_wpct, w183_mass, w_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'W-183')
w184_n = n_density_solver(w184_wpct, w184_mass, w_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'W-184')
w186_n = n_density_solver(w186_wpct, w186_mass, w_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'W-186')
v50_n = n_density_solver(v50_wpct, v50_mass, v_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'V-50')
v51_n = n_density_solver(v51_wpct, v51_mass, v_wfrac, p_clad, alpha_clad, delta_t_clad, A_num, 'V-51')

# Solve coolant number density
na23_n = ['NA23_7     ', 'N/A', A_num*p_sodium/na23_mass]


##################################################################################
####                          Part 2: Volume fractions                        ####
##################################################################################


# SS316 grid plate data
p_plate = 7.93
alpha_plate = 0.0000107322282132897


# Expansion factors
expand_plate = 1+alpha_plate*(t_in-20)
expand_clad = 1+alpha_clad*(t_clad-20)
expand_fuel = 1.05


#################################################################
####                    Assembly dimensions                  ####
#################################################################
# Constants
gap = 0.4318
hex_pitch = 16.142
wall_thickness = 0.3937
pitch = 0.9934
pin_diameter = 0.8306
pin_to_diameter = pitch/pin_diameter
clad_thickness = 0.0559
wire_diameter = 0.1575
wire_pitch = 20.32


# Calculated dimension values at room temp and operating temp
hexagon_pitch = ['Hexagon Pitch', hex_pitch, hex_pitch*expand_plate]
duct_wall_thickness = ['Duct Wall Thickness', wall_thickness, wall_thickness*(1+alpha_clad*delta_t_clad)]
duct_outside_flat2flat = ['Duct Outside Flat-to-Flat', hexagon_pitch[1]-gap, (hexagon_pitch[1]-gap)*expand_clad]
interassembly_gap = ['Interassembly Gap', gap, hexagon_pitch[2]-duct_outside_flat2flat[2]]
duct_inside_flat2flat = ['Duct Inside Flat-to-Flat', duct_outside_flat2flat[1]-2*duct_wall_thickness[1],\
                         (duct_outside_flat2flat[1]-2*duct_wall_thickness[1])*expand_clad]
pin_diameter = ['Pin Diameter', pin_diameter, pin_diameter*expand_clad]
pin_to_diameter_ratio = ['Pin to Diameter Ratio', pin_to_diameter,pin_to_diameter]
pin_pitch = ['Pin Pitch', pitch, pin_to_diameter_ratio[2]*pin_diameter[2]]
cladding_thickness = ['Cladding Thickness', clad_thickness, clad_thickness*expand_clad]
cladding_inner_diameter = ['Cladding Inner Diameter', pin_diameter[1]-2*cladding_thickness[1],\
                           (pin_diameter[1]-2*cladding_thickness[1])*expand_clad]
smeared_density = ['Smeared Density', sd_roomtemperature, sd_hot]
fuel_slug_diameter = ['Fuel Slug Diameter', cladding_inner_diameter[1]*math.sqrt(smeared_density[1]*0.01),\
                      cladding_inner_diameter[2]*math.sqrt(smeared_density[2]*0.01)]
wirewrap_diameter = ['Wirewrap Diameter', wire_diameter, wire_diameter*expand_clad]
helical_wirewrap_pitch = ['Helical Wirewrap Pitch', wire_pitch, wire_pitch*expand_clad]
num_pins = ['Number of Pins', 217, 217]
roomtemp_wire_tilt = helical_wirewrap_pitch[1]/math.sqrt((math.pi*(pin_diameter[1]+wirewrap_diameter[1]))**2+\
                                                         helical_wirewrap_pitch[1]**2)
hot_wire_tilt = helical_wirewrap_pitch[2]/math.sqrt((math.pi*(pin_diameter[2]+wirewrap_diameter[2]))**2+\
                                                         helical_wirewrap_pitch[2]**2)
wire_tilt = ['Wirewrap Tilt', roomtemp_wire_tilt, hot_wire_tilt]


# Areas
assembly_area = ['Assembly Area', math.sqrt(3)*hexagon_pitch[1]**2*0.5, math.sqrt(3)*hexagon_pitch[2]**2*0.5]
fuel_slug_area = ['Fuel Slug Area', num_pins[1]*math.pi*fuel_slug_diameter[1]**2*0.25,\
                  num_pins[2]*math.pi*fuel_slug_diameter[2]**2*0.25]
bond_sodium_area = ['Bond Sodium Area', num_pins[1]*math.pi*(cladding_inner_diameter[1]**2-fuel_slug_diameter[1]**2)*0.25,\
                   num_pins[2]*math.pi*(cladding_inner_diameter[2]**2-fuel_slug_diameter[2]**2)*0.25]
cladding_area = ['Cladding Area', num_pins[1]*math.pi*(pin_diameter[1]**2-cladding_inner_diameter[1]**2)*0.25,\
                   num_pins[2]*math.pi*(pin_diameter[2]**2-cladding_inner_diameter[2]**2)*0.25]
roomtemp_wire_area = num_pins[1]*math.pi*(wirewrap_diameter[1]*0.5)**2/wire_tilt[1]
hot_wire_area = num_pins[2]*math.pi*(wirewrap_diameter[2]*0.5)**2/wire_tilt[2]
wire_area = ['Wirewrap Area', roomtemp_wire_area, hot_wire_area]
duct_area = ['Duct Area', math.sqrt(3)*(duct_outside_flat2flat[1]**2-duct_inside_flat2flat[1]**2)*0.5,\
             math.sqrt(3)*(duct_outside_flat2flat[2]**2-duct_inside_flat2flat[2]**2)*0.5]
inside_duct_area = ['Inside the Duct Area',\
                    math.sqrt(3)*duct_outside_flat2flat[1]**2*0.5-fuel_slug_area[1]-bond_sodium_area[1]-cladding_area[1]-\
                   wire_area[1]-duct_area[1],\
                   math.sqrt(3)*duct_outside_flat2flat[2]**2*0.5-fuel_slug_area[2]-bond_sodium_area[2]-cladding_area[2]-\
                   wire_area[2]-duct_area[2]]
inter_assembly_gap_area = ['Inter-Assembly Gap Area', math.sqrt(3)*(hexagon_pitch[1]**2-duct_outside_flat2flat[1]**2)*0.5\
                          , math.sqrt(3)*(hexagon_pitch[2]**2-duct_outside_flat2flat[2]**2)*0.5]

# Volume fractions
roomtemp_fuel_slug_vfrac = fuel_slug_area[1]/assembly_area[1]
hot_fuel_slug_vfrac = fuel_slug_area[2]/assembly_area[2]
roomtemp_bond_sodium_vfrac = bond_sodium_area[1]/assembly_area[1]
hot_bond_sodium_vfrac = bond_sodium_area[2]/assembly_area[2]
fuel_slug_vfrac = ['Fuel Slug Volume Fraction', roomtemp_fuel_slug_vfrac, hot_fuel_slug_vfrac]
bond_sodium_vfrac = ['Bond Sodium Volume Fraction', roomtemp_bond_sodium_vfrac, hot_bond_sodium_vfrac]
total_fuel_vfrac = ['Total Fuel Volume Fraction', fuel_slug_vfrac[1]+bond_sodium_vfrac[1],\
                    fuel_slug_vfrac[2]+bond_sodium_vfrac[2]]

roomtemp_cladding_vfrac = cladding_area[1]/assembly_area[1]
hot_cladding_vfrac = cladding_area[2]/assembly_area[2]
roomtemp_wire_vfrac = wire_area[1]/assembly_area[1]
hot_wire_vfrac = wire_area[2]/assembly_area[2]
roomtemp_duct_vfrac = duct_area[1]/assembly_area[1]
hot_duct_vfrac = duct_area[2]/assembly_area[2]
cladding_vfrac = ['Cladding Volume Fraction', roomtemp_cladding_vfrac, hot_cladding_vfrac]
wire_vfrac = ['Wire Volume Fraction', roomtemp_wire_vfrac, hot_wire_vfrac]
duct_vfrac = ['Duct Volume Fraction', roomtemp_duct_vfrac, hot_duct_vfrac]
total_structure_vfrac = ['Total Structure Volume Fraction', cladding_vfrac[1]+wire_vfrac[1]+duct_vfrac[1],\
                        cladding_vfrac[2]+wire_vfrac[2]+duct_vfrac[2]]

roomtemp_inside_duct_vfrac = inside_duct_area[1]/assembly_area[1]
hot_inside_duct_vfrac = inside_duct_area[2]/assembly_area[2]
roomtemp_inter_assembly_gap_vfrac = inter_assembly_gap_area[1]/assembly_area[1]
hot_inter_assembly_gap_vfrac = inter_assembly_gap_area[2]/assembly_area[2]
inside_duct_vfrac = ['Inside Duct Volume Fraction', roomtemp_inside_duct_vfrac, hot_inside_duct_vfrac]
inter_assembly_gap_vfrac = ['Inter-Assembly Gap  VolumeFraction', roomtemp_inter_assembly_gap_vfrac,\
                            hot_inter_assembly_gap_vfrac]
total_coolant_vfrac = ['Total Coolant Volume Fraction', inside_duct_vfrac[1]+inter_assembly_gap_vfrac[1],\
                      inside_duct_vfrac[2]+inter_assembly_gap_vfrac[2]]

###########################################################################
####                    Homogenized Nuclide Densities                  ####
###########################################################################

# Make vector of all fuel nuclide densities
fuel_nuclide_densities = [u235_n, u238_n, zr90_n, zr91_n, zr92_n, zr94_n, zr96_n]
fuel_homogenized_nuclide_densities = [[[] for i in range(0,4)] for i in range(len(fuel_nuclide_densities))]
# Populate homoginized nuclide densities vector
for i in range(len(fuel_nuclide_densities)):
    fuel_homogenized_nuclide_densities[i][0] = fuel_nuclide_densities[i][0]
    fuel_homogenized_nuclide_densities[i][1] = fuel_nuclide_densities[i][1]*fuel_slug_vfrac[1]
    fuel_homogenized_nuclide_densities[i][2] = fuel_homogenized_nuclide_densities[i][1]*assembly_area[1]/\
    (assembly_area[2]*(1+alpha_fuel*delta_t_fuel))
    fuel_homogenized_nuclide_densities[i][3] = fuel_homogenized_nuclide_densities[i][2]/expand_fuel
    
# Make vector of all cladding nuclide densities
clad_nuclide_densities = [fe54_n,  fe56_n, fe57_n, fe58_n, ni58_n, ni60_n, ni61_n, ni62_n, ni64_n, cr50_n, cr52_n,\
                          cr53_n, cr54_n, mo92_n, mo94_n, mo95_n, mo96_n, mo97_n, mo98_n, mo100_n, mn55_n, c_n,  si28_n,\
                          si29_n, si30_n, w180_n, w182_n, w183_n, w184_n, w186_n, v50_n, v51_n]
clad_homogenized_nuclide_densities = [[[] for i in range(0,4)] for i in range(len(clad_nuclide_densities))]
# Populate homoginized nuclide densities vector
for i in range(len(clad_nuclide_densities)):
    clad_homogenized_nuclide_densities[i][0] = clad_nuclide_densities[i][0]
    clad_homogenized_nuclide_densities[i][1] = clad_nuclide_densities[i][1]*total_structure_vfrac[1]
    clad_homogenized_nuclide_densities[i][2] = clad_homogenized_nuclide_densities[i][1]*assembly_area[1]/\
    (assembly_area[2]*(1+alpha_clad*delta_t_clad))
    clad_homogenized_nuclide_densities[i][3] = clad_homogenized_nuclide_densities[i][2]/expand_clad

# Make vector for coolant nuclide densities
coolant_homogenized_nuclide_densities = [[] for i in range(0,4)]
coolant_homogenized_nuclide_densities[0] = na23_n[0]
coolant_homogenized_nuclide_densities[2] = na23_n[2]*total_coolant_vfrac[2]

###########################################################
####                    Printing Data                  ####
###########################################################

# Make one big vector
homogenized_nuclide_densities = [[] for i in range(len(fuel_homogenized_nuclide_densities)+\
                                                   len(clad_homogenized_nuclide_densities)+\
                                       len(coolant_homogenized_nuclide_densities)-3)]

for i in range(len(fuel_nuclide_densities)):
    homogenized_nuclide_densities[i] = fuel_homogenized_nuclide_densities[i]
for i in range(len(fuel_nuclide_densities),len(fuel_nuclide_densities)+len(clad_nuclide_densities)):
    homogenized_nuclide_densities[i] = clad_homogenized_nuclide_densities[i-len(fuel_nuclide_densities)]
homogenized_nuclide_densities[-1] = coolant_homogenized_nuclide_densities


# Print results
for i in range(len(homogenized_nuclide_densities)):
    print(homogenized_nuclide_densities[i])


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
out_file.write('Enrichment Weight Percent = '+str(enr)+' \n')
out_file.write('Enrichment Weight Percent Split (Outer/Inner) = '+str(1)+' \n')
out_file.write('Smeared Density Percent = '+str(sd_roomtemperature)+' \n')
out_file.write('U-10Zr alloy density = '+str(round(p_fuel,5))+' g/cm^3 \n')
out_file.write('MCC3ID     Density [#/barn-cm]\n')
for i in range(28):
  out_file.write(homogenized_nuclide_densities[i][0]+\
                 str('{:1.5E}'.format(homogenized_nuclide_densities[i][3]*10.**-24.))+\
                 '\n')
out_file.write(homogenized_nuclide_densities[-1][0]+\
               str('{:1.5E}'.format(homogenized_nuclide_densities[-1][2]*10.**-24.))+\
               '\n')

out_file.close()

print('Nuclide Densities Printed To '+out_file_name)