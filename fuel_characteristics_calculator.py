'''
    This function takes in the weight enrichment for the fuel, and proceeds to find the nuclide densities for
    U-235, U-238, and Zr. The fuel assumed is a metallic U-10Zr, an alloy of uranium and zirconium that is 10 percent
    Zr by weight.
'''
import math
import matplotlib.pyplot as plt
import numpy as np
import locale

#######################################################################################

def nuclide_solver(w_enrichment):
    w_enrichment = w_enrichment/100
    
    # Atomic masses (amu)
    Na = 6.022*pow(10,23) # Avagadro's number
    m_u235 = 235.0439299
    m_u238 = 238.0507883
    m_zr90 = 89.9047044
    m_zr91 = 90.9056458
    m_zr92 = 91.9050408
    m_zr94 = 93.9063152
    m_zr96 = 95.9082734
    
    # Mass fraction of alloy that is U and Zr
    w_u = 0.9;
    w_zr = 0.1;
    
    # Find densities of U-235 and U-238 metal by first finding the number density of uranium metal
    N_u = 19.1*6.022*pow(10,23)/(238.02891)
    p_u235 = m_u235*N_u/Na
    p_u238 = m_u238*N_u/Na
        
    # Find densities of Zr isotopes by first finding the number density of zirconium metal
    N_zr = 6.51*6.022*pow(10,23)/(91.224)
    p_zr90 = m_zr90*N_zr/Na
    p_zr91 = m_zr91*N_zr/Na
    p_zr92 = m_zr92*N_zr/Na
    p_zr94 = m_zr94*N_zr/Na
    p_zr96 = m_zr96*N_zr/Na

    # Find atomic percent of uranium that is U-235
    a_enrichment = (w_enrichment*m_u238)/(((1-w_enrichment)*m_u235)+(w_enrichment*m_u238))
    # Atomic percentages of natural zirconium that is each isotope
    frac_zr90 = 0.5145
    frac_zr91 = 0.1122
    frac_zr92 = 0.1715
    frac_zr94 = 0.1738
    frac_zr96 = 0.0280

    # Find atomic percent of alloy that is U and Zr
    x_u = (w_u/(a_enrichment*m_u235+(1-a_enrichment)*m_u238))/\
    (w_zr/(frac_zr90*m_zr90+frac_zr91*m_zr91+frac_zr92*m_zr92+frac_zr94*m_zr94+frac_zr96*m_zr96)\
     +w_u/(a_enrichment*m_u235+(1-a_enrichment)*m_u238))
    
    x_zr = (w_zr/(frac_zr90*m_zr90+frac_zr91*m_zr91+frac_zr92*m_zr92+frac_zr94*m_zr94+frac_zr96*m_zr96))/\
    (w_zr/(frac_zr90*m_zr90+frac_zr91*m_zr91+frac_zr92*m_zr92+frac_zr94*m_zr94+frac_zr96*m_zr96)+\
     w_u/(a_enrichment*m_u235+(1-a_enrichment)*m_u238))

    # Find atomic percent of alloy that is each isotope
    x_u235 = a_enrichment*x_u
    x_u238 = (1-a_enrichment)*x_u
    x_zr90 = frac_zr90*x_zr
    x_zr91 = frac_zr91*x_zr
    x_zr92 = frac_zr92*x_zr
    x_zr94 = frac_zr94*x_zr
    x_zr96 = frac_zr96*x_zr

    # Find mass-density of alloy, molar mass of alloy, and then number density of alloy
    p_alloy = x_zr90*p_zr90+x_zr91*p_zr91+x_zr92*p_zr92+x_zr94*p_zr94+x_zr96*p_zr96+x_u235*p_u235+x_u238*p_u238
    m_alloy = x_zr90*m_zr90+x_zr91*m_zr91+x_zr92*m_zr92+x_zr94*m_zr94+x_zr96*m_zr96+x_u235*m_u235+x_u238*m_u238
    n_alloy = p_alloy*Na/m_alloy

    # Find number densities of fuel nuclides in barn-cm
    n_u235 = n_alloy*x_u235*pow(10,-24)
    n_u238 = n_alloy*x_u238*pow(10,-24)
    n_zr90 = n_alloy*x_zr90*pow(10,-24)
    n_zr91 = n_alloy*x_zr91*pow(10,-24)
    n_zr92 = n_alloy*x_zr92*pow(10,-24)
    n_zr94 = n_alloy*x_zr94*pow(10,-24)
    n_zr96 = n_alloy*x_zr96*pow(10,-24)
    
    return (round(p_alloy,5), round(n_u235,5), round(n_u238,5), round(n_zr90,5), round(n_zr91,5), round(n_zr92,5),\
           round(n_zr94,5), round(n_zr96,5))

#######################################################################################

def mass_solver(p_fuel, sd, num_assemblies):
    # Find number of fuel pins
    num_pins = num_assemblies*271
    # Find total volume of fuel
    clad_inner_diameter = 0.755-0.056 # (cm)
    fuel_diameter = clad_inner_diameter*sd/100
    fuel_height = 81.3
    v_fuel = num_pins*((math.pi/4)*pow(fuel_diameter,2)*fuel_height)
    
    # Find mass of fuel in kg
    fuel_mass = v_fuel*p_fuel/1000 
    
    return fuel_mass

def cost_solver(w_enrichment, fuel_mass):
    feedX = 0.00711 # Feed enrichment
    productX = w_enrichment/100
    tailX = 0.00228 # Tailing enrichment
    
    # Prices in dollars per kilo
    SWUP = 115.42
    convP = 12
    mineP = 85.56
    
    # Find masses at different points in the process
    feed_mass = fuel_mass*(productX-tailX)/(feedX-tailX)
    tail_mass = feed_mass-fuel_mass
    conv_mass = 1.05*feed_mass
    mine_mass = conv_mass*1.1792499
    
    # Calculate SWU and costs
    SWU = fuel_mass*(value_function(productX)-value_function(tailX))-\
    feed_mass*(value_function(feedX)-value_function(tailX))
    
    SWU_cost = SWU*SWUP
    conv_cost = conv_mass*convP
    mine_cost = mine_mass*mineP

    total_cost = round(SWU_cost+conv_cost+mine_cost,2)
    
    return total_cost

#######################################################################################

# This is the value function used in SWU calculations
def value_function(x):
    return (1-2*x)*math.log((1-x)/x)

#######################################################################################
# Run all functions in conjunction to produce a full core analysis

def solve_core_compositionANDcosts():
    locale.setlocale(locale.LC_ALL, 'en_US')
    # Poll for inputs
    enrichment_inner = float(input('Enter fuel\'s % weight enrichment for the inner core: '))
    enrichment_outer = float(input('Enter fuel\'s % weight enrichment for the outer core: '))
    sd_inner = float(input('Enter fuel\'s smeared density in percent: '))
    sd_outer = sd_inner
    sodium_fraction = float(input('Enter percent sodium area fraction in pin cell: '))/100
    fuel_fraction = float(input('Enter percent fuel area fraction in pin cell: '))/100
    
    # Find densities for inner and outer core
    (p_fuel_inner, n_u235_inner, n_u238_inner, n_zr90_inner, n_zr91_inner, n_zr92_inner, n_zr94_inner, n_zr96_inner)\
    = nuclide_solver(enrichment_inner)
    (p_fuel_outer, n_u235_outer, n_u238_outer, n_zr90_outer, n_zr91_outer, n_zr92_outer, n_zr94_outer, n_zr96_outer)\
    = nuclide_solver(enrichment_inner)
    
    # Find fuel masses for inner and outer core
    inner_mass = mass_solver(p_fuel_inner, sd_inner, 78)
    outer_mass = mass_solver(p_fuel_outer, sd_outer, 102)
    
    # Find costs for inner and outer core
    inner_cost = cost_solver(enrichment_inner, inner_mass)
    outer_cost = cost_solver(enrichment_outer, outer_mass)
    
    # Find adjusted number densities for MC^2
    sodium_density = 0.927
    sodium_mass = 22.989769
    n_sodium = sodium_density*6.022*pow(10,23)/sodium_mass
    n_sodium_adjusted = n_sodium*sodium_fraction*pow(10,-24)
    n_u235_inner_adjusted = round(n_u235_inner*fuel_fraction,5)
    n_u238_inner_adjusted = round(n_u238_inner*fuel_fraction,5)
    n_zr90_inner_adjusted = round(n_zr90_inner*fuel_fraction,5)
    n_zr91_inner_adjusted = round(n_zr91_inner*fuel_fraction,5)
    n_zr92_inner_adjusted = round(n_zr92_inner*fuel_fraction,5)
    n_zr94_inner_adjusted = round(n_zr94_inner*fuel_fraction,5)
    n_zr96_inner_adjusted = round(n_zr96_inner*fuel_fraction,5)
    n_u235_outer_adjusted = round(n_u235_outer*fuel_fraction,5)
    n_u238_outer_adjusted = round(n_u238_outer*fuel_fraction,5)
    n_zr90_outer_adjusted = round(n_zr90_outer*fuel_fraction,5)
    n_zr91_outer_adjusted = round(n_zr91_outer*fuel_fraction,5)
    n_zr92_outer_adjusted = round(n_zr92_outer*fuel_fraction,5)
    n_zr94_outer_adjusted = round(n_zr94_outer*fuel_fraction,5)
    n_zr96_outer_adjusted = round(n_zr96_outer*fuel_fraction,5)
    
    print('\n******************************************************************')
    print('INNER CORE DATA\n')
    print('Density of U-10Zr fuel (inner core):\t',p_fuel_inner,'\t(g/cm3)\n')
    print('Number density of U-235 (inner core):\t',n_u235_inner,'\t(1/b-cm)')
    print('Number density of U-238 (inner core):\t',n_u238_inner,'\t(1/b-cm)')
    print('Number density of Zr-90 (inner core):\t',n_zr90_inner,'\t(1/b-cm)')
    print('Number density of Zr-91 (inner core):\t',n_zr91_inner,'\t(1/b-cm)')
    print('Number density of Zr-92 (inner core):\t',n_zr92_inner,'\t(1/b-cm)')
    print('Number density of Zr-94 (inner core):\t',n_zr94_inner,'\t(1/b-cm)')
    print('Number density of Zr-96 (inner core):\t',n_zr96_inner,'\t(1/b-cm)\n')
    print('Inner core mass:\t',round(inner_mass,3),'\t(kg)\n')
    print('Inner core cost:\t$',locale.format_string("%d", inner_cost, grouping=True),'\n')
    print('******************************************************************')
    print('OUTER CORE DATA\n')
    print('Density of U-10Zr fuel (outer core):\t',p_fuel_outer,'\t(g/cm3)\n')
    print('Number density of U-235 (outer core):\t',n_u235_outer,'\t(1/b-cm)')
    print('Number density of U-238 (outer core):\t',n_u238_outer,'\t(1/b-cm)')
    print('Number density of Zr-90 (outer core):\t',n_zr90_outer,'\t(1/b-cm)')
    print('Number density of Zr-91 (outer core):\t',n_zr91_outer,'\t(1/b-cm)')
    print('Number density of Zr-92 (outer core):\t',n_zr92_outer,'\t(1/b-cm)')
    print('Number density of Zr-94 (outer core):\t',n_zr94_outer,'\t(1/b-cm)')
    print('Number density of Zr-96 (outer core):\t',n_zr96_outer,'\t(1/b-cm)\n')
    print('Outer core mass:\t',round(outer_mass,3),'\t(kg)\n')
    print('Outer core cost:\t$',locale.format_string("%d", outer_cost, grouping=True),'\n')
    print('******************************************************************')
    print('TOTAL CORE DATA\n')
    print('Total core mass:\t',round(inner_mass+outer_mass,3),'\t(kg)\n')
    print('Total core cost:\t$',locale.format_string("%d", inner_cost+outer_cost, grouping=True),'\n')    
    print('******************************************************************')
    print('ADJUSTED NUMBER DENSITIES FOR MC2\n')
    print('Number density of sodium:\t',round(n_sodium_adjusted,5),'\t(1/b-cm)\n')
    
    print('Number density of U-235 (inner core):\t',n_u235_inner_adjusted,'\t(1/b-cm)')
    print('Number density of U-238 (inner core):\t',n_u238_inner_adjusted,'\t(1/b-cm)')
    print('Number density of Zr-90 (inner core):\t',n_zr90_inner_adjusted,'\t(1/b-cm)')
    print('Number density of Zr-91 (inner core):\t',n_zr91_inner_adjusted,'\t(1/b-cm)')
    print('Number density of Zr-92 (inner core):\t',n_zr92_inner_adjusted,'\t(1/b-cm)')
    print('Number density of Zr-94 (inner core):\t',n_zr94_inner_adjusted,'\t(1/b-cm)')
    print('Number density of Zr-96 (inner core):\t',n_zr96_inner_adjusted,'\t(1/b-cm)\n')

    print('Number density of U-235 (outer core):\t',n_u235_outer_adjusted,'\t(1/b-cm)')
    print('Number density of U-238 (outer core):\t',n_u238_outer_adjusted,'\t(1/b-cm)')
    print('Number density of Zr-90 (outer core):\t',n_zr90_outer_adjusted,'\t(1/b-cm)')
    print('Number density of Zr-91 (outer core):\t',n_zr91_outer_adjusted,'\t(1/b-cm)')
    print('Number density of Zr-92 (outer core):\t',n_zr92_outer_adjusted,'\t(1/b-cm)')
    print('Number density of Zr-94 (outer core):\t',n_zr94_outer_adjusted,'\t(1/b-cm)')
    print('Number density of Zr-96 (outer core):\t',n_zr96_outer_adjusted,'\t(1/b-cm)\n')
    
    
    
    
solve_core_compositionANDcosts()