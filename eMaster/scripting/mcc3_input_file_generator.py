#!/usr/bin/python
'''
    This script takes in a run_densities.txt file and a unperturbed MCC3 input
    file. It outputs a MCC3 input file with all number densities replaced 
    specified in the run_densities.txt file replacing whatever previous values
    the original MCC3 input had.
    Additionally, the output file will have all fission products at infinite 
    dilution unless a FP_densities.txt file is passed. Given a FP_densities.txt
    file fission products will similiarly be replaced.
'''
# Helper functions and maps ####################################################

# For formating
indent = '                      '

# List of nuclides in HT9 Cladding - Densities are never changed
clad_list = ['FE54_7', 'FE56_7', 'FE57_7', 'FE58_7', 'CR50_7', 'CR52_7', 
             'CR53_7', 'CR54_7', 'MN55_7', 'NI58_7', 'NI60_7', 'NI61_7',
             'NI62_7', 'NI64_7', 'MO92_7', 'MO94_7', 'MO95_7', 'MO96_7',
             'MO97_7', 'MO98_7', 'MO1007']

coolant_dens;

# Read in densities from run_densities
def read_dens_to_dict(densities):
  line = densities.readline()
  while('MCC3ID     Density [#/barn-cm]' not in line):
    line = densities.readline()
  line = densities.readline()
  dens_dict = {}
  while line:
    split_line = line.split()
    if split_line[0] == 'NA23_7':
      coolant_dens = split_line[1]
    else:
      dens_dict[split_line[0]] = split_line[1]
    line = densities.readline()
  return dens_dict

# Helpers to check id's
def check_is_fuel(dens_dict, id):
  return (id in dens_dict)

def check_is_coolant(id):
  return (id == 'NA23_7')

def check_is_clad(id):
  for i in range(len(clad_list)):
    if id == clad_list[i]:
      return True
  return False

# Generate New MCC3 Input ######################################################

import sys
from math import *
import os.path
from os import path

# Get original mcc3 input file
original_inp = sys.argv[1] 
if (not path.exists(original_inp)):
  original_inp = str(input('MCC3 input not found! Enter MCC3 input file: \n'))

# Get num_densities file
dens_file_name = sys.argv[2]
if (not path.exists(dens_file_name)):
  dens_file_name = str(input('Input not found! Enter run_densities file: \n'))

# Get fission product densities file
try:
  fp_file_name = input('Enter FP_densities file: \n')
except SyntaxError:
  fp_file_name = None
if fp_file_name is not None and (not path.exists(str(fp_file_name))):
  fp_file_name = str(input('Input not found! Enter FP_densities file: \n'))
else:
  print('No fission product densities file passed, I won\'t mess with them.\n')

# Open files
original_handle=open(original_inp,'r')
dens_handle = open(dens_file_name, 'r')
# Temporary ##########################################################################
if fp_file_name:
  print('Only infinite dilution functionality currently implemented')

# Read in run densities from input
dens_dict = read_dens_to_dict(dens_handle)
dens_handle.close()

# Open output file
print_file_name = 'mcc3_modified_input.n'
i = 1
while path.exists(print_file_name):
  print_file_name = 'mcc3_modified_input_'+str(i)+'.n'
  i = i+1
print_file = open(print_file_name,'w')

line = original_handle.readline()
while line:
  fission_products_zone = False
  # Read through to a block that needs changing
  while 'START' not in line and line:
    print_file.write(line)
    line = original_handle.readline()
  print_file.write(line)
  line = original_handle.readline() # Junk START line

  # Read through block while making edits
  while ('END' not in line) and line:
    if '!' not in line: # Skip comments
      split_line = line.split()
      if ( (not fission_products_zone) and check_is_clad(split_line[0]) ):
        print_file.write(indent+ split_line[0]+'  '+split_line[1]+'  '+ \
                         split_line[2]+'  '+split_line[3]+'\n')
      elif ( (not fission_products_zone) and 
             check_is_fuel(dens_dict, split_line[0]) ):
        print_file.write(indent+ split_line[0]+'  '+split_line[1]+'  '+\
                         dens_dict[split_line[0]]+'  '+split_line[3]+'\n')
      elif ( (not fission_products_zone) and check_is_coolant(split_line[0]) ): 
        print_file.write(indent+split_line[0]+'  '+split_line[1]+'  '+ \
                         coolant_dens+'  '+split_line[3]+'\n')  
        #print_file.write(indent+split_line[0]+'  '+split_line[1]+'  '+ \
        #                split_line[2]+'  '+split_line[3]+'\n')      
      else: # Must be fission product
        print_file.write(line)
        #print_file.write(indent+split_line[0]+'  '+split_line[1]+ \
        #                 '  1.00000E-20  '+split_line[3]+'\n')
    elif '! fission products' in line:
      fission_products_zone = True
      print_file.write(line)
    else:
      print_file.write(line)
    # End if / else
    line = original_handle.readline()
# End reading 
print_file.close()

print('Modified MCC3 input printed to '+str(print_file_name))




