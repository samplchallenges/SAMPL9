#!/usr/bin/env python

# Credit:
# This adapted from Andrea Rizzi's file of the same name which he wrote for SAMPL6
# at https://github.com/samplchallenges/SAMPL6/blob/master/host_guest/Analysis/ExperimentalMeasurements/generate_tables.py

# =============================================================================
# GLOBAL IMPORTS
# =============================================================================

import os
import math
import csv
import json
from collections import OrderedDict

import numpy as np
from simtk import unit as u


# =============================================================================
# CONSTANTS
# =============================================================================

T = 298 * u.kelvin
R = u.MOLAR_GAS_CONSTANT_R
RELATIVE_TITRANT_CONC_ERROR = 0.03

WP6_GUESTS_SMILES_PATH = '/home/amezcum1/SAMPL9/host_guest/WP6/guest_files/WP6_guest_smiles.txt'
#WP6_GUESTS_NAMES_PATH = '/home/amezcum1/SAMPL9/host_guest/WP6/guest_files/WP6_guest_names.txt'
#CD_GUESTS_NAMES_PATH = '/home/amezcum1/SAMPL9/host_guest/bCD/guest_files/bCD_guest_names.txt'
CD_GUESTS_SMILES_PATH = '/home/amezcum1/SAMPL9/host_guest/bCD/guest_files/guest_smiles.txt'
CD_HOST_NAMES = ['bCD', 'HbCD']

# Experimental results as provided by the Isaacs and Gilson groups.
# The error is relative. None means that the error is <1%.
EXPERIMENTAL_DATA = OrderedDict([

    ('bCD-PMZ', OrderedDict([
        ('DG', -20.81 * u.kilojoules_per_mole), ('dDG', None * u.kilojoules_per_mole),
        ('DH', -24.76 * u.kilojoules_per_mole), ('dDH', None * u.kilojoules_per_mole),
        ('TDS', 3.95 * u.kilojoules_per_mole), ('dTDS', None * u.kilojoule_per_mole),
        ('n', 1.09)
    ])),
    ('bCD-PMT', OrderedDict([
        ('DG', -18.73 * u.kilojoules_per_mole), ('dDG', None * u.kilojoules_per_mole),
        ('DH', -16.45 * u.kilojoules_per_mole), ('dDH', None * u.kilojoules_per_mole),
        ('TDS', -2.28 * u.kilojoules_per_mole), ('dTDS', None * u.kilojoule_per_mole),
        ('n', 0.94)
    ])),
    ('bCD-CPZ', OrderedDict([
        ('DG', -22.66 * u.kilojoules_per_mole), ('dDG', None * u.kilojoules_per_mole),
        ('DH', -26.62 * u.kilojoules_per_mole), ('dDH', None * u.kilojoules_per_mole),
        ('TDS', 3.96 * u.kilojoules_per_mole), ('dTDS', None * u.kilojoule_per_mole),
        ('n', 0.77)
    ])),
    ('bCD-TDZ', OrderedDict([
        ('DG', -23.86 * u.kilojoules_per_mole), ('dDG', None * u.kilojoules_per_mole),
        ('DH', -20.62 * u.kilojoules_per_mole), ('dDH', None * u.kilojoules_per_mole),
        ('TDS', -3.24 * u.kilojoules_per_mole), ('dTDS', None * u.kilojoule_per_mole),
        ('n', 1.14)
    ])),
    ('bCD-TFP', OrderedDict([
        ('DG', -21.18 * u.kilojoules_per_mole), ('dDG', None * u.kilojoules_per_mole),
        ('DH', -16.31 * u.kilojoules_per_mole), ('dDH', None * u.kilojoules_per_mole),
        ('TDS', -4.87 * u.kilojoules_per_mole), ('dTDS', None * u.kilojoule_per_mole),
        ('n', 1.18)
    ])),
    ('HbCD-PMZ', OrderedDict([
        ('DG', -21.15 * u.kilojoules_per_mole), ('dDG', None * u.kilojoules_per_mole),
        ('DH', -21.41 * u.kilojoules_per_mole), ('dDH', None * u.kilojoules_per_mole),
        ('TDS', 0.26 * u.kilojoules_per_mole), ('dTDS', None * u.kilojoule_per_mole),
        ('n', 0.99)
    ])),
    ('HbCD-PMT', OrderedDict([
        ('DG', -22.42 * u.kilojoules_per_mole), ('dDG', None * u.kilojoules_per_mole),
        ('DH', -16.96 * u.kilojoules_per_mole), ('dDH', None * u.kilojoules_per_mole),
        ('TDS', -5.46 * u.kilojoules_per_mole), ('dTDS', None * u.kilojoule_per_mole),
        ('n', 1.06)
    ])),
    ('HbCD-CPZ', OrderedDict([
        ('DG', -22.60 * u.kilojoules_per_mole), ('dDG', None * u.kilojoules_per_mole),
        ('DH', -24.95 * u.kilojoules_per_mole), ('dDH', None * u.kilojoules_per_mole),
        ('TDS', 2.35 * u.kilojoules_per_mole), ('dTDS', None * u.kilojoule_per_mole),
        ('n', 0.77)
    ])),
    ('HbCD-TDZ', OrderedDict([
        ('DG', -27.03 * u.kilojoules_per_mole), ('dDG', None * u.kilojoules_per_mole),
        ('DH', -38.34 * u.kilojoules_per_mole), ('dDH', None * u.kilojoules_per_mole),
        ('TDS', 11.31 * u.kilojoules_per_mole), ('dTDS', None * u.kilojoule_per_mole),
        ('n', 0.87)
    ])),
    ('HbCD-TFP', OrderedDict([
        ('DG', -23.17 * u.kilojoules_per_mole), ('dDG', None * u.kilojoules_per_mole),
        ('DH', -31.05 * u.kilojoules_per_mole), ('dDH', None * u.kilojoules_per_mole),
        ('TDS', 7.88 * u.kilojoules_per_mole), ('dTDS', None * u.kilojoule_per_mole),
        ('n', 0.56)
    ])),
    ('WP6-G1', OrderedDict([
        ('Ka_1', 52900 / u.molar), ('dKa_1', 700 / u.molar),
        ('Ka_2', 51500 / u.molar), ('dKa_2', 1600 / u.molar),
        ('DH_1', -8.08 * u.kilocalories_per_mole), ('dDH_1', 0.02 * u.kilocalories_per_mole),
        ('DH_2', -7.87 * u.kilocalories_per_mole), ('dDH_2', 0.03 * u.kilocalories_per_mole),
        #('TDS_1', 8.09 * u.kilojoules_per_mole), ('dTDS_1', 0.45 * u.kilojoules_per_mole),
        #('TDS_2', 8.56 * u.kilojoules_per_mole), ('dTDS_2', 0.47 * u.kilojoules_per_mole),
        ('TDS', None), ('dTDS', None),
        ('n', 1.00)
    ])),
    ('WP6-G2', OrderedDict([
        ('Ka_1', 45900000 / u.molar), ('dKa_1', 3500000 / u.molar),
        ('Ka_2', 44800000 / u.molar), ('dKa_2', 2100000 / u.molar),
        ('DH_1', -6.10 * u.kilocalories_per_mole), ('dDH_1', 0.02 * u.kilocalories_per_mole),
        ('DH_2', -6.15 * u.kilocalories_per_mole), ('dDH_2', 0.02 * u.kilocalories_per_mole),
        #('TDS_1', -17.04 * u.kilojoules_per_mole), ('dTDS_1', 1.75 * u.kilojoules_per_mole),
        #('TDS_2', -17.84 * u.kilojoules_per_mole), ('dTDS_2', 1.76 * u.kilojoules_per_mole),
        ('TDS', None), ('dTDS', None),
        ('n', 1.00)
    ])),
    ('WP6-G3', OrderedDict([
        ('Ka_1', 645000 / u.molar), ('dKa_1', 18000 / u.molar),
        ('Ka_2', 625000 / u.molar), ('dKa_2', 30000 / u.molar),
        ('DH_1', -4.75 * u.kilocalories_per_mole), ('dDH_1', 0.02 * u.kilocalories_per_mole),
        ('DH_2', -4.68 * u.kilocalories_per_mole), ('dDH_2', 0.02 * u.kilocalories_per_mole),
        #('TDS_1', 4.54 * u.kilojoules_per_mole), ('dTDS_1', 1.40 * u.kilojoules_per_mole),
        #('TDS_2', 8.23 * u.kilojoules_per_mole), ('dTDS_2', 0.53 * u.kilojoules_per_mole),
        ('TDS', None), ('dTDS', None),
        ('n', 1.00)
    ])),
    ('WP6-G4', OrderedDict([
        ('Ka_1', 50800 / u.molar), ('dKa_1', 1100 / u.molar),
        ('Ka_2', 49700 / u.molar), ('dKa_2', 1600 / u.molar),
        ('DH_1', -4.15 * u.kilocalories_per_mole), ('dDH_1', 0.02 * u.kilocalories_per_mole),
        ('DH_2', -3.98 * u.kilocalories_per_mole), ('dDH_2', 0.02 * u.kilocalories_per_mole),
        #('TDS_1', -10.78 * u.kilojoules_per_mole), ('dTDS_1', 1.19 * u.kilojoules_per_mole),
        #('TDS_2', -13.56 * u.kilojoules_per_mole), ('dTDS_2', 1.47 * u.kilojoules_per_mole),
        ('TDS', None), ('dTDS', None),
        ('n', 1.00)
    ])),
    ('WP6-G5', OrderedDict([
        ('Ka_1', 9010 / u.molar), ('dKa_1', 230 / u.molar),
        ('Ka_2', 9310 / u.molar), ('dKa_2', 220 / u.molar),
        ('DH_1', -3.95 * u.kilocalories_per_mole), ('dDH_1', 0.03 * u.kilocalories_per_mole),
        ('DH_2', -4.15 * u.kilocalories_per_mole), ('dDH_2', 0.03 * u.kilocalories_per_mole),
        #('TDS_1', 2.98 * u.kilojoules_per_mole), ('dTDS_1', 1.80 * u.kilojoules_per_mole),
        #('TDS_2', 0.81 * u.kilojoules_per_mole), ('dTDS_2', 2.66 * u.kilojoules_per_mole),
        ('TDS', None), ('dTDS', None),
        ('n', 1.00)
    ])),
    ('WP6-G6', OrderedDict([
        ('Ka_1', 709000 / u.molar), ('dKa_1', 44000 / u.molar),
        ('Ka_2', 726000 / u.molar), ('dKa_2', 44000 / u.molar),
        ('DH_1', -6.90 * u.kilocalories_per_mole), ('dDH_1', 0.07 * u.kilocalories_per_mole),
        ('DH_2', -7.10 * u.kilocalories_per_mole), ('dDH_2', 0.04 * u.kilocalories_per_mole),
        #('TDS_1', -21.82 * u.kilojoules_per_mole), ('dTDS_1', 1.73 * u.kilojoules_per_mole),
        #('TDS_2', -21.21 * u.kilojoules_per_mole), ('dTDS_2', 5.73 * u.kilojoules_per_mole),
        ('TDS', None), ('dTDS', None),
        ('n', 1.00)
    ])),
    ('WP6-G7', OrderedDict([
        ('Ka_1', 131000 / u.molar), ('dKa_1', 5000 / u.molar),
        ('Ka_2', 129000 / u.molar), ('dKa_2', 5000 / u.molar),
        ('DH_1', -3.18 * u.kilocalories_per_mole), ('dDH_1', 0.02 * u.kilocalories_per_mole),
        ('DH_2', -3.17 * u.kilocalories_per_mole), ('dDH_2', 0.02 * u.kilocalories_per_mole),
        #('TDS_1', 4.05 * u.kilojoules_per_mole), ('dTDS_1', 0.73 * u.kilojoules_per_mole),
        #('TDS_2', 5.51 * u.kilojoules_per_mole), ('dTDS_2', 0.53 * u.kilojoules_per_mole),
        ('TDS', None), ('dTDS', None),
        ('n', 1.00)
    ])),
    ('WP6-G8', OrderedDict([
        ('Ka_1', 23500 / u.molar), ('dKa_1', 400 / u.molar),
        ('Ka_2', 23000 / u.molar), ('dKa_2', 600 / u.molar),
        ('DH_1', -9.55 * u.kilocalories_per_mole), ('dDH_1', 0.05 * u.kilocalories_per_mole),
        ('DH_2', -9.51 * u.kilocalories_per_mole), ('dDH_2', 0.05 * u.kilocalories_per_mole),
        #('TDS_1', -33.19 * u.kilojoules_per_mole), ('dTDS_1', 3.40 * u.kilojoules_per_mole),
        #('TDS_2', -34.21 * u.kilojoules_per_mole), ('dTDS_2', 2.75 * u.kilojoules_per_mole),
        ('TDS', None), ('dTDS', None),
        ('n', 1.00)
    ])),
    ('WP6-G9', OrderedDict([
            ('Ka_1', 37500 / u.molar), ('dKa_1', 3100 / u.molar),
            ('Ka_2', 36800 / u.molar), ('dKa_2', 3700 / u.molar),
            ('DH_1', -5.31 * u.kilocalories_per_mole), ('dDH_1', 0.08 * u.kilocalories_per_mole),
            ('DH_2', -5.21 * u.kilocalories_per_mole), ('dDH_2', 0.08 * u.kilocalories_per_mole),
            #('TDS_1', 3.76 * u.kilojoules_per_mole), ('dTDS_1', 1.68 * u.kilojoules_per_mole),
            #('TDS_2', 0.33 * u.kilojoules_per_mole), ('dTDS_2', 4.18 * u.kilojoules_per_mole),
            ('TDS', None), ('dTDS', None),
            ('n', 1.00)
        ])),
    ('WP6-G10', OrderedDict([
        ('Ka_1', 16100000 / u.molar), ('dKa_1', 800000 / u.molar),
        ('Ka_2', 16000000 / u.molar), ('dKa_2', 1200000 / u.molar),
        ('DH_1', -6.23 * u.kilocalories_per_mole), ('dDH_1', 0.02 * u.kilocalories_per_mole),
        ('DH_2', -6.11 * u.kilocalories_per_mole), ('dDH_2', 0.02 * u.kilocalories_per_mole),
        #('TDS_1', -13.96 * u.kilojoules_per_mole), ('dTDS_1', 0.48 * u.kilojoules_per_mole),
        #('TDS_2', -13.14 * u.kilojoules_per_mole), ('dTDS_2', 1.33 * u.kilojoules_per_mole),
        ('TDS', None), ('dTDS', None),
        ('n', 1.00)
    ])),
    ('WP6-G11', OrderedDict([
        ('Ka_1', 33700 / u.molar), ('dKa_1', 500 / u.molar),
        ('Ka_2', 33100 / u.molar), ('dKa_2', 500 / u.molar),
        ('DH_1', -5.61 * u.kilocalories_per_mole), ('dDH_1', 0.02 * u.kilocalories_per_mole),
        ('DH_2', -5.45 * u.kilocalories_per_mole), ('dDH_2', 0.03 * u.kilocalories_per_mole),
        #('TDS_1', 6.28 * u.kilojoules_per_mole), ('dTDS_1', 0.84 * u.kilojoules_per_mole),
        #('TDS_2', 8.23 * u.kilojoules_per_mole), ('dTDS_2', 0.54 * u.kilojoules_per_mole),
        ('TDS', None), ('dTDS', None),
        ('n', 1.00)
    ])),
    ('WP6-G12', OrderedDict([
        ('Ka_1', 94300000 / u.molar), ('dKa_1', 3100000 / u.molar),
        ('Ka_2', 84000000 / u.molar), ('dKa_2', 5100000 / u.molar),
        ('DH_1', -7.45 * u.kilocalories_per_mole), ('dDH_1', 0.02 * u.kilocalories_per_mole),
        ('DH_2', -7.41 * u.kilocalories_per_mole), ('dDH_2', 0.02 * u.kilocalories_per_mole),
        #('TDS_1', -17.85 * u.kilojoules_per_mole), ('dTDS_1', 1.53 * u.kilojoules_per_mole),
        #('TDS_2', -19.59 * u.kilojoules_per_mole), ('dTDS_2', 1.56 * u.kilojoules_per_mole),
        ('TDS', None), ('dTDS', None),
        ('n', 1.00)
    ])),
    ('WP6-G13', OrderedDict([
        ('Ka_1', 1630000 / u.molar), ('dKa_1', 110000 / u.molar),
        ('Ka_2', 1610000 / u.molar), ('dKa_2', 110000 / u.molar),
        ('DH_1', -4.98 * u.kilocalories_per_mole), ('dDH_1', 0.04 * u.kilocalories_per_mole),
        ('DH_2', -4.96 * u.kilocalories_per_mole), ('dDH_2', 0.03 * u.kilocalories_per_mole),
        #('TDS_1', 0.70 * u.kilojoules_per_mole), ('dTDS_1', 0.74 * u.kilojoules_per_mole),
        #('TDS_2', 0.63 * u.kilojoules_per_mole), ('dTDS_2', 0.72 * u.kilojoules_per_mole),
        ('TDS', None), ('dTDS', None),
        ('n', 1.00)
    ])),
])



# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def load_smiles(file_path):
    """Return the list of guests IDs and SMILES."""
    guests = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            smiles, gid = line.split(';', 1)
            guests.append([smiles.strip(), gid.strip()])
    return guests

def load_names(file_path):
    """Return the list of guests IDs and names."""
    guests = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            name, gid = line.split(';', 1)
            guests.append([name.strip(), gid.strip()])
    return guests


def compute_DG(Ka, dKa):
    """Compute the free energy from the association constant.

    Parameters
    ----------
    Ka : simtk.Quantity
        Association constant.
    dKa : simtk.Quantity
        Association constant uncertainty.

    Returns
    -------
    DG : simtk.Quantity
        Binding free energy.
    dDG : simtk.Quantity
        Binding free energy uncertainty.

    """
    concentration_unit = 1 / Ka.unit
    DG = -R * T * np.log(Ka*concentration_unit)
    # Propagate error.
    if dKa is None:
        dDG = None
    else:
        dDGdKa = -R * T / Ka  # Derivative dDG(Ka)/dKa.
        # Have to use u.sqrt to avoid bug with simtk.unit
        dDG = u.sqrt(dDGdKa**2 * dKa**2)
    return DG, dDG

def compute_Ka(DG, dDG):
    """Compute the association constant from the free energy.

    Parameters
    ----------
    DG : simtk.Quantity
        Free energy
    dDG : simtk.Quantity
        Uncertainty in free energy

    Returns
    -------
    Ka : simtk.Quantity
        Association constant.
    dKa : simtk.Quantity
        Association constant uncertainty.

    """
    concentration_unit = u.molar
    Ka = np.exp(-DG/(R*T))*1/concentration_unit
    # Propagate error.
    if dDG is None:
        dKa = None
    else:
        dKadDG = - Ka / (R*T)  # Derivative dKa(DG)/dDG.
        dKa = u.sqrt(dKadDG**2 * dDG**2)

    return Ka, dKa


def compute_TDS(DG, dDG, DH, dDH):
    """Compute the entropy from free energy and enthalpy.

    Parameters
    ----------
    DG : simtk.Quantity
        Free energy.
    dDG : simtk.Quantity
        Free energy uncertainty.
    DH : simtk.Quantity
        Enthalpy.
    dDH : simtk.Quantity
        Enthalpy uncertainty.

    Returns
    -------
    TDS : simtk.Quantity
        Entrop.
    dTDS : simtk.Quantity
        Binding free energy uncertainty.

    """
    TDS = DH - DG
    dTDS = u.sqrt(dDH**2 + dDG**2)
    return TDS, dTDS


def strip_units(quantities):
    for k, v in quantities.items():
        if isinstance(v, u.Quantity):
            # We only have energies and association and dissociation constants.
            if 'Ka' in k:
                quantities[k] = v.value_in_unit(v.unit)
            elif 'Kd' in k:
                quantities[k] = v.value_in_unit(v.unit)
            else:
                quantities[k] = v.value_in_unit(u.kilocalories_per_mole)


def reduce_to_first_significant_digit(quantity, uncertainty):
    # If strings, just return
    if isinstance(quantity, str) or isinstance(uncertainty, str):
        return quantity, uncertainty
    first_significant_digit = math.floor(math.log10(abs(uncertainty)))
    quantity = round(quantity, -first_significant_digit)
    uncertainty = round(uncertainty, -first_significant_digit)
    return quantity, uncertainty


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    # Load names and SMILES of guests.
    molecule_names = {}

    smiles_by_host = {
        'WP6': load_smiles(WP6_GUESTS_SMILES_PATH),
        'bCD' : load_smiles(CD_GUESTS_SMILES_PATH),
        'HbCD' : load_smiles(CD_GUESTS_SMILES_PATH),
    }

    #smiles_by_host = {
    #        'WP6': load_smiles(WP6_GUESTS_SMILES_PATH)
    #            }

    names_by_host = {
        'WP6': load_names(WP6_GUESTS_NAMES_PATH),
        'bCD' : load_names(CD_GUESTS_NAMES_PATH),
        'HbCD' : load_names(CD_GUESTS_NAMES_PATH),
    }

    #names_by_host = {
    #    'WP6': load_names(WP6_GUESTS_NAMES_PATH)
    #    }

    for host in CD_HOST_NAMES:
        smiles_by_host[host] = load_smiles(CD_GUESTS_SMILES_PATH)
        names_by_host[host] = load_names(CD_GUESTS_NAMES_PATH)

    #for host in ['WP6']:
    for host in ['WP6']+CD_HOST_NAMES:
        molecule_names[host] = {}
        for smi, gid in smiles_by_host[host]:
            for name, gid2 in names_by_host[host]:
                if gid==gid2:
                    molecule_names[host][gid] = smi, name

    output_dict = OrderedDict()
    upper_bound_molecules = dict(Ka=set(), DH=set(), TDS=set())

    for system_name, system_data in EXPERIMENTAL_DATA.items():
        host_name, gid = system_name.split('-')

        # Load SMILES and common name of the molecule.
        molecule_smiles, molecule_name = molecule_names[host_name][gid]

        # Create entry in the output dictionary.
        output_dict[system_name] = OrderedDict([
            ('name', molecule_name),
            ('SMILES', molecule_smiles),
        ])
        output_dict[system_name].update(system_data)
        system_data = output_dict[system_name]  # Shortcut.

        # If this data has two values, combine: First deal with measured values
        # Note that Kd/Ka values should be combined as free energies (since that's the normally distributed quantity)
        for data_type in ['Kd', 'Ka', 'DH']:
            if data_type+'_1' in system_data:
                if 'DH' in data_type: #just take mean
                    final_val = np.mean( [system_data[data_type+'_1'], system_data[data_type+'_2']])
                    system_data[data_type] = final_val

                    # Also compute uncertainty -- the larger of the propagated uncertainty and the standard error in the mean
                    final_unc = u.sqrt( system_data['d'+data_type+'_1']**2 + system_data['d'+data_type+'_2']**2 )
                    std_err = u.sqrt( 0.5*( (system_data[data_type+'_1']-final_val)**2 + (system_data[data_type+'_2']-final_val)**2) )
                    if std_err > final_unc:
                        final_unc = std_err
                    system_data['d'+data_type] = final_unc
                #Otherwise first convert to free energy then take mean
                elif 'Kd' in data_type: #First convert to free energy then take mean
                    # If we have Kd data instead of Ka data, convert
                    if 'Kd_1' in system_data and not 'Ka_1' in system_data:
                        system_data['Ka_1'] = 1./system_data['Kd_1']
                        # Handle uncertainty -- 1/Kd^2 * dKd
                        system_data['dKa_1'] = system_data['dKd_1']/(system_data['Kd_1']**2)
                        system_data['Ka_2'] = 1./system_data['Kd_2']
                        # Handle uncertainty -- 1/Kd^2 * dKd
                        system_data['dKa_2'] = system_data['dKd_2']/(system_data['Kd_2']**2)
                elif 'Ka' in data_type:
                    if 'Ka_1' in system_data and not 'Ka' in system_data:
                        # Now convert to free energy
                        DG_1, dDG_1 = compute_DG(system_data['Ka_1'], system_data['dKa_1'])
                        DG_2, dDG_2 = compute_DG(system_data['Ka_2'], system_data['dKa_2'])
                        # Take mean
                        DG = (DG_1+DG_2)/2.
                        # Compute uncertainty
                        final_unc = u.sqrt( (dDG_1)**2 + (dDG_2)**2)
                        std_err = u.sqrt( 0.5*( (DG_1-DG)**2 + (DG_2-DG)**2) )
                        if std_err > final_unc:
                            final_unc = std_err
                        # Convert back to Ka and store
                        Ka, dKa = compute_Ka( DG, final_unc)
                        system_data['Ka'] = Ka
                        system_data['dKa'] = dKa


        # Incorporate the relative concentration uncertainties into quantities.
        # Skip this for the Gibb data, where concentration errors are already accounted for and we already have free energy
        TDS = None
        dTDS = None
        DG = None
        dDG = None
        if not 'OA' in system_name:
            for k in ['Ka', 'DH']:
                quantity = system_data[k]
                # Compute relative uncertainty
                relative_uncertainty = system_data['d' + k]/quantity
                # Use upper-bound of 1% if <1% is reported. Keep track of these molecules.
                if relative_uncertainty is None:
                    upper_bound_molecules[k].add(system_name)
                    relative_uncertainty = 0.01
                # Incorporate the relative concentration uncertainties into quantities.
                relative_uncertainty = u.sqrt( relative_uncertainty**2 + RELATIVE_TITRANT_CONC_ERROR**2)
                # Convert relative to absolute errors.
                system_data['d' + k] = abs(quantity * relative_uncertainty)

            # Propagate Ka and DH error into DG and TDS.
            DG, dDG = compute_DG(system_data['Ka'], system_data['dKa'])
            system_data['DG'] = DG
            system_data['dDG'] = dDG
            TDS, dTDS = compute_TDS(system_data['DG'], system_data['dDG'],
                                    system_data['DH'], system_data['dDH'])
            system_data['TDS'] = TDS
            system_data['dTDS'] = dTDS

        # If we have a free energy but not a Ka, compute Ka
        if not 'Ka' in system_data:
            try:
                system_data['Ka'], system_data['dKa'] = compute_Ka(system_data['DG'], system_data['dDG'])
            except TypeError:
                if system_data['DG']=='NaN':
                    system_data['Ka']='NaN'
                    system_data['dKa']='NaN'

        # Strip units.
        strip_units(system_data)

        # Consistency checks.
        if system_data['TDS']!='NaN' and system_data['DG']!='NaN' and system_data['DH']!='NaN':
            assert np.isclose(system_data['DG'], system_data['DH'] - system_data['TDS'], atol=0.10000000000001, rtol=0.0)

            if DG is not None:
                computed_DG = DG.value_in_unit(u.kilocalories_per_mole)
                assert np.isclose(np.around(computed_DG, decimals=2), system_data['DG'], atol=0.0200000000000001, rtol=0.0)
            if TDS is not None:
                computed_TDS = TDS.value_in_unit(u.kilocalories_per_mole)
                assert np.isclose(np.around(computed_TDS, decimals=2), system_data['TDS'], atol=0.0200000000000001, rtol=0.0)


        # Report only error most significant digit.
        for k in ['Ka', 'DH', 'TDS', 'DG']:
            if k in system_data:
                quantity, uncertainty = system_data[k], system_data['d' + k]
                if uncertainty is not None:
                    system_data[k], system_data['d' + k] = reduce_to_first_significant_digit(quantity, uncertainty)

    # Create output JSON file.
    with open('experimental_measurements.json', 'w') as f:
        json.dump(output_dict, f)

    # Create output CSV file.
    # Convert single dict to list of dicts.
    csv_dicts = []
    for system_id, system_data in output_dict.items():
        csv_dict = OrderedDict([('ID', system_id)])
        csv_dict.update(system_data)
        csv_dicts.append(csv_dict)
    with open('experimental_measurements.csv', 'w') as f:
        writer = csv.DictWriter(f, csv_dicts[0].keys(), delimiter=';')
        writer.writeheader()
        writer.writerows(csv_dicts)

    # Create a LaTex table.
    os.makedirs('PDFTable', exist_ok=True)
    old_host = ''
    with open('PDFTable/experimental_measurements.tex', 'w', encoding='utf-8') as f:
        f.write('\\documentclass{article}\n'
                '\\usepackage[a4paper,margin=0.4in,tmargin=0.5in,landscape]{geometry}\n'
                '\\usepackage{tabu}\n'
                '\\pagenumbering{gobble}\n'
                '\\begin{document}\n'
                '\\begin{center}\n'
                '\\footnotesize\n'
                '\\begin{tabu}')

        # Cell alignment.
        field_names = ['ID', 'name', '$K_a$ (M$^{-1}$)', '$\\Delta G$ (kcal/mol) $^{(a)}$', '$\\Delta H$ (kcal/mol)', '$T\\Delta S$ (kcal/mol) $^{(b)}$', '$n$']
        f.write('{| ' + ' | '.join(['c' for _ in range(len(field_names))]) + ' |}\n')

        # Table header.
        f.write('\\hline\n')
        f.write('\\rowfont{\\bfseries} ' + ' & '.join(field_names) + ' \\\\\n')
        f.write('\\hline\n')

        # Print lines.
        for csv_dict in csv_dicts:

            # Separate hosts with a double horizontal line.
            host_name = csv_dict['ID'].split('-')[0]
            if host_name != old_host:
                f.write('\\hline\n')
                old_host = host_name

            #if csv_dict['ID']=='clip-g10' or 'OA-g5' in csv_dict['ID']:
            #    # One name can't be dealt with; reformat
            #    csv_dict['name'] = "Can't format in LaTeX"

            row = '{ID} & {name}'

            # add a superscript to reflect the different experiments
            superscript = ''
            if csv_dict['ID'] == 'WP6-G1' or csv_dict['ID'] == 'WP6-G3':
                superscript += '$^a$'
            elif csv_dict['ID'] == 'WP6-G2' or csv_dict['ID'] == 'WP6-G10':
                superscript += '$^b$'
            elif csv_dict['ID'] == 'WP6-G6' or csv_dict['ID'] == 'WP6-G13':
                superscript += '$^c$'
            elif csv_dict['ID'] == 'WP6-G7':
                superscript += '$^d$'
            elif csv_dict['ID'] == 'WP6-G4' or csv_dict['ID'] == 'WP6-G8':
                superscript += '$^e$'
            elif csv_dict['ID'] == 'WP6-G5' or csv_dict['ID'] == 'WP6-G9' or csv_dict['ID'] == 'WP6-G11':
                superscript += '$^f$'
            elif csv_dict['ID'] == 'WP6-G12':
                superscript += '$^g$'


            for k in ['Ka', 'DG', 'DH', 'TDS']:
                row += ' & '

                # Report Ka in scientific notation.
                if k == 'Ka':
                    if not isinstance(k, str):
                        first_significant_digit = math.floor(math.log10(abs(csv_dict['d' + k])))
                        csv_dict['d' + k] /= 10**first_significant_digit
                        csv_dict[k] /= 10**first_significant_digit
                        row += '('
                row += '{' + k + '} +- {d' + k + '}'
                if k == 'Ka':
                    if not isinstance(k, str):
                        row += ') $\\times$ 10'
                        if first_significant_digit != 1:
                            row += '$^{{{{{}}}}}$'.format(first_significant_digit)

                ## Check if we used the upperbound.
                #superscript = ''
                ## if k != 'DG' and csv_dict['ID'] in upper_bound_molecules[k]:
                ##     superscript += 'a'
                #if k == 'Ka':
                #    if csv_dict['n'] == 0.33:
                #        superscript += 'd'
                #    elif csv_dict['n'] == 0.5 or csv_dict['n'] == 2:
                #        superscript += 'c'
                #if superscript != '':
                #    row += ' $^{{(' + superscript + ')}}$'

            row += (' & {n: .2f} \\\\\n'
                    '\\hline\n')

            row = row.format(**csv_dict)

            # Escape underscores for latex formatting
            row = row.replace('_','\_')

            # Write
            f.write(row)

        f.write('\end{tabu}\end{center}\\vspace{5mm}\n'
                'All quantities are reported as point estimate +- statistical error from the ITC data fitting procedure. '
                '$^a$ Measure directly by ITC during the tiration of WP6 (0.1 mM) in the cell with guest (1 mM) in '
                'the syring. $^b$ Measured by the ITC competition titration by mixing G7 (0.2 mM) and WP6 (0.1 mM) '
                'in the cell with guest (1 mM) in the syringe. $^c$ Measured directly by the ITC during the titration'
                'of WP6(0.05 mM) in the cell with guest (0.5 mM) in the syringe. $^d$ Measured directly by the ITC during '
                'the titration of WP6 (0.2 mM) in the cell with guest (2 mM) in the syringe. $^e$ Measured directly by the'
                'ITC during the tiration of WP6 (0.5 mM) in the cell with guest (5 mM) in the syringe. $^f$ Measured '
                'directly by the ITC during the tiration of WP6 (1 mM) in the cell with guest (10 mM) in the syringe. '
                '$^g$ Measured by the ITC competition titration by mixing G7 (0.5 mM) and WP6 (0.1 mM) in the cell with '
                'guest (1 mM) in the syringe.\\\\\n'
                '($^a$) Statistical errors were propagated from the $K_a$ measurements. \\\\\n'
                '($^b$) All experiments were performed at 298 K. \\\\\n'
                '($^c$) Units of M$^{-2}$. \\\\\n'
                '($^d$) Units of M$^{-3}$.\n'
                '\end{document}\n')
