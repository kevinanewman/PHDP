"""
    Physical constants

"""

from phdp import *

constants = {
    'KH2O-gas': 3.50,
    'Mair_g/mol': 28.96559,
    'MC_g/mol': 12.0107,
    'MCO2_g/mol': 44.0095,
    'MO_g/mol': 15.9994,
    'MS_g/mol': 32.0655,
    'MN_g/mol': 14.0067,
    'MCO_g/mol': 28.0101,
    'MTHC_g/mol': 13.87539,
    'MCH4_g/mol': 16.0425,
    'MNMHC_g/mol': 13.87539,
    'MNOx_g/mol': 46.0055,
    'MN2O_g/mol': 44.0128,
    'MH_g/mol': 1.00794,
    'MCH4N2O_g/mol': 60.05526,
    'wCH4N2O_Mass Fraction of urea in DEF': 0.325,
    'EmfuelCref	MJ/kg': 49.3112,
    'wCref': 0.874,
    'MH2O_g/mol': 18.01528,
}


def update_constants():
    # DEF constants
    constants['Mdef'] = constants['MC_g/mol'] + constants['MH_g/mol'] * 17.85 + constants['MO_g/mol'] * 7.92 + \
                    constants['MN_g/mol'] * 2

    constants['wcDef'] = constants['MC_g/mol'] / constants['Mdef']
    constants['whDef'] = constants['MH_g/mol'] * 17.85 / constants['Mdef']
    constants['woDef'] = constants['MO_g/mol'] * 7.92 / constants['Mdef']
    constants['wnDef'] = constants['MN_g/mol'] * 2 / constants['Mdef']

    # Test Fuel constants
    constants['FuelHTCRAT_ratio'] = phdp_globals.test_data['Header']['FuelHTCRAT_ratio'].item()
    constants['Mfuel_g/mol'] = constants['MC_g/mol'] + constants['FuelHTCRAT_ratio'] * constants['MH_g/mol']
    constants['wcFuel'] = constants['MC_g/mol'] / constants['Mfuel_g/mol']
    constants['whFuel'] = constants['FuelHTCRAT_ratio'] * constants['MH_g/mol'] / constants['Mfuel_g/mol']

    print_dict(constants)
