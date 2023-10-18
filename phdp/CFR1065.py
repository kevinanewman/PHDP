"""

CFR 1065 code

"""

import numpy as np
from common.phdp_functions import ASTM_round
from constants import constants


def vapor_pressure_of_water_kPa(Tsat_K):
    """
    CFR 1065.645-1
    https://www.ecfr.gov/current/title-40/chapter-I/subchapter-U/part-1065/subpart-G/section-1065.645

    Calculate the vapor pressure of water for a given saturation temperature condition, Tsat, as follows,
    or use good engineering judgment to use a different relationship of the vapor pressure of water to a given
    saturation temperature condition:

    (1) For humidity measurements made at ambient temperatures from (0 to 100) °C, or for humidity measurements made
     over super-cooled water at ambient temperatures from (−50 to 0) °C, use the equation below.

    Args:
        Tsat_K (float): saturation temperature of water at measured conditions, K

    Returns:
        Vapor pressure of water at saturation temperature condition, kPa

    """
    return 10 ** ((10.79574 * (1 - (273.16 / Tsat_K))) -
                  (5.028 * np.log10(Tsat_K / 273.16)) +
                  (0.000150475 * (1 - (10 ** (-8.2969 * ((Tsat_K / 273.16) - 1)))))
                  + (0.00042873 * (10 ** (4.76955 * (1 - (273.16 / Tsat_K))) - 1))
                  - 0.2138602)


def relative_humidity(RH_pct, pH2O_kPa, pABS_kPa):
    """
    CFR 1065.645-4
    https://www.ecfr.gov/current/title-40/chapter-I/subchapter-U/part-1065/subpart-G/section-1065.645

    If you measure humidity as a relative humidity, RH,
    determine the amount of water in an ideal gas, xH2O, as follows.

    Args:
        RH_pct (float): percent relative humidity, e.g. ``48.89``
        pH2O (float): vapor pressure of water at measured conditions
        pABS_kPa (float): wet static absolute pressure at measured conditions

    Returns:
        xH2O mol/mol, relative humidity

    """
    RH = RH_pct / 100

    return RH * pH2O_kPa / pABS_kPa


def dewpoint_temp_K(pH2O_scaled):
    """
    CFR 1065.645-5
    https://www.ecfr.gov/current/title-40/chapter-I/subchapter-U/part-1065/subpart-G/section-1065.645

    Dewpoint determination from relative humidity and dry bulb temperature.  This paragraph (d) describes how to
    calculate dewpoint temperature from relative humidity, RH. This is based on “ITS–90 Formulations for Vapor Pressure,
    Frostpoint Temperature, Dewpoint Temperature, and Enhancement Factors in the Range −100 to + 100 °C”
    (Hardy, B., The Proceedings of the Third International Symposium on Humidity & Moisture, Teddington, London,
    England, April 1998). Calculate pH20sat as described in paragraph (a) of this section based on setting Tsat equal
    to Tamb. Calculate pH20scaled by multiplying pH20sat by RH. Calculate the dewpoint, Tdew, from pH20
    using the equation below.

    Args:
        pH2O_scaled (float): water vapor pressure scaled to the relative humidity at the location of the
            relative humidity measurement, Tsat = Tamb

    Returns:
        The dewpoint temperature, K

    """
    return (207.98233 - 20.156028 * np.log(pH2O_scaled) + 0.46778925 *
            (np.log(pH2O_scaled)) ** 2 - 9.2288067 * 10 ** (-6) * (np.log(pH2O_scaled)) ** 3) / \
        (1 - 0.13319669 * np.log(pH2O_scaled) +
         5.6577518 * 10 ** (-3) * (np.log(pH2O_scaled)) ** 2
         - 7.5172865 * 10 ** (-5) * (np.log(pH2O_scaled)) ** 3)


def alpha(qmFuel_gph, DEFMassFlowRate_Avg_gph):
    """
    CFR 1065.655-20
    https://www.ecfr.gov/current/title-40/chapter-I/subchapter-U/part-1065/subpart-G/section-1065.655

    Calculate alpha, atomic hydrogen-to-carbon ratio

    Args:
        qmFuel_gph (Series): fuel mass flow rate, grams per hour
        DEFMassFlowRate_Avg_gph (Series): DEF mass flow rate, grams per hour

    Returns:
        alpha, atomic hydrogen-to-carbon ratio

    """
    fuel_rate_negative = qmFuel_gph < 0

    alpha = (constants['MC_g/mol'] / constants['MH_g/mol'] *
       (qmFuel_gph * constants['whFuel'] + DEFMassFlowRate_Avg_gph * constants['whDef']) /
       (qmFuel_gph * constants['wcFuel'] + DEFMassFlowRate_Avg_gph * constants['wcDef']))

    alpha.loc[fuel_rate_negative] = 1.86  # for now, eventually probably: alpha.loc[~fuel_rate_negative].mean()

    return alpha


def beta(qmFuel_gph, DEFMassFlowRate_Avg_gph):
    """
    CFR 1065.655-20
    https://www.ecfr.gov/current/title-40/chapter-I/subchapter-U/part-1065/subpart-G/section-1065.655

    Calculate beta, atomic oxygen-to-carbon ratio

    Args:
        qmFuel_gph (Series): fuel mass flow rate, grams per hour
        DEFMassFlowRate_Avg_gph (Series): DEF mass flow rate, grams per hour

    Returns:
        beta, atomic oxygen-to-carbon ratio

    """
    fuel_rate_negative = qmFuel_gph < 0

    beta = (constants['MC_g/mol'] / constants['MO_g/mol'] *
            (DEFMassFlowRate_Avg_gph * constants['woDef']) /
            (qmFuel_gph * constants['wcFuel'] + DEFMassFlowRate_Avg_gph * constants['wcDef']))

    beta.loc[fuel_rate_negative] = 0.03  # for now, eventually probably: beta.loc[~fuel_rate_negative].mean()

    return beta


def delta(qmFuel_gph, DEFMassFlowRate_Avg_gph):
    """
    CFR 1065.655-20
    https://www.ecfr.gov/current/title-40/chapter-I/subchapter-U/part-1065/subpart-G/section-1065.655

    Calculate delta, atomic nitrogen-to-carbon ratio

    Args:
        qmFuel_gph (Series): fuel mass flow rate, grams per hour
        DEFMassFlowRate_Avg_gph (Series): DEF mass flow rate, grams per hour

    Returns:
        delta, atomic nitrogen-to-carbon ratio

    """
    fuel_rate_negative = qmFuel_gph < 0

    delta = (constants['MC_g/mol'] / constants['MN_g/mol'] *
        (DEFMassFlowRate_Avg_gph * constants['wnDef'])/
             (qmFuel_gph * constants['wcFuel'] + DEFMassFlowRate_Avg_gph * constants['wcDef']))

    delta.loc[fuel_rate_negative] = 0.007  # for now, eventually probably: delta.loc[~fuel_rate_negative].mean()

    return delta


if __name__ == '__main__':
    test_pass = True

    pH2O = ASTM_round(vapor_pressure_of_water_kPa(282.65), 6)
    test_pass = test_pass and pH2O == 1.186581

    xH2O = ASTM_round(relative_humidity(50.77, 2.3371, 99.980), 6)
    test_pass = test_pass and xH2O == 0.011868

    Tdew_K = ASTM_round(dewpoint_temp_K(925.717), 6)
    test_pass = test_pass and Tdew_K == 279.001023

    if not test_pass:
        os._exit(-1)
