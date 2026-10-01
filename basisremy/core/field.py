####################################################################################################
#                                             field.py                                             #
####################################################################################################
#                                                                                                  #
# Purpose: One rule for the static field. The spectrometer frequency is the primary quantity       #
#          when a header carries it (REMY rounds B0 to two decimals, so a 123.25 MHz Siemens       #
#          system would otherwise simulate at 2.89 T = 123.05 MHz); the field is derived from it   #
#          unrounded. Without a frequency the field gives the frequency. Every backend and the     #
#          exporters then see a consistent pair.                                                   #
#                                                                                                  #
####################################################################################################

# Gyromagnetic ratios in MHz/T (gamma / 2 pi). 1H uses the value the simulation
# engines use (FID-A, MRSCloud, Spinach adapters: 42.577), so B0 -> MHz -> B0
# round-trips exactly through them.
GAMMA_MHZ_PER_T = {
    '1H': 42.577, '2H': 6.536, '13C': 10.7084, '15N': -4.316, '17O': -5.772,
    '19F': 40.078, '23Na': 11.262, '31P': 17.235,
}


def gamma_mhz_per_t(nucleus):
    """Gyromagnetic ratio for a nucleus label such as '1H' (default 1H)."""
    return GAMMA_MHZ_PER_T.get(str(nucleus).strip() if nucleus else '1H',
                               GAMMA_MHZ_PER_T['1H'])


def _number(value):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if v > 0 else None


def reconcile_field(info):
    """Return the ('Center Freq' [MHz], 'B0' [T]) pair implied by *info*.

    The frequency wins when both are present; the result is empty when
    neither is a number.
    """
    gamma = gamma_mhz_per_t(info.get('Nucleus'))
    freq = _number(info.get('Center Freq'))
    b0 = _number(info.get('B0'))
    if freq is not None:
        return {'Center Freq': freq, 'B0': freq / gamma}
    if b0 is not None:
        return {'Center Freq': b0 * gamma, 'B0': b0}
    return {}
