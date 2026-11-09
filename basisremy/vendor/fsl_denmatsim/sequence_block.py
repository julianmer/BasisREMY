"""Sequence block representation and RF/gradient unit helpers."""

import numpy as np

from .constants import GYRO_MAG_RATIO

__all__ = [
    "SequenceBlock",
    "active_gradient_axes",
    "gradient_multiplier",
    "rf_multiplier",
    "simBlock",
]


def gradient_multiplier(units, nucleus):
    """Return multiplier to convert gradient units to T/m."""
    units = units.lower()
    if units == "mt":
        return 1e-3
    elif units == "t":
        return 1
    elif units == "hz":
        return 1 / (GYRO_MAG_RATIO[nucleus] * 1e6)
    else:
        raise ValueError(f'Unknown GradUnits {units}. must be one of "mT", "T" or "Hz".')


def rf_multiplier(units, nucleus):
    """Return multiplier to convert RF amplitude units to Hz."""
    units = units.lower()
    if units == "hz":
        return 1
    elif units == "t":
        return GYRO_MAG_RATIO[nucleus] * 1e6
    elif units == "mt":
        return GYRO_MAG_RATIO[nucleus] * 1e3
    elif units == "ut":
        return GYRO_MAG_RATIO[nucleus]
    else:
        raise ValueError(f"Unknown RFUnits {units}.")


def active_gradient_axes(gradient):
    """Return per-axis gradient magnitudes and active axes."""
    if gradient.ndim == 1:
        gradient_magnitude = gradient
    else:
        gradient_magnitude = np.sum(np.abs(gradient), axis=1)

    axes = np.nonzero(gradient_magnitude)[0]
    if axes.size == 0:
        axes = np.asarray([None])

    return gradient_magnitude, axes


class SequenceBlock:
    """One sequence segment: RF pulse, delay/rephase period, and optional coherence filter."""

    def __init__(self, RF, delay, rephase, cfilter, RFUnits="Hz", GUnits="T", nucleus="1H"):
        """Create a sequence block with normalised RF and gradient units."""
        self.nucleus = nucleus
        self.rf_units = RFUnits
        self.gradient_units = GUnits
        self.delay_duration = delay
        self.coherence_filter = cfilter

        grad_multiplier = gradient_multiplier(GUnits, nucleus)
        self.rephase_gradient = grad_multiplier * np.array(rephase)

        if RF is None:
            self.complex_pulse = None
            self.amplitude_pulse = None
            self.phase_pulse = None
            self.time_step = 0.0
            self.frequency_offset = 0.0
            self.phase_offset = 0.0
            self.slice_gradient = np.array([0.0, 0.0, 0.0])
            self.gradient_magnitude = self.slice_gradient
            self.gradient_axis = np.asarray([None])
            return

        pulse_multiplier = rf_multiplier(RFUnits, nucleus)
        if "ampScale" in RF.keys():
            pulse_multiplier = pulse_multiplier * RF["ampScale"]

        if "cmplx" in RF.keys():
            self.complex_pulse = pulse_multiplier * np.array(RF["cmplx"])
            npoints = len(self.complex_pulse)
            self.amplitude_pulse = None
            self.phase_pulse = None
        elif "amp" in RF.keys():
            self.amplitude_pulse = pulse_multiplier * np.array(RF["amp"])
            npoints = len(self.amplitude_pulse)
            if "phase" in RF.keys():
                self.phase_pulse = np.array(RF["phase"])
            else:
                self.phase_pulse = 0.0
            self.complex_pulse = None
        else:
            raise ValueError('RF block must contain either "cmplx" or "amp".')

        self.time_step = RF["time"] / npoints
        self.frequency_offset = float(RF["frequencyOffset"])
        self.phase_offset = float(RF["phaseOffset"])
        self.slice_gradient = grad_multiplier * np.array(RF["grad"])
        self.gradient_magnitude, self.gradient_axis = active_gradient_axes(self.slice_gradient)

    def get_pulse(self):
        """Return complex RF pulse contained in block."""
        if (self.complex_pulse is None) and (self.amplitude_pulse is None):
            return None
        elif self.complex_pulse is None:
            rf = self.amplitude_pulse * np.exp(1j * self.phase_pulse)
        else:
            rf = self.complex_pulse

        return rf * np.exp(1j * self.phase_offset)

    def getPulse(self):
        """Backward-compatible alias for get_pulse."""
        return self.get_pulse()

    @property
    def cmplxPulse(self):
        return self.complex_pulse

    @cmplxPulse.setter
    def cmplxPulse(self, value):
        self.complex_pulse = value

    @property
    def ampPulse(self):
        return self.amplitude_pulse

    @ampPulse.setter
    def ampPulse(self, value):
        self.amplitude_pulse = value

    @property
    def phsPulse(self):
        return self.phase_pulse

    @phsPulse.setter
    def phsPulse(self, value):
        self.phase_pulse = value

    @property
    def tStep(self):
        return self.time_step

    @tStep.setter
    def tStep(self, value):
        self.time_step = value

    @property
    def offset(self):
        return self.frequency_offset

    @offset.setter
    def offset(self, value):
        self.frequency_offset = value

    @property
    def phase(self):
        return self.phase_offset

    @phase.setter
    def phase(self, value):
        self.phase_offset = value

    @property
    def sliceGrad(self):
        return self.slice_gradient

    @sliceGrad.setter
    def sliceGrad(self, value):
        self.slice_gradient = value
        self.gradient_magnitude, self.gradient_axis = active_gradient_axes(self.slice_gradient)

    @property
    def sum_grads(self):
        return self.gradient_magnitude

    @sum_grads.setter
    def sum_grads(self, value):
        self.gradient_magnitude = value

    @property
    def axis(self):
        return self.gradient_axis

    @axis.setter
    def axis(self, value):
        if value is None:
            self.gradient_axis = np.asarray([None])
        elif np.isscalar(value):
            self.gradient_axis = np.asarray([value])
        else:
            self.gradient_axis = np.asarray(value)

    @property
    def rephase(self):
        return self.rephase_gradient

    @rephase.setter
    def rephase(self, value):
        self.rephase_gradient = value

    @property
    def delay(self):
        return self.delay_duration

    @delay.setter
    def delay(self, value):
        self.delay_duration = value

    @property
    def cfilter(self):
        return self.coherence_filter

    @cfilter.setter
    def cfilter(self, value):
        self.coherence_filter = value


# Backward-compatible name for existing callers.
simBlock = SequenceBlock
