"""Propagation routines for sequence simulation modes."""

from collections import Counter

import numpy as np

from .simulator import simulator

MAX_FULL_PROPAGATOR_CACHE_UNIQUE = 256
MAX_DELAY_PROPAGATOR_CACHE_UNIQUE = 256
MAX_REPEATED_PULSE_PROPAGATOR_CACHE_UNIQUE = 256


# Convert (rephase) gradient areas to amplitudes for a given duration
def areaToAmp(areas, durations):
    """Convert gradient area to amplitude given duration."""
    return areas / durations


def rephase_amplitude_for_delay(block):
    """Return rephase gradient amplitude, or None for a zero-duration no-op."""
    if block.delay == 0:
        if np.any(block.rephase != 0):
            raise ValueError("Cannot apply a non-zero rephase gradient with zero delay.")
        return None

    return areaToAmp(block.rephase, block.delay)


def apply_rephase_delay(obj, p, block, currR, rephase_amplitude, propagator_cache=None, cache_key=None):
    """Apply delay/rephase propagation when the block has a non-zero delay."""
    if rephase_amplitude is None:
        return p

    if propagator_cache is None:
        return obj.applyGrad(rephase_amplitude, currR, block.delay, p=p)

    prop = get_delay_rephase_propagator(obj, block, currR, rephase_amplitude, propagator_cache, cache_key)
    return obj.applyPropagator(prop, p=p)


def apply_coherence_filter(obj, p, block):
    """Apply a block coherence filter if present."""
    if block.cfilter is None:
        return p

    return obj.selectCoherence(block.cfilter, p=p)


def build_full_spatial_coordinates(spatialDims):
    """Return flattened x/y/z coordinates for full spatial simulation."""
    X, Y, Z = np.meshgrid(*spatialDims, indexing="ij")
    return np.stack((X, Y, Z), axis=-1).reshape(-1, 3)


def effective_gradient_values(gradient, coordinate):
    """Return the gradient-position product for each gradient time point."""
    if gradient.ndim == 1:
        return np.asarray([gradient @ coordinate])

    return gradient.T @ coordinate


def effective_gradient_key(gradient, coordinate, decimals=12):
    """Return a rounded cache key for the effective gradient offset."""
    return tuple(np.round(effective_gradient_values(gradient, coordinate), decimals=decimals))


def array_value_key(values):
    """Return a hashable value key for a numpy-compatible array."""
    if values is None:
        return None

    values = np.ascontiguousarray(values)
    return values.shape, values.dtype.str, values.tobytes()


def pulse_signature(block):
    """Return the parts of a block that determine its RF Hamiltonian."""
    pulse = block.getPulse()
    if pulse is None or block.tStep == 0:
        return None

    return array_value_key(pulse), float(block.offset), float(block.tStep)


class RepeatedPulsePropagatorCache:
    """Cache RF+slice propagators for pulse blocks repeated within one simulation."""

    def __init__(self, seqBlocks, max_unique=MAX_REPEATED_PULSE_PROPAGATOR_CACHE_UNIQUE):
        self._block_signatures = {id(block): pulse_signature(block) for block in seqBlocks}
        signature_counts = Counter(signature for signature in self._block_signatures.values() if signature is not None)
        self._repeated_signatures = {signature for signature, count in signature_counts.items() if count > 1}
        self._cache = {}
        self._disabled_signatures = set()
        self._max_unique = max_unique

    def is_reusable(self, block):
        """Return True if this block has a repeated pulse signature worth caching."""
        signature = self._block_signatures.get(id(block))
        return signature in self._repeated_signatures and signature not in self._disabled_signatures

    def get(self, block, gradient_key):
        """Return a cached propagator for block and gradient key, if available."""
        signature = self._block_signatures.get(id(block))
        if signature not in self._repeated_signatures or signature in self._disabled_signatures:
            return None

        return self._cache.get(signature, {}).get(gradient_key)

    def store(self, block, gradient_key, propagator):
        """Store a propagator if this repeated pulse still fits in the cache."""
        signature = self._block_signatures.get(id(block))
        if signature not in self._repeated_signatures or signature in self._disabled_signatures:
            return

        signature_cache = self._cache.setdefault(signature, {})
        if gradient_key not in signature_cache and len(signature_cache) >= self._max_unique:
            self._disabled_signatures.add(signature)
            self._cache.pop(signature, None)
            return

        signature_cache[gradient_key] = propagator


def reusable_propagator_keys(gradient, coordinates, max_unique):
    """Return reusable propagator keys, or None when caching is not worthwhile."""
    keys = [effective_gradient_key(gradient, currR) for currR in coordinates]
    unique_count = len(set(keys))
    if unique_count == len(keys):
        return None
    if unique_count > max_unique:
        return None
    if unique_count > (len(keys) // 2):
        return None

    return keys


def reusable_full_propagator_keys(block, full_coordinates):
    """Return reusable RF+slice propagator keys, or None when caching is not worthwhile."""
    return reusable_propagator_keys(block.sliceGrad, full_coordinates, MAX_FULL_PROPAGATOR_CACHE_UNIQUE)


def reusable_delay_rephase_keys(rephase_amplitude, coordinates):
    """Return reusable delay/rephase propagator keys, or None when caching is not worthwhile."""
    if rephase_amplitude is None:
        return None

    return reusable_propagator_keys(rephase_amplitude, coordinates, MAX_DELAY_PROPAGATOR_CACHE_UNIQUE)


def get_rf_slice_propagator(obj, block, Hconst, currR, propagator_cache=None, cache_key=None, pulse_cache=None):
    """Return a cached or newly-created RF+slice propagator for one coordinate."""
    if pulse_cache is not None and pulse_cache.is_reusable(block) and cache_key is None:
        cache_key = effective_gradient_key(block.sliceGrad, currR)

    if pulse_cache is not None and cache_key is not None:
        prop = pulse_cache.get(block, cache_key)
        if prop is not None:
            return prop

    if propagator_cache is not None and cache_key is not None:
        prop = propagator_cache.get(cache_key)
        if prop is not None:
            if pulse_cache is not None:
                pulse_cache.store(block, cache_key, prop)
            return prop

    H = Hconst + obj.getHGrad(block.sliceGrad, currR)
    prop = simulator.createPropagator(H, block.tStep)
    if propagator_cache is not None and cache_key is not None:
        propagator_cache[cache_key] = prop
    if pulse_cache is not None and cache_key is not None:
        pulse_cache.store(block, cache_key, prop)

    return prop


def get_full_rf_slice_propagator(obj, block, Hconst, currR, propagator_cache, cache_key):
    """Return a cached or newly-created RF+slice propagator for one full-mode coordinate."""
    return get_rf_slice_propagator(obj, block, Hconst, currR, propagator_cache, cache_key)


def get_delay_rephase_propagator(obj, block, currR, rephase_amplitude, propagator_cache, cache_key):
    """Return a cached or newly-created delay/rephase propagator."""
    if cache_key is None:
        H = obj.Hcs + obj.Hj + obj.getHGrad(rephase_amplitude, currR)
        return simulator.createPropagator(H, block.delay)

    prop = propagator_cache.get(cache_key)
    if prop is None:
        H = obj.Hcs + obj.Hj + obj.getHGrad(rephase_amplitude, currR)
        prop = simulator.createPropagator(H, block.delay)
        propagator_cache[cache_key] = prop

    return prop


def copy_density_to_grid(p, grid_shape):
    """Return density matrix p copied across a spatial grid."""
    return np.broadcast_to(p, tuple(grid_shape) + p.shape).astype(np.complex128, copy=True)


def copy_density_to_points(p, npoints):
    """Return density matrix p copied across a flat point dimension."""
    return np.broadcast_to(p, (npoints,) + p.shape).astype(np.complex128, copy=True)


def build_axis_coordinates(spatial_axis, axis):
    """Return position vectors with values on one active axis."""
    coordinates = np.zeros((spatial_axis.size, 3))
    coordinates[:, axis] = spatial_axis
    return coordinates


def iter_interleaved_coordinates(spatialDims, innerDim):
    """Yield density-matrix indices, coordinates, and propagator indices for interleaved loops."""
    outerdims = np.arange(0, 3)
    outerdims = outerdims[outerdims != innerDim]

    p_index = [0, 0, 0]
    currR = np.zeros(3)
    for iDx, r1 in enumerate(spatialDims[outerdims[0]]):
        p_index[outerdims[0]] = iDx
        currR[outerdims[0]] = r1
        for jDx, r2 in enumerate(spatialDims[outerdims[1]]):
            p_index[outerdims[1]] = jDx
            currR[outerdims[1]] = r2
            for kDx, r3 in enumerate(spatialDims[innerDim]):
                p_index[innerDim] = kDx
                currR[innerDim] = r3
                yield tuple(p_index), currR, kDx


def propagate_1d_sequence(obj, seqBlocks, axisOrder, spatialDims, pulse_cache=None):
    """Run sequence blocks using the 1D projection method."""
    if pulse_cache is None:
        pulse_cache = RepeatedPulsePropagatorCache(seqBlocks)

    currentP = obj.thermalEq()
    previousAx = None
    coordinate_cache = {}
    for ax, block in zip(axisOrder, seqBlocks):
        if ax is None:
            ax = 0
        if ax not in coordinate_cache:
            coordinate_cache[ax] = build_axis_coordinates(spatialDims[ax], ax)

        # Loop over spatial dimension relevant to this pulse
        if currentP.ndim == 2:  # Single density mat
            currentP = copy_density_to_grid(currentP, spatialDims[ax].shape)
        # If this axis is not the same as the last then take the mean and redistribute over new axis.
        elif previousAx != ax:
            currentP = copy_density_to_grid(np.mean(currentP, axis=0), spatialDims[ax].shape)

        # Precalculate the RF Hamiltonian so that it isn't repeated unnecessarily for each point
        # import pdb; pdb.set_trace()
        Hrf = obj.getHRF(block.getPulse(), block.offset)
        Hconst = Hrf + obj.Hcs + obj.Hj
        rph = rephase_amplitude_for_delay(block)
        delay_keys = reusable_delay_rephase_keys(rph, coordinate_cache[ax])
        delay_cache = None if delay_keys is None else {}
        for iDx, currR in enumerate(coordinate_cache[ax]):
            # Pulse
            prop = get_rf_slice_propagator(obj, block, Hconst, currR, pulse_cache=pulse_cache)
            currentP[iDx] = obj.applyPropagator(prop, p=currentP[iDx])

            # Delay with rephase
            delay_key = None if delay_keys is None else delay_keys[iDx]
            currentP[iDx] = apply_rephase_delay(
                obj,
                currentP[iDx],
                block,
                currR,
                rph,
                delay_cache,
                delay_key,
            )

            # Coherence selection
            currentP[iDx] = apply_coherence_filter(obj, currentP[iDx], block)
        previousAx = ax

    return currentP


def propagate_interleaved_sequence(obj, seqBlocks, axisOrder, spatialDims, pulse_cache=None):
    """Run sequence blocks using the interleaved spatial simulation method."""
    if pulse_cache is None:
        pulse_cache = RepeatedPulsePropagatorCache(seqBlocks)

    x, y, z = spatialDims
    currentP = obj.thermalEq()
    axis_coordinate_cache = {}
    for ax, block in zip(axisOrder, seqBlocks):
        if ax is None:
            ax = 0
        if ax not in axis_coordinate_cache:
            axis_coordinate_cache[ax] = build_axis_coordinates(spatialDims[ax], ax)

        # Hybrid method. We can't use 1d projection method but we know there aren't
        # any concurrent slice select gradients.
        # We can therefore precalculate the pulse propagators for a column of isochromats.
        # Then repeatedly apply this column to all values of the second and third dimensions.

        # Make an array of density matrices of the correct shape
        if currentP.ndim == 2:
            fullshape = (x.size, y.size, z.size)
            currentP = copy_density_to_grid(currentP, fullshape)

        # Make propagators for rf pulse + slice select gradient over current axis
        currProp = []
        Hrf = obj.getHRF(block.getPulse(), block.offset)
        Hconst = Hrf + obj.Hcs + obj.Hj
        for currR in axis_coordinate_cache[ax]:
            currProp.append(get_rf_slice_propagator(obj, block, Hconst, currR, pulse_cache=pulse_cache))

        # This next step isn't pretty! But I couldn't figure out how to vectorise the step in the numpy framework
        rph = rephase_amplitude_for_delay(block)
        delay_keys = reusable_delay_rephase_keys(
            rph,
            (currR for _, currR, _ in iter_interleaved_coordinates(spatialDims, ax)),
        )
        delay_cache = None if delay_keys is None else {}
        for delay_idx, (p_index, currR, prop_index) in enumerate(iter_interleaved_coordinates(spatialDims, ax)):
            iDx, jDx, kDx = p_index
            tmpP = currentP[iDx, jDx, kDx, :, :]
            # Pulse
            tmpP = obj.applyPropagator(currProp[prop_index], p=tmpP)

            # Rephase+delay
            delay_key = None if delay_keys is None else delay_keys[delay_idx]
            tmpP = apply_rephase_delay(obj, tmpP, block, currR, rph, delay_cache, delay_key)

            # Coherence selection
            tmpP = apply_coherence_filter(obj, tmpP, block)

            currentP[iDx, jDx, kDx, :, :] = tmpP

    return currentP


def propagate_full_sequence(obj, seqBlocks, axisOrder, spatialDims, pulse_cache=None):
    """Run sequence blocks using the full spatial simulation method."""
    if pulse_cache is None:
        pulse_cache = RepeatedPulsePropagatorCache(seqBlocks)

    full_coordinates = build_full_spatial_coordinates(spatialDims)
    currentP = obj.thermalEq()
    for _, block in zip(axisOrder, seqBlocks):
        # Loop over mesh grid of all spatial points
        # This could take a very long time!
        if currentP.ndim == 2:
            currentP = copy_density_to_points(currentP, full_coordinates.shape[0])

        Hrf = obj.getHRF(block.getPulse(), block.offset)
        Hconst = Hrf + obj.Hcs + obj.Hj
        rph = rephase_amplitude_for_delay(block)
        propagator_keys = reusable_full_propagator_keys(block, full_coordinates)
        propagator_cache = {}
        delay_keys = reusable_delay_rephase_keys(rph, full_coordinates)
        delay_cache = None if delay_keys is None else {}
        for iDx, currR in enumerate(full_coordinates):
            # Pulse
            cache_key = None if propagator_keys is None else propagator_keys[iDx]
            prop = get_rf_slice_propagator(
                obj,
                block,
                Hconst,
                currR,
                propagator_cache,
                cache_key,
                pulse_cache,
            )
            currentP[iDx] = obj.applyPropagator(prop, p=currentP[iDx])

            # Delay with rephase
            delay_key = None if delay_keys is None else delay_keys[iDx]
            currentP[iDx] = apply_rephase_delay(
                obj,
                currentP[iDx],
                block,
                currR,
                rph,
                delay_cache,
                delay_key,
            )

            # Coherence selection
            currentP[iDx] = apply_coherence_filter(obj, currentP[iDx], block)

    return currentP


def propagate_sequence(obj, seqBlocks, axisOrder, spatialDims, method):
    """Run sequence blocks and return the final density matrix array before readout."""
    pulse_cache = RepeatedPulsePropagatorCache(seqBlocks)
    if method == "1d":
        return propagate_1d_sequence(obj, seqBlocks, axisOrder, spatialDims, pulse_cache)
    elif method == "interleaved":
        return propagate_interleaved_sequence(obj, seqBlocks, axisOrder, spatialDims, pulse_cache)
    elif method == "full":
        return propagate_full_sequence(obj, seqBlocks, axisOrder, spatialDims, pulse_cache)
    else:
        raise ValueError(f"Unknown simulation method {method}.")
