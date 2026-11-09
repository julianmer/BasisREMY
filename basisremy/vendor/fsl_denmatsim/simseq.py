import numpy as np
from copy import deepcopy

from .constants import GYRO_MAG_RATIO
from . import propagation
from .sequence_block import SequenceBlock, simBlock as simBlock
from . import simulator as sim

# Backward-compatible names for callers that reached these helpers via simseq.
areaToAmp = propagation.areaToAmp
apply_rephase_delay = propagation.apply_rephase_delay
build_full_spatial_coordinates = propagation.build_full_spatial_coordinates
propagate_sequence = propagation.propagate_sequence
rephase_amplitude_for_delay = propagation.rephase_amplitude_for_delay


def normalise_spin_system(spinsys):
    """Return a validated copy of a spin system with numpy array fields."""
    out = dict(spinsys)

    if "shifts" in out:
        if isinstance(out["shifts"], (float, int)):
            out["shifts"] = np.array([out["shifts"]], dtype=float)
        else:
            out["shifts"] = np.array(out["shifts"], dtype=float)
    else:
        raise ValueError("spinsys must have a field shifts.")

    if "J" in out:
        out["j"] = out.pop("J")
    if "j" in out:
        if isinstance(out["j"], (float, int)):
            out["j"] = np.array([[out["j"]]], dtype=float)
        else:
            out["j"] = np.array(out["j"], dtype=float)
    else:
        raise ValueError("spinsys must have a field j.")

    return out


def normalise_resolution(resolution):
    """Return resolution as a positive three-element integer array."""
    resolution = np.asarray(resolution)
    if resolution.ndim == 0:
        resolution = np.repeat(resolution, 3)

    if resolution.shape != (3,):
        raise ValueError(f"resolution must be a scalar or a three-element sequence, got {resolution}.")

    if not np.issubdtype(resolution.dtype, np.integer):
        if not np.all(np.equal(resolution, np.asarray(resolution, dtype=int))):
            raise ValueError(f"resolution entries must be integers, got {resolution}.")

    resolution = np.asarray(resolution, dtype=int)
    if (resolution < 1).any():
        raise ValueError(f"resolution entries must be >= 1, got {resolution}.")

    return resolution


def make_spatial_axis(bounds, resolution, units):
    """Return a spatial axis in metres, sampling the centre for singleton axes."""
    bounds = np.asarray(bounds, dtype=float)
    if bounds.shape != (2,):
        raise ValueError(f"Spatial bounds must be a two-element sequence, got {bounds}.")

    resolution = int(resolution)
    if resolution < 1:
        raise ValueError(f"Spatial resolution must be >= 1, got {resolution}.")

    if units.lower() == "mm":
        spaceScale = 1000
    elif units.lower() == "cm":
        spaceScale = 100
    elif units.lower() == "m":
        spaceScale = 1
    else:
        raise ValueError(f"Unknown spatial units {units}. Use mm, cm or m.")

    if resolution == 1:
        axis = np.array([np.mean(bounds)])
    else:
        axis = np.linspace(bounds[0], bounds[1], resolution)

    return axis / spaceScale


def normalise_sequence_params(parameterdict):
    """Return a validated copy of sequence parameters."""
    out = deepcopy(parameterdict)
    out = grumble(out)
    out["resolution"] = normalise_resolution(out["resolution"])

    if "Nucleus" in out:
        if out["Nucleus"] not in GYRO_MAG_RATIO:
            raise ValueError(
                f"Unknown nucleus {out['Nucleus']}. Only {', '.join(list(GYRO_MAG_RATIO.keys()))} are supported."
            )
    else:
        out["Nucleus"] = "1H"

    if "method" in out:
        out["method"] = out["method"].lower()
        if out["method"] not in {"1d", "interleaved", "full"}:
            raise ValueError('method must be one of "1d", "interleaved", or "full".')

    return out


# Grumble over any inconsistencies in the passed parameters
def grumble(parameterdict):
    """Function to handle and clean sequence parameters passed to seqsim."""
    pulses = parameterdict["RF"]
    R = np.array(parameterdict["rephaseAreas"])
    delays = np.array(parameterdict["delays"])
    coherenceFilter = parameterdict["CoherenceFilter"]

    # Check if number of RF, delays, gradients,rephase areas etc are consistent.
    # And otherwise have consistent shapes (3 axes)
    checkList = [R, delays, coherenceFilter]
    checkListNames = ["rephaseAreas", "delays", "CoherenceFilter"]
    simSteps = len(pulses)
    for param, pName in zip(checkList, checkListNames):
        currentLength = len(param)
        if currentLength != simSteps:
            raise ValueError(
                f"{simSteps} RF pulses were specified but {pName} has {currentLength} entries ({param}). Aborting."
            )

        if any(x in pName for x in ["rephaseAreas"]):
            if param.shape[1] != 3:
                raise ValueError(
                    f"{pName} must be a NRFx3 ({simSteps}x3) matrix. It has shape {param.shape}. Aborting."
                )

    # Check RF has required fields if not then add those that are possible
    for iDx, rf in enumerate(pulses):
        # Throw errors
        if "time" not in rf.keys():
            raise ValueError(f"RF {iDx} does not contain a time field. Aborting.")
        if "grad" not in rf.keys():
            raise ValueError(f"RF {iDx} does not contain a grad field. Aborting.")

        # Add default values
        if "frequencyOffset" not in rf.keys():
            parameterdict["RF"][iDx].update({"frequencyOffset": 0.0})
        if "phaseOffset" not in rf.keys():
            parameterdict["RF"][iDx].update({"phaseOffset": 0.0})

        # Check possible times
        if (np.array(rf["time"]) < 0).any():
            raise ValueError(f"pulseTime in RF {iDx} must have values >= 0.0.")

    # If some of the optional parameters aren't passed then generate them here
    # Optional parameters with defaults are:
    optParams = ["spatiallyResolve", "centralShift", "RFUnits", "GradUnits", "spaceUnits"]
    defaults = [False, 0.0, "Hz", "T", "m"]
    for optP, val in zip(optParams, defaults):
        if optP not in parameterdict.keys():
            parameterdict.update({optP: val})

    # Check time parameters are not less than zero
    if (np.array(delays) < 0).any():
        raise ValueError(f"Delays must have values >= 0.0. Delays = {delays}.")

    # Check delays don't have 0 time if there is a rephase gradient
    for r, d in zip(R, delays):
        if (r != 0).any() and (d == 0):
            raise ValueError(f"Delay has a value of 0 when there is a rephase gradient in same block (rephase = {r}).")

    return parameterdict


def grumblespins(spinsys):
    """Function to handle and clean spin system parameters passed to seqsim."""
    return normalise_spin_system(spinsys)


# This function checks the order of the entries of the slice-select gradient variable and determines
# what mode the simulator runs in.
# If the gradients occur on axes in an order without reoccurrence e.g 2,1,1,3 then 1D projection may be used
# If they occur in an interleaved order e.g. 1,2,1,3 then use an intermediate 'interleaved' mode.
# If gradients occur simultaneously across 2 or more axes then use the full simulation method. Slow!
def gradOrderChecker(seqBlocks):
    """Determine simulation mode from gradient order in sequence blocks."""
    fullSim = False
    disable1D = False
    axisUsed = []
    for block in seqBlocks:
        if block.axis.size > 1:
            print(
                "More than one gradient axis on per pulse (G).\n "
                "Setting to full spatial simulation mode.\n "
                "This might take a while!"
            )
            fullSim = True
            break

        currentAxis = block.axis

        if currentAxis[0] is None:
            axisUsed.append(currentAxis[0])
            continue
        if (axisUsed == currentAxis).any():  # Same axis used again, not a problem yet
            # Get the last axis used and guard against Nones (or all Nones)
            last_axis_used = None
            counter = 0
            while last_axis_used is None:
                counter -= 1
                if np.abs(counter) > len(axisUsed):
                    continue
                last_axis_used = axisUsed[counter]
            if (last_axis_used != currentAxis).any():  # Not the last used axis, now it's a problem
                disable1D = True
                print("Slice select gradients not in linear order, 1D projection method disabled.")

        axisUsed.append(currentAxis[0])

    axisUsed = np.array(axisUsed)
    if fullSim:
        return "full", axisUsed
    elif disable1D:
        return "interleaved", axisUsed
    else:
        # 1D is fine
        return "1d", axisUsed


# In the 1D projection method the simulator can't cope with rephase gradients applied on another axis.
# If this is the case (e.g. 2nd rephase in STEAM) then insert a new entry into all of the looping variables.
# This new entry has zeros for most entries, has an entry for the rephase in the different axis
# and shares the delay time with the previous axis.
def splitConcurrentRephase(seqBlock):
    """Split sequence blocks with multiple gradient directions occurring in the rephase section"""
    seqBlockOut = []
    axisOrderOut = []
    for iDx, block in enumerate(seqBlock):  # Loop over blocks
        numAxisReph = np.sum(block.rephase != 0)  # Number of elements in row of R that aren't 0
        dims = np.arange(0, 3)
        dims = dims[dims != block.axis]
        insertNeeded = (block.rephase[dims] != 0).any()
        # if numAxisReph>1: # If more than 1
        if insertNeeded:
            currentLength = len(seqBlockOut)  # Store current length
            stableBlock = deepcopy(block)
            for axis, rephG in enumerate(stableBlock.rephase):  # Loop over elements in row
                if block.sum_grads[axis] != 0:  # If on the same axis as slice select then insert at currentLength
                    # Modify values in block
                    toInsert = np.zeros(3)
                    toInsert[axis] = rephG
                    block.rephase = toInsert
                    block.delay = stableBlock.delay / numAxisReph
                    # Insert
                    seqBlockOut.insert(currentLength, block)
                    axisOrderOut.insert(currentLength, block.axis[0])

                elif rephG == 0:  # If zero and not on slice select axis just skip
                    continue
                else:  # Otherwise make new line and append
                    # import pdb; pdb.set_trace()
                    toInsert = np.zeros(3)
                    toInsert[axis] = rephG
                    tmpBlock = SequenceBlock(None, stableBlock.delay / numAxisReph, toInsert, None)
                    tmpBlock.axis = axis
                    seqBlockOut.append(tmpBlock)
                    axisOrderOut.append(axis)

        else:  # If just one value of rephase keep as is.
            seqBlockOut.append(block)
            if block.axis[0] is None:
                axisOrderOut.append(None)
            else:
                axisOrderOut.append(int(block.axis[0]))

    return seqBlockOut, np.asarray(axisOrderOut)


def build_spatial_dimensions(parameterdict):
    """Return x, y, z spatial axes in metres."""
    x = make_spatial_axis(parameterdict["x"], parameterdict["resolution"][0], parameterdict["spaceUnits"])
    y = make_spatial_axis(parameterdict["y"], parameterdict["resolution"][1], parameterdict["spaceUnits"])
    z = make_spatial_axis(parameterdict["z"], parameterdict["resolution"][2], parameterdict["spaceUnits"])
    return [x, y, z]


def build_sequence_blocks(parameterdict):
    """Return normalised sequence blocks from sequence parameters."""
    seqBlocks = []
    for rf, delays, reph, cf in zip(
        parameterdict["RF"], parameterdict["delays"], parameterdict["rephaseAreas"], parameterdict["CoherenceFilter"]
    ):
        seqBlocks.append(
            SequenceBlock(
                rf,
                delays,
                reph,
                cf,
                parameterdict["RFUnits"],
                parameterdict["GradUnits"],
                nucleus=parameterdict["Nucleus"],
            )
        )

    return seqBlocks


def choose_simulation_method(parameterdict, seqBlocks):
    """Return simulation method and gradient axis order."""
    if "method" in parameterdict.keys():
        methodPoss, axisOrder = gradOrderChecker(seqBlocks)
        method = parameterdict["method"]
        if (method == "1d") & ((methodPoss == "interleaved") | (methodPoss == "full")):
            print(f"1D method requested but not possible, running with {methodPoss}")
            method = methodPoss
        elif (method == "interleaved") & (methodPoss == "full"):
            print(f"interleaved method requested but not possible, running with {methodPoss}")
            method = methodPoss
    else:
        method, axisOrder = gradOrderChecker(seqBlocks)

    return method, axisOrder


def configure_spatial_resolution(parameterdict, method):
    """Return method and spatially-resolved output flag."""
    resolve = bool(parameterdict.get("spatiallyResolve", False))
    if resolve and method == "1d":
        method = "interleaved"

    return method, resolve


def prepare_sequence_blocks_for_method(method, seqBlocks, axisOrder):
    """Perform method-specific block preparation."""
    if method == "1d":
        return splitConcurrentRephase(seqBlocks)

    return seqBlocks, axisOrder


def create_simulator(spinsys, parameterdict):
    """Create simulator object from normalised inputs."""
    nucleus = parameterdict["Nucleus"]
    return sim.simulator(
        spinsys, B0=parameterdict["B0"], centralShift=parameterdict["centralShift"], gamma=GYRO_MAG_RATIO[nucleus] * 1e6
    )


def readout_sequence(obj, currentP, method, resolve, points, dwellTime, lw, spatialDims):
    """Generate final density matrix and FID output."""
    x, y, z = spatialDims
    if method == "1d":
        finalP = np.mean(currentP, axis=0)
        FID, ax = obj.readout(points, dwellTime, lw, p=finalP)
    elif method == "interleaved":
        if resolve:
            outShape = (points, x.size, y.size, z.size)
            currentP = currentP.reshape((-1,) + currentP.shape[-2:])
            FID = np.zeros(outShape, dtype=np.complex128)
            FID = FID.reshape((points, currentP.shape[0]))
            for iDx, p in enumerate(currentP):
                FID[:, iDx], ax = obj.readout(points, dwellTime, lw, p=p)
            FID = FID.reshape(outShape)
            finalP = currentP  # Just for returning
        else:
            finalP = np.mean(currentP.reshape((-1,) + currentP.shape[-2:]), axis=0)
            FID, ax = obj.readout(points, dwellTime, lw, p=finalP)
    elif method == "full":
        if resolve:
            outShape = (points, x.size, y.size, z.size)
            FID = np.zeros(outShape, dtype=np.complex128)
            FID = FID.reshape((points, currentP.shape[0]))
            for iDx, p in enumerate(currentP):
                FID[:, iDx], ax = obj.readout(points, dwellTime, lw, p=p)
            FID = FID.reshape(outShape)
            finalP = currentP  # Just for returning
        else:
            finalP = np.mean(currentP, axis=0)
            FID, ax = obj.readout(points, dwellTime, lw, p=finalP)

    return FID, ax, finalP


# Main sequence simulation function
def simseq(spinsys, parameterdict, verbose=True):
    """Sequence density matrix simulations.

    Args:
        spinsys (dict): Dict describing spin system to simulate.
        parameterdict (dict): Dict describing sequence and simulation parameters.

    Returns:
        FID (ndarray): Complex FID signal from readout.
        ax (dict): Dict of time, frequency and ppm axes.
        finalP (ndarray): Final state of the density matrix before readout.
    """
    # Check inputs
    parameterdict = normalise_sequence_params(parameterdict)
    spinsys = normalise_spin_system(spinsys)

    # Readout
    lw = parameterdict["Rx_LW"]
    dwellTime = 1 / parameterdict["Rx_SW"]
    points = parameterdict["Rx_Points"]
    receiverPhs = parameterdict["Rx_Phase"]

    # Spatial dimension
    spatialDims = build_spatial_dimensions(parameterdict)

    # Main sequence blocks comprising the RF, delays, rephase areas and coherence filters
    seqBlocks = build_sequence_blocks(parameterdict)

    # Determine simulation method
    # Options are 1D projection, out of order or full simulation.
    method, axisOrder = choose_simulation_method(parameterdict, seqBlocks)

    # Handle option for resolving spatial response
    method, resolve = configure_spatial_resolution(parameterdict, method)

    if verbose:
        print(f"Simulation running using mode {method}. Axis order = {axisOrder}.")

    # Deal with rephase gradients on different axes
    # Split the delay into more blocks, inserting zeros into pulses etc
    seqBlocks, axisOrder = prepare_sequence_blocks_for_method(method, seqBlocks, axisOrder)

    # Initialize the class object
    obj = create_simulator(spinsys, parameterdict)
    currentP = propagate_sequence(obj, seqBlocks, axisOrder, spatialDims, method)

    # Generate final result
    FID, ax, finalP = readout_sequence(obj, currentP, method, resolve, points, dwellTime, lw, spatialDims)

    FID = FID * np.exp(-1j * receiverPhs)  # *= doesn't work for some reason.

    return FID, ax, finalP
