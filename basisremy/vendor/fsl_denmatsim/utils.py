import json
import os
from . import simulator as sim
from . import simseq
from .constants import GYRO_MAG_RATIO
import numpy as np
# File for odds and ends functions


# Read in the spin systems defined in the spinSystems json file distributed with the package
def readBuiltInSpins():
    spinSysFile = os.path.join(os.path.dirname(__file__), "spinSystems.json")
    with open(spinSysFile, "r") as jsonFile:
        jsonString = jsonFile.read()
        spinsys = json.loads(jsonString)

    return spinsys


# Reconstruct a FID from a saved density matrix (probably from a results file)
# at a given time resolution and number of points, avoids interpolation.
def FIDFromDensityMat(p, spinsys, B0, points, dwellTime, lw, offset=0, receiverPhs=0, nucleus="1H"):
    if isinstance(spinsys, list):
        combinedScaledFID = []
        for s, pp in zip(spinsys, p):
            s = simseq.grumblespins(s)
            obj = sim.simulator(s, B0, gamma=GYRO_MAG_RATIO[nucleus] * 1e6, centralShift=offset)
            FID, _ = obj.readout(points, dwellTime, lw, p=pp)
            FID = FID * np.exp(-1j * receiverPhs)
            combinedScaledFID.append(FID * s["scaleFactor"])
        combinedScaledFID = np.sum(combinedScaledFID, axis=0)
    else:
        spinsys = simseq.grumblespins(spinsys)
        obj = sim.simulator(spinsys, B0, gamma=GYRO_MAG_RATIO[nucleus] * 1e6, centralShift=offset)
        FID, _ = obj.readout(points, dwellTime, lw, p=p)
        FID = FID * np.exp(-1j * receiverPhs)
        combinedScaledFID = FID * spinsys["scaleFactor"]

    return combinedScaledFID
