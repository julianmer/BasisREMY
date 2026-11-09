import numpy as np
import scipy.linalg as lin


### Define some maths functions ###
# Kronecker product of list of matrices
# Applied in order ((0x1)x2)x3...
def kronInOrder(matList):
    if len(matList) == 1:
        return matList[0]

    out = np.kron(matList[0], matList[1])
    for m in matList[2:]:
        out = np.kron(out, m)
    return out


# See page 138 of Levitt Spin Dynamics.
# exp(A) is not equal to the exponential of each element of A!
def diagOpExp(op):
    # calculate the exponential of a diagonal operator
    return np.diag(np.exp(np.diagonal(op)))


def diagonaliseMat(op):
    w, v = np.linalg.eig(op)
    D = np.diag(w)
    X = v
    if np.array_equal(op, op.conj().T) & (np.linalg.norm(v, axis=0) == 1).all():
        Xp = X.conj().T
    else:
        Xp = np.linalg.inv(X)
    return X, D, Xp


def opExp(op):
    # if np.count_nonzero(op - np.diag(np.diagonal(op)))==0:
    #     #print('is diagonal')
    #     out = diagOpExp(op)
    # else:
    #     #print('is not diagonal')
    #     # Diagonalise the matrix first
    #     X,D,Xp = diagonaliseMat(op)
    #     out = X@diagOpExp(D)@Xp

    return lin.expm(op)


### Define basic angular momentum operators ###
def singleIz():
    return 0.5 * np.array([[1, 0], [0, -1]])


def singleIx():
    return 0.5 * np.array([[0, 1], [1, 0]])


def singleIy():
    return (1 / (2j)) * np.array([[0, 1], [-1, 0]])


def createIz(nspins):
    # Take
    AllIz = []
    for iDx in range(0, nspins):
        currentOps = []
        for jDx in range(0, nspins):
            if iDx == jDx:
                currentOps.append(singleIz())
            else:
                currentOps.append(np.eye(2))
        AllIz.append(kronInOrder(currentOps))

    return AllIz


def createIx(nspins):
    # Take
    AllIx = []
    for iDx in range(0, nspins):
        currentOps = []
        for jDx in range(0, nspins):
            if iDx == jDx:
                currentOps.append(singleIx())
            else:
                currentOps.append(np.eye(2))
        AllIx.append(kronInOrder(currentOps))

    return AllIx


def createIy(nspins):
    # Take
    AllIy = []
    for iDx in range(0, nspins):
        currentOps = []
        for jDx in range(0, nspins):
            if iDx == jDx:
                currentOps.append(singleIy())
            else:
                currentOps.append(np.eye(2))
        AllIy.append(kronInOrder(currentOps))

    return AllIy
