"""System Calculations Module.

This module contains functions for calculations on multi-coil systems, like length
and resistance of serie associated coils and magnetic field for only multiple coils 
systems with same current using numerical methods (e.g., Biot-Savart Law).
"""
from numpy import sqrt, vstack, ndarray
from .mathematics.constants import BX,BY,BZ
def calculateMultipleCoilsLength(coilList):
    '''
        Calculates the sum of lengths from multiple coils.

        :param coilList list|numpy.ndarray: a list of ordered coil paths.
        :param invertRAxis bool: (optional) whether each coil path is in the format [[x1,y1,z1],...] rather than [[X], [Y], [Z]].

        :returns float: the sum of the lengths from the coils in coilList, in the same unit as the inumpyut coordinates.
    '''
    length = 0
    for coil in coilList:
        length += coil.length
    return length


def calculateMultipleCoilsResistance(coilList, invertRAxis:bool=False):
    '''
        Calculates the resistance of multiple coils as if they are in series and have uniform resistivity and cross sectional area.
        All the units should be consistent, so if the coordinates of the path are in meters, the cross sectional area should be in squared meters.

        :param coilList list|numpy.ndarray: a list or array of ordered coil paths.
        :param crossSectionalArea float: (optional) the area of a cross section of the wire, the default is unitary, in meters.
        :param resistivity float: (optional) the resistivity of the wire material, the default is copper.
        :param invertRAxis bool: (optional) whether each coil path is in the format [[x1,y1,z1],...] rather than [[X], [Y], [Z]].

        :returns float: the resistance of the coils from the list in series.
    '''

    length = calculateMultipleCoilsLength(coilList)
    resistance = length * coilList[0].resistivity / coilList[0]._crossSectionalArea
    return resistance

def calculateMultipleCoils3D(coilList, pointsList: ndarray, I:float=1, invertRAxis:bool=False,
                         invertPAxis:bool=False, calculateB:bool=True, verbose:bool=False):
    '''
     Calculates the Biot-Savart law for multiple coil paths with the same current.

        :param coilList list|numpy.ndarray: a list of paths that each coil takes, in meters.
        :param pointsList numpy.ndarray: an array of points where the magnetic field will be calculated at, in meters.
        :param I float: (optional) the current going through the coil, in amperes.
        :param invertRAxis bool: (optional) whether the coilPath is in the format of [[X1,Y1,Z1],...] rather than [[X], [Y], [Z]].
        :param invertPAxis bool: (optional) whether the pointsList is in the format of [[X1,Y1,Z1],...] rather than [[X], [Y], [Z]].
        :param calculateB bool: (optional) whether or not the calculation of the magnetic field modulus should take place.
        :param verbose bool: (optional) makes the function print the progress of the calculations, usefull for long lists of points.

        :returns numpy.ndarray: a list of coordinates and the respective magnetic field values for each point caused by the coilList.
            The format is in the same shape as pointsList, beign either [[x1,y1,z1,bx1,by1,bz1*,b*],...] for [[X],[Y],[Z],[Bx],[By],[Bz]*,[B]*].
    '''
    nCoils = len(coilList)

    if verbose:
        print(f"\nCalculating coil 1 out of {nCoils:d}")

    # Calculates the first coil separately for convenience.
    returnal = coilList[0].biotSavart3d(pointsList, invertPAxis = invertPAxis)
    for i in range(1, nCoils):
        if verbose:
            print(f"\nCalculating coil {i+1:d} out of {nCoils:d}")
    returnal [BX:BZ+1] += coilList[i].biotSavart3d( pointsList, invertPAxis = invertPAxis)[BX:BZ+1]

    # Calculates the modulus of the magnetic field for the points
    if calculateB:
        Bfield = sqrt( returnal[BX]**2 + returnal[BY]**2 + returnal[BZ]**2 )
    returnal = vstack((returnal, Bfield))

    return returnal
