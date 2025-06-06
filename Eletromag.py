import numpy as np
'''
    This package summarizes various electromagnetism and utility calculations.
'''

#Defining constants for clarity's sake
X, Y, Z, BX, BY, BZ, B = 0, 1, 2, 3, 4, 5, 6

# This is MU0/4pi in T.m/A
MU0_PRIME = 1e-7

MU0 = MU0_PRIME * 4 * np.pi


def circleArea(radius):
    circleArea = 2 * np.pi * radius
    return circleArea
def rectangleArea(width,length):
    rectangleArea = width * length
    return rectangleArea
def squareArea(side):
    squareArea = side**2
    return squareArea


def crossSectionalArea(fill_ratio:float=1,*,radius:float=None,width:float=None,length:float=None,side:float=None):
    if side != None:
        return squareArea(side) * fill_ratio
    elif width and length != None:
        return rectangleArea(width,length) * fill_ratio
    elif radius != None:
        return circleArea(radius) * fill_ratio
    else: 
        raise TypeError('Invalid paramers')
      
def createLine(Pa, Pb,*, max_seg_len:float = 1,n_points:int = None):
    """
        Calculates the list of coordinates (coil path) between two different points in 3D space.

        :Pa| list or np.array: Coordinates of the initial point.
        :Pb| list or np.array: Coordinates of the final point.
        :max_seg_len| float: (optional) Maximum length of each segment. Default is 1.
        :n_points| int: (optional) The number of points in the coil path. Default is None.


        Returns:
            list: A list of coordinates representing the segmented line (coil path).
    """
    # Transform all the objects in float if they are compactible, otherwise it will raise a value error
    try:
        Pa = [float(coordinate) for coordinate in Pa]
        Pb = [float(coordinate) for coordinate in Pb]
        max_seg_len = float(max_seg_len)

    except:
        raise ValueError("All elements in the input must be numbers, lists or arrays.")
    
    assert Pa != Pb, "The initial and final points must be different."
    assert max_seg_len > 0,'The maximun segment legth must be an positive number'
    assert type(n_points) == 'NoneType' or int
        
    projection = [b - a for a, b in zip(Pa, Pb)]
    length = sum(p**2 for p in projection) ** 0.5
    
    assert length >= max_seg_len,'The distance between the points must be equal or higher than the maximum segment length'
    
    max_seg_len = length/n_points if n_points != None else max_seg_len

    ratio = length / max_seg_len
    coilPath = [
        [projection[0] / ratio * i + Pa[0],
         projection[1] / ratio * i + Pa[1],
         projection[2] / ratio * i + Pa[2]]
        for i in range(int(ratio) + 1)
    ]

    if Pb not in coilPath:
        coilPath.append(Pb)

    return coilPath

def createArch(center, radius, start_angle, angle, max_seg_len, n_points=None, anticlockwise: bool = False):
    """
    Calculates the list of coordinates (coil path) in a specific arch in 3D space.

    :center| list or np.array: Coordinates of the initial point.
    :radius| float: Radius of the arch.
    :start_angle| float: Starting angle (radians).
    :angle| float: Total angle (radians) to sweep.
    :max_seg_len| float: (optional) Maximum length of each segment. Default is 1.
    :n_points| int: (optional) Number of points. If None, computed from max_seg_len.
    :anticlockwise| bool: If True, the arc is swept in the negative angular direction.

    Returns:
        list: A list of coordinates representing the segmented arch (coil path).
    """
    try:
        center = [float(coordinate) for coordinate in center]
        radius = float(radius)
        angle = float(angle)
        start_angle = float(start_angle)
        max_seg_len = float(max_seg_len)
    except:
        raise ValueError("All elements in the input must be numbers, lists or arrays.")

    assert radius > 0 and angle > 0, 'The radius and angle must be positive numbers'
    assert n_points is None or isinstance(n_points, int), 'n_points must be None or an integer'

    length = angle * radius
    n_points = int(length / max_seg_len) + 1 if n_points is None else n_points
    theta = angle / n_points

    if not anticlockwise:
        coilPath = [
            [
                center[0] + radius * np.cos(start_angle + i * theta),
                center[1] + radius * np.sin(start_angle + i * theta),
                center[2]
            ]
            for i in range(n_points + 1)
        ]
    else:
        coilPath = [
            [
                center[0] + radius * np.cos(start_angle - i * theta),
                center[1] + radius * np.sin(start_angle - i * theta),
                center[2]
            ]
            for i in range(n_points + 1)
        ]
    return coilPath

class Coil:
    
    
    def __init__(self, coilPath:np.ndarray, invertRAxis:bool=False, *, crossSectionalArea:float = 1 , resistivity:float=1.7e-8):
        '''
            initialize a instance from Coil class
        

            :param coilPath np.ndarray: the ordered array of points where the path goes trough.
            :param invertRAxis bol: (optional) used if the coilPath is in the format of [[X1,Y1,Z1],...] rather than [[X], [Y], [Z]].
            :param crossSectionalArea float: (optional) the area of a cross section of the wire, the default is unitary.
            :param resistivity float: (optional) the resistivity of the wire material, the default is copper.

        '''
        if isinstance(coilPath, str):
            self.coilPath = np.loadtxt(coilPath)
        else:
            self.coilPath = coilPath
        
        if not invertRAxis:
            self.coilPath = np.moveaxis(self.coilPath, 0, 1)
        if crossSectionalArea == None:
            raise ValueError('Insert a valid cross sectional area') 
        length = self.__calculateCoilLength()
        self._length = length
        self._resistivity = resistivity
        self._crossSectionalArea = crossSectionalArea
        self._resistance = self.__calculateCoilResistance()
        
    @property 
    def resistivity(self):
        return self._resistivity
    
    @resistivity.setter
    def resistivity (self, new_resistivity):
        if new_resistivity <= 0:
            raise ValueError("The value must be positive")
        self._resistivity = new_resistivity
        self._resistance = self.__calculateCoilResistance()
        
    @property 
    def length(self):
        return self._length
    
    @property 
    def resistance(self):
        return self._resistance
    
    @resistance.setter
    def resistance (self, new_resistance):
        if new_resistance <= 0:
            raise ValueError("The value must be positive")
        self._resistance = new_resistance
        
    def __calculateCoilLength(self):
        '''
            Calculates the length of a coil
        '''
        dl = self.coilPath[1:] - self.coilPath[:-1]

        return np.sum(np.linalg.norm( dl, axis=1 ))

    def __calculateCoilResistance(self):
        '''
            Calculates the resistance of a coil
        '''
  
        resistance = self._length * self._resistivity / self._crossSectionalArea
        return resistance
    
    def __BiotSavart1pDimensionless(self, r0:np.ndarray):
        '''
            Calculates the Biot-Savart integral for the point at r0.
            Neither the current nor the vacuum permissivity / 4pi are taken into account. This makes the process of aclculating multiple points faster.

            :param coilPath np.ndarray: the ordered array of points where the path goes trough.
            :param r0 np.ndarray(float): the point to check for the magnetic field.

            :returns float: the integral part of the Biot-Savart law for the coil path and the point r0.
        '''
        avgR = .5 * ( self.coilPath[:-1] + self.coilPath[1:])

        rPrime = r0-avgR

        rMod = np.linalg.norm(rPrime, axis=1)
        rMod = np.abs(rMod)[:, np.newaxis]

        dl = self.coilPath[1:] - self.coilPath[:-1]

        integ = np.sum( np.cross(dl, rPrime) / (rMod**3), axis=0 )
        return integ 

    def biotSavart1p(self, r0:np.ndarray, I:float):
        '''
            Calculates the magnetic field at point r0 by using Biot-Savart and assuming constant current.

            :param coilPath np.ndarray: the array of points where the path goes trough.
            :param r0 np.ndarray(float): the point to check for the magnetic field.
            :param I float: (optional) the current going through the coil.

            :returns np.ndarray: a list of the magnetic field components in r0 because of the current going through coilPath.
         '''
        MU0_PRIME = 1e-7

        outsideValue = I*MU0_PRIME
        
        # Calculates the integral.
        integ = self.__BiotSavart1pDimensionless(r0)
        return integ*outsideValue
    
    def biotSavart3d( self, pointsList:np.ndarray, I:float = 1, invertPAxis:bool=False ):
        '''
            Calculates the magnetic fields for an array of points in space by using Biot-Savart and assuming constant current.

            :param coilPath np.ndarray: the array of points where the path goes through, in meters.
            :param pointsList np.ndarray: the array of points to check for the magnetic field, in meters.
            :param I float: (optional) the current going through the coil, in Amperes.
            :param invertPAxis bool: (optional) whether the pointsList is in the format of [[x1,y1,z1],...] rather than [[X], [Y], [Z]].

            :returns np.ndarray: a list of coordinates and the respective magnetic field values for each point caused by the coilPath.
                The format is in the same shape as pointsList, beign either [[x1,y1,z1,bx1,by1,bz1],...] for [[X],[Y],[Z],[Bx],[By],[Bz]].
        '''

        MU0_PRIME = 1e-7

        if not invertPAxis:
            pointsList = np.moveaxis(pointsList, 0, 1)

        pListLen = np.shape(pointsList)[0]
        results = []
        for i in range(pListLen):
            singleResult = self.__BiotSavart1pDimensionless(pointsList[i])
            results.append(singleResult)

          # Multiplies the integrals by the outside factor. We do this all at the same time to save computational time.
        results = np.array(results) * I * MU0_PRIME

        # Returns the lists to the default orientation and concatenates the new data to each position.
        results = np.moveaxis(results, 0, 1)
        pointsList = np.moveaxis(pointsList, 0, 1)
        # returnal ends up looking like [[X], [Y], [Z], [Bx], [By], [Bz]]
        returnal = np.concatenate((pointsList, results))

        # Returns the new data in the same orientation that pointsList was given.
        if invertPAxis:
            returnal = np.moveaxis(returnal, 0, 1)
        return returnal
    
    def dissipationPotency (self, I):
        '''
        Calculates the dissipated potency of the coil

        :param I float: the current running through the coils.
        :param 

        '''
        dissipationPotency = self._resistance * I**2
        return dissipationPotency
    
 
def calculateMultipleCoils3D(coilList, pointsList:np.ndarray, I:float=1, invertRAxis:bool=False,
                         invertPAxis:bool=False, calculateB:bool=True, verbose:bool=False):

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
        Bfield = np.sqrt( returnal[BX]**2 + returnal[BY]**2 + returnal[BZ]**2 )
    returnal = np.vstack((returnal, Bfield))

    return returnal

def calculateMultipleCoilsLength(coilList):

    length = 0
    for coil in coilList:
        length += coil.length
    return length
def calculateMultipleCoilsResistance(coilList, invertRAxis:bool=False):
   
    length = calculateMultipleCoilsLength(coilList)
    resistance = length * coilList[0].resistivity / coilList[0]._crossSectionalArea
    return resistance