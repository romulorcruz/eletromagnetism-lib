import numpy as np
import matplotlib.pyplot as plt
BX = 0
BY = 1
BZ = 2
def circleArea(radius):
    return 2 * np.pi * radius
def rectangleArea(width,length):
    return width * length
def squareArea(side):
    return side**2

class CoilList:
    def __init__(self, list):
        for i in list:
            if type(i) != Coil:
                raise TypeError('The objects of CoilList must be Coil objects')
            
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
def createArch(center, radius, angle,initial_angle ,max_seg_len,n_points = None, anticlockwise:bool = False):
    """
        Calculates the list of coordinates (coil path) in a specific arch in 3D space.

        :center| list or np.array: Coordinates of the initial point.
        :radius| float: Radius of the arch
        :angle| float: Angle of the arch.
        :max_seg_len| float: (optional) Maximum length of each segment. Default is 1.
        :n_points| int: (optional) The number of points in the coil path. Default is None.

        Returns:
            list: A list of coordinates representing the segmented arch (coil path).
    """
    # Transform all the objects in the center list to float if they are compactible, otherwise it will raise a value error
    try:
        center = [float(coordinate) for coordinate in center]
        radius = float(radius)
        angle = float(angle)
        max_seg_len = float(max_seg_len)
    except:
        raise ValueError("All elements in the input must be numbers, lists or arrays.")
        
    assert radius and angle > 0,'The radius and angle must be positive numbers'
    assert type(n_points) == 'NoneType' or int
    
    length = angle * radius
    n_points = int(length / max_seg_len) + 1 if n_points == None else n_points
    theta = angle / n_points


    if anticlockwise == False:
        coilPath_ = [[
            center[0] + radius * np.cos((-i) * theta + initial_angle),
            center[1] + radius * np.sin((-i) * theta + initial_angle),
            center[2]
            ]
            for i in range(n_points + 1)
        ]

    else:
        coilPath_ = [
        [
            center[0] + radius * np.cos((-i) * theta + initial_angle),
            center[1] + radius * np.sin((-i) * theta + initial_angle),
            center[2]
        ]
        for i in range(n_points + 1)
    ]
    coilPath = []
    [coilPath.append(point) for point in coilPath_ if point not in coilPath]
    return coilPath

    
class Coil:
    
    
    def __init__(self, coilPath:np.ndarray, invertRAxis:bool=False, *, crossSectionalArea:float, resistivity:float=1.7e-8):
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
            self.coilPath = np.array(coilPath)
        
        if invertRAxis:
            self.coilPath = np.moveaxis(self.coilPath, 0, 1)
        if crossSectionalArea == None:
            raise ValueError('Insert a valid cross sectional area') 
        length = self.__calculateCoilLength()
        self._length = length
        self._resistivity = resistivity
        self._crossSectionalArea = crossSectionalArea
        self.resistance = self.__calculateCoilResistance()
        
    @property 
    def resistivity(self):
        return self._resistivity
    
    @resistivity.setter
    def resistivity (self, new_resistivity):
        if new_resistivity <= 0:
            raise ValueError("The value must be positive")
        self._resistivity = new_resistivity
        self.resistance = self.__calculateCoilResistance()
        
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
    
    def biotSavart3d(self, pointsList:np.ndarray, I:float = 1, invertPAxis:bool=False ):
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
            singleResult = self.__BiotSavart1pDimensionless(pointsList[i, :])
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
        dissipationPotency = self.resistance * I**2
        return dissipationPotency
    
    def plotCoil(self):
        '''Tá funcionando só em dois eixos'''
        x = [point[0] for point in self.coilPath]
        y = [point[1] for point in self.coilPath]
        return plt.plot(x,y)
    
 
def calculateMultipleCoils3D(coilList, pointsList:np.ndarray, I:float=1, invertRAxis:bool=False,
                         invertPAxis:bool=False, calculateB:bool=True, verbose:bool=False):

    nCoils = len(coilList)

    if verbose:
        print(f"\nCalculating coil 1 out of {nCoils:d}")

    # Calculates the first coil separately for convenience.
    returnal = Coil.biotSavart3d(pointsList, I, invertRAxis, invertPAxis)
    for i in range(1, nCoils):
        if verbose:
            print(f"\nCalculating coil {i+1:d} out of {nCoils:d}")
        returnal[BX:BZ+1] += Coil.biotSavart3d(pointsList, I, invertRAxis, invertPAxis)[BX:BZ+1]

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
    resistance = length * coilList[0].resistivity / coilList[0].crossSectionalArea
    return resistance

def crossSectionalArea(fill_ratio:float=1,*,radius:float=None,width:float=None,length:float=None,side:float=None):
    if side != None:
        return squareArea(side) * fill_ratio
    elif width and length != None:
        return rectangleArea(width,length) * fill_ratio
    elif radius != None:
        return circleArea(radius) * fill_ratio
    else: 
        raise TypeError('Invalid paramers')