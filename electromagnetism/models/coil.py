"""Electromagnetic Coil Module.

This module defines the `Coil` class, encapsulating the properties of a 
single electromagnetic coil. It provides functionalities for calculating
intrinsic properties like length and resistance, and computing the 
magnetic field it generates at various points in space using numerical 
methods (e.g., Biot-Savart Law).
"""

from numpy import array, ndarray, loadtxt, moveaxis, newaxis, cross, concatenate, shape
from numpy.linalg import norm
import numpy as np
from ..mathematics.constants import MU0_PRIME
from electromagnetism.mathematics.geometry import helicoid
from electromagnetism.mathematics.geometry import racetrack3d
import plotly.express as px

class Coil:
    """Coil class.
    This class encapsulate the properties of a single electromagnetic coil. It provides 
    functionalities for calculating intrinsic properties like length and resistance, and 
    computing the magnetic field it generates at various points in space using numerical 
    methods (e.g., Biot-Savart Law)."""
    def __init__(self, coilPath:ndarray, invertRAxis:bool=False,
                 *, crossSectionalArea:float = 1 , resistivity:float=1.7e-8):
        '''
            initialize a instance from Coil class
        

            :param coilPath numpy.ndarray: the ordered array of points where the path goes trough.
            :param invertRAxis bol: (optional) used if the coilPath is in the format of 
            [[X1,Y1,Z1],...] rather than [[X], [Y], [Z]].
            :param crossSectionalArea float: (optional) the area of a cross section of the wire in 
            squared meters, the default is unitary.
            :param resistivity float: (optional) the resistivity of the wire material in Ohm meter, 
            the default is copper.

        '''
        if isinstance(coilPath, str):
            self.coilPath = loadtxt(coilPath)
        else:
            self.coilPath = coilPath

        if not invertRAxis:
            self.coilPath = moveaxis(self.coilPath, 0, 1)
        if crossSectionalArea is None:
            raise ValueError('Insert a valid cross sectional area')
        length = self.__calculateCoilLength()
        self._length = length
        self._resistivity = resistivity
        self._crossSectionalArea = crossSectionalArea
        self._resistance = self.__calculateCoilResistance()

    @property
    def resistivity(self):
        '''
            Returns the resistivity of the wire material.
        '''
        return self._resistivity

    @resistivity.setter
    def resistivity (self, new_resistivity):
        '''
            Sets the resistivity value of the wire material.

            :param new_resistivity float: A new value for the resistivity that must be positive.

            :raises ValueError: if new_resistance is not positive.
        '''
        if new_resistivity <= 0:
            raise ValueError("The value must be positive")
        self._resistivity = new_resistivity
        self._resistance = self.__calculateCoilResistance()

    @property
    def length(self):
        '''
            returns the length of the coil.
        '''
        return self._length

    @property
    def resistance(self):
        '''
            Returns the resistance of the coil.
        '''
        return self._resistance

    @resistance.setter
    def resistance (self, new_resistance):
        '''
            Manually sets the resistance value of the coil.

            :param new_resistance float: A new value for the resistance that must be positive.

            :raises ValueError: if new_resistance is not positive.
        '''
        if new_resistance <= 0:
            raise ValueError("The value must be positive")
        self._resistance = new_resistance

    def __calculateCoilLength(self):
        '''
            Calculates the length of a coil in meters
        '''
        dl = self.coilPath[1:] - self.coilPath[:-1]
        return np.sum(norm(dl, axis=1 ))


    def __calculateCoilResistance(self):
        '''
            Calculates the resistance of a coil in Ohms
        '''

        resistance = self._length * self._resistivity / self._crossSectionalArea
        return resistance

    def __BiotSavart1pDimensionless(self, r0: ndarray, integration_method: str = 'Simpson'):
        '''
            Calculates the Biot-Savart integral for the point at r0.
            Neither the current nor the vacuum permissivity / 4pi are taken into account. 
            This makes the process of aclculating multiple points faster.

            :param r0 numpy.ndarray(float): the point to check for the magnetic field.

            :returns float: the integral part of the Biot-Savart law for the coil path 
            and the point r0.
        '''
        if integration_method == 'Riemann':
            avgR = 0.5 * (self.coilPath[:-1] + self.coilPath[1:])
            dl = self.coilPath[1:] - self.coilPath[:-1]
                        
            rPrime = r0 - avgR
            rMod = norm(rPrime, axis=1)

            integrand = np.cross(dl, rPrime) / rMod[:, np.newaxis]**3
            integ = np.sum(integrand, axis=0)
            return integ

        elif integration_method == 'Simpson':

            points = self.coilPath
            rPrime = r0 - self.coilPath
            rMod = norm(rPrime, axis=1)
            epsilon = 1e-12
            rMod = np.maximum(rMod, epsilon)
            dl = np.gradient(self.coilPath, axis=0)
            integrand = np.cross(dl, rPrime) / rMod[:, np.newaxis]**3

            h = norm(self.coilPath[1] - self.coilPath[0])

            weights = 2 * np.ones(len(self.coilPath))
            weights[1::2] = 4
            weights[0], weights[-1] = 1, 1

            integ = (h/3) * np.sum(integrand * weights[:, np.newaxis], axis=0)
            return integ



    def biotSavart1p(self, r0:ndarray, I:float):
        '''
            Calculates the magnetic field at point r0 by using Biot-Savart
            and assuming constant current.

            :param r0 numpy.ndarray(float): the point to check for the magnetic field.
            :param I float: (optional) the current going through the coil in meters.

            :returns numpy.ndarray: a list of the magnetic field components in r0 
            because of the current going through coilPath.
         '''
        outsideValue = I*MU0_PRIME

        # Calculates the integral.
        integ = self.__BiotSavart1pDimensionless(r0)
        return integ*outsideValue

    def biotSavart3d( self, pointsList:ndarray, I:float = 1, invertPAxis:bool=False ):
        '''
            Calculates the magnetic fields for an array of points in space by using Biot-Savart
            and assuming constant current.

            :param pointsList numpy.ndarray: the array of points to check for the magnetic field.
            :param I float: (optional) the current going through the coil, in Amperes.
            :param invertPAxis bool: (optional) whether the pointsList is in the format of 
            [[x1,y1,z1],...] rather than [[X], [Y], [Z]].

            :returns numpy.ndarray: a list of coordinates and the respective magnetic field values 
            for each point caused by the coilPath.
                The format is in the same shape as pointsList, beign either 
                [[x1,y1,z1,bx1,by1,bz1],...] for [[X],[Y],[Z],[Bx],[By],[Bz]].
        '''

        if invertPAxis:
            pointsList = moveaxis(pointsList, 0, 1)

        pListLen = shape(pointsList)[0]
        results = []
        for i in range(pListLen):
            singleResult = self.__BiotSavart1pDimensionless(pointsList[i])
            results.append(singleResult)

          # Multiplies the integrals by the outside factor.
        results = array(results) * I * MU0_PRIME

        # Returns the lists to the default orientation and concatenates the
        # new data to each position.
        results = moveaxis(results, 0, 1)
        pointsList = moveaxis(pointsList, 0, 1)
        # returnal ends up looking like [[X], [Y], [Z], [Bx], [By], [Bz]]
        returnal = concatenate((pointsList, results))

        # Returns the new data in the same orientation that pointsList was given.
        if invertPAxis:
            returnal = moveaxis(returnal, 0, 1)
        return returnal

    def dissipationPotency (self, I):
        '''
            Calculates the dissipated potency of the coil

            :param I float: the current running through the coils in meters.
            
            Returns the dissipated potency of the coil in Watts.

        '''
        dissipationPotency = self._resistance * I**2
        return dissipationPotency


    def cloud(self, n):
        min = np.min(self.coilPath)
        max = np.max(self.coilPath)

        xs = np.linspace(min, max, n)
        ys = np.linspace(min, max, n)
        zs = np.linspace(min, max, n)

        space = np.moveaxis(np.array([xs,ys,zs]), 1,0)
        return space
    
    
    def plot(self):
        fig = px.line_3d(self.coilPath, x = self.coilPath[:,0], y = self.coilPath[:,1], z = self.coilPath[:,2])
        fig.show()
    
    
class Solenoid(Coil):
    def __init__(self, n_turns: int, z_initial_point: float, z_final_point: float, radius: float, max_seg_len: float,
                 *, crossSectionalArea: float = 1.0, resistivity: float = 1.7e-8,invertRAxis: bool = False):

        if n_turns <= 0:
            raise ValueError("Invalid parameters.")
        if radius <= 0:
            raise ValueError("Invalid parameters.")
        if max_seg_len <= 0:
            raise ValueError("Invalid parameters.")
        self.n_turns = n_turns
        self.z_initial_point = z_initial_point
        self.z_final_point = z_final_point
        self.radius = radius
        self.max_seg_len = max_seg_len
        self._resistivity = resistivity
        self._crossSectionalArea = crossSectionalArea
        self.coilPath = np.array(helicoid(self.n_turns, [0,0,self.z_initial_point], [0,0,self.z_final_point], self.radius, self.max_seg_len))
 
 
        super().__init__(self.coilPath, invertRAxis=invertRAxis, crossSectionalArea=crossSectionalArea, resistivity=resistivity)

class Racetrack(Coil):
    
    def __init__(self, center, inwidth: float, inlength: float, max_seg_len: float, int_radius: float, thickness: float,
                 *, crossSectionalArea: float = 1.0, resistivity: float = 1.7e-8,invertRAxis: bool = False):
        if int_radius <= 0:
            raise ValueError("Invalid parameters.")
        if max_seg_len <= 0:
            raise ValueError("Invalid parameters.")
        
        self.center = center
        self.inwidth = inwidth
        self.inlength = inlength
        self.max_seg_len = max_seg_len
        self.int_radius = int_radius
        self.thickness = thickness
        self._resistivity = resistivity
        self._crossSectionalArea = crossSectionalArea
        self.coilPath = np.array(racetrack3d(self.center, self.inwidth, self.inlength, self.max_seg_len,
                                             self.int_radius, self.thickness))
        
        super().__init__(self.coilPath, invertRAxis=invertRAxis, crossSectionalArea=crossSectionalArea, resistivity=resistivity)
    