"""Electromagnetic Coil Module.

This module defines the `Coil` class, encapsulating the properties of a 
single electromagnetic coil. It provides functionalities for calculating
intrinsic properties like length and resistance, and computing the 
magnetic field it generates at various space in space using numerical 
methods (e.g., Biot-Savart Law).
"""

from numpy import array, ndarray, loadtxt, moveaxis, newaxis, cross, concatenate, shape
from numpy.linalg import norm
import numpy as np
from ..mathematics.constants import MU0_PRIME
from electromagnetism.mathematics.geometry import helicoid
from electromagnetism.mathematics.geometry import racetrack3d
import plotly.express as px
import plotly.graph_objects as go
from scipy.integrate import simpson
from plotly.io import show
class Coil:
    def __init__(self, coilPath:ndarray, invertRAxis:bool=False,
                 *, crossSectionalArea:float = 1 , resistivity:float=1.7e-8):
        if isinstance(coilPath, str):
            self.coilPath = loadtxt(coilPath)
        else:
            self.coilPath = coilPath

        if invertRAxis:
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
        '''
        if new_resistance <= 0:
            raise ValueError("The value must be positive")
        self._resistance = new_resistance
    
    @property
    def crossSectionalArea(self):
        '''
            Returns the Cross Sectional Area of the wire.
        '''
        return self._crossSectionalArea
 
    @crossSectionalArea.setter
    def crossSectionalArea(self, new_Cross_sec_area):
        '''
            Sets the cross Sectional Area value of the wire .
 
            :param new_Cross_sec_area float: A new value for the cross sectional area that must be positive.
 
            :raises ValueError: if new_Cross_sec_area is not positive.
        '''
        if new_Cross_sec_area <= 0:
            raise ValueError("The value must be positive")
        self._crossSectionalArea = new_Cross_sec_area
        self._resistance = self.__calculateCoilResistance()

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
            This makes the process of aclculating multiple space faster.

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
            avgR = 0.5 * (self.coilPath[:-1] + self.coilPath[1:])
            dl = self.coilPath[1:] - self.coilPath[:-1]
            rPrime = r0 - avgR

            integrand = np.cross(dl, rPrime) / norm(rPrime, axis=1)[:, np.newaxis]**3
            integ = simpson(integrand, axis=0)
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

        integ = self.__BiotSavart1pDimensionless(r0)
        return integ*outsideValue

    def biotSavart3d( self, spaceList:ndarray, I:float = 1, invertPAxis:bool=False, * ,integration_method: str = 'Simpson' ):
        '''
            Calculates the magnetic fields for an array of space in space by using Biot-Savart
            and assuming constant current.

            :param spaceList numpy.ndarray: the array of space to check for the magnetic field.
            :param I float: (optional) the current going through the coil, in Amperes.
            :param invertPAxis bool: (optional) whether the spaceList is in the format of 
            [[x1,y1,z1],...] rather than [[X], [Y], [Z]].

            :returns numpy.ndarray: a list of coordinates and the respective magnetic field values 
            for each point caused by the coilPath.
                The format is in the same shape as spaceList, beign either 
                [[x1,y1,z1,bx1,by1,bz1],...] for [[X],[Y],[Z],[Bx],[By],[Bz]].
        '''

        if invertPAxis:
            spaceList = moveaxis(spaceList, 0, 1)

        pListLen = shape(spaceList)[0]
        results = []
        for i in range(pListLen):
            singleResult = self.__BiotSavart1pDimensionless(spaceList[i], integration_method=integration_method)
            results.append(singleResult)

        # Multiplies the integrals by the outside factor.
        results = array(results) * I * MU0_PRIME
        # returnal ends up looking like [[X], [Y], [Z], [Bx], [By], [Bz]]
        results = moveaxis(results, 0, 1)
        spaceList = moveaxis(spaceList, 0, 1)
        returnal = concatenate((spaceList, results))

        # Returns the new data in the same orientation that spaceList was given.
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


    def cloud(self, padding: float = 0.0,n:int = 10,i = 1, integration_method:float='Simpson', plane_axis=None, plane_value='mid', plane_thickness=0.0, show=False):
            '''
                Generate a cloud of points by a mesh grid and calculates the magnetic field at each point.
                :param padding float: space between the maximum and minimum point of the object and the start and end of meshgrid to every axis
                :param n int: number of points in every axis before the generation of mesh grid
                :param i float: the value of current in coil
                :param integration_method str: integration method for Biot-Savart law
                
                :returns fig,b numpy.array: returns the figure if plot is true and the vector b with the coordinates of the points and the magnetic field
                at x,y and z and the total field by the norm of x,y and z components.
            '''
            x = np.linspace(self.coilPath[:,0].min()-padding, self.coilPath[:,0].max()+padding,n)
            y = np.linspace(self.coilPath[:,1].min()-padding, self.coilPath[:,1].max()+padding,n)
            z = np.linspace(self.coilPath[:,2].min()-padding, self.coilPath[:,2].max()+padding,n)
    
            xx,yy, zz = np.meshgrid(x,y,z)
            space = np.zeros((x.shape[0] * y.shape[0] * z.shape[0], 3))
            space[:, 0] = xx.flatten()
            space[:, 1] = yy.flatten()
            space[:, 2] = zz.flatten()
    
            b = self.biotSavart3d(space,integration_method=integration_method, I= i)
            b_t = np.linalg.norm((b[3],b[4],b[5]), axis=0)
            b = np.concatenate((b,[b_t]))
    
            fig = go.Figure()
            fig.add_trace(go.Scatter3d(
            x=space[:, 0],
            y=space[:, 1],
            z=space[:, 2],
            mode="markers",
            marker=dict(
                size=5,                    
                color=b_t,
                colorscale="Plasma",
                opacity=0.5,
                colorbar=dict(title="|B|"),
            ),
            name="Field"
            ))
    
            if plane_axis is not None:
                axis_map = {'x': 0, 'y': 1, 'z': 2}
                ax = axis_map.get(plane_axis.lower())
    
                if ax is not None:
                    coord = space[:, ax]
    
                    if plane_value == 'mid':
                        val = 0.5 * (coord.min() + coord.max())
                    else:
                        val = float(plane_value)
    
                    if plane_thickness <= 0.0:
                        if plane_axis.lower() == 'x':
                            dz = x[1] - x[0]
                        elif plane_axis.lower() == 'y':
                            dz = y[1] - y[0]
                        else:
                            dz = z[1] - z[0]
                        plane_thickness_eff = abs(dz)
                    else:
                        plane_thickness_eff = float(plane_thickness)
    
                    mask = np.abs(coord - val) <= plane_thickness_eff / 2.0
    
                    if mask.any():
                        fig.add_trace(go.Scatter3d(
                            x=space[mask, 0],
                            y=space[mask, 1],
                            z=space[mask, 2],
                            mode="markers",
                            marker=dict(
                                size=4,              
                                color=b_t[mask],
                                colorscale="Plasma",
                                opacity=0.9
                            ),
                            name=f"Field on {plane_axis} = {val:.3g}"
                        ))
            fig.update_layout(
                scene=dict(
                    xaxis_title="x",
                    yaxis_title="y",
                    zaxis_title="z",
                    aspectmode="data",
                ),
                margin=dict(l=0, r=0, b=0, t=30),
                title="Magnetic field cloud"
            )
    
            if show:
                fig.show()
        
            return fig, b
    
    def plot(self):
        pts = self.coilPath  # (N, 3)
        x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
 
        copper_metallic = [
            [0.00, "#2b1306"],  # marrom bem escuro, sombra
            [0.20, "#5a2610"],  # marrom avermelhado
            [0.40, "#8c3f1c"],  # cobre escuro
            [0.60, "#c7632a"],  # cobre médio
            [0.80, "#e9924a"],  # cobre mais claro/brilho
            [1.00, "#ffe0b3"],  # highlight quase dourado
        ]
 
        fig = go.Figure(
            data=[
                go.Scatter3d(
                    x=x,
                    y=y,
                    z=z,
                    mode="lines",
                    line=dict(
                        width=6,
                        color=z,                 # usa z como “valor” de cor
                        colorscale=copper_metallic,     # ou "Viridis", "Turbo", etc.
                        cmin=float(z.min()),
                        cmax=float(z.max())
                    )
                )
            ]
        )
 
        fig.update_layout(
            scene=dict(aspectmode="data"),
            coloraxis_colorbar=dict(title="z")
        )
 
        if show:
            fig.show()
        return fig

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
        self.coilPath = np.array(helicoid(self.n_turns, [0,0,self.z_initial_point],
            self.z_initial_point + self.z_final_point, self.radius, self.max_seg_len))
 
 
        super().__init__(self.coilPath, invertRAxis=invertRAxis, crossSectionalArea=crossSectionalArea, resistivity=resistivity)

class Racetrack(Coil):
    
    def __init__(self, center, inwidth: float, inlength: float, max_seg_len: float, int_radius: float, thickness: float,
                 height: float,
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
        self.height = height
        self.coilPath = np.array(racetrack3d(self.center, self.inwidth, self.inlength, self.max_seg_len,
                                             self.int_radius, self.thickness, self.height))
        
        super().__init__(self.coilPath, invertRAxis=invertRAxis, crossSectionalArea=crossSectionalArea, resistivity=resistivity)
    