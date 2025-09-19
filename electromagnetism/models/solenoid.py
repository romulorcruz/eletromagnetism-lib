from numpy import array, ndarray, loadtxt, moveaxis, newaxis, cross, concatenate, shape
import numpy as np
from numpy.linalg import norm
from electromagnetism.mathematics.geometry import helicoid

class Solenoid(Coil):
    def __init__(self, n_turns: int, initial_point, final_point, radius: float, max_seg_len: float,
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
        self.coilPath = np.array(helicoid(self.n_turns, self.z_initial_point, self.z_final_point, self.radius, self.max_seg_len))
 
 
        super().__init__(self.coilPath, invertRAxis=invertRAxis, crossSectionalArea=crossSectionalArea, resistivity=resistivity)

    