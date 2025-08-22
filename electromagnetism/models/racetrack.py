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