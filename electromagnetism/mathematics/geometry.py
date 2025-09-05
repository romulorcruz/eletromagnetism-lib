"""Geometry Module.

This module contains all geometric tools that can be useful for
resistance calculation and coil path generations.

Contains:
Area calculations: circle, rectangle, square.
Geometrics figure coordinates generators: arch, line, racetrack.
"""
from numpy import sin, cos, linspace, array, moveaxis
from numpy.linalg import norm
from .constants import pi

origin = [0,0,0]

def circleArea(radius: float):
    '''
        Calculates the area of a circle based on its radius.
        
        :param radius float: the radius of the circle in meters.
        
        :returns: float: The area of the circle squared meters.

     '''
    circleArea = pi * radius **2
    return circleArea


def rectangleArea(width: float,length: float):
    '''
        Calculates the area of a rectangle based on its width and legnth.
        
        :param width float: the width of the rectangle in meters.
        :param length float: the length of the rectangle in meters.

        :returns: float: The area of the rectangle in squared meters.

     '''
    rectangleArea = width * length
    return rectangleArea


def squareArea(side: float):
    '''
        Calculates the area of a square using teh length of its side.
        
        :param side float: the length of the square's side in meters.
        
        :returns: float: The area of the square in squared meters.

     '''    
    squareArea = side**2
    return squareArea


def crossSectionalArea(fill_ratio:float=1,*,radius:float=None,width:float=None,length:float=None,side:float=None):
    '''
        Calculates the cross-sectional area of a conductor based on its shape and a fill ratio.

        Only one shape parameters (radius, width+length or side) must be provided. 
        The function multiplies the area obtained by the fill_ratio.

        :param fill_ratio float: (optional) Ratio representing the fraction of the shape area occupied by conductor (default is 1).
        :param radius: float, optional
            The radius of a circular cross-section.
        :param width: float, optional
            The width of a rectangular cross-section.
        :param length: float, optional
            The length of a rectangular cross-section.
        :param side: float, optional
            The side length of a square cross-section.

        :returns: float
            The calculated cross-sectional area adjusted by the fill_ratio.

        :raises TypeError: if the parameters do not define a valid shape.
    '''
    if side is not None:
        return squareArea(side) * fill_ratio
    elif width and length is not None:
        return rectangleArea(width,length) * fill_ratio
    elif radius is not None:
        return circleArea(radius) * fill_ratio
    else: 
        raise TypeError('Invalid paramers')
      

def createLine(Pa, Pb,*, max_seg_len:float = 1,n_points:int = None):
    """
        Calculates the list of coordinates (coil path) between two different points in 3D space.

        :Pa| list or numpy.array: Coordinates of the initial point.
        :Pb| list or numpy.array: Coordinates of the final point.
        :max_seg_len| float: (optional) Maximum length of each segment. Default is 1.
        :n_points| int: (optional) The number of points in the coil path. Default is None.


        Returns:
            list: A list of coordinates representing the segmented line (coil path).
    """
    # Transform all the objects in float if they are compactible, otherwise it will raise a value error
    try:
        Pa = array([float(coordinate) for coordinate in Pa])
        Pb = array([float(coordinate) for coordinate in Pb])
        max_seg_len = float(max_seg_len)

    except:
        raise ValueError("All elements in the input must be numbers, lists or arrays.")
    
    assert Pa.all() != Pb.all(), "The initial and final points must be different."
    assert max_seg_len > 0,'The maximun segment legth must be an positive number'
    #assert isinstance(n_points, ('None', 'int'))
    length = norm(Pa - Pb)
    assert length >= max_seg_len,'The distance between the points must be equal or higher than the maximum segment length'
    if n_points == None:
        n_points = int(round(length/max_seg_len))
    path = [
        linspace(Pa[0],Pb[0],n_points),
        linspace(Pa[1],Pb[1],n_points),
        linspace(Pa[2],Pb[2],n_points)]
    path = moveaxis(path, 0, 1)

    return path


def createArch(center: list, radius: float, start_angle: float, angle: float,
                max_seg_len: float, n_points=None, anticlockwise: bool = False):
    """
    Calculates the list of coordinates (coil path) in a specific arch in 3D space.

    :center| list or numpy.array: Coordinates of the initial point.
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
    assert n_points is None or isinstance(n_points, int), "n_points must be None or an integer"

    length = angle * radius
    n_points = int(length / max_seg_len) + 1 if n_points is None else n_points
    theta = angle / n_points

    if not anticlockwise:
        coilPath = [
            [
                center[0] + radius * sin(start_angle + i * theta),
                center[1] + radius * cos(start_angle + i * theta),
                center[2]
            ]
            for i in range(n_points + 1)
        ]
    else:
        coilPath = [
            [
                center[0] + radius * sin(start_angle - i * theta),
                center[1] + radius * cos(start_angle - i * theta),
                center[2]
            ]
            for i in range(n_points + 1)
        ]
    return coilPath

def helicoid(n: int, Pa,Pb,r:float, max_seg_len:float) :
    mold = createArch([0,0,0], radius=r, start_angle=0,angle=2*pi,max_seg_len=max_seg_len)

    x = [coordinate[0] for coordinate in mold]
    y = [coordinate[1] for coordinate in mold]
    z = [coordinate[2] for coordinate in createLine(Pa=Pa,Pb=Pb, max_seg_len=max_seg_len)]

    path = [[x_,y_,z_] for x_,y_,z_ in zip(x,y,z)]
    return path

def race_track(center,width: float, length: float, max_seg_len:float, int_radius:float):
    x_0,y_0,z_0 = center[0],center[1],center[2]
    C_1 = [x_0+width/2,y_0+length/2,z_0]
    C_2 = [x_0+width/2,y_0-length/2,z_0]
    C_3 = [x_0-width/2,y_0-length/2,z_0]
    C_4 = [x_0-width/2,y_0+length/2,z_0]

    arch_1 = createArch(C_1,int_radius,0,pi/2,max_seg_len=max_seg_len)
    arch_2 = createArch(C_2,int_radius,pi/2,pi/2,max_seg_len=max_seg_len)
    arch_3 = createArch(C_3,int_radius,pi,pi/2,max_seg_len=max_seg_len)
    arch_4 = createArch(C_4,int_radius,pi+pi/2,pi/2,max_seg_len=max_seg_len)

    line_1 = createLine(arch_1[-1],arch_2[1],max_seg_len=max_seg_len)
    line_2 = createLine(arch_2[-1],arch_3[1],max_seg_len=max_seg_len)
    line_3 = createLine(arch_3[-1],arch_4[1],max_seg_len=max_seg_len)
    line_4 = createLine(arch_4[-1],arch_1[1],max_seg_len=max_seg_len)

    path =  arch_1 + line_1 + arch_2 + line_2 + arch_3 + line_3 + arch_4 + line_4

    return path
def racetrack3d(center, inwidth, inlength, max_seg_len, int_radius, thickness):
    
    N_coils_h = int(thickness/max_seg_len)
    inwidth_ = inwidth
    inlength_ = inlength
    int_radius_ = int_radius
    
    Total_racetrack = []
    
    for i in range (N_coils_h):
        inwidth_ += 0.01
        inlength_ +=0.01
        int_radius_ += 0.01

        Total_racetrack += race_track(center, inwidth_, inlength_, max_seg_len, int_radius_)

    return Total_racetrack
