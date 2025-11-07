"""Geometry Module.

This module contains all geometric tools that can be useful for
resistance calculation and coil path generations.

Contains:
Area calculations: circle, rectangle, square.
Geometrics figure coordinates generators: arch, line, racetrack.
"""
import numpy as np
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
      

def line(Pa, Pb,*, max_seg_len:float = 1,n_points:int = None):
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
        Pa = np.array([float(coordinate) for coordinate in Pa])
        Pb = np.array([float(coordinate) for coordinate in Pb])
        max_seg_len = float(max_seg_len)

    except:
        raise ValueError("All elements in the input must be numbers, lists or arrays.")
    
    assert not np.array_equal(Pa, Pb), "The initial and final points must be different."
    assert max_seg_len > 0,'The maximun segment legth must be an positive number'
    #assert isinstance(n_points, ('None', 'int'))
    length = norm(Pa - Pb)
    assert length >= max_seg_len,'The distance between the points must be equal or higher than the maximum segment length'
    if n_points is None:
        n_points = int(np.ceil(length / max_seg_len)) + 1
        path = np.moveaxis(np.linspace(Pa, Pb, n_points),0,1)
    else:
        path = np.moveaxis(np.linspace(Pa, Pb, n_points),0,1)

    return path

def arc(center: list, radius: float, start_angle: float, angle: float,
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

    if n_points is not None:
        angles = np.linspace(start_angle, angle, n_points)
    else:
        arc_length = abs(angle) * radius
        n_segments = int(arc_length / max_seg_len) + 1
        angles = np.linspace(start_angle, start_angle + angle, n_segments)
    if not anticlockwise:
        x = center[0] + radius * np.sin(angles)
        y = center[1] + radius * np.cos(angles)
    else:
        x = center[0] - radius * np.sin(angles)
        y = center[1] - radius * np.cos(angles)
    z = np.full_like(angles, center[2])
    path = np.stack([x, y, z], axis=1)
    return path

def helicoid(n: int, Pa, z_len, r: float, max_seg_len: float):
    """Returns an helicoid path

    :n|int: number of turns of the helicoid
    :Pa|ArrayLike: origin of the helicoid object
    :z_len|float: the length of the helicoid trought the dimension that is not arch
    :max_seg_len|float: the max distance between the points

    Returns:
        path of the helicoid
    """
    total_angle = n * 2 * np.pi
    
    spiral_mold = np.array(arc(Pa, radius=r, start_angle=0,
                               angle=total_angle,
                               max_seg_len=max_seg_len))
    x = spiral_mold[:, 0]
    y = spiral_mold[:, 1]
    z = np.linspace(Pa[2], Pa[2] + z_len, num=x.shape[0])

    path = np.stack([x, y, z], axis=1)
    return path

def race_track(center,width: float, length: float, max_seg_len:float, int_radius:float):
    x_0, y_0, z_0 = center[0], center[1], center[2]

    centers = np.array([
        [x_0 + width / 2, y_0 + length / 2, z_0],
        [x_0 + width / 2, y_0 - length / 2, z_0],
        [x_0 - width / 2, y_0 - length / 2, z_0],
        [x_0 - width / 2, y_0 + length / 2, z_0]
    ])

    arc_angles = np.array([0, np.pi / 2, np.pi, 3 * np.pi / 2])
    arcs = [arc(centers[i], int_radius, arc_angles[i], np.pi/2, max_seg_len) for i in range(4)]
    lines = [
        line(arcs[0][-1], arcs[1][1], max_seg_len=max_seg_len),
        line(arcs[1][-1], arcs[2][1], max_seg_len=max_seg_len),
        line(arcs[2][-1], arcs[3][1], max_seg_len=max_seg_len),
        line(arcs[3][-1], arcs[0][1], max_seg_len=max_seg_len)
    ]
    path = np.concatenate([arcs[0], lines[0], arcs[1], lines[1], arcs[2], lines[2], arcs[3], lines[3]], axis=0)

    return path

def racetrack2d(center, inwidth, inlength, max_seg_len, int_radius, thickness):
    N_coils_h = int(thickness / max_seg_len)

    inwidths = inwidth + max_seg_len * np.arange(N_coils_h)
    inlengths = inlength + max_seg_len * np.arange(N_coils_h)
    int_radii = int_radius + max_seg_len * np.arange(N_coils_h)

    racetracks = [race_track(center, w, l, max_seg_len, r) for w, l, r in zip(inwidths, inlengths, int_radii)]
    return np.concatenate(racetracks, axis=0)

def racetrack3d(center, inwidth, inlength, max_seg_len, int_radius, thickness, height):
    
    if height == 0:
        return racetrack2d(center, inwidth, inlength, max_seg_len, int_radius, thickness)
    
    N_layers = int(height / max_seg_len)                                            
    z_space = np.arange(N_layers) * max_seg_len                           
    z_coords = [center[2] + z for z in z_space]
    all_layers = [
        racetrack2d(np.array([center[0], center[1], z]), inwidth, inlength, max_seg_len, int_radius, thickness)
        for z in z_coords
    ]

    return np.concatenate(all_layers, axis=0)