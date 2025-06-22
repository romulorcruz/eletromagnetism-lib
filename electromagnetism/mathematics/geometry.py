"""Geometry Module.

This module contains all geometric tools that can be useful for
resistance calculation and coil path generations.

Contains:
Area calculations: circle, rectangle, square.
Geometrics figure coordinates generators: arch, line, racetrack.
"""
from numpy import sin, cos, pi, array, linspace

def circle_area(radius: float) -> float:
    """Calculates the area of a circle based on its radius.
    
    Args:
        radius (float): The radius of the circle in meters.
        
    Returns:
        float: The area of the circle in squared meters.
    """
    if radius <= 0:
        raise ValueError("The radius must be a positive number.")
    return pi * radius ** 2


def rectangle_area(width: float,length: float) -> float:
    """Calculates the area of a rectangle based on its width and length.
    
    Args:
        width (float): The width of the rectangle in meters.
        length (float): The length of the rectangle in meters.
        
    Returns:
        float: The area of the rectangle in squared meters.
    """
    if width <= 0:
        raise ValueError("The width must be a positive number.") 
    
    if length <= 0:
        raise ValueError("The length must be a positive number.") 
    
    return width * length

def square_area(side: float) -> float:
    """Calculates the area of a square based on its side.
    
    Args:
        side (float): The side of the square in meters.
        
    Returns:
        float: The area of the square in squared meters.
    """
    if side <= 0:
        raise ValueError("The side must be a positive number") 
    return side**2

def create_line(Pa: float, Pb: float,*, max_seg_len:float = 1,n_points:int = None) -> array:
    """Calculates the list of coordinates (coil path) between two different points in 3D space.

        Args:
            Pa (list or numpy.array): Coordinates of the initial point.
            Pb (list or numpy.array): Coordinates of the final point.
            max_seg_len (float): (optional) Maximum length of each segment. Default is 1.
            n_points (int): (optional) The number of points in the coil path. Default is None.


        Returns:
            path: A numpy.array of coordinates representing the segmented line.
    """
    # Transform all the objects in float if they are compactible, otherwise it will raise a value error
    try:
        Pa = [float(coordinate) for coordinate in Pa]
        Pb = [float(coordinate) for coordinate in Pb]
        max_seg_len = float(max_seg_len)

    except:
        raise ValueError("All elements in the input must be numbers, lists or arrays.")
    if Pa == Pb:
        raise ValueError("The initial and final points must be different.")
    if max_seg_len <= 0:
        raise ValueError("The maximun segment legth must be an positive number")
    if n_points is not None and isinstance(n_points, int):
        raise TypeError("The number of points must be a positive integer.")
    if n_points is not None and n_points <= 0:
        raise ValueError("The number of points must be a positive integer.")

    projection = [b - a for a, b in zip(Pa, Pb)]
    length = sum(p**2 for p in projection) ** 0.5

    if max_seg_len >= length:
        raise ValueError("The distance between the points must be equal or higher than the maximum segment length")

    max_seg_len = length/n_points if n_points != None else max_seg_len

    ratio = length / max_seg_len
    path = [
        [projection[0] / ratio * i + Pa[0],
         projection[1] / ratio * i + Pa[1],
         projection[2] / ratio * i + Pa[2]]
        for i in range(int(ratio) + 1)
    ]

    if Pb not in path:
        path.append(Pb)

    return array(path)


def create_arch(center: list, radius: float, start_angle: float, angle: float,
                max_seg_len: float, n_points=None, anticlockwise: bool = False):
    """Calculates the list of coordinates (coil path) in a specific arch in 3D space.

    center (list or numpy.array): Coordinates of the initial point.
    radius (float): Radius of the arch.
    start_angle (float): Starting angle (radians).
    angle (float): Total angle (radians) to sweep.
    max_seg_len (float): (optional) Maximum length of each segment. Default is 1.
    n_points (int): (optional) Number of points. If None, computed from max_seg_len.
    anticlockwise (bool): If True, the arc is swept in the negative angular direction.

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
        raise TypeError("All elements in the input must be numbers, lists or arrays.")
    
    if radius <= 0 or angle <= 0:
        raise ValueError("The radius and angle must be positive numbers.")
    if n_points is not None and not isinstance(n_points, int):
        raise TypeError("The number of points must be a positive integer.")
    if n_points is not None and n_points <= 0:
        raise ValueError("The number of points must be a positive integer.")

    length = angle * radius
    n_points = int(length / max_seg_len) + 1 if n_points is None else n_points
    theta = angle / n_points

    if not anticlockwise:
        path = [
            [
                center[0] + radius * cos(start_angle + i * theta),
                center[1] + radius * sin(start_angle + i * theta),
                center[2]
            ]
            for i in range(n_points + 1)
        ]

    else:
        path = [
            [
                center[0] + radius * cos(start_angle - i * theta),
                center[1] + radius * sin(start_angle - i * theta),
                center[2]
            ]
            for i in range(n_points + 1)
        ]
    return path
