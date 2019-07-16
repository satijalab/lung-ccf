"""
    Drawing 3D - Library for moving a drawing cursor in a 3D numpy volume.
    by Tommaso Biancalani <tbiancal@broadinstitute.org>
"""

import numpy as np
import xml.etree.ElementTree as ET
import logging
from tqdm import tqdm

# comment line below to disable output debug
# logging.basicConfig(level=logging.ERROR)
logging.basicConfig(level=logging.INFO)


class Volume3D:
    """
    Allocate a 3D numpy array and provides primitives to draw on it.
    """
    def __init__(self, shape=(512, 512, 512)):
        """
        Allocate volumes with dimension `(x_pixels, y_pixels, z_pixels)`.
        """
        self.shape = shape
        self.voxels = np.zeros(self.shape, dtype=np.uint8)
        # self.branch_points = np.zeros(self.shape)
        # self.end_points = np.zeros(self.shape)
        logging.debug('Volume allocated with shape ' + str(self.shape))

    def _get_sphere(self, x, y, z, r):
        """
        Return coordinates of a sphere with center in `(x, y, z)` and radius `r.`
        :param x:
        :param y:
        :param z:
        :param r:
        :return:
        """
        # Pts should lie within allowed volume.
        x_min = np.max([0, x - r])
        x_max = np.min([x + r, self.shape[0]])
        y_min = np.max([0, y - r])
        y_max = np.min([y + r, self.shape[1]])
        z_min = np.max([0, z - r])
        z_max = np.min([z + r, self.shape[2]])

        # Get cube of points within allowed volume.
        x_grid, y_grid, z_grid = np.meshgrid(np.arange(x_min, x_max),
                                             np.arange(y_min, y_max),
                                             np.arange(z_min, z_max))

        # test all points to see which ones are within the radius of line thickness
        ixs = np.sqrt(np.power(x_grid - x, 2) + np.power(y_grid - y, 2) + np.power(z_grid - z, 2)) <= r

        return x_grid[ixs], y_grid[ixs], z_grid[ixs]

    def draw_sphere(self, x, y, z, r):
        """
        Draw a sphere centered in `(x, y, z)` and radius `r`.
        :param x:
        :param y:
        :param z:
        :param r:
        :return:
        """
        xs, ys, zs = self._get_sphere(x, y, z, r)
        for x, y, z in zip(xs, ys, zs):
            # TODO stammerda puo' essere ottimizzata
            x = int(round(x))
            y = int(round(y))
            z = int(round(z))
            x_oob_flag = x < 0 or x >= self.shape[0]
            y_oob_flag = y < 0 or y >= self.shape[1]
            z_oob_flag = z < 0 or z >= self.shape[2]
            if x_oob_flag or y_oob_flag or z_oob_flag:
                break
            self.voxels[z, y, x] = 255  # z first!

    def draw_line(self, start_position, end_position, radius):
        start_position = np.asarray(start_position)
        end_position = np.asarray(end_position)
        distance_vector = end_position - start_position
        distance_norm = np.linalg.norm(distance_vector)
        distance_versor = distance_vector / distance_norm

        distance_norm_int = round(int(distance_norm))
        if distance_norm_int == 0:
            logging.error('draw_line: vector is zero')
        logging.debug('Drawing line...')
        for t in tqdm(range(distance_norm_int)):
            center = start_position + t * distance_versor
            self.draw_sphere(*center, r=radius)


class Point:
    """
        A point in 3D space. Variables follow convention in figure.
        Importantly, `x, y, z` are always rounded to positive integers.
    """
    def __init__(self, **kwargs):
        """
        Define point by supplying either x, y, z or r, phi, theta.
        x,y,z should be positive integers (rounded to them otherwise).
        r is positive float
        phi (deg) between 0 - 360
        theta (deg) between 0 - 180
        :param kwargs:
        """
        if 'x' in kwargs and 'y' in kwargs and 'z' in kwargs:
            self.x = int(round(kwargs.get('x')))
            self.y = int(round(kwargs.get('y')))
            self.z = int(round(kwargs.get('z')))
            polar_dict = self.from_cartesian_to_polar(x=self.x, y=self.y, z=self.z)
            self.theta = polar_dict['theta']
            self.phi = polar_dict['phi']
            self.r = polar_dict['r']
        elif 'r' in kwargs and 'theta' in kwargs and 'phi' in kwargs:
            # Allocate polar coordinates
            self.phi = kwargs.get('phi') % 360
            self.theta = kwargs.get('theta') % 360
            if self.theta > 180 or self.theta < 0:
                logging.error('theta must be in between 0 and 180')
                raise NotImplementedError
            self.r = kwargs.get('r')
            # Convert coordinates to Cartesian and round them to closest integer
            cartesian_dict = self.from_polar_to_cartesian(r=self.r,
                                                          theta=self.theta,
                                                          phi=self.phi)
            self.x = int(round(cartesian_dict['x']))
            self.y = int(round(cartesian_dict['y']))
            self.z = int(round(cartesian_dict['z']))
            # Now reconvert them to polar
            polar_dict = self.from_cartesian_to_polar(x=self.x, y=self.y, z=self.z)
            self.theta = polar_dict['theta']
            self.phi = polar_dict['phi']
            self.r = polar_dict['r']
        else:
            logging.error('Arguments must be either Cartesian (x, y, z) or polar (r, theta, phi).')

    def get_cartesian_coordinates(self):
        return {'x': self.x, 'y': self.y, 'z': self.z}

    def get_polar_coordinates(self):
        return {'r': self.r, 'theta': self.theta, 'phi': self.phi}

    @staticmethod
    def from_cartesian_to_polar(*, x, y, z):
        """
        Convert xyz to polar coordinate. Return r, theta, phi.
        Angles are returned in degs.
        :param x:
        :param y:
        :param z:
        :return:
        """
        r = np.sqrt(x**2 + y**2 + z**2)
        if r == 0:
            return {'r': 0, 'theta': 0, 'phi': 0}
        else:
            theta = np.arccos(z/r)
            theta = np.rad2deg(theta)
            phi = np.arctan2(y, x)
            # if x == 0 and y < 0:
            #     phi = 3 * np.pi / 2
            # elif x == 0 and y >= 0:
            #     phi = np.pi / 2
            # elif x > 0:
            #     if y >= 0:
            #         phi = np.arctan(y/x)
            #     else:  # y < 0
            #         phi = np.arctan(-y/x) + (3*np.pi/2)
            # elif x < 0:
            #     if y >= 0:
            #         phi = np.arctan(-y/x) + (np.pi / 2)
            #     else:
            #         phi = np.arctan(y/x) + np.pi
            # else:
            #     raise NotImplementedError
            phi = np.rad2deg(phi)
            phi %= 360
            theta %= 360
            if theta < 0 or theta > 180:
                raise ArithmeticError
            return {'r': r, 'theta': theta, 'phi': phi}

    @staticmethod
    def from_polar_to_cartesian(*, r, theta, phi, round_off=True):
        """
        Convert r, theta, phi into Cartesian coordinates.
        Angles are provided in deg.
        :param r:
        :param theta:
        :param phi:
        :param round_off:
        Boolean flag. Whether to round off the cartesian coordinates.
        :return:
        """
        # Constraint polar angles into their domains
        theta %= 360
        phi %= 360
        if theta < 0 or theta > 180:
            logging.error('theta must be in between 0 and 180.')
            raise ValueError
        # Convert to radians
        phi = np.deg2rad(phi)
        theta = np.deg2rad(theta)
        logging.debug('phi: {}, theta: {}'.format(phi, theta))
        # Change coordinates
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)
        # Round numbers to closest integer
        if round_off is True:
            x = int(round(x))
            y = int(round(y))
            z = int(round(z))
        return {'x': x, 'y': y, 'z': z}


class Direction(Point):
    """
    A versor with application point that lies on a plane_versor.
    Inherits `Point` class to define application point.
    Add two angles to define versor direction.
    Define plane_versor by using orthogonal versor to plane_versor.
    """
    # Threshold below which floating numbers are considered zero
    threshold = 1e-14

    def __init__(self, **kwargs):
        """
        Bla bla bla. Use `n_x, n_y, n_z` for plane_versor.
        :param kwargs:
        """
        # Init application point.
        super().__init__(**kwargs)
        # Init facing direction
        if 'altitude' in kwargs and 'azimuth' in kwargs:
            self.altitude = kwargs.get('altitude')
            if self.altitude < 0 or self.altitude > 180:
                logging.error('altitude must be in between 0 and 180.')
            self.azimuth = kwargs.get('azimuth')
            self.facing_versor = None
            self._set_facing_versor()
        else:
            logging.error('Need `altitude` and `azimuth` arguments')
        # Init plane_versor using user-defined otherwise default
        if 'plane_versor' in kwargs:
            self.plane_versor = kwargs.get('plane_versor')
        else:
            self.plane_versor = (1, 0, 0)  # default

        # if `check_quantity` not zero, facing versor doesn't lie in branching plane
        check_quantity = np.abs(np.dot(self.plane_versor, self.facing_versor))
        check_quantity = self._chop(check_quantity)
        if check_quantity != 0:
            logging.error('Facing versor not in branching plane.')
            logging.error('Plane versor {}. Facing versor {}.'.format(self.plane_versor,
                                                                      self.facing_versor))
        # Init skeleton
        # self.current_node = SkeletonNode(x=x, y=y, z=z)  # .
        # self.skeleton = ET.ElementTree(self.current_node)

    def __repr__(self):
        cart_dict = self.get_cartesian_coordinates()
        str_x = 'x: ' + str(cart_dict['x'])
        str_y = ' y: ' + str(cart_dict['y'])
        str_z = ' z: ' + str(cart_dict['z'])
        str_pos = 'Position: ' + str_x + str_y + str_z + '\n'
        dir_dict = self.get_direction()
        str_az = ' azimuth: ' + str(dir_dict['azimuth'])
        str_al = ' altitude: ' + str(dir_dict['altitude'])
        str_dir = 'Facing angles: ' + str_az + str_al + '\n'
        str_dir_vec = 'Facing versor: ' + str(self.facing_versor) + '\n'
        str_plane = 'Branching plane versor :' + str(self.plane_versor) + '\n'
        return str_pos + str_dir + str_dir_vec + str_plane

    def get_plane_versor(self):
        nx, ny, nz = self.plane_versor
        return {'nx': nx, 'ny': ny, 'nz': nz}

    def _set_facing_versor(self):
        """
        Setter for `facing_versor` which gives the versor corresponding to angles
        azimuth and altitude (w/o accounting for application point).
        """
        cart_dict = self.from_polar_to_cartesian(r=1,
                                                 theta=self.altitude,
                                                 phi=self.azimuth,
                                                 round_off=False)
        x, y, z = cart_dict['x'], cart_dict['y'], cart_dict['z']
        x = self._chop(x)
        y = self._chop(y)
        z = self._chop(z)
        self.facing_versor = (x, y, z)
        logging.debug('`set_facing_versor`: x {}, y {} z {}'.format(x, y, z))

    def _chop(self, number_float):
        """
        Similar to `Chop` in Mathematica. If less than threshold, returns zero.
        :return:
        """
        threshold = self.__class__.threshold
        if np.abs(number_float) < threshold:
            return 0.
        else:
            return number_float

    def _set_azimuth_and_altitude(self):
        """
        Setter for attributes `azimuth` and `altitude`. Use attribute `facing_versor`
        to find the theta.
        :return:
        """
        x, y, z = self.facing_versor
        polar_dict = self.from_cartesian_to_polar(x=x, y=y, z=z)
        self.azimuth = polar_dict['phi']
        self.altitude = polar_dict['theta']

    def get_direction(self):
        return {'altitude': self.altitude, 'azimuth': self.azimuth}

    def _set_direction(self, *, azimuth=None, altitude=None):
        """
        Update direction to `new_azimuth` and/or `altitude` (deg). Do
        not update `plane_versor` accordingly therefore be careful.
        :param azimuth:
        :param altitude:
        :return:
        """
        altitude = self.altitude if altitude is None else altitude
        azimuth = self.azimuth if azimuth is None else azimuth
        if altitude < 0 or altitude > 180:
            logging.error('Altitude must be in between zero and 180.')
        self.altitude = altitude
        self.azimuth = azimuth
        self._set_facing_versor()
        if np.dot(self.plane_versor, self.facing_versor) > 10e-15:
            logging.error('Facing versor does not lie on branching plane!')

    def reverse_direction(self):
        """
        Reverse direction faced by cursor
        :return:
        """
        cartesian_dict = self.from_polar_to_cartesian(r=1,
                                                      theta=self.altitude,
                                                      phi=self.azimuth,
                                                      round_off=False,)
        x = -cartesian_dict['x']
        y = -cartesian_dict['y']
        z = -cartesian_dict['z']
        logging.debug('reverse_direction x {} y {} z {}'.format(x, y, z))
        polar_dict = self.from_cartesian_to_polar(x=x, y=y, z=z)
        logging.debug('radius after reverse_direction {}.'.format(polar_dict['r']))
        new_azimuth = polar_dict['phi']
        new_altitude = polar_dict['theta']
        self._set_direction(azimuth=new_azimuth, altitude=new_altitude)

    def move_forward(self, increment):
        """
        Move application point forward by `increment` along versor direction.
        If increment is too small to move cursor returns False.
        :param increment:
        :return:
        """
        assert increment > 0
        altitude = np.deg2rad(self.altitude)
        azimuth = np.deg2rad(self.azimuth)
        delta_x = increment * np.sin(altitude) * np.cos(azimuth)
        delta_y = increment * np.sin(altitude) * np.sin(azimuth)
        delta_z = increment * np.cos(altitude)

        old_position = (self.x, self.y, self.z)
        self.x += delta_x
        self.y += delta_y
        self.z += delta_z

        self.x = int(round(self.x))
        self.y = int(round(self.y))
        self.z = int(round(self.z))

        if (self.x, self.y, self.z) == old_position:
            logging.error('move_forward: position hasnt changed, increment is too small')
            return False
        else:
            return True

    @staticmethod
    def _get_rotation_matrix_from_axis_angle(*, axis, angle_rad):
        """
        Get rotation matrix from axis and theta from:
        https://en.wikipedia.org/wiki/Rotation_matrix
        The matrix computes the rotation around `axis` by and theta `angle_rad`.
        :param axis:
        tuple of three element defining plane versor
        :param angle_rad:
        CAREFUL: theta is in rads
        :return:
        """
        cos = np.cos(angle_rad)
        sin = np.sin(angle_rad)
        one_m_cos = 1. - cos
        ux, uy, uz = axis
        R_xx = cos + ((ux ** 2) * one_m_cos)
        R_xy = (ux * uy * one_m_cos) - (uz * sin)
        R_xz = (ux * uz * one_m_cos) + (uy * sin)
        R_yx = (uy * ux * one_m_cos) + (uz * sin)
        R_yy = cos + ((uy ** 2) * one_m_cos)
        R_yz = (uy * uz * one_m_cos) - (ux * sin)
        R_zx = (uz * ux * one_m_cos) - (uy * sin)
        R_zy = (uz * uy * one_m_cos) + (ux * sin)
        R_zz = cos + ((uz ** 2) * one_m_cos)

        R = np.array([[R_xx, R_xy, R_xz],
                      [R_yx, R_yy, R_yz],
                      [R_zx, R_zy, R_zz]], dtype='float64')
        return R

    def turn(self, angle):
        """
        Turn facing direction by `theta` by remaining onto branching plane_versor.
        Positive angles turn right (clockwise), negative angles turn left, with
        versor normal to plane pointing north.
        :param angle:
        Angle in deg in between -180 and 180.
        :return:
        """
        if angle > 180 or angle < -180:
            logging.error('Angle should be in between -180 and 180')
            raise ValueError
        angle = np.deg2rad(angle)
        R = self._get_rotation_matrix_from_axis_angle(axis=self.plane_versor, angle_rad=angle)
        logging.debug('Rotation matrix is {}'.format(R))

        self._set_facing_versor()  # TODO it is safe to remove if code is bugless
        ax, ay, az = np.dot(R, self.facing_versor)
        ax = self._chop(ax)
        ay = self._chop(ay)
        az = self._chop(az)
        self.facing_versor = (ax, ay, az)
        logging.debug('turn: new facing versor {}'.format(self.facing_versor))
        self._set_azimuth_and_altitude()

    def rotate_plane(self, angle):
        """
        Rotate branching plane by `theta` (deg). Rotation occurs around direction defined
        by `facing_versor` applied to current cursor position. Rotation follows
        right-hand rule.
        :param angle:
        :return:
        """
        angle_rad = np.deg2rad(angle)
        self._set_facing_versor()  # TODO safe to remove is code is bugless
        R = self._get_rotation_matrix_from_axis_angle(axis=self.facing_versor,
                                                      angle_rad=angle_rad)
        vx, vy, vz = np.dot(R, self.plane_versor)
        vx = self._chop(vx)
        vy = self._chop(vy)
        vz = self._chop(vz)
        self.plane_versor = (vx, vy, vz)


class Cursor(Direction):
    """
    Manipulate cursor in 3D numpy volume.
    """
    def __init__(self,
                 volume=None,  # pointer to existing Volume3D
                 volume_size=(512, 512, 512),
                 cursor_position=(256, 256, 510),
                 plane_versor=(1, 0, 0),
                 azimuth=0,
                 altitude=180,
                 ):
        """
        Initialize cursor in 3D volume, with  position `x, y, z` and
        facing direction given by `azimuth, altitude` angles (deg).
        :param x:
        positive integer
        :param y:
        positive integer (coordinate along y axis)
        :param z:
        positive integer
        :param azimuth:
        :param altitude:
        """

        x, y, z = cursor_position
        # Init `Direction` using Cartesian coordinates
        if x is None or x < 0:
            logging.error('x coordinate must be a positive integer')
            raise ValueError
        if y is None or y < 0:
            logging.error('y coordinate must be a positive integer')
            raise ValueError
        if z is None or z < 0:
            logging.error('z coordinate must be a positive integer')
            raise ValueError
        if altitude is None or altitude < 0 or altitude > 180:
            logging.error('altitude must be in between 0 and 180.')
            raise ValueError
        if azimuth is None or azimuth < 0 or azimuth > 360:
            logging.error('azimuth must be in between 0 and 360.')
            raise ValueError
        super().__init__(x=x, y=y, z=z, altitude=altitude, azimuth=azimuth, plane_versor=plane_versor)

        if volume is None:
            # Allocate new 3D Volume
            self.volume = Volume3D(shape=volume_size)
        else:  # use existing volume
            self.volume = volume

    def draw_forward(self, increment, radius):
        """
        Move cursor forward by `increment` along facing direction,
        and draw corresponding line on numpy volume.
        :param increment:
        :param radius:
        :return:
        """
        # TODO program skeleton
        # add node to skeleton and update skeleton position
        # self.current_node = SkeletonNode(*end_coord, parent_node=self.current_node)

        # update cursor
        start_position = self.x, self.y, self.z
        self.move_forward(increment)
        end_position = self.x, self.y, self.z

        self.volume.draw_line(start_position=start_position, end_position=end_position, radius=radius)
        logging.debug('draw_forward: position {}'.format((self.x, self.y, self.z)))

