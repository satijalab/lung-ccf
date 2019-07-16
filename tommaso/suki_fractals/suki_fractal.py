"""
    Library to create fractals from Bela Suki paper.
    by Tommaso Biancalani <tbiancal@broadinstitute.org>
"""

from drawing_3d import *


class BranchPoint:
    """
    A point in a branch. Typically an extremity, but class can be inherited to
    subclass branching points (see `SkeletonPoint`).
    """
    def __init__(self, *, cursor_position, branch_diameter=None,
                 theta=None, azimuth=None, altitude=None, generation_number=0):
        """

        :param branch_diameter:
        Diameter of the branch at this point (it changes within a `FractalBranch`).
        :param cursor_position:
        :param theta:
        If child point, the corresponding theta.
        :param azimuth:
        azimuth of facing cursor (unneeded if parent point).
        :param altitude:
        altitude of facing cursor (unneeded if parent point).
        :param generation_number:
        All `BranchingPoint` has a generation number which should be the same
        of the corresponding `SkeletonPoint`.
        """
        self.branch_diameter = branch_diameter
        self.cursor_position = cursor_position
        self.theta = theta
        self.azimuth = azimuth  # Facing direction of cursor
        self.altitude = altitude
        self.generation_number = generation_number

    def is_out(self, volume_size):
        """
        Return True if point is out of volume
        :return:
        """
        x, y, z = self.cursor_position
        if x < 0 or y < 0 or z < 0:
            return True
        if x > volume_size[0] or y > volume_size[1] or z > volume_size[2]:
            return True
        return False


class SkeletonPoint(ET.Element, BranchPoint):
    """
    Branching point that constitutes the skeleton of fractal.
    """
    def __init__(self, cursor_position, branch_diameter=None, theta=None,
                 azimuth=None, altitude=None, parent_point=None, generation_number=0):
        """
        :param cursor_position:
        :param parent_point:
        Branching point that is parent of this one.
        """
        if parent_point is None:
            name = 'root'
        else:
            name = 'node'

        # it is too confusing to use `super` here
        ET.Element.__init__(self, name)
        BranchPoint.__init__(self,
                             cursor_position=cursor_position,
                             branch_diameter=branch_diameter,
                             theta=theta,
                             azimuth=azimuth,
                             altitude=altitude,
                             generation_number=generation_number,
                             )

        # Save point labels.
        # attributes need to be string to perform XML serialization
        self.set('x', str(cursor_position[0]))
        self.set('y', str(cursor_position[1]))
        self.set('z', str(cursor_position[2]))

        # Set branch generation
        if parent_point is not None:
            # append as child of parent
            parent_point.append(self)
            self.generation_number = generation_number
        self.set('generation_number', str(self.generation_number))
        logging.debug('Creating {} at {} with depth {}'.format(name, cursor_position, self.generation_number))
        if parent_point is not None:
            logging.debug('Parent node has generation {}'.format(parent_point.generation_number))


class FractalBranch:
    """
    Planar fractal branch consisting of a parent point
    and two children points.
    """
    def __init__(self, *, point_0: BranchPoint,
                 point_1: BranchPoint,
                 point_2: BranchPoint,
                 plane_versor: tuple,):

        # Branch points on branch
        self.point_0 = point_0
        self.point_1 = point_1
        self.point_2 = point_2

        # self.generation = generation
        self.plane_versor = plane_versor


class SukiFractal:
    """
    Implement fractal from paper
    "A three-dimensional model of the human airway tree" - 1999
    """
    def __init__(self, n=2.8, r=.4, d0=100, max_generation=2, volume_size=(100, 100, 100)):

        # Parameters from paper
        self.n = n  # diameter exponent
        self.r = r  # volume ratio parameter

        # Numpy volume onto which to write
        self.volume = Volume3D(shape=volume_size)

        # Draw branches recursively until `max_generation` is reached
        self.max_generation = max_generation
        max_x, max_y, max_z = volume_size

        # Initial cursor position
        cursor_position = (max_x//2, max_y//2, max_z-10)

        # draw parent node of skeleton
        self.current_node = SkeletonPoint(cursor_position=cursor_position,
                                          branch_diameter=d0,
                                          azimuth=90,
                                          altitude=180,
                                          )
        self.skeleton = ET.ElementTree(self.current_node)

        self.draw_lung_recursively(point=self.current_node,
                                   plane_versor=(1, 0, 0),
                                   )

    def draw_lung_recursively(self, point: BranchPoint, plane_versor):
        """
        Draw branches starting from parent `point` having `plane_versor`.
        :param point:
        :param plane_versor:
        :return:
        """

        if point.generation_number == self.max_generation:
            # quit if max generation is achieved
            logging.info('Maximum generation achieved.')
            return
        else:
            # Print debug info
            plane_versor_rounded = [round(x, 2) for x in plane_versor]
            logging.info('Drawing branch on plane {}.'.format(plane_versor_rounded))

            # Draw branch on plane
            branch = self.draw_branch(point_0=point,
                                      plane_versor=plane_versor,
                                      )
            if branch is None:
                # Reach the end of volume
                logging.info('Branch reached end of volume.')
                return

            # Take facing versor from point one
            facing_versor_1 = self._get_facing_versor(azimuth=branch.point_1.azimuth,
                                                      altitude=branch.point_1.altitude)
            # Get rotated plane versor for point one
            plane_versor_1 = self.__class__._get_rotated_plane_versor(angle=90,
                                                                      facing_versor=facing_versor_1,
                                                                      plane_versor=branch.plane_versor)

            # Draw branch into point one
            self.draw_lung_recursively(point=branch.point_1,
                                       plane_versor=plane_versor_1,
                                       )
            logging.debug('Now drawing pt2')
            # Take facing versor from point two
            facing_versor_2 = self._get_facing_versor(azimuth=branch.point_2.azimuth,
                                                      altitude=branch.point_2.altitude)
            # Get rotated plane versor for point two
            plane_versor_2 = self.__class__._get_rotated_plane_versor(angle=90,
                                                                      facing_versor=facing_versor_2,
                                                                      plane_versor=branch.plane_versor)
            self.draw_lung_recursively(point=branch.point_2,
                                       plane_versor=plane_versor_2,
                                       )

    def draw_branch(self, point_0: BranchPoint, plane_versor):
        """
        Draw a planar `FractalBranch` starting from `point_0`.
        :param point_0:
        :param plane_versor:
        :return:
        """
        current_generation_number = self.current_node.generation_number

        # logging.info('Drawing branch with diameter {}.'.format(point_0.branch_diameter))
        # Init Cursor to draw on branch
        cursor = Cursor(volume=self.volume,
                        cursor_position=point_0.cursor_position,
                        plane_versor=plane_versor,
                        azimuth=point_0.azimuth,
                        altitude=point_0.altitude,
                        )

        # Draw from parent point to bifurcation point
        length = point_0.branch_diameter
        cursor.draw_forward(increment=length, radius=point_0.branch_diameter // 2)
        
        # If new cursor position is out of bounds, end branch by returning None
        if self._is_out(cursor_position=(cursor.x, cursor.y, cursor.z),
                        volume_size=self.volume.shape) is True:
            return None
            

        # update skeleton
        self.current_node = SkeletonPoint(cursor_position=(cursor.x, cursor.y, cursor.z),
                                          generation_number=current_generation_number+1,
                                          parent_point=self.current_node)

        # Compute diameters and angles of children points
        diameter_1, diameter_2 = self.get_children_diameters(d0=point_0.branch_diameter,
                                                             r=self.r,
                                                             n=self.n)
        angle_1, angle_2 = self.get_branching_angles(n=self.n, r=self.r)

        # Move to child point 1 and store it
        cursor.turn(angle_1)
        length_1 = 2 * diameter_1
        cursor.draw_forward(increment=length_1, radius=diameter_1//2)

        # Exit if next pt falls outside volume
        if self._is_out(cursor_position=(cursor.x, cursor.y, cursor.z),
                        volume_size=self.volume.shape) is True:
            return None

        point_1 = BranchPoint(branch_diameter=diameter_1,
                              cursor_position=(cursor.x, cursor.y, cursor.z),
                              theta=angle_1,
                              altitude=cursor.altitude,
                              azimuth=cursor.azimuth,
                              generation_number=current_generation_number+1,
                              )

        # Come back to bifurcation point. Move to point 2 and save it.
        cursor.reverse_direction()
        cursor.move_forward(length_1)
        cursor.reverse_direction()
        cursor.turn(-angle_1)
        cursor.turn(-angle_2)
        length_2 = 2 * diameter_2
        cursor.draw_forward(increment=length_2, radius=diameter_2//2)

        if self._is_out(cursor_position=(cursor.x, cursor.y, cursor.z),
                        volume_size=self.volume.shape) is True:
            return None

        point_2 = BranchPoint(branch_diameter=diameter_2,
                              cursor_position=(cursor.x, cursor.y, cursor.z),
                              theta=angle_2,
                              altitude=cursor.altitude,
                              azimuth=cursor.azimuth,
                              generation_number=current_generation_number+1,
                              )

        # Save fractal branch
        branch = FractalBranch(point_0=point_0, point_1=point_1, point_2=point_2,
                               plane_versor=plane_versor)

        return branch

    @staticmethod
    def _is_out(cursor_position, volume_size):
        """
        Return True if `cursor_position` falls outside `volume_size`.
        :param cursor_position:
        :param volume_size:
        :return:
        """
        x, y, z = cursor_position
        if x < 0 or y < 0 or z < 0:
            return True
        if x > (volume_size[0] - 1) or \
           y > (volume_size[1] - 1) or \
           z > (volume_size[2] - 1):
            return True
        
        return False

    @staticmethod
    def _get_facing_versor(*, azimuth, altitude):
        """
        Get `facing_versor` from  azimuth and altitude
        (w/o accounting for application point).
        """
        cart_dict = Point.from_polar_to_cartesian(r=1,
                                                 theta=altitude,
                                                 phi=azimuth,
                                                 round_off=False)
        x, y, z = cart_dict['x'], cart_dict['y'], cart_dict['z']
        facing_versor = (x, y, z)
        return facing_versor

    @staticmethod
    def get_children_diameters(*, d0, r, n):
        """
        Eq. (4) - Return diameters as integers
        :param d0:
        :param r:
        :param n:
        :return:
        """
        d1 = d0 * np.power(r, 1/n)
        d2 = d0 * np.power((1.-r), 1/n)
        d1 = int(round(d1))
        d2 = int(round(d2))
        return d1, d2

    @staticmethod
    def get_branching_angles(*, r, n):
        """
        Eq. (6) and returns angles in deg
        :param r:
        :param n:
        :return:
        """
        r4n = np.power(r, 4/n)
        one_m_r4n = np.power((1-r), 4/n)
        r2n = np.power(r, 2/n)
        one_m_r2n = np.power((1-r), 2/n)
        cos_theta_1 = (1. + r4n - one_m_r4n) / (2 * r2n)
        cos_theta_2 = (1. + one_m_r4n - r4n) / (2 * one_m_r2n)
        theta_1 = np.rad2deg(np.arccos(cos_theta_1))
        theta_2 = np.rad2deg(np.arccos(cos_theta_2))
        return theta_1, theta_2

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

    @classmethod
    def _get_rotated_plane_versor(cls, angle, facing_versor, plane_versor):
        """
        Rotate branching plane (ie `plane_versor` by `theta` (deg).
        Rotation occurs around direction defined by `facing_versor`
        applied to current cursor position. Rotation follows
        right-hand rule. (Method taken from `Cursor` class).
        :param angle:
        in deg
        :param facing_versor:
        tuple
        :param plane_versor:
        :return:
        """
        angle_rad = np.deg2rad(angle)
        rotation_matrix = cls._get_rotation_matrix_from_axis_angle(axis=facing_versor,
                                                                    angle_rad=angle_rad)
        vx, vy, vz = np.dot(rotation_matrix, plane_versor)
        vx /= np.linalg.norm((vx, vy, vz))  # TODO should be already normalized
        vy /= np.linalg.norm((vx, vy, vz))
        vz /= np.linalg.norm((vx, vy, vz))
        return vx, vy, vz
