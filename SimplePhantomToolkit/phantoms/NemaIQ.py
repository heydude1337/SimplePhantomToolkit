import numpy as np
from copy import copy
import pandas as pd
from SimplePhantomToolkit import Phantoms, Sphere, Cylinder, Phantom, geometry, NemaIQPlot

class NemaIQ(Phantom):

# =============================================================================
#     Physical dimensions of the phantom, measured on CT
# =============================================================================
    # length of active compartment
    PHANTOM_HEIGHT = 195
    VOLUME = 9700 # volume of large compartment

    # spheres
    SPHERE_DIAMETERS = (37, 28, 22, 17, 13, 10)
    # angle for each sphere location, largest sphere at 0 degrees
    SPHERE_POSTION_ANGLE = (0, 60, 120, 180, 240, 300)
        # radius of sphere along which spheres are positions
    SPHERE_LARGE_DIAMETER = 115

    SPHERE_OFFSET = 70 # z distance spheres from top lid of phantom
    CLOCKWISE = 0
    ANTI_CLOCKWISE = 1

    # Lung insert
    LUNG_INSERT_DIAMETER = 45
    LUNG_OFFSET_LID = 5

    # Uniform regions
    UNIFORM_DIAMETER = 50
    UNIFORM_LENGTH = 75
    # offset from the sphere plane, equal to largest sphere, also offset from
    # bottom
    UNIFORM_OFFSET = 37
    UNIFORM_DIAMETER = 50
    UNIFORM_LARGE_DIAMETER = 140

    # orientation of phantom in the z direction
    UP = 0
    DOWN = 1

    _LOG_LEVEL = Phantoms.LEVEL_DEBUG

    def __init__(self, image = None, threshold = 0.5,
                 up_down = 0,
                 rotation = 0,
                 anchor = (0, 0, 0),
                 sphere_angle = 0):

        """ Nema Image Quality Phantom. Phantom should by aligned along the
        z-axis, rotation should be 0 or 1 defining if the spheres are oriented
        in a clockwise direction or anti clockwise direction. anchor is the
        center around which the spheres are located. Angle is the angle at
        which the largest sphere is located. All other elements are located
        relative to this position. """


        super(NemaIQ, self).__init__(anchor = anchor, image = image)

        self.up_down = up_down
        self.rotation = rotation
        self.sphere_angle = sphere_angle


        # cache properties that are expensive to calculate
        self._phantom = None
        self._spheres = None
        self._uniform_regions = None
        self._lung_insert = None

        msg = 'NemaIQ phantom with orientation %s, rotation %s at position %s'
        self.logger.debug(msg, str(up_down), str(rotation), str(anchor))

    def __str__(self):
        text = super(NemaIQ, self).__str__()
        text += '\norientation: {0}\n rotation: {1}'.format(self.up_down,
                               self.rotation)
        return text

    def plot(self):
        return NemaIQPlot(self).plot()

    @property
    def sphere_angle(self):
        """ The angle in the x -y planeat which the largest sphere can
        be located. Zero angle means the largest spheres lies on the positive
        x-axis """
        return self._sphere_angle

    @sphere_angle.setter
    def sphere_angle(self, value):
        self._sphere_angle = value
        self.reset()

    @property
    def up_down(self):
        """ Orientation of the phantom along the z-axis. up_down is 0 if the
        spheres are positioned above the centroid of the phantom and 0 if the
        speres are placed below the centroid of the phantom (= upside down)
        """
        return self._up_down

    @up_down.setter
    def up_down(self, value):
        self._up_down = value
        self.reset()

    @property
    def rotation(self):
        """ rotation is 0 for a clockwise orientation of the spheres in the
        x-y plane and 1 for a anti clockwise orientation. """
        return self._rotation

    @rotation.setter
    def rotation(self, value):
        self._rotation = value
        self.reset()

    @property
    def lung_insert(self):
        """ Return a cylinder phantom that represents the lung insert """
        if self._lung_insert is None:
            length =  NemaIQ.PHANTOM_HEIGHT - self.LUNG_OFFSET_LID

            z_top = NemaIQ.SPHERE_OFFSET - self.LUNG_OFFSET_LID
            z_bottom = z_top - length

            if self.up_down == NemaIQ.DOWN:
                z_top *= -1
                z_bottom *= -1

            cylinder = Cylinder(image = self.image,
                                axis = [[0,0, z_top], [0,0, z_bottom]],
                                radius = NemaIQ.LUNG_INSERT_DIAMETER / 2,
                                anchor = self.anchor)
            self.lung_insert = cylinder
        return self.lung_insert

    @property
    def spheres(self):
        """ Return phantom elements for the spheres that are derived from a
        region growing algorithm (within SimpleITK). """
        if self._spheres is None:
            anchors = NemaIQ._sphere_positions(self.anchor,
                                               self.sphere_angle,
                                               self.rotation,
                                               self.up_down)


            spheres = Phantoms()
            max_d = np.max(self.SPHERE_DIAMETERS)
            for anchor in anchors:
                sphere = GrowRegionFromSphere(image = self.image,
                                              anchor = anchor,
                                              radius = max_d /2,
                                              background = self.mean_background)
                spheres.append(sphere)
            self._spheres = spheres

        return self._spheres

    @property
    def uniform_regions(self):
        """ Return cylindrical phantom elements representing uniform regions
        in the background """
        if self._uniform_regions is None:
            axes, anchors = NemaIQ._uniform_positions(self.anchor,
                                                      self.sphere_angle,
                                                      up_down = self.up_down)

            cylinders = Phantoms()

            for axis, anchor in zip(axes, anchors):
                cylinder = Cylinder(image = self.image,
                                    axis = axis,
                                    anchor = anchor,
                                    radius = self.UNIFORM_DIAMETER/2)
                cylinders.append(cylinder)
            self._uniform_regions = cylinders
        return self._uniform_regions

    @property
    def phantom(self):
        """ Return all phantom elements (spheres, background and lung insert).
        """
        if self._phantom is None:
            self.logger.info('Expensive Phantom Generation')
            _phantom = self.spheres + self.uniform_regions #+ self.lung_insert
            self._phantom = _phantom
        return self._phantom

    @property
    def mask(self):
        return self.phantom.mask

    @property
    def mean_background(self):
        """ Returns the mean of the mean values in the uniform regions """
        return np.mean(self.uniform_regions.analyze()['mean'])

    def reset(self):
        self.logger.info('Reset')
        super(NemaIQ, self).reset()
        self._phantom = None
        self._spheres = None
        self._uniform_regions = None
        self._lung_insert = None


    def analyze(self):
        results = self.uniform_regions.analyze()
        for diameter, sphere in zip(self.SPHERE_DIAMETERS, self.spheres):
            result = sphere.analyze()
            result['diameter'] = diameter
            result['physical_volume_ml'] = 4/3 * np.pi * (diameter/2) ** 3 / 1000
            result['measured_volume_ml'] = np.prod(self.image.GetSpacing()) * result['count'] / 1000
            results = pd.concat((results, result))
        return results


    @staticmethod
    def _uniform_positions(anchor, angle0, up_down = 0):
        """ Position cylinders in the x-y plane. Anchor is the point around all
        spheres are located. Cylinders will be placed in the uniform
        background. Angle0 is the angle that locates the largest
        sphere relative to the positive x-axis.  """

        angles = np.asarray(NemaIQ.SPHERE_POSTION_ANGLE) / 180 * np.pi + angle0
        radius = NemaIQ.UNIFORM_LARGE_DIAMETER /2

        anchors = np.ones((len(angles), 3))
        anchors[:, 0] = radius * np.cos(angles) + anchor[0]

        anchors[:, 1] = radius * np.sin(angles) + anchor[1]

        anchors[:, 2] = anchor[2]

        if up_down == NemaIQ.UP:
            z_top = -1 * NemaIQ.UNIFORM_OFFSET
            z_bottom = z_top - NemaIQ.UNIFORM_LENGTH

        elif up_down == NemaIQ.DOWN:
            z_top = NemaIQ.UNIFORM_OFFSET
            z_bottom = z_top + NemaIQ.UNIFORM_LENGTH

        axes = []
        for anchor in anchors:
            axis = [[0, 0, z_top], [0,0, z_bottom]]
            axes.append(axis)

        return np.asarray(axes), anchors

    @staticmethod
    def _sphere_positions(anchor, angle0=0, rotation = 0, up_down = 0):
        """ Position spheres in the x-y plane. Anchor is the point around all
        spheres will be located. Angle0 is the angle that locates the largest
        sphere relative to the positive x-axis. """

        angles = np.asarray(NemaIQ.SPHERE_POSTION_ANGLE) / 180 * np.pi
        radius = NemaIQ.SPHERE_LARGE_DIAMETER /2

        if rotation == NemaIQ.CLOCKWISE:
            angles *= -1

        if up_down == NemaIQ.DOWN:
            angles *= -1

        angles += angle0

        anchors = np.ones((len(angles), 3))

        anchors[:, 0] = radius * np.cos(angles) + anchor[0]
        anchors[:, 1] = radius * np.sin(angles) + anchor[1]
        anchors[:, 2] = anchor[2]

        return anchors

class NemaIQAuto(NemaIQ):
    def __init__(self, image):
        """ Automatically find hot spheres and position all phantom elements.
        """
        debug = self.logger.level == self.LEVEL_DEBUG

        anchor, angle, rotation = NemaIQAuto._orient_spheres_NemaIQ(image,
                                                                    debug = debug)

        # If the center of gravity of the entire phantom lies below the sphere
        # plane the orientation of the phantom is UP otherwise DOWN.

        if geometry.centroid(image)[2] < anchor[2]:
            orientation = self.UP
        else:
            orientation = self.DOWN
            # if phantom is flipped upside down, rotation inverses as well
            rotation = not rotation


        super(NemaIQAuto, self).__init__(image=image, anchor=anchor,
             sphere_angle = angle, rotation = rotation, up_down = orientation)

    @staticmethod
    def _orient_spheres_NemaIQ(image, nspheres = 4, threshold = 0.5, debug = False):
        """ Automatically orient Nema IQ phantom in image.
        nspheres defines number of spheres that are visible above background.
        treshold is set and used to delineate spheres from background.

        returns center, angle and orientation

        center is the geometric center around which the spheres are positioned.
        angle the rotation around z-axis for the spheres.
        orientation is clockwise (0) or anti clockwise (1)."""

        im = copy(image) # leave original image as is

        anchors = [] # collect the physical point location

        radius_sphere = 40 # use as temporary mask for sphere

        mean_values = [] # collect the summed voxel value for each sphere

        # finde spheres above background
        for i in range(0, nspheres):

           # global maximum value is inside a sphere
           max_value = geometry.sitk_max(im)

           # get voxel location for maximum value
           locI = geometry.voxel_index_for_value(im, value = max_value)[0]

           # get physical location for maximum value
           loc = im.TransformIndexToPhysicalPoint(locI)

           # create a mask to mask out this sphere to find next one
           sphere = Sphere(im, anchor = loc, radius = radius_sphere)

           # centroid of the masked image is a the estimate of sphere location
           anchors += [geometry.centroid(im, mask = sphere.mask)]

           # mean value should be indicative of sphere size
           mean_values += [geometry.sitk_mean(geometry.apply_mask(sphere.image,
                                                                sphere.mask))]

           # remove sphere from image in order to find next maximum in sphere
           im = geometry.apply_mask(im, (1-sphere.mask))

        # make numpy
        anchors = np.array(anchors)
        mean_values = np.array(mean_values)

        # sort anchors by region size
        index = np.argsort(-1* mean_values)
        mean_values = mean_values[index]
        anchors = anchors[index]

        # fit circle to points to find the center along which they are positioned
        x = anchors[:, 0]
        y = anchors[:, 1]

        # find center and radius od a circle trough the center of all spheres
        xc, yc, radius = geometry.fit_circle(x, y)
        zc = np.mean(anchors[:, 2])
        center = (xc, yc, zc)

        if debug:
            from matplotlib import pyplot as plt
            import SimpleITK as sitk
            plt.figure()
            anchors = [image.TransformPhysicalPointToContinuousIndex(ai)
                       for ai in anchors]

            anchorsI = np.array(anchors)


            zI = int(np.round(np.mean(anchorsI[:,2])))
            plt.imshow(sitk.GetArrayFromImage(image)[zI,:, :])
            for pI, color in zip(anchorsI, ['r', 'g', 'b', 'y']):
                plt.plot(pI[0], pI[1], color + 'd')

            centerI = image.TransformPhysicalPointToContinuousIndex((xc, yc, 0))
            plt.plot(*centerI[0:2], 'ro')

        # find polar angles to obtain clockwise or anti clockwise orientation
        theta = np.arctan2(y-yc, x-xc)

        if theta[1] > theta[0]:
            orientation = NemaIQ.ANTI_CLOCKWISE
        else:
            orientation = NemaIQ.CLOCKWISE

        return center, theta[0], orientation

class GrowRegion(Phantom):
    """ Phantom element that creates a mask by a region growing algorithm. """

    def __init__(self, image, seed = (0,0,0), threshold = 0.5, background = 0):
        """ Specify a threshold for the region growing algorithm and the seed
        point. """
        super(GrowRegion, self).__init__(image, coordinates_phantom = seed)

        self.anchor = seed
        self.threshold = threshold
        self._background = background

    @property
    def threshold(self):
        """
        Set treshold for the region growing algorithm as fraction of max value.
        """
        return self._threshold

    @threshold.setter
    def threshold(self, threshold):
        self._threshold = threshold
        self.reset()

    @property
    def background(self):
        """
        Set the background activity that will be added to the max voxel value,
        before tresholding, default = 0.
        """
        return self._background

    @background.setter
    def background(self, value):
        self._background = value
        self.reset()

    @property
    def seed(self):
        """ Seed point in physical units [mm] """
        return self.anchor

    @seed.setter
    def seed(self, seed):
        self.anchor = seed
        self.reset()

    @property
    def seedI(self):
        """ Seed point in voxel indices [i, j, k] """
        return self.transform_world_to_index(self.seed)

    @seedI.setter
    def seedI(self, value):
        self.seed = self.transform_index_to_world(value)

    def _make_sitk_mask(self):
        #seed = self.coordinates_index
        dbg_msg = 'Growing region from seed {0} with treshold {1} and background {2}'
        self.logger.debug(dbg_msg.format(self.seedI, self.threshold, self.background))
        mask = geometry.region_growing(self.image, seed = self.seedI,
                                       threshold = self.threshold,
                                       background = self.background)

        return mask

    def analyze(self):
        result = super(GrowRegion, self).analyze()
        result['threshold'] = self.threshold
        return result


class GrowRegionFromSphere(GrowRegion):
    """ Extension to the GrowRegion class, seed is found by looking within a
    sphere with a certain radius around the defined initial seed. The new seed is
    used for region growing. """

    def __init__(self, image = None, anchor = (0,0,0),
                 threshold = 0.5, radius = 10, background = 0):



        super(GrowRegionFromSphere, self).__init__(image = image,
             seed = anchor, threshold = threshold)

        self._sphere = None
        self._radius = radius
        self._seedI = None
        self.background = background

    @property
    def radius(self):
        """
        Defines the radius of the spherical search region for max voxel value.
        """
        return self._radius

    @radius.setter
    def radius(self, value):
        self._radius = value
        self.reset()

    @property
    def sphere(self):
        """ Returns a phantom sphere. Within this sphere the hottest voxel will
        be used as seed point. """
        if self._sphere is None:
            sphere = Sphere(image = self.image,
                            radius = self.radius,
                            anchor = self.anchor)
            self._sphere = sphere
        return self._sphere

    @sphere.setter
    def sphere(self, value):
        self._sphere = value

    @property
    def seed(self):
        return self.transform_index_to_world(self.seedI)

    @property
    def seedI(self):
        if self._seedI is None:
            max_value = geometry.sitk_max(self.image, self.sphere.mask)
            masked_im = geometry.apply_mask(self.image, self.sphere.mask)
            self._seedI = geometry.voxel_index_for_value(masked_im, max_value)[0]
        return self._seedI

    def reset(self):
        super(GrowRegionFromSphere, self).reset()
        self.sphere = None
        self._seedI = None






