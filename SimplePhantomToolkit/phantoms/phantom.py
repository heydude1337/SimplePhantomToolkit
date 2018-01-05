# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 11:07:21 2016

@author: Marcel
"""

import SimpleITK as sitk
import pandas as pd
import numpy as np


from SimplePhantomToolkit import Logger, PhantomPlot, geometry

from SimplePhantomToolkit.analysis import sitk_mask_analysis



class CoordinateSystem(Logger):
    """
    Creates a coordinate system with appropiate transformation to transform
    physical coordinates, phantom coordinates and index coordinates to each other.

    index coordinates:   Location specified as voxel location i, j, (k)
    world coordinates:   Coordinates in physical units (mm) in the dicom image
    phantom coordinates: Coordinates in physical units relative to the origin,
                         describing the phantom.
    """


    _LOG_LEVEL = Logger.LEVEL_INFO

    def __init__(self, angle = None, anchor = None, coordinates_phantom=None):
        """
        Create coordinate system, angle specifies the rotation angle between world
        and phantom coordinates. Translation is specified by anchor. The origin of
        the phantom is moved to the anchor to obtain the world coordinates of the
        phantom."""

        self.angle = angle
        self.anchor = anchor
        self.coordinates_phantom = coordinates_phantom
        self._p_to_w_transform = None

    @property
    def anchor(self):
        """ Origin in world coordinates in mm. Anchor and angle define a rigid
        transform that is applied to phantom coordinates to obtain the world
        coordinates."""
        return self._anchor

    @anchor.setter
    def anchor(self, anchor):
        self._set_anchor(anchor)

    def _set_anchor(self, anchor):
        anchor = np.asarray(anchor, dtype = float)
        if not(hasattr(self, 'anchor')):
            self._anchor = anchor
            self.reset()
        if not(np.all(self.anchor == anchor)):
            self._anchor = anchor
            self.reset()

    @property
    def angle(self):
        """ Euler angle(s) in radians to define rotation of phantom
        coordinates. Anchor and angle define a rigid transform that is applied
        to phantom coordinates to obtain the world coordinates."""
        return self._angle

    @angle.setter
    def angle(self, angle):
          self._set_angle(angle)

    def _set_angle(self, angle):
        angle = np.asarray(angle, dtype = float)
        if not(hasattr(self, '_angle')):
            self._angle = angle


        if not(np.all(angle == self.angle)):
            self._angle = angle
            self.logger.debug('Angle set at %s', angle)
            self.reset()

    def move(self, vector):
        """ Translate phantom in world."""
        self.logger.debug('Moving %s mm', str(vector))
        self.anchor += vector

    def rotate(self, angle = None):
        """ Rotate phantom in world. """
        self.logger.debug('Rotating %s radians', str(angle))
        self.angle += angle

    def reset(self):
        """ Reset calculated values. Call this function to force recalculation.
        """
        if self._p_to_w_transform is not None:
            self.logger.debug('3D transform reset')
            self._p_to_w_transform = None

        self.logger.debug('Phantom has moved!')

# =============================================================================
# Coordinates for phantom, world and index
# =============================================================================

    @property
    def coordinates_phantom(self):
        """ Returns coordinate(s) of the phantom in mm relative to the anchor.
        """
        return self._coordinates_phantom

    @coordinates_phantom.setter
    def coordinates_phantom(self, coordinates):
        coordinates = np.asarray(coordinates, dtype=float)
        self._coordinates_phantom = coordinates
        self.reset()

    @property
    def coordinates_world(self):
        """ Returns coordinate(s) of phantom in mmm in world. """
        T = self.transform_phantom_to_world
        return T(self.coordinates_phantom)

    @property
    def coordinates_continuous_index(self):
        """ Returns coordinate(s) of phantom in image space continuous indices.
        """
        T = self.transform_phantom_to_continuous_index
        return  T(self.coordinates_phantom.tolist())

    @property
    def coordinates_index(self):
        """ Returns coordinate(s) of phantom in image space continuous indices.
        """
        return  np.round(self.coordinates_continuous_index)


# =============================================================================
# Transformation between each coordinate system
# =============================================================================

    @property
    def transform_phantom_to_world(self):
        """
        Rigid tranformation function to transform phantom to world coordinates.
        """
        T = self._phantom_to_world_transform.TransformPoint
        return  lambda p: self._apply_transform(T, p)

    @property
    def transform_world_to_phantom(self):
        """
        Rigid Tranformation function to transform world to phantom coordinates.
        """
        T = self._phantom_to_world_transform.GetInverse().TransformPoint
        return  lambda p: self._apply_transform(T, p)

    @property
    def transform_index_to_world(self):
        """
        Tranformation function to transform index to world coordinates.
        """
        T = self.image.TransformContinuousIndexToPhysicalPoint
        return  lambda p: self._apply_transform(T, p)

    @property
    def transform_world_to_continuous_index(self):
        """
        Tranformation function to transform world to continuous index coordinates.
        """
        T = self.image.TransformPhysicalPointToContinuousIndex
        return lambda p: self._apply_transform(T, p)

    # derived transformations
    @property
    def transform_phantom_to_continuous_index(self):
        """
        Tranformation function to transform phantom to continuous index coordinates.
        """
        T1 = self.transform_phantom_to_world

        T2 = self.transform_world_to_continuous_index

        return  lambda p: T2(T1(p))

    @property
    def transform_phantom_to_index(self):
        """
        Tranformation function to transform phantom to index coordinates.
        """
        T1 = self.phantom_to_continuous_index
        T2 = np.round

        return  lambda p: T2(T1(p)).astype(int).tolist()

    @property
    def transform_world_to_index(self):
        """
        Tranformation function to transform world to index coordinates.
        """
        T1 = self.transform_world_to_continuous_index
        T2 = np.round
        return  lambda p: T2(T1(p)).astype(int).tolist()

    @property
    def transform_index_to_phantom(self):
        """
        Tranformation function to transform index to phantom coordinates.
        """
        T1 = self.transform_index_to_world
        T2 = self.transform_world_to_phantom

        return lambda p: T2(T1(p))

# =============================================================================
# Calculate and store rigid body between phantom and world
# =============================================================================
    @property
    def _phantom_to_world_transform(self):
        # retrieve stored transform. If not stored calculate and store
        # note _p_to_w_transform is set to None when reset is called
        if self._p_to_w_transform  is None:
            self._p_to_w_transform = self._get_phantom_to_world_transform
        return self._p_to_w_transform

    @property
    def _get_phantom_to_world_transform(self):
        # calculate rigid body transform between phantom and world

        anchor = self.anchor.tolist()

        angle = self.angle.tolist()

        if self.image.GetDimension() == 2:
            self.logger.debug('Setting 2D transform')
            transform = sitk.Euler2DTransform()
            transform.SetAngle(angle)
            transform.SetTranslation(anchor)
        elif self.image.GetDimension() == 3:
            self.logger.debug('Setting 3D transform')
            parameters = (*angle, *anchor)
            transform = sitk.Euler3DTransform()
            transform.SetParameters(parameters)
        return transform

# =============================================================================
# Help functions
# =============================================================================

    @staticmethod
    def _apply_transform(transform, coordinates):
        # Apply a transfor to coordinate(s). SimpleITK does not accept numpy
        # scalars or arrays as input. This is a work around to apply a SimpleITK
        # transfor on numpy coordinates

        coordinates = np.asarray(coordinates, dtype=float)

        if coordinates.ndim == 1:
            transformed =  np.asarray(transform(coordinates.tolist()), dtype=float)
        elif coordinates.ndim == 2:
            # recursive call
            transformed = []
            for ci in coordinates:
                transformed += [CoordinateSystem._apply_transform(transform, ci)]

            transformed = np.stack(transformed)
        else:
            error_msg = ('Coordinates must be 1 or 2 dimensional. '
                         '{0} dimension found.')
            raise IndexError(error_msg.format(coordinates.ndim))
        return transformed


class Phantom(CoordinateSystem):
    """ Base class from which all phantoms are derived. Phantoms are defined by
        a point, or set of points.
    """

    def __init__(self, image = None,
                       coordinates_phantom=None,
                       anchor=None,
                       angle=None):

        """ Create phantom object on a specified image. Phantom location is
        defined by (a set of) coodinate(s). Anchor defines to point relative to
        which the coordinate(s) is/are defined. Angle defines the relative
        rotation between the world and the phantom """

        self._image = image
        self._mask = None
        self._p_to_w_transform = None

        super(Phantom, self).__init__(angle=angle, anchor=anchor,
             coordinates_phantom=coordinates_phantom)

#==============================================================================
#   Override setters to set default values
#==============================================================================

    def _set_angle(self, angle):
        if angle is None:
            if self.image.GetDimension() == 2:
                angle = 0
            elif self.image.GetDimension() == 3:
                angle = (0, 0, 0)

        super(Phantom, self)._set_angle(angle)

    def _set_anchor(self, anchor):
        if anchor is None:
            anchor = [0, 0, 0]
            anchor =  anchor[0:self.image.GetDimension()]
        super(Phantom, self)._set_anchor(anchor)


#==============================================================================
#   Properties (Read/Write)
#==============================================================================

    @property
    def image(self):
        """ Returns/sets the image. Image is of type SimpleITK.Image. """
        return self._image

    @image.setter
    def image(self, image):
        if image is None:
            image = sitk.Image((10,10,10), sitk.sitkInt8)
        elif isinstance(image, np.ndarray):
            image = sitk.GetImageFromArray(image)

        if image is not self._image:
            self._image = image
            self.reset()


    @property
    def mask(self):
        """ Calculate and return a SimpleITK mask/label image for the phantom
        """
        if self._mask is None:
            self._mask = self._make_sitk_mask()
        return self._mask

    def _make_sitk_mask(self):
        error_msg = '_make_sitk_mask must be implemented by sublass of Phantom.'
        raise NotImplementedError(error_msg)


#==============================================================================
#   Methods (Calculations)
#==============================================================================
    def analyze(self):
        """ Obtain mean, max, min etc. for the image mask. """
        result = sitk_mask_analysis(self.image, label_mask=self.mask, ignore0=True)
        result = pd.DataFrame(result)
        result['location'] = str(self.coordinates_world)
        result['type'] = self.__class__.__name__
        return result

    def reset(self):
        super(Phantom, self).reset()
        if self._mask is not None:
            self.logger.debug('phantom mask reset')
            self._mask = None


#==============================================================================
#   Display functions
#==============================================================================

    def display(self):
        """ display using ImageJ """
        mask = sitk.Cast(self.mask, self.image.GetPixelIDValue())
        sitk.Show(mask)
        sitk.Show(self.image)
        sitk.Show(self.image * mask)

    def plot(self):
        """ Use matplotlib to plot the phantom. """
        return PhantomPlot(self).plot()


    def __str__(self):
        s = ('Phantom Type:\t {0} \nPhantom:\t {1} \nWorld Location:\t {2} '
             '\nRotation Angle\t {3}')
        return s.format(type(self).__name__,
                        str(self.coordinates_phantom),
                        self.anchor,
                        self.angle)


    def __repr__(self):
        return self.__str__()

    def __add__(self, rhs):
        if isinstance(rhs, Phantom):
            result = Phantoms([self, rhs])
        elif isinstance(rhs, Phantoms):
            result = rhs.append(self)
        else:
            raise ValueError('Cannot add type {0} to Phantoms'.format(type(rhs)))
        return result

    def __radd__(self, other):
        return self.__add__(other)


class Phantoms(list, Logger):
    _LOG_LEVEL = Phantom.LEVEL_DEBUG

    """ Combine multiple phantoms and generate single output """
    def __init__(self, *args):
        super(Phantoms, self).__init__(*args)

        # caching of calculations
        self._mask = None
        self._masks = None

    def __getitem__(self, key):
        # ensure a Phantoms object is returned for slices
        result = super(Phantoms, self).__getitem__(key)

        # for a sice a list is returned, convert to Phantoms
        if type(result) is list:
            result = Phantoms(result)
        return result

    def __add__(self, rhs):
        # Join multiple phantoms to a Phantoms object
        self.logger.debug('Adding type %s to Phantoms', str(type(rhs)))
        if isinstance(rhs, Phantom):
            self.append(rhs)
            result = self
        elif isinstance(rhs, Phantoms):
            result = super(Phantoms, self).__add__(rhs)
        else:
            raise ValueError('Cannot add type {0} to Phantoms'.format(type(rhs)))

        if type(result) is list:
            result = Phantoms(result)
        return result

    def __radd__(self, other):
        return self.__add__(other)

    @property
    def image(self):
        """ Returns/sets the image. Image is of type SimpleITK.Image """
        if len(self) > 0:
            return self[0].image
        else:
            return self._image

    @image.setter
    def image(self, image):
        self._image = image
        for phantom in self:
            phantom.image = image

    @property
    def mask(self):
        """ Calculate and return a SimpleITK mask/label image for the phantom.
        """
        if self._mask is None:
            self._mask = sum(self.masks)
        return self._mask

    @property
    def masks(self):
        """ Calculate and return a SimpleITK mask/label image for each phantom
        element.
        """
        if self._masks is None:
            self.logger.info('Expensive mask calculation')
            self._masks = [phantom.mask for phantom in self]
        return self._masks


    def analyze(self, masks = None, index = []):
        """ Obtain mean, max, min etc. for the image mask. """
        analysis = [phantom.analyze() for phantom in self]
        analysis = pd.concat(analysis)
        return analysis


    def reset(self):
        self.logger.debug('Phantoms reset')
        self._mask = None
        self._masks = None
        for element in self:
            element.reset()


    def plot(self):
        """ Use matplotlib to plot the phantom. """
        return PhantomPlot(self).plot()

    def move(self, vector):
        """ Translate phantom in world."""
        self.logger.debug('Moving %s mm', str(vector))
        for phantom in self:
            phantom.anchor += np.asarray(vector)
        self.reset()


    def rotate(self, angle, point):
        """ rotate each phantom element around a point """
        self.logger.debug('Rotating %s around %s', str(angle), str(point))
        for phantom in self:
            phantom.rotate(angle)
            phantom.anchor = geometry.rotate(phantom.anchor, angle, point)
        self.reset()


