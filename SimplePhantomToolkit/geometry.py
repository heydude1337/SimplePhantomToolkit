# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 12:59:21 2016

@author: Marcel
"""
import numpy as np
import SimpleITK as sitk
from scipy import linalg

def grid_from_image(image):
    """ Similar to numpy.meshgrid using sitk. Grids will be in world (physical)
    space. """
    imsize = image.GetSize()
    spacing = image.GetSpacing()
    origin = image.GetOrigin()
    direction = image.GetDirection()
    grid = sitk.PhysicalPointSource(sitk.sitkVectorFloat64, imsize,
                                    origin, spacing,  direction)

    dim =  image.GetDimension()
    grid = [sitk.VectorIndexSelectionCast(grid, i) \
                         for i in range(0, dim)]

    for gi in grid:
        gi.CopyInformation(image)
    return grid

def rotate_image(image, center, angle, interpolator = sitk.sitkNearestNeighbor):
    """ Rotates an image around center the given angle. For 3D images angle
    must have three components defining the euler angles. """

    if all([ai==0 for ai in angle]): return image


    if image.GetDimension() == 2:

      transform = sitk.Euler2DTransform()
      transform.SetAngle(angle)
      transform.SetTranslation(center)
    elif image.GetDimension() == 3:

      parameters = (*angle,0,0,0)
      transform = sitk.Euler3DTransform()
      transform.SetParameters(parameters)
      transform.SetCenter(center)

    image = sitk.Resample(image, transform, interpolator)
    return image


# =============================================================================
# Geometry functions that work on both sitk an numpy
# =============================================================================
def dot(v, w):
    """ Dot product between two vectors. """
    if sitk.Image in (type(v[0]), type(w[0])):
        r = 0
        for vi, wi in zip(v, w):
            r += vi * wi
    else:
        r = v.dot(w)
    return r

def norm(v):
    """ norm of vector """
    if isinstance(v[0], sitk.Image):
        return sitk.Sqrt(sum([vi**2 for vi in v]))
    else:
        return float(linalg.norm(v))

def normalize(v):
    """ Normalize vector """
    if isinstance(v, sitk.Image):
        m = norm(v)
        nv = [vi/m for vi in v]
    else:
        nv = v / norm(v)
    return nv

def cross(u, v):
    """ Cross product between two vectors. """
    if isinstance(v[0], sitk.Image) or isinstance(u[0], sitk.Image):
        s1 = (u[1]*v[2] - u[2] * v[1])
        s2 = (u[2]*v[0] - u[0] * v[2])
        s3 = (u[0]*v[1] - u[1] * v[0])
        w = (s1, s2, s3)
    else:
        w = np.cross(u,v)
    return w

# =============================================================================
#  Geometric masks
# =============================================================================

def sitk_cylinder(image, axis=((10,10,10),(90, 90, 90)), radius = 20):
    """ Make a cylinder wiht a center, radius and heiht in a matrix of size
       imsize"""

    p = np.array(axis[0])
    q = np.array(axis[1])

    n = p - q

    # make image grid
    x, y, z = grid_from_image(image)

    u = (float(p[0]) -x, float(p[1])-y, float(p[2])-z)

    dxyz = cross(n.tolist(), u)

    d = sitk.Sqrt(dxyz[0]**2 + dxyz[1]**2 + dxyz[2]**2) / norm(n) < radius

    side1 = dot(n.tolist(), (x-float(p[0]), y-float(p[1]), z-float(p[2]))) <= 0
    side2 = dot(n.tolist(), (x-float(q[0]), y-float(q[1]), z-float(q[2]))) >= 0

    mask = d * side1 * side2
    return mask

def sitk_circle(image, center = (16, 16), radius = 8):
    """ Make a cirlcle at a center with specified radius. """
    gridx, gridy = grid_from_image(image)
    mask = ((gridx - center[0])**2 + (gridy - center[1])**2) <= (radius**2)
    return mask

def sitk_sphere(image = None, center = (16, 16, 16), radius = 1):
    """ Make a sphere at a center with specified radius. """

    gridx, gridy, gridz = grid_from_image(image)

    mask = sitk.Image(image.GetSize(), sitk.sitkInt8)
    mask.CopyInformation(image)

    mask = ((gridx - center[0]) ** 2) +  \
           ((gridy - center[1]) ** 2) +  \
           ((gridz - center[2]) ** 2) < (radius ** 2)

    return mask


def region_growing(image, seed = None, threshold = 0.5, background = 0):
    """ Use a region growing algorithm to obtain a binary mask. Seed is the
    seedpoint and threshold the fraction of the voxel value of the seed that
    defines the lowest value to be included. """

    upper = image[seed] + background
    lower = threshold * upper


    f = sitk.ConnectedThresholdImageFilter()
    f.SetLower(lower)
    f.SetUpper(upper)
    f.SetSeed(seed)

    mask = f.Execute(image)

    return mask

def apply_mask(image, mask):
    return image * sitk.Cast(mask, image.GetPixelIDValue())


def centroid(image, mask = None):
  """ Calculate the centroid (center of gravity of an image. """
  image = sitk.Cast(image, sitk.sitkFloat64)
  if mask is not None:
      image = apply_mask(image, mask)

  grid = grid_from_image(image)

  w = [g * image for g in grid]

  sum_im = sitk_sum(image)

  centroidI = []
  for wi in w:
   wi_sum = sitk_sum(wi)
   centroidI+=[wi_sum/sum_im]

  return centroidI

def sitk_min(image, mask = None):
    sitk_filter = _statistics_filter(image, mask = mask)
    if mask is None:
        value = sitk_filter.GetMinimum()
    if mask is not None:
        value = sitk_filter.GetMinimum(1)
    return value

def sitk_max(image, mask = None):
    sitk_filter = _statistics_filter(image, mask = mask)
    if mask is None:
        value = sitk_filter.GetMaximum()
    if mask is not None:
        value = sitk_filter.GetMaximum(1)
    return value

def sitk_mean(image, mask = None):
    sitk_filter = _statistics_filter(image, mask = mask)
    if mask is None:
        value = sitk_filter.GetMean()
    if mask is not None:
        value = sitk_filter.GetMean(1)
    return value

def sitk_std(image, mask = None):
    sitk_filter = _statistics_filter(image, mask = mask)
    if mask is None:
        value = sitk_filter.GetSigma()
    if mask is not None:
        value = sitk_filter.GetSigma(1)
    return value

def sitk_sum(image, mask = None):
    sitk_filter = _statistics_filter(image, mask = mask)
    if mask is None:
        value = sitk_filter.GetSum()
    if mask is not None:
        value = sitk_filter.GetSum(1)
    return value


def _statistics_filter(image, mask = None):

    if mask is None:
        sitk_filter = sitk.StatisticsImageFilter()
        sitk_filter.Execute(image)

    if mask is not None:
        sitk_filter = sitk.LabelStatisticsImageFilter()
        sitk_filter.Execute(image, mask)

    return sitk_filter

def voxel_index_for_value(image, value = 1):

    im_array = sitk.GetArrayFromImage(image)

    if value == 'min':
        value = np.min(im_array)
    if value == 'max':
        value = np.max(im_array)

    index = np.where(im_array == value)

    z, y, x  = ([*index[0]], [*index[1]], [*index[2]])

    x = [int(xi) for xi in x]
    y = [int(yi) for yi in y]
    z = [int(zi) for zi in z]

    return tuple(zip(x, y, z))



def fit_circle(x, y):
    # source: http://scipy-cookbook.readthedocs.io/items/Least_Squares_Circle.html
    #! python
    # == METHOD 1 ==
    # method_1 = 'algebraic'

    # coordinates of the barycenter
    x_m = np.mean(x)
    y_m = np.mean(y)

    # calculation of the reduced coordinates
    u = x - x_m
    v = y - y_m

    # linear system defining the center (uc, vc) in reduced coordinates:
    #    Suu * uc +  Suv * vc = (Suuu + Suvv)/2
    #    Suv * uc +  Svv * vc = (Suuv + Svvv)/2
    Suv  = np.sum(u*v)
    Suu  = np.sum(u**2)
    Svv  = np.sum(v**2)
    Suuv = np.sum(u**2 * v)
    Suvv = np.sum(u * v**2)
    Suuu = np.sum(u**3)
    Svvv = np.sum(v**3)

    # Solving the linear system
    A = np.array([ [ Suu, Suv ], [Suv, Svv]])
    B = np.array([ Suuu + Suvv, Svvv + Suuv ])/2.0
    uc, vc = linalg.solve(A, B)

    xc_1 = x_m + uc
    yc_1 = y_m + vc

    # Calcul des distances au centre (xc_1, yc_1)
    Ri_1     = np.sqrt((x-xc_1)**2 + (y-yc_1)**2)
    R_1      = np.mean(Ri_1)
    # residu_1 = sum((Ri_1-R_1)**2)

    return xc_1, yc_1, R_1

def rotate(point, angle, center = None):
    fcn = rotate_fcn(angle, center)
    if type(point[0]) in (list, tuple):
        rpoint = list(map(fcn, point))
    else:
        rpoint = fcn(point)
    return rpoint

def rotate_fcn(angle, center = None):
    if type(angle) in (list, tuple):
        T = sitk.Euler3DTransform()
        angle_fcn = T.SetRotation
        if center is None:
            center = (0,0,0)
    else:
        T = sitk.Euler2DTransform()
        angle_fcn = T.SetAngle
        if center is None:
            center = (0,0)

    angle_fcn(angle)
    T.SetCenter(center)
    return T.TransformPoint


