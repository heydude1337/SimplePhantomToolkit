import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
from SimplePhantomToolkit import Phantom, Phantoms, geometry

class SphereBase(Phantom):
    """ Sphere defined by two opposite points on the surface (point1, point2)
    """
    def __init__(self, image=None, points=((1,0,0),(-1,0,0)),
                 anchor=(0,0,0)):

        super(SphereBase, self).__init__(image = image,
                                     coordinates_phantom = points,
                                     anchor = anchor)

    def _make_sitk_mask(self):
        p1, p2 = self.coordinates_world

        center = 0.5 * (p1 + p2)
        radius = geometry.norm(p2 - p1) / 2

        mask = geometry.sitk_sphere(image=self.image,
                                    center=center,
                                    radius=radius)
        return mask


class Sphere(SphereBase):
    def __init__(self, image=None, radius=1, anchor = (0,0,0)):
        p1 = (0, 0, -radius)
        p2 = (0, 0, radius)

        super(Sphere, self).__init__(image=image,
                                      points = (p1, p2),
                                      anchor=anchor)


    @property
    def radius(self):
        p1, p2 = self.coordinates_phantom
        return geometry.norm(p2 - p1) / 2


    @radius.setter
    def radius(self, radius):
        p1 = [0, 0, -radius]
        p2 = [0, 0, radius]
        self.coordinates_phantom = (p1, p2)
        self.reset()


class CylinderBase(Phantom):
    """ Cylinder defined by three points. point1 and point2 are central the
      circularn edges of the cylinder. point3 is located at the boundary of
      the circle around point 1. The length of point3 - point1 equals the
      radius of the cylinder. """

    def __init__(self, image=None,
                     points=((0, 0, -1), (0,0, 1), (0,1, -1)),
                     anchor = (0, 0, 0)):



        super(CylinderBase, self).__init__(image=image,
                                            coordinates_phantom=points,
                                            anchor = anchor)


    def _make_sitk_mask(self):
        p1, p2, p3 = self.coordinates_world
        radius = geometry.norm(p3 - p1)
        axis = (p1, p2)

        self.logger.debug('making mask at {0} with radius {1}'.format(axis,radius))
        mask = geometry.sitk_cylinder(self.image, axis = axis, radius = radius)

        return mask

class Cylinder(CylinderBase):
    def __init__(self, image=None, axis = ((0,0,0), (0,0,1)),
                 radius = 1, anchor = (0,0,0)):
        msg = 'Creating cylinder with axis %s and radius %s at location %s'
        self.logger.debug(msg, str(axis), str(radius), str(anchor))
        axis = np.asarray(axis)

        points = self._points_for_axis_and_size(axis, radius)

        super(Cylinder, self).__init__(image=image, points = points,
              anchor = anchor)


    @staticmethod
    def _points_for_axis_and_size(axis, radius):
        # calculate a point on the edge of the circular side of the cylinder
        # (p3). The norm between p1 and p3 is then the radius.
        # find point 3 perpendicular to p1, p2 at distance p3 -p1

        p1, p2 = axis

        # define axis as vector
        v = p2 - p1

        # rule out both points are equal
        if geometry.norm(v) == 0:
            raise ValueError('Axis should be defined by two unique points')

        # random point that does not lie on the axis
        ok = False
        while not ok:
            q = np.random.rand(3)
            if geometry.dot(q,v) != 0:
                ok = True

        # vector from p1 to q
        w = q - p1

        # defines a plane on which the axis lies
        n = geometry.normalize(geometry.cross(v, w))

        p3 = p1 + (n * radius)

        return p1, p2, p3


    @property
    def radius(self):
        """ Get/Set the radius of the cylinder """
        p1, _, p3 = self.coordinates_phantom
        return geometry.norm(p3 -p1)

    @radius.setter
    def radius(self, radius):
        p1, _, p3 = self.coordinates_phantom

        n = np.normalize(p3-p1)
        p3 = p1 + (n * radius)
        self.coordinates_phantom[2] = p3
        self.reset()

    @property
    def axis(self):
        """ Get/Set the axis of the cylinder. Axis is defined by two points.
        """
        return self.coordinates_phantom[0:2]

    @axis.setter
    def axis(self, axis):
        axis = np.array(axis)
        self.coordinates_phantom[0:2] = axis
        self.reset()

if __name__ == "__main__":


    image = sitk.GetImageFromArray(np.random.rand(128,128, 128))


    image.SetOrigin((-64, -64, -64))
    cylinder = CylinderBase(image,
               anchor = (0, 0, 0),
               points = ((0, 0, -64), (0, 64, 64), (0, 10, -64)))

    plt.figure()
    cylinder.plot()

    cylinder = Cylinder(image,
               anchor = (0, 0, 0),
               axis = ((0, 0, -64), (0, 64, 64)),
               radius = 10)

    plt.figure()
    cylinder.plot()




    sphere = Sphere(image=image,
                    anchor = (0,0,0),
                    radius = 50)
    sphere.plot()

    N = 10
    centers = np.random.rand(N,3) * 128 - 64
    radius = np.random.rand(N) * 30

    spheres = Phantoms([Sphere(image=image, radius = ri,
                               anchor= ci) for ri, ci in zip(radius, centers)])
    spheres.plot()

    radius = np.random.rand(N) * 10
    axis = np.random.rand(N, 2, 3) * 128 - 64

    cylinders = Phantoms([Cylinder(image=image, radius = ri, axis = ai) for ai, ri in zip(axis, radius)])

    cylinders.plot()

    pass