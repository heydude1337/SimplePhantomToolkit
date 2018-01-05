from SimplePhantomToolkit import Phantom, Phantoms, geometry
import numpy as np


class CircleBase(Phantom):
    def __init__(self, image = None, point1=(-1,0 ), point2=(1,0), 
                 anchor = (0,0)):
        """ Define a 2-Dimensional circle phantom from two oppsite points on the
        circumference of the circle.

           image:   SimpleITK image or numpy array, 2D
           point1, point2: Oppisite points
           anchor: See Phantom Documentation
    
        """
        
        super(CircleBase, self).__init__(image = image,
             coordinates_phantom = (point1, point2),
             anchor = anchor)
  

    def _make_sitk_mask(self):
        
        
        c = self.coordinates_world
        radius = 0.5 * geometry.norm(c[1] - c[0]) 
        center = 0.5 * (c[0] + c[1])
        
        self.logger.debug('Making mask at point: {0} with radius: {1}'.format(center, radius))
        
            
        mask = geometry.sitk_circle(self.image, center = center, radius = radius)
        return mask


class Circle(CircleBase):
    def __init__(self, image = None, center = (0,0), radius = 1, 
                 anchor = (0,0)):
        
        center = np.array(center)
        point1 = center - [radius, 0]
        point2 = center + [radius, 0]
        
        super(Circle, self).__init__(image=image, point1=point1, 
             point2 = point2, anchor=anchor)
        
        
if __name__ == '__main__':
    SIZE = 500
    import SimpleITK as sitk
    image = sitk.Image((SIZE, SIZE), sitk.sitkFloat32)
    phantom1 = CircleBase(image = image, point1 = (0.1*SIZE,0.1*SIZE), 
                          point2 = (0.4 * SIZE, 0.4 * SIZE))
    phantom1.plot()

    phantom2 = Circle(image=image, center = (0.3 * SIZE, 0.3 * SIZE), radius = 0.15 * SIZE)
    phantom2.plot()
    
    centers = np.random.rand(10,2) * SIZE
    radius = np.random.rand(10) * 0.2 * SIZE
    
    phantom3 = Phantoms([Circle(image, ci, ri) for ci, ri in zip(centers, radius)])
    phantom3.plot()
    
    
    
    

   