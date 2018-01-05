# SimplePhantomToolkit
Analyze phantom experiments (PET, CT) using the SimpleITK framework

## Sample Code

The SimplePhantomToolkit works on SimpleIKT image objects. Let's assume 
'image' is a 3D SimpleIKT image.

    from SimplePhantomToolkit import Cylinder
    
    myphantom = Cylinder(image, 
                         axis = (0,0, -5), (0,0,5), 
                         radius = 5,
                         anchor = (64, 64, 64))
                         
    result = myphantom.analyze()
    

A Cylinder phantom is defined by the 'axis', the two points that are the center
points on the circular faces of the cylinder and by the 'radius'. The 
argument anchor defines the physical point in the image that corresponds to 
the origin of the phantom (0,0,0). In this example the cylinder is placed in 
the image with the axis starting at (64,64,59) and ending at (64,64,69).

The last line of code returns the mean, max, std etc. within the defined 
cylinder.


    myphantom.mask
    
Returns an SimpleITK image mask (voxel are 1 within the cylinder and 0 outside).

    myphantom.plot()
    
Plots the phantom using matplotlib on three orthogonal slices through the 
center of the phantom.

   myphantom.display()
   
Calls the SimpleITK Show method for the image and the mask. By default these
are displayed using ImageJ

   myphantom.rotate((theta_x, theta_y, theta_z))
   
Rotates the phantom using the Euler angles '(theta_x, theta_y, theta_z)'.
In the 2D case just call 'myphantom.rotate(theta)' with theta a scalar.
                        
## Available Phantoms

* 2D Phantoms
  * Circle
* 3D Phantoms
  * Sphere
  * Cylinder
  * PET Phantom (Siemens Biograph mCT, homogenous cylinder)
  * Nema Image Quality Phantom
  
See help for description of each phantom
                        