# -*- coding: utf-8 -*-

import SimpleITK as sitk
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pylab


from SimplePhantomToolkit import Logger, geometry


try:
    import seaborn as sns
    sns.set_style("whitegrid")
except:
    pass




params = {'legend.fontsize': 'x-large',
          'figure.figsize': (16, 9),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}

pylab.rcParams.update(params)



AXES_TITLES_3D = ['transversal [y-x]', 'axial 1 [z-x]', 'axial 2 [z-y]']

class PhantomPlot(Logger):


  ALPHA = 0.5
  COLORMAP_IMAGE = plt.cm.gray
  COLORMAP_OVERLAY = plt.cm.jet
  _LOG_LEVEL = Logger.LEVEL_ERROR

  def __init__(self, phantom = None, center = None):

    self.phantom = phantom
    self.image = phantom.image
    self.dimensions = self.phantom.image.GetDimension()

    self.logger.debug('Init plot for phantom {0} with dimensions {1}' \
                 .format(type(self.phantom).__name__, self.dimensions))

    self._center = None
    self._figure = None

  @property
  def figure(self):
      if self._figure is None:
          self._figure = plt.figure()
      return self._figure

  @property
  def center(self):
      if self._center is None:
          # take centroid of mask as default center
          centroid = geometry.centroid(self.phantom.mask)
          centroidI = self.image.TransformPhysicalPointToIndex(centroid)
          self._center = np.asarray(list(centroidI))
      return self._center

  @center.setter
  def center(self, value):

      self._center = np.asarray(value)


  def plot(self):
    if self.dimensions == 2:
      axes = self._plot_image2D()
      self._plot_overlay2D(axes)

    elif self.dimensions == 3:
       axes = self._plot_image3D()
       self._plot_overlay3D(axes)

    else:
      ERROR_MSG = 'Cannot init plot with image dimensions {0}'
      self.logger.error(ERROR_MSG.format(self.dimensions))
      return

    return axes

  def _plot_image2D(self):
    self.logger.debug('new 2D plot')
    image = sitk.GetArrayFromImage(self.phantom.image)
    #axes = plt.axes()
    axes = self.figure.add_subplot(111)
    axes.imshow(image, cmap = self.COLORMAP_IMAGE)

    return axes

  def _plot_overlay2D(self, axes):
    # plot mask as overlay with alpha blending
    self.logger.debug('adding mask overlay to axes: {0}'.format(axes))
    np_mask = sitk.GetArrayFromImage(self.phantom.mask)
    axes.imshow(np_mask>0, alpha = self.ALPHA, cmap=self.COLORMAP_OVERLAY)
    return axes


  def _plot_image3D(self, axtitles = AXES_TITLES_3D):
    # Plot saggital, coronal and transvers slices
    self.logger.debug('new 3D plot at {0}'.format(self.center))
    axes = make_axes(col_names = axtitles, figure = self.figure)

    image = self.phantom.image
    im1 = image[int(self.center[0]),:,:]
    im2 = image[:,int(self.center[1]),:]
    im3 = image[:,:,int(self.center[2])]

    for (title, ax), im in zip(axes.items(), (im1, im2, im3)):
      im=sitk.GetArrayFromImage(im)
      ax.imshow(im, cmap = self.COLORMAP_IMAGE)
      ax.set_title(title)
    return axes

  def _plot_overlay3D(self, axes):
    # plot mask as overlay with alpha blending
    self.logger.debug('new 3D overlay at {0}'.format(self.center))
    mask = self.phantom.mask

    m1 = mask[int(self.center[0]),:,:]
    m2 = mask[:,int(self.center[1]),:]
    m3 = mask[:,:,int(self.center[2])]

    for (title, ax), m in zip(axes.items(), (m1, m2, m3)):
      self.logger.debug('add overlay {0}'.format(title))
      m = sitk.GetArrayFromImage(m)
      ax.imshow(m > 0, alpha = self.ALPHA, cmap=self.COLORMAP_OVERLAY)




class NemaIQPlot(PhantomPlot):
    def __init__(self, phantom):

        centers = [sphere.seedI for sphere in phantom.spheres]
        self.centers = np.stack(centers)
        center = np.round(np.mean(self.centers, 0))
        self.NemaIQ = phantom

        super(NemaIQPlot, self).__init__(phantom = phantom.spheres,
                                         center = center)

    def plot(self, figure = None):
        axes = super(NemaIQPlot, self).plot()

        ax3 = axes[AXES_TITLES_3D[2]]
        ax3.plot(self.centers[:,0], self.centers[:, 1], 'x')

        label_xy = [center - self.center for center in self.centers]
        for (i, center), label_xy in zip(enumerate(self.centers), label_xy):
            # source: https://stackoverflow.com/questions/5147112/how-to-put-individual-tags-for-a-scatter-plot
            plt.annotate(str(i), xy=center[0:2], xytext=[label_xy[0]*3, -1 * label_xy[1]*3],
                         textcoords='offset points', ha='right', va='bottom',
                         bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                         arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))


def make_axes(col_names=None, row_names=None, figure = None):
  """ Create subplot and return a dict to access them with (col_name, row_name)
      as key, values are axes objects """
  axes = {}

  if figure is None:
      figure = plt.figure()

  if row_names is None:
    ncols = len(col_names)
    for c, col_name in enumerate(col_names):
      axes[col_name] = figure.add_subplot(1, ncols, c+1)
  elif col_names is None:
    nrows = len(row_names)
    for r, row_name in enumerate(row_names):
      axes[row_name] = figure.add_subplot(nrows, 1, r+1)
  else:
    ncols = len(col_names)
    nrows = len(row_names)

    for c, col_name in enumerate(col_names):
        for r, row_name in enumerate(row_names):
          n = r*ncols + c + 1
          axes[(col_name, row_name)] = figure.add_subplot(nrows, ncols, n)

  return axes




def show(phantom):
  """ Display phantom in 3D using imageJ """
  mask = phantom.mask
  sitk.Show(phantom.image)
  sitk.Show(mask)




