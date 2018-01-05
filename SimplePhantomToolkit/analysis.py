
import SimpleITK as sitk
#==============================================================================
# Helper functions for using the phantom object
#==============================================================================

def sitk_mask_analysis(image, label_mask, ignore0 = True):
  """ Obtain image statistics for voxels defined in label_mask, a list is
      returned with the results for each value of the mask (1, 2 3, 4).
      A value of 0 inside the mask is usually the background, which is
      omitted if ignore0 is True"""

  f = sitk.LabelStatisticsImageFilter()
  f.Execute(image, label_mask)
  return [sitk_result(f, label) for label in f.GetLabels() \
          if label > 0 or not(ignore0)]

def sitk_result(sitk_statistics, label):
  """Convert sitk statistics to a dictionary"""

  result = {}
  result['max']     =  sitk_statistics.GetMaximum(label)
  result['min']     =  sitk_statistics.GetMinimum(label)
  result['std']     =  sitk_statistics.GetSigma(label)
  result['mean']    =  sitk_statistics.GetMean(label)
  result['sum']     =  sitk_statistics.GetSum(label)
  result['count']   =  sitk_statistics.GetCount(label)
  result['label']   =  label
  return result