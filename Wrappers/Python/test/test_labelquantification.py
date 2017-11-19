'''
Unit tests for instrument classes 
'''

import unittest
from ccpi.quantification.LabelQuantification import LabelQuantification
import numpy as np
import math
from tifffile import TiffFile    
class TestLabelQuantification(unittest.TestCase):    

    def test_compute_label_quantification(self):
        img = TiffFile('test/FoamData.tif')        
        data = img.asarray()
        voxel_size = np.ones(3,dtype=np.float32)
        origin = np.zeros(3,dtype=np.float32)
        [featureNames, featureValues] = LabelQuantification(data, origin, voxel_size, float(np.amin(data)), float(np.amax(data)), 100.0)
        self.assertEqual(featureNames[0],b'Bbox_diag')
        self.assertAlmostEqual(featureValues[0,0],320.52924, 5)
        self.assertEqual(featureNames[2],b'Label')
        self.assertAlmostEqual(featureValues[0,2],1.0)        
        self.assertAlmostEqual(featureValues[1,2],2.0)            
        img.close()
    
     
        
if __name__ == "__main__":
    unittest.main()