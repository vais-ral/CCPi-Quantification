'''
Unit tests for instrument classes 
'''

import unittest
from ccpi.quantification import LabelQuantificationUShort
import numpy as np
import math
from tifffile import TiffFile    
class TestLabelQuantification(unittest.TestCase):    

    def test_compute_label_quantification(self):
        img = TiffFile('test/FoamData.tif')        
        data = img.asarray()
        voxel_size = np.ones(3,dtype=float)
        origin = np.zeros(3,dtype=int)
        lqs = LabelQuantificationUShort(data, origin, voxel_size, float(np.amin(data)), float(np.amax(data)), 100.0)
        lqs.compute()
        output = lqs.getOutput()
        print(output)
        img.close()
    
     
        
if __name__ == "__main__":
    unittest.main()