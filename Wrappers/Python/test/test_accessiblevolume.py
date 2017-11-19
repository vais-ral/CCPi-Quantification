'''
Unit tests for instrument classes 
'''

import unittest
from ccpi.quantification.AccessibleVolume import AccessibleVolume
import numpy as np
import math
from tifffile import TiffFile    
class TestAccessibleVolume(unittest.TestCase):    
        
    def test_compute_accessible_volume(self):
        data = TiffFile('test/Data128.tif')
        mask = TiffFile('test/DataMask128.tif')
        voxel_size = np.ones(3,dtype=np.float32)
        origin = np.zeros(3,dtype=np.float32)
        av = AccessibleVolume(data.asarray(), mask.asarray(), origin, voxel_size, math.log(80.0), math.log(600.0), 11, 9.0)
        self.assertAlmostEqual(av[0,0], 8.00000000e+01)        
        self.assertAlmostEqual(av[0,1], 9.69123542e-01)
        self.assertAlmostEqual(av[10,0], 6.00000061e+02)        
        self.assertAlmostEqual(av[10,1], 0)        
        data.close()
        mask.close() 
    
     
        
if __name__ == "__main__":
    unittest.main()