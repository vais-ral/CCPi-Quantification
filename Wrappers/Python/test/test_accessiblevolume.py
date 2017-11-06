'''
Unit tests for instrument classes 
'''

import unittest
from ccpi.quantification import AccessibleVolumeInput, AccessibleVolume
import numpy as np
import math
from tifffile import TiffFile    
class TestAccessibleVolume(unittest.TestCase):    
    def test_read_volumeinput(self):
        data = TiffFile('test/Data128.tif')
        mask = TiffFile('test/DataMask128.tif')
        voxel_size = np.ones(3,dtype=float)
        origin = np.zeros(3,dtype=int)
        input = AccessibleVolumeInput(voxel_size, origin, data.asarray(), mask.asarray())
        print(input.getScafoldVolume())
        print(input.getScafoldPorosity())
        print(input.getDimensions())
        print(input.getVoxelSize())
        print(input.getOrigin())        
        data.close()
        mask.close()
        
    def test_compute_accessible_volume(self):
        data = TiffFile('test/Data128.tif')
        mask = TiffFile('test/DataMask128.tif')
        voxel_size = np.ones(3,dtype=float)
        origin = np.zeros(3,dtype=int)
        input = AccessibleVolumeInput(voxel_size, origin, data.asarray(), mask.asarray())
        av = AccessibleVolume(input, math.log(80.0), math.log(600.0), 11, 9.0)
        av.compute()
        output = av.getAccessibleVolume()
        print(output)
        data.close()
        mask.close() 
    
     
        
if __name__ == "__main__":
    unittest.main()