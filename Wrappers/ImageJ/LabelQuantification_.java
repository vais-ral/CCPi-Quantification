
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import java.util.HashMap;

/*
 * This implements LabelQuantification_
 * based on the code given by Shen (MXIF)
 */

/**
 *
 * @author Mr. Srikanth Nagella
 */
public class LabelQuantification_ implements PlugInFilter{
    ImagePlus imagePlus;
    static {
        System.loadLibrary("CCPiImageJPlugin_Cpp");
    }
    
    @Override
    public int setup(String string, ImagePlus ip) {
        if (string.equals("about")) {
            showAbout();
            return DONE;
        }
        imagePlus = ip;
        return DOES_8G/*+DOES_16+DOES_32*/+STACK_REQUIRED;
    }    
    

    @Override
    public void run(ImageProcessor ip) { 
        if(imagePlus.getBitDepth()==8)
            Process(Short.class,imagePlus.getImageStack());
    }
    
    <T> void Process(Class<T> type,ImageStack inputImage)
    {
	
      ij.gui.GenericDialog gd = new ij.gui.GenericDialog("Label Quantification");
	  float minFeatureSize = 100;
	  float vSize = 1;
      gd.addNumericField("Minimum Feature Size: ", minFeatureSize, 0);
      gd.addNumericField("Voxel Size: ", vSize, 0);
      gd.showDialog();
      if (gd.wasCanceled()) return;
      minFeatureSize = (float)gd.getNextNumber();
      vSize = (float)gd.getNextNumber();	
        int[] dims = new int[3];
        dims[0] = inputImage.getWidth();
        dims[1] = inputImage.getHeight();
        dims[2] = imagePlus.getNSlices();  
        float[] voxelSize = new float[3];
        voxelSize[0] = vSize;
        voxelSize[1] = vSize;
        voxelSize[2] = vSize;
        float[] origin = new float[3];
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;        
		//TODO: Find the values
        float minLabelValue = (float)inputImage.getProcessor(1).getMin();
        float maxLabelValue = (float)inputImage.getProcessor(1).getMax();
        short[] data =new short[dims[0]*dims[1]*dims[2]];//readRaw("Data128.raw",dims);
        for(int sliceNumber=0;sliceNumber < inputImage.getSize();sliceNumber++){
            byte[] inputImageArray = (byte[]) inputImage.getPixels(sliceNumber+1);    
            int sindex = sliceNumber*(inputImage.getHeight())*(inputImage.getWidth());
			if(inputImage.getProcessor(sliceNumber+1).getMin() < minLabelValue)
				minLabelValue = (float)inputImage.getProcessor(sliceNumber+1).getMin();
			if(inputImage.getProcessor(sliceNumber+1).getMax() > maxLabelValue)
				maxLabelValue = (float)inputImage.getProcessor(sliceNumber+1).getMax();
            for(int h=0; h< inputImage.getHeight(); h++){
                int hindex = h * (inputImage.getWidth());
                for(int w=0;w<inputImage.getWidth();w++){
                    data[sindex+hindex+w]= inputImageArray[hindex+w];
                }
            }
        }        
		int vtkDataType = 3;// for unsigned char TODO: change this to accept more types
        CCPiImageDataUnsignedChar image = new CCPiImageDataUnsignedChar(data, dims, true);
		CCPiImageJUserInterface ui = new CCPiImageJUserInterface();
        CCPiLabelQuantificationITKImplUnsignedChar filter = new CCPiLabelQuantificationITKImplUnsignedChar(image, ui, origin, dims, voxelSize, minLabelValue, maxLabelValue, minFeatureSize, vtkDataType);
        filter.Compute();
		CCPiLabelQuantificationResult output = filter.GetOutput();
		ij.measure.ResultsTable resultsTable = new ij.measure.ResultsTable();
		StringVector columnHeaders = output.GetQuantityNames();
		IntList		 rowIndexes = output.GetLabelIndexes();
		for(int columnIdx = 0;columnIdx<columnHeaders.size();columnIdx++){
			for(int rowIdx =0;  rowIdx<rowIndexes.size(); rowIdx++) {
				double value = output.GetValue(columnHeaders.get(columnIdx), rowIndexes.get(rowIdx));
				resultsTable.setValue(columnHeaders.get(columnIdx), /*rowIndexes.get(rowIdx)*/ rowIdx, value);
			}
		}
		resultsTable.show("Label Quantification Results");
    }
    
    void showAbout() {
        IJ.showMessage("LabelQuantification (CCPi)...", " Generates foam quantities based on code provided by Shen (MXIF).");
    }    
    
}
