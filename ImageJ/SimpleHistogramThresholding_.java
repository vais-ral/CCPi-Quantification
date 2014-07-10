
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import java.util.HashMap;

/*
 * This implements SimpleHistogramThresholding_
 */

/**
 *
 * @author Mr. Srikanth Nagella
 */
public class SimpleHistogramThresholding_ implements PlugInFilter{
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
        int[] dims = new int[3];
        dims[0] = inputImage.getWidth();
        dims[1] = inputImage.getHeight();
        dims[2] = imagePlus.getNSlices();  
        float[] voxelSize = new float[3];
        voxelSize[0] = 1;
        voxelSize[1] = 1;
        voxelSize[2] = 1;
        float[] origin = new float[3];
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;        
        float minIntensity = 0;
        float maxIntensity = Short.MAX_VALUE;
        short[] data =new short[dims[0]*dims[1]*dims[2]];//readRaw("Data128.raw",dims);
        for(int sliceNumber=0;sliceNumber < inputImage.getSize();sliceNumber++){
            byte[] inputImageArray = (byte[]) inputImage.getPixels(sliceNumber+1);            
            int sindex = sliceNumber*(inputImage.getHeight())*(inputImage.getWidth());
            for(int h=0; h< inputImage.getHeight(); h++){
                int hindex = h * (inputImage.getWidth());
                for(int w=0;w<inputImage.getWidth();w++){
                    data[sindex+hindex+w]= inputImageArray[hindex+w];
                }
            }
        }        
        CCPiImageDataUnsignedChar image = new CCPiImageDataUnsignedChar(data, dims, true);
        CCPiSimpleHistogramThresholdingITKImplUnsignedChar filter = new CCPiSimpleHistogramThresholdingITKImplUnsignedChar(image, dims, voxelSize, origin, minIntensity, maxIntensity);
        filter.Compute();
        short[] outputImageData =  filter.GetOutputImage().GetImage();
        ImagePlus outputImageIP =  IJ.createHyperStack("SimpleHistograThresholding", inputImage.getWidth(), inputImage.getHeight(),imagePlus.getNChannels(), imagePlus.getNSlices(), imagePlus.getNFrames(), imagePlus.getBitDepth());
        int count=0;
        for(int sliceNumber=0;sliceNumber < inputImage.getSize();sliceNumber++){
            ImageProcessor op = outputImageIP.getStack().getProcessor(sliceNumber+1);
            byte[] outputImageArray = (byte[]) op.getPixels();
            int sindex = sliceNumber*(inputImage.getHeight())*(inputImage.getWidth());
            for(int h=0; h< inputImage.getHeight(); h++){
                int hindex = h * (inputImage.getWidth());
                for(int w=0;w<inputImage.getWidth();w++){
                 //   data[sindex+hindex+w]= inputImageArray[hindex+w];
                     outputImageArray[hindex+w] = (byte) outputImageData[sindex+hindex+w];
                }
            }
        }        
        outputImageIP.show();
        outputImageIP.updateAndDraw();                
    }
    
    void showAbout() {
        IJ.showMessage("SimpleHistogramThresholding (CCPi)...", " Finds the two peaks in the Image and binary Thresholds half way between the two peaks.");
    }    
    
}
