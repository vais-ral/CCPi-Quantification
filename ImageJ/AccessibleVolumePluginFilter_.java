
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import javax.swing.JOptionPane;
/*
 * To change this template, choose Tools | Templates and open the template in
 * the editor.
 */
/**
 *
 * @author sn65
 */
public class AccessibleVolumePluginFilter_ implements PlugInFilter {

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
        return DOES_8G+STACK_REQUIRED;
    }

    @Override
    public void run(ImageProcessor ip) { 
        HashMap<Integer,String> imageIDTitleMap = GenerateImageIDAndTitleHashMap();
        CCPiAccessibleVolumeInputDialog accessibleVolumeInputDialog = new CCPiAccessibleVolumeInputDialog(null,true);
        String[] maskFileNames = (String[])imageIDTitleMap.values().toArray(new String[0]);   
        for(String maskFileName : maskFileNames)
            accessibleVolumeInputDialog.addMaskFileName(maskFileName);
        accessibleVolumeInputDialog.setVisible(true);           
        if(!accessibleVolumeInputDialog.isCancelButtonPressed())        
            Process(imagePlus.getImageStack(),WindowManager.getImage(accessibleVolumeInputDialog.getMaskFileName()).getImageStack(),accessibleVolumeInputDialog.getOutputFileName(),accessibleVolumeInputDialog.getMinSphereDiameter(),accessibleVolumeInputDialog.getMaxSphereDiameter(),accessibleVolumeInputDialog.getNumberOfSpheres(),accessibleVolumeInputDialog.getImageResolution());
    }

    void Process(ImageStack inputImage, ImageStack maskImage,String outputFilename,double minSphereDiameter,double maxSphereDiameter,int numberOfSpheres,double imageResolution)
    {
        int[] dims = new int[3];
        dims[0] = inputImage.getWidth();
        dims[1] = inputImage.getHeight();
        dims[2] = imagePlus.getNSlices();
        float[] voxelSize = new float[3];
        voxelSize[0] = (float)imageResolution;
        voxelSize[1] = (float)imageResolution;
        voxelSize[2] = (float)imageResolution;
        float[] origin = new float[3];
        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;
        
        short[] data =new short[dims[0]*dims[1]*dims[2]];//readRaw("Data128.raw",dims);
        short[] maskData = new short[dims[0]*dims[1]*dims[2]];//readRaw("Mask128.raw",dims);    
        short[] outputImageData = new short[dims[0]*dims[1]*dims[2]];

        for(int sliceNumber=0;sliceNumber < inputImage.getSize();sliceNumber++){
            byte[] inputImageArray = (byte[]) inputImage.getPixels(sliceNumber+1);
            byte[] inputMaskImageArray = (byte[]) maskImage.getPixels(sliceNumber+1);            
            int sindex = sliceNumber*(inputImage.getHeight())*(inputImage.getWidth());
            for(int h=0; h< inputImage.getHeight(); h++){
                int hindex = h * (inputImage.getWidth());
                for(int w=0;w<inputImage.getWidth();w++){
                    data[sindex+hindex+w]= inputImageArray[hindex+w];
                    maskData[sindex+hindex+w] = inputMaskImageArray[hindex+w];
                }
            }
        }
        
        CCPiImageJUserInterface ui = new CCPiImageJUserInterface();
		CCPiImageDataUnsignedChar imgData = new CCPiImageDataUnsignedChar(data, dims, true);
		CCPiImageDataUnsignedChar imgMaskData = new CCPiImageDataUnsignedChar(maskData, dims, true);
		CCPiImageDataUnsignedChar imgOutputData = new CCPiImageDataUnsignedChar(outputImageData, dims, true);
        CCPiAccessibleVolumeInputImages input = new CCPiAccessibleVolumeInputImages(dims,voxelSize,origin,imgData,imgMaskData);
        CCPiAccessibleVolumeITKImpl filter = new CCPiAccessibleVolumeITKImpl(input, ui, imgOutputData, (float)Math.log(minSphereDiameter), (float)Math.log(maxSphereDiameter), numberOfSpheres, (float)imageResolution);
        filter.Compute();

        Map<Double,Double> result = filter.GetAccessibleVolume();
        WriteAccessibleVolumeResultToCSVFile(outputFilename,result);
        outputImageData = filter.GetOutputImage().GetImage();
        ImagePlus outputImage =  IJ.createHyperStack("AccessibleVolume", inputImage.getWidth(), inputImage.getHeight(),imagePlus.getNChannels(), imagePlus.getNSlices(), imagePlus.getNFrames(), imagePlus.getBitDepth());
        int count=0;
        for(int sliceNumber=0;sliceNumber < inputImage.getSize();sliceNumber++){
            ImageProcessor op = outputImage.getStack().getProcessor(sliceNumber+1);
            byte[] outputImageArray = (byte[]) op.getPixels();
            int sindex = sliceNumber*(inputImage.getHeight())*(inputImage.getWidth());
            for(int h=0; h< inputImage.getHeight(); h++){
                int hindex = h * (inputImage.getWidth());
                for(int w=0;w<inputImage.getWidth();w++){
                 //   data[sindex+hindex+w]= inputImageArray[hindex+w];
                     outputImageArray[hindex+w] = (byte) outputImageData[sindex+hindex+w];
                     if(outputImageData[sindex+hindex+w]!=0) count++;
                }
            }
        }        
        outputImage.show();
        outputImage.updateAndDraw();
    }
    
    HashMap<Integer,String> GenerateImageIDAndTitleHashMap()
    {
         HashMap<Integer,String> result = new HashMap<Integer,String>();
         //Get Image ids
         int[] imageIDArray = WindowManager.getIDList();
         for(int index=0;index<imageIDArray.length;index++)
         {
             result.put(imageIDArray[index], WindowManager.getImage(imageIDArray[index]).getTitle());
         }
         return result;
    }
    
    void WriteAccessibleVolumeResultToCSVFile(String filename, Map<Double,Double> result)
    {
        try{
        File outputFile = new File(filename);        
        if(!outputFile.exists())
            outputFile.createNewFile();
        FileWriter fw = new FileWriter(outputFile.getAbsoluteFile());        
        BufferedWriter bfw = new BufferedWriter(fw);
		for(Map.Entry<Double,Double> entry: result.entrySet()){
            bfw.write(entry.getKey()+","+entry.getValue());
            bfw.write("\n");
        }
        bfw.close();
        }catch(IOException ex){        
        }
    }
    void showAbout() {
        IJ.showMessage("Accessible Volume (CCPi)...", " Uses the algorithm described in Appendix A of PhD thesis \"Non-destructive quantification of tissue scaffolds and augmentation implants using X-ray microtomography\" by Sheng Yue, Department of Materials, Imperial College London.");
    }
}
