
import ij.IJ;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author sn65
 */
public class CCPiImageJUserInterface extends CCPiUserApplicationInterface{
    @Override
  public void LogMessage(String message) {
      IJ.getTextPanel().append(message);
  }

    @Override
  public void SetStatusMessage(String message) {
      IJ.showStatus(message);
  }

    @Override
  public void SetProgressValue(float value) {
      IJ.showProgress(value);
  }

    @Override
  public boolean isCancel() {
      return false;
  }   
}
