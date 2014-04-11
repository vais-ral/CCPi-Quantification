/** 
 * @file Header file for module to track particles between 3D frames.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date May 2013
 */

#ifndef CCPIPARTICLETRACKING_H
#define CCPIPARTICLETRACKING_H

#include <vector>

using namespace std;

#include <mclib/McHandle.h>  // Smart pointer template class

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortFloatTextN.h>
#include <hxcore/HxPortIntSlider.h>

#include <mclib/McDArray.h>

#include <Inventor/nodes/SoSeparator.h> // Open Inventor scene node

#include "CCPiParticleTracker.h"

#include "api.h"

 // Forward declarations
class HxUniformScalarField3;
class HxConnection;

/** 
 * A module that calculates the tracks between particles in successive
 * 3D volume frames
 * 
 * The algorithm is described in detail in 
 * Feature point tracking and trajectory analysis for video imaging in 
 * cell biology, I.F. Sbalzarini, P. Koumoutsakos, 
 * Journal of Structural Biology 151(2005) 182-195
 * 
 */
class CCPI_API CCPiParticleTracking : public HxCompModule
{
    HX_HEADER(CCPiParticleTracking);

  public:

/** 
 * Default constructor. Creates ports and sets default values
 */
    CCPiParticleTracking();
    ~CCPiParticleTracking();

    HxPortDoIt portAction;

    /** Port providing float text input field for displacement */
    HxPortFloatTextN portDisplacement;
    /** Port for total number of frames including the first frame.
     * Affects size of m_frames vector */
    HxPortIntSlider portNumFrames;

    /** Perform the calculation in this module. Called by Avizo. */
    virtual void compute();
    
    /** Update connection objects in module GUI according to number
     * of frames entered. */
    virtual void update();

  private:
  
    /// The object that tracks particles
    CCPiParticleTracker m_tracker;
    
    /// Open Inventor scene to draw track lines
    McHandle<SoSeparator> m_scene;
    
    /// Array of connections to the other frames
    McDArray<HxConnection*> m_frames;

    /** 
     * Add particles to the tracker. Particles are the centroids of objects
     * found in the 3D volume.
     * @param field The 3D volume data.
     * @param frameId The number of the frame.
     */
    void addParticles(HxUniformScalarField3 *field, size_t frameId);
    
    /**
     * Plot the tracks of particles
     */
    void plotTracks();
    
    /**
     * Create the connection objects for the other frames with particles
     * we wish to track. 
     * @param numConnections The number of extra frame connections
     */
    void createFrameConnections(const int& numConnections);
    
    /**
     * Colours for plotting track lines
     */
    static unsigned int m_colours[3];
    
    /**
     * Colours for plotting spheres representing particles
     */
    static float m_sphereColours[][3];
};

#endif // CCPIPARTICLETRACKING_H
