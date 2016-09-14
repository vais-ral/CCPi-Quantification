/** 
 * @file Implementation file for module to track particles between 3D frames.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date May 2013
 */

#include <string>

using namespace std;

#include <QApplication>

#include <hxcore/HxConnection.h>
#include <hxcore/HxMessage.h>
#include <hxfield/HxUniformScalarField3.h>

#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/nodes/SoVertexProperty.h>

#include "CCPiParticleTracking.h"

HX_INIT_CLASS(CCPiParticleTracking,HxCompModule)

// Colours for track lines
unsigned int CCPiParticleTracking::m_colours[3] = {0xFF0000FF, 0x00FF00FF, 0x0000FFFF};
// Colours for spheres representing particles
float CCPiParticleTracking::m_sphereColours[][3] = { {1,0,0}, {0,1,0}, {0,0,1} };

/** 
 * Default constructor. Creates ports and sets default values
 */
CCPiParticleTracking::CCPiParticleTracking() :
    HxCompModule(HxUniformScalarField3::getClassTypeId()),
    portAction(this,"action",QApplication::translate("CCPiParticleTracking", "Action")),
    portDisplacement(this,"displacement",QApplication::translate("CCPiParticleTracking", "Displacement"),1),
    portNumFrames(this,"numFrames",QApplication::translate("CCPiParticleTracking", "Number of frames"))
{
    portAction.setLabel(0,"DoIt");
    
    // Define default values for numeric ports
    portDisplacement.setMinMax(1.0,1000.0);
    portDisplacement.setValue(30.0);
    portNumFrames.setMinMax(0, 15);
    portNumFrames.setValue(2);
    
    // Create the connections for frames. We have one connection as 
    // a module object (to the first frame) so we create one less connection
    // than the number of frames.
    createFrameConnections( portNumFrames.getValue()-1 );

    // Create the Open Inventor scene
    m_scene = new SoSeparator;
    
}

/**
 * Destructor
 */
CCPiParticleTracking::~CCPiParticleTracking()
{
    // Hide the scene I created
    hideGeom(m_scene);
    
    // Delete dynamic connections.
    for ( mclong i = 0; i < m_frames.size(); ++i )
        delete m_frames[i];
    m_frames.resize( 0 );
}

/**
 * Create the connection objects for the other frames with particles
 * we wish to track. 
 * @param numConnections The number of extra frame connections
 */
void CCPiParticleTracking::createFrameConnections(const int& numberOfConnections)
{
    // Delete extra connections.
    while ( m_frames.size() > numberOfConnections )
    {
        delete m_frames.last();
        m_frames.removeLast();
    }

    // Create new connections.
    while ( m_frames.size() < numberOfConnections )
    {
        QString name = QString("frame%1").arg(m_frames.size());
        QString label = QApplication::translate("CCPiParticleTracking", "Frame %1").arg(m_frames.size()+2);
        HxConnection* connection = new HxConnection(this, name, label, 
          HxUniformScalarField3::getClassTypeId() );
        connection->show(); // To display the connection within the Properties Area.
        m_frames.append( connection );
    }    
}

/** Update connection objects in module GUI according to number
 * of frames entered. */
void CCPiParticleTracking::update()
{
    if ( portNumFrames.isNew() )
    {
        // Number of dynamic connections has changed: create new connections / delete extra ones.
        createFrameConnections( portNumFrames.getValue()-1 );
    }
}

/** Perform the calculation in this module. Called by Avizo. */
void CCPiParticleTracking::compute()
{
    // Check whether the action port button was clicked
    if (!portAction.wasHit()) return;
    
    // Access the input data object. The member portData, which is of type
    // HxConnection, is inherited from HxModule.
    HxUniformScalarField3 *field = (HxUniformScalarField3*) portData.getSource();
    
    // Clear the tracker of any previous particles and tracks
    m_tracker.clear();
    
    // Add particles from master connection - the first frame
    addParticles(field, 0);
    
    // Assume fields are labelled for now so we need to add particles from each frame
    // to the tracker
    for (int i = 0; i < m_frames.size(); i++) {
        field = (HxUniformScalarField3*) m_frames[i]->getSource();
        if (field) addParticles(field, i+1);
    }
    
    // Link the particles
    int linkrange = 2;
    float displacement = portDisplacement.getValue();
    m_tracker.linkParticles(linkrange, displacement);
    
    // Calculate the tracks now particles are linked
    m_tracker.calculateTracks();
    
    theMsg->stream() << "There are " << m_tracker.getNumberOfTracks() << " tracks" << endl;
    
    // Print the tracks
    m_tracker.printTracks(theMsg->stream());
    
    // Plot the tracks
    plotTracks();
}

/** Add particles to the tracker. Particles are the centroids of objects
 * found in the 3D volume.
 * @param field The 3D volume data.
 * @param frameId The number of the frame.
 */
void CCPiParticleTracking::addParticles(HxUniformScalarField3 *field, size_t frameId)
{
    // Get the min and max values of input data 
    float min = 0.0, max = 0.0;
    field->getRange(min, max);
    
    // The maximum number of particles 
    int numParticles = (int)ceil(max)-(int)floor(min)+1;
    
    // Voxel count and centroid coordinates for each particle
    vector<int> voxelCount(numParticles);
    vector<float> x(numParticles);
    vector<float> y(numParticles);
    vector<float> z(numParticles);
    for (int i = 0; i < numParticles; i++) {
        x[i] = y[i] = z[i] = 0.0;
        voxelCount[i] = 0;
    }
    
    // Size of data volume
	McDim3l dims = field->lattice().getDims();
    
    // Loop over field to find voxels for each label and add coordinate
    for (int k = 0; k < dims[2]; k++) {
        for (int j = 0; j < dims[1]; j++) {
            for (int i = 0; i < dims[0]; i++) {
                int value = field->evalReg(i,j,k);
                // Ignore zero values since that's the background
                if (value != 0) {
                    x[value] += i;
                    y[value] += j;
                    z[value] += k;
                    voxelCount[value]++;
                }
            }
        }
    }
    
    // Calculate centroids. Ignore i = 0 as that's the background
    for (int i = 1; i < numParticles; i++) {
        x[i] /= (float)voxelCount[i];
        y[i] /= (float)voxelCount[i];
        z[i] /= (float)voxelCount[i];
        m_tracker.addParticle(x[i], y[i], z[i], frameId);
    }
}

/**
 * Plot the tracks of particles
 */
void CCPiParticleTracking::plotTracks()
{
    // Clear the existing tracks
    m_scene->removeAllChildren();
    
    // Loop over all tracks
    for (size_t i = 0; i < m_tracker.getNumberOfTracks(); i++) {
        vector< vector<float> > trackPoints = m_tracker.getTrack(i);
        
        // Use the technique from the Open Inventor guide at
        // http://www.vsg3d.com/support/oiv_doc/TechDoc/PG-GettingStarted.pdf
        // Pages numbered 41 and 42.
        
        // Loop over points for this track and set up coords data structure.
        // Also draw a sphere for each point
        size_t numTrackPoints = trackPoints.size();
        // Set up vertex property node for scene. Holds coords of vertices for
        // track lines
        SoVertexProperty *vProp = new SoVertexProperty;
        for (size_t j = 0; j < numTrackPoints; j++) {
            // Add sphere - use translation to get position then have to 
            // translate back again.
            SoSphere *sphere = new SoSphere;
            sphere->radius = 2;
            SoTranslation *trans = new SoTranslation;
            trans->translation.setValue(trackPoints[j][0], trackPoints[j][1], trackPoints[j][2]);
            SoTranslation *backTrans = new SoTranslation;
            backTrans->translation.setValue(-trackPoints[j][0], -trackPoints[j][1], -trackPoints[j][2]);
            SoMaterial *material = new SoMaterial;
            material->diffuseColor.setValue(m_sphereColours[i%3]);
            m_scene->addChild(material);
            m_scene->addChild(trans);
            m_scene->addChild(sphere);
            m_scene->addChild(backTrans);
            // Add this track point as a vertex
            vProp->vertex.set1Value((int)j, trackPoints[j][0], trackPoints[j][1], trackPoints[j][2]);
        }
        // Set colour of this track
		vProp->orderedRGBA.setValue( m_colours[i%3] );
        
        // Now draw lines
        int numPoints[] = {(int)trackPoints.size()};
        SoLineSet *lineSet = new SoLineSet;
        lineSet->numVertices.setValues(0, 1, numPoints);
        lineSet->vertexProperty.setValue(vProp);
        
        // Add to scene
        m_scene->addChild(lineSet);

    }
    
    showGeom(m_scene);
}
