/** 
 * @file Header file for class that is a particle track in a set of frames.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date May 2013
 */

#ifndef CCPITRACK_H
#define CCPITRACK_H

#include "CCPiDefines.h"
#include <iostream>
#include <vector>

using namespace std;

#include "CCPiParticle.h"

/**
 * Defines a track that is basically an array of sequential particles.
 */
 
class CCPI_EXPORT CCPiTrack {

  public:
  
    /**
     * Constructor.
     * Constructs a trajectory from the given particle array.
     * Sets its length according to information of the first and last particles
     * 
     * @param particles the array containing all the particles defining this track
     */
    CCPiTrack(vector<CCPiParticle> particles);
    
    /**
     * Get the particles
     * @return Vector of particles
     */
    vector<CCPiParticle> getParticles() {
        return m_existing_particles;
    }
    
    /**
     * Set the serial number 
     * @param serialNumber The serial number for this track
     */
    void setSerialNumber(int serialNumber) {
        m_serialNumber = serialNumber;
    }
    
    /**
     * Get the serial number 
     * @return The serial number for this track
     */
    int getSerialNumber() {
        return m_serialNumber;
    }

    /**
     * Populates the gaps vector.
     * Each entry represents a gap using the indices of the particles that
     * have a gap between them. Two sequential particles that are more 
     * then 1 frame apart give a gap.
     */
    void populateGaps();
    
    /**
     * Get the gaps in this track.
     * @return The vector of gaps. Each gap is a pair of indices - the first
     * is to the the particle at the start of the gap and the second is to
     * the particle at the end of the gap.
     */
    vector< vector<size_t> > getGaps() {
        return m_gaps;
    }
     
    /**
     * Get the number of gaps
     * @return The number of gaps.
     */
    int getNumGaps() {
        return m_numGaps;
    }
    
    /**
     * Print the track to an output stream.
     * @param out The output stream.
     */
    void printTrack(ostream& out);
    
  private:
  
    /// All particles in this track in order
    vector<CCPiParticle> m_existing_particles;
    /// Number of frames this track spans
    int m_length;

    /** 
     * Array of 2 indices of particles in m_existing_particles vector marking
     * the start and end points of a gap in this track. two sequential 
     * particles that are more then 1 frame apart give a gap.
     */
    vector< vector<size_t> > m_gaps; 	
    int m_numGaps;

    /// Serial number of this trajectory (for report and display)
    int m_serialNumber;



};

#endif // CCPITRACK_H
