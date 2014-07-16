/** 
 * @file Header file for class that carries out particle tracking of multiple
 * particles in a set of frames.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date May 2013
 */

#ifndef CCPIPARTICLETRACKER_H
#define CCPIPARTICLETRACKER_H

#include "CCPiDefines.h"
#include <functional>
#include <iostream>
#include <list>
#include <vector>

using namespace std;

#include "CCPiParticle.h"
#include "CCPiTrack.h"



/**
 * Class for tracking particles over a number of frames.
 * 
 * Use as follows
 * ...
 */
 
class CCPI_EXPORT CCPiParticleTracker {

  public:
  
    /// Default constructor
    CCPiParticleTracker() : m_linkrange(2) {}
    
    /**
     * Add a particle to a frame
     * @param x CCPiParticle x position
     * @param y CCPiParticle y position
     * @param z CCPiParticle z position
     * @param frameId Frame that particle is to be added to. First frame
     * is frame 0.
     */
    void addParticle(double x, double y, double z, size_t frameId);
     
    /**
     * Compute the tracks
     */
    void linkParticles(int linkrange, float displacement);
    
    /**
     * Generates tracks according to the information available in 
     * each frame and particle. 
     */
    void calculateTracks();
    
    /**
     * Clears all the data so there are no particles and no tracks. 
     */
    void clear();
    
    /** 
     * Get the number of tracks calculated
     * @return The number of tracks
     */
    size_t getNumberOfTracks() {
        return m_tracks.size();
    }
    
    /**
     * Get a track. Data given as a vector of the 3D points for the particles
     * in the track.
     * @return Vector of points in the track. Empty if the track doesn't exist.
     */
    vector< vector<float> > getTrack(size_t trackId);
    
    /** 
     * Print a track
     * @param out Stream to print to
     * @param trackId The track to print. Will print a warning message if 
     * this track number is too big.
     * @return Output stream after printing
     */
    void printTrack(ostream& out, size_t trackId);
    
    /** 
     * Print all the tracks
     * @param out Stream to print to
     * @return Output stream after printing
     */
    void printTracks(ostream& out);
    
  private:
  
    /// The frames
    vector< vector<CCPiParticle> > m_frames;
    
    /// The tracks
    vector< CCPiTrack > m_tracks;
    
    /// The range for calculating links
    size_t m_linkrange;
       
};

#endif // CCPIPARTICLETRACKER_H
