/** 
 * @file Header file for class representing a particle in particle tracking module.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date May 2013
 */

#ifndef CCPIPARTICLE_H
#define CCPIPARTICLE_H

#include "CCPiDefines.h"
#include <iostream>
#include <vector>

using namespace std;

/**
 * Defines a particle that holds all the relevant info for it.
 */
class CCPI_EXPORT CCPiParticle {

  public:
  
    /// Original coordinates. Can be refined
	double m_x, m_y, m_z;
    /// Original coordinates. Not to be changed.
	double m_original_x, m_original_y, m_original_z;
    
    /// Number of the frame this particle belonges to - 0 based
	size_t m_frame;
    
    /// Flag that is used while detecting and linking particles
	bool m_special;
    
    /// Array that holds in position i the particle number in frame i
    /// that this particle is linked to 
	vector<int> m_next;
     
	/// Intensity moments. Not used yet but will be when we get to particle
    /// detection from images.
	float m_m0;						
	float m_m1, m_m2, m_m3, m_m4;
    
    /// Non-particle discrimination score
	float m_score;
	
    /// Number of frames this particle is linked over
	size_t m_linkrange;	

    // For debugging?
	int m_nbIterations;

	/**
	 * Constructor. 
	 * @param x Original x coordinates
	 * @param y Original y coordinates
	 * @param z Original z coordinates
     * @param frame_num Number of the frame this particle belonges to
     * @param linkrange
	 */
	CCPiParticle(double x, double y, double z, size_t frame_num, size_t linkrange);

    /**
     * Set the frame number of this particle
     * @param frame The frame number
     */
	void setFrame(size_t frame) {
		m_frame = frame;
	}

    /**
     * Get the frame number of this particle
     * @return The frame number of this particle
     */
    size_t getFrame() {
		return m_frame;
	}
    
    /**
     * Set the linkrange for this particle. This is number of subsequent
     * frames we (will) have linked particles for.
     * @param linkrange The link range
     */
    void setLinkrange(size_t linkrange) {
        m_linkrange = linkrange;
        m_next.resize(linkrange);
    }
	
    /**
     * Get the x,y,z coordinates of the particle
     * @return Coordinates as vector. Entries are x,y,z.
     */
	vector<double> getPosition() {
		vector<double> result(3);
        result[0] = m_x;
        result[1] = m_y;
        result[2] = m_z;
		return result;
	}
    
    /**
     * Equality operator
     * Check all the components are identical
     */
    bool operator==(CCPiParticle const& p);
    
};

/// Utility function to print a particle
inline ostream& operator<<(ostream& out, const CCPiParticle& p)
{
    return out << '(' << p.m_x << ',' << p.m_y << ',' << p.m_z << ')';
}

#endif // CCPIPARTICLE_H
