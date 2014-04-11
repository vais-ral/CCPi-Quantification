/** 
 * @file Implementation file for class that is a particle track in a set of frames.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date May 2013
 */
 
#include "CCPiTrack.h"

/**
 * Constructor.
 * Constructs a trajectory from the given particle array.
 * Sets its length according to information of the first and last particles
 * 
 * @param particles the array containing all the particles defining this track
 */
CCPiTrack::CCPiTrack(vector<CCPiParticle> particles) {

    m_existing_particles = particles;
    // The length is the last trjectory frame - the first frame (first frame can be 0) 
    m_length = (int)(m_existing_particles.back().getFrame() - 
        m_existing_particles[0].getFrame());
    m_numGaps = 0;
}

void CCPiTrack::populateGaps() {

    for (size_t i = 0; i < m_existing_particles.size()-1; i++) {
        // If two sequential particles are more then 1 frame apart - GAP 
        if (m_existing_particles[i+1].getFrame() - m_existing_particles[i].getFrame() > 1) {
            vector<size_t> gap(2);
            gap[0] = i;
            gap[1] = i+1;
            m_gaps.push_back(gap);
            m_numGaps++;
        }
    }
}

/**
 * Print the track to an output stream.
 * @param out The output stream.
 */
void CCPiTrack::printTrack(ostream& out) {
    
    for (size_t i = 0; i < m_existing_particles.size(); i++) {
        out << m_existing_particles[i] << " -> ";
    }
    out << endl;
}


