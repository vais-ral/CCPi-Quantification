/** 
 * @file Implementation file for class representing Particle in particle 
 * tracking module.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date May 2013
 */
 
#include "CCPiParticle.h"

/**
 * Constructor. 
 * @param x Original x coordinates
 * @param y Original y coordinates
 * @param z Original z coordinates
 * @param frame_num Number of the frame this particle belonges to
 * @param linkrange
 */
CCPiParticle::CCPiParticle (double x, double y, double z, size_t frame_num, size_t linkrange) {
    m_x = x;
    m_original_x = x;
    m_y = y;
    m_original_y = y;
    m_z = z;
    m_original_z = z;
    m_special = true;
    setFrame(frame_num);
    setLinkrange(linkrange);
    m_nbIterations = 0;
    m_m0 = m_m1 = m_m2 = m_m3 = m_m4 = 0.0;
    m_score = 0.0;
}

/**
 * Equality operator
 * Check all the components are identical
 */
bool CCPiParticle::operator==(CCPiParticle const& p) {
    
    bool equal = (m_x == p.m_x) && (m_y == p.m_y) && (m_z == p.m_z);
    
    equal = equal && (m_original_x == p.m_original_x) &&
            (m_original_y == p.m_original_y) && (m_original_z == p.m_original_z);
            
    equal = equal && (m_special == p.m_special) && (m_score == p.m_score);
    
    equal = equal && (m_frame == p.m_frame) && (m_next == p.m_next) &&
            (m_linkrange == p.m_linkrange);
    
    equal = equal && (m_nbIterations == p.m_nbIterations) &&
            (m_m0 == p.m_m0) && (m_m1 == p.m_m1) && (m_m2 == p.m_m2) && 
            (m_m3 == p.m_m3) && (m_m4 == p.m_m4);
    
    return equal;
}
