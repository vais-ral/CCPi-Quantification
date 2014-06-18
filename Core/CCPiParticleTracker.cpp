/** 
 * @file Implementation file for class that carries out particle tracking of multiple
 * particles in a set of frames.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date May 2013
 */
 
#include "CCPiParticleTracker.h"

/**
 * Add a particle to a frame
 * @param x CCPiParticle x position
 * @param y CCPiParticle y position
 * @param z CCPiParticle z position
 * @param frameId Frame that particle is to be added to. First frame
 * is frame 0.
 */
void CCPiParticleTracker::addParticle(double x, double y, double z, 
                                      size_t frameId) 
{

    // Check the number of frames we already have to see if frameId means
    // we need to resize the frames vector
    if (frameId >= m_frames.size()) {
        m_frames.resize(frameId+1);
    }
    
    // Create the particle object and add it to the frame
    CCPiParticle p(x, y, z, frameId, 2);
    m_frames[frameId].push_back(p);
 }
 
void CCPiParticleTracker::linkParticles(int linkrange, float displacement) {
    
    m_linkrange = linkrange;
    size_t numFrames = m_frames.size();

    // set the length of the particles next array according to the linkrange
    // it is done now since link range can be modified after first run
    for (size_t fr = 0; fr < numFrames; fr++) {
        for (size_t pr = 0; pr < m_frames[fr].size(); pr++) {
            m_frames[fr].at(pr).setLinkrange(linkrange);
        }
    }
    int curr_linkrange = linkrange;

    // If the linkrange is too big, set it the right value 
    if ((int)numFrames < (curr_linkrange + 1))
        curr_linkrange = (int)numFrames - 1;

    // Loop over frames
    for (size_t m = 0; m < numFrames - curr_linkrange; m++) {
        //cout << "Linking Frame " << m+1 << endl;
        
        // Set up data for particles in this frame
        size_t numParticles = m_frames[m].size();
        for (size_t i = 0; i < numParticles; i++) {
            m_frames[m].at(i).m_special = false;
            for (int n = 0; n < linkrange; n++)
                m_frames[m].at(i).m_next[n] = -1;
        }

        // Loop forward from current frame
        for (int n = 0; n < curr_linkrange; n++) {
            float max_cost = (float)(n + 1) * displacement * (float)(n + 1) * displacement;

            size_t numParticlesNext = m_frames[m + (n + 1)].size();

            // Set up the cost matrix
            vector< vector<float> > cost(numParticles+1);
			for (size_t i = 0; i < numParticles+1; i++) cost[i].resize(numParticlesNext+1);
            
            // Set up the association matrix
            vector< vector<bool> > g(numParticles+1);
			for (size_t i = 0; i < numParticles+1; i++)g[i].resize(numParticlesNext+1);

            // g_x stores the index of the currently associated particle for 
            // particles in the next frame
            vector<int> g_x(numParticlesNext+1);
            // g_y stores the index of the currently associated particle for 
            // particles in current frame
            vector<int> g_y(numParticles+1);
            
            // okv is a helper vector in the initialization phase. It keeps track of the empty columns
            // in g.
            vector<bool> okv(numParticlesNext+1);
            for (size_t i = 0; i < numParticlesNext+1; i++) okv[i] = true;
            
            vector<CCPiParticle> p1 = m_frames[m];
            vector<CCPiParticle> p2 = m_frames[m + (n + 1)];

            // Fill in the costs
            for (size_t i = 0; i < numParticles; i++) {
                for (size_t j = 0; j < numParticlesNext; j++) {
                    cost[i][j] = 
                        (p1[i].m_x - p2[j].m_x)*(p1[i].m_x - p2[j].m_x) + 
                        (p1[i].m_y - p2[j].m_y)*(p1[i].m_y - p2[j].m_y) + 
                        (p1[i].m_z - p2[j].m_z)*(p1[i].m_z - p2[j].m_z) + 
                        (p1[i].m_m0 - p2[j].m_m0)*(p1[i].m_m0 - p2[j].m_m0) + 
                        (p1[i].m_m2 - p2[j].m_m2)*(p1[i].m_m2 - p2[j].m_m2);
                    g[i][j] = false;
                }
            }

            for (size_t i = 0; i < numParticles + 1; i++) {
                cost[i][numParticlesNext] = max_cost;
                g[i][numParticlesNext] = false;
                g_y[i] = false;
            }
            for (size_t j = 0; j < numParticlesNext + 1; j++) {
                cost[numParticles][j] = max_cost;
                g[numParticles][j] = false;
                g_x[j] = false;
            }
            cost[numParticles][numParticlesNext] = 0.0f;
            g[numParticles][numParticlesNext] = false;

            // Initialize the relation matrix

            for (size_t i = 0; i < numParticles; i++) { // Loop over the x-axis
                //cout << "Linking Frame "  << m+1 << ": Initializing Relation matrix" << endl;

                float min = max_cost;
                int prev = -1;
                for (size_t j = 0; j < numParticlesNext; j++) { // Loop over the y-axis without the dummy

                    // Let's see if we can use this coordinate
                    if (okv[j] && min > cost[i][j]) {
                        min = cost[i][j];
                        if (prev >= 0) {
                            okv[prev] = true;
                            g[i][prev] = false;
                        }

                        okv[j] = false;
                        g[i][j] = true;
                        
                        prev = (int)j;
                    }
                }

                // Check if we have a dummy particle
                if (min == max_cost) {
                    if (prev >= 0) {
                        okv[prev] = true;
                        g[i][prev] = false;
                    }
                    g[i][numParticlesNext] = true;
                    okv[numParticlesNext] = false;
                }
            }
            
            // Look for columns that are zero
            for (size_t j = 0; j < numParticlesNext; j++) {
                int ok = 1;
                for (size_t i = 0; i < numParticles + 1; i++) {
                    if(g[i][j])
                        ok = 0;
                }
                if (ok == 1)
                    g[numParticles][j] = true;
            }

            // Build g_x and g_y, a speedup for g */
            for (size_t i = 0; i < numParticles + 1; i++) {
                for (size_t j = 0; j < numParticlesNext + 1; j++) {
                    if (g[i][j]) {
                        g_x[j] = (int)i;
                        g_y[i] = (int)j;
                    }
                }
            }
            g_x[numParticlesNext] = (int)numParticles;
            g_y[numParticles] = (int)numParticlesNext;
                    
            
            // Now the relation matrix needs to be optimized
            //cout << "Linking Frame " << (m+1) << ": Optimizing Relation matrix" << endl;
            float min = -1.0;
            while (min < 0.0) {
                min = 0.0;
                int prev_i = 0, prev_j = 0, prev_x = 0, prev_y = 0;
                for (size_t i = 0; i < numParticles + 1; i++) {
                    for (size_t j = 0; j < numParticlesNext + 1; j++) {
                        if (i == numParticles && j == numParticlesNext)
                            continue;

                        if ((g[i][j] == false) && (cost[i][j] <= max_cost)) {
                            // Calculate the reduced cost
                            
                            // Look along the x-axis, including
                            // the dummy particles
                            int x = g_x[j];

                            // Look along the y-axis, including
                            // the dummy particles
                            int y = g_y[i];
                            
                            
                            // z is the reduced cost
                            float z = cost[i][j] + 
                                        cost[x][y] - 
                                        cost[i][y] - 
                                        cost[x][j];
                            
                            if (z > -1.0e-10)
                                z = 0.0;
                            if (z < min) {
                                min = z;
                                prev_i = (int)i;
                                prev_j = (int)j;
                                prev_x = (int)x;
                                prev_y = (int)y;
                            }
                        }
                    }
                }

                if (min < 0.0) {
                    g[prev_i][prev_j] = true;
                    g_x[prev_j] = prev_i;
                    g_y[prev_i] = prev_j;
                    g[prev_x][prev_y] = true;
                    g_x[prev_y] = prev_x;
                    g_y[prev_x] = prev_y;
                    g[prev_i][prev_y] = false;
                    g[prev_x][prev_j] = false;
                    
                    // ensure the dummies still map to each other
                    g_x[numParticlesNext] = (int)numParticles;
                    g_y[numParticles] = (int)numParticlesNext;
                }
            }

            // After optimization, the particles needs to be linked
            for (size_t i = 0; i < numParticles; i++) {
                for (size_t j = 0; j < numParticlesNext; j++) {
                    if (g[i][j] == true) {
                        m_frames[m].at(i).m_next[n] = (int)j;
                    }
                }
            }
        }

        if (m == (numFrames - curr_linkrange - 1) && curr_linkrange > 1)
            curr_linkrange--;
    }

    // At the last frame all trajectories end
    for (size_t i = 0; i < m_frames[numFrames-1].size(); i++) {
        m_frames[numFrames-1].at(i).m_special = false;
        for (int n = 0; n < linkrange; n++)
            m_frames[numFrames-1].at(i).m_next[n] = -1;
    }
    
}

/**
 * Generates tracks according to the information available in 
 * each frame and particle. 
 */
void CCPiParticleTracker::calculateTracks() {

    size_t numFrames = m_frames.size();

    // Vector to hold particles for current trajctory
    vector<CCPiParticle> curr_traj_particles;
    curr_traj_particles.reserve(numFrames);		
	
    for (size_t i = 0; i < numFrames; i++) {
        for (size_t j = 0; j < m_frames[i].size(); j++) {
            
            if (!m_frames[i].at(j).m_special) {
                m_frames[i].at(j).m_special = true;
                int found = -1;
                // Go over all particles that this particle (particles[j]) is linked to
                size_t n; // Used when we break from loop so declare here
                for (n = 0; n < m_linkrange; n++) {
                    // If it is NOT a dummy particle - stop looking
                    if (m_frames[i].at(j).m_next[n] != -1) {
                        found = (int)n;
                        break;
                    }
                }
                // If this particle is not linked to any other
                // go to next particle and dont add a track
                if (found == -1)
                    continue;

                // If this particle is linked to a "real" particle that was already linked
                // break the trajectory and start again from the next particle. dont add a trajectory
                if (m_frames[i + n + 1].at(m_frames[i].at(j).m_next[n]).m_special) 
                    continue;

                // This particle is linked to another "real" particle that is not already linked
                // so we have a trajectory					
                curr_traj_particles.push_back(m_frames[i].at(j));
                size_t k = i;
                int m = (int)j;
                do {
                    int found = -1;
                    for (size_t n = 0; n < m_linkrange; n++) {
                        if (m_frames[k].at(m).m_next[n] != -1) {
                            // If this particle is linked to a "real" particle that
                            // that is NOT already linked, continue with building the trajectory
                            if (m_frames[k + n + 1].at(m_frames[k].at(m).m_next[n]).m_special == false) {
                                found = (int)n;
                                break;
                            } else {									
                                break;
                            }
                        }
                    }
                    if (found == -1)
                        break;
                    m = m_frames[k].at(m).m_next[found];
                    k += (found + 1);
                    curr_traj_particles.push_back(m_frames[k].at(m));
                    m_frames[k].at(m).m_special = true;
                } while (m != -1);					

                // Create the current track
                CCPiTrack curr_traj(curr_traj_particles);

                // Set current track parameters
                curr_traj.setSerialNumber((int)m_tracks.size());
                curr_traj.populateGaps();
                // Add current track to the vector
                m_tracks.push_back(curr_traj);
                // Clear temporary vector
                curr_traj_particles.clear();
                curr_traj_particles.reserve(numFrames);		
            }				
        }
    }		
}

/**
 * Clears all the data so there are no particles and no tracks. 
 */
void CCPiParticleTracker::clear()
{
    m_frames.clear();
    m_tracks.clear();
    m_linkrange = 2;
}

/**
 * Get a track. Data given as a vector of the 3D points for the particles
 * in the track.
 * @return Vector of points in the track. Empty if the track doesn't exist.
 */
vector< vector<float> > CCPiParticleTracker::getTrack(size_t trackId)
{
    // The vector to return
    vector< vector<float> > trackPoints;

    // Check trackId actually exists
    if (trackId >= m_tracks.size()) {
        return trackPoints;
    }
    
    // Gather the points from the track
    vector<CCPiParticle> particles = m_tracks[trackId].getParticles();
    
    // Iterate over particles and get the coordinates for return object
    for (size_t i = 0; i < particles.size(); i++) {
        vector<float> coords(3);
        coords[0] = particles[i].m_x;
        coords[1] = particles[i].m_y;
        coords[2] = particles[i].m_z;
        trackPoints.push_back(coords);
    }
 
    return trackPoints;
    
}


/** 
 * Print a track
 * @param out Stream to print to
 * @param trackId The track to print. Will print a warning message if 
 * this track number is too big.
 * @return Output stream after printing
 */
void CCPiParticleTracker::printTrack(ostream& out, size_t trackId) 
{
    // Check trackId actually exists
    if (trackId >= m_tracks.size()) {
        out << "**WARNING** This track does not exist. There are only " <<
            m_tracks.size() << "tracks" << endl;
        return;
    }
            
    // Print the track
    m_tracks[trackId].printTrack(out);
}

/** 
 * Print all the tracks
 * @param out Stream to print to
 * @return Output stream after printing
 */
void CCPiParticleTracker::printTracks(ostream& out) 
{
    // Print the tracks
    for (size_t i = 0; i < m_tracks.size(); i++) {
        m_tracks[i].printTrack(out);
    }
}

