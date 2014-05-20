package comp559.cloth;

import java.util.ArrayList;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * Particle class for 599 assignment 2
 * @author kry
 */
public class Particle {

    /** Pinned property for this particle.  Particle should not move if true. */
    boolean pinned = false;
       
    /** Mass of the particle in kg */
    double mass = 1;
    
    /** Index of this particle in the array of particles */
    int index;
        
    /** Position of the particle */
    Point3d p = new Point3d();
    
    /** Velocity of the particle */
    Vector3d v = new Vector3d();
    
    /** Initial position of the particle */
    Point3d p0 = new Point3d();
    
    /** Initial velocity of the particle */
    Vector3d v0 = new Vector3d();
    
    /** 
     * Forces acting on the particle.  Accumulate them here, but
     * remember to reset them to zero before starting.
     */
    Vector3d f = new Vector3d();
    
    /**
     * A list of springs to which this particle is attached.  This
     * is probably not necessary unless you want to remove a particle
     * from your system (i.e., remove the particle, and delete the 
     * springs that were attached to it).
     */
    ArrayList<Spring> springs = new ArrayList<Spring>();
    
    /**
     * Creates a new particle with the given initial position and velocity
     * @param p0 the initial position of the particle
     * @param v0 the initial velocity of the particle
     * @param index an index for looking this particle up in the particles array
     */
    public Particle( Point3d p0, Vector3d v0, int index ) {
        this.p0.set(p0);
        this.v0.set(v0);
        this.index = index;
        reset();
    }
    
    /**
     * Resets the position of this particle
     */
    public void reset() {
        p.set(p0);
        v.set(v0);
        f.set(0,0,0);
    }
    
    /**
     * Adds the given force to this particle
     * @param force
     */
    public void addForce( Vector3d force ) {
        f.add(force);
    }
   
    /**
     * Clears all forces acting on this particle
     */
    public void clearForce() {
        f.set(0,0,0);
    }
    
    /**
     * Computes the distance of a point to this particle
     * @param x
     * @param y
     * @param z
     * @return the distance
     */
    public double distance( double x, double y, double z ) {
        Point3d tmp = new Point3d( x, y, z );
        return tmp.distance(p);
    }
}
