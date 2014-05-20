package comp559.cloth;

import javax.vecmath.Vector3d;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

/**
 * Spring class for 599 assignment 2
 * @author kry
 */
public class Spring {

    Particle p1 = null;
    Particle p2 = null;
    
    static double k = 1;
    static double c = 1;
    
    double l0 = 0;
    
    /**
     * Creates a spring between the two given particles.
     * @param p1
     * @param p2
     */
    public Spring( Particle p1, Particle p2 ) {
        this.p1 = p1;
        this.p2 = p2;
        recomputeRestLength();
        p1.springs.add(this);
        p2.springs.add(this);
    }
    
    /**
     * Recomputes the rest length between two particles
     * from their initial positions.
     */
    public void recomputeRestLength() {
        l0 = p1.p0.distance( p2.p0 );
    }
    
    /**
     * Apply equal and opposite spring force and spring damping force 
     * to the two particles that this spring connects.
     */
    public void apply() {
        Vector3d force = new Vector3d();
        
        force.sub( p2.p, p1.p );
        double l = force.length();
        force.normalize();
        force.scale( (l-l0)*k );
        p1.addForce(force);
        force.scale(-1);
        p2.addForce(force);
        
        force.sub( p2.p, p1.p );
        force.normalize();
        Vector3d v = new Vector3d();
        v.sub(p2.v, p1.v);
        double rv = force.dot(v);
        force.scale(c*rv);
        p1.addForce(force);
        force.scale(-1);
        p2.addForce(force);        
    }
    
    /**
     * Computes the force and adds it to the appropriate components of the force vector.
     * @param f
     */
    public void addForce(Vector f) {
    	Vector3d force = new Vector3d();
        
        force.sub( p2.p, p1.p );
        double l = force.length();
        force.normalize();
        force.scale( (l-l0)*k );
        f.add(p1.index*3,force.x);
    	f.add(p1.index*3+1,force.y);
    	f.add(p1.index*3+2,force.z);
        force.scale(-1);
        f.add(p2.index*3,force.x);
    	f.add(p2.index*3+1,force.y);
    	f.add(p2.index*3+2,force.z);
        
        force.sub( p2.p, p1.p );
        force.normalize();
        Vector3d v = new Vector3d();
        v.sub(p2.v, p1.v);
        double rv = force.dot(v);
        force.scale(c*rv);
        f.add(p1.index*3,force.x);
    	f.add(p1.index*3+1,force.y);
    	f.add(p1.index*3+2,force.z);
        force.scale(-1);
        f.add(p2.index*3,force.x);
    	f.add(p2.index*3+1,force.y);
    	f.add(p2.index*3+2,force.z); 
    }
    
    /** 
     * TODO: Compute this spring's contribution to the gradient matrix
     * and modifies the provided matrix accordingly.
     * 
     * @param K the matrix
     */
    public void gradient( Matrix K ) {
        // There are naturally lots of ways to do this, but one 
        // simple way of computing the stiffness matrix is to have 
        // each spring add its contribution, as we saw in class.  Note
        // that each particle has an index or id number which can
        // be useful for knowing which bits you should be modifying.

        // TODO: complete this function!
        
    	// dfdx = -k * ((1 - l0/|l|)(I - l * l^T) + l * l^T)
    	
    	DenseVector tmpV = new DenseVector(3);
    	DenseMatrix tmpM = new DenseMatrix(3,3);
    	DenseMatrix I = new DenseMatrix(3,3);
    	
    	Vector3d l = new Vector3d(p1.p.x - p2.p.x, p1.p.y - p2.p.y, p1.p.z - p2.p.z);
    	
    	tmpV.set(0,l.x/l.length());
    	tmpV.set(1,l.y/l.length());
    	tmpV.set(2,l.z/l.length());
    	
    	tmpM.rank1(tmpV);
    	I.set(tmpM);
    	I.scale(-1);
    	I.add(0, 0, 1.0);
    	I.add(1, 1, 1.0);
    	I.add(2, 2, 1.0);
    	I.scale(1 - l0/l.length());
    	
    	tmpM.add(I);
    	tmpM.scale(-k);
    	
    	K.add(p1.index*3, p1.index*3, tmpM.get(0, 0));
    	K.add(p1.index*3, p1.index*3+1, tmpM.get(0, 1));    	
    	K.add(p1.index*3, p1.index*3+2, tmpM.get(0, 2));
    	K.add(p1.index*3+1, p1.index*3, tmpM.get(1, 0));
    	K.add(p1.index*3+1, p1.index*3+1, tmpM.get(1, 1));
    	K.add(p1.index*3+1, p1.index*3+2, tmpM.get(1, 2));
    	K.add(p1.index*3+2, p1.index*3, tmpM.get(2, 0));
    	K.add(p1.index*3+2, p1.index*3+1, tmpM.get(2, 1));
    	K.add(p1.index*3+2, p1.index*3+2, tmpM.get(2, 2));
    	
    	K.add(p2.index*3, p2.index*3, tmpM.get(0, 0));
    	K.add(p2.index*3, p2.index*3+1, tmpM.get(0, 1));    	
    	K.add(p2.index*3, p2.index*3+2, tmpM.get(0, 2));
    	K.add(p2.index*3+1, p2.index*3, tmpM.get(1, 0));
    	K.add(p2.index*3+1, p2.index*3+1, tmpM.get(1, 1));
    	K.add(p2.index*3+1, p2.index*3+2, tmpM.get(1, 2));
    	K.add(p2.index*3+2, p2.index*3, tmpM.get(2, 0));
    	K.add(p2.index*3+2, p2.index*3+1, tmpM.get(2, 1));
    	K.add(p2.index*3+2, p2.index*3+2, tmpM.get(2, 2));
    	
    	tmpM.scale(-1);
    	
    	K.add(p1.index*3, p2.index*3, tmpM.get(0, 0));
    	K.add(p1.index*3, p2.index*3+1, tmpM.get(0, 1));
    	K.add(p1.index*3, p2.index*3+2, tmpM.get(0, 2));
    	K.add(p1.index*3+1, p2.index*3, tmpM.get(1, 0));
    	K.add(p1.index*3+1, p2.index*3+1, tmpM.get(1, 1));
    	K.add(p1.index*3+1, p2.index*3+2, tmpM.get(1, 2));
    	K.add(p1.index*3+2, p2.index*3, tmpM.get(2, 0));
    	K.add(p1.index*3+2, p2.index*3+1, tmpM.get(2, 1));
    	K.add(p1.index*3+2, p2.index*3+2, tmpM.get(2, 2));
    	
    	K.add(p2.index*3, p1.index*3, tmpM.get(0, 0));
    	K.add(p2.index*3, p1.index*3+1, tmpM.get(0, 1));
    	K.add(p2.index*3, p1.index*3+2, tmpM.get(0, 2));
    	K.add(p2.index*3+1, p1.index*3, tmpM.get(1, 0));
    	K.add(p2.index*3+1, p1.index*3+1, tmpM.get(1, 1));
    	K.add(p2.index*3+1, p1.index*3+2, tmpM.get(1, 2));
    	K.add(p2.index*3+2, p1.index*3, tmpM.get(2, 0));
    	K.add(p2.index*3+2, p1.index*3+1, tmpM.get(2, 1));
    	K.add(p2.index*3+2, p1.index*3+2, tmpM.get(2, 2));
    }
    
    /** 
     * TODO: Compute this spring's contribution to the gradient matrix
     * and modifies the provided matrix accordingly.
     * 
     * @param K the matrix
     */
    public void gradient2( Matrix K ) {
        // There are naturally lots of ways to do this, but one 
        // simple way of computing the stiffness matrix is to have 
        // each spring add its contribution, as we saw in class.  Note
        // that each particle has an index or id number which can
        // be useful for knowing which bits you should be modifying.

        // TODO: complete this function!
    	
    	// dfdv = -c * l * l^T
    	
    	DenseVector tmpV = new DenseVector(3);
    	DenseMatrix tmpM = new DenseMatrix(3,3);
    	
    	Vector3d l = new Vector3d(p1.p.x - p2.p.x, p1.p.y - p2.p.y, p1.p.z - p2.p.z);
    	
    	tmpV.set(0,l.x/l.length());
    	tmpV.set(1,l.y/l.length());
    	tmpV.set(2,l.z/l.length());
    	
    	tmpM.rank1(tmpV);
    	tmpM.scale(-c);
    	
    	K.add(p1.index*3, p1.index*3, tmpM.get(0, 0));
    	K.add(p1.index*3, p1.index*3+1, tmpM.get(0, 1));    	
    	K.add(p1.index*3, p1.index*3+2, tmpM.get(0, 2));
    	K.add(p1.index*3+1, p1.index*3, tmpM.get(1, 0));
    	K.add(p1.index*3+1, p1.index*3+1, tmpM.get(1, 1));
    	K.add(p1.index*3+1, p1.index*3+2, tmpM.get(1, 2));
    	K.add(p1.index*3+2, p1.index*3, tmpM.get(2, 0));
    	K.add(p1.index*3+2, p1.index*3+1, tmpM.get(2, 1));
    	K.add(p1.index*3+2, p1.index*3+2, tmpM.get(2, 2));
    	
    	K.add(p2.index*3, p2.index*3, tmpM.get(0, 0));
    	K.add(p2.index*3, p2.index*3+1, tmpM.get(0, 1));    	
    	K.add(p2.index*3, p2.index*3+2, tmpM.get(0, 2));
    	K.add(p2.index*3+1, p2.index*3, tmpM.get(1, 0));
    	K.add(p2.index*3+1, p2.index*3+1, tmpM.get(1, 1));
    	K.add(p2.index*3+1, p2.index*3+2, tmpM.get(1, 2));
    	K.add(p2.index*3+2, p2.index*3, tmpM.get(2, 0));
    	K.add(p2.index*3+2, p2.index*3+1, tmpM.get(2, 1));
    	K.add(p2.index*3+2, p2.index*3+2, tmpM.get(2, 2));
    	
    	tmpM.scale(-1);
    	
    	K.add(p1.index*3, p2.index*3, tmpM.get(0, 0));
    	K.add(p1.index*3, p2.index*3+1, tmpM.get(0, 1));
    	K.add(p1.index*3, p2.index*3+2, tmpM.get(0, 2));
    	K.add(p1.index*3+1, p2.index*3, tmpM.get(1, 0));
    	K.add(p1.index*3+1, p2.index*3+1, tmpM.get(1, 1));
    	K.add(p1.index*3+1, p2.index*3+2, tmpM.get(1, 2));
    	K.add(p1.index*3+2, p2.index*3, tmpM.get(2, 0));
    	K.add(p1.index*3+2, p2.index*3+1, tmpM.get(2, 1));
    	K.add(p1.index*3+2, p2.index*3+2, tmpM.get(2, 2));
    	
    	K.add(p2.index*3, p1.index*3, tmpM.get(0, 0));
    	K.add(p2.index*3, p1.index*3+1, tmpM.get(0, 1));
    	K.add(p2.index*3, p1.index*3+2, tmpM.get(0, 2));
    	K.add(p2.index*3+1, p1.index*3, tmpM.get(1, 0));
    	K.add(p2.index*3+1, p1.index*3+1, tmpM.get(1, 1));
    	K.add(p2.index*3+1, p1.index*3+2, tmpM.get(1, 2));
    	K.add(p2.index*3+2, p1.index*3, tmpM.get(2, 0));
    	K.add(p2.index*3+2, p1.index*3+1, tmpM.get(2, 1));
    	K.add(p2.index*3+2, p1.index*3+2, tmpM.get(2, 2));
    }
}
