/*
 * Created on May 30, 2006
 */
package comp559.cloth;

import comp559.cloth.Filter;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

/**
 * Implementation of conjugate gradient that allows the user to 
 * easily choose a maximum number of iterations or a maximum 
 * relative error.
 * @author Paul Kry
 */
public class ConjugateGradientMTJ {
    
    DenseVector r;
    DenseVector p;
    DenseVector Ap;
    double rmag;
    
    int iterations;    
    
    private Filter filter = null;
    
    /**
     * Creates a new conjugate gradient solver
     * @param n the size of the system
     */
    public ConjugateGradientMTJ( int n ) {
        r = new DenseVector( n );
        p = new DenseVector( n );
        Ap = new DenseVector( n );
    }
    
    /**
     * Sets the filter to use during solves 
     * @param f
     */
    public void setFilter(Filter f) {
        filter = f;
    }
    
    /** 
     * Performs conjugate gradient for the given number of iterations.
     * Note the comments in the code as to places where one might want
     * to filter values to satisfy constraints.
     * @param A
     * @param b
     * @param x
     * @param numIts
     */
    public void solve( Matrix A, Vector b, Vector x, int numIts ) {        
        A.mult(x, r);
        r.scale(-1);
        r.add(b);        
        filter.filter( r );        
        for ( int i = 0; i < numIts; i++ ) {            
            if ( i== 0 ) {
                p.set(r);
                rmag = r.dot(r);
            } else {
                double rmagnew = r.dot(r);
                double betak = rmagnew / rmag;
                rmag = rmagnew;
                p.scale( betak );
                p.add( r );
            }
            A.mult(p,Ap);            
            filter.filter( Ap );            
            double pAp = p.dot(Ap);
            if ( Math.abs(pAp) < 1e-12 ) return;
            if ( rmag < 1e-12 ) return; // close enough!            
            double alphak = rmag / pAp;                        
            x.add( alphak, p );            
            filter.filter( x );            
            r.add( -alphak, Ap );
            // note: might also want to keep track of the residuals to 
            // see how the solution is progressing.
            //double resid = r.norm( Norm.Two );
            //residuals[i+1] = resid;
            //if ( resid/residuals[0] < relerr ) break;
        }
    }
    
}
