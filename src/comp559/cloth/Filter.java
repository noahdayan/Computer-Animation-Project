package comp559.cloth;

import no.uib.cipr.matrix.Vector;

/**
 * Velocity filter to use with a conjugate gradients solve
 * @author kry
 */
public interface Filter {

    /**
     * removes disallowed parts of v by projection
     * @param v
     */
    public void filter( Vector v );
    
}
