package comp559.cloth;

/**
 * Interface for a numerical integration method
 * 
 * Note that the state is stored in a variable called y, but other good variable names 
 * would have been x as used in class, or q as used in class for a system with changed
 * variables or some other representation of state.  
 * 
 * The reason that y is used is because you can draw a graph with time in the horizontal
 * axis, and let y be your normal vertical axis.  This is also the conventional variable
 * name in the numerical recipes book.
 * 
 * y (state)
 * ^
 * |
 * | ,--._/ trajectory
 * |/
 * |
 * +-------> t (time)
 *
 * See the NR book online, chapter 17, for additional information on integration of ODEs:
 * http://www.nrbook.com/nr3/
 * 
 * @author kry
 */
public interface Integrator {

    /** 
     * @return the name of this numerical integration method
     */
    public String getName();
    
    /** 
     * Advances the system at t by h 
     * @param y The state at time h
     * @param n The dimension of the state (i.e., y.length)
     * @param t The current time (in case the derivs function is time dependent)
     * @param h The step size
     * @param yout  The state of the system at time t+h
     * @param derivs The object which computes the derivative of the system state
     */
    public void step(double[] y, int n, double t, double h, double[] yout, Function derivs); 
}
