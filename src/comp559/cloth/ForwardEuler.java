package comp559.cloth;

/**
 * TODO: finish this class!
 * @author kry
 */
public class ForwardEuler implements Integrator {
    
    @Override
    public String getName() {
        return "Forward Euler";
    }
    
    /** 
     * Advances the system at t by h 
     * @param y The state at time h
     * @param n The dimension of the state (i.e., y.length)
     * @param t The current time (in case the derivs function is time dependent)
     * @param h The step size
     * @param yout  The state of the system at time t+h
     * @param derivs The object which computes the derivative of the system state
     */
    @Override
    public void step(double[] y, int n, double t, double h, double[] yout, Function derivs) {
        // TODO: implement this method
        
    	double[] dydt = new double[n];
    	
    	derivs.derivs(t, y, dydt);
    	
    	for(int i = 0; i < n; i++)
    	{
    		yout[i] = y[i] + h*dydt[i];
    	}
    }

}
