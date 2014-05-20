package comp559.cloth;

public class Midpoint implements Integrator {

    @Override
    public String getName() {
        return "midpoint";
    }

    @Override
    public void step(double[] y, int n, double t, double h, double[] yout, Function derivs) {
        // TODO Auto-generated method stub

    	double[] dydt = new double[n];
    	double[] midpoint = new double[n];
    	
    	derivs.derivs(t, y, dydt);
    	
    	for(int i = 0; i < n; i++)
    	{
    		midpoint[i] = y[i] + h*dydt[i]/2.0;
    	}
    	
    	derivs.derivs(t+h/2.0, midpoint, midpoint);
    	
    	for(int i = 0; i < n; i++)
    	{
    		yout[i] = y[i] + h*midpoint[i];
    	}
    }

}
