package comp559.cloth;

public class SymplecticEuler implements Integrator {

    @Override
    public String getName() {
        return "symplectic Euler";
    }

    @Override
    public void step(double[] y, int n, double t, double h, double[] yout, Function derivs) {
        // TODO Auto-generated method stub

    	double[] dydt = new double[n];
    	
    	derivs.derivs(t, y, dydt);
    	
    	for(int i = 0; i < n; i+=6)
    	{
    		yout[i+3] = y[i+3] + h*dydt[i+3];
    		yout[i+4] = y[i+4] + h*dydt[i+4];
    		yout[i+5] = y[i+5] + h*dydt[i+5];
    		yout[i] = y[i] + h*yout[i+3];
    		yout[i+1] = y[i+1] + h*yout[i+4];
    		yout[i+2] = y[i+2] + h*yout[i+5];
    	}
    }

}
