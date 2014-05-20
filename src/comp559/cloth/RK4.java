package comp559.cloth;

public class RK4 implements Integrator {
    
    @Override
    public String getName() {
        return "RK4";
    }

    @Override
    public void step(double[] y, int n, double t, double h, double[] yout, Function derivs) {
        // TODO Auto-generated method stub
        
    	double[] k1 = new double[n];
    	double[] k2 = new double[n];
    	double[] k3 = new double[n];
    	double[] k4 = new double[n];
    	
    	derivs.derivs(t, y, k1);
    	
    	for(int i = 0; i < n; i++)
    	{
    		k2[i] = y[i] + h*k1[i]/2.0;
    	}
    	
    	derivs.derivs(t+h/2.0, k2, k2);
    	
    	for(int i = 0; i < n; i++)
    	{
    		k3[i] = y[i] + h*k2[i]/2.0;
    	}
    	
    	derivs.derivs(t+h/2.0, k3, k3);
    	
    	for(int i = 0; i < n; i++)
    	{
    		k4[i] = y[i] + h*k3[i];
    	}
    	
    	derivs.derivs(t+h, k4, k4);
    	
    	for(int i = 0; i < n; i++)
    	{
    		yout[i] = y[i] + (1.0/6.0)*h*k1[i] + (1.0/3.0)*h*k2[i] + (1.0/3.0)*h*k3[i] + (1.0/6.0)*h*k4[i];
    	}
    }
}
