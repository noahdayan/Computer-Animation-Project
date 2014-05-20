package comp559.cloth;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import javax.media.opengl.GL;
import javax.media.opengl.GL2;
import javax.media.opengl.GLAutoDrawable;
import javax.swing.JPanel;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import comp559.cloth.ConjugateGradientMTJ;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;
import mintools.viewer.EasyViewer;
import mintools.viewer.SceneGraphNode;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;

/**
 * Implementation of a simple particle system
 * @author kry
 */
public class ParticleSystem implements SceneGraphNode, Function, Filter {

    /** 
     * The list of all particles. 
     */
    private List<Particle> particles = new LinkedList<Particle>();
    
    /** 
     * The set of springs connecting the particles.  By using a set,
     * it will be easier to remove springs if ever we are deleting particles.
     */
    private Set<Spring> springs = new HashSet<Spring>();
    
    private List<Particle[]> triangles = new LinkedList<Particle[]>();
    
    /** time elapsed in the simulation */
    double time = 0;
    
    /**
     * Creates a particle system
     */
    public ParticleSystem() {
        createSystem();
    }

    /**
     * A grid of particles to approximate a piece of cloth
     */
    private ParticleGrid cloth = new ParticleGrid();
    
    // TODO: note that you may want to define other matrices and vectors.  Just
    // the same, you can probably get by with just one matrix by setting values 
    // in the appropriate places.
    /**
     * The stiffness matrix.  Note that this will actually be a sparse matrix,
     * either a FlexCompRowMatrix or a CompRowMatrix.
     */
    private Matrix K;
    
    public Integrator integrator;
    
    public double[] state = new double[1];
	public double[] stateOut = new double[1];
    
    private ConjugateGradientMTJ CG;
    
    private FlexCompRowMatrix dfdx;
    private FlexCompRowMatrix dfdv;
    private DenseVector deltaxdot;
    private DenseVector b;
    private DenseVector f;
    
    private DenseVector x0;
	private DenseVector v0;
	private DenseVector xNew;
	private DenseVector vNew;
	private FlexCompRowMatrix I;
	private DenseVector tmp;
    
    /**
     * Creates one of a number of simple test systems.
     */
    public void createSystem() {
        
        particles.clear();
        springs.clear();
        triangles.clear();
        
        // use this to test!
        // Particle p1 = new Particle( new Point3d(0, 0, 0), new Vector3d(0, 0, 0), 0 );
        // Particle p2 = new Particle( new Point3d(1, 0, 0), new Vector3d(0, 0, 0), 1 );
        // particles.add( p1 );
        // particles.add( p2 );
        // p1.pinned = true;
        // springs.add( new Spring( p1, p2 ) );
        
        /* // triangle test
        Particle p1 = new Particle( new Point3d(0, -1, 1), new Vector3d(0, 0, 0), 0 );
        Particle p2 = new Particle( new Point3d(-1, -1, -1), new Vector3d(0, 0, 0), 1 );
        Particle p3 = new Particle( new Point3d(1, -1, -1), new Vector3d(0, 0, 0), 2 );
        Particle p4 = new Particle( new Point3d(0, 0, 0), new Vector3d(0, 0, 0), 3 );
        particles.add( p1 );
        particles.add( p2 );
        particles.add( p3 );
        particles.add( p4 );
        p1.pinned = true;
        p2.pinned = true;
        p3.pinned = true;
        springs.add( new Spring( p1, p2 ) );
        springs.add( new Spring( p2, p3 ) );
        springs.add( new Spring( p3, p1 ) );
        Particle[] t = new Particle[3];
        t[0] = p1;
        t[1] = p2;
        t[2] = p3;
        triangles.add(t); */
        
        cloth.create( particles, springs, triangles );

        // TODO: You probably want to initialize your matrices here!
        // For efficiency, you might want to make it as a
        // FlexCompRowMatrix, but then build a CompRowMatrix
        // that has a fixed sparsity structure.  Certainly, you do
        // not want to build the matrix on every call to update particles!
        
        init();
    }
    
    /**
     * Initializes the system 
     * Allocates the arrays and vectors necessary for the solve of the full system
     */
    public void init() {
    	int N = particles.size();
        
        K = new FlexCompRowMatrix(3*N, 3*N);
        
        integrator = new RK4();
        
        CG = new ConjugateGradientMTJ(3*N);
        CG.setFilter(this);
        
        dfdx = new FlexCompRowMatrix(3*N, 3*N);
        dfdv = new FlexCompRowMatrix(3*N, 3*N);
        deltaxdot = new DenseVector(3*N);
        b = new DenseVector(3*N);
        f = new DenseVector(3*N);
    }
    
    /**
     * Gets the particles in the system
     * @return the particle set
     */
    public List<Particle> getParticles() {
        return particles;
    }
    
    /**
     * Resets the positions of all particles
     */
    public void resetParticles() {
        for ( Particle p : particles ) {
            p.reset();
        }
        time = 0;
    }
    
    /**
     * Deletes all particles
     */
    public void clearParticles() {        
        particles.clear();        
        springs.clear();
        cloth.clear();
    }

    /**
     * Updates the position of all particles
     * @param elapsed
     */
    public void updateParticles( double elapsed ) {

        Spring.k = k.getValue();
        Spring.c = c.getValue();
                           
        // TODO: Here again is a crappy forward Euler hack so that the system 
        // does something rather than nothing.  You may want to adapt your 
        // integrators from the previous assignment to try them out here, but
        // ultimately, the goal is to do a backward Euler step here.
        // You need to build the stiffness matrix for the current configuration,
        // combine it with the mass, compute values for the right hand side, and 
        // then solve (probably with conjugate gradient).  There three options
        // for dealing with constraints.  You can include lagrange multipliers
        // in your system, or you can add a filter to the conjugate gradient 
        // code, or you can effectively treat pinned particles as if they were 
        // not part of the system (i.e., treat the system as a lower number of
        // degrees of freedom, though this approach will not let you process the
        // spherical obstacle easily).
        
        if(useImplicit.getValue())
        {
        	int n = getPhaseSpaceDim();
        	
        	if(n != state.length)
            {
            	state = new double[n];
            	stateOut = new double[n];
            }
        	
        	getPhaseSpace(state);
        	backwardEuler(state, n, elapsed, stateOut);
        	setPhaseSpace(stateOut);
        }
        else if(useExplicit.getValue())
        {
        	int n = getPhaseSpaceDim();
        	
        	if(n != state.length)
            {
            	state = new double[n];
            	stateOut = new double[n];
            }
        	
        	getPhaseSpace(state);
        	integrator.step(state, n, time, elapsed, stateOut, this);
        	setPhaseSpace(stateOut);
        }
        else
        {
	        Vector3d tmp = new Vector3d();
	        double h = elapsed / 100;
	        for ( int substeps = 0; substeps < 100; substeps++ ) {
	            for ( Particle p : particles ) {
	                p.f.set( 0, useg.getValue() ? p.mass*g.getValue() : 0, 0 );
	                tmp.scale(-drag.getValue(), p.v);
	                p.f.add(tmp);
	            }
	            for ( Spring s : springs ) {
	                s.apply();
	            }
	            for ( Particle p : particles ) {
	                if ( p.pinned ) continue;
	                
	                tmp.scale( h / p.mass, p.f );
	                p.v.add( tmp );
	                tmp.scale( h, p.v );
	                p.p.add( tmp );
	                p.f.set(0,0,0);
	            }
	        }
        }
        time = time + elapsed;
        
        if(usePenaltyForces.getValue())
        {
        	penaltyForces(elapsed);
        }
        
        if(usePostStepFix.getValue())
        {
        	postStepFixBox();
        	
        	if(useObstacle.getValue())
        	{
        		postStepFixSphere();
        	}
        }
    }
    
    /**
     * Creates a new spring between two particles and adds it to the system.
     * @param p1
     * @param p2
     * @return the new spring
     */
    public Spring createSpring( Particle p1, Particle p2 ) {
        Spring s = new Spring( p1, p2 ); 
        springs.add( s );         
        return s;
    }
    
    @Override
    public void init(GLAutoDrawable drawable) {
        // do nothing
    }

    @Override
    public void display(GLAutoDrawable drawable) {
        GL2 gl = drawable.getGL().getGL2();

        gl.glDisable( GL2.GL_LIGHTING );
        if ( drawPoints.getValue() ) {
            gl.glPointSize( 10 );
            gl.glBegin( GL.GL_POINTS );
            final double[] particleColour = { 0,.95, 0, 0.5 };
            gl.glColor4dv( particleColour, 0 );
            for ( Particle p : particles ) {                
                gl.glVertex3d( p.p.x, p.p.y, p.p.z );
            }
            gl.glEnd();
        }
        gl.glEnable( GL2.GL_LIGHTING );
        
        cloth.display( drawable );
             
        
        if ( useObstacle.getValue() ) {
            float[] frontColour = { .3f, .2f, .2f, 1};
            float[] specColour  = { .8f, .8f, .8f, 1};
            gl.glMaterialfv( GL.GL_FRONT, GL2.GL_DIFFUSE, frontColour, 0 );
            gl.glMaterialfv( GL.GL_FRONT, GL2.GL_SPECULAR, specColour, 0 );
            gl.glMateriali( GL.GL_FRONT, GL2.GL_SHININESS, 92 );
            gl.glPushMatrix();
            gl.glTranslated( xpos.getValue(), ypos.getValue(), zpos.getValue() );
            EasyViewer.glut.glutSolidSphere( radius.getValue()*radiusratio.getValue(), 128, 64 );
            gl.glPopMatrix();
        }
                
        gl.glDisable( GL2.GL_LIGHTING );
       
        if ( drawSprings.getValue() ) {
            gl.glColor4d(0,.5,.5,.5);
            gl.glLineWidth(2f);
            gl.glBegin( GL.GL_LINES );
            for (Spring s : springs) {
                gl.glVertex3d( s.p1.p.x, s.p1.p.y, s.p1.p.z );
                gl.glVertex3d( s.p2.p.x, s.p2.p.y, s.p2.p.z );
            }
            gl.glEnd();
        }
    }
    
    private BooleanParameter useg = new BooleanParameter( "use gravity", true );
    private DoubleParameter g = new DoubleParameter( "gravity", -9.8, -20, 20 );
    private DoubleParameter k = new DoubleParameter( "spring constant", 20, 0, 200 );
    private DoubleParameter c = new DoubleParameter( "spring damping", 0.05, 0, 10 );
    private DoubleParameter drag = new DoubleParameter( "viscous damping", .001, 0, .1 );
    private DoubleParameter restitution = new DoubleParameter( "restitution", 0, 0, 1 );
    private DoubleParameter H = new DoubleParameter( "thickness", 0.1, 0, 1 );
    private IntParameter iterations = new IntParameter( "iterations", 100, 1, 100 );
    
    private BooleanParameter useObstacle = new BooleanParameter( "use obstacle" , true );
    private BooleanParameter useImplicit = new BooleanParameter( "use implicit integration" , true );
    private BooleanParameter useExplicit = new BooleanParameter( "use explicit integration" , false );
    private BooleanParameter usePostStepFix = new BooleanParameter( "use post step fix" , true );
    private BooleanParameter usePenaltyForces = new BooleanParameter( "use penalty forces" , true );
    
    private BooleanParameter drawSprings = new BooleanParameter( "draw springs", false );
    private BooleanParameter drawPoints = new BooleanParameter( "draw points", false );
    
    /** parameters for controlling the position and size of an optional obstacle */
    private DoubleParameter xpos = new DoubleParameter( "x position", .5, -2, 2 );
    private DoubleParameter ypos = new DoubleParameter( "y position", 0, -2, 2 );
    private DoubleParameter zpos = new DoubleParameter( "z position", -.5, -2, 2 );
    private DoubleParameter radius = new DoubleParameter( "radius", .5, 0, 2 );
    private DoubleParameter radiusratio = new DoubleParameter( "radius display ratio", .95, .9, 1 );
    
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        
        vfp.add( cloth.getControls() );
        
        vfp.add( useg.getControls() );
        vfp.add( useObstacle.getControls() );
        vfp.add( useImplicit.getControls() );
        vfp.add( useExplicit.getControls() );
        vfp.add( usePostStepFix.getControls() );
        vfp.add( usePenaltyForces.getControls() );
        
        vfp.add( g.getSliderControls(false) );
        vfp.add( k.getSliderControls(false) );
        vfp.add( c.getSliderControls(false) );
        vfp.add( drag.getSliderControls(false) );
        vfp.add( restitution.getSliderControls(false) );
        vfp.add( H.getSliderControls(false) );
        vfp.add( iterations.getSliderControls() );
        
        vfp.add( xpos.getSliderControls(false) );
        vfp.add( ypos.getSliderControls(false) );
        vfp.add( zpos.getSliderControls(false) );
        vfp.add( radius.getSliderControls(false) );
        vfp.add( radiusratio.getSliderControls(false) );
        
        vfp.add( drawSprings.getControls() );
        vfp.add( drawPoints.getControls() );
                
        return vfp.getPanel();        
    }
    
    @Override
    public String toString() {
        return "particles = " + particles.size() + "\n" +  
               "k = " + k.getValue() + "\n" +
               "b = " + c.getValue() + "\n" +
               "time = " + time;
    }
    
    @Override
    public void filter(Vector v) {
        for ( Particle p : particles ) {
            if ( !p.pinned ) continue;
            v.set( p.index*3+0, 0 );
            v.set( p.index*3+1, 0 );
            v.set( p.index*3+2, 0 );
        }
    }
    
    /**
     *  Evaluates derivatives for ODE integration.
     * @param t time 
     * @param y state
     * @param dydt to be filled with the derivative
     */
    @Override
    public void derivs(double t, double[] y, double[] dydt) {
        // set particle positions to given values
        setPhaseSpace( y );
        
        // TODO: compute forces, and accelerations, and set dydt
        
        for(Particle p: particles)
        {
        	p.clearForce();
        	
        	if(useg.getValue())
            {
                p.f.y = g.getValue() * p.mass;
            }
        	
        	p.f.x -= drag.getValue() * p.v.x;
        	p.f.y -= drag.getValue() * p.v.y;
        	p.f.z -= drag.getValue() * p.v.z;
        }
        
        for(Spring s: springs)
        {
        	s.apply();
        }
        
        int i = 0;
        for(Particle p: particles)
        {
        	dydt[i++] = p.v.x;
        	dydt[i++] = p.v.y;
        	dydt[i++] = p.v.z;
        	dydt[i++] = p.f.x / p.mass;
        	dydt[i++] = p.f.y / p.mass;
        	dydt[i++] = p.f.z / p.mass;
        }
    }
    
    private void backwardEuler(double[] state, int n, double h, double[] stateOut)
    {
    	x0 = new DenseVector(3*particles.size());
    	v0 = new DenseVector(3*particles.size());
    	xNew = new DenseVector(3*particles.size());
    	vNew = new DenseVector(3*particles.size());
    	I = new FlexCompRowMatrix(3*particles.size(), 3*particles.size());
    	tmp = new DenseVector(3*particles.size());
    	
    	for(int i = 0; i < 3*particles.size(); i++)
    	{
    		I.set(i, i, 1.0);
    	}
    	
    	for(int i = 0; i < state.length; i+=6)
    	{
    		int j = i/2;
    		x0.set(j, state[i]);
    		x0.set(j+1, state[i+1]);
    		x0.set(j+2, state[i+2]);
    		v0.set(j, state[i+3]);
    		v0.set(j+1, state[i+4]);
    		v0.set(j+2, state[i+5]);
    	}
    	
    	for(Spring s: springs)
        {
    		s.addForce(f);
    		s.gradient(dfdx);
    		s.gradient2(dfdv);
        }
    	
    	for(Particle p: particles)
    	{
    		f.add(p.index*3, -drag.getValue() * p.v.x);
    		f.add(p.index*3+1, -drag.getValue() * p.v.y);
    		f.add(p.index*3+2, -drag.getValue() * p.v.z);
    		f.add(p.index*3+1, g.getValue() * p.mass);
    	}
    	
    	dfdx.mult(v0, tmp);
    	K.set(I.add(dfdv.scale(-h)).add(dfdx.scale(-h*h)));
    	b.set(f.add(tmp.scale(h)).scale(h));
    	
    	CG.solve(K, b, deltaxdot, iterations.getValue());
    	
    	vNew.set(v0.add(deltaxdot));
    	tmp.set(vNew);
    	xNew.set(x0.add(tmp.scale(h)));
    	
    	for(int i = 0; i < stateOut.length; i+=6)
    	{
    		int j = i/2;
    		stateOut[i] = xNew.get(j);
    		stateOut[i+1] = xNew.get(j+1);
    		stateOut[i+2] = xNew.get(j+2);
    		stateOut[i+3] = vNew.get(j);
    		stateOut[i+4] = vNew.get(j+1);
    		stateOut[i+5] = vNew.get(j+2);
    	}
    }
    
    /**
     * Gets the phase space state of the particle system
     * @param y
     */
    public void getPhaseSpace( double[] y ) {
        int count = 0;
        for ( Particle p : particles ) {
            y[count++] = p.p.x;
            y[count++] = p.p.y;
            y[count++] = p.p.z;
            y[count++] = p.v.x;
            y[count++] = p.v.y;
            y[count++] = p.v.z;
        }
    }
    
    /**
     * Gets the dimension of the phase space state
     * (particles * 3 dimensions * 2 for velocity and position)
     * @return dimension
     */
    public int getPhaseSpaceDim() {        
        return particles.size() * 6;
    }
    
    /**
     * Sets the phase space state of the particle system
     * @param y
     */
    public void setPhaseSpace( double[] y ) {
        int count = 0;
        for ( Particle p : particles ) {
            if ( p.pinned ) {
                count += 6;
            } else {
                p.p.x = y[count++];
                p.p.y = y[count++];
                p.p.z = y[count++];
                p.v.x = y[count++];
                p.v.y = y[count++];
                p.v.z = y[count++];
            }
        }
    }
    
    /**
     * Fixes positions and velocities after a step to deal with wall collisions 
     */
    public void postStepFixBox() {
        for ( Particle p : particles ) {
            if ( p.pinned ) {
                p.v.set(0,0,0);
            }
        }
        // do wall collisions
        double r = restitution.getValue();
        int width = 2;
        int height = 2;
        for ( Particle p : particles ) {            
            if ( p.p.x <= -width ) {
                p.p.x = -width;
                if ( p.v.x < 0 ) p.v.x = - p.v.x * r;
                if ( p.f.x < 0 ) p.f.x = 0;                
            }
            if ( p.p.x >= width ) {
                p.p.x = width;
                if (p.v.x > 0 ) p.v.x = - p.v.x * r;
                if (p.f.x > 0 ) p.f.x = 0;
            } 
            
            if ( p.p.y >= height ) {
                p.p.y = height;
                if ( p.v.y > 0 ) p.v.y = - p.v.y * r;
                if ( p.f.y > 0 ) p.f.y = 0;
            } 
            if ( p.p.y <= -height ) {
                p.p.y = -height;
                if ( p.v.y < 0 ) p.v.y = - p.v.y * r;
                if ( p.f.y < 0 ) p.f.y = 0;
            }
            
            if ( p.p.z <= -width ) {
                p.p.z = -width;
                if ( p.v.z < 0 ) p.v.z = - p.v.z * r;
                if ( p.f.z < 0 ) p.f.z = 0;                
            }
            if ( p.p.z >= width ) {
                p.p.z = width;
                if (p.v.z > 0 ) p.v.z = - p.v.z * r;
                if (p.f.z > 0 ) p.f.z = 0;
            }
        }
    }
    
    /**
     * Fixes positions and velocities after a step to deal with sphere collisions 
     */
    public void postStepFixSphere() {
        for ( Particle p : particles ) {
            if ( p.pinned ) {
                p.v.set(0,0,0);
            }
        }
        // do sphere collisions
        double r = restitution.getValue();
        for ( Particle p : particles ) {            
            if(p.distance(xpos.getValue(), ypos.getValue(), zpos.getValue()) <= radius.getValue())
            {
            	Vector3d v = new Vector3d(p.p.x - xpos.getValue(), p.p.y - ypos.getValue(), p.p.z - zpos.getValue());
            	v.normalize();
            	v.scale(radius.getValue());
            	p.p = new Point3d(xpos.getValue() + v.x, ypos.getValue() + v.y, zpos.getValue() + v.z);
            	
            	if(p.v.length() > 0)
            	{
            		p.v.scale(r);
            		p.v.negate();
            	}
            	if(p.f.length() > 0)
            	{
            		p.f.set(0,0,0);
            	}
            }
        }
    }
    
    private void penaltyForces(double h)
    {
    	for(Particle p: particles)
    	{
    		for(Particle[] t: triangles)
    		{
    			if(p != t[0] && p != t[1] && p != t[2])
    			{
	    			Point3d p1 = t[0].p;
	    			Point3d p2 = t[1].p;
	    			Point3d p3 = t[2].p;
	    			
	    			Vector3d v1 = new Vector3d();
	                Vector3d v2 = new Vector3d();
	    			Vector3d normal = new Vector3d();
	    			
	    			v1.sub(p2,p1);
	                v2.sub(p3,p1);
	                normal.cross(v1,v2);
	                normal.normalize();
	                
	                Vector3d v = new Vector3d();
	                v.sub(p.p,p1);
	                
	                double dot = normal.dot(v);
	                
	                if(Math.abs(dot) < H.getValue())
	                {
	                	double l = v.dot(normal);
	                	normal.scale(l);
	                	
	                	Point3d project = new Point3d();
	                	project.sub(p.p, normal);
	                	
	                	Vector3d v0 = new Vector3d();
	                	v0.sub(project,p1);
	                	
	                	Vector3d normalA = new Vector3d();
	                	Vector3d normalB = new Vector3d();
	                	Vector3d normalC = new Vector3d();
	                	
	                	Vector3d v3 = new Vector3d();
		                Vector3d v4 = new Vector3d();
		                Vector3d v5 = new Vector3d();
		                Vector3d v6 = new Vector3d();
	                	
		                v3.sub(p3,p2);
		                v4.sub(project,p2);
		                v5.sub(p1,p3);
		                v6.sub(project,p3);
		                
		                normalA.cross(v3, v4);
	                	normalB.cross(v5, v6);
	                	normalC.cross(v1, v0);
	                	
	                	normal.cross(v1,v2);
	                	
	                	double a = normal.dot(normalA)/normal.lengthSquared();
	                	double b = normal.dot(normalB)/normal.lengthSquared();
	                	double c = normal.dot(normalC)/normal.lengthSquared();
	                	
	                	normal.normalize();
	                	
	                	if(0 < a && a < 1 && 0 < b && b < 1 && 0 < c && c < 1)
	                	{
	                		Vector3d temp = new Vector3d();
	                		temp.sub(project,p.p);
	                		
	                		if(normal.dot(p.v) < 0)
	                		{
	                			normal.negate();
	                		}
	                		
	                		double d = H.getValue() - normal.dot(temp);
	                		
	                		Vector3d vRel = new Vector3d();
	                		Vector3d vTriangle = new Vector3d(a*t[0].v.x + b*t[1].v.x + c*t[2].v.x, a*t[0].v.y + b*t[1].v.y + c*t[2].v.y, a*t[0].v.z + b*t[1].v.z + c*t[2].v.z);
	                		
	                		vRel.sub(p.v, vTriangle);
	                		
	                		if(vRel.dot(normal) >= 0.1*d/h)
	    					{
	    						double j = -Math.min(h*Spring.k*d, p.mass*(0.1*d/h - vRel.dot(normal)));
	    						double i = 2.0*j / (1.0+a*a+b*b+c*c);
	    						
	    						if(!p.pinned)
	    						{
	    							p.v.add(new Vector3d(-i*normal.x/p.mass, -i*normal.y/p.mass, -i*normal.z/p.mass));
	    						}
	    						if(!t[0].pinned)
	    						{
	    							t[0].v.add(new Vector3d(i*a*normal.x/t[0].mass, i*a*normal.y/t[0].mass, i*a*normal.z/t[0].mass));
	    						}
	    						if(!t[1].pinned)
	    						{
	    							t[1].v.add(new Vector3d(i*b*normal.x/t[1].mass, i*b*normal.y/t[1].mass, i*b*normal.z/t[1].mass));
	    						}
	    						if(!t[2].pinned)
	    						{
	    							t[2].v.add(new Vector3d(i*c*normal.x/t[2].mass, i*c*normal.y/t[2].mass, i*c*normal.z/t[2].mass));
	    						}
	    					}
	                	}
	                }
    			}
    		}
    	}
    }
}

    