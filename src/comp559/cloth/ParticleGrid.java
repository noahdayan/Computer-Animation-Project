package comp559.cloth;

import java.util.List;
import java.util.Set;

import javax.media.opengl.GL;
import javax.media.opengl.GL2;
import javax.media.opengl.GLAutoDrawable;
import javax.swing.JPanel;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import mintools.parameters.BooleanParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;
import mintools.viewer.SceneGraphNode;

/**
 * A container for a grid of particles.
 * 
 * Note that this grid of particles will look like 
 * a piece of cloth, but the spring based bending 
 * forces are not a fantastic approximation.
 * @author kry
 */
public class ParticleGrid implements SceneGraphNode {

    /**
     * Grid of particles.  This double array helps us keep track of these
     * particles for the purposes of rendering.
     */
    Particle[][] grid = null;
    
    /** 
     * Per particle normals for our grid, since we want to render them 
     * as a surface.
     */
    private Vector3d[][] gridNormals = null;

    /**
     * Creates a new particle grid container, but does not actually
     * create any particles.
     */
    public ParticleGrid() {
        // do nothing
    }
    
    /**
     * Clears the grid of particles
     */
    public void clear() {
        grid = null;
    }
    
    /**
     * Creates a grid of particles and connects them with springs.
     * Two corners of the grid are pinned.
     * The particles and springs created are added to the list and set 
     * objects which are passed as parameters.
     * @param particles
     * @param springs
     */
    public void create( List<Particle> particles, Set<Spring> springs, List<Particle[]> triangles ) {

        Point3d p = new Point3d();
        Vector3d zero = new Vector3d( 0,0,0 );
                    
        int N = (int) gridDimension.getValue();
        double particleMass = 1.0 / (N*N);
        grid = new Particle[N][N];        
        gridNormals = new Vector3d[N][N];
        
        for ( int i = 0; i < N; i++ ) {
            for ( int j = 0; j < N; j++ ) {
                p.x = -1 + 2.0 * i/(N-1);
                p.y = 1;
                p.z = -1 + 2.0 * j/(N-1);                
                grid[i][j] = new Particle( p, zero, particles.size() );
                particles.add(grid[i][j]);                
                grid[i][j].mass = particleMass;                    
                gridNormals[i][j] = new Vector3d();
            }
        }
        
        // create springs connecting the grid of particles
        for ( int i = 0; i < N; i++ ) {
            for ( int j = 0; j < N; j++ ) {
                if ( i+1 < N ) springs.add( new Spring( grid[i][j], grid[i+1][j]) );
                if ( j+1 < N ) springs.add( new Spring( grid[i][j], grid[i][j+1]) );
                if ( i+1 < N && j+1 < N ) {
                    // add cross stitches... even though they get rendered as triangles.
                    // this way the grid forces are at least symmetric
                    springs.add( new Spring( grid[i+1][j], grid[i][j+1]) );
                    springs.add( new Spring( grid[i][j], grid[i+1][j+1]) );                     
                }
            }                
        }
        
        for ( int i = 0; i < grid.length-1; i++ ) {
            for ( int j = 0; j < grid[0].length-1; j++ ) {
            	
            	Particle[] triangle = new Particle[3];
            	
            	triangle[0] = grid[i][j];
            	triangle[1] = grid[i+1][j];
            	triangle[2] = grid[i+1][j+1];
                
            	triangles.add(triangle);
            	
            	triangle[0] = grid[i][j];
            	triangle[1] = grid[i+1][j+1];
            	triangle[2] = grid[i][j+1];
            	
            	triangles.add(triangle);
            	
            }
        }
        
        // pinned at two corners
        grid[0][0].pinned = true;
        grid[N-1][0].pinned = true;
    }
 
    public void init(GLAutoDrawable drawable) {
        // nothing to do!
    }
    
    /**
     * This display call will draw the grid of particles as a 
     * smooth sheet, either with flat shaded triangles, or smooth
     * shaded triangles where per vertex (i.e., per particle) normals
     * are computed by looking at the neighbouring vertices.
     */
    public void display( GLAutoDrawable drawable ) {
        GL2 gl = drawable.getGL().getGL2();
        
        if ( grid != null && drawGrid.getValue() ) {
            Vector3d v1 = new Vector3d();
            Vector3d v2 = new Vector3d();
            
            // draw the cloth using the double sided lighting model and 
            // different materials for front and back.
            gl.glDisable( GL.GL_CULL_FACE );
            gl.glLightModeli( GL2.GL_LIGHT_MODEL_TWO_SIDE, GL.GL_TRUE );            
            float[] frontColour = { .71f, .71f, 0,    1};
            float[] backColour  = { .71f,  0,   .71f, 1};
            gl.glMaterialfv( GL.GL_FRONT, GL2.GL_DIFFUSE, frontColour, 0 );
            gl.glMaterialfv( GL.GL_BACK, GL2.GL_DIFFUSE, backColour, 0 );

            if ( !drawFlatShaded.getValue() ) {    
                // note that this could be made much faster by setting
                // up index and vertex buffers and using glDrawElements
                // but this is just a rendering efficiency issue.
                Vector3d n;
                // compute per vertex normals!                
                for ( int i = 0; i < grid.length; i++ ) {
                    for ( int j = 0; j < grid[0].length; j++ ) {
                        if ( i > 0 && i < grid.length-1) {
                            v1.sub( grid[i+1][j].p, grid[i-1][j].p );                        
                            v1.scale(0.5);
                        } else if ( i > 0 ) {
                            v1.sub( grid[i][j].p, grid[i-1][j].p );                        
                        } else {
                            v1.sub( grid[i+1][j].p, grid[i][j].p );
                        }
                        
                        if ( j > 0 && j < grid[0].length-1) {
                            v2.sub( grid[i][j+1].p, grid[i][j-1].p );                        
                            v2.scale(0.5);
                        } else if ( j > 0 ) {
                            v2.sub( grid[i][j].p, grid[i][j-1].p );                        
                        } else {
                            v2.sub( grid[i][j+1].p, grid[i][j].p );
                        }
                        
                        n = gridNormals[i][j];
                        n.cross(v1,v2);
                        n.normalize();                        
                    }
                }
                gl.glBegin(GL.GL_TRIANGLES);
                for ( int i = 0; i < grid.length-1; i++ ) {
                    for ( int j = 0; j < grid[0].length-1; j++ ) {
                        {
                            Point3d p1 = grid[i][j].p;
                            Point3d p2 = grid[i+1][j].p;
                            Point3d p3 = grid[i+1][j+1].p;
                            n = gridNormals[i][j]; gl.glNormal3d( n.x,n.y,n.z);
                            gl.glVertex3d(p1.x,p1.y,p1.z);
                            n = gridNormals[i+1][j]; gl.glNormal3d( n.x,n.y,n.z);
                            gl.glVertex3d(p2.x,p2.y,p2.z);
                            n = gridNormals[i+1][j+1]; gl.glNormal3d( n.x,n.y,n.z);
                            gl.glVertex3d(p3.x,p3.y,p3.z);
                        }
                        {                    
                            Point3d p1 = grid[i][j].p;
                            Point3d p2 = grid[i+1][j+1].p;
                            Point3d p3 = grid[i][j+1].p;                            
                            n = gridNormals[i][j]; gl.glNormal3d( n.x,n.y,n.z);
                            gl.glVertex3d(p1.x,p1.y,p1.z);
                            n = gridNormals[i+1][j+1]; gl.glNormal3d( n.x,n.y,n.z);
                            gl.glVertex3d(p2.x,p2.y,p2.z);
                            n = gridNormals[i][j+1]; gl.glNormal3d( n.x,n.y,n.z);
                            gl.glVertex3d(p3.x,p3.y,p3.z);
                        }                    
                    }
                }
                gl.glEnd();
            } else {                
                Vector3d n = new Vector3d();                
                gl.glBegin(GL.GL_TRIANGLES);
                for ( int i = 0; i < grid.length-1; i++ ) {
                    for ( int j = 0; j < grid[0].length-1; j++ ) {
                        {
                            Point3d p1 = grid[i][j].p;
                            Point3d p2 = grid[i+1][j].p;
                            Point3d p3 = grid[i+1][j+1].p;
                            v1.sub(p2,p1);
                            v2.sub(p3,p1);
                            n.cross(v1,v2);
                            n.normalize();
                            gl.glNormal3d(n.x,n.y,n.z);
                            gl.glVertex3d(p1.x,p1.y,p1.z);
                            gl.glVertex3d(p2.x,p2.y,p2.z);
                            gl.glVertex3d(p3.x,p3.y,p3.z);
                        }                        
                        {                    
                            Point3d p1 = grid[i][j].p;
                            Point3d p2 = grid[i+1][j+1].p;
                            Point3d p3 = grid[i][j+1].p;
                            v1.sub(p2,p1);
                            v2.sub(p3,p1);
                            n.cross(v1,v2);
                            n.normalize();
                            gl.glNormal3d(n.x,n.y,n.z);
                            gl.glVertex3d(p1.x,p1.y,p1.z);
                            gl.glVertex3d(p2.x,p2.y,p2.z);
                            gl.glVertex3d(p3.x,p3.y,p3.z);
                        }                    
                    }
                }
                gl.glEnd();
            }            
            gl.glLightModeli( GL2.GL_LIGHT_MODEL_TWO_SIDE, GL.GL_FALSE );
            gl.glEnable( GL.GL_CULL_FACE );
        }
    }
    
    private BooleanParameter drawGrid = new BooleanParameter( "draw grid", true );
    private IntParameter gridDimension = new IntParameter ("grid dimension", 15, 2, 40 );
    private BooleanParameter drawFlatShaded = new BooleanParameter( "flat shaded", false );
    
    @Override
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        vfp.add( gridDimension.getControls() );
        vfp.add( drawGrid.getControls() );
        vfp.add( drawFlatShaded.getControls() );
        return vfp.getPanel();
    }
    
}
