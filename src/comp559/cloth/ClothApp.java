package comp559.cloth;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import javax.media.opengl.GL;
import javax.media.opengl.GL2;
import javax.media.opengl.GLAutoDrawable;
import javax.swing.JButton;
import javax.swing.JPanel;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.swing.VerticalFlowPanel;
import mintools.viewer.EasyViewer;
import mintools.viewer.FancyAxis;
import mintools.viewer.Interactor;
import mintools.viewer.SceneGraphNode;

/**
 * Main program for 599 assignment number 2
 * @author kry
 */
public class ClothApp implements SceneGraphNode, Interactor {

    private EasyViewer ev;
    
    private ParticleSystem system;
    
    /**
     * Entry point for application
     * @param args
     */
    public static void main(String[] args) {
        new ClothApp();        
    }
    
    /**
     * Creates the application / scene instance
     */
    public ClothApp() {
        system = new ParticleSystem();
        ev = new EasyViewer( "Cloth", this, new Dimension(640,480), new Dimension(640,480) );
        ev.addInteractor(this);
    }
     
    @Override
    public void init(GLAutoDrawable drawable) {
        GL2 gl = drawable.getGL().getGL2();
        gl.glEnable( GL.GL_BLEND );
        gl.glBlendFunc( GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA );
        gl.glEnable( GL.GL_LINE_SMOOTH );
        gl.glEnable( GL2.GL_POINT_SMOOTH );
        system.init(drawable);
    }
    
    private FancyAxis fa = new FancyAxis(0.5);
       
    private BoxRoom boxRoom = new BoxRoom();
    
    @Override
    public void display(GLAutoDrawable drawable) {
        GL2 gl = drawable.getGL().getGL2();
        
        if ( run.getValue() ) {   
            system.updateParticles( stepsize.getValue() );                
        }

        if ( drawAxis.getValue() ) { fa.draw(gl); }
        boxRoom.display(drawable);
        system.display( drawable );
               
        EasyViewer.beginOverlay( drawable );
        String text = system.toString() + "\n" + 
                      "h = " + stepsize.getValue() + "\n";
        EasyViewer.printTextLines( drawable, text );
        EasyViewer.endOverlay(drawable);    

        if ( run.getValue() || stepped ) {
            stepped = false;        
            if ( record.getValue() ) {
                // write the frame
                File file = new File( "stills/" + dumpName + format.format(nextFrameNum) + ".png" );                                             
                nextFrameNum++;
                file = new File(file.getAbsolutePath().trim());
                ev.snapshot(drawable, file);
            }
        }
    }

    private BooleanParameter record = new BooleanParameter( "record (press ENTER in canvas to toggle)", false );
    
    /** 
     * boolean to signal that the system was stepped and that a 
     * frame should be recorded if recording is enabled
     */
    private boolean stepped = false;
        
    private String dumpName = "dump";
    
    private int nextFrameNum = 0;
    
    private NumberFormat format = new DecimalFormat("00000");
    
    private BooleanParameter run = new BooleanParameter( "simulate", false );

    private DoubleParameter stepsize = new DoubleParameter( "step size", 0.05, 1e-5, 1 );
        
    private BooleanParameter drawAxis = new BooleanParameter( "draw axis (to help get your bearings in 3D)", false );
    
    @Override
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        JButton create1 = new JButton("create test system 1");
        vfp.add( create1 );
        create1.addActionListener( new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                system.createSystem();
            }
        });
        
        vfp.add( record.getControls() );
        vfp.add( run.getControls() );
        vfp.add( stepsize.getSliderControls(true) );
        vfp.add( system.getControls() );
        vfp.add( drawAxis.getControls());
                
        return vfp.getPanel();
    }
    
    @Override
    public void attach(Component component) {

        component.addKeyListener( new KeyAdapter() {
            @Override
            public void keyPressed(KeyEvent e) {
                if ( e.getKeyCode() == KeyEvent.VK_SPACE ) {
                    run.setValue( ! run.getValue() ); 
                } else if ( e.getKeyCode() == KeyEvent.VK_S ) {
                    system.updateParticles( stepsize.getValue() );                
                    stepped = true;
                } else if ( e.getKeyCode() == KeyEvent.VK_R ) {
                    system.resetParticles();                    
                } else if ( e.getKeyCode() == KeyEvent.VK_C ) {                   
                    system.clearParticles();
                } else if ( e.getKeyCode() == KeyEvent.VK_ESCAPE ) {
                    ev.stop();
                } else if ( e.getKeyCode() == KeyEvent.VK_ENTER ) {
                    record.setValue( ! record.getValue() );
                }
                if ( e.getKeyCode() != KeyEvent.VK_ESCAPE ) ev.redisplay();
            }
        } );
    }
 
}
