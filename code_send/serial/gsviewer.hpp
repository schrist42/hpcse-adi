#ifndef GS_VIEWER_HPP
#define GS_VIEWER_HPP


#include "grayscott.hpp"


/**
 * Viewer for the Gray-Scott Reaction Diffusion simulation.
 * 
 * This class implements a viewer for the Gray-Scott Reaction Diffusion simulation.
 * It visualizes the field U as the simulation is running using OpenGL.
 */
class GSViewer
{
public:
    /**
     * Construct a Gray-Scott viewer.
     * 
     * The constructor initializes glut.
     * 
     * @param argc  argument from the main function
     * @param argv  argument from the main function
     */
    GSViewer(int argc, char* argv[]);
    
    /**
     * Run the visualization.
     * 
     * The visualization function will run the the simulation and visualize the field U.
     * It initializes the window and will run the glutMainLoop.
     * 
     * @param simulation    the simulation that has to be run and visualized
     */
    void visualize(GrayScott *simulation);
    
private:
    /**
     * glutReshapeFunc
     * 
     * This function handles the reshaping of the window.
     * It gets called by glutReshapeFunc.
     * The function has to be static for glut.
     * 
     * @param w     new width of the window
     * @param h     new height of the window
     */
    static void resize(int w, int h);
    
    /**
     * glutDisplayFunc
     * 
     * This function handles all the visualization as well as invoces a step in the simulation.
     * It gets called by glutDisplayFunc in the glutMainLoop.
     * The function has to be static for glut.
     */
    static void display();
    
    /**
     * glutKeyboardFunc
     * 
     * Keyboard callback functions.
     */
    static void keyDown(unsigned char key, int x, int y);
    static void keyUp(unsigned char key, int x, int y);
};


#endif // GS_VIEWER_HPP
