//#include <GL/glut.h>
#include <GL/freeglut.h>
#include <iostream>
#include <algorithm>

#include "gsviewer.hpp"


int width = 256*2;
int height = 256*2;

GrayScott *simulation;

bool pauseFlag = true;


// Helper function declaration
void writeStep();
void drawText(std::string text, double x, double y);
void updateTitle();

double red( double gray );
double green( double gray );
double blue( double gray );



GSViewer::GSViewer(int argc, char* argv[])
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
}


void GSViewer::visualize(GrayScott *sim)
{
    simulation = sim;
    
    glutInitWindowSize(width, height);
    glutInitWindowPosition(10,10);
    
    glutCreateWindow("Gray-Scott Reaction Diffusion - PAUSED, press space"); // paused on start
    
    glutReshapeFunc(resize);
    glutDisplayFunc(display);
    glutIdleFunc(display);
    
    // register GLUT callbacks (e.g. functions to be executed when looping)
	glutKeyboardFunc(keyDown);
	glutKeyboardUpFunc(keyUp);
    
    glClearColor(0,0,0,1);
    
    // Note: glutSetOption is only available with freeglut
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);
    
    glutMainLoop();
}


void GSViewer::resize(int w, int h)
{
    width = w;
    height = h;
    
    const float ar = (float) width / (float) height; 
	
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-ar, ar, -0.25, 1.25, -0.25, 100.0); 
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}



#define F(x,y) field[(x) + (y)*N]

void GSViewer::display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // perform simulation step
    if (!pauseFlag) {
        simulation->step();
    }
    
    int N = simulation->size(); // number of cells in each direction
    
    std::vector<double> field = simulation->getU();
    
    
    glPushMatrix();
    glScalef(2./(double)width, 2./(double)height, 1);
    glTranslatef(-width/2.,-height/2.,0);
    
    
    // draw squares
    
    // square width and height respectively
    int sw = width/N;
    int sh = height/N;

    // loop over all cells
    for (int i=0; i<N; ++i) {
        for (int j=0; j<N; ++j) {
            
            double color = F(i,j);

            glColor3f(red(color), green(color), blue(color));
            
            glRecti(i*sw, j*sh, (i+1)*sw, (j+1)*sh);
        }
    }

    glPopMatrix();
    
    writeStep();
    
    glutSwapBuffers();
}


// Keyboard events
void GSViewer::keyDown(unsigned char key, int x, int y)
{
	glutIgnoreKeyRepeat(1);
	switch(key) 
	{
		// quit
		case 'q':
			//exit(0);
			glutLeaveMainLoop();
			break;
		
		case 27: // esc key
		    //exit(0);
		    glutLeaveMainLoop();
		    break;

		// pause
		case ' ': 
			pauseFlag = !pauseFlag;
			updateTitle();
			break;
	}

	glutPostRedisplay();
}


void GSViewer::keyUp(unsigned char key, int x, int y)
{
	glutIgnoreKeyRepeat(1);
	switch(key)
	{}

	glutPostRedisplay();
}






void writeStep()
{
    double x = (double)width/2.-(double)width/10.;
	double y = (double)height-20.;

    char text[50];
    sprintf(text, "Time = %.2fs", simulation->getTime());
    drawText(text,x,y);
}


void drawText(std::string text, double x, double y)
{
	// to set position in screen coordinates, setup projection and modelview matrices
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0.0, width, 0.0, height);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glColor3f(1,1,1);
	glRasterPos2f(x,y);

	for (int i=0; i<text.length(); i++) { 
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, text.data()[i]);
	}

	// restore back the matrices
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
}


void updateTitle()
{
	std::string title = "Gray-Scott Reaction Diffusion";
	if (pauseFlag)
		title += " - PAUSED, press space";
	glutSetWindowTitle(title.c_str());
}



/////////////////// HELPER FUNCTIONS FOR COLOR /////////////////////////////////
// Source: http://stackoverflow.com/a/7706668

double red( double gray ) {
    if (gray < 0.6)
        return (637.5 * std::max(gray-0.2, 0.)) / 255.;
    else
        return 1;
}
double green( double gray ) {
    if (gray < 0.6)
        return ( 637.5 * std::max(gray-0.2, 0.) ) / 255.;
    else
        return ( -637.5 * (gray-1.) ) / 255.;
}
double blue( double gray ) {
    if (gray < 0.6)
        return ( -637.5 * (gray-0.6) ) / 255.;
    else
        return 0.;
}





