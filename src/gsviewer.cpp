//#include <GL/glut.h>
#include <GL/freeglut.h>
#include <iostream>
#include <algorithm>

#include "gsviewer.hpp"


#define TIMERSECS 1000

int width = 1000;
int height = 1000;

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
//    glutDisplayFunc(display);
//    glutIdleFunc(display);
    glutTimerFunc(TIMERSECS, display, 0);
    
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

void GSViewer::display(int value)
{
    glutTimerFunc(10, &display, 0);
    
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // perform simulation step
    if (!pauseFlag) {
        simulation->step();
    }
    
    int N = simulation->size(); // number of cells in each direction
    
    std::vector<double> field = simulation->getU();
    
    double min = *std::min_element(field.begin(), field.end());
    double max = *std::max_element(field.begin(), field.end());
    
    glPushMatrix();
    glScalef(2./(double)width, 2./(double)height, 1);
    glTranslatef(-width/2.,-height/2.,0);
    
    
    // to test coordinate system
//    glColor3f(1,1,1);
//    glRecti(0.25*width,0.25*height,width,height);

//    // bottom left corner is (0,0), top right corner is (width,height)
//    glPointSize(50);
//    glBegin( GL_POINTS );
//        glColor3d(1,0,0);
//        glVertex2d(0,0);
//        
//        glColor3d(0,1,0);
//        glVertex2d(width,height);
//    glEnd();
    
    
    // draw squares
    
    // square width and height respectively
    int sw = width/N;
    int sh = height/N;

    // loop over all cells
    for (int i=0; i<N; ++i) {
        for (int j=0; j<N; ++j) {
            
            double color = F(i,j);
            
            // map color from [min,max] to [-1,1]
            //
            // To map [A, B] --> [a, b]
            // use this formula: (val - A)*(b-a)/(B-A) + a
            
            color = (color-min) * (1.-(-1.)) / (max-min) + (-1.);
            
//            std::cout << color << "\n";
//            std::cout << red(color) << " " <<  green(color) << " " << blue(color) << "\n";
            
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
	{
	}

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

double interpolate( double val, double y0, double x0, double y1, double x1 ) {
    return (val-x0)*(y1-y0)/(x1-x0) + y0;
}

double base( double val ) {
    if ( val <= -0.75 ) return 0;
    else if ( val <= -0.25 ) return interpolate( val, 0.0, -0.75, 1.0, -0.25 );
    else if ( val <= 0.25 ) return 1.0;
    else if ( val <= 0.75 ) return interpolate( val, 1.0, 0.25, 0.0, 0.75 );
    else return 0.0;
}

double red( double gray ) {
    return base( gray + 0.5 );
}
double green( double gray ) {
    return base( gray );
}
double blue( double gray ) {
    return base( gray - 0.5 );
}





