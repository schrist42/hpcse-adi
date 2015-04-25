#ifndef GS_VIEWER_HPP
#define GS_VIEWER_HPP


#include "grayscott.hpp"


class GSViewer
{
public:
    GSViewer(int argc, char* argv[]);
    
    void visualize(GrayScott *simulation);
    
private:
    static void resize(int w, int h);
    static void display();
};


#endif
