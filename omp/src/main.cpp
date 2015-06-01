#include <iostream>
#include <boost/program_options.hpp>

#include "grayscott.hpp"
#include "gsviewer.hpp"


//GrayScott *simulation;


bool process_command_line(int argc, char** argv,
                          int&    N,
                          double& L,
                          double& dt,
                          double& Du,
                          double& Dv,
                          double& F,
                          double& k,
                          int&    nSteps,
                          bool&   visualize,
                          std::string& pngName,
                          unsigned int& nThreads,
                          bool&   benchmark)
{
	// Define and parse the program options
	namespace po = boost::program_options; 
	po::options_description desc("Options"); 
	po::variables_map vm;
	
	try {
		desc.add_options() 
			("help,h",                                                        "Print help message"                   ) 
			("ncells,N", po::value<int>(&N)->default_value(256),              "Number of cells in one dimension"     )
			("lenght,L", po::value<double>(&L)->default_value(2.),            "Length of the domain in one dimension")
			("dt",       po::value<double>(&dt)->default_value(1),            "Time step"                            )
			("du,u",     po::value<double>(&Du)->default_value(2e-5,"2e-5"),  "Diffusion coefficient for u"          ) 
			("dv,v",     po::value<double>(&Dv)->default_value(1e-5,"1e-5"),  "Diffusion coefficient for v"          ) 	
			(",F",       po::value<double>(&F)->default_value(0.007,"0.007"), "Model parameter 1"                    )
			(",k",       po::value<double>(&k)->default_value(0.046,"0.046"), "Model parameter 2"                    )
			("nsteps,s", po::value<int>(&nSteps)->default_value(1000),        "Number of steps"                      )
			("pngname",  po::value<std::string>(&pngName)->default_value("alpha"), "Name for output png"             )
			("nthreads,t", po::value<unsigned int>(&nThreads)->default_value(2,"2"), "Number of threads"             )
			("benchmark",                                                     "Set to perform benchmark"             );

		po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

		// parse help option
		if (vm.count("help")) {
			std::cout << "\n" << desc << "\n"; 
			return false;
		}

		po::notify(vm);
	}
	catch( std::exception & e ) {
		std::cerr << e.what() << "\n\n";
		std::cout << desc;
		return false;
	}
	catch( ... ) {
		std::cerr << "Unknown error!" << "\n\n";
		std::cout << desc;
        return false;
	}
	
	
	// parse options without value
	if (vm.count("visualize")) {
		visualize = true;
	}
	if (vm.count("benchmark")) {
		benchmark = true;
	}

	return true; // everything worked correctly
}



int main(int argc, char* argv[])
{
    int    N;
    double L;
    double dt;
    double Du;
    double Dv;
    double F;
    double k;
    int    nSteps;
	bool   visualize = false;
	std::string pngname;
	unsigned int nThreads;
	bool   benchmark = false;
	
	bool result = process_command_line(argc, argv, N, L, dt, Du, Dv, F, k, nSteps, visualize, pngname, nThreads, benchmark);
	if (!result)
	    return 1;
	    

    GrayScott* simulation = new GrayScott(N, L, dt, Du, Dv, F, k, nSteps, pngname, nThreads);
    
    if (benchmark) {
        simulation->benchmark();
    }
    else {
        simulation->run();
    }

    delete simulation;
    
    return 0;
}




