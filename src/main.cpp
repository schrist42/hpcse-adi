#include <iostream>
#include <boost/program_options.hpp>

#include "grayscott.hpp"


bool process_command_line(int argc, char** argv,
                          int& N,
                          double& L,
                          double& dt,
                          double& Du,
                          double& Dv,
                          double& F,
                          double& k,
                          int& nSteps)
{
	// Define and parse the program options
	namespace po = boost::program_options; 
	po::options_description desc("Options"); 
	po::variables_map vm;
	
	try {
		desc.add_options() 
			("help,h", "Print help message") 
			("ncells,N", po::value<int>(&N)->default_value(100), "Number of cells in one dimension")
			("lenght,L", po::value<double>(&L)->default_value(1.), "Length of the domain in one dimension")
			("dt", po::value<double>(&dt)->default_value(0.1,"0.1"), "Function number to use")
			("du,u", po::value<double>(&Du)->default_value(2e-5,"2e-5"), "Diffusion coefficient for u") 
			("dv,v", po::value<double>(&Dv)->default_value(1e-5,"1e-5"), "Diffusion coefficient for v") 	
			(",F", po::value<double>(&F)->default_value(0.03,"0.03"), "Model parameter 1")
			(",k", po::value<double>(&k)->default_value(0.062,"0.062"), "Model parameter 2")
			("nSteps,s", po::value<int>(&nSteps)->default_value(1e3,"1e3"), "Number of steps");

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
	
	
//	// parse options without value
//	if (vm.count("msp")) {
//		stats.mean_squared_disp_ = true;
//	}

	return true; // everything worked correctly
}



int main(int argc, char* argv[])
{
    int N;
    double L;
    double dt;
    double Du;
    double Dv;
    double F;
    double k;
    int nSteps;
	
	
	bool result = process_command_line(argc, argv, N, L, dt, Du, Dv, F, k, nSteps);
	if (!result)
	    return 1;
	    

    GrayScott gs(N, L, dt, Du, Dv, F, k, nSteps);
    gs.run();
    
    return 0;
}




