#include <iostream>
#include <boost/program_options.hpp>
#include <mpi.h>

#include "grayscott.hpp"
//#include "gsviewer.hpp"


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
                          std::string& pngName)
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
			("nsteps,s", po::value<int>(&nSteps)->default_value(5000),        "Number of steps"                      )
			("pngname",  po::value<std::string>(&pngName)->default_value("alpha"), "Name for output png"             );

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
//	if (vm.count("visualize")) {
//		visualize = true;
//	}

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
	std::string pngname;
	
	// set/read parameters
	bool result = process_command_line(argc, argv, N, L, dt, Du, Dv, F, k, nSteps, pngname);
	if (!result)
	    return 1;
	
	
	/// initialize MPI domain
    MPI_Init(&argc, &argv);
    
    world_info world;
    
    MPI_Comm_size(MPI_COMM_WORLD, &world.size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world.rank);
    
    
    
    int dims[2] = {0,1};
    MPI_Dims_create(world.size, 2, dims);
    world.dims_x = dims[0];
    world.dims_y = dims[1];
    
    if (world.rank == 0)
        std::cout
        << "processes: " << world.size << "\n"
        << "dims_x: " << world.dims_x << "\n"
        << "dims_y: " << world.dims_y << "\n"
        << std::endl;
    
        
    /// make gridpoints multiple of procs in y-direction
    const int tmp = N;
    N = tmp % world.dims_y == 0 ? tmp : tmp + (world.dims_y - tmp % world.dims_y);
    assert(N % world.dims_y == 0);
    
    
    GrayScott* simulation = new GrayScott(N, -1., 1., dt, Du, Dv, F, k, nSteps, pngname, world);
    
    MPI_Barrier(MPI_COMM_WORLD);

    simulation->run();

    MPI_Barrier(MPI_COMM_WORLD);

    delete simulation;
    
    
    MPI_Finalize();
    
    return 0;
}




