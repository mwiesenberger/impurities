#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#ifdef WITH_MPI
#include <mpi.h>
#endif // WITH_MPI

#ifdef WITH_GLFW
#include "draw/host_window.h"
#endif // WITH_GLFW

#include "dg/algorithm.h"
#include "dg/file/file.h"


#include "impurities.h"
#include "init.h"
#include "diag.h"

int main( int argc, char* argv[])
{
#ifdef WITH_MPI
    dg::mpi_init( argc, argv);
    MPI_Comm comm;
    dg::mpi_init2d( dg::DIR, dg::PER, comm, std::cin, true);
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
#endif //WITH_MPI
    dg::file::WrappedJsonValue js( dg::file::error::is_throw);
    if( argc != 3 )
    {
        DG_RANK0 std::cerr << "ERROR: Wrong number of arguments!\nUsage: "
                << argv[0]<<" [input.json] [output.nc]\n \n"<<std::endl;
        dg::abort_program();
    }
    const impurities::Parameters p;
    try{
        dg::file::file2Json( argv[1], js.asJson(),
                dg::file::comments::are_discarded, dg::file::error::is_throw);
        p = impurities::Parameters( js);
    } catch( std::exception& e) {
        DG_RANK0 std::cerr << "ERROR in input file "<<argv[1]<<std::endl;
        DG_RANK0 std::cerr << e.what()<<std::endl;
        dg::abort_program();
    }
    DG_RANK0 std::cout << js.asJson() << std::endl;

    //Construct grid
    dg::x::CartesianGrid2d grid( p.x[0], p.x[1], p.y[0], p.y[1], p.n, p.Nx, p.Ny,
            dg::DIR, dg::PER
        #ifdef WITH_MPI
        , comm
        #endif //WITH_MPI
        );
    //Construct equations
    impurities::Equations<dg::x::CartesianGrid2d, dg::x::DMatrix, dg::x::DVec>
        eqs( grid, p);
    // Construct initial condition
    std::map<std::string, dg::x::DVec> y0;
    try{
        y0 = impurities::initial_conditions(grid, eqs, p, js);
    }catch ( std::exception& error){
        DG_RANK0 std::cerr << "Error in input file " << argv[1]<< std::endl;
        DG_RANK0 std::cerr << error.what() << std::endl;
        dg::abort_program();
    }
    // Construct timestepper
    std::string tableau;
    double rtol, atol, time = 0.;
    try{
        rtol = js["timestepper"].get("rtol", 1e-5).asDouble();
        atol = js["timestepper"].get("atol", 1e-5).asDouble();
        tableau = js[ "timestepper"].get( "tableau",
                "Bogacki-Shampine-4-2-3").asString();
    }catch ( std::exception& error){
        DG_RANK0 std::cerr << "Error in input file " << argv[1]<< std::endl;
        DG_RANK0 std::cerr << error.what() << std::endl;
        dg::abort_program();
    }
    dg::Adaptive< dg::ERKStep< dg::x::DVec>> adapt(tableau, omega);
    dg::AdaptiveTimeloop<dg::x::DVec> timeloop( adapt, eqs,
                        dg::pid_control, dg::l2norm, rtol, atol);

    ////////////////////////////////////////////////////////////////////

    return 0;

}
