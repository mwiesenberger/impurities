#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>

#include "draw/host_window.h"
//#include "draw/device_window.cuh"

#include "dg/file/json_utilities.h"

#include "toeflI.cuh"
#include "parameters.h"

int main( int argc, char* argv[])
{
    //Parameter initialisation
    std::stringstream title;
    Json::Value js;
    if( argc == 1)
        dg::file::file2Json( "input.json", js, dg::file::comments::are_discarded);
    else if( argc == 2)
        dg::file::file2Json( argv[1], js, dg::file::comments::are_discarded);
    else
    {
        std::cerr << "ERROR: Too many arguments!\nUsage: "<< argv[0]<<" [filename]\n";
        return -1;
    }
    const imp::Parameters p( js);
    p.display( std::cout);
    /////////glfw initialisation ////////////////////////////////////////////
    dg::file::file2Json( "window_params.json", js, dg::file::comments::are_discarded);
    GLFWwindow* w = draw::glfwInitAndCreateWindow( js["width"].asDouble(), js["height"].asDouble(), "");
    draw::RenderHostData render(js["rows"].asDouble(), js["cols"].asDouble());
    /////////////////////////////////////////////////////////////////////////
    ////////////////////////////////set up computations///////////////////////////
    dg::CartesianGrid2d grid( 0, p.lx, 0, p.ly, p.n, p.Nx, p.Ny, p.bc_x, p.bc_y);
    //create RHS 
    dg::Diffusion< dg::CartesianGrid2d, dg::DMatrix, dg::DVec > diffusion( grid, p); 
    dg::ToeflI< dg::CartesianGrid2d, dg::DMatrix, dg::DVec > toeflI( grid, p); 

    //create initial vector
    dg::Gaussian gaussian( p.posX*grid.lx(), p.posY*grid.ly(), p.sigma, p.sigma, p.amp); //gaussian width is in absolute values
    std::vector<dg::DVec> y0(3, dg::DVec( grid.size()) );
    dg::Helmholtz<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> & gamma = toeflI.gamma();

    //////////////////initialisation of timestepper and first step///////////////////
    dg::DefaultSolver<std::vector<dg::DVec> > solver( diffusion, y0, y0[0].size(), p.eps_time);
    dg::ImExMultistep< std::vector<dg::DVec> > karniadakis( "ImEx-BDF-3-3", y0);
    karniadakis.init( std::tie( toeflI, diffusion, solver), 0, y0, p.dt);


    dg::DVec dvisual( grid.size(), 0.);
    dg::HVec hvisual( grid.size(), 0.), visual(hvisual);
    dg::IHMatrix equi = dg::create::backscatter( grid);
    draw::ColorMapRedBlueExt colors( 1.);
    //create timer
    dg::Timer t;
    double time = 0;
    const double mass_blob0 = toeflI.mass();
    double E0 = toeflI.energy(), energy0 = E0, E1 = 0, diff = 0;
    std::cout << "Begin computation \n";
    std::cout << std::scientific << std::setprecision( 2);
    unsigned step = 0;
    dg::Elliptic<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> laplacianM( grid,  dg::centered);
    while ( !glfwWindowShouldClose( w ))
    {
        //transform field to an equidistant grid
        title << std::setprecision(2) << std::scientific;
        for( unsigned i=0; i<y0.size(); i++)
        {
            dg::assign( y0[i], hvisual);
            dg::blas2::gemv( equi, hvisual, visual);
            //compute the color scale
            colors.scale() =  (float)thrust::reduce( visual.begin(), visual.end(), 0., dg::AbsMax<double>() );
            if( colors.scale() == 0) 
                colors.scale() = 1.;
            //draw ions
            title <<"n / "<<colors.scale()<<"\t";
            render.renderQuad( visual, grid.n()*grid.Nx(), grid.n()*grid.Ny(), colors);
        }
        //transform phi
        dg::blas2::gemv( laplacianM, toeflI.potential()[0], y0[1]);
        dg::assign( y0[1], hvisual);
        dg::blas2::gemv( equi, hvisual, visual);
        //compute the color scale
        colors.scale() =  (float)thrust::reduce( visual.begin(), visual.end(), 0., dg::AbsMax<double>() );
        //draw phi and swap buffers
        title <<"omega / "<<colors.scale()<<"\t";
        title << std::fixed; 
        title << " &&   time = "<<time;
        render.renderQuad( visual, grid.n()*grid.Nx(), grid.n()*grid.Ny(), colors);
        glfwSetWindowTitle(w,title.str().c_str());
        title.str("");
        glfwPollEvents();
        glfwSwapBuffers( w);

        //step 
#ifdef DG_BENCHMARK
        t.tic();
#endif//DG_BENCHMARK
        for( unsigned i=0; i<p.itstp; i++)
        {
            step++;
            std::cout << "(m_tot-m_0)/m_0: "<< (toeflI.mass()-mass_blob0)/mass_blob0<<"\t";
            E0 = E1;
            E1 = toeflI.energy();
            diff = (E1 - E0)/p.dt;
            double diss = toeflI.energy_diffusion( );
            std::cout << "(E_tot-E_0)/E_0: "<< (E1-energy0)/energy0<<"\t";
            std::cout << diff << " "<<diss<<"\t";
            std::cout << "Accuracy: "<< 2.*fabs((diff-diss)/(diff+diss))<<"\n";

            try{ karniadakis.step( std::tie( toeflI, diffusion, solver), time, y0);}
            catch( dg::Fail& fail) {
                std::cerr << "CG failed to converge to "<<fail.epsilon()<<"\n";
                std::cerr << "Does Simulation respect CFL condition?\n";
                glfwSetWindowShouldClose( w, GL_TRUE);
                break;
            }
        }
#ifdef DG_BENCHMARK
        t.toc();
        std::cout << "\n\t Step "<<step;
        std::cout << "\n\t Average time for one step: "<<t.diff()/(double)p.itstp<<"s\n\n";
#endif//DG_BENCHMARK
    }
    glfwTerminate();
    ////////////////////////////////////////////////////////////////////

    return 0;

}
