#pragma once

#include <json/json.h>
#include "dg/file/json_utilities.h"

namespace imp
{
/**
 * @brief Provide a mapping between input file and named parameters
 */
struct Parameters
{
    unsigned n, Nx, Ny;
    unsigned n_out, Nx_out, Ny_out;

    double eps_pol, eps_gamma;

    double lx, ly;
    enum dg::bc bc_x, bc_y;

    double nu, kappa;

    std::vector<double> a, mu, tau;

    unsigned num_species;

    unsigned vorticity;
    unsigned mode;
    double wall_pos, wall_amp, wall_sigma;
    /**
     * @brief constructor to make a const object
     *
     * @param js json object
     */
    Parameters( const dg::file::WrappedJsonValue& js) {
        num_species = 2;
        n  = js["grid"].get("n",3).asUInt();
        Nx = js["grid"].get("Nx", 100).asUInt();
        Ny = js["grid"].get("Ny", 100).asUInt();
        n_out  = js["output"].get("n",3).asUInt();
        Nx_out = js["output"].get("Nx",100).asUInt();
        Ny_out = js["output"].get("Ny",100).asUInt();

        eps_pol = js["eps_pol"].asDouble();
        eps_gamma = js["eps_gamma"].asDouble();
        eps_time = js["eps_time"].asDouble();
        kappa = js["curvature"].asDouble();
        nu = js["nu_perp"].asDouble();
        lx = js["grid"].get("lx",200).asDouble();
        ly = js["grid"].get("ly",200).asDouble();
        bc_x = dg::str2bc(js["bc_x"].asString());
        bc_y = dg::str2bc(js["bc_y"].asString());
        for ( unsigned u=0; i<num_species; u++)
        {
            a[i] = js["physical"]["a"].get( i, 1.).asDouble();
            tau[i] = js["physical"]["tau"].get(i, 1.).asDouble();
            mu[i]  = js["physical"]["mu"].get(i,1.).asDouble();
        }

        //a[0] = -1, a[1] = 1-a[2];
        //mu[0] = 0, mu[1] = 1;
        //tau[0] = -1;
    }

    void display( std::ostream& os = std::cout ) const
    {
        os << "Physical parameters are: \n"
            <<"    Viscosity:       = "<<nu<<"\n"
            <<"    Curvature_y:     = "<<kappa<<"\n"
            <<"    Ion-temperature: = "<<tau[1]<<"\n";
        os << "Number of species "<<num_sepcies<<"\n";
        for( unsigned u=0; u<num_species; u++)
        {
            os <<"    a_"<<u<<"   = "<<a[i]<<"\n"
               <<"    mu_"<<u<<"  = "<<mu[i]<<"\n"
               <<"    tau_"<<u<<" = "<<tau[i]<<"\n";
        }
        os << "Boundary parameters are: \n"
            <<"    lx = "<<lx<<"\n"
            <<"    ly = "<<ly<<"\n";
        os << "Boundary conditions in x are: \n"
            <<"    "<<bc2str(bc_x)<<"\n";
        os << "Boundary conditions in y are: \n"
            <<"    "<<bc2str(bc_y)<<"\n";
        os << "Algorithmic parameters are: \n"
            <<"    n  = "<<n<<"\n"
            <<"    Nx = "<<Nx<<"\n"
            <<"    Ny = "<<Ny<<"\n"
            <<"    dt = "<<dt<<"\n";
        os << "Stopping for CG:         "<<eps_pol<<"\n"
            <<"Stopping for Gamma CG:   "<<eps_gamma<<std::endl;
        //the endl is for the implicit flush
    }
};
}//namespace imp
