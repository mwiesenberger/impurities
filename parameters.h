#pragma once

#include <json/json.h>
#include "dg/file/json_utilities.h"

namespace impurities
{
/**
 * @brief Provide a mapping between input file and named parameters
 */
struct Parameters
{
    unsigned n, Nx, Ny;

    unsigned num_stages;
    std::vector<double> eps_pol, eps_gamma;
    enum dg::direction pol_dir, diff_dir;

    double x[2], y[2];
    std::map<std::string, enum dg::bc> bcx, bcy;

    double kappa;
    double epsilon_D;

    std::map<std::string, double> a, mu, tau, nu_perp;

    std::vector<std::string> species;
    unsigned num_species;

    Parameters(){}

    /**
     * @brief constructor to make a const object
     *
     * @param js json object
     */
    Parameters( const dg::file::WrappedJsonValue& js) {
        n  = js["grid"].get("n",3).asUInt();
        Nx = js["grid"].get("Nx", 100).asUInt();
        Ny = js["grid"].get("Ny", 100).asUInt();

        num_stages = js["elliptic"]["stages"].asUInt();
        eps_pol.resize(num_stages);
        eps_gamma.resize(num_stages);
        eps_pol[0] = js["elliptic"]["eps_pol"][0].asDouble();
        eps_gamma[0] = js["elliptic"]["eps_gamma"][0].asDouble();
        for( unsigned u=1;u<num_stages; u++)
        {
            eps_pol[u] = js["elliptic"][ "eps_pol"][u].asDouble();
            eps_gamma[u] = js["elliptic"]["eps_gamma"][u].asDouble();
            eps_pol[u]*=eps_pol[0];
            eps_gamma[u]*=eps_gamma[0];
        }
        pol_dir =  dg::str2direction( js["elliptic"]["direction"].asString());
        diff_dir = dg::centered;
        kappa = js["curvature"].asDouble();
        epsilon_D = js["potential"]["epsilon_D"].asDouble();
        x[0] = js["grid"]["x"].get(0,0.).asDouble();
        x[1] = js["grid"]["x"].get(1,200.).asDouble();
        y[0] = js["grid"]["y"].get(0,0.).asDouble();
        y[1] = js["grid"]["y"].get(1,200.).asDouble();
        num_species = js["species"].size();
        double suma = 0.;
        for( unsigned u=0; u<num_species; u++)
        {
            std::string name = js["species"][u]["name"].asString();
            species.push_back( name);
            a[name] = js["species"][u]["a"].asDouble();
            suma += a[name];
            tau[name] = js["species"][u]["tau"].asDouble();
            mu[name] = js["species"][u]["mu"].asDouble();
            nu_perp[name] = js["species"][u]["nu_perp"].asDouble();
            bcx[name] = dg::str2bc(js["species"][u]["bc"][0].asString());
            bcy[name] = dg::str2bc(js["species"][u]["bc"][1].asString());
        }
        bcx["potential"] = dg::str2bc(js["potential"]["bc"][0].asString());
        bcy["potential"] = dg::str2bc(js["potential"]["bc"][1].asString());
        if( !dg::is_same( suma , 0, 1e-15))
            throw dg::Error( dg::Message(_ping_)<<" The sum of a`s is not zero but "<<suma<<"\n");
    }

    void display( std::ostream& os = std::cout ) const
    {
        os << "Physical parameters are: \n"
            <<"    Curvature_y:     = "<<kappa<<"\n";
        os << "Number of species "<<num_species<<"\n";
        for( unsigned u=0; u<num_species; u++)
        {
            std::string s = species[u];
            os <<"     name_"<<s<<"  = "<<s<<"\n"
               <<"        a_"<<s<<"  = "<<a.at(s)<<"\n"
               <<"       mu_"<<s<<"  = "<<mu.at(s)<<"\n"
               <<"      tau_"<<s<<"  = "<<tau.at(s)<<"\n"
               <<"       nu_"<<s<<"  = "<<nu_perp.at(s)<<"\n"
               <<"      bcx_"<<s<<"  = "<<dg::bc2str(bcx.at(s))<<"\n"
               <<"      bcy_"<<s<<"  = "<<dg::bc2str(bcy.at(s))<<"\n";
        }
        os << "Potential\n"
           <<"      eps_D= "<<epsilon_D<<"\n"
           <<"      bcx  = "<<dg::bc2str(bcx.at("potential"))<<"\n"
           <<"      bcy  = "<<dg::bc2str(bcy.at("potential"))<<"\n";

        os << "Boundary parameters are: \n"
            <<"    x = "<<x[0]<<" x "<<x[1]<<"\n"
            <<"    y = "<<y[0]<<" x "<<y[1]<<"\n";
        os << "Algorithmic parameters are: \n"
            <<"    n  = "<<n<<"\n"
            <<"    Nx = "<<Nx<<"\n"
            <<"    Ny = "<<Ny<<"\n";
        os << "Stopping for CG:         "<<eps_pol[0]<<"\n"
           << "    stages =             "<<num_stages<<"\n"
            <<"Stopping for Gamma CG:   "<<eps_gamma[0]<<std::endl;
        //the endl is for the implicit flush
    }
};
}//namespace imp
