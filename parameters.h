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
    std::string timestepper, tableau;

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

        timestepper = js["timestepper"].get("type", "multistep").asString();
        tableau     = js["timestepper"].get("tableau", "TVB-3-3").asString();

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

};
}//namespace imp
