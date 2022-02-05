#pragma once
// diag.h
#pragma once
#include "equations.h"
#include "parameters.h"

namespace impurities
{

struct Variables
{
    Equations<dg::x::CartesianGrid2d, dg::x::DMatrix, dg::x::DVec>& rhs;
    const dg::x::CartesianGrid2d& grid;
    const Parameters& p,
    const std::map<std::string, dg::x::DVec>& ns;
};

struct Record
{
    std::string name; // variable name in the output file
    std::string long_name; // longer description as an attribute
    std::function<void(dg::x::DVec&, Variables&, std::string)> function;
    // function that generates the data points for the variable
};

// time - independent output (only called once)
std::vector<Record> diagnostics2d_static_list = {
    { "xc", "x-coordinate in Cartesian coordinate system",
        []( dg::x::DVec& result, Variables& v, std::string s ) {
            result = dg::evaluate( dg::cooX2d, v.grid);
        }
    },
    { "yc", "y-coordinate in Cartesian coordinate system",
        []( dg::x::DVec& result, Variables& v, std::string s ) {
            result = dg::evaluate( dg::cooY2d, v.grid);
        }
    },
    { "weights", "Gaussian Integration weights",
        []( dg::x::DVec& result, Variables& v, std::string s ) {
            result = dg::create::weights( v.grid);
        }
    }
};

// time - dependent output (called periodically)
std::vector<Record> diagnostics2d_s_list = {
    {"n", "Real density in 2d",
        []( dg::x::DVec& result, Variables& v, std::string s) {
             dg::blas1::copy(v.rhs.real_n(s), result);
        }
    },
    {"psi", "stream function",
        []( dg::x::DVec& result, Variables& v, std::string s) {
             dg::blas1::copy(v.rhs.potential.at(s), result);
        }
    },
    {"S", " a_s tau_s N ln N",
        []( dg::x::DVec& result, Variables& v, std::string s) {
             double astaus = v.p.a.at(s)*v.p.tau.at(s);
             dg::blas1::evaluate( result, dg::equals(), [&]DG_DEVICE(
                         double ns){
                     return astaus*v.ns.at(s)*ln(v.ns.at(s)); },
                         v.ns.at(s));
        }
    },
    {"U", " 0.5 a_s tau_s N u_E^2",
        []( dg::x::DVec& result, Variables& v, std::string s) {
             dg::blas1::pointwiseDot( 0.5*v.p.a[s]*v.p.mu.at(s), v.ns[s],
                     v.eqs.uE2(), 0., result);
        }
    }

//        diff_ = -p.nu*blas2::dot( 1., w2d, lapy[0]);
//        double Gi[3];
//        Gi[0] = - blas2::dot( 1., w2d, lapy[0]) - blas2::dot( lapy[0], w2d, lny[0]); // minus 
//        for( unsigned i=1; i<3; i++)
//            Gi[i] = - p.tau[i]*(blas2::dot( 1., w2d, lapy[i]) + blas2::dot( lapy[i], w2d, lny[i])); // minus 
//        double Gphi[3];
//        for( unsigned i=0; i<3; i++)
//            Gphi[i] = -blas2::dot( phi[i], w2d, lapy[i]);
//        //std::cout << "ge "<<Ge<<" gi "<<Gi<<" gphi "<<Gphi<<" gpsi "<<Gpsi<<"\n";
//        ediff_ = p.nu*( Gi[0] - Gphi[0] + p.a[1]*(Gi[1] + Gphi[1]) + p.a[2]*( Gi[2] + Gphi[2]));
};


} //namespace impurities

