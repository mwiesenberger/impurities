#pragma once

#include "dg/algorithm.h"
#include "parameters.h"

namespace impurities
{

template< class Geometry, class Matrix, class Container >
struct Equations
{
    using value_type = dg::get_value_type<Container>;

    Equations( ) = default;
    Equations( const Geometry& g, Parameters p);

    const Container& psi(std::string name ) const { return m_psi.at(name);}
    const Container& potential( ) const { return m_phi;}

    dg::Helmholtz<Geometry, Matrix, Container>& gamma_inv(std::string name) {
        return m_multi_helmN.at(name)[0];
    }
    const std::function<void( const Container&, Container&)>& gamma(std::string
            name) {return m_inv_helmN.at(name);
    }
    const Container& uE2() const {return m_UE2;}
    void compute_real_n( std::string s, const Container& n, Container& real_n)
    {
        dg::blas1::pointwiseDot( m_binv, m_binv, m_temp);
        dg::blas1::pointwiseDot( m_p.a.at(s)*m_p.mu.at(s), n, m_temp, 0., m_chi);
        m_multi_pol[0].set_chi( m_chi);
        dg::apply( m_multi_pol[0], m_phi, m_omega);
        // m_multi_pol is negative !
        dg::blas1::axpbypgz( -1., m_omega, 1., m_gamma_n.at(s), 0., real_n);
    }
    void compute_lapPhi( Container& lapP)
    {
        m_multi_pol[0].set_chi( 1.);
        dg::blas2::symv( -1., m_multi_pol[0], m_phi, 0., lapP);
    }
    void compute_lapPsi( std::string s, Container& lapP)
    {
        m_multi_pol[0].set_chi( 1.);
        dg::blas2::symv( -1., m_multi_pol[0], m_psi.at(s), 0., lapP);
    }
    void compute_lapN( std::string s, const Container& ns,  Container& lapN)
    {
        dg::blas2::symv( -1., m_laplacianM.at(s), ns, 0., lapN);
    }


    void operator()(double t, const std::map<std::string, Container>& y, std::map<std::string, Container>& yp);
    unsigned ncalls()const {return m_ncalls;}

private:
    Equations( const Equations& ) = delete; // because of lambdas holding this
    Equations( Equations&& ) = delete;      // because of lambdas holding this
    void compute_phi( double t, const std::map<std::string, Container>& y);
    void compute_psi( double t, const Container& potential);

    Container m_binv; //magnetic field

    Container m_chi, m_omega;
    Container m_tilden, m_temp;
    Container m_phi;
    Container m_dyn, m_dxpsi, m_dypsi;
    std::array<Container,2> m_v;
    // output relevant quantities
    Container m_UE2;
    std::map<std::string,Container> m_psi, m_gamma_n;

    std::vector<Container> m_multi_chi;

    dg::MultigridCG2d<Geometry, Matrix, Container> m_multigrid;
    std::vector<dg::PCG<Container> > m_multi_pcg;

    std::map<std::string, std::vector<dg::Helmholtz<Geometry, Matrix,
        Container>>> m_multi_helmN;
    std::vector<dg::Helmholtz<Geometry, Matrix, Container >> m_multi_helmP;
    std::vector<dg::Elliptic< Geometry, Matrix, Container >> m_multi_pol;

    std::map<std::string, std::function <void( const Container&, Container&)>> m_inv_helmN, m_inv_helmP;
    std::function< void( const Container&, Container&)> m_inv_pol;

    // initial guess
    dg::Extrapolation<Container> m_extra_phi;
    std::map<std::string, dg::Extrapolation<Container>> m_extra_n, m_extra_psi;

    // derivatives for each species
    std::map<std::string, dg::Elliptic< Geometry, Matrix, Container >> m_laplacianM;
    std::map<std::string, dg::Advection< Geometry, Matrix, Container>> m_adv;
    std::array<Matrix,2> m_centered;
    Parameters m_p;
    unsigned m_ncalls = 0;
    Container m_weights;
};

template< class Geometry, class Matrix, class Container>
Equations< Geometry, Matrix, Container>::Equations( const Geometry& grid, Parameters p) :
    m_binv( dg::evaluate( dg::LinearX( p.kappa, 1.), grid)),
    m_chi( dg::evaluate( dg::zero, grid )), m_omega(m_chi),
    m_tilden(m_chi), m_temp(m_chi),
    m_phi(m_chi),
    m_dyn(m_chi), m_dxpsi(m_chi), m_dypsi(m_chi),
    m_v({m_chi,m_chi}), m_UE2(m_chi),
    m_extra_phi( 2, m_omega),
    m_p(p)
{
    m_weights = dg::create::volume(grid);
    m_multigrid.construct( grid, m_p.num_stages);
    m_multi_chi = m_multigrid.project(m_chi);
    for( unsigned u=0; u<m_p.num_stages; u++)
    {
        m_multi_pcg.push_back({ m_multi_chi[u], 10000});
    }
    for( auto s : m_p.species)
    {
        m_psi[s] = m_chi;
        m_gamma_n[s] = m_chi;
        m_extra_n[s] = m_extra_psi[s] = m_extra_phi;
        m_laplacianM[s] = { grid, m_p.bcx.at(s), m_p.bcy.at(s), m_p.diff_dir};
        m_adv[s] = { grid,  m_p.bcx.at(s), m_p.bcy.at(s)};
        // construct inversion for Gamma operators
        for( unsigned u=0; u<m_p.num_stages; u++)
        {
            m_multi_helmN[s].push_back({ -0.5*m_p.tau.at(s)*m_p.mu.at(s), {m_multigrid.grid(u), m_p.bcx.at(s),
                    m_p.bcy.at(s), m_p.pol_dir}});
        }
        if( m_p.mu.at(s) == 0 || m_p.tau.at(s) == 0)
        {
            m_inv_helmN[s] = [](const auto& y, auto& x){ dg::blas1::copy( y, x);};
            m_inv_helmP[s] = [](const auto& y, auto& x){ dg::blas1::copy( y, x);};
        }
        else
        {
            m_inv_helmN[s] = [this, s=s, tilde=m_chi] ( const auto& y, auto& x) mutable
            {
                m_multigrid.set_benchmark( true, "Gamma N "+s);
                dg::blas1::copy( y, tilde);
                if( m_p.bcx.at(s) != dg::PER || m_p.bcy.at(s) != dg::PER)
                {
                    dg::blas1::plus( tilde, -1.);
                    dg::blas1::plus( x, -1.);
                }
                m_multigrid.solve( m_multi_helmN.at(s), x, tilde, m_p.eps_gamma);
                if( m_p.bcx.at(s) != dg::PER || m_p.bcy.at(s) != dg::PER)
                    dg::blas1::plus( x, +1.);
            };

            m_inv_helmP[s] = [this, s=s] ( const auto& y, auto& x)
            {
                m_multigrid.set_benchmark( true, "Gamma P "+s);
                for( unsigned u=0; u<m_p.num_stages; u++)
                    m_multi_helmP[u].alpha() = -0.5*m_p.tau.at(s)*m_p.mu.at(s);
                m_multigrid.solve( m_multi_helmP, x, y, m_p.eps_gamma);
            };
            //! now we store this as a variable! Delete copy and move!
        }
    }
    std::string s = "potential";
    m_centered = {dg::create::dx( grid, m_p.bcx.at(s)),
                  dg::create::dy( grid, m_p.bcy.at(s))};
    for( unsigned u=0; u<m_p.num_stages; u++)
    {
        m_multi_helmP.push_back( {0., {m_multigrid.grid(u), m_p.bcx.at(s),
                m_p.bcy.at(s), m_p.pol_dir}});
        m_multi_pol.push_back( {m_multigrid.grid(u), m_p.bcx.at(s), m_p.bcy.at(s),
                m_p.pol_dir});
    }
    m_inv_pol = [this] (const auto& y, auto& x)
    {
        m_multigrid.set_benchmark( true, "Polarisation");
        m_multigrid.solve( m_multi_pol, x, y, m_p.eps_pol);

    };
}


//idx is impurity species one or two
template< class G, class M, class Container>
void Equations<G, M, Container>::compute_psi( double t, const Container& potential)
{
    m_multi_pol[0].variation(m_binv, potential, m_UE2); // u_E^2
    for( auto s : m_p.species)
    {
        m_extra_psi.at(s).extrapolate( t, m_psi.at(s));
        dg::apply( m_inv_helmP.at(s), potential, m_psi.at(s));
        m_extra_psi.at(s).update( t, m_psi.at(s));

        dg::blas1::axpby( 1., m_psi.at(s), -0.5*m_p.mu.at(s), m_UE2, m_psi.at(s));
    }
}


template<class G, class Matrix, class Container>
void Equations<G, Matrix, Container>::compute_phi( double t, const std::map<std::string, Container>& y)
{
    // Compute chi
    dg::blas1::copy( m_p.epsilon_D, m_chi);
    dg::blas1::pointwiseDot( m_binv, m_binv, m_temp);
    for( auto s : m_p.species)
        dg::blas1::pointwiseDot( m_p.a.at(s)*m_p.mu.at(s), y.at(s), m_temp, 1., m_chi);

    // update multi_pol object
    m_multigrid.project( m_chi, m_multi_chi);
    for( unsigned u=0; u<m_p.num_stages; u++)
        m_multi_pol[u].set_chi( m_multi_chi[u]);

    //Compute rhs: \sum_s a_s Gamma_s N_s
    dg::blas1::copy( 0., m_omega);
    for( auto s : m_p.species)
    {
        // Solve Gamma N
        m_extra_n.at(s).extrapolate( t, m_gamma_n.at(s));
        dg::apply( m_inv_helmN.at(s), y.at(s), m_gamma_n.at(s));
        m_extra_n.at(s).update( t, m_gamma_n.at(s));
        dg::blas1::axpby( m_p.a.at(s), m_gamma_n.at(s), 1., m_omega);
    }

    // Solve for potential
    m_extra_phi.extrapolate( t, m_phi);
    dg::apply( m_inv_pol, m_omega, m_phi);
    m_extra_phi.update( t, m_phi);
}

template< class G, class M, class Container>
void Equations< G, M, Container>::operator()(double time, const std::map<std::string, Container>& y, std::map<std::string, Container>& yp)
{
    m_ncalls ++;
    // y.at(s) == n_s
    dg::Timer t;
    t.tic();


    // compute m_phi
    compute_phi( time, y);
    // compute all m_psi
    compute_psi( time, m_phi);
    t.toc();
    std::cout << "# Computing  phi and psi took [s] "<<t.diff()<<"\n";
    t.tic();


    for( auto s : m_p.species)
    {
        // we need to subtract 1 because dg expects 0 boundary conditions
        dg::blas1::copy( y.at(s), m_tilden);
        if( m_p.bcx.at(s) != dg::PER || m_p.bcy.at(s) != dg::PER)
            dg::blas1::plus( m_tilden, -1.);

        // ExB + Curv advection with updwind scheme
        dg::blas2::symv( m_centered[0], m_psi.at(s), m_dxpsi);
        dg::blas2::symv( m_centered[1], m_psi.at(s), m_dypsi);
        dg::blas1::pointwiseDot( -1., m_binv, m_dypsi, 0., m_v[0]);
        dg::blas1::pointwiseDot( +1., m_binv, m_dxpsi, 0., m_v[1]);
        dg::blas1::plus( m_v[1], -m_p.tau.at(s)*m_p.kappa);
        m_adv.at(s).upwind( -1., m_v[0], m_v[1], m_tilden, 0., yp.at(s));

        // Div ExB velocity
        dg::blas1::pointwiseDot( m_p.kappa, y.at(s), m_dypsi, 1., yp.at(s));

        // Add diffusion
        if( m_p.nu_perp.at(s) != 0)
        {
            dg::blas2::gemv( m_laplacianM.at(s), m_tilden, m_temp);
            dg::blas2::gemv( -m_p.nu_perp.at(s), m_laplacianM.at(s), m_temp, 1., yp.at(s));
        }
    }
    t.toc();
    std::cout << "# Computing rhs [s] took "<<t.diff()<<"\n";
}

}//namespace impurities

