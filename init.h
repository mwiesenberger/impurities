#pragma once

namespace imp
{

std::vector<dg::x::DVec> initial_conditions(
    const dg::x::CartesianGrid2d& grid,
    dg::file::WrappedJsonValue init)
{
    std::vector<dg::x::HVec> y0(3, dg::x::HVec( grid.size()) );
    std::string initial = init.get("type", "lamb").asString();
    double amp, sigma, posX, posY;
        vorticity = js["vorticity"].asDouble();
        mode = js["mode"].asUInt();
        wall_pos = js["wall_pos"].asDouble();
        wall_amp = js["wall_amp"].asDouble();
        wall_sigma = js["wall_sigma"].asDouble();
        amp = js["amplitude"].asDouble();
        sigma = js["sigma"].asDouble();
        posX = js["posX"].asDouble();
        posY = js["posY"].asDouble();
    if( initial == "zero")
    {
        omega = dg::evaluate( dg::zero, grid);
    }

    if( "gaussian" == initial)
    {
        if( p.vorticity == 0)
        {
            gamma.alpha() = -0.5*p.tau[1];
            y0[0] = dg::evaluate( gaussian, grid);
            dg::blas2::symv( gamma, y0[0], y0[1]); // n_e = \Gamma_i n_i -> n_i = ( 1+alphaDelta) n_e' + 1 

            dg::blas1::scal( y0[1], 1./p.a[1]); //n_i ~1./a_i n_e
            y0[2] = dg::evaluate( dg::zero, grid);
        }
        else
        {
            y0[1] = y0[0] = dg::evaluate( gaussian, grid);
            dg::blas1::scal( y0[1], 1/p.a[1]);
            y0[2] = dg::evaluate( dg::zero, grid);
        }
    }
    if( initial == "wall")
    {
        //init wall in y0[2]
        dg::GaussianX wall( p.wall_pos*grid.lx(), p.wall_sigma, p.wall_amp); 
        dg::DVec wallv = dg::evaluate( wall, grid);
        gamma.alpha() = -0.5*p.tau[2]*p.mu[2];
        dg::blas2::symv( gamma, wallv, y0[2]); 
        if( p.a[2] != 0.)
            dg::blas1::scal( y0[2], 1./p.a[2]); //n_z ~1./a_z

        //init blob in y0[1]
        gamma.alpha() = -0.5*p.tau[1];
        y0[0] = dg::evaluate( gaussian, grid);
        dg::blas2::symv( gamma, y0[0], y0[1]); 
        if( p.a[2] == 1)
        {
            std::cerr << "No blob with trace ions possible!\n";
            return -1;
        }
        dg::blas1::scal( y0[1], 1./p.a[1]); //n_i ~1./a_i n_e

        //sum up
        if( p.a[2] != 0)
            dg::blas1::axpby( 1., wallv, 1., y0[0]); //add wall to blob in n_e
    }
    if( init == "blob") 
    {
        gamma.alpha() = -0.5*p.tau[2]*p.mu[2];
        y0[0] = dg::evaluate( gaussian, grid);
        dg::blas2::symv( gamma, y0[0], y0[2]); 
        if( p.a[2] == 0)
        {
            std::cerr << "No impurity blob with trace impurities possible!\n";
            return -1;
        }
        dg::blas1::axpby( 1./p.a[2], y0[2], 0., y0[2]); //n_z ~1./a_z n_e
        y0[1] = dg::evaluate( dg::zero, grid);
    }
    throw dg::Error( dg::Message() << "Initial condition "
                    <<initial<<" not recognized!");
    return std::vector<dg::x::DVec>>(omega);
}

} // namespace imp
