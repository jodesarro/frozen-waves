#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <complex>
#include <future>
#include <vector>

// A library to evaluate Bessel functions available at https://github.com/jodesarro/bessel-library.
#include "bessel-library.hpp" 

// A library to handle MTX files available at https://github.com/jodesarro/data-file-library.
#include "mtxdat-library.hpp" 

// A library to handle WL files available at https://github.com/jodesarro/data-file-library.
#include "wldat-library.hpp"

using namespace std;

// Imaginary unity.
static const complex<double> I = complex<double>(0., 1.);



// Function to calculate the field given x.
void field_for_x( int ix, int xpoints, int ypoints, double zmin, double dz, double zpoints, double * z0_p,
    complex<double> * A_p, complex<double> * b_p, complex<double> * bj0, int inmax, int immax, int ipmax, complex<double> * field )
{
    for ( int iy=0; iy<ypoints; iy++ )
    {
        for ( int iz=0; iz<zpoints; iz++ )
        {
            double z = zmin + static_cast<double>(iz)*dz;
            int ip = ipmax-1;
            complex<double> psi_p = complex<double> (0., 0.);
            for ( int im=0; im<immax; im++ )
            {
                for ( int in=0; in<inmax; in++ )
                {
                    psi_p += A_p[in + inmax*im + inmax*immax*ip]*bj0[in + inmax*im + inmax*immax*ix + inmax*immax*xpoints*iy]*exp( I*b_p[in + inmax*ip]*(z-z0_p[ip]) );
                }
            }
            field[ix + xpoints*iy + xpoints*ypoints*iz] = psi_p;
        }
    }
}


// Function to calculate Bessel functions given x.
void bj0_for_x( int ix, double x, int xpoints, double * x0, double ymin, double dy, int ypoints, double * y0,
    complex<double> * h_p, int inmax, int immax, int ipmax, complex<double> * bj0 )
{
    for ( int iy=0; iy<ypoints; iy++ )
    {
        double y = ymin + static_cast<double>(iy)*dy;
        for ( int im=0; im<immax; im++ )
        {
            double square_of_diff_x_x0 = (x-x0[im])*(x-x0[im]);
            double square_of_diff_y_y0 = (y-y0[im])*(y-y0[im]);
            double varrho = sqrt( square_of_diff_x_x0 + square_of_diff_y_y0 );
            for ( int in=0; in<inmax; in++ )
            {
                int ip = ipmax-1;
                bj0[in + inmax*im + inmax*immax*ix + inmax*immax*xpoints*iy] = bessel::cyl_j(0, h_p[in+inmax*ip]*varrho);
            }
        }
    }
}


int main()
{
    // Starting time counter.
    long int start_time = time( nullptr );
    

    // Redirecting all cout to a text file
    ofstream cout("output.txt");


    // Printing information.
    cout << "| LSLS-SCALAR-RESTRICTED-P-4f | 07 APR 2024 | @JODESARRO |" << endl;
    cout << "A C++ routine to evaluate the scalar field of longitudinally structured light sheets, restricted to the case where Qm=Q, Lm=L and Nm=N, after P 4f systems." << endl;
    cout << "https://github.com/jodesarro/frozen-waves\n" << endl;


    // !EXPRESSIONS NOT CHECKED FOR COMPLEX REFRACTIVE INDEX!

    
    // Constants definition
    const double none = 0.;


    // Specification of general parameters of frozen waves.
    cout << "Loading general parameters..." << endl;
    complex<double> nref = complex<double> ( 1. , 0. ); 
    double l0 = 514.0e-9;
    double k0 = 2.0*M_PI/l0;
    static const int N = 27; // It is assumed Nm=N.
    static const int M = 1;
    double z0 = 0.1;
    cout << "Chosen value for the refractive index: " << nref << endl;
    cout << "Chosen value for the wavelength (in vacuum): " << l0 << " m" << endl;
    cout << "Calculated value for the angular frequency (in vacuum): " << k0*299792458 << " rad/s" << endl;
    cout << "Chosen value for N parameter: " << N << endl;
    cout << "Chosen value for M parameter: " << M << endl;
    cout << "Chosen value for z0 parameter: " << z0 << endl;   


    // Calculation of maximum values for indexes n and m of frozen waves.
    //  n->[in] : n=in-N, 0<=in<inmax
    //  m->[im] : m=im+1, 0<=im<immax
    static const int inmax = N*2+1;
    static const int immax = M;


    // Specification of L parameter.
    //  It is assumed Lm=L.
    double L = 2.;
    cout << "Chosen value for L parameter: " << L << " m" << endl;


    // Specification of 4f systems parameters

    // Number of 4f-systems
    static const int P = 3;

    // Calculation of maximum values for indexes n and m of frozen waves.
    //  p->[ip] : p=ip, 0<=ip<ipmax
    static const int ipmax = P+1;

    // Lenses parameters
    //  In the absense of 4f systems, set f1_p[p]=1., f2_p[p]=-f1_p[p] and d_p[p]=0 for all p.
    double d_p [ipmax] = { 0.01, 0.01, 0.01, none }; // d_p[ipmax-1] may be any number, it will not be used.
    double f1_p [ipmax] = { none, 0.15, 0.15, 0.15 }; // f1_p[0] may be any number, it will not be used.
    double f2_p [ipmax] = { none, 0.025, 0.05, 0.025  }; // f2_p[0] may be any number, it will not be used.

    double z0_p [ipmax]; z0_p[0] = z0;
    double S_p [ipmax]; S_p[0] = none; // S_p[0] will not be used.
    double M_p [ipmax];  M_p[0] = none; // M_p[0] will not be used.
    for (int ip=1; ip<ipmax; ip++)
    {
        S_p[ip] = f2_p[ip] + f1_p[ip];
        M_p[ip] = -f2_p[ip]/f1_p[ip];
        z0_p[ip] = (z0_p[ip-1]-d_p[ip-1])*M_p[ip]*M_p[ip]-S_p[ip]*M_p[ip];
    }

    complex<double> Phi_p [ipmax]; Phi_p[0] = none; // Phi_p[0] will not be used.
    double total_reescalation_transverse = 1.;
    double total_reescalation_longitudinal = 1.;
    double total_reescalation_amplitude = 1.;
    double total_longitudinal_displacement = z0_p[ipmax-1];
    for (int ip=1; ip<ipmax; ip++)
    {
        Phi_p[ip] = k0*nref*(S_p[ip]*(1.-1./M_p[ip])+z0_p[ip]*(1.-1./(M_p[ip]*M_p[ip])));
        total_reescalation_transverse *= abs(M_p[ip]);
        total_reescalation_longitudinal *= M_p[ip]*M_p[ip];
        total_reescalation_amplitude *= (1./M_p[ip])*(1./M_p[ip]);
        total_longitudinal_displacement += (S_p[ip] + d_p[ip-1]);
    }

    cout << "Calculated value for total transverse reescalation: " << total_reescalation_transverse << endl;
    cout << "Calculated value for total longitudinal reescalation: " << total_reescalation_longitudinal << endl;
    cout << "Calculated value for total amplitude reescalation: " << total_reescalation_amplitude << endl;
    cout << "Calculated value for total longitudinal displacement: " << total_longitudinal_displacement << endl;
    cout << "Calculated value for longitudinal displacement from last lens: " << z0_p[ipmax-1] << endl;


    // Specification of Q parameter.
    //  It is assumed Qm=Q.

    //  Use the next 5 lines only to calculate Q with a chosen spot radius for the final beam.
    //double chosen_final_spot_radius_fw = 1.75258e-05;
    //cout << "Chosen final spot radius of a FW: ";
    //cout << chosen_final_spot_radius_fw << " m" << endl;
    //double a = 1.0 - 0.5*pow( (2.4048*l0)/( 2.0*M_PI*real(nref)*chosen_final_spot_radius_fw/total_reescalation_transverse ), 2.0 );
    //cout << "Calculated value for the a parameter: " << a << endl;

    //  Use the next 2 lines only to specify a fixed value for Q.
    double a = 0.999993;
    cout << "Chosen value for the a parameter: " << a << endl;

    double Q = a*k0*real(nref);
    cout << "Calculated value for Q parameter: " << Q << " 1/m" << endl;


    // Specification of positions x'm and y'm relative to the origin.
    double dx0 = 0.;
    double dy0 = 0.;
    double x0 [immax];
    double y0 [immax];
    for ( int im=0; im<immax; im++ )
    {
        double m = static_cast<double>(im + 1);
        x0[im] = (m-1.0)*dx0;
        y0[im] = 0.;
    }
    cout << "Separation between consecutives FWs: (" << dx0 << ", " << dy0 << ") m" << endl;


    // Specification of data calculation parameters.
    double xmin = -total_reescalation_transverse*150.e-6;
    double xmax = -xmin;
    static const int xpoints = 201;
    double ymin = xmin;
    double ymax = xmax;
    static const int ypoints = 201;
    double zmin = z0_p[ipmax-1]; // z-axis With origin at the last lens.
    double zmax = z0_p[ipmax-1] + 0.35*total_reescalation_longitudinal;
    static const int zpoints = 201;
    cout << "x E (" << xmin << " m, " << xmax << " m, " << xpoints << " pt)" << endl;
    cout << "y E (" << ymin << " m, " << ymax << " m, " << ypoints << " pt)" << endl;
    cout << "z E (" << zmin << " m, " << zmax << " m, " << zpoints << " pt)" << endl;


    // Specification of async or lazy calculation of xy planes.
    bool is_async = true;
    auto launch_policy = (is_async) ? (launch::async) : (launch::deferred) ;


    // Usually there is no need to change values in all of the following.


    // Calculation of wavenumbers.
    //  The hnm=hn and βnm=βn are consequences of Nm=N, Qm=Q and Lm=L.
    cout << "Calculating wavenumbers..." << endl;
    complex<double> k = k0*nref;
    complex<double> h_p [inmax*ipmax];
    complex<double> b_p [inmax*ipmax];
    for ( int in=0; in<inmax; in++ )
    {
        double n = static_cast<double>( in - N );
        double bnr = Q + 2.0*M_PI*n/L;
        b_p[in + inmax*0] = complex<double> ( bnr, bnr*imag(nref)/real(nref) );
        h_p[in + inmax*0] = sqrt(2.)*k*sqrt(1.-b_p[in+inmax*0]/k);
        for (int ip=1; ip<ipmax; ip++)
        {
            h_p[in + inmax*ip] = h_p[in + inmax*(ip-1)]/M_p[ip];
            b_p[in + inmax*ip] = k-h_p[in + inmax*ip]*h_p[in + inmax*ip]/(2.*k);
        }
    }


    // Resistant attenuation or compensation parameter.
    //  The βI0m=βI0 is a consequence of Qm=Q.
    double bi0 = imag(b_p[N + inmax*0]);


    // Calculation of maximum value possible for N by requiring 0≤β0≤kr.
    //  The Nmaxm=Nmax is consequence of Nm=N, Qm=Q and Lm=L.
    int Nmax = floor( L * min( k0*real(nref) - Q, Q ) / (2.0*M_PI) );
    cout << "Calculated maximum value possible for N: " << Nmax << endl;


    // Calculation of axicon angles for generating incident Bessel beams inside the considered media.
    //  The θnm=θn is a consequence of βnm=βn.
    double theta_pN = acos( real(b_p[(inmax-1) + inmax*0])/(k0*real(nref)) );
    cout << "Calculated incident Bessel beam axicon angle for n=+N: ";
    cout << theta_pN << " rad = " << (360./(2.0*M_PI))*theta_pN << " degree" << endl;
    double theta_mN = acos( real(b_p[0 + inmax*0])/(k0*real(nref)) );
    cout << "Calculated incident Bessel beam axicon angle for n=-N: ";
    cout << theta_mN << " rad = " << (360./(2.0*M_PI))*theta_mN << " degree" << endl;


    // Calculation of aperture radius for experimental generation of incident Bessel beams
    //  inside the considered media.
    //  The Rnm=Rn is a consequence of θnm=θn.
    double calculated_aperture_bb_pN = L*tan(theta_pN);
    cout << "Calculated incident Bessel beam aperture radius for n=+N: ";
    cout << calculated_aperture_bb_pN << " m" << endl;
    double calculated_aperture_bb_mN = L*tan(theta_mN);
    cout << "Calculated incident Bessel beam aperture radius for n=-N: ";
    cout << calculated_aperture_bb_mN << " m" << endl;


    // Calculation of spot radius of incident, ideal linear frozen wave.
    //  That is, spot radius of the ideal Bessel beam with n=0.
    double calculated_spot_radius_fw = 2.4048/real(h_p[N + inmax*0]);
    cout << "Calculated spot radius of a incident FW: ";
    cout << calculated_spot_radius_fw << " m" << endl;


    // Calculation of spot radius of the incident, ideal linear frozen wave by asymptotic expansion.
    //  That is, spot radius of the ideal Bessel beam with n=0 by asymptotic expansion.
    double calculated_asymptotic_spot_radius_fw = 0.25*3.0*M_PI/real(h_p[N + inmax*0]);
    cout << "Calculated spot radius by asymptotic expansion of a incident FW: ";
    cout << calculated_asymptotic_spot_radius_fw << " m" << endl;


    // Calculation of spot radius of final, ideal linear frozen wave.
    //  That is, spot radius of the ideal Bessel beam with n=0.
    double calculated_final_spot_radius_fw_ = calculated_spot_radius_fw*total_reescalation_transverse;
    cout << "Calculated spot radius of a final FW: ";
    cout << calculated_final_spot_radius_fw_ << " m" << endl;


    // Calculation of parameters related to complex refractive index.
    /*if ( imag(nref) != 0. )
    {
        // Calculation of Δ parameter.
        //  The Δm=Δ is a consequence of Nm=N, Qm=Q and Lm=L.
        double delta = abs( imag(b[inmax-1]) - imag(b[0]) )/imag(b[N]);
        cout << "Calculated delta parameter: " << delta << endl;


        // Calculation of penetration depth of Bessel beams.
        //  The δnm=δn is a consequence of βnm=βn.
        double del_pN = 0.5/imag(b[inmax-1]);
        cout << "Calculated penetration depth for n=+N: " << del_pN << " m" << endl;
        double del_mN = 0.5/imag(b[0]);
        cout << "Calculated penetration depth for n=-N: " << del_mN << " m" << endl;


        if (!is_purely_real_transverse_wavenumber)
        {
            // Calculation of maximum Bessel beam aperture radius possible for an experimental generation
            //  inside the considered media due to the complex transversal wavenumber.
            //  The Rnm=Rn is a consequence of hnm=hn.
            double calculated_max_aperture_possible_bb_pN = 0.5/imag(h[inmax-1]);
            cout << "Calculated maximum aperture radius possible, due to complex transversal wavenumber, for n=+N: ";
            cout << calculated_max_aperture_possible_bb_pN << " m" << endl;
            double calculated_max_aperture_possible_bb_mN = 0.5/imag(h[0]);
            cout << "Calculated maximum aperture radius possible, due to complex transversal wavenumber, for n=-N: ";
            cout << calculated_max_aperture_possible_bb_mN << " m" << endl;


            // Calculation of minimum Bessel beam aperture radius possible for an experimental generation
            //  inside the considered media due to the complex transversal wavenumber.
            //  The Rnm=Rn is a consequence of hnm=hn.
            double calculated_min_aperture_possible_bb_pN = 0.25*3.0*M_PI/real(h[inmax-1]);
            cout << "Calculated minimum aperture radius possible, due to complex transversal wavenumber, for n=+N: ";
            cout << calculated_min_aperture_possible_bb_pN << " m" << endl;
            double calculated_min_aperture_possible_bb_mN = 0.25*3.0*M_PI/real(h[0]);
            cout << "Calculated minimum aperture radius possible, due to complex transversal wavenumber, for n=-N: ";
            cout << calculated_min_aperture_possible_bb_mN << " m" << endl;


            // Calculation of ratio nI/nR<2/3π that ensure finite behavior of Bessel beams of complex argument.
            double calculated_ratio_ni_nr = imag(nref)/real(nref);
            cout << "Calculated ratio nI/nR to ensure finite behavior of Bessel beams of complex argument: " << calculated_ratio_ni_nr << " < " << 2.0/(3.0*M_PI) << endl;
        }
    }*/

   
    // The start and end positions of the stepfunctions of F(z).
    double bound_l1 = 0.10;
    double bound_l2 = 0.15;
    double bound_l3 = 0.20;
    double bound_l4 = 0.25;


    // Calculation of A coefficients by means of an approximation of the integral
    //  by trapezoidal method for equally spaced z.
    //  The Anm = An is a consequence of Lm=L, Qm=Q, Nm=N.
    cout << "Calculating the A coefficients..." << endl;
    complex<double> A_p [inmax*immax*ipmax];
    for ( int im=0; im<immax; im++ )
    {
        for ( int in=0; in<inmax; in++ )
        {
            complex<double> int_l2l1;
            complex<double> int_l4l3;
            if ( in == N )
            {
                if ( bi0 == 0. )
                {
                    int_l2l1 = bound_l2-bound_l1;
                    int_l4l3 = bound_l4-bound_l3;
                }
                else
                {
                    int_l2l1 = ( exp(bi0*bound_l2) - exp(bi0*bound_l1) )/( bi0 );
                    int_l4l3 = ( exp(bi0*bound_l4) - exp(bi0*bound_l3) )/( bi0 );
                }
            }
            else
            {
                double n = static_cast<double>(in - N);
                int_l2l1 = ( exp((bi0-I*2.0*M_PI*n/L)*bound_l2) - exp((bi0-I*2.0*M_PI*n/L)*bound_l1) )/( bi0-2.0*I*M_PI*n/L );
                int_l4l3 = ( exp((bi0-I*2.0*M_PI*n/L)*bound_l4) - exp((bi0-I*2.0*M_PI*n/L)*bound_l3) )/( bi0-2.0*I*M_PI*n/L );
            }
            A_p[in + inmax*im + inmax*immax*0] = (int_l2l1+int_l4l3)/L;
            for ( int ip=1; ip<ipmax; ip++)
            {
                A_p[in + inmax*im + inmax*immax*ip] = exp(-I*Phi_p[ip])*A_p[in + inmax*im + inmax*immax*(ip-1)]/(M_p[ip]);
            }
        }
    }


    // Calculation of spatial increments.
    double dx, dy, dz;
    dx = (xpoints==1) ? (0.) : (xmax-xmin)/(xpoints-1) ;
    dy = (ypoints==1) ? (0.) : (ymax-ymin)/(ypoints-1) ;
    dz = (zpoints==1) ? (0.) : (zmax-zmin)/(zpoints-1) ;


    // Calculation of Bessel functions J0.
    cout << "Calculating Bessel functions..." << endl;
    complex<double> * bj0 = new complex<double> [inmax*immax*xpoints*ypoints];
    std::vector<std::future<void>> future_bj0_for_x;
    for ( int ix=0; ix<xpoints; ix++ )
    {
        double x = xmin + static_cast<double>(ix)*dx;
        future_bj0_for_x.push_back( async(launch_policy, bj0_for_x,
            ix, x, xpoints, x0, ymin, dy, ypoints, y0, h_p, inmax, immax, ipmax, bj0) );
    }
    for(auto &e : future_bj0_for_x)
    {
        e.get();
    }
    

    // Calculation of scalar field.
    cout << "Calculating scalar field..." << endl;
    complex<double> * field = new complex<double> [xpoints*ypoints*zpoints];
    std::vector<std::future<void>> future_field_for_x;
    for ( int ix=0; ix<xpoints; ix++ )
    {
        double x = xmin + static_cast<double>(ix)*dx;
        future_field_for_x.push_back( async( launch_policy, field_for_x,
            ix, xpoints, ypoints, zmin, dz, zpoints, z0_p, A_p, b_p, bj0, inmax, immax, ipmax, field ) );
    }
    for(auto &e : future_field_for_x)
    {
        e.get();
    }
    delete[] bj0;


    // Writing field to WL file.
    cout << "Exporting scalar field..." << endl;
    wldat_export( "field.wl", field, xpoints, ypoints, zpoints, 16, true );
    delete[] field;


    // Printing time counted.
    cout << "Done! Elapsed time: " << ( time(nullptr) - start_time ) << " seconds.";


    puts("\a");
    system("pause");
}