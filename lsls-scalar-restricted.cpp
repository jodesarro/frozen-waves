#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <complex>
#include <future>
#include <vector>

// A library to evaluate Bessel functions available at https://github.com/jodesarro/bessel-library.
#include "../bessel-library/bessel-library.hpp" 

// A library to handle MTX files available at https://github.com/jodesarro/mtxdat-library.
#include "../mtxdat-library/mtxdat-library.hpp" 

// A library to handle WL files available at https://github.com/jodesarro/wldat-library.
#include "../wldat-library/wldat-library.hpp"

using namespace std;

// Imaginary unity.
static const complex<double> I = complex<double>(0., 1.);


// Function to calculate the field given x.
void field_for_x( int ix, int xpoints, int ypoints, double zmin, double dz, double zpoints,
    complex<double> * A, complex<double> * b, complex<double> * bj0, int inmax, int immax, complex<double> * field )
{
    for ( int iy=0; iy<ypoints; iy++ )
    {
        for ( int iz=0; iz<zpoints; iz++ )
        {
            double z = zmin + static_cast<double>(iz)*dz;
            complex<double> psi = complex<double> (0., 0.);
            for ( int im=0; im<immax; im++ )
            {
                for ( int in=0; in<inmax; in++ )
                {
                    psi += A[in + inmax*im]*bj0[in+inmax*im+inmax*immax*ix+inmax*immax*xpoints*iy]*exp( I*b[in]*z );
                }
            }
            field[ix + xpoints*iy + xpoints*ypoints*iz] = psi;
        }
    }
}


// Function to calculate Bessel functions given x.
void bj0_for_x( int ix, double x, int xpoints, double * x0, double ymin, double dy, int ypoints, double * y0,
    complex<double> * h, int inmax, int immax, complex<double> * bj0 )
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
                bj0[in+inmax*im+inmax*immax*ix+inmax*immax*xpoints*iy] = bessel::cyl_j(0, h[in]*varrho);
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
    cout << "| LSLS-SCALAR-RESTRICTED | 05 SEP 2023 | @JODESARRO |" << endl;
    cout << "A C++ routine to evaluate the scalar field of longitudinally structured light sheets restricted to the case where Qm=Q, Lm=L and Nm=N." << endl;
    cout << "https://github.com/jodesarro/frozen-waves\n" << endl;


    // Specification of general parameters of frozen waves.
    cout << "Loading general parameters..." << endl;
    complex<double> nref = complex<double> ( 1. , 0. );
    //double l0 = 451.0e-9; // B
    //double l0 = 532.0e-9; // G
    //double l0 = 635.0e-9; // R
    double l0 = 632.0e-9;
    double k0 = 2.0*M_PI/l0;
    static const int N = 20; // It is assumed Nm=N.
    static const int M = 101;
    cout << "Chosen value for the refractive index: " << nref << endl;
    cout << "Chosen value for the wavelength (in vacuum): " << l0 << " m" << endl;
    cout << "Calculated value for the angular frequency (in vacuum): " << k0*299792458 << " rad/s" << endl;
    cout << "Chosen value for N parameter: " << N << endl;
    cout << "Chosen value for M parameter: " << M << endl;   


    // Calculation of maximum values for indexes n and m of frozen waves.
    //  n->[in] : n=in-N, 0<=in<inmax
    //  m->[im] : m=im+1, 0<=im<immax
    static const int inmax = N*2+1;
    static const int immax = M;


    // Specification of L parameter.
    //  It is assumed Lm=L.
    double L = 0.060;
    cout << "Chosen value for L parameter: " << L << " m" << endl;   


    // Specification of Q parameter.
    //  It is assumed Qm=Q.

    //  Use the next 2 lines only to calculate Q with a chosen spot radius.
    //double chosen_spot_radius_fw = 9.0e-6;
    //double Q = sqrt( 1. - pow(2.4048*l0/(chosen_spot_radius_fw*2.*M_PI), 2.) )*2.*M_PI/l0;
    //double a = Q/real(nref)*k0;
    
    //  Use the next line only to specify a fixed value for Q.
    double a = 0.9986;
    double Q = a*real(nref)*k0;
    
    cout << "Chosen value for the a parameter: " << a << endl;   
    cout << "Calculated value for Q parameter: " << Q << " 1/m" << endl;


    // Specification of positions x'm and y'm relative to the origin.
    double dx0 = 50.e-6;
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


    // Specification of whether to ensure purely real transverse wavenumber.
    //  Set true to ensure purely real transverse wavenumber or false otherwise.
    bool is_purely_real_transverse_wavenumber = false;


    // Specification of data calculation parameters.
    double xmin = 0.;
    double xmax = double(M-1)*dx0;
    static const int xpoints = 1001;
    double ymin = 0.;
    double ymax = 0.;
    static const int ypoints = 1;
    double zmin = 0.;
    double zmax = L;
    static const int zpoints = 151;
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
    complex<double> h [inmax];
    complex<double> b [inmax];
    if ( is_purely_real_transverse_wavenumber)
    {
        for ( int in=0; in<inmax; in++ )
        {
            double n = static_cast<double>( in - N );
            double bnr = Q + 2.0*M_PI*n/L;
            double bni = k0*k0*real(nref)*imag(nref)/bnr;
            b[in] = complex<double> ( bnr, bni );
            h[in] = sqrt( (real(nref)*real(nref)-imag(nref)*imag(nref))*k0*k0 - bnr*bnr + bni*bni );
        }
    }
    else
    {
        for ( int in=0; in<inmax; in++ )
        {
            double n = static_cast<double>( in - N );
            double bnr = Q + 2.0*M_PI*n/L;
            b[in] = complex<double> ( bnr, bnr*imag(nref)/real(nref) );
            h[in] = sqrt( k*k - b[in]*b[in] );
        }
    }


    // Resistant attenuation or compensation parameter.
    //  The βI0m=βI0 is a consequence of Qm=Q.
    double bi0 = imag(b[N]);


    // Calculation of maximum value possible for N by requiring 0≤β0≤kr.
    //  The Nmaxm=Nmax is consequence of Nm=N, Qm=Q and Lm=L.
    int Nmax = floor( L * min( k0*real(nref) - Q, Q ) / (2.0*M_PI) );
    cout << "Calculated maximum value possible for N: " << Nmax << endl;


    // Calculation of axicon angles for generating Bessel beams inside the considered media.
    //  The θnm=θn is a consequence of βnm=βn.
    double theta_pN = acos( real(b[inmax-1])/(k0*real(nref)) );
    cout << "Calculated Bessel beam axicon angle for n=+N: ";
    cout << theta_pN << " rad = " << (360./(2.0*M_PI))*theta_pN << " degree" << endl;
    double theta_mN = acos( real(b[0])/(k0*real(nref)) );
    cout << "Calculated Bessel beam axicon angle for n=-N: ";
    cout << theta_mN << " rad = " << (360./(2.0*M_PI))*theta_mN << " degree" << endl;


    // Calculation of aperture radius for experimental generation of Bessel beams
    //  inside the considered media.
    //  The Rnm=Rn is a consequence of θnm=θn.
    double calculated_aperture_bb_pN = L*tan(theta_pN);
    cout << "Calculated Bessel beam aperture radius for n=+N: ";
    cout << calculated_aperture_bb_pN << " m" << endl;
    double calculated_aperture_bb_mN = L*tan(theta_mN);
    cout << "Calculated Bessel beam aperture radius for n=-N: ";
    cout << calculated_aperture_bb_mN << " m" << endl;


    // Calculation of spot radius of a ideal linear frozen wave.
    //  That is, spot radius of the ideal Bessel beam with n=0.
    double calculated_spot_radius_fw = 2.405/real(h[N]);
    cout << "Calculated single FW spot radius: ";
    cout << calculated_spot_radius_fw << " m" << endl;


    // Calculation of spot radius of a ideal linear frozen wave by asymptotic expansion.
    //  That is, spot radius of the ideal Bessel beam with n=0 by asymptotic expansion.
    double calculated_asymptotic_spot_radius_fw = 0.25*3.0*M_PI/real(h[N]);
    cout << "Calculated single FW spot radius by asymptotic expansion: ";
    cout << calculated_asymptotic_spot_radius_fw << " m" << endl;


    // Calculation of parameters related to complex refractive index.
    if ( imag(nref) != 0. )
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
    }


    // Reading MTX function-f file to a function F[ii+iimax*ik].
    //  x->[ii] : 0<=ii<iimax
    //  z->[ik] : 0<=ik<ikmax
    cout << "Loading file containing the function F..." << endl;
    int iimax;
    int ikmax;
    mtxdat_getsize( "function-f.mtx", iimax, ikmax );
    double * F = new double [iimax*ikmax];
    mtxdat_import( "function-f.mtx", F, iimax, ikmax );


    // Calculation of A coefficients by means of an approximation of the integral
    //  by trapezoidal method for equally spaced z.
    //  The Anm = An is a consequence of Lm=L, Qm=Q, Nm=N.
    cout << "Calculating the A coefficients..." << endl;
    complex<double> A [inmax*immax];
    for ( int im=0; im<immax; im++ )
    {
        for ( int in=0; in<inmax; in++ )
        {
            double n = static_cast<double>(in - N);
            int ii = floor( static_cast<double>(im)*static_cast<double>(iimax-1)/static_cast<double>(immax-1) );
            complex<double> aux = 0.5*( F[ii + iimax*(ikmax-1)]*exp( bi0*L ) + F[ii + 0] );
            for ( int ik=1; ik<ikmax-1; ik++ )
            {
                double z = L*static_cast<double>(ik)/static_cast<double>(ikmax-1);
                aux += F[ii + iimax*ik]*exp( -I*2.0*M_PI*n*z/L )*exp( bi0*z );
            }
            A[in + inmax*im] = aux / static_cast<double>(ikmax-1);
        }
    }
    delete[] F;


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
            ix, x, xpoints, x0, ymin, dy, ypoints, y0, h, inmax, immax, bj0) );
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
            ix, xpoints, ypoints, zmin, dz, zpoints, A, b, bj0, inmax, immax, field ) );
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
