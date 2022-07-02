// Cited references are available at https://github.com/jodesarro/frozen-waves.

#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <complex>

// A library to evaluate Bessel functions available at https://github.com/jodesarro/bessel-library.
#include "../bessel-library/bessel-library.hpp"

// A library to handle MTX files available at https://github.com/jodesarro/mtxdat-library.
#include "../mtxdat-library/mtxdat-library.hpp"

// A library to handle WL files available at https://github.com/jodesarro/wldat-library.
#include "../wldat-library/wldat-library.hpp"

using namespace std;

int main()
{
    // Starting time counter.
    long int start_time = time( nullptr );


    // Printing information.
    puts("| LFW-SCALAR-TRADITIONAL | 02 JUL 2022 | @JODESARRO |");
    puts("A C++ routine to evaluate the scalar field of linear frozen waves using the traditional method.");
    puts("https://github.com/jodesarro/frozen-waves\n");


    // Specification of general parameters of frozen waves.
    puts("Loading general parameters...");
    complex<double> nref = complex<double> ( 1.5 , 0.49e-6 );
    double l0 = 308.0e-9;
    double k0 = 2.0*M_PI/l0;
    complex<double> k = k0*nref;
    static const int N = 20;


    // Calculation of maximum value for index n of frozen waves.
    //  n->[in] : n=in-N, 0<=in<inmax
    static const int inmax = N*2+1;


    // Specification of L parameter.
    double L = 0.33;


    // Specification of Q parameter.
    double Q = 0.9999*real(nref)*k0;
    cout << "A chosen value for Q parameter: " << Q << " 1/m" << endl;


    // Calculation of maximum value possible for N by requiring 0≤β0≤kr.
    int Nmax = floor( L * min( k0*real(nref) - Q, Q ) / (2.0*M_PI) );
    cout << "Calculated maximum value possible for N: " << Nmax << endl;


    // Calculation of wavenumbers.
    puts("Calculating wavenumbers...");
    complex<double> h [inmax];
    complex<double> b [inmax];
    for ( int in=0; in<inmax; in++ )
    {
        double n = static_cast<double>( in - N );
        double bnr = Q + 2.0*M_PI*n/L;
        b[in] = complex<double> ( bnr, bnr*imag(nref)/real(nref) );
        h[in] = sqrt( k*k - b[in]*b[in] );
    }


    // Calculation of Δ parameter.
    double delta = abs( imag(b[inmax-1]) - imag(b[0]) )/imag(b[N]);
    cout << "Calculated delta parameter: " << delta << endl;


    // Calculation of penetration depth of Bessel beams.
    double del_min = 0.5/imag(b[inmax-1]);
    cout << "Calculated penetration depth for n=+N: " << del_min << " m" << endl;
    double del_max = 0.5/imag(b[0]);
    cout << "Calculated penetration depth for n=-N: " << del_max << " m" << endl;


    // Calculation of axicon angles for generating Bessel beams inside the considered media.
    double theta_min = acos( real(b[inmax-1])/(k0*real(nref)) );
    cout << "Calculated axicon angle for n=+N: ";
    cout << theta_min << " rad = " << (360./(2.0*M_PI))*theta_min << " degree" << endl;
    double theta_max = acos( real(b[0])/(k0*real(nref)) );
    cout << "Calculated axicon angle for n=-N: ";
    cout << theta_max << " rad = " << (360./(2.0*M_PI))*theta_max << " degree" << endl;


    // Calculation of aperture radius for experimental generation of Bessel beams
    //  inside the considered media.
    double calculated_aperture_bb_min = L*tan(theta_min);
    cout << "Calculated aperture radius for n=+N: ";
    cout << calculated_aperture_bb_min << " m" << endl;
    double calculated_aperture_bb_max = L*tan(theta_max);
    cout << "Calculated aperture radius for n=-N: ";
    cout << calculated_aperture_bb_max << " m" << endl;


    // Calculation of parameters related to Bessel functions of complex argument.
    if ( imag(nref) != 0. )
    {
        // Calculation of maximum Bessel beam aperture radius possible for an experimental generation
        //  inside the considered media due to the complex transversal wavenumber.
        double calculated_max_aperture_possible_bb = 0.5/imag(h[0]);
        cout << "Calculated maximum aperture radius possible: ";
        cout << calculated_max_aperture_possible_bb << " m" << endl;


        // Calculation of ratio nI/nR<2/3π that ensure finite behavior of Bessel beams of complex argument.
        double calculated_ratio_ni_nr = imag(nref)/real(nref);
        cout << "Calculated ratio nI/nR: ";
        cout << calculated_ratio_ni_nr << " < " << 2.0/(3.0*M_PI) << endl;
    }


    // Calculation of spot radius of a ideal linear frozen wave.
    //  That is, spot radius of the ideal Bessel beam with n=0.
    double calculated_spot_radius_lfw = 2.405/real(h[N]);
    cout << "Calculated LFW spot radius: ";
    cout << calculated_spot_radius_lfw << " m" << endl;


    // Calculation of spot radius of a ideal linear frozen wave by asymptotic expansion.
    //  That is, spot radius of the ideal Bessel beam with n=0 by asymptotic expansion.
    double calculated_asymptotic_spot_radius_lfw = 0.25*3.0*M_PI/real(h[N]);
    cout << "Calculated LFW spot radius by asymptotic expansion: ";
    cout << calculated_asymptotic_spot_radius_lfw << " m" << endl;


    // Specification of resistant (attenuation or compensation) parameter as defined in [1].
    //  Set bi0=imag(b[N]) for attenuation-resistant method
    //  or bi0=0. for non-attenuation-resistant method.
    double bi0 = imag(b[N]);


    // Specification of positions x0 and y0 relative to the origin.
    //  In [1] the authors assumed x0=y0=0, i.e, the beam is trasversally located at the origin.
    //  Note that, in cylindrical coordinates, ρ0²=x0²+y0² and φ0=atan(y0/x0).
    double x0 = 0.;
    double y0 = 0.;


    // Specification of data calculation parameters.
    double xmin = -0.004;
    double xmax = 0.004;
    static const int xpoints = 151;
    double ymin = 0.;
    double ymax = 0.;
    static const int ypoints = 1;
    double zmin = 0.;
    double zmax = L;
    static const int zpoints = 150;
    cout << "x E (" << xmin << " m, " << xmax << " m, " << xpoints << " pt)" << endl;
    cout << "y E (" << ymin << " m, " << ymax << " m, " << ypoints << " pt)" << endl;
    cout << "z E (" << zmin << " m, " << zmax << " m, " << zpoints << " pt)" << endl;


    // Usually there is no need to change values in all of the following.


    // Reading MTX intensity profile file to function F[ik].
    //  z->[ik] : 0<=ik<ikmax
    puts("Loading file containing the intensity profile...");
    int iimax; // iimax==1.
    int ikmax;

    mtxdat_getsize( "intensity-profile.mtx", iimax, ikmax );
    double * F = new double [iimax*ikmax];
    mtxdat_import( "intensity-profile.mtx", F, iimax, ikmax );


    // Calculation of A coefficients by means of an approximation of the integral
    //  by trapezoidal method for equally spaced z.
    puts("Calculating the A coefficients...");
    complex<double> A [inmax];
    for ( int in=0; in<inmax; in++ )
    {
        double n = static_cast<double>(in - N);
        complex<double> aux = 0.5*( F[ikmax-1]*exp( bi0*L ) + F[0] );
        for ( int ik=1; ik<ikmax-1; ik++ )
        {
            double z = L*static_cast<double>(ik)/static_cast<double>(ikmax-1);
            aux += F[ik]*exp( complex<double>(0., -2.0*M_PI*n*z/L) )*exp( bi0*z );
        }
        A[in] = aux / static_cast<double>(ikmax-1);
    }
    delete[] F;


    // Calculation of scalar field.
    puts("Calculating scalar field...");
    complex<double> * field = new complex<double> [xpoints*ypoints*zpoints];
    double dx, dy, dz;
    dx = (xpoints==1) ? (0.) : (xmax-xmin)/(xpoints-1) ;
    dy = (ypoints==1) ? (0.) : (ymax-ymin)/(ypoints-1) ;
    dz = (zpoints==1) ? (0.) : (zmax-zmin)/(zpoints-1) ;
    for ( int iz=0; iz<zpoints; iz++ )
    {
        // Printing progress.
        cout << fixed << setprecision(1) << "\b\b\b\b\b\b\b\b\b\b";
        cout << 100.*static_cast<double>(iz)/static_cast<double>(zpoints) << "%";

        double z = zmin + static_cast<double>(iz)*dz;

        for ( int iy=0; iy<ypoints; iy++ )
        {
            double y = ymin + static_cast<double>(iy)*dy;
            double square_of_diff_y_y0 = (y-y0)*(y-y0);
            for ( int ix=0; ix<xpoints; ix++ )
            {
                double x = xmin + static_cast<double>(ix)*dx;
                double square_of_diff_x_x0 = (x-x0)*(x-x0);
                complex<double> psi = complex<double> (0., 0.);
                double varrho = sqrt( square_of_diff_x_x0 + square_of_diff_y_y0 );
                for ( int in=0; in<inmax; in++ )
                {
                    psi += A[in]*bessel::cyl_j0(h[in]*varrho)*exp( b[in]*complex<double>(0., z) );
                }
                field[ix + xpoints*iy + xpoints*ypoints*iz] = psi;
            }
        }
    }


    // Clearing progress.
    cout << "\b\b\b\b\b\b\b\b\b\b";


    // Writing field to WL file.
    puts("Exporting scalar field...");
    wldat_export( "field.wl", field, xpoints, ypoints, zpoints, 16, true );
    delete[] field;


    // Printing time counted.
    cout << "Done! Elapsed time: " << ( time(nullptr) - start_time ) << " seconds." << endl;


    system("pause");
}