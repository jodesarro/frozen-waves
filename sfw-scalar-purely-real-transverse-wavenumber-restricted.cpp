// Cited references available at README.md.

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
    puts("| SFW-SCALAR-PURELY-REAL-TRANSVERSE-WAVENUMBER-RESTRICTED | 02 JUL 2022 | @JODESARRO |");
    puts("A C++ routine to evaluate the scalar field of surface frozen waves ensuring purely real transverse wavenumber and restricted to the case where Qm=Q, Lm=L and Nm=N.");
    puts("https://github.com/jodesarro/frozen-waves\n");


    // Specification of general parameters of frozen waves.
    puts("Loading general parameters...");
    complex<double> nref = complex<double> ( 1.4 , 2.56e-6 );
    double l0 = 632.0e-9;
    double k0 = 2.0*M_PI/l0;
    complex<double> k = k0*nref;
    static const int N = 15; // It is assumed Nm=N.
    static const int M = 400;


    // Calculation of maximum values for indexes n and m of frozen waves.
    //  n->[in] : n=in-N, 0<=in<inmax
    //  m->[im] : m=im+1, 0<=im<immax
    static const int inmax = N*2+1;
    static const int immax = M;


    // Specification of L parameter.
    //  It is assumed Lm=L.
    double L = 0.05;


    // Specification of Q parameter.
    //  It is assumed Qm=Q.
    double Q = 0.999*real(nref)*k0;
    cout << "A chosen value for Q parameter: " << Q << " 1/m" << endl;


    // Calculation of maximum value possible for N, Equation (7) of [2].
    int Nmax = floor( L * min( k0*sqrt(0.5*( real(nref)*real(nref) - imag(nref)*imag(nref)
            + sqrt( pow(real(nref)*real(nref) - imag(nref)*imag(nref),2.0)
            + pow(2.0*real(nref)*imag(nref),2.0) ))) - Q, Q ) / (2.0*M_PI) );
    cout << "Calculated maximum value possible for N: " << Nmax << endl;


    // Calculation of wavenumbers, Equations (4) and (5) of [2].
    //  The hnm=hn and βnm=βn are consequences of Nm=N, Qm=Q and Lm=L.
    puts("Calculating wavenumbers...");
    complex<double> h [inmax];
    complex<double> b [inmax];
    for ( int in=0; in<inmax; in++ )
    {
        double n = static_cast<double>( in - N );
        double bnr = Q + 2.0*M_PI*n/L;
        double bni = k0*k0*real(nref)*imag(nref)/bnr;
        b[in] = complex<double> ( bnr, bni );
        h[in] = sqrt( (real(nref)*real(nref)-imag(nref)*imag(nref))*k0*k0 - bnr*bnr + bni*bni );
    }


    // Calculation of Δ parameter.
    //  The Δm=Δ is a consequence of Nm=N, Qm=Q and Lm=L.
    double delta = abs( imag(b[inmax-1]) - imag(b[0]) )/imag(b[N]);
    cout << "Calculated delta parameter: " << delta << endl;


    // Calculation of penetration depth of Bessel beams using the idea of Chapter II of [2].
    //  The δnm=δn is a consequence of βnm=βn.
    double del_min = 0.5/imag(b[0]);
    cout << "Calculated penetration depth for n=-N: " << del_min << " m" << endl;
    double del_max = 0.5/imag(b[inmax-1]);
    cout << "Calculated penetration depth for n=+N: " << del_max << " m" << endl;


    // Calculation of axicon angles for generating Bessel beams outside the considered media.
    // Although the complex axicon angles inside the considered media were not discussed in Chapter II of [2],
    //  it is possible to calculate real axicon angles outside that media. Below is the calculation of the
    //  axicon angles in a dielectric media of unitary refractive index (air or vacuum)
    //  for z≤0- where the beams are supposed to be generated.
    //  The θnm=θn is a consequence of βnm=βn.
    double theta_min = asin( real(h[inmax-1])/k0 );
    cout << "Calculated axicon angle in vacuum (z<=0-) for n=+N: ";
    cout << theta_min << " rad = " << (360./(2.0*M_PI))*theta_min << " degree" << endl;
    double theta_max = asin( real(h[0])/k0 );
    cout << "Calculated axicon angle in vacuum (z<=0-) for n=-N: ";
    cout << theta_max << " rad = " << (360./(2.0*M_PI))*theta_max << " degree" << endl;


    // Calculation of aperture radius for experimental generation of Bessel beams inside a media of
    //  unitary refractive index (air or vacuum) in z≤0- where the beams is supposed to be generated.
    //  The Rnm=Rn is a consequence of θnm=θn.
    double calculated_aperture_bb_min = L*tan(theta_min);
    cout << "Calculated aperture radius in vacuum (z<=0-) for n=+N: ";
    cout << calculated_aperture_bb_min << " m" << endl;
    double calculated_aperture_bb_max = L*tan(theta_max);
    cout << "Calculated aperture radius in vacuum (z<=0-) for n=-N: ";
    cout << calculated_aperture_bb_max << " m" << endl;


    // Calculation of maximum value possible for N for experimental generation of Bessel beams inside a media of
    //  unitary refractive index (air or vacuum) in z≤0- where the beams is supposed to be generated.
    int Nmax2 = floor( L * min( Q - k0*sqrt(0.5*( real(nref)*real(nref) - imag(nref)*imag(nref) - 1.0
        + sqrt( pow(real(nref)*real(nref) - imag(nref)*imag(nref) - 1.0,2.0)
        + pow(2.0*real(nref)*imag(nref),2.0) ))), Q ) / (2.0*M_PI) );
    cout << "Calculated maximum value possible for N to generate the beam in vacuum (z<=0-): " << Nmax2 << endl;


    // Calculation of spot radius of a ideal linear frozen wave.
    //  That is, spot radius of the ideal Bessel beam with n=0 [3,4].
    double calculated_spot_radius_lfw = 2.405/real(h[N]);
    cout << "Calculated LFW spot radius: ";
    cout << calculated_spot_radius_lfw << " m" << endl;


    // Calculation of spot radius of a ideal linear frozen wave by asymptotic expansion.
    //  That is, spot radius of the ideal Bessel beam with n=0 by asymptotic expansion.
    double calculated_asymptotic_spot_radius_lfw = 0.25*3.0*M_PI/real(h[N]);
    cout << "Calculated LFW spot radius by asymptotic expansion: ";
    cout << calculated_asymptotic_spot_radius_lfw << " m" << endl;


    // Specification of resistant (attenuation or compensation) parameter as defined in [4].
    //  The βI0m=βI0 is a consequence of Qm=Q.
    //  Set bi0=imag(b[N]) for attenuation-resistant method
    //  or bi0=0. for non-attenuation-resistant method.
    double bi0 = imag(b[N]);


    // Specification of positions x0m and y0m relative to the origin.
    //  There are two possibilities for a surface frozen wave [3]:
    //   (a) a cartesian surface (plane) inclined by a constant and m-independent angle φ0=atan(y0m/x0m);
    //   (b) a cylindrical surface of constant and m-independent radius ρ0=sqrt(x0m²+y0m²);
    //  Below it was assumed (a) with φ0=0, y0m=0 and x0m=(m-1)*dx0, where dx0 is the separation Δx0 in x, as assumed in [4].
    double dx0 = 4.0*calculated_spot_radius_lfw;
    double x0 [immax];
    double x0x0 [immax];
    for ( int im=0; im<immax; im++ )
    {
        double m = static_cast<double>(im + 1);
        x0[im] = (m-1.0)*dx0;
    }
    cout << "Separation between consecutives LFWs: " << dx0 << " m" << endl;


    // Specification of data calculation parameters.
    double xmin = 0.;
    double xmax = x0[immax-1];
    static const int xpoints = 150;
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


    // Reading MTX intensity profile file to function F[ii+iimax*ik].
    //  x->[ii] : 0<=ii<iimax
    //  z->[ik] : 0<=ik<ikmax
    puts("Loading file containing the intensity profile...");
    int iimax;
    int ikmax;

    mtxdat_getsize( "intensity-profile.mtx", iimax, ikmax );
    double * F = new double [iimax*ikmax];
    mtxdat_import( "intensity-profile.mtx", F, iimax, ikmax );


    // Calculation of A coefficients by means of an approximation of the integral
    //  by trapezoidal method for equally spaced z using Equation (4) of [4].
    //  The Anm = An is a consequence of Lm=L, Qm=Q, Nm=N.
    puts("Calculating the A coefficients...");
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
                aux += F[ii + iimax*ik]*exp( complex<double>(0., -2.0*M_PI*n*z/L) )*exp( bi0*z );
            }
            A[in + inmax*im] = aux / static_cast<double>(ikmax-1);
        }
    }
    delete[] F;


    // Calculation of scalar field by Equation (2) of [4].
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
            double square_of_diff_y_y0 = (y-0.)*(y-0.);
            for ( int ix=0; ix<xpoints; ix++ )
            {
                double x = xmin + static_cast<double>(ix)*dx;
                complex<double> psi = complex<double> (0., 0.);
                for ( int im=0; im<immax; im++ )
                {
                    double square_of_diff_x_x0 = (x-x0[im])*(x-x0[im]);
                    double varrho = sqrt( square_of_diff_x_x0 + square_of_diff_y_y0 );
                    for ( int in=0; in<inmax; in++ )
                    {
                        psi += A[in + inmax*im]*bessel::cyl_j0(h[in]*varrho)*exp( b[in]*complex<double>(0., z) );
                    }
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