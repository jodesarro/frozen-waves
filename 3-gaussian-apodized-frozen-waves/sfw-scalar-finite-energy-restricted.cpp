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
void field_for_x( int ix, int xmin, double dx, int xpoints, int ymin, double dy, int ypoints, double zmin, double dz, double zpoints,
    double * x0, double * y0, complex<double> * b, complex<double> * h,  complex<double> k, double w0w0,  complex<double> * A,
    int inmax, int immax, complex<double> * field )
{
    double x = xmin + static_cast<double>(ix)*dx;
    for ( int iy=0; iy<ypoints; iy++ )
    {
        double y = ymin + static_cast<double>(iy)*dy;
        for ( int iz=0; iz<zpoints; iz++ )
        {
            double z = zmin + static_cast<double>(iz)*dz;
            complex<double> mu = 1.0 + complex<double>(0., 2.)*z/( w0w0*k );
            complex<double> w0w0mu = w0w0*mu;
            complex<double> psi = complex<double> (0., 0.);
            for ( int im=0; im<immax; im++ )
            {
                double square_of_diff_x_x0 = (x-x0[im])*(x-x0[im]);
                double square_of_diff_y_y0 = (y-y0[im])*(y-y0[im]);
                double varrho = sqrt( square_of_diff_x_x0 + square_of_diff_y_y0 );
                double varrhovarrho = varrho*varrho;
                for ( int in=0; in<inmax; in++ )
                {
                    psi += A[in + inmax*im]*bessel::cyl_j0(h[in]*varrho/mu)
                            *exp( + I*b[in]*z/mu - I*z*k/mu + I*z*k - varrhovarrho/w0w0mu )/mu;
                }
            }
            field[ix + xpoints*iy + xpoints*ypoints*iz] = psi;
        }
    }
}


int main()
{
    // Starting time counter.
    long int start_time = time( nullptr );


    // Printing information.
    puts("| SFW-SCALAR-FINITE-ENERGY-RESTRICTED | 21 APR 2023 | @JODESARRO |");
    puts("A C++ routine to evaluate the scalar field of surface frozen waves of finite energy and restricted to the case where Qm=Q, Lm=L, Nm=N and w0m=w0.");
    puts("https://github.com/jodesarro/frozen-waves\n");


    // Specification of general parameters of frozen waves.
    puts("Loading general parameters...");
    complex<double> nref = complex<double> ( 1.0 , 0. );
    //double l0 = 451.0e-9; // B
    //double l0 = 532.0e-9; // G
    //double l0 = 635.0e-9; // R
    double l0 = 451.0e-9;
    double k0 = 2.0*M_PI/l0;
    static const int N = 19; // It is assumed Nm=N.
    static const int M = 25;


    // Calculation of maximum values for indexes n and m of frozen waves.
    //  n->[in] : n=in-N, 0<=in<inmax
    //  m->[im] : m=im+1, 0<=im<immax
    static const int inmax = N*2+1;
    static const int immax = M;


    // Specification of L parameter.
    //  It is assumed Lm=L.
    double L = 0.05;


    // Longitudinal space of nulls to add in the end of F.
    //  The DLm = DL is consequence of Lm=L.
    double DL = 0.;


    // Specification of Q parameter.
    //  It is assumed Qm=Q.

    //  Use the next 2 lines only to calculate Q with a chosen spot radius.
    //double chosen_spot_radius_fw = 9.0e-6;
    //double Q = ( 1.0 - 0.5*(2.405/(chosen_spot_radius_fw*abs(k)))*(2.405/(chosen_spot_radius_fw*abs(k))) ) * real(k);
    
    //  Use the next line only to specify a fixed value for Q.
    double Q = 0.9986*real(nref)*k0;
    
    cout << "Calculated value for Q parameter: " << Q << " 1/m" << endl;


    // Specification of gaussian waist radius (z=0) to apodize each frozen wave.
    //  The Gaussian waist radius is the inverse of the q parameter: w0m=1/qm.
    //  It is assumed w0m=w0.
    double w0 = 1000. * M_SQRT2 * (2.405/(M_SQRT2*sqrt(1.0-Q/(real(nref)*k0))*abs(k0*nref))); // n times spot radius
    //double w0 = L * (M_SQRT2*sqrt(1.0-Q/(real(nref)*k0))*abs(k0*nref)) / ( 2.0*real(nref)*k0 );
    cout << "Gaussian waist radius (z=0): " << w0 << " m" << endl;
    double w0w0 = w0*w0;
   

    // Specification of positions x'm and y'm relative to the origin.
    double dx0 = L/214.;
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


    // Specification of whether to compensate for apodization with G.
    //  Set true to compensate or false to not compensate.
    bool is_apodization_compensation = true;


    // Specification of data calculation parameters.
    double xmin = 0.;
    double xmax = double(M-1)*dx0;
    static const int xpoints = 401;
    double ymin = 0.;
    double ymax = 0.;
    static const int ypoints = 1;
    double zmin = 0.;
    double zmax = 2.*L;
    static const int zpoints = 429;
    cout << "x E (" << xmin << " m, " << xmax << " m, " << xpoints << " pt)" << endl;
    cout << "y E (" << ymin << " m, " << ymax << " m, " << ypoints << " pt)" << endl;
    cout << "z E (" << zmin << " m, " << zmax << " m, " << zpoints << " pt)" << endl;


    // Specification of async or lazy calculation of xy planes.
    bool is_async = true;
    auto launch_policy = (is_async) ? (launch::async) : (launch::deferred) ;


    // Usually there is no need to change values in all of the following.


    // Calculation of wavenumbers.
    //  The hnm=hn and βnm=βn are consequences of Nm=N, Qm=Q and Lm=L.
    puts("Calculating wavenumbers...");
    complex<double> k = k0*nref;
    complex<double> h [inmax];
    complex<double> b [inmax];
    for ( int in=0; in<inmax; in++ )
    {
        double n = static_cast<double>( in - N );
        double bnr = Q + (2.0*M_PI)*n/L;
        double bni = imag(k)*( 2.0 - bnr/real(k) );
        b[in] = complex<double> ( bnr, bni );
        h[in] = M_SQRT2*sqrt( 1.0 - bnr/real(k) )*abs(k);
    }

    
    // Resistant attenuation or compensation parameter.
    //  The βI0m=βI0 is a consequence of Qm=Q.
    double bi0 = imag(b[N]);


    // Calculation of maximum value possible for N.
    //  Is is assumed Nm=N.
    int Nmax = floor( L * min( k0*real(nref) - Q, Q ) / (2.0*M_PI) );
    cout << "Calculated maximum value possible for N: " << Nmax << endl;


    // Calculation of Gaussian spot radius (transverse intensity width) apodization at z=0.
    //  The ΔρGm=ΔρG is a consequence of w0m=w0.
    double calculated_spot_radius_gb = w0 / M_SQRT2;
    cout << "Calculated Gaussian spot radius (transverse intensity width) apodization at z=0: " << calculated_spot_radius_gb << " m" << endl;


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
    cout << "Calculated Bessem beam aperture radius for n=+N: ";
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
    }


    // Reading MTX function-f file to a function F[ii+iimax*ik].
    //  Importing F and compressing it to 0≤z≤L-DL and adding zeros in L-DL<z≤L
    //  x->[ii] : 0<=ii<iimax
    //  z->[ik] : 0<=ik<ikmax
    puts("Loading file containing the function F...");
    int iimax;
    int ikmax;

    mtxdat_getsize( "function-f.mtx", iimax, ikmax );
    double ik_per_meter = (ikmax - 1.)/(L-DL - 0.);
    int ikmax_old = ikmax;
    ikmax += floor( ik_per_meter*(DL-0.) );
    double * F = new double [iimax*ikmax];
    mtxdat_import( "function-f.mtx", F, iimax, ikmax_old );
    for ( int i = iimax*ikmax_old; i < iimax*ikmax; i++ )
    {
        F[i] = double(0);
    }


    // Calculation of A coefficients by means of an approximation of the integral
    //  by trapezoidal method for equally spaced z.
    //  The Anm = An is a consequence of Lm=L, Qm=Q, Nm=N.
    puts("Calculating the A coefficients...");
    complex<double> A [inmax*immax];
    if (is_apodization_compensation)
    {
        // Using finite energy method.    
        complex<double> G0 = 1.0;
        complex<double> GL = exp( - imag(k)*L - h[N]*h[N]*L*L/( w0w0*real(k)*real(k) ) );
        for ( int im=0; im<immax; im++ )
        {
            for ( int in=0; in<inmax; in++ )
            {
                double n = static_cast<double>(in - N);
                int ii = floor( static_cast<double>(im)*static_cast<double>(iimax-1)/static_cast<double>(immax-1) );
                complex<double> aux = 0.5*( F[ii + iimax*(ikmax-1)]/GL + F[ii + 0]/G0 );
                for ( int ik=1; ik<ikmax-1; ik++ )
                {
                    double z = L*static_cast<double>(ik)/static_cast<double>(ikmax-1);
                    complex<double> Gz = exp( - imag(k)*z - h[N]*h[N]*z*z/( w0w0*real(k)*real(k) ) );
                    aux += ( F[ii + iimax*ik]/Gz )*exp( -I*2.0*M_PI*n*z/L );
                }
                A[in + inmax*im] = aux / static_cast<double>(ikmax-1);
            }
        }
    }
    else
    {
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
    }
    delete[] F;


    // Calculation of spatial increments.
    double dx, dy, dz;
    dx = (xpoints==1) ? (0.) : (xmax-xmin)/(xpoints-1) ;
    dy = (ypoints==1) ? (0.) : (ymax-ymin)/(ypoints-1) ;
    dz = (zpoints==1) ? (0.) : (zmax-zmin)/(zpoints-1) ;


    // Calculation of scalar field.
    puts("Calculating scalar field...");
    complex<double> * field = new complex<double> [xpoints*ypoints*zpoints];
    std::vector<std::future<void>> future_field_for_x;
    for ( int ix=0; ix<xpoints; ix++ )
    {
        double x = xmin + static_cast<double>(ix)*dx;
        future_field_for_x.push_back( async( launch_policy, field_for_x,
            ix, xmin, dx, xpoints, ymin, dy, ypoints, zmin, dz, zpoints, x0, y0, b, h, k, w0w0, A, inmax, immax, field ) );
    }
    for(auto &e : future_field_for_x)
    {
        e.get();
    }


    // Writing field to WL file.
    puts("Exporting scalar field...");
    wldat_export( "field.wl", field, xpoints, ypoints, zpoints, 16, true );
    delete[] field;


    // Printing time counted.
    cout << "Done! Elapsed time: " << ( time(nullptr) - start_time ) << " seconds." << endl;


    system("pause");
}