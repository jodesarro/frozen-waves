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
void field_for_x( int ix, int xpoints, int ypoints, double zmin, double dz, double zpoints, int ismax, double * d,
    complex<double> * B, complex<double> * b, complex<double> * bj0, complex<double> * tau, complex<double> * gamma,
    int inmax, int immax, complex<double> * field )
{
    for ( int iy=0; iy<ypoints; iy++ )
    {
        for ( int iz=0; iz<zpoints; iz++ )
        {
            double z = zmin + static_cast<double>(iz)*dz;

            // Identifying the medium s (the number is).
            int is;
            for ( is=0; is<ismax-1; is++ )
            {
                if ( z < d[is] )
                {
                    break;
                }
            }

            // Calculationg field
            complex<double> psi = complex<double> (0., 0.);
            for ( int im=0; im<immax; im++ )
            {
                for ( int in=0; in<inmax; in++ )
                {
                    psi += B[in + inmax*im]*bj0[in+inmax*im+inmax*immax*ix+inmax*immax*xpoints*iy]*
                            ( tau[in+inmax*is]*exp(I*b[in+inmax*is]*z) + gamma[in+inmax*is]*exp(-I*b[in+inmax*is]*z) );
                }
            }
            field[ix + xpoints*iy + xpoints*ypoints*iz] = psi;
        }
    }
}


// Function to calculate Bessel functions given x.
void bj0_for_x( int ix, double x, int xpoints, double * x0, double ymin, double dy, int ypoints, double * y0,
    double * h, int inmax, int immax, complex<double> * bj0 )
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
    cout << "| LSLS-SCALAR-STRATIFIED-RESTRICTED | 15 NOV 2023 | @JODESARRO |" << endl;
    cout << "A C++ routine to evaluate the scalar field of longitudinally structured light sheets in stratified structure. The beam is structured in a lossy or lossless medium after passing through several lossless media and is restricted to the case where Qm=Q, Lm=L and Nm=N." << endl;
    cout << "https://github.com/jodesarro/frozen-waves\n" << endl;
    

    // Specification of general parameters of frozen waves.
    cout << "Loading general parameters..." << endl;
    double l0 = 632.0e-9;
    double k0 = 2.0*M_PI/l0;
    static const int N = 60; // It is assumed Nm=N.
    static const int M = 401;
    cout << "Chosen value for the wavelength (in vacuum): " << l0 << " m" << endl;
    cout << "Calculated value for the angular frequency (in vacuum): " << k0*299792458 << " rad/s" << endl;
    cout << "Chosen value for N parameter: " << N << endl;
    cout << "Chosen value for M parameter: " << M << endl;   


    // Calculation of maximum values for indexes n and m of frozen waves.
    //  n->[in] : n=in-N, 0<=in<inmax
    //  m->[im] : m=im+1, 0<=im<immax
    static const int inmax = N*2+1;
    static const int immax = M;


    // Parameters related to the stratified structure.

    // Number S>1 of media. The media are represented by index is.
    //  s->[is] : is=s-1, 0<=is<ismax
    static const int ismax = 1+2+1;

    double d [ismax-1] = { 0., 0.015, 0.02 };
    complex<double> nref [ismax] = { 1., 2.5, 2., 1.5 };

    // Solution of sugar and water
    /*double Dd = 0.001;
    double * d = new double [ismax-1];
    d[0] = 0.;
    complex<double> * nref = new complex<double> [ismax];
    nref[0] = 1.;
    for (int is=1; is<ismax-1; is++)
    {
        d[is] = d[is-1] + Dd;
        double z = ( double(is) - 0.5 )*Dd;
        nref[is] = 1.335 + 0.03657/( 1.+0.0006778*exp(z*134.6) );
    }
    nref[ismax-1]=1.;*/
   
    cout << "Number S of layers: " << ismax << endl;    
    cout << "Index of refraction for s=1: " << nref[0] << endl;
    cout << "Index of refraction for s=S: " << nref[ismax-1] << endl;
    cout << "Longitudinal position of the last interface (s=S-1): " << d[ismax-2] << " m" << endl;


    // Specification of L parameter.
    //  It is assumed Lm=L.
    double L = 0.04 + d[ismax-2];
    cout << "Chosen value for L parameter: " << L << " m" << endl;   


    // Longitudinal space of nulls to add in the end of F.
    //  The DLm = DL is consequence of Lm=L.
    double DL = d[ismax-2];


    // Specification of Q parameter.
    //  It is assumed Qm=Q.

    //  Use the next 10 lines only to calculate Q with a chosen spot radius.
    //double chosen_spot_radius_fw = 9.0e-6;
    //double a;
    //if ( imag(nref[ismax-1]) != 0. )
    //{
    //    //a = ? //not calculated yet.
    //}
    //else
    //{
    //    a = sqrt( 1.0 - pow(l0*2.4048/(2.0*M_PI*real(nref[ismax-1])*chosen_spot_radius_fw), 2.0 );
    //}
    
    //  Use the next line only to specify a fixed value for Q.
    double a = 0.9986;
    
    double Q = a*k0*real(nref[ismax-1]);
    
    cout << "Chosen value for the a parameter: " << a << endl;   
    cout << "Calculated value for Q parameter: " << Q << " 1/m" << endl;


    // Specification of positions x'm and y'm relative to the origin.
    double dx0 = 2.*4.57325e-06;
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


    // Specification of reflection compensation method.
    //  Set generate_with_B=true to use the coefficient Bnm to generate the beam inside the first lossless medium (Bnm=Anm/Tnm)
    //  or generate_with_B=false to use the traditional coefficient Anm to generate the beam inside the first lossless medium.
    bool generate_with_B = true;


    // Specification of data calculation parameters.
    double xmin = 0.;
    double xmax = x0[immax-1];
    static const int xpoints = 801;
    double ymin = 0.;
    double ymax = 0.;
    static const int ypoints = 1;
    double zmin = 0.13;
    double zmax = L;
    static const int zpoints = 201;
    cout << "x E (" << xmin << " m, " << xmax << " m, " << xpoints << " pt)" << endl;
    cout << "y E (" << ymin << " m, " << ymax << " m, " << ypoints << " pt)" << endl;
    cout << "z E (" << zmin << " m, " << zmax << " m, " << zpoints << " pt)" << endl;


    // Specification of async or lazy calculation of xy planes.
    bool is_async = true;
    auto launch_policy = (is_async) ? (launch::async) : (launch::deferred) ;


    // Usually there is no need to change values in all of the following.


    // Calculation of wavenumbers.
    //  The hnm=hn and βnms=βns are consequences of Nm=N, Qm=Q and Lm=L.
    cout << "Calculating wavenumbers..." << endl;
    complex<double> * b = new complex<double> [inmax*ismax];
    double h [inmax];
    for ( int in=0; in<inmax; in++ )
    {
        double n = static_cast<double>( in - N );
        double brnS = Q + 2.0*M_PI*n/L;
        double binS = k0*k0*real(nref[ismax-1])*imag(nref[ismax-1])/brnS;
        b[in + inmax*(ismax-1)] = complex<double> ( brnS, binS );
        h[in] = sqrt( real( k0*nref[ismax-1]*k0*nref[ismax-1] - b[in + inmax*(ismax-1)]*b[in + inmax*(ismax-1)] ) );
        for (int is=0; is<ismax-1; is++)
        {
            b[in + inmax*is] = sqrt( k0*nref[is]*k0*nref[is] - h[in]*h[in] );
        }
    }


    // Resistant attenuation or compensation parameter.
    //  The βI0mS=βI0S is a consequence of Qm=Q.
    double bi0S = imag(b[N+inmax*(ismax-1)]);


    // Calculation of minimum value possible for the refractive indices smaller than nref[ismax-1].
    //  The nrefminm=nrefmin is consequence of am=a.
    double nrefmin = sqrt( real(nref[ismax-1])*real(nref[ismax-1])*(1.0-a*a) + imag(nref[ismax-1])*imag(nref[ismax-1])*(1.0/(a*a)-1.0) );
    cout << "Calculated minimum value possible for the refractive indices smaller than the last one: " << nrefmin << endl;


    // Calculation of maximum value possible for N.
    //  The Nmaxm=Nmax is consequence of Nm=N, Qm=Q and Lm=L.
    double * kprime = new double [ismax];
    for ( int is=0; is<ismax; is++ )
    {
        if ( imag(nref[ismax-1]) == 0. && real(nref[ismax-1]) <= real(nref[is]) )
        {
            kprime[is] = 0.;
        }
        else
        {
            double temp = real(nref[ismax-1])*real(nref[ismax-1]) - imag(nref[ismax-1])*imag(nref[ismax-1]) - real(nref[is])*real(nref[is]);
            kprime[is] = k0*sqrt(0.5*( temp + sqrt( temp*temp + pow(2.0*real(nref[ismax-1])*imag(nref[ismax-1]), 2.0) )));
        }
        
    }
    int Nmax = floor( L * min( Q - *max_element(kprime, kprime + ismax), k0*real(nref[ismax-1]) - Q )/(2.0*M_PI) );
    cout << "Calculated maximum value possible for N: " << Nmax << endl;


    // Calculation of reflection and transmission coefficients.
    // The Tnm=Tn, Rnm=Rn, τnms=τns and Γnms=Γns are consequences of Nm=N, Qm=Q and Lm=L.
    cout << "Calculation of reflection and transmission coefficients..." << endl;
    complex<double> * tau = new complex<double> [inmax*ismax];
    complex<double> * gamma = new complex<double> [inmax*ismax];
    complex<double> T [inmax];
    complex<double> R [inmax];
    complex<double> * M11 = new complex<double> [inmax*ismax];
    complex<double> * M12 = new complex<double> [inmax*ismax];
    complex<double> * M21 = new complex<double> [inmax*ismax];
    complex<double> * M22 = new complex<double> [inmax*ismax];
    for ( int in=0; in<inmax; in++ )
    {
        // Matrix M for s=1 (for is=0) is not defined.

        // Matrix M for s=2 (for is=1).
        M11[in + inmax*1] = 1.;
        M12[in + inmax*1] = 0.;
        M21[in + inmax*1] = 0.;
        M22[in + inmax*1] = 1.;

        // Matrix M for s=3 to s=S (for 2<=is<ismax).
        for (int is=2; is<ismax; is++)
        {
            // Calculation of matrix N for s-1, i.e., for s=2 to S-1 (for 1<=is<ismax-1).
            double deltadsm1 = d[is-1]-d[is-2];
            complex<double> N11sm1 = cos( b[in + inmax*(is-1)] * deltadsm1 );
            complex<double> N12sm1 = sin( b[in + inmax*(is-1)] * deltadsm1 )/b[in + inmax*(is-1)];
            complex<double> N21sm1 = -b[in + inmax*(is-1)]*sin( b[in + inmax*(is-1)] * deltadsm1 );
            complex<double> N22sm1 = cos( b[in + inmax*(is-1)] * deltadsm1 );

            // Calculation of matrix M for s, i.e., for s=3 to S (for 2<=is<ismax).
            M11[in + inmax*is] = N11sm1*M11[in + inmax*(is-1)] + N12sm1*M21[in + inmax*(is-1)];
            M12[in + inmax*is] = N11sm1*M12[in + inmax*(is-1)] + N12sm1*M22[in + inmax*(is-1)];
            M21[in + inmax*is] = N21sm1*M11[in + inmax*(is-1)] + N22sm1*M21[in + inmax*(is-1)];
            M22[in + inmax*is] = N21sm1*M12[in + inmax*(is-1)] + N22sm1*M22[in + inmax*(is-1)];
        }

        complex<double> M11S = M11[in + inmax*(ismax-1)];
        complex<double> M12S = M12[in + inmax*(ismax-1)];
        complex<double> M21S = M21[in + inmax*(ismax-1)];
        complex<double> M22S = M22[in + inmax*(ismax-1)];
        complex<double> bn1 = b[in + inmax*0];
        complex<double> bnS = b[in + inmax*(ismax-1)];

        // Coefficient of transmission to the last medium.
        T[in] = ( I*2.*bn1 )/
                ( -M21S + bn1*bnS*M12S + I*( bn1*M22S + bnS*M11S ) );

        // Coefficient of reflection to the first medium.
        R[in] = ( +M21S + bn1*bnS*M12S + I*( bn1*M22S - bnS*M11S ) )/
                ( -M21S + bn1*bnS*M12S + I*( bn1*M22S + bnS*M11S ) );

        // Transmission and reflection coefficients in the first medium s=1 (for is=0).
        tau[in + inmax*0] = 1.;
        gamma[in + inmax*0] = R[in];

        // Transmission and reflection coefficients in a medium s for s=2 to S-1 (for 1<=is<ismax-1).
        for (int is=1; is<ismax-1; is++)
        {
            tau[in + inmax*is] = ( (1.0+R[in])*(b[in + inmax*is]*M11[in + inmax*is]-I*M21[in + inmax*is])
                    + (1.0-R[in])*(bn1*M22[in + inmax*is] + I*bn1*b[in + inmax*is]*M12[in + inmax*is]) )
                    /(2.0*b[in + inmax*is]*exp(I*b[in + inmax*is]*d[is-1] ));

            gamma[in + inmax*is] = ( (1.0+R[in])*(b[in + inmax*is]*M11[in + inmax*is]+I*M21[in + inmax*is])
                    - (1.0-R[in])*(bn1*M22[in + inmax*is] - I*bn1*b[in + inmax*is]*M12[in + inmax*is]) )
                    /(2.0*b[in + inmax*is]*exp(-I*b[in + inmax*is]*d[is-1] ));
        }

        // Transmission and reflection coefficients in the last medium s=S (for is=ismax-1).
        tau[in + inmax*(ismax-1)] = T[in] * exp(-I*bnS*d[ismax-2] );
        gamma[in + inmax*(ismax-1)] = 0.;

    }


    // Calculation of axicon angles for generating Bessel beams inside the first lossless medium
    //  in z≤0 where the beams are supposed to be generated.
    //  The θnm=θn is a consequence of hnm=hn.
    double theta_pN = asin( real(h[inmax-1])/(k0*real(nref[0])) );
    cout << "Calculated axicon angle inside first lossless medium (z<=0) for n=+N: ";
    cout << theta_pN << " rad = " << (360./(2.0*M_PI))*theta_pN << " degree" << endl;
    double theta_mN = asin( real(h[0])/(k0*real(nref[0])) );
    cout << "Calculated axicon angle inside first lossless medium (z<=0) for n=-N: ";
    cout << theta_mN << " rad = " << (360./(2.0*M_PI))*theta_mN << " degree" << endl;


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


    // Calculation of parameters related to complex refractive index of the last medium.
    if ( imag(nref[ismax-1]) != 0. )
    {
        // Calculation of Δ parameter.
        //  The Δm=Δ is a consequence of Nm=N, Qm=Q and Lm=L.
        double delta = abs( imag(b[inmax-1 + inmax*(ismax-1)]) - imag(b[0 + inmax*(ismax-1)]) )/imag(b[N + inmax*(ismax-1)]);
        cout << "Calculated delta parameter: " << delta << endl;


        // Calculation of penetration depth of Bessel beams inside the last lossy medium (z≥d[ismax-2]).
        //  The δnm=δn is a consequence of βnm=βn.
        double del_pN = 0.5/imag(b[inmax-1 + inmax*(ismax-1)]);
        cout << "Calculated penetration depth inside last lossy medium for n=+N: " << del_pN << " m" << endl;
        double del_mN = 0.5/imag(b[0 + inmax*(ismax-1)]);
        cout << "Calculated penetration depth inside last lossy medium for n=-N: " << del_mN << " m" << endl;
    }


    // Calculation of total transmission coefficient.
    //  The Tnm=Tn is a consequence of Nm=N, Qm=Q and Lm=L.
    complex<double> calculated_transmission_coef_pN = T[inmax-1];
    cout << "Calculated total transmission coefficient for n=+N: " << calculated_transmission_coef_pN << " = ";
    cout << abs(calculated_transmission_coef_pN) << " exp(" << arg(calculated_transmission_coef_pN) << " i)" << endl;
    complex<double> calculated_transmission_coef_0 = T[N];
    cout << "Calculated total transmission coefficient for n=0: " << calculated_transmission_coef_0 << " = ";
    cout << abs(calculated_transmission_coef_0) << " exp(" << arg(calculated_transmission_coef_0) << " i)" << endl;
    complex<double> calculated_transmission_coef_mN = T[0];
    cout << "Calculated total transmission coefficient for n=-N: " << calculated_transmission_coef_mN << " = ";
    cout << abs(calculated_transmission_coef_mN) << " exp(" << arg(calculated_transmission_coef_mN) << " i)" << endl;


    // Reading MTX function-f file to a function F[ii+iimax*ik].
    //  Importing F and compressing it to 0≤z≤L-DL and adding zeros in L-DL<z≤L
    //  x->[ii] : 0<=ii<iimax
    //  z->[ik] : 0<=ik<ikmax
    cout << "Loading file containing the function F..." << endl;
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


    // Calculation of A and B coefficients by means of an approximation of the integral
    //  by trapezoidal method for equally spaced z.
    //  The Anm = An and Bnm=Bn is a consequence of Lm=L, Qm=Q, Nm=N.
    cout << "Calculating the A and B coefficients..." << endl;
    complex<double> A [inmax*immax];
    complex<double> B [inmax*immax];
    for ( int im=0; im<immax; im++ )
    {
        for ( int in=0; in<inmax; in++ )
        {
            double n = static_cast<double>(in - N);
            int ii = floor( static_cast<double>(im)*static_cast<double>(iimax-1)/static_cast<double>(immax-1) );
            complex<double> aux = 0.5*( F[ii + iimax*(ikmax-1)]*exp( bi0S*L ) + F[ii + 0] );
            for ( int ik=1; ik<ikmax-1; ik++ )
            {
                double z = L*static_cast<double>(ik)/static_cast<double>(ikmax-1);
                aux += F[ii + iimax*ik]*exp( -I*2.0*M_PI*n*z/L )*exp( bi0S*z );
            }
    
            A[in + inmax*im] = aux / static_cast<double>(ikmax-1);

            if ( generate_with_B )
            {
                B[in + inmax*im] = A[in + inmax*im]/T[in];
            }
            else
            {
                B[in + inmax*im] = A[in + inmax*im];
            }

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
            ix, xpoints, ypoints, zmin, dz, zpoints, ismax, d, B, b, bj0, tau, gamma, inmax, immax, field ) );
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