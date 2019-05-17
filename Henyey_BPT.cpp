#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>


using namespace std;

/*
 This is a first try at using the Henyey method.  So I'll start out solving some equations which are a bit simpler than the stuff that I'll use in the oscillation equations case.  Mainly avoiding having to read any input data at this stage.
 
 Testing the different recurrence relations, this is:
 
 ************
 ------------
 * CASE IIa *
 ------------
 ************
 
 */


// Here I set out the functions which I define after main()

int CompMult(double* xr, double* xi, double* yr, double* yi, double* zr, double* zi);

int CompDiv(double* xr, double* xi, double* yr, double* yi, double* zr, double* zi);

int MatrixMult(double X[][2][2], double Y[][2][2], double Z[][2][2], int i, int j, int k);

int MatrixInv(double X[][2][2], double Z[][2][2], int i);

int VectorMult(double X[][2][2], double Y[][2][1], double Z[][2][1], int i, int j, int k);

int CMatrixMult(double Xr[][2][2], double Xi[][2][2] , double Yr[][2][2], double Yi[][2][2] , double Zr[][2][2], double Zi[][2][2], int i, int j, int k);

int CMatrixInv(double Xr[][2][2], double Xi[][2][2] , double Zr[][2][2], double Zi[][2][2], int i, int j);

int CVectorMult(double Xr[][2][2], double Xi[][2][2] , double Yr[][2][1], double Yi[][2][1] , double Zr[][2][1], double Zi[][2][1], int i, int j, int k);

int CMatrixDiagInv(double Xr[][2][2], double Xi[][2][2] , double Zr[][2][2], double Zi[][2][2], int i, int j);

int CompSQRT(double* ar, double* ai , double* br, double* bi);

int FifthOrderExtrap(double var[], double r[], double C);

int MeasureGrid(int* J_new_add, double* ratio_max_add);

int MakeGrid(int J_new, double ratio_max);

int FunctionF(double* xadd, double* f);

int FunctionG(double* xadd, double* g);

int FunctionY(double* xadd, double* y);


int main()
{
    
    int J,k,z;
    double C, parameter_nonad, location_nonad, prefactor;
    
    
    // J must be defined and given a value before any arrays which need it are defined, or else you'll get a segmentation fault because the arrays won't know how big they are.
    // So J is given a temporary value here until some MESA data is actually read in and use to define the size of the vectors instead.
    
    J=17500; // At the moment, it seems that this produces a seg fault when J >= 17500, which seems to be a result of filling up a memory limit, potentially the RAM?  But that would seem unlikely... It was helped by changing the dummy matrices to be [1][2][2] instead of [J][2][2], so it is a total memory issue rather than the memory used by any given array.
    z=J-1;
    
    
    // This section is to do with reading input from a file
    // If the mass of the star changes in the input file, you need to make sure that Mstar and D are changed appropriately
    
    ifstream infile;
    infile.open("Input/profiles_Henyey_conv_Ubuntu_1_0_alpha2H_conv_vel.txt");
    //"Input/profiles_Henyey_conv_Ubuntu_1_0_alpha2Hs_conv_vel.txt"
    //"Input/profiles_Henyey_conv_Ubuntu_1_4.txt" (beware of it not having conv_vel in it)
    
    k=0;
    
    double input;
    int no_of_lines;
    string line;
    
    no_of_lines = 0;
    
    ifstream linefile;
    linefile.open("Input/profiles_Henyey_conv_Ubuntu_1_0_alpha2H_conv_vel.txt");
    
    getline(linefile,line);
    while (linefile) {
        if (line != "") {
            no_of_lines = no_of_lines + 1;
        }
        
        getline(linefile,line);
        
    }
    
    linefile.close();
    
    cout << "The number of non-empty lines in this file is: " << no_of_lines << "\n" ;
    
    
    
    
    J = no_of_lines;
    
    
    
    
    z = J-1;
    
    
    // Here the input arrays are defined, and J cannot be changed again
    
    double zone[J], lnT[J], lnRho[J], grav[J], radius_cm[J], rmid_cm[J];
    double temperature[J], rho[J], pressure[J], grada[J], cp[J], chiRho[J], chiT[J];
    double opacity[J], dkap_dlnrho_face[J], dkap_dlnT_face[J], flux[J], brunt_A[J];
    double K[J], rho_face[J], scale_height_cm[J], pressure_scale_height_cm[J];
    double cv[J], gamma1[J], gamma3[J], conv_L_div_L[J], flux_tot[J], conv_vel[J];
    
    
    
    cout << "FLAG - about to take the input\n";
    
    
    // Get the first line's input
    infile >> input;
    
    // Then process the first line's input, then try to get another line.
    // This works because infile==true if the last thing it tried to read successfully gave it something.
    // Therefore, the while loop will only end when the condition is tested after the read has failed.
    // Therefore we need to keep the reading of line m and the processing of line m separated by a test of the condition
    while (infile) {
        
        zone[z] = input;
        
        infile >> input;
        
        lnT[z] = input*2.30258509299404568401799; // We need ln instead of log, so we use ln(x)= ln(10) log(x), where ln(10) = 2.30258509299404568401799
        
        infile >> input;
        
        lnRho[z] = input*2.30258509299404568401799;
        
        infile >> input;
        
        grav[z] = input;
        
        infile >> input;
        
        //5
        radius_cm[z] = input;
        
        infile >> input;
        
        rmid_cm[z] = input;
        
        infile >> input;
        
        temperature[z] = input;
        
        infile >> input;
        
        rho[z] = input;
        
        infile >> input;
        
        pressure[z] = input;
        
        infile >> input;
        
        //10
        grada[z] = input;
        
        infile >> input;
        
        cp[z] = input;
        
        infile >> input;
        
        chiRho[z] = input;
        
        infile >> input;
        
        chiT[z] = input;
        
        infile >> input;
        
        opacity[z] = input;
        
        infile >> input;
        
        //15
        dkap_dlnrho_face[z] = input;
        
        infile >> input;
        
        dkap_dlnT_face[z] = input;
        
        infile >> input;
        
        // This is to see if messing with the background flux makes much of a difference
        flux[z] = input; // *(1.0+0.0009*sin(10000*(radius_cm[z]/radius_cm[J-1])*(radius_cm[z]/radius_cm[J-1])*(radius_cm[z]/radius_cm[J-1])*(radius_cm[z]/radius_cm[J-1])*(radius_cm[z]/radius_cm[J-1])*(radius_cm[z]/radius_cm[J-1])));
        
        infile >> input;
        
        brunt_A[z] = input;
        
        infile >> input;
        
        // The following variables don't need to be interpolated, as they are just used for adding in extra cells for when non-adiabaticity becomes important
        
        scale_height_cm[z] = input;
        
        infile >> input;
        
        //20
        pressure_scale_height_cm[z] = input;
        
        infile >> input;
        
        cv[z] = input;
        
        infile >> input;
        
        gamma1[z] = input;
        
        infile >> input;
        
        gamma3[z] = input;
        
        infile >> input;
        
        conv_L_div_L[z] = input;
        
        infile >> input;
        
        conv_vel[z] = input;
        
        infile >> input; // This is to get the input for the zone for the next loop (and therefore to decide whether or not to do the next loop, too)
        
        
        
        K[z] = 4.0 * 7.565767e-15 * 2.99792458e10 * temperature[z]*temperature[z]*temperature[z] / (3.0 * opacity[z] * rho[z]);
        
        // Here we record the total flux (both convective and radiative)
        
        flux_tot[z] = flux[z];
        
        // Here flux is modified in order to just account for the radiative flux, not ALL of the flux
        
        flux[z] = flux[z] * (1.0 - conv_L_div_L[z]);
        
        
        z = z - 1;
    }
    
    /*
     // Here I add in the centre-most zone
     
     C = 100000; // This is a measure of how many times close to x=0 the outside edge of the centre-most zone will be (for reference, 100 wil be the standard size)
     
     radius_cm[0] = radius_cm[1]/C;
     rmid_cm[0] = rmid_cm[1]/C;
     zone[0] = zone[1] + 1;
     
     
     FifthOrderExtrap(lnT, rmid_cm, C);
     FifthOrderExtrap(lnRho, rmid_cm, C);
     FifthOrderExtrap(grav, rmid_cm, C);
     FifthOrderExtrap(temperature, rmid_cm, C);
     FifthOrderExtrap(rho, rmid_cm, C);
     FifthOrderExtrap(pressure, rmid_cm, C);
     FifthOrderExtrap(grada, rmid_cm, C);
     FifthOrderExtrap(cp, rmid_cm, C);
     FifthOrderExtrap(chiRho, rmid_cm, C);
     FifthOrderExtrap(chiT, rmid_cm, C);
     FifthOrderExtrap(opacity, rmid_cm, C);
     FifthOrderExtrap(brunt_A, rmid_cm, C);
     FifthOrderExtrap(K, rmid_cm, C);
     
     FifthOrderExtrap(dkap_dlnrho_face, radius_cm, C);
     FifthOrderExtrap(dkap_dlnT_face, radius_cm, C);
     FifthOrderExtrap(flux, radius_cm, C);
     FifthOrderExtrap(rho_face, radius_cm, C);
     
     */
    
    
    
    
    
    
    
    
    
    
    infile.close();
    
    double r_a, r_b, r_c;
    
    for (k=0; k < J-1; k = k+1) {
        
        // Here r_a and r_b are briefly being used as stand-ins for delta[0] and delta[1] respectively
        
        r_a = (rmid_cm[k+1] - radius_cm[k])/(rmid_cm[k+1] - rmid_cm[k]);
        r_b = (radius_cm[k] - rmid_cm[k])/(rmid_cm[k+1] - rmid_cm[k]);
        
        rho_face[k] = r_a*rho[k] + r_b*rho[k+1];
        
        
    }
    
    // rho_face
    /* This extrapolates to get the density at the outermost face, using:
     
     x = ( r_c*a - r_a*c + (r_a/r_b)*((r_c - r_a)/(r_c - r_b))*( r_b*c - r_c*b ) ) / ( r_c - r_a - (r_a/r_b)*((r_c - r_a)/(r_c - r_b))*( r_c - r_b ) )
     
     where:
     
     a = rho[J-1]
     b = rho[J-2]
     c = rho[J-3]
     
     r_a = radius_cm[J-1] - rmid_cm[J-1]
     r_b = radius_cm[J-1] - rmid_cm[J-2]
     r_c = radius_cm[J-1] - rmid_cm[J-3]
     
     */
    
    r_a = radius_cm[J-1] - rmid_cm[J-1];
    r_b = radius_cm[J-1] - rmid_cm[J-2];
    r_c = radius_cm[J-1] - rmid_cm[J-3];
    
    rho_face[J-1] = ( r_c*rho[J-1] - r_a*rho[J-3] + (r_a/r_b)*((r_c - r_a)/(r_c - r_b))*( r_b*rho[J-3] - r_c*rho[J-2] ) ) / ( r_c - r_a - (r_a/r_b)*((r_c - r_a)/(r_c - r_b))*( r_c - r_b ) );
    
    
    
    
    
    
    
    
    /*
     
     Here are some constants, as used in MESA, which can be found in mesa-r9575/const/public/const_def.f90
     Just thought they could be handy to have around.
     
     
     standard_cgrav ! = 6.67428d-8
     ! gravitational constant (g^-1 cm^3 s^-2)
     planck_h ! = 6.62606896D-27
     ! Planck's constant (erg s)
     hbar ! = planck_h / (2*pi)
     qe ! = 4.80320440D-10
     ! electron charge (esu == (g cm^3 s^-2)^(1/2))
     avo ! = 6.02214179d23
     ! Avogadro's constant (mole^-1)
     clight ! = 2.99792458d10
     ! speed of light in vacuum (cm s^1)
     kerg ! = 1.3806504D-16
     ! Boltzmann's constant (erg K^-1)
     
     boltz_sigma ! = 5.670400D-5
     ! boltzmann's sigma = crad*clight/4 (erg cm^-2 K^-4 s^-1)
     crad ! = boltz_sigma*4/clight = 7.565767d-15 (erg cm^-3 K^-4)
     ! radiation density constant, a (erg cm^-3 K^-4); Prad ! = crad * T^4 / 3
     
     
     
     ! astronomical constants
     ! solar age, L, and R values from Bahcall et al, ApJ 618 (2005) 1049-1056.
     msol ! = 1.9892d33  ! solar mass (g)
     rsol ! = 6.9598d10 ! solar radius (cm)
     lsol ! = 3.8418d33  ! solar luminosity (erg s^-1)
     agesol ! = 4.57d9  ! solar age (years)
     Msun ! = msol
     Rsun ! = rsol
     Lsun ! = lsol
     Msun33 ! = msol*1d-33
     Rsun11 ! = rsol*1d-11
     Lsun33 ! = lsol*1d-33
     Teffsol ! = 5777d0 ! temperature (k)
     loggsol ! = 4.4378893534131256d0 ! log surface gravity ! log(g/(cm s^-2))
     teffsun ! = teffsol
     loggsun ! = loggsol
     mbolsun ! = 4.746 ! Bolometric magnitude of the Sun
     ly ! = 9.460528d17 ! light year (cm)
     pc ! = 3.261633d0 * ly ! parsec (cm)
     secyer ! 3.1558149984d7 ! seconds per year
     dayyer ! 365.25 ! days per year
     m_earth ! = 5.9764d27 ! earth mass (g)
     r_earth ! = 6.37d8 ! earth radius (cm)
     au ! = 1.495978921d13 ! astronomical unit (cm)
     m_jupiter ! = 1.8986d30 ! jupiter mass (g)
     r_jupiter ! = 6.9911d9 ! jupiter mean radius (cm)
     semimajor_axis_jupiter ! = 7.7857d13 ! jupiter semimajor axis (cm)
     
     
     */
    
    
    double R, m, l, f, G, mp, D, omega, Mstar, d_one, d_two, dp_dr_BC, dlnT_dr_BC,prop,dlnRho_dr_BC,flux_BC,Period;
    
    // This is set by the tidal forcing
    m = 2.0;
    l = 2.0;
    
    R = radius_cm[J-1];
    
    flux_BC = flux[J-1];
    
    
    
    G = 6.67428e-8;
    
    mp = 1.0; // Planetary mass in terms of Jupiter masses
    mp = 1.8986e30 * mp; // converted into g
    
    Mstar = 1.0; // Stellar mass in solar masses
    Mstar = Mstar * 1.9892e33; // Converted into g
    
    
    
    
    D = 0.0512; //*1.118587; // Orbital radius of planet in AU, 1.58740105197 = 4^(1/3), which acts to double the period of the orbit
    
    
    // Here I modify D in au
    
    //D = D*5;
    
    
    D = 1.495978921e13 * D; // Converted into cm
    
    // Asserting D in cm directly because I calculated it
    // 1M non-res = 9.47699e11
    // 1M resonant= 4.6495e11
    //1.6M non-res= 1.13104e12
    //1.6Mresonant= 5.5097e11
    //D = 9.47699e11;
    
    f = -(G * mp) / (4.0 * D * D * D);
    
    cout << "f = " << f << "\n";
    
   
    
    
    
    
    omega = sqrt(G * ((Mstar + mp)/(D*D*D)));
    
    
    // Here I artificially modify JUST omega, and nothing else -- this is an unphysical situation because it leaves the orbital distance and masses unchanged, but alters the frequency in isolation
    
    omega = omega*1.0000;
    
    
    
    //f = - (mp * omega * omega) / (4.0 * Mstar);
    
    //cout << "*\n*\n*\n*\n*\n*\n*\n" ;
    //cout << "f: " << f << "\n";
    //cout << "*\n*\n*\n*\n*\n*\n*\n" ;
    
    
    
    
    
    cout << "omega = " << omega << "\n\n";
    
    cout << "\n\n" << "for JP, f = " << sqrt(4.0*omega*omega*R*R*R/(G*Mstar)) << "\n\n";
    
    cout << "equilibrium xi / R at surface = " << mp*R*R*R/(4.0*Mstar*D*D*D) << "\n";
    
    d_one = rmid_cm[J-1] - rmid_cm[J-2];
    
    d_two = rmid_cm[J-1] - rmid_cm[J-3];
    
    dp_dr_BC = (d_one + d_two)*pressure[J-1]/(d_one*d_two) - pressure[J-2]*d_two / ( d_one * (d_two - d_one) ) + pressure[J-3] * d_one / (d_two * (d_two - d_one) );
    
    dlnT_dr_BC = (d_one + d_two)*lnT[J-1]/(d_one*d_two) - lnT[J-2]*d_two / ( d_one * (d_two - d_one) ) + lnT[J-3] * d_one / (d_two * (d_two - d_one) );
    
    dlnRho_dr_BC = (d_one + d_two)*lnRho[J-1]/(d_one*d_two) - lnRho[J-2]*d_two / ( d_one * (d_two - d_one) ) + lnRho[J-3] * d_one / (d_two * (d_two - d_one) );
    
    cout << "dp_dr_BC = " << dp_dr_BC << "\n";
    
    cout << "dlnT_dr_BC = " << dlnT_dr_BC << "\n";
    
    cout << "dlnRho_dr_BC = " << dlnRho_dr_BC << "\n";
    
    
    
    /*
     This is to try to increase resolution in the centre-most zones.
     
     L = number of zones per group
     N = max number of sub-zones per zone (e.g. the highest resolution increase)
     E = extra number of cells
     
     
     
     E = L SUM_{k=1}^{k = Nmax - 1}(k) = ( L * Nmax * ( Nmax - 1 ) ) / 2
     
     Then the rest are all left untouched.
     
     Because of all of the extra zones, J needs to be re-defined, and this could cause troubles with the other variables when you try to upscale this.  ***********BEWARE!!!***************
     
     
     
     
     */
    
    
    // J is redefined here, and the old value is stored as Jold
    
    
    
    // count[kold] = the number of cells used to get to the outer edge of cell kold (including that edge).  Therefore J = count[Jold-1], and the extra number of cells = J - Jold
    int N[J], count[J], Jold, k_start, k_end;
    
    int U, E, L, Nmax, n;
    // these variables pertain to the cell size choice
    double dx_max_general, dx_max_surface, dx_max, delta_surface;
    
    int dx_max_surface_tracker;
    
    Jold = J;
    
    count[0] = 1;
    
    dx_max_surface_tracker = 0;
    
    
    // This changes how much non-adiabaticity comes into play: 1.0 = non-adiabatic, larger makes it less so (approaches the adiabatic limit)
    prop = 1.0;
    
    
    
    
    // This bit opens the file to write the dx_max_surface data into
    ofstream dx_max_surface_file;
    dx_max_surface_file.open("Output/dx_max_surface_Ubuntu.dat", ios::out);
    
    // This sets the precision at which values are printed at to the named file output
    //matrixfile.precision(10);
    
    
    
    
    
    // This sets N[0] and N[Jold-1] =0, as the first and last cells will be left alone
    N[0] = 0;
    N[Jold-1] = 0;
    
    
    // This sums over the extra cells needed to achieve at least the minimum resolution, given by the maximum cell size: dx_max
    // It only goes up to k < Jold - 2 because we don't want to split the final cell, that is the cell from radius_cm[Jold-2] to radius_cm[Jold-1]
    
    for (k=1; k < Jold-1; k = k + 1) {
        
        // delta_surface is the cell width for the non-adiabaticity to be important
        
        delta_surface = sqrt(( (gamma3[k]-1) * flux[k] * scale_height_cm[k] )/( pressure[k] * omega * prop));
        
        // by defining dx_max inside this loop, we can make it a function of radius_cm, for instance
        
        dx_max_general = 0.000015 + 0.0003*(radius_cm[k]/R)*(radius_cm[k]/R);
        
        dx_max_surface = (delta_surface/R) * 1.0;
        
        
        if ( (dx_max_general < dx_max_surface) || ((radius_cm[k]/R)<0.9985)) {
            
            dx_max = dx_max_general;
            
        } else {
            
            dx_max_surface_tracker = dx_max_surface_tracker + 1;
            
            dx_max = dx_max_surface;
            
        }
        
        N[k] = ((radius_cm[k] - radius_cm[k-1])/(R*dx_max));
        
        count[k] = count[k-1] + 1 + N[k];
        
        //cout << "At k = " << k << ", N = " << N << "\n";
        
        dx_max_surface_file << radius_cm[k]/R << "\t" << dx_max_surface << "\n";
        
        
    }
    
    
    dx_max_surface_file.close();
    
    
    
    /*
     count[Jold-1] = count[Jold-2] + 1;
     
     cout << "\n\ncount[Jold-4] = " << count[Jold-4] << "\n\n";
     
     cout << "\n\ncount[Jold-3] = " << count[Jold-3] << "\n\n";
     
     cout << "\n\ncount[Jold-2] = " << count[Jold-2] << "\n\n";
     
     cout << "\n\ncount[Jold-1] = " << count[Jold-1] << "\n\n";
     
     // This gives the new value of J, by taking into account all of the extra cells needed
     J = count[Jold-1];
     */
    
    
    // This section calculates the locations of the new set of grid points according to your choice of functions
    // Currently both functions must be functions of x, but this may be changed in the future
    int J_new;
    double ratio_max, function_f, function_g, y, x, r_grid_cm;
    
    x = 0.7;
    
    FunctionF(&x,&function_f);
    FunctionG(&x,&function_g);
    FunctionY(&x,&y);
    
    cout << "x = " << x << "\n";
    cout << "f = " << function_f << "\n";
    cout << "g = " << function_g << "\n";
    cout << "y = " << y << "\n";
    
    ratio_max = 1.1;
    
    MeasureGrid(&J_new, &ratio_max);
    
    cout << "J_new = " << J_new << "\n";
    
    
    
    MakeGrid(J_new, ratio_max);
    
    
    
    
    // Here I override the previous redefinition of J, and define it according to whatever functions I have chosen to use
    J = J_new;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    // This sets the precision at which values are printed at (potentially for both the terminal and files, that is not yet clear)
    cout.precision(15);
    
    
    
    cout << "FLAG - about to introduce largest arrays\n";
    
    
    /*
     
     In order to reduce the burden on the memory, all other than alpha, gamma, u, v, RECu, RECv and RECc will only be defined for one value of i at a time, so everything used to calculate alphas and gammas must be in the loop which involves A, B etc being defined.  And you'll need to be careful with the outer BCs, although if they are defined for J-1 last, then it should be okay as they will remain defined as such outside of the for loop.
     
     */
    
    
    
    // These are the matrices in the oscillation equations
    double Ar[1][2][2],Ai[1][2][2],Cr[1][2][2],Ci[1][2][2],Dr[1][2][2],Di[1][2][2];
    double Er[1][2][2],Ei[1][2][2],Fr[1][2][2],Fi[1][2][2],Hr[1][2][2],Hi[1][2][2];
    
    // These are the vectors in the oscillation equations
    double ur[J][2][1],ui[J][2][1],vr[J][2][1],vi[J][2][1],Mr[1][2][1],Mi[1][2][1],Nr[1][2][1],Ni[1][2][1];
    
    // This is the matrix involved in the u-v relation
    double alphar[J][2][2],alphai[J][2][2];
    
    // This is the vector involved in the u-v relation
    double gammar[J][2][1],gammai[J][2][1];
    
    // In order to work back from u_{i+1} and v_{i+1} to v_{i}, 3 other arrays will need to retain their total memory, to give us: v_{i} = RECu u_{i+1} + RECv v_{i+1} + RECc
    double RECur[J][2][2], RECui[J][2][2], RECvr[J][2][2], RECvi[J][2][2], RECcr[J][2][1], RECci[J][2][1];
    
    
    cout << "FLAG - part way through the introductions\n";
    
    // Here some dummy matrices and vectors must be introduced in order to manage complex calculations later on, including some named and defined ones, and some total dummies
    
    // These are the defined ones
    double Pr[1][2][2], Pi[1][2][2], Qr[1][2][2], Qi[1][2][2], Rr[1][2][2], Ri[1][2][2];
    
    cout << "FLAG - done all the J x thing x thing arrays now\n"; // This flag is not reached - therefore I think it's the arrays just above which are currently pushing things over the edge - can I rearrange it so that these are just 1x2x2 matrices, by defining them within the loop?  (Need to check that nothing crosses over to another cell, but I think it should be all ok.)
    
    // These are the total dummies, with the naming convention as: dum = dummy; M / V = matrix / vector.  Beware when using these, as they will come in with values already attached, so make sure that the first thing that you do is to define them according to your particular purpose at that time.
    double dumMAr[1][2][2], dumMAi[1][2][2], dumMBr[1][2][2], dumMBi[1][2][2], dumMCr[1][2][2], dumMCi[1][2][2], dumMDr[1][2][2], dumMDi[1][2][2], dumMEr[1][2][2], dumMEi[1][2][2];
    double dumVAr[1][2][1], dumVAi[1][2][1], dumVBr[1][2][1], dumVBi[1][2][1], dumVCr[1][2][1], dumVCi[1][2][1], dumVDr[1][2][1], dumVDi[1][2][1], dumVEr[1][2][1], dumVEi[1][2][1];
    
    cout << "Most of the way through the introductions\n";
    
    
    // These are matrices involved with the outer boundary conditions, and as they only apply at the boundary they don't need to be the same size as the general arrays, but only need to include one 2x2 matrix each, or one 2D vector.  The matrices and vectors are kept separate for clarity.  The first row is the matrices in the boundary condition equation, the second is a series of defined matrices to simplify the expression for v_{J-1}, and the third row is the vector from the boundary condition equation.
    double etar[1][2][2], etai[1][2][2], mur[1][2][2], mui[1][2][2], nur[1][2][2], nui[1][2][2];
    double BCar[1][2][2], BCai[1][2][2], BCbr[1][2][2], BCbi[1][2][2], BCcr[1][2][2], BCci[1][2][2], BCdr[1][2][2], BCdi[1][2][2], BCer[1][2][2], BCei[1][2][2], BCfr[1][2][2], BCfi[1][2][2];
    double xr[1][2][1], xi[1][2][1];
    
    cout << "FLAG - large arrays all introduced\n";
    
    // These are the rescaling parameters, and to keep things general, they are going to be complex.
    double T_ar, T_ai, T_br, T_bi, T_cr, T_ci, T_dr, T_di, dummyr, dummyi;
    
    // This is just another dummy index
    int i;
    
    
    /*
     
     The basic equations are as follows:
     
     A_{i,i+1} u_{i} + C_{i,i+1} u_{i+1} + D_{i,i+1} v_{i+1} = M_{i,i+1}
     
     E_{i,i+1} u_{i} + F_{i,i+1} v_{i}   + H_{i,i+1} v_{i+1} = N_{i,i+1}
     
     u_{i} + alpha_{i} v_{i} + gamma_{i} = 0
     
     The equations are arranged as such because of where the quantities are defined within the cell.
     u = (\xi_{r} , F_{r}), both of which are defined at the outer edge of the cell, whereas v = (p ,T), both of which are defined in the middle of the cell.
     As such, the first vector equation is valid at the boundary between region i and region i+1, so can only contain /xi_{i} and F_{i}, but the gradients of p and T, which will be valid at this point (as they will be of the form (p_{i+1} - p_{i})/(rmid_{i+1} - rmid_{i}) ).
     The second vector equation is valid in the middle of cell i+1, so can involve p_{i+1}, T_{i+1} and the gradients of \xi and F.
     
     Using the above equations to elimate u_{i} and v_{i}, and comparing coefficents with the equation u_{i+1} + alpha_{i+1} v_{i+1} + gamma_{i+1} = 0 gives the following relations (CASE II):
     
     alpha_{i+1} = C^{-1} [ D - A ( alpha_{i} F^{-1} E - 1 )^{-1} alpha_{i} F_{-1} H ]
     
     gamma_{i+1} = C^{-1} [ A ( alpha_{i} F^{-1} E - 1 )^{-1} ( gamma_{i} + alpha_{i} F^{-1} N ) - M ]
     
     where a lack of subscript implies _{i, i+1}
     
     
     Using the top two equations to eliminate u_{i} gives an expression for v_{i} in terms of u_{i+1} and v_{i+1} as (CASE a):
     
     u_{i} = A^{-1} [ M - C u_{i+1} - D v_{i+1} ]3
     
     Which can be recast as:
     
     u_{i} = RECu u_{i+1} + RECv v_{i+1} + RECc
     
     where
     
     RECu = - A^{-1} C
     RECv = - A^{-1} D
     RECc = A^{-1} M
     
     
     
     
     To minimise needlessly repeated calculations, these equations will be re-written in terms of recurring blocks of matrices, such as:
     
     P = E A^{-1}
     
     Q = alpha_{i} F^{-1}
     
     R = [ alpha_{i} F^{-1} E - 1 ]^{-1} = [ Q E - 1 ]^{-1}
     
     
     Using these expression in the equations for alpha_{i+1}, gamma_{i+1} and v_{i} gives us:
     
     alpha_{i+1} = C^{-1} [ D - A R Q H ]
     gamma_{i+1} = C^{-1} [ A R ( gamma_{i} + Q N ) - M ]
     v_{i} = F^{-1} [ N + P ( C u_{i+1} - M ) + ( P D - H ) v_{i+1} ]
     
     
     
     To help with evaluating gradients to do with p and T in the centre of zone J-1, the following are used:
     
     dp_dr_BC = (dp/dr)[J-1] = (d_one + d_two)*pressure[J-1]/(d_one*d_two) - pressure[J-1]*d_two / ( d_one * (d_two - d_one) ) + pressure[J-3] * d_one / (d_two * (d_two - d_one) )
     
     dlnT_dr_BC = (dlnT/dr)[J-1] = (d_one + d_two)*lnT[J-1]/(d_one*d_two) - lnT[J-1]*d_two / ( d_one * (d_two - d_one) ) + lnT[J-3] * d_one / (d_two * (d_two - d_one) )
     
     where
     
     d_one = rmid_cm[J-1] - rmid_cm[J-2]
     
     d_two = rmid_cm[J-1] - rmid_cm[J-3]
     
     
     */
    
    /*
     
     
     (1/r^2)*d(r^2*a)/dr+c+d   =
     a+(1/r^2)*d(r^2*b)/dr+c+d =
     b+c+d+d'   =
     a+c+c'+d   =
     
     These are built around giving the solution:
     a = (r/R)*sin(r/R)
     b = sqrt(r/R)*sin(sqrt(200.0)*r/R)
     c = 50.0*sin(sqrt(17.0)*r/R)/(sqrt(17.0)*r/R)
     d = cosh(r/R) + 3.0*cos(sqrt(23.0)*r/R)
     
     (1/r^2)*d(r^2*a)/dr = (3.0/R)*sin(r/R) + (r/(R*R))*cos(r/R)
     (1/r^2)*d(r^2*b)/dr = (5.0/(2.0*R*sqrt(r)))*sin(sqrt(200)*r/R) + (sqrt(200*r)/(R*R))*cos(sqrt(200)*r/R)
     dp/dr = (50.0/r)*( cos(sqrt(17.0)*r/R) - (R/(sqrt(17.0)*r))*sin(sqrt(17.0)*r/R) )
     dT/dr = sinh(r/R)/r - (sqrt(207)/R)*sin(sqrt(23)*r/R)
     
     
     
     The BCs are:
     a+c   = sin(1) + cos(1)
     a+b+d = -cos(1)
     
     */
    
    // Before iterating through to get all the alpha matrices and gamma vectors, we must use the boundary conditions at the centre: alpha_{0} = 0; gamma_{0} = 0
    
    alphar[0][0][0] = 0;
    alphar[0][0][1] = 0;
    alphar[0][1][0] = 0;
    alphar[0][1][1] = 0;
    
    gammar[0][0][0] = 0;
    gammar[0][1][0] = 0;
    
    alphai[0][0][0] = 0;
    alphai[0][0][1] = 0;
    alphai[0][1][0] = 0;
    alphai[0][1][1] = 0;
    
    gammai[0][0][0] = 0;
    gammai[0][1][0] = 0;
    
    cout << "alpha_{0} and gamma_{0} have been set to 0, which is equivalent to using the inner BCs \n";
    
    int xcount, labelminx, labelmink, kold, kmin;
    
    double count_doub,N_doub,k_doub,x_convection;
    
    // These are defined to keep only two variables at any one time - k and k+1, which will be re-written for each new zone
    double lnT_HR[2], lnRho_HR[2], grav_HR[2], radius_cm_HR[2], rmid_cm_HR[2];
    double temperature_HR[2], rho_HR[2], pressure_HR[2], grada_HR[2], cp_HR[2], chiRho_HR[2], chiT_HR[2];
    double opacity_HR[2], dkap_dlnrho_face_HR[2], dkap_dlnT_face_HR[2], flux_HR[2], brunt_A_HR[2];
    double K_HR[2], rho_face_HR[2], gamma1_HR[2], flux_tot_HR[2], conv_L_div_L_HR[2], conv_vel_HR[2], scale_height_cm_HR[2];
    
    // This is defined for the k and k+1 values of delta which help with linear interpolation: F_{at rmid[1]} = delta[0]*flux[0] + delta[1]*flux[1]
    // The third one is for interpolating cell-central variables, because you need to involve k+2 variables for that
    double delta[3];
    
    // This is an exception, as this needs to be output at the end for plotting purposes
    double radius_cm_HR_output[J],rmid_cm_HR_output[J],flux_HR_output[J],pressure_HR_output[J],temperature_HR_output[J],test[J],rho_HR_output[J],eta_Terquem_output[J],grav_HR_output[J],rho_face_HR_output[J],Fprime_for_T[J],Fprime_for_p[J],opacity_HR_output[J],K_HR_output[J],chiRho_HR_output[J],chiT_HR_output[J],brunt_A_HR_output[J],cp_HR_output[J],grada_HR_output[J],conv_L_div_L_HR_output[J],flux_tot_HR_output[J],dkap_dlnT_face_HR_output[J],dkap_dlnrho_face_HR_output[J],gamma1_HR_output[J];
    
    
    
    
    
    
    // For testing, these variables are introduced.
    double a[2], b[2], c[2], d[2];
    
    
    int sum;
    
    sum = 0;
    
    
    
    
    
    // This bit opens the file to write the Matrix data into
    ofstream matrixfile;
    matrixfile.open("Output/Matrices_rescaled_Voscish.dat", ios::out);
    
    // This sets the precision at which values are printed at to the named file output
    matrixfile.precision(10);
    
    matrixfile << "radius_cm_HR_output[k]/R" << "\t\t\t" << "Ar[0][0][0]" << "\t\t" << "Ar[0][0][1]" << "\t\t" << "Ar[0][1][0]" << "\t\t" << "Ar[0][1][1]" << "\t\t" << "Cr[0][0][0]" << "\t\t" << "Cr[0][0][1]" << "\t\t" << "Cr[0][1][0]" << "\t\t" << "Cr[0][1][1]" << "\t\t" << "Dr[0][0][0]" << "\t\t" << "Dr[0][0][1]" << "\t\t" << "Dr[0][1][0]" << "\t\t" << "Dr[0][1][1]" << "\t\t" << "Er[0][0][0]" << "\t\t" << "Er[0][0][1]" << "\t\t" << "Er[0][1][0]" << "\t\t" << "Er[0][1][1]" << "\t\t" << "Fr[0][0][0]" << "\t\t" << "Fr[0][0][1]" << "\t\t" << "Fr[0][1][0]" << "\t\t" << "Fr[0][1][1]" << "\t\t" << "Hr[0][0][0]" << "\t\t" << "Hr[0][0][1]" << "\t\t" << "Hr[0][1][0]" << "\t\t" << "Hr[0][1][1]" << "\t\t" << "Mr[0][0][0]" << "\t\t" << "Mr[0][1][0]" << "\t\t" << "Nr[0][0][0]" << "\t\t" << "Nr[0][1][0]" << "\t\t\t" << "detA" << "\t\t" << "detC" << "\t\t" << "detD" << "\t\t" << "detE" << "\t\t" << "detF" << "\t\t" << "detH" << "\n";
    
    
    
    
    
    // This initialises all of the array_HR variables (as the k=0 case is unaltered)
    radius_cm_HR_output[0] = radius_cm[0];
    rmid_cm_HR_output[0] = rmid_cm[0];
    
    radius_cm_HR[0] = radius_cm[0];
    rmid_cm_HR[0] = rmid_cm[0];
    flux_HR[0] = flux[0];
    dkap_dlnrho_face_HR[0] = dkap_dlnrho_face[0];
    dkap_dlnT_face_HR[0] = dkap_dlnT_face[0];
    opacity_HR[0] = opacity[0];
    lnRho_HR[0] = lnRho[0];
    rho_HR[0] = rho[0];
    cp_HR[0] = cp[0];
    temperature_HR[0] = temperature[0];
    lnT_HR[0] = lnT[0];
    pressure_HR[0] = pressure[0];
    grav_HR[0] = grav[0];
    K_HR[0] = K[0];
    chiRho_HR[0] = chiRho[0];
    chiT_HR[0] = chiT[0];
    grada_HR[0] = grada[0];
    brunt_A_HR[0] = brunt_A[0];
    rho_face_HR[0] = rho_face[0];
    gamma1_HR[0] = gamma1[0];
    conv_L_div_L_HR[0] = conv_L_div_L[0];
    conv_vel_HR[0] = conv_vel[0];
    scale_height_cm_HR[0] = scale_height_cm[0];
    flux_tot_HR[0] = flux_tot[0];
    
    
    
    
    
    
    // Here I  initialise the array r_grid_cm[3], which will contain r_grid_cm at k-1, k and k+1. But if k=0, then k-1 doesn't work, so it is just 0.0;
    // This information is being read from the file which has the locations listed in order, from the centre of the star
    
    ifstream gridfile;
    gridfile.open("Memory/grid_new_centre.txt");
    gridfile.precision(15);
    
    // This takes the first thing in the file, which is for k=0
    gridfile >> input;
    
    double r_grid_array[3];
    
    
    // We imagine we are starting with k=-1, so we have the oth and 1st elements = 0.0, then element 2 has the input
    
    r_grid_array[0] = 0.0;
    
    r_grid_array[1] = 0.0;
    
    r_grid_array[2] = input;
    
    cout << "r_grid_array is:\n";
    cout << r_grid_array[0] << "\n";
    cout << r_grid_array[1] << "\n";
    cout << r_grid_array[2] << "\n";
    
    
    //
    //
    //
    //
    //
    //
    //
    //
    //
    // This is the start of the big for loop in which all of the matrix stuff is setup and first used
    //
    //
    //
    //
    //
    //
    
    for (k=0; k < J-1; k=k+1) {
        
        
        
        // Here we shuffle along the variable arrays for this new value of k
        
        // This takes in the next line of the gridfile, which will be the k+1 value
        gridfile >> input;
        
        r_grid_array[0] = r_grid_array[1];
        
        r_grid_array[1] = r_grid_array[2];
        
        r_grid_array[2] = input;
        
        
        
        
        // This defines the variables for k, which will then all need to repeated because we also need to know the variables for k+1 at the same time
        
        // Here I define the location of the grid point as r_grid_cm
        r_grid_cm = R*r_grid_array[1];
        
        if ( r_grid_cm < radius_cm[0] ) {
            // This sets us to use cells 0 and 1 of the old grid if we need to get anything smaller than radius_cm[0]
            // Note that this requires extrapolation, sadly.
            kold = 1;
            
        } else {
            
            if ( r_grid_cm >= radius_cm[Jold-1] ) {
                // This sets us to use the outer cell, and nothing beyond it (prevents segfaults)
                // Whilst we shouldn't need to extrapolate at the surface for the cell-surface variables, this will be needed for the cell-mid variables
                kold = Jold-1;
                
            } else {
                
                for (kold = 0; radius_cm[kold] < r_grid_cm; kold = kold + 1) {
                    
                    // I don't need anything in here
                    // NB - kold is then the first point at which r_grid_cm has been exceeded
                    // therefore the interpolation is done between kold-1 and kold
                    
                }
                
            }
            
        }
        
        
        // Now we define the values for delta[0] and delta[1], for the cell-face variables
        // q[k] = delta[0]*q0[kold-1] + delta[1]*q0[kold]
        
        delta[0] = (radius_cm[kold] - r_grid_cm)/(radius_cm[kold] - radius_cm[kold-1]);
        
        delta[1] = (r_grid_cm - radius_cm[kold-1])/(radius_cm[kold] - radius_cm[kold-1]);
        
        
        // Here we then define the cell-face variables, because we've already taken care of all of the possible dodgy cases
        
        radius_cm_HR_output[k] = delta[0]*radius_cm[kold-1] + delta[1]*radius_cm[kold];
        
        
        
        radius_cm_HR[0] = delta[0]*radius_cm[kold-1] + delta[1]*radius_cm[kold];
        
        flux_HR[0] = delta[0]*flux[kold-1] + delta[1]*flux[kold];
        
        dkap_dlnrho_face_HR[0] = delta[0]*dkap_dlnrho_face[kold-1] + delta[1]*dkap_dlnrho_face[kold];
        
        dkap_dlnT_face_HR[0] = delta[0]*dkap_dlnT_face[kold-1] + delta[1]*dkap_dlnT_face[kold];
        
        opacity_HR[0] = delta[0]*opacity[kold-1] + delta[1]*opacity[kold];
        
        rho_face_HR[0] = delta[0]*rho_face[kold-1] + delta[1]*rho_face[kold];
        
        flux_tot_HR[0] = delta[0]*flux_tot[kold-1] + delta[1]*flux_tot[kold];
        
        conv_L_div_L_HR[0] = delta[0]*conv_L_div_L[kold-1] + delta[1]*conv_L_div_L[kold];
        
        conv_vel_HR[0] = delta[0]*conv_vel[kold-1] + delta[1]*conv_vel[kold];
        
        scale_height_cm_HR[0] = delta[0]*scale_height_cm[kold-1] + delta[1]*scale_height_cm[kold];
        
        //
        //
        //
        //
        //
        // Now we move on to the calculation for the cell-mid variables.
        //
        //
        //
        //
        // We first need to find where the middle of the cell is, and we reuse r_grid_new to store it
        
        
        
        r_grid_cm = 0.5*R*( r_grid_array[1] + r_grid_array[0] );
        
        
        
        
        // STILL NEED TO DEFINE rmid_cm_HR_output[k]
        rmid_cm_HR_output[k] = r_grid_cm;
        
        // Now we need to find kold for this cell-centre, according to the cell-centres of the old grid
        
        if ( r_grid_cm < rmid_cm[0] ) {
            // This sets us to use cells 0 and 1 of the old grid if we need to get anything smaller than radius_cm[0]
            // Note that this requires extrapolation, sadly.
            kold = 1;
            
        } else {
            
            if ( r_grid_cm >= rmid_cm[Jold-1] ) {
                // This sets us to use the outer cell, and nothing beyond it (prevents segfaults)
                // Whilst we shouldn't need to extrapolate at the surface for the cell-surface variables, this will be needed for the cell-mid variables
                kold = Jold-1;
                
            } else {
                
                for (kold = 0; rmid_cm[kold] < r_grid_cm; kold = kold + 1) {
                    
                    // I don't need anything in here
                    // NB - kold is then the first point at which r_grid_cm has been exceeded
                    // therefore the interpolation is done between kold-1 and kold
                    
                }
                
            }
            
        }
        
        // Now we define the values for delta[0] and delta[1], for the cell-centre variables
        // q[k] = delta[0]*q0[kold-1] + delta[1]*q0[kold]
        
        delta[0] = (rmid_cm[kold] - r_grid_cm)/(rmid_cm[kold] - rmid_cm[kold-1]);
        
        delta[1] = (r_grid_cm - rmid_cm[kold-1])/(rmid_cm[kold] - rmid_cm[kold-1]);
        
        
        // Here we define the cell-centre variables
        
        rmid_cm_HR[0] = delta[0]*rmid_cm[kold-1] + delta[1]*rmid_cm[kold];
        lnRho_HR[0] = delta[0]*lnRho[kold-1] + delta[1]*lnRho[kold];
        rho_HR[0] = delta[0]*rho[kold-1] + delta[1]*rho[kold];
        cp_HR[0] = delta[0]*cp[kold-1] + delta[1]*cp[kold];
        temperature_HR[0] = delta[0]*temperature[kold-1] + delta[1]*temperature[kold];
        
        lnT_HR[0] = delta[0]*lnT[kold-1] + delta[1]*lnT[kold];
        pressure_HR[0] = delta[0]*pressure[kold-1] + delta[1]*pressure[kold];
        grav_HR[0] = delta[0]*grav[kold-1] + delta[1]*grav[kold];
        K_HR[0] = delta[0]*K[kold-1] + delta[1]*K[kold];
        chiRho_HR[0] = delta[0]*chiRho[kold-1] + delta[1]*chiRho[kold];
        chiT_HR[0] = delta[0]*chiT[kold-1] + delta[1]*chiT[kold];
        grada_HR[0] = delta[0]*grada[kold-1] + delta[1]*grada[kold];
        brunt_A_HR[0] = delta[0]*brunt_A[kold-1] + delta[1]*brunt_A[kold];
        gamma1_HR[0] = delta[0]*gamma1[kold-1] + delta[1]*gamma1[kold];
        
        
        
        
        // Here we define the variables for k+1
        // It's the same process, except we use i=k+1, and replace all the "k"s with "i"s, and the "0"s with "1"s in the variables
        
        i = k+1;
        
        // Here I define the location of the grid point as r_grid_cm
        r_grid_cm = R*r_grid_array[2];
        
        if ( r_grid_cm < radius_cm[0] ) {
            // This sets us to use cells 0 and 1 of the old grid if we need to get anything smaller than radius_cm[0]
            // Note that this requires extrapolation, sadly.
            kold = 1;
            
        } else {
            
            if ( r_grid_cm >= radius_cm[Jold-1] ) {
                // This sets us to use the outer cell, and nothing beyond it (prevents segfaults)
                // Whilst we shouldn't need to extrapolate at the surface for the cell-surface variables, this will be needed for the cell-mid variables
                kold = Jold-1;
                
            } else {
                
                for (kold = 0; radius_cm[kold] < r_grid_cm; kold = kold + 1) {
                    
                    // I don't need anything in here
                    // NB - kold is then the first point at which r_grid_cm has been exceeded
                    // therefore the interpolation is done between kold-1 and kold
                    
                }
                
            }
            
        }
        
        
        // Now we define the values for delta[0] and delta[1], for the cell-face variables
        // q[k] = delta[0]*q0[kold-1] + delta[1]*q0[kold]
        
        delta[0] = (radius_cm[kold] - r_grid_cm)/(radius_cm[kold] - radius_cm[kold-1]);
        
        delta[1] = (r_grid_cm - radius_cm[kold-1])/(radius_cm[kold] - radius_cm[kold-1]);
        
        
        // Here we then define the cell-face variables, because we've already taken care of all of the possible dodgy cases
        
        radius_cm_HR_output[i] = delta[0]*radius_cm[kold-1] + delta[1]*radius_cm[kold];
        
        radius_cm_HR[1] = delta[0]*radius_cm[kold-1] + delta[1]*radius_cm[kold];
        flux_HR[1] = delta[0]*flux[kold-1] + delta[1]*flux[kold];
        dkap_dlnrho_face_HR[1] = delta[0]*dkap_dlnrho_face[kold-1] + delta[1]*dkap_dlnrho_face[kold];
        dkap_dlnT_face_HR[1] = delta[0]*dkap_dlnT_face[kold-1] + delta[1]*dkap_dlnT_face[kold];
        opacity_HR[1] = delta[0]*opacity[kold-1] + delta[1]*opacity[kold];
        rho_face_HR[1] = delta[0]*rho_face[kold-1] + delta[1]*rho_face[kold];
        flux_tot_HR[1] = delta[0]*flux_tot[kold-1] + delta[1]*flux_tot[kold];
        conv_L_div_L_HR[1] = delta[0]*conv_L_div_L[kold-1] + delta[1]*conv_L_div_L[kold];
        conv_vel_HR[1] = delta[0]*conv_vel[kold-1] + delta[1]*conv_vel[kold];
        scale_height_cm_HR[1] = delta[0]*scale_height_cm[kold-1] + delta[1]*scale_height_cm[kold];
        
        
        
        
        
        //
        //
        //
        //
        //
        // Now we move on to the calculation for the cell-mid variables.
        //
        //
        //
        //
        // We first need to find where the middle of the cell is, and we reuse r_grid_new to store it
        
        
        
        r_grid_cm = 0.5*R*( r_grid_array[2] + r_grid_array[1] );
        
        
        
        
        // STILL NEED TO DEFINE rmid_cm_HR_output[i]
        rmid_cm_HR_output[i] = r_grid_cm;
        
        // Now we need to find kold for this cell-centre, according to the cell-centres of the old grid
        
        if ( r_grid_cm < rmid_cm[0] ) {
            // This sets us to use cells 0 and 1 of the old grid if we need to get anything smaller than radius_cm[0]
            // Note that this requires extrapolation, sadly.
            kold = 1;
            
        } else {
            
            if ( r_grid_cm >= rmid_cm[Jold-1] ) {
                // This sets us to use the outer cell, and nothing beyond it (prevents segfaults)
                // Whilst we shouldn't need to extrapolate at the surface for the cell-surface variables, this will be needed for the cell-mid variables
                kold = Jold-1;
                
            } else {
                
                for (kold = 0; rmid_cm[kold] < r_grid_cm; kold = kold + 1) {
                    
                    // I don't need anything in here
                    // NB - kold is then the first point at which r_grid_cm has been exceeded
                    // therefore the interpolation is done between kold-1 and kold
                    
                }
                
            }
            
        }
        
        // Now we define the values for delta[0] and delta[1], for the cell-centre variables
        // q[i] = delta[0]*q0[kold-1] + delta[1]*q0[kold]
        
        delta[0] = (rmid_cm[kold] - r_grid_cm)/(rmid_cm[kold] - rmid_cm[kold-1]);
        
        delta[1] = (r_grid_cm - rmid_cm[kold-1])/(rmid_cm[kold] - rmid_cm[kold-1]);
        
        
        // Here we define the cell-centre variables
        
        rmid_cm_HR[1] = delta[0]*rmid_cm[kold-1] + delta[1]*rmid_cm[kold];
        lnRho_HR[1] = delta[0]*lnRho[kold-1] + delta[1]*lnRho[kold];
        rho_HR[1] = delta[0]*rho[kold-1] + delta[1]*rho[kold];
        cp_HR[1] = delta[0]*cp[kold-1] + delta[1]*cp[kold];
        temperature_HR[1] = delta[0]*temperature[kold-1] + delta[1]*temperature[kold];
        lnT_HR[1] = delta[0]*lnT[kold-1] + delta[1]*lnT[kold];
        pressure_HR[1] = delta[0]*pressure[kold-1] + delta[1]*pressure[kold];
        grav_HR[1] = delta[0]*grav[kold-1] + delta[1]*grav[kold];
        K_HR[1] = delta[0]*K[kold-1] + delta[1]*K[kold];
        chiRho_HR[1] = delta[0]*chiRho[kold-1] + delta[1]*chiRho[kold];
        chiT_HR[1] = delta[0]*chiT[kold-1] + delta[1]*chiT[kold];
        grada_HR[1] = delta[0]*grada[kold-1] + delta[1]*grada[kold];
        brunt_A_HR[1] = delta[0]*brunt_A[kold-1] + delta[1]*brunt_A[kold];
        gamma1_HR[1] = delta[0]*gamma1[kold-1] + delta[1]*gamma1[kold];
        
        
        
        
        
        
        
        
        
        
        
        // cout << "Flag with information: Jold=" << Jold << "   J=" << J << "   k=" << k << "   i=" << i << "   radius_cm_HR_output[k]/R=" << radius_cm_HR_output[k]/R << "   radius_cm_HR_output[i]/R=" << radius_cm_HR_output[i]/R << "   kold=" << kold << "    N[kold] = " << N[kold] << "\n";
        
        
        
        
        
        
        //cout << "Past the interpolation bit for k = " << k << "\n";
        
        
        
        /*
         
         This is meant to do the interpolation.
         
         It first evaluates whether the point is in the higher-resolution region, if not, it just changes the value the it reads its input from; if yes, it then works out what bit of the higher-resolution region it is from.
         
         First we calculate the value of x it is at, which narrows down the possible values of k.
         
         Then we calculate the value of k it is at.
         
         Then we get the value of n that it is at.
         
         Then we've got all of the interpolation parameters, and we're all set to go.
         
         Note that we'll need to repeat this for the k+1 cell because we involve both k and k+1 in the equations at the same time.
         
         
         
         Equations
         _________
         
         xcount = floor( (kold - 1)/L )
         
         label = ((xcount * L) / 2)*((2*Nmax) + 1 - xcount) + ( kold - (xcount*L) - 1)*(Nmax - xcount) + n + 1
         
         labelminx = minimum value of label for a given value of x
         labelminx = ((x * L) / 2)*((2*Nmax) + 1 - x) + 1
         
         labelmink = minimum value of label for a given value of kold
         labelmink = ((xcount * L) / 2)*((2*Nmax) + 1 - xcount) + ( kold - (xcount*L) - 1)*(Nmax - xcount) + 1
         
         kmin = minimum value of old k for a given value of x
         kmin = (L*xcount) + 1
         
         */
        
        
        
        
        
        
        prefactor = ( flux_tot_HR[0] * conv_L_div_L_HR[0] * ( rmid_cm_HR[1] - rmid_cm_HR[0] ) * ( (delta[0]/K_HR[0]) + (delta[1]/K_HR[1]) ) / ( temperature_HR[1] - temperature_HR[0] ) ) / ( (lnT_HR[1] - lnT_HR[0]) - (delta[0]*gamma1_HR[0] + delta[1]*gamma1_HR[1])*( log(pressure_HR[1]) - log(pressure_HR[0]) ) );
        
        
        // This checks the interpolation of any given variable array
        test[k] = brunt_A_HR[0]*grav_HR[0]/(radius_cm_HR[0]*m*m*omega*omega); //prefactor; //conv_L_div_L_HR[0]; //gamma1_HR[0]; //brunt_A_HR[0]*grav_HR[0]/(radius_cm_HR[0]*m*m*omega*omega); //brunt_A_HR[0]; // sqrt((flux_HR[0] * radius_cm_HR[0]/(pressure_HR[0]*omega*brunt_A_HR[0]))*(grada_HR[0]))/R; //(radius_cm_HR[1] - radius_cm_HR[0])/R; //((temperature_HR[1]-temperature_HR[0])/(rmid_cm_HR[1]-rmid_cm_HR[0]))*( - 4.0*7.5657e-15*2.99792458e10*temperature_HR[0]*temperature_HR[0]*temperature_HR[0]/( 3.0*opacity_HR[0]*rho_HR[0] ) ); //brunt_A_HR[0]*grav_HR[0]/radius_cm_HR[0];
        test[k+1] = brunt_A_HR[1]*grav_HR[1]/(radius_cm_HR[1]*m*m*omega*omega); //prefactor; //conv_L_div_L_HR[1]; //gamma1_HR[1]; //brunt_A_HR[1]*grav_HR[1]/(radius_cm_HR[1]*m*m*omega*omega); //brunt_A_HR[1] ; // sqrt((flux_HR[1] * radius_cm_HR[1]/(pressure_HR[1]*omega*brunt_A_HR[1]))*(grada_HR[1]))/R; //(radius_cm_HR[1] - radius_cm_HR[0])/R; //((temperature_HR[1]-temperature_HR[0])/(rmid_cm_HR[1]-rmid_cm_HR[0]))*( - 4.0*7.5657e-15*2.99792458e10*temperature_HR[0]*temperature_HR[0]*temperature_HR[0]/( 3.0*opacity_HR[0]*rho_HR[0] ) ); //brunt_A_HR[1]*grav_HR[1]/radius_cm_HR[1];
        
        
        
        
        
        
        
        
        flux_HR_output[k] = flux_HR[0];
        flux_HR_output[k+1] = flux_HR[1];
        
        pressure_HR_output[k] = pressure_HR[0];
        pressure_HR_output[k+1] = pressure_HR[1];
        
        temperature_HR_output[k] = temperature_HR[0];
        temperature_HR_output[k+1] = temperature_HR[1];
        
        rho_HR_output[k] = rho_HR[0];
        rho_HR_output[k+1] = rho_HR[1];
        
        
        rho_face_HR_output[k] = rho_face_HR[0];
        rho_face_HR_output[k+1] = rho_face_HR[1];
        
        
        opacity_HR_output[k] = opacity_HR[0];
        opacity_HR_output[k+1] = opacity_HR[1];
        
        
        dkap_dlnT_face_HR_output[k] = dkap_dlnT_face_HR[0];
        dkap_dlnT_face_HR_output[k+1] = dkap_dlnT_face_HR[1];
        
        dkap_dlnrho_face_HR_output[k] = dkap_dlnrho_face_HR[0];
        dkap_dlnrho_face_HR_output[k+1] = dkap_dlnrho_face_HR[1];
        
        
        gamma1_HR_output[k] = gamma1_HR[0];
        gamma1_HR_output[k+1] = gamma1_HR[1];
        
        
        
        K_HR_output[k] = K_HR[0];
        K_HR_output[k+1] = K_HR[1];
        
        
        
        grav_HR_output[k] = grav_HR[0];
        grav_HR_output[k+1] = grav_HR[1];
        
        
        chiRho_HR_output[k] = chiRho_HR[0];
        chiRho_HR_output[k+1] = chiRho_HR[1];
        
        chiT_HR_output[k] = chiT_HR[0];
        chiT_HR_output[k+1] = chiT_HR[1];
        
        brunt_A_HR_output[k] = brunt_A_HR[0];
        brunt_A_HR_output[k+1] = brunt_A_HR[1];
        
        cp_HR_output[k] = cp_HR[0];
        cp_HR_output[k+1] = cp_HR[1];
        
        grada_HR_output[k] = grada_HR[0];
        grada_HR_output[k+1] = grada_HR[1];
        
        conv_L_div_L_HR_output[k] = conv_L_div_L_HR[0];
        conv_L_div_L_HR_output[k+1] = conv_L_div_L_HR[1];
        
        flux_tot_HR_output[k] = flux_tot_HR[0];
        flux_tot_HR_output[k+1] = flux_tot_HR[1];
        
        
        // Made such that F' = dT0/dr * ( Fprime_for_T * T'/T0   +   dT'/dT0   -   dxi/dr   +   Fprime_for_p * p'/p0 )
        // The strange prefactor is therefore: -(4*a*c*T^3)/(3*opacity*rho)
        Fprime_for_T[k] = ( flux_HR[0] ) * ( 3.0 - dkap_dlnT_face_HR[0]/(0.5*(opacity_HR[0]+opacity_HR[1])) + 0.5*((chiT_HR[0]/chiRho_HR[0])+(chiT_HR[1]/chiRho_HR[1]))*( 1.0 + dkap_dlnrho_face_HR[0]/(0.5*(opacity_HR[0]+opacity_HR[1])) ) );
        Fprime_for_T[k+1] = ( flux_HR[1] ) * ( 3.0 - dkap_dlnT_face_HR[1]/opacity_HR[1] + (chiT_HR[1]/chiRho_HR[1])*( 1.0 + dkap_dlnrho_face_HR[1]/opacity_HR[1] ) );
        
        Fprime_for_p[k] = ( flux_HR[0] ) * (-1.0/(0.5*(chiRho_HR[0]+chiRho_HR[1])))*( 1.0 + dkap_dlnrho_face_HR[0]/(0.5*(opacity_HR[0]+opacity_HR[1])) );
        Fprime_for_p[k+1] = ( flux_HR[1] ) * (-1.0/chiRho_HR[1])*( 1.0 + dkap_dlnrho_face_HR[1]/opacity_HR[1] );
        
        // These versions have a self-calculated value for the radiative flux, which has major issues at the surface
        //Fprime_for_T[k] = ( - 4.0*7.5657e-15*2.99792458e10*temperature_HR[0]*temperature_HR[0]*temperature_HR[0]/( 3.0*opacity_HR[0]*rho_HR[0] ) ) * ( 3.0 - dkap_dlnT_face_HR[0]/opacity_HR[0] + (chiT_HR[0]/chiRho_HR[0])*( 1.0 + dkap_dlnrho_face_HR[0]/opacity_HR[0] ) );
        //Fprime_for_T[k+1] = ( - 4.0*7.5657e-15*2.99792458e10*temperature_HR[1]*temperature_HR[1]*temperature_HR[1]/( 3.0*opacity_HR[1]*rho_HR[1] ) ) * ( 3.0 - dkap_dlnT_face_HR[1]/opacity_HR[1] + (chiT_HR[1]/chiRho_HR[1])*( 1.0 + dkap_dlnrho_face_HR[1]/opacity_HR[1] ) );
        
        //Fprime_for_p[k] = ( - 4.0*7.5657e-15*2.99792458e10*temperature_HR[0]*temperature_HR[0]*temperature_HR[0]/( 3.0*opacity_HR[0]*rho_HR[0] ) ) * (-1.0/chiRho_HR[0])*( 1.0 + dkap_dlnrho_face_HR[0]/opacity_HR[0] );
        //Fprime_for_p[k+1] = ( - 4.0*7.5657e-15*2.99792458e10*temperature_HR[1]*temperature_HR[1]*temperature_HR[1]/( 3.0*opacity_HR[1]*rho_HR[1] ) ) * (-1.0/chiRho_HR[1])*( 1.0 + dkap_dlnrho_face_HR[1]/opacity_HR[1] );
        
        
        
        
        if (k==J-1) {
            // If we are at the very, very surface, we just use the previous value because it will be very probably pretty similar
            eta_Terquem_output[k] = eta_Terquem_output[k-1];
            
        } else {
            
            eta_Terquem_output[k] = (0.5*(brunt_A_HR[0] + brunt_A_HR[1])) * ( ( pressure_HR[1] + pressure_HR[0] ) / ( pressure_HR[1] - pressure_HR[0] ) ) * ( ( rmid_cm_HR[1] - rmid_cm_HR[0] ) / ( rmid_cm_HR[1] + rmid_cm_HR[0] ) );
            
        }
        
        
        
        
        // This defines the delta coefficients for interpolation
        delta[0] = (rmid_cm_HR[1] - radius_cm_HR[0])/(rmid_cm_HR[1] - rmid_cm_HR[0]);
        delta[1] = (radius_cm_HR[0] - rmid_cm_HR[0])/(rmid_cm_HR[1] - rmid_cm_HR[0]);
        
        /*
         
         This is the case with the equations used for the actual oscillations.
         
         
         (1/(rho * r * r))*(d/dr)*( r * r * r * rho * a ) + ( (1/chiRho) - (p/rho) * ( l(l+1)/(r*r*omega*omega) ) * c - ( chiT / chiRho ) * d  ===>  3*a + r * ( d(lnRho)/dr ) * a + r * d(a)/dr + ( (1/chiRho) - (p/rho) * ( l(l+1)/(r*r*omega*omega) ) * c - ( chiT / chiRho ) * d
         
         
         i * A * ( chiRho / chiT ) * a + (1.0/( omega * rho * T * cp * r * r ))*( 2*r*F*b + r*r*(d(F)/dr)*b + r*r*F*(d(b)/dr) ) - i*grada*c + ( i + ( (l*(l+1.0)/(r*r)) * ( K/(omega*rho*cp) ) ) )*d
         // NB - A here is defined in MESA as A = N*N*r/g, whereas Terquem defines it to be - N*N/g   =====>>>  The difference is a factor of -r
         
         
         (d(r)/dT)*(F/K)*b - (1.0/chiRho)*( 1.0 + (dkap_dlnRho / kappa) )*c + ( 4.0 - dkap_dlnT / kappa + (chiT/chiRho)*(1.0 + dkap_dlnRho / kappa) )*d + d(d)/dlnT
         
         
         ( m*m*omega*omega*r/g ) * a - ( 1.0/chiRho )*c - (1.0/(g * rho))*d( p*c )/dr + ( chiT/chiRho )*d
         
         
         
         
         BC:
         a + c
         a + b + d
         
         
         This gives the case where everything = x
         
         Mr[0][0][0] = 3.0 + 2.0*(rmid_cm_HR[1]/R);
         Mr[0][1][0] = 3.0 + (rmid_cm_HR[1]/R);
         
         Mi[0][0][0] = 0;
         Mi[0][1][0] = 0;
         
         
         Nr[0][0][0] = 1.0 + 3.0*(radius_cm_HR[0]/R);
         Nr[0][1][0] = 1.0 + 3.0*(radius_cm_HR[0]/R);
         
         Ni[0][0][0] = 0;
         Ni[0][1][0] = 0;
         
         
         xr[0][0][0] = 2.0*(rmid_cm_HR[1]/R);
         xr[0][1][0] = 3.0*(rmid_cm_HR[1]/R);
         
         xi[0][0][0] = 0;
         xi[0][1][0] = 0;
         
         
         
         
         This gives the case where:
         a = x^2
         b = sin(x)
         c = cos(3.6*x)
         d = cosh(x)
         
         Mr[0][0][0] = 4.0*(rmid_cm_HR[1]/R) + cos(3.6*rmid_cm_HR[1]/R) + cosh(rmid_cm_HR[1]/R);
         Mr[0][1][0] = (rmid_cm_HR[1]/R)*(rmid_cm_HR[1]/R) + 2.0*(rmid_cm_HR[1]/R)/(rmid_cm_HR[1]/R) + cos(rmid_cm_HR[1]/R) + cos(3.6*rmid_cm_HR[1]/R) - cosh(rmid_cm_HR[1]/R);
         
         Mi[0][0][0] = 0;
         Mi[0][1][0] = 0;
         
         Nr[0][0][0] = sin(radius_cm_HR[0]/R) + cos(3.6*radius_cm_HR[0]/R) + cosh(radius_cm_HR[0]/R) + sinh(radius_cm_HR[0]/R);
         Nr[0][1][0] = (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R) + cos(3.6*radius_cm_HR[0]/R) - 3.6*sin(3.6*radius_cm_HR[0]/R) + cosh(radius_cm_HR[0]/R);
         
         Ni[0][0][0] = 0;
         Ni[0][1][0] = 0;
         
         xr[0][0][0] = (rmid_cm_HR[1]/R)*(rmid_cm_HR[1]/R) + cos(3.6*rmid_cm_HR[1]/R);
         xr[0][1][0] = (rmid_cm_HR[1]/R)*(rmid_cm_HR[1]/R) + sin(rmid_cm_HR[1]/R) + cosh(rmid_cm_HR[1]/R);
         
         xi[0][0][0] = 0;
         xi[0][1][0] = 0;
         
         
         
         
         
         
         
         This gives the case where:
         a = (r/R)*sin(r/R)
         b = sqrt(r/R)*sin(sqrt(200.0)*r/R)
         c = 50.0*sin(sqrt(17.0)*r/R)/(sqrt(17.0)*r/R)
         d = cosh(r/R) + 3.0*cos(sqrt(23.0)*r/R)
         
         NB: in order to avoid the code producing a pole in c', we use the expanded form of c' as: (-(850 x)/3 + (1445 x^3)/3 - (24565 x^5)/84 + (417605 x^7)/4536 - (7099285 x^9)/399168 + (24137569 x^11)/10378368)
         This approximation does break down a bit as you approach 1, but not as badly as when you get a pole at x=0
         
         c' approx = (-(850 * (radius_cm_HR[0]/R))/3 + (1445 * (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R))/3 - (24565 * (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R))/84 + (417605 * (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R))/4536 - (7099285 * (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R))/399168 + (24137569 * (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R))/10378368)
         
         M00 = sin(x) + x*cos(x) + sqrt(12800*x)*cos(sqrt(200)*x) + (4/sqrt(x))*sin(sqrt(200)*x) + 50*sin(sqrt(17)*x)/(sqrt(17)*x) + cosh(x) + 3*cos(sqrt(23)*x)
         M10 = sin(x) + x*cos(x) + sqrt(200*x)*cos(sqrt(200)*x) + (0.5/sqrt(x))*sin(sqrt(200)*x) + 150*sin(sqrt(17)*x)/(sqrt(17)*x) - 2*cosh(x) - 6*cos(sqrt(23)*x)
         
         N00 = x*sin(x) + 4*sqrt(x)*sin(sqrt(200)*x) + 6*sinh(x) - 6*sqrt(207)*sin(sqrt(23)*x) - (-(850 x)/3 + (1445 x^3)/3 - (24565 x^5)/84 + (417605 x^7)/4536 - (7099285 x^9)/399168 + (24137569 x^11)/10378368)
         N10 = x*sin(x) + 7*sqrt(x)*sin(sqrt(200)*x) - sinh(x) + sqrt(207)*sin(sqrt(23)*x) - 5*(-(850*x)/3 + (1445*x**3)/3 - (24565*x**5)/84 + (417605*x**7)/4536 - (7099285*x**9)/399168 + (24137569*x**11)/10378368)
         
         Mr[0][0][0] = sin(rmid_cm_HR[1]/R) + (rmid_cm_HR[1]/R)*cos(rmid_cm_HR[1]/R) + sqrt(12800*rmid_cm_HR[1]/R)*cos(sqrt(200)*rmid_cm_HR[1]/R) + (4.0/sqrt(rmid_cm_HR[1]/R))*sin(sqrt(200)*rmid_cm_HR[1]/R) + 50.0*sin(sqrt(17.0)*rmid_cm_HR[1]/R)/(sqrt(17.0)*rmid_cm_HR[1]/R) + cosh(rmid_cm_HR[1]/R) + 3.0*cos(sqrt(23.0)*rmid_cm_HR[1]/R);
         Mr[0][1][0] = sin(rmid_cm_HR[1]/R) + (rmid_cm_HR[1]/R)*cos(rmid_cm_HR[1]/R) + sqrt(200*rmid_cm_HR[1]/R)*cos(sqrt(200)*rmid_cm_HR[1]/R) + (0.5/sqrt(rmid_cm_HR[1]/R))*sin(sqrt(200)*rmid_cm_HR[1]/R) + 3.0*50.0*sin(sqrt(17.0)*rmid_cm_HR[1]/R)/(sqrt(17.0)*rmid_cm_HR[1]/R) - 2.0*cosh(rmid_cm_HR[1]/R) - 2.0*3.0*cos(sqrt(23.0)*rmid_cm_HR[1]/R);
         
         Mi[0][0][0] = 0;
         Mi[0][1][0] = 0;
         
         if (radius_cm_HR[0]/R < 0.3) {
         
         Nr[0][0][0] = (radius_cm_HR[0]/R)*sin(radius_cm_HR[0]/R) + 4.0*sqrt(radius_cm_HR[0]/R)*sin(sqrt(200)*radius_cm_HR[0]/R) + 6.0*sinh(radius_cm_HR[0]/R) - 6.0*sqrt(207.0)*sin(sqrt(23.0)*radius_cm_HR[0]/R) - (-(850 * (radius_cm_HR[0]/R))/3 + (1445 * (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R))/3 - (24565 * (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R))/84 + (417605 * (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R))/4536 - (7099285 * (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R))/399168 + (24137569 * (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R))/10378368);
         Nr[0][1][0] = (radius_cm_HR[0]/R)*sin(radius_cm_HR[0]/R) + 7.0*sqrt(radius_cm_HR[0]/R)*sin(sqrt(200)*radius_cm_HR[0]/R) - sinh(radius_cm_HR[0]/R) + sqrt(207.0)*sin(sqrt(23.0)*radius_cm_HR[0]/R) - 5.0*(-(850 * (radius_cm_HR[0]/R))/3 + (1445 * (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R))/3 - (24565 * (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R))/84 + (417605 * (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R))/4536 - (7099285 * (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R))/399168 + (24137569 * (radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)*(radius_cm_HR[0]/R))/10378368);
         
         } else {
         
         Nr[0][0][0] = (radius_cm_HR[0]/R)*sin(radius_cm_HR[0]/R) + 4.0*sqrt(radius_cm_HR[0]/R)*sin(sqrt(200)*radius_cm_HR[0]/R) + 6.0*sinh(radius_cm_HR[0]/R) - 6.0*sqrt(207.0)*sin(sqrt(23.0)*radius_cm_HR[0]/R) - (50.0/sqrt(17))*( ((sqrt(17)*(radius_cm_HR[0]/R)*cos(sqrt(17)*radius_cm_HR[0]/R) - sin(sqrt(17)*radius_cm_HR[0]/R)))/((radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)) );
         Nr[0][1][0] = (radius_cm_HR[0]/R)*sin(radius_cm_HR[0]/R) + 7.0*sqrt(radius_cm_HR[0]/R)*sin(sqrt(200)*radius_cm_HR[0]/R) - sinh(radius_cm_HR[0]/R) + sqrt(207.0)*sin(sqrt(23.0)*radius_cm_HR[0]/R) - 5.0*(50.0/sqrt(17))*( ((sqrt(17)*(radius_cm_HR[0]/R)*cos(sqrt(17)*radius_cm_HR[0]/R) - sin(sqrt(17)*radius_cm_HR[0]/R)))/((radius_cm_HR[0]/R)*(radius_cm_HR[0]/R)) );
         
         
         }
         
         
         Ni[0][0][0] = 0;
         Ni[0][1][0] = 0;
         
         xr[0][0][0] = (rmid_cm_HR[1]/R)*sin(rmid_cm_HR[1]/R) + 50.0*sin(sqrt(17.0)*rmid_cm_HR[1]/R)/(sqrt(17.0)*rmid_cm_HR[1]/R);
         xr[0][1][0] = (rmid_cm_HR[1]/R)*sin(rmid_cm_HR[1]/R) + sqrt(rmid_cm_HR[1]/R)*sin(sqrt(200)*rmid_cm_HR[1]/R) + cosh(rmid_cm_HR[1]/R) + 3.0*cos(sqrt(23.0)*rmid_cm_HR[1]/R);
         
         xi[0][0][0] = 0;
         xi[0][1][0] = 0;
         
         
         
         */
        
        //        cout << 1/(radius_cm[k+1] - radius_cm[k]) << "\n";
        
        
        // Option 1 (be aware of the need to modify for rho_face in Ar[0][0][0] and Cr[0][0][0])
        /*
         Ar[0][0][0] = - radius_cm_HR[0]*radius_cm_HR[0]*rho_HR[1] / ( R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) );
         Ar[0][0][1] = 0;
         Ar[0][1][0] = 0;
         Ar[0][1][1] = -radius_cm_HR[0]*radius_cm_HR[0]/( R*R*(radius_cm_HR[1] - radius_cm_HR[0]) );
         
         Ai[0][0][0] = 0;
         Ai[0][0][1] = 0;
         Ai[0][1][0] = 0; // m*omega*rho_HR[1]*temperature_HR[1]*brunt_A_HR[1]*chiRho_HR[1]*cp_HR[1]*radius_cm_HR[0]*0.5 / (R*R*chiT_HR[1]);
         Ai[0][1][1] = 0;
         
         Cr[0][0][0] = radius_cm_HR[1]*radius_cm_HR[1]*rho_HR[1] / ( R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) );
         Cr[0][0][1] = 0;
         Cr[0][1][0] = 0;
         Cr[0][1][1] = radius_cm_HR[1]*radius_cm_HR[1]/( R*R*(radius_cm_HR[1] - radius_cm_HR[0]) );
         
         Ci[0][0][0] = 0;
         Ci[0][0][1] = 0;
         Ci[0][1][0] = 0; // m*omega*rho_HR[1]*temperature_HR[1]*brunt_A_HR[1]*chiRho_HR[1]*cp_HR[1]*radius_cm_HR[1]*0.5 / (R*R*chiT_HR[1]);
         Ci[0][1][1] = 0;
         
         Dr[0][0][0] = ( rho_HR[1]*rmid_cm_HR[1]*rmid_cm_HR[1] / (chiRho_HR[1] * pressure_HR[1]*R*R) )  -  ( ( l*(l + 1.0) )/( m*m*omega*omega*R*R ) );
         Dr[0][0][1] = - rho_HR[1] * rmid_cm_HR[1] * rmid_cm_HR[1] * chiT_HR[1] / ( temperature_HR[1] * R * R * chiRho_HR[1] );
         Dr[0][1][0] = 0;
         Dr[0][1][1] = l*(l + 1.0)*K_HR[1] / (R*R);
         
         Di[0][0][0] = 0;
         Di[0][0][1] = 0;
         Di[0][1][0] = 0; // - m*omega*rho_HR[1]*temperature_HR[1]*cp_HR[1]*grada_HR[1]*rmid_cm_HR[1]*rmid_cm_HR[1] / ( R*R*pressure_HR[1] );
         Di[0][1][1] = 0; // m*omega*rho_HR[1]*cp_HR[1]*rmid_cm_HR[1]*rmid_cm_HR[1];
         
         Er[0][0][0] = 0;
         Er[0][0][1] = - ( (delta[0]/K_HR[0]) + (delta[1]/K_HR[1]) );
         Er[0][1][0] = -m*m*omega*omega*( (delta[0]*rho_HR[0]) + (delta[1]*rho_HR[1]) );
         Er[0][1][1] = 0;
         
         Ei[0][0][0] = 0;
         Ei[0][0][1] = 0;
         Ei[0][1][0] = 0;
         Ei[0][1][1] = 0;
         
         Fr[0][0][0] = (( temperature_HR[1] - temperature_HR[0] )/( rmid_cm_HR[1] - rmid_cm_HR[0] ))*( delta[0]*( 1.0 + (dkap_dlnrho_face_HR[0]/opacity_HR[0]) )/( pressure_HR[0]*chiRho_HR[0] ) );
         Fr[0][0][1] = (1.0/( rmid_cm_HR[1] - rmid_cm_HR[0] ))  +  ( (lnT_HR[1] - lnT_HR[0])/(rmid_cm_HR[1] - rmid_cm_HR[0]) )*delta[0]*( -3.0 + (dkap_dlnT_face_HR[0]/opacity_HR[0])  -  (chiT_HR[0]/chiRho_HR[0])*( 1.0 + (dkap_dlnrho_face_HR[0]/opacity_HR[0]) ) );
         Fr[0][1][0] = (-1.0/(rmid_cm_HR[1] - rmid_cm_HR[0]))  +  (grav_HR[0]*rho_HR[0]*delta[0])/(chiRho_HR[0]*pressure_HR[0]);
         Fr[0][1][1] = -delta[0]*( chiT_HR[0]/chiRho_HR[0] )*(grav_HR[0]*rho_HR[0]/temperature_HR[0]);
         
         Fi[0][0][0] = 0;
         Fi[0][0][1] = 0;
         Fi[0][1][0] = 0;
         Fi[0][1][1] = 0;
         
         Hr[0][0][0] = (( temperature_HR[1] - temperature_HR[0] )/( rmid_cm_HR[1] - rmid_cm_HR[0] ))*( delta[1]*( 1.0 + (dkap_dlnrho_face_HR[0]/opacity_HR[1]) )/( pressure_HR[1]*chiRho_HR[1] ) );
         Hr[0][0][1] = (-1.0/( rmid_cm_HR[1] - rmid_cm_HR[0] ))  +  ( (lnT_HR[1] - lnT_HR[0])/(rmid_cm_HR[1] - rmid_cm_HR[0]) )*delta[1]*( -3.0 + (dkap_dlnT_face_HR[0]/opacity_HR[1])  -  (chiT_HR[1]/chiRho_HR[1])*( 1.0 + (dkap_dlnrho_face_HR[0]/opacity_HR[1]) ) );
         Hr[0][1][0] = (1.0/(rmid_cm_HR[1] - rmid_cm_HR[0]))  +  (grav_HR[1]*rho_HR[1]*delta[1])/(chiRho_HR[1]*pressure_HR[1]);
         Hr[0][1][1] = -delta[1]*( chiT_HR[1]/chiRho_HR[1] )*(grav_HR[1]*rho_HR[1]/temperature_HR[1]);
         
         Hi[0][0][0] = 0;
         Hi[0][0][1] = 0;
         Hi[0][1][0] = 0;
         Hi[0][1][1] = 0;
         */
        
        
        // Option 2  --  NB: I don't think this is currently correct
        /*
         Ar[0][0][0] = - radius_cm_HR[0]*radius_cm_HR[0]*rho_HR[0] / ( R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) );
         Ar[0][0][1] = 0;
         Ar[0][1][0] = 0;
         Ar[0][1][1] = -radius_cm_HR[0]*radius_cm_HR[0]/( R*R*(radius_cm_HR[1] - radius_cm_HR[0]) );
         
         Ai[0][0][0] = 0;
         Ai[0][0][1] = 0;
         Ai[0][1][0] = 0; // m*omega*rho_HR[1]*temperature_HR[1]*brunt_A_HR[1]*chiRho_HR[1]*cp_HR[1]*radius_cm_HR[0]*0.5 / (R*R*chiT_HR[1]);
         Ai[0][1][1] = 0;
         
         Cr[0][0][0] = radius_cm_HR[1]*radius_cm_HR[1]*rho_HR[1] / ( R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) );
         Cr[0][0][1] = 0;
         Cr[0][1][0] = 0;
         Cr[0][1][1] = radius_cm_HR[1]*radius_cm_HR[1]/( R*R*(radius_cm_HR[1] - radius_cm_HR[0]) );
         
         Ci[0][0][0] = 0;
         Ci[0][0][1] = 0;
         Ci[0][1][0] = 0; // m*omega*rho_HR[1]*temperature_HR[1]*brunt_A_HR[1]*chiRho_HR[1]*cp_HR[1]*radius_cm_HR[1]*0.5 / (R*R*chiT_HR[1]);
         Ci[0][1][1] = 0;
         
         Dr[0][0][0] = ( rho_HR[1]*rmid_cm_HR[1]*rmid_cm_HR[1] / (chiRho_HR[1] * pressure_HR[1]*R*R) )  -  ( ( l*(l + 1.0) )/( m*m*omega*omega*R*R ) );
         Dr[0][0][1] = - rho_HR[1] * rmid_cm_HR[1] * rmid_cm_HR[1] * chiT_HR[1] / ( temperature_HR[1] * R * R * chiRho_HR[1] );
         Dr[0][1][0] = 0;
         Dr[0][1][1] = l*(l + 1.0)*K_HR[1] / (R*R);
         
         Di[0][0][0] = 0;
         Di[0][0][1] = 0;
         Di[0][1][0] = 0; // - m*omega*rho_HR[1]*temperature_HR[1]*cp_HR[1]*grada_HR[1]*rmid_cm_HR[1]*rmid_cm_HR[1] / ( R*R*pressure_HR[1] );
         Di[0][1][1] = 0; // m*omega*rho_HR[1]*cp_HR[1]*rmid_cm_HR[1]*rmid_cm_HR[1]/(R*R);
         
         Er[0][0][0] = 0;
         Er[0][0][1] = - ( (delta[0]/K_HR[0]) + (delta[1]/K_HR[1]) );
         Er[0][1][0] = -m*m*omega*omega*( (delta[0]*rho_HR[0]) + (delta[1]*rho_HR[1]) );
         Er[0][1][1] = 0;
         
         Ei[0][0][0] = 0;
         Ei[0][0][1] = 0;
         Ei[0][1][0] = 0;
         Ei[0][1][1] = 0;
         
         Fr[0][0][0] = (( temperature_HR[1] - temperature_HR[0] )/( rmid_cm_HR[1] - rmid_cm_HR[0] ))*( delta[0]*( 1.0 + (dkap_dlnrho_face_HR[0]/opacity_HR[0]) )/( pressure_HR[0]*chiRho_HR[0] ) );
         Fr[0][0][1] = (1.0/( rmid_cm_HR[1] - rmid_cm_HR[0] ))  +  ( (lnT_HR[1] - lnT_HR[0])/(rmid_cm_HR[1] - rmid_cm_HR[0]) )*delta[0]*( -3.0 + (dkap_dlnT_face_HR[0]/opacity_HR[0])  -  (chiT_HR[0]/chiRho_HR[0])*( 1.0 + (dkap_dlnrho_face_HR[0]/opacity_HR[0]) ) );
         Fr[0][1][0] = (-1.0/(rmid_cm_HR[1] - rmid_cm_HR[0]))  +  (grav_HR[0]*rho_HR[0]*delta[0])/(chiRho_HR[0]*pressure_HR[0]);
         Fr[0][1][1] = -delta[0]*( chiT_HR[0]/chiRho_HR[0] )*(grav_HR[0]*rho_HR[0]/temperature_HR[0]);
         
         Fi[0][0][0] = 0;
         Fi[0][0][1] = 0;
         Fi[0][1][0] = 0;
         Fi[0][1][1] = 0;
         
         Hr[0][0][0] = (( temperature_HR[1] - temperature_HR[0] )/( rmid_cm_HR[1] - rmid_cm_HR[0] ))*( delta[1]*( 1.0 + (dkap_dlnrho_face_HR[0]/opacity_HR[1]) )/( pressure_HR[1]*chiRho_HR[1] ) );
         Hr[0][0][1] = (-1.0/( rmid_cm_HR[1] - rmid_cm_HR[0] ))  +  ( (lnT_HR[1] - lnT_HR[0])/(rmid_cm_HR[1] - rmid_cm_HR[0]) )*delta[1]*( -3.0 + (dkap_dlnT_face_HR[0]/opacity_HR[1])  -  (chiT_HR[1]/chiRho_HR[1])*( 1.0 + (dkap_dlnrho_face_HR[0]/opacity_HR[1]) ) );
         Hr[0][1][0] = (1.0/(rmid_cm_HR[1] - rmid_cm_HR[0]))  +  (grav_HR[1]*rho_HR[1]*delta[1])/(chiRho_HR[1]*pressure_HR[1]);
         Hr[0][1][1] = -delta[1]*( chiT_HR[1]/chiRho_HR[1] )*(grav_HR[1]*rho_HR[1]/temperature_HR[1]);
         
         Hi[0][0][0] = 0;
         Hi[0][0][1] = 0;
         Hi[0][1][0] = 0;
         Hi[0][1][1] = 0;
         */
        
        
        
        // Option 2.i
        /*
         Ar[0][0][0] = - radius_cm_HR[0]*radius_cm_HR[0]*rho_face_HR[0] / ( rho_HR[1]*R*( radius_cm_HR[1] - radius_cm_HR[0] ) ); // - radius_cm_HR[0]*radius_cm_HR[0]*radius_cm_HR[0]*rho_face_HR[0] / ( rho_HR[1]*R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) );
         Ar[0][0][1] = 0;
         Ar[0][1][0] = 0;
         Ar[0][1][1] = - ( radius_cm_HR[0]*radius_cm_HR[0]*flux_BC )/( m*omega*rho_HR[1]*temperature_HR[1]*cp_HR[1]*R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) ); // - ( radius_cm_HR[0]*radius_cm_HR[0]*flux_HR[0] )/( m*omega*rho_HR[1]*temperature_HR[1]*cp_HR[1]*R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) );
         
         Ai[0][0][0] = 0.0;
         Ai[0][0][1] = 0.0;
         Ai[0][1][0] = 0.0; // ( brunt_A_HR[1]/( R ) ) * ( chiRho_HR[1] / chiT_HR[1] ) * 0.5 * radius_cm_HR[0]; // ( brunt_A_HR[1]/( R*R ) ) * ( chiRho_HR[1] / chiT_HR[1] ) * 0.5 * radius_cm_HR[0] * radius_cm_HR[0];
         Ai[0][1][1] = 0.0;
         
         Cr[0][0][0] = radius_cm_HR[1]*radius_cm_HR[1]*rho_face_HR[1] / ( rho_HR[1]*R*( radius_cm_HR[1] - radius_cm_HR[0] ) ); // radius_cm_HR[1]*radius_cm_HR[1]*radius_cm_HR[1]*rho_face_HR[1] / ( rho_HR[1]*R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) );
         Cr[0][0][1] = 0.0;
         Cr[0][1][0] = 0.0;
         Cr[0][1][1] = ( radius_cm_HR[1]*radius_cm_HR[1]*flux_BC )/( m*omega*rho_HR[1]*temperature_HR[1]*cp_HR[1]*R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) ); // ( radius_cm_HR[1]*radius_cm_HR[1]*flux_HR[1] )/( m*omega*rho_HR[1]*temperature_HR[1]*cp_HR[1]*R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) );
         
         Ci[0][0][0] = 0.0;
         Ci[0][0][1] = 0.0;
         Ci[0][1][0] = 0.0; // ( brunt_A_HR[1]/( R ) ) * ( chiRho_HR[1] / chiT_HR[1] ) * 0.5 * radius_cm_HR[1]; // ( brunt_A_HR[1]/( R*R ) ) * ( chiRho_HR[1] / chiT_HR[1] ) * 0.5 * radius_cm_HR[1] * radius_cm_HR[1];
         Ci[0][1][1] = 0.0;
         
         Dr[0][0][0] = ( rmid_cm_HR[1]*rmid_cm_HR[1] / ( chiRho_HR[1]*R*R ) )  -  ( (l*(l+1.0)*pressure_HR[1])/(m*m*omega*omega*R*R*rho_HR[1]) );
         Dr[0][0][1] = - (chiT_HR[1] / chiRho_HR[1]) * ( rmid_cm_HR[1]*rmid_cm_HR[1] / (R*R) );
         Dr[0][1][0] = 0.0;
         Dr[0][1][1] = ( l * ( l + 1.0 ) * K_HR[1] )/( m * omega * rho_HR[1] * cp_HR[1] * R * R );
         
         Di[0][0][0] = 0.0;
         Di[0][0][1] = 0.0;
         Di[0][1][0] = 0.0; // - grada_HR[1] * rmid_cm_HR[1] * rmid_cm_HR[1] / ( R*R );
         Di[0][1][1] = 0.0; // rmid_cm_HR[1] * rmid_cm_HR[1] / ( R * R );
         
         Er[0][0][0] = 0.0;
         Er[0][0][1] = - flux_BC * ( rmid_cm_HR[1] - rmid_cm_HR[0] ) * ( (delta[0]/K_HR[0]) + (delta[1]/K_HR[1]) ) / ( temperature_HR[1] - temperature_HR[0] ); // - flux_HR[0] * ( rmid_cm_HR[1] - rmid_cm_HR[0] ) * ( (delta[0]/K_HR[0]) + (delta[1]/K_HR[1]) ) / ( temperature_HR[1] - temperature_HR[0] );
         Er[0][1][0] = - 1.0; // - radius_cm_HR[0] / R;
         Er[0][1][1] = 0.0;
         
         Ei[0][0][0] = 0.0;
         Ei[0][0][1] = 0.0;
         Ei[0][1][0] = 0.0;
         Ei[0][1][1] = 0.0;
         
         Fr[0][0][0] = ( delta[0] / chiRho_HR[0] ) * ( 1.0 + ( dkap_dlnrho_face_HR[0] / opacity_HR[0] ) );
         Fr[0][0][1] = delta[0] * ( -4.0 + ( dkap_dlnT_face_HR[0] / opacity_HR[0] ) - ( chiT_HR[0]/chiRho_HR[0] )*( 1.0 + ( dkap_dlnrho_face_HR[0] / opacity_HR[0] ) ) )   +   ( (delta[0]*temperature_HR[0]) + (delta[1]*temperature_HR[1]) )/( temperature_HR[1] - temperature_HR[0] );
         Fr[0][1][0] = (1.0/( m*m*omega*omega*R )) * (   (- pressure_HR[0] * ( (delta[0]/rho_HR[0]) + (delta[1]/rho_HR[1]) ))/( rmid_cm_HR[1] - rmid_cm_HR[0] )   +   ( delta[0] * grav_HR[0] / chiRho_HR[0] )   );
         Fr[0][1][1] = - ( delta[0] * chiT_HR[0] * grav_HR[0] )/( m*m*omega*omega*R*chiRho_HR[0] );
         
         Fi[0][0][0] = 0.0;
         Fi[0][0][1] = 0.0;
         Fi[0][1][0] = 0.0;
         Fi[0][1][1] = 0.0;
         
         Hr[0][0][0] = ( delta[1] / chiRho_HR[1] ) * ( 1.0 + ( dkap_dlnrho_face_HR[0] / opacity_HR[1] ) );
         Hr[0][0][1] = delta[1] * ( -4.0 + ( dkap_dlnT_face_HR[0] / opacity_HR[1] ) - ( chiT_HR[1]/chiRho_HR[1] )*( 1.0 + ( dkap_dlnrho_face_HR[0] / opacity_HR[1] ) ) )   -   ( (delta[0]*temperature_HR[0]) + (delta[1]*temperature_HR[1]) )/( temperature_HR[1] - temperature_HR[0] );
         Hr[0][1][0] = (1.0/( m*m*omega*omega*R )) * (   ( pressure_HR[1] * ( (delta[0]/rho_HR[0]) + (delta[1]/rho_HR[1]) ))/( rmid_cm_HR[1] - rmid_cm_HR[0] )   +   ( delta[1] * grav_HR[1] / chiRho_HR[1] )   );
         Hr[0][1][1] = - ( delta[1] * chiT_HR[1] * grav_HR[1] )/( m*m*omega*omega*R*chiRho_HR[1] );
         
         Hi[0][0][0] = 0.0;
         Hi[0][0][1] = 0.0;
         Hi[0][1][0] = 0.0;
         Hi[0][1][1] = 0.0;
         */
        
        
        
        
        
        Ar[0][0][0] = - radius_cm_HR[0]*radius_cm_HR[0]*rho_face_HR[0] / ( rho_HR[1]*R*( radius_cm_HR[1] - radius_cm_HR[0] ) ); // - radius_cm_HR[0]*radius_cm_HR[0]*radius_cm_HR[0]*rho_face_HR[0] / ( rho_HR[1]*R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) );
        Ar[0][0][1] = 0;
        Ar[0][1][0] = 0;
        Ar[0][1][1] = - ( radius_cm_HR[0]*radius_cm_HR[0]*flux_BC )/( m*omega*rho_HR[1]*temperature_HR[1]*cp_HR[1]*R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) ); // - ( radius_cm_HR[0]*radius_cm_HR[0]*flux_HR[0] )/( m*omega*rho_HR[1]*temperature_HR[1]*cp_HR[1]*R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) );
        
        Ai[0][0][0] = 0.0;
        Ai[0][0][1] = 0.0;
        Ai[0][1][0] = prop*( brunt_A_HR[1]/( R ) ) * ( chiRho_HR[1] / chiT_HR[1] ) * 0.5 * radius_cm_HR[0]; // ( brunt_A_HR[1]/( R*R ) ) * ( chiRho_HR[1] / chiT_HR[1] ) * 0.5 * radius_cm_HR[0] * radius_cm_HR[0];
        Ai[0][1][1] = 0.0;
        
        Cr[0][0][0] = radius_cm_HR[1]*radius_cm_HR[1]*rho_face_HR[1] / ( rho_HR[1]*R*( radius_cm_HR[1] - radius_cm_HR[0] ) ); // radius_cm_HR[1]*radius_cm_HR[1]*radius_cm_HR[1]*rho_face_HR[1] / ( rho_HR[1]*R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) );
        Cr[0][0][1] = 0.0;
        Cr[0][1][0] = 0.0;
        Cr[0][1][1] = ( radius_cm_HR[1]*radius_cm_HR[1]*flux_BC )/( m*omega*rho_HR[1]*temperature_HR[1]*cp_HR[1]*R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) ); // ( radius_cm_HR[1]*radius_cm_HR[1]*flux_HR[1] )/( m*omega*rho_HR[1]*temperature_HR[1]*cp_HR[1]*R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) );
        
        Ci[0][0][0] = 0.0;
        Ci[0][0][1] = 0.0;
        Ci[0][1][0] = prop*( brunt_A_HR[1]/( R ) ) * ( chiRho_HR[1] / chiT_HR[1] ) * 0.5 * radius_cm_HR[1]; // ( brunt_A_HR[1]/( R*R ) ) * ( chiRho_HR[1] / chiT_HR[1] ) * 0.5 * radius_cm_HR[1] * radius_cm_HR[1];
        Ci[0][1][1] = 0.0;
        
        Dr[0][0][0] = ( rmid_cm_HR[1]*rmid_cm_HR[1] / ( chiRho_HR[1]*R*R ) )  -  ( (l*(l+1.0)*pressure_HR[1])/(m*m*omega*omega*R*R*rho_HR[1]) );
        Dr[0][0][1] = - (chiT_HR[1] / chiRho_HR[1]) * ( rmid_cm_HR[1]*rmid_cm_HR[1] / (R*R) );
        Dr[0][1][0] = 0.0;
        Dr[0][1][1] = ( l * ( l + 1.0 ) * K_HR[1] )/( m * omega * rho_HR[1] * cp_HR[1] * R * R );
        
        Di[0][0][0] = 0.0;
        Di[0][0][1] = 0.0;
        Di[0][1][0] = prop*(- grada_HR[1]) * rmid_cm_HR[1] * rmid_cm_HR[1] / ( R*R );
        Di[0][1][1] = prop*rmid_cm_HR[1] * rmid_cm_HR[1] / ( R * R );
        
        Er[0][0][0] = 0.0;
        Er[0][0][1] = - flux_BC * ( rmid_cm_HR[1] - rmid_cm_HR[0] ) * ( (delta[0]/K_HR[0]) + (delta[1]/K_HR[1]) ) / ( temperature_HR[1] - temperature_HR[0] ); // - flux_HR[0] * ( rmid_cm_HR[1] - rmid_cm_HR[0] ) * ( (delta[0]/K_HR[0]) + (delta[1]/K_HR[1]) ) / ( temperature_HR[1] - temperature_HR[0] );
        Er[0][1][0] = - 1.0; // - radius_cm_HR[0] / R;
        Er[0][1][1] = 0.0;
        
        Ei[0][0][0] = 0.0;
        Ei[0][0][1] = 0.0;
        Ei[0][1][0] = 0.0;
        Ei[0][1][1] = 0.0;
        
        Fr[0][0][0] = ( delta[0] / chiRho_HR[0] ) * ( 1.0 + ( dkap_dlnrho_face_HR[0] / opacity_HR[0] ) );
        Fr[0][0][1] = delta[0] * ( -4.0 + ( dkap_dlnT_face_HR[0] / opacity_HR[0] ) - ( chiT_HR[0]/chiRho_HR[0] )*( 1.0 + ( dkap_dlnrho_face_HR[0] / opacity_HR[0] ) ) )   +   ( (delta[0]*temperature_HR[0]) + (delta[1]*temperature_HR[1]) )/( temperature_HR[1] - temperature_HR[0] );
        Fr[0][1][0] = (1.0/( m*m*omega*omega*R )) * (   (- pressure_HR[0] * ( (delta[0]/rho_HR[0]) + (delta[1]/rho_HR[1]) ))/( rmid_cm_HR[1] - rmid_cm_HR[0] )   +   ( delta[0] * grav_HR[0] / chiRho_HR[0] )   );
        Fr[0][1][1] = - ( delta[0] * chiT_HR[0] * grav_HR[0] )/( m*m*omega*omega*R*chiRho_HR[0] );
        
        Fi[0][0][0] = 0.0;
        Fi[0][0][1] = 0.0;
        Fi[0][1][0] = 0.0;
        Fi[0][1][1] = 0.0;
        
        Hr[0][0][0] = ( delta[1] / chiRho_HR[1] ) * ( 1.0 + ( dkap_dlnrho_face_HR[0] / opacity_HR[1] ) );
        Hr[0][0][1] = delta[1] * ( -4.0 + ( dkap_dlnT_face_HR[0] / opacity_HR[1] ) - ( chiT_HR[1]/chiRho_HR[1] )*( 1.0 + ( dkap_dlnrho_face_HR[0] / opacity_HR[1] ) ) )   -   ( (delta[0]*temperature_HR[0]) + (delta[1]*temperature_HR[1]) )/( temperature_HR[1] - temperature_HR[0] );
        Hr[0][1][0] = (1.0/( m*m*omega*omega*R )) * (   ( pressure_HR[1] * ( (delta[0]/rho_HR[0]) + (delta[1]/rho_HR[1]) ))/( rmid_cm_HR[1] - rmid_cm_HR[0] )   +   ( delta[1] * grav_HR[1] / chiRho_HR[1] )   );
        Hr[0][1][1] = - ( delta[1] * chiT_HR[1] * grav_HR[1] )/( m*m*omega*omega*R*chiRho_HR[1] );
        
        Hi[0][0][0] = 0.0;
        Hi[0][0][1] = 0.0;
        Hi[0][1][0] = 0.0;
        Hi[0][1][1] = 0.0;
        
        
        
        

        
        
        
        
        
        
        
        
        
        
        
        /*
         
         *
         *
         *
         *
         *
         *
         *
         *
         *
         *
         This is where I amend the matrices if there is convection going on, as a rough test.
         *
         *
         *
         *
         *
         *
         *
         *
         *
         *
         *
         
         */
        
        //test[k] = 0.0;
        
        //if (brunt_A_HR[0] <= 0.0) {
            
            // This chooses which is smaller - either F_{rad} / F_{tot} OR John Papaloizou's suggestion for looking at the convective efficiency
            /*
            if ( abs( gamma1_HR[0] * (brunt_A_HR[0] / radius_cm[0]) * (- 0.5) * ( pressure_HR[1] + pressure_HR[0] ) * ( rmid_cm_HR[1] - rmid_cm_HR[0] ) / ( pressure_HR[1] - pressure_HR[0] )) <= (flux_HR[0]/flux_BC) * ( (radius_cm_HR[0]*radius_cm_HR[0])/(R*R) ) ) {
                
                parameter_nonad = abs( gamma1_HR[0] * (brunt_A_HR[0] / radius_cm[0]) * (- 0.5) * ( pressure_HR[1] + pressure_HR[0] ) * ( rmid_cm_HR[1] - rmid_cm_HR[0] ) / ( pressure_HR[1] - pressure_HR[0] ));
                
            } else {
                
                parameter_nonad = (flux_HR[0]/flux_BC) * ( (radius_cm_HR[0]*radius_cm_HR[0])/(R*R) );
                
            }
            */
            // parameter_nonad may be referred to as lambda or \lambda in my emails - it's about suppressing the div(F') term
            //parameter_nonad = abs( gamma1_HR[0] * (brunt_A_HR[0] / radius_cm[0]) * (- 0.5) * ( pressure_HR[1] + pressure_HR[0] ) * ( rmid_cm_HR[1] - rmid_cm_HR[0] ) / ( pressure_HR[1] - pressure_HR[0] )); //(flux_HR[0]/flux_BC) * ( (radius_cm_HR[0]*radius_cm_HR[0])/(R*R) ); //0.0000000001; //0.000001; //0.0000000001;
            
            /*
            test[k] = parameter_nonad;
            
            location_nonad = radius_cm_HR[0] / R;
            */
            /*
            // This is an override to actually match what John is asking for
            if (radius_cm_HR[0] / R >= 0.9995) {
                
                parameter_nonad = 1.0;
                
            }
             */
            /*
            // Imaginary part is already 0
            Ar[0][1][1] = parameter_nonad * Ar[0][1][1];
            
            // Imaginary part is already 0
            Cr[0][1][1] = parameter_nonad * Cr[0][1][1];
            
            // Imaginary part contains the coefficient for \tilde{p}, which is not affected, so is left alone
            Dr[0][1][1] = parameter_nonad * Dr[0][1][1];
            */
        //}
        
        
        /*
         // This sets the equations to be adiabatic for x > 0.7 (that is, once the convective envelope starts)
         if ( 0.7 < radius_cm_HR[0]/R  && radius_cm_HR[0]/R < 1.1) {
         
         prop = 0.0;
         
         if ( 0.9 < radius_cm_HR[0]/R  && radius_cm_HR[0]/R < 0.95 ) {
         
         prop = 1.0 - ( ( (radius_cm_HR[0]/R) - 0.9) / 0.05);
         
         }
         
         
         
         
         if ( 0.85 < radius_cm_HR[0]/R  && radius_cm_HR[0]/R < 0.95 ) {
         
         prop = ( ( (radius_cm_HR[0]/R) - 0.85) / 0.1);
         
         }
         
         
         
         
         
         Ai[0][0][0] = prop*Ai[0][0][0];
         Ai[0][0][1] = prop*Ai[0][0][1];
         Ai[0][1][0] = prop*Ai[0][1][0];
         Ai[0][1][1] = prop*Ai[0][1][1];
         
         Ci[0][0][0] = prop*Ci[0][0][0];
         Ci[0][0][1] = prop*Ci[0][0][1];
         Ci[0][1][0] = prop*Ci[0][1][0];
         Ci[0][1][1] = prop*Ci[0][1][1];
         
         Di[0][0][0] = prop*Di[0][0][0];
         Di[0][0][1] = prop*Di[0][0][1];
         Di[0][1][0] = prop*Di[0][1][0];
         Di[0][1][1] = prop*Di[0][1][1];
         
         
         
         }
         */
        
        
        
        
        
        a[0] = sinh(4.0*(radius_cm_HR[0]/R))*sin(sinh(4.0*(radius_cm_HR[0]/R))) + 0.01*sin(1.0/(0.00001+(radius_cm_HR[0]/R)));
        a[1] = sinh(4.0*(radius_cm_HR[1]/R))*sin(sinh(4.0*(radius_cm_HR[1]/R))) + 0.01*sin(1.0/(0.00001+(radius_cm_HR[1]/R)));
        
        b[0] = (radius_cm_HR[0]/R)*sinh(radius_cm_HR[0]/R);
        b[1] = (radius_cm_HR[1]/R)*sinh(radius_cm_HR[1]/R);
        
        c[0] = 200.0 - 300.0*rmid_cm_HR[0]/R;
        c[1] = 200.0 - 300.0*rmid_cm_HR[1]/R;
        
        d[0] = 10.0 - 9.0*rmid_cm_HR[0]/R;
        d[1] = 10.0 - 9.0*rmid_cm_HR[1]/R;
        
        
        
        /*
         Mr[0][0][0] = (l*(l + 1.0)/(m*m*omega*omega*R*R))*rho_HR[1] * f * rmid_cm_HR[1] * rmid_cm_HR[1];
         Mr[0][1][0] = 0.0;
         
         Mi[0][0][0] = 0.0;
         Mi[0][1][0] = 0.0;
         
         
         Nr[0][0][0] = 0.0;
         Nr[0][1][0] = - 2.0*f*radius_cm_HR[0]*( (delta[0]*rho_HR[0]) + (delta[1]*rho_HR[1]) );
         
         
         Ni[0][0][0] = 0;
         Ni[0][1][0] = 0;
         */
        
        /*
         Mr[0][0][0] = Ar[0][0][0]*a[0] + Ar[0][0][1]*b[0] + Cr[0][0][0]*a[1] + Cr[0][0][1]*b[1] + Dr[0][0][0]*c[1] + Dr[0][0][1]*d[1];
         Mr[0][1][0] = Ar[0][1][0]*a[0] + Ar[0][1][1]*b[0] + Cr[0][1][0]*a[1] + Cr[0][1][1]*b[1] + Dr[0][1][0]*c[1] + Dr[0][1][1]*d[1];
         
         Mi[0][0][0] = 0.0;
         Mi[0][1][0] = 0.0;
         
         
         Nr[0][0][0] = Er[0][0][0]*a[0] + Er[0][0][1]*b[0] + Fr[0][0][0]*c[0] + Fr[0][0][1]*d[0] + Hr[0][0][0]*c[1] + Hr[0][0][1]*d[1];
         Nr[0][1][0] = Er[0][1][0]*a[0] + Er[0][1][1]*b[0] + Fr[0][1][0]*c[0] + Fr[0][1][1]*d[0] + Hr[0][1][0]*c[1] + Hr[0][1][1]*d[1];
         
         
         Ni[0][0][0] = 0;
         Ni[0][1][0] = 0;
         */
        
        
        /*
         Mr[0][0][0] = ( ( ( l*(l+1.0)*rmid_cm_HR[1]*rmid_cm_HR[1]*f )/( m*m*omega*omega*R*R ) ) )*( 1.0 + prop*sin(radius_cm_HR[0]*rmid_cm_HR[1]*rho_HR[0]*pressure_HR[1])); // Ar[0][0][0]*a[0] + Ar[0][0][1]*b[0] + Cr[0][0][0]*a[1] + Cr[0][0][1]*b[1] + Dr[0][0][0]*c[1] + Dr[0][0][1]*d[1];  // ( ( l*(l+1.0)*rmid_cm_HR[1]*rmid_cm_HR[1]*f )/( m*m*omega*omega*R*R ) );
         Mr[0][1][0] = 0.0; // Ar[0][1][0]*a[0] + Ar[0][1][1]*b[0] + Cr[0][1][0]*a[1] + Cr[0][1][1]*b[1] + Dr[0][1][0]*c[1] + Dr[0][1][1]*d[1];  // 0.0;
         
         Mi[0][0][0] = 0.0; // Ai[0][0][0]*a[0] + Ai[0][0][1]*b[0] + Ci[0][0][0]*a[1] + Ci[0][0][1]*b[1] + Di[0][0][0]*c[1] + Di[0][0][1]*d[1];  // 0.0;
         Mi[0][1][0] = 0.0; // Ai[0][1][0]*a[0] + Ai[0][1][1]*b[0] + Ci[0][1][0]*a[1] + Ci[0][1][1]*b[1] + Di[0][1][0]*c[1] + Di[0][1][1]*d[1];  // 0.0;
         
         
         Nr[0][0][0] = 0.0; // Er[0][0][0]*a[0] + Er[0][0][1]*b[0] + Fr[0][0][0]*c[0] + Fr[0][0][1]*d[0] + Hr[0][0][0]*c[1] + Hr[0][0][1]*d[1];  // 0.0;
         Nr[0][1][0] = ( - 2.0 * f * radius_cm_HR[0] / ( m*m*omega*omega*R ) ) * ( 1.0 + prop*sin(temperature_HR[1]*rho_HR[0]*radius_cm_HR[0]*temperature_HR[0]*rho_HR[1])); // Er[0][1][0]*a[0] + Er[0][1][1]*b[0] + Fr[0][1][0]*c[0] + Fr[0][1][1]*d[0] + Hr[0][1][0]*c[1] + Hr[0][1][1]*d[1];  // - 2.0 * f * radius_cm_HR[0] / ( m*m*omega*omega*R );
         
         Ni[0][0][0] = 0.0; // Ei[0][0][0]*a[0] + Ei[0][0][1]*b[0] + Fi[0][0][0]*c[0] + Fi[0][0][1]*d[0] + Hi[0][0][0]*c[1] + Hi[0][0][1]*d[1];  // 0.0;
         Ni[0][1][0] = 0.0; // Ei[0][1][0]*a[0] + Ei[0][1][1]*b[0] + Fi[0][1][0]*c[0] + Fi[0][1][1]*d[0] + Hi[0][1][0]*c[1] + Hi[0][1][1]*d[1];  // 0.0;
         */
        
        
        
        
        
        Mr[0][0][0] = ( ( ( l*(l+1.0)*rmid_cm_HR[1]*rmid_cm_HR[1]*f )/( m*m*omega*omega*R*R ) ) ); //*( 1.0 + prop*sin(radius_cm_HR[0]*rmid_cm_HR[1]*rho_HR[0]*pressure_HR[1])); // Ar[0][0][0]*a[0] + Ar[0][0][1]*b[0] + Cr[0][0][0]*a[1] + Cr[0][0][1]*b[1] + Dr[0][0][0]*c[1] + Dr[0][0][1]*d[1];  // ( ( l*(l+1.0)*rmid_cm_HR[1]*rmid_cm_HR[1]*f )/( m*m*omega*omega*R*R ) );
        Mr[0][1][0] = 0.0; // Ar[0][1][0]*a[0] + Ar[0][1][1]*b[0] + Cr[0][1][0]*a[1] + Cr[0][1][1]*b[1] + Dr[0][1][0]*c[1] + Dr[0][1][1]*d[1];  // 0.0;
        
        Mi[0][0][0] = 0.0; // Ai[0][0][0]*a[0] + Ai[0][0][1]*b[0] + Ci[0][0][0]*a[1] + Ci[0][0][1]*b[1] + Di[0][0][0]*c[1] + Di[0][0][1]*d[1];  // 0.0;
        Mi[0][1][0] = 0.0; // Ai[0][1][0]*a[0] + Ai[0][1][1]*b[0] + Ci[0][1][0]*a[1] + Ci[0][1][1]*b[1] + Di[0][1][0]*c[1] + Di[0][1][1]*d[1];  // 0.0;
        
        
        Nr[0][0][0] = 0.0; // Er[0][0][0]*a[0] + Er[0][0][1]*b[0] + Fr[0][0][0]*c[0] + Fr[0][0][1]*d[0] + Hr[0][0][0]*c[1] + Hr[0][0][1]*d[1];  // 0.0;
        Nr[0][1][0] = ( - 2.0 * f * radius_cm_HR[0] / ( m*m*omega*omega*R ) ); // * ( 1.0 + prop*sin(temperature_HR[1]*rho_HR[0]*radius_cm_HR[0]*temperature_HR[0]*rho_HR[1])); // Er[0][1][0]*a[0] + Er[0][1][1]*b[0] + Fr[0][1][0]*c[0] + Fr[0][1][1]*d[0] + Hr[0][1][0]*c[1] + Hr[0][1][1]*d[1];  // - 2.0 * f * radius_cm_HR[0] / ( m*m*omega*omega*R );
        
        Ni[0][0][0] = 0.0; // Ei[0][0][0]*a[0] + Ei[0][0][1]*b[0] + Fi[0][0][0]*c[0] + Fi[0][0][1]*d[0] + Hi[0][0][0]*c[1] + Hi[0][0][1]*d[1];  // 0.0;
        Ni[0][1][0] = 0.0; // Ei[0][1][0]*a[0] + Ei[0][1][1]*b[0] + Fi[0][1][0]*c[0] + Fi[0][1][1]*d[0] + Hi[0][1][0]*c[1] + Hi[0][1][1]*d[1];  // 0.0;
        
        
        
        
        
        
        
        
        
        
        
        
        // Here I amend the matrix elements to account for convection, done in nice easy stages.
        
        // Stage 1: only varying entropy
        
        // prefactor gives  (dr/dT0) * (F_conv / K0) / ( ( dlnT - grada dlnp ) )
        // Note that I have cancelled out the (rmid[1] - rmid[0]) bits, to streamline things a bit
        
        //prefactor = ( flux_tot_HR[0] * conv_L_div_L_HR[0] * ( rmid_cm_HR[1] - rmid_cm_HR[0] ) * ( (delta[0]/K_HR[0]) + (delta[1]/K_HR[1]) ) / ( temperature_HR[1] - temperature_HR[0] ) ) / ( (lnT_HR[1] - lnT_HR[0]) - (delta[0]*grada_HR[0] + delta[1]*grada_HR[1])*( log(pressure_HR[1]) - log(pressure_HR[0]) ) );
        
        prefactor = ( flux_tot_HR[0] * conv_L_div_L_HR[0] * ( (delta[0]/(K_HR[0] * cp_HR[0] * (brunt_A_HR[0]/rmid_cm_HR[0]) * (chiRho_HR[0]/chiT_HR[0]) )) + (delta[1]/(K_HR[1] * cp_HR[1] * (brunt_A_HR[1]/rmid_cm_HR[1]) * (chiRho_HR[1]/chiT_HR[1]) ))) / ( temperature_HR[1] - temperature_HR[0] ) );
        
        
        x_convection = 24.0 * 5.67e-5 * (delta[0]*temperature_HR[0] + delta[1]*temperature_HR[1]) * (delta[0]*temperature_HR[0] + delta[1]*temperature_HR[1]) * (delta[0]*temperature_HR[0] + delta[1]*temperature_HR[1]) / (   (delta[0]*rho_HR[0] + delta[1]*rho_HR[1])*(delta[0]*rho_HR[0] + delta[1]*rho_HR[1]) * (delta[0]*cp_HR[0] + delta[1]*cp_HR[1]) * 2.0*scale_height_cm_HR[0] * conv_vel_HR[0] * (delta[0]*opacity_HR[0] + delta[1]*opacity_HR[1])   );
        
        
        
        if( brunt_A_HR[0] <= 0.0 ) {
            /*
             Fr[0][0][0] = Fr[0][0][0] + ( prefactor * ( - delta[0]*grada_HR[0] - delta[1]*grada_HR[1] ) * ( ((-1.0) * pressure_HR[0] * ( (delta[0]/pressure_HR[0]) + (delta[1]/pressure_HR[1]) )) - ((delta[0]/pressure_HR[0])*( pressure_HR[1] - pressure_HR[0] )) ) );
             
             Fr[0][0][1] = Fr[0][0][1] + ( prefactor * ( ((-1.0) * temperature_HR[0] * ( (delta[0]/temperature_HR[0]) + (delta[1]/temperature_HR[1]) )) - ((delta[0]/temperature_HR[0])*( temperature_HR[1] - temperature_HR[0] )) ) );
             
             Hr[0][0][0] = Hr[0][0][0] + ( prefactor * ( - delta[0]*grada_HR[0] - delta[1]*grada_HR[1] ) * ( ((+1.0) * pressure_HR[1] * ( (delta[0]/pressure_HR[0]) + (delta[1]/pressure_HR[1]) )) - ((delta[1]/pressure_HR[1])*( pressure_HR[1] - pressure_HR[0] )) ) );
             
             Hr[0][0][1] = Hr[0][0][1] + ( prefactor * ( ((+1.0) * temperature_HR[1] * ( (delta[0]/temperature_HR[0]) + (delta[1]/temperature_HR[1]) )) - ((delta[1]/temperature_HR[1])*( temperature_HR[1] - temperature_HR[0] )) ) );
             */
            
            // ds'/ds0
            
            Fr[0][0][0] = Fr[0][0][0] + ( prefactor * cp_HR[0] * grada_HR[0] );
            
            Fr[0][0][1] = Fr[0][0][1] - ( prefactor * cp_HR[0] );
            
            Hr[0][0][0] = Hr[0][0][0] - ( prefactor * cp_HR[1] * grada_HR[1] );
            
            Hr[0][0][1] = Hr[0][0][1] + ( prefactor * cp_HR[1] );
            
            
            /*
            // this adds in the other terms (that is, T'/T_{0}, rho'/rho_{0}, kappa'/kappa_{0}, but NOT v'_{c},v_{c,0})
            
             Fr[0][0][0] = Fr[0][0][0] + ( (  ( (radius_cm_HR[1] - radius_cm_HR[0]) / ( temperature_HR[1] - temperature_HR[0] ) ) * flux_tot_HR[0] * conv_L_div_L_HR[0] * ( (delta[0]/K_HR[0]) + (delta[1]/K_HR[1]) ) / ( 1.0 + x_convection )  ) * (          ( 1.0 + ( 3.0 + dkap_dlnrho_face_HR[0]/opacity_HR[0] )*x_convection )*(  delta[0] / chiRho_HR[0]  )          ) );
             
             Fr[0][0][1] = Fr[0][0][1] + ( (  ( (radius_cm_HR[1] - radius_cm_HR[0]) / ( temperature_HR[1] - temperature_HR[0] ) ) * flux_tot_HR[0] * conv_L_div_L_HR[0] * ( (delta[0]/K_HR[0]) + (delta[1]/K_HR[1]) ) / ( 1.0 + x_convection )  ) * (          (  delta[0]  )*( 1.0 - 2.0*x_convection - (1.0 + 3.0*x_convection)*(chiT_HR[0]/chiRho_HR[0])  +  (x_convection/opacity_HR[0])*( dkap_dlnT_face_HR[0] - (chiT_HR[0]/chiRho_HR[0])*dkap_dlnrho_face_HR[0] )           ) ) );
             
             Hr[0][0][0] = Hr[0][0][0] + ( (  ( (radius_cm_HR[1] - radius_cm_HR[0]) / ( temperature_HR[1] - temperature_HR[0] ) ) * flux_tot_HR[0] * conv_L_div_L_HR[0] * ( (delta[0]/K_HR[0]) + (delta[1]/K_HR[1]) ) / ( 1.0 + x_convection )  ) * (          ( 1.0 + ( 3.0 + dkap_dlnrho_face_HR[0]/opacity_HR[0] )*x_convection )*(  delta[1] / chiRho_HR[1]  )          ) );
             
             Hr[0][0][1] = Hr[0][0][1] + ( (  ( (radius_cm_HR[1] - radius_cm_HR[0]) / ( temperature_HR[1] - temperature_HR[0] ) ) * flux_tot_HR[0] * conv_L_div_L_HR[0] * ( (delta[0]/K_HR[0]) + (delta[1]/K_HR[1]) ) / ( 1.0 + x_convection )  ) * (          (  delta[1]  )*( 1.0 - 2.0*x_convection - (1.0 + 3.0*x_convection)*(chiT_HR[1]/chiRho_HR[1])  +  (x_convection/opacity_HR[0])*( dkap_dlnT_face_HR[0] - (chiT_HR[1]/chiRho_HR[1])*dkap_dlnrho_face_HR[0] )           ) ) );
             
            
            
            // this adds in the bit for v'_{c}
            
            Fr[0][0][0] = Fr[0][0][0] + (   (1.0 + 2.0*x_convection)/(2.0*(1.0 + x_convection))   )*( prefactor * cp_HR[0] * grada_HR[0] );
            
            Fr[0][0][1] = Fr[0][0][1] - (   (1.0 + 2.0*x_convection)/(2.0*(1.0 + x_convection))   )*( prefactor * cp_HR[0] );
            
            Hr[0][0][0] = Hr[0][0][0] - (   (1.0 + 2.0*x_convection)/(2.0*(1.0 + x_convection))   )*( prefactor * cp_HR[1] * grada_HR[1] );
            
            Hr[0][0][1] = Hr[0][0][1] + (   (1.0 + 2.0*x_convection)/(2.0*(1.0 + x_convection))   )*( prefactor * cp_HR[1] );
            
            Nr[0][0][0] =  ( (  ( (radius_cm_HR[1] - radius_cm_HR[0]) / ( temperature_HR[1] - temperature_HR[0] ) ) * flux_tot_HR[0] * conv_L_div_L_HR[0] * ( (delta[0]/K_HR[0]) + (delta[1]/K_HR[1]) ) / ( 1.0 + x_convection )  ) * (          f * radius_cm_HR[0] * ( 1.0/grav_HR[0] )          ) );
            */
            
            
            
            
            
        }
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        if (k == 0) {
            
            // This makes dumMA = C^{-1}
            CMatrixInv(Cr,Ci,dumMAr,dumMAi,0,0);
            
            // This makes dumMB = alpha_{0}
            CMatrixMult(dumMAr,dumMAi,Dr,Di,dumMBr,dumMBi,0,0,0);
            
            // This makes dumVA = - M
            dumVAr[0][0][0] = - Mr[0][0][0];
            dumVAr[0][1][0] = - Mr[0][1][0];
            
            dumVAi[0][0][0] = - Mi[0][0][0];
            dumVAi[0][1][0] = - Mi[0][1][0];
            
            // This makes dumVB = gamma_{0}
            CVectorMult(dumMAr,dumMAi,dumVAr,dumVAi,dumVBr,dumVBi,0,0,0);
            
            T_ar = dumVBr[0][0][0];
            
            CompDiv(&dumMBr[0][1][1], &dumMBi[0][1][1], &dumMBr[0][0][1], &dumMBi[0][0][1], &dummyr, &dummyi);
            // This assigns T_b
            CompMult(&dummyr, &dummyi, &dumVBr[0][0][0], &dumVBi[0][0][0], &T_br, &T_bi);
            
            // This assigns T_c
            CompDiv(&dumVBr[0][0][0], &dumVBi[0][0][0], &dumMBr[0][0][0], &dumMBi[0][0][0], &T_cr, &T_ci);
            
            // This assigns T_d
            CompDiv(&dumVBr[0][0][0], &dumVBi[0][0][0], &dumMBr[0][0][1], &dumMBi[0][0][1], &T_dr, &T_di);
            
            cout << "T_a = " << T_ar << " + " << T_ai << "i\n";
            cout << "T_b = " << T_br << " + " << T_bi << "i\n";
            cout << "T_c = " << T_cr << " + " << T_ci << "i\n";
            cout << "T_d = " << T_dr << " + " << T_di << "i\n";
            
            
            // If this section is active, then this prevents any rescaling
            T_ar = 1.0;
            T_ai = 0.0;
            T_br = 1.0;
            T_bi = 0.0;
            T_cr = 1.0;
            T_ci = 0.0;
            T_dr = 1.0;
            T_di = 0.0;
            
            
        }
        
        
        // This rescales the A - H matrices
        
        for (i = 0; i < 2; i = i + 1) {
            
            
            // This rescales the T_a bits
            
            CompMult(&Ar[0][i][0], &Ai[0][i][0], &T_ar, &T_ai, &dummyr, &dummyi);
            
            Ar[0][i][0] = dummyr;
            Ai[0][i][0] = dummyi;
            
            
            CompMult(&Cr[0][i][0], &Ci[0][i][0], &T_ar, &T_ai, &dummyr, &dummyi);
            
            Cr[0][i][0] = dummyr;
            Ci[0][i][0] = dummyi;
            
            
            CompMult(&Er[0][i][0], &Ei[0][i][0], &T_ar, &T_ai, &dummyr, &dummyi);
            
            Er[0][i][0] = dummyr;
            Ei[0][i][0] = dummyi;
            
            
            // This rescales the T_b bits
            
            CompMult(&Ar[0][i][1], &Ai[0][i][1], &T_br, &T_bi, &dummyr, &dummyi);
            
            Ar[0][i][1] = dummyr;
            Ai[0][i][1] = dummyi;
            
            
            CompMult(&Cr[0][i][1], &Ci[0][i][1], &T_br, &T_bi, &dummyr, &dummyi);
            
            Cr[0][i][1] = dummyr;
            Ci[0][i][1] = dummyi;
            
            
            CompMult(&Er[0][i][1], &Ei[0][i][1], &T_br, &T_bi, &dummyr, &dummyi);
            
            Er[0][i][1] = dummyr;
            Ei[0][i][1] = dummyi;
            
            
            
            // This rescales the T_c bits
            
            CompMult(&Dr[0][i][0], &Di[0][i][0], &T_cr, &T_ci, &dummyr, &dummyi);
            
            Dr[0][i][0] = dummyr;
            Di[0][i][0] = dummyi;
            
            
            CompMult(&Fr[0][i][0], &Fi[0][i][0], &T_cr, &T_ci, &dummyr, &dummyi);
            
            Fr[0][i][0] = dummyr;
            Fi[0][i][0] = dummyi;
            
            
            CompMult(&Hr[0][i][0], &Hi[0][i][0], &T_cr, &T_ci, &dummyr, &dummyi);
            
            Hr[0][i][0] = dummyr;
            Hi[0][i][0] = dummyi;
            
            
            // This rescales the T_d bits
            
            CompMult(&Dr[0][i][1], &Di[0][i][1], &T_dr, &T_di, &dummyr, &dummyi);
            
            Dr[0][i][1] = dummyr;
            Di[0][i][1] = dummyi;
            
            
            CompMult(&Fr[0][i][1], &Fi[0][i][1], &T_dr, &T_di, &dummyr, &dummyi);
            
            Fr[0][i][1] = dummyr;
            Fi[0][i][1] = dummyi;
            
            
            CompMult(&Hr[0][i][1], &Hi[0][i][1], &T_dr, &T_di, &dummyr, &dummyi);
            
            Hr[0][i][1] = dummyr;
            Hi[0][i][1] = dummyi;
            
            
            
            
            
            
        }
        
        
        
        
        
        
        
        // This sets the boundary condition matrices
        
        if (k == J-2) {
            
            // Option 1
            /*
             etar[0][0][0] = 0.5*dp_dr_BC;
             etar[0][0][1] = 0;
             etar[0][1][0] = 0.5*( ( (log(flux_HR[1]) - log(flux_HR[0]))/(radius_cm_HR[1] - radius_cm_HR[0]) )  -  4.0*dlnT_dr_BC );
             etar[0][1][1] = 0.5/flux_HR[0];
             
             etai[0][0][0] = 0;
             etai[0][0][1] = 0;
             etai[0][1][0] = 0;
             etai[0][1][1] = 0;
             
             mur[0][0][0] = 0.5*dp_dr_BC;
             mur[0][0][1] = 0;
             mur[0][1][0] = 0.5*( ( (log(flux_HR[1]) - log(flux_HR[0]))/(radius_cm_HR[1] - radius_cm_HR[0]) )  -  4.0*dlnT_dr_BC );
             mur[0][1][1] = 0.5/flux_HR[1];
             
             mui[0][0][0] = 0;
             mui[0][0][1] = 0;
             mui[0][1][0] = 0;
             mui[0][1][1] = 0;
             
             nur[0][0][0] = 1.0;
             nur[0][0][1] = 0;
             nur[0][1][0] = 0;
             nur[0][1][1] = -4.0/temperature_HR[1];
             
             nui[0][0][0] = 0;
             nui[0][0][1] = 0;
             nui[0][1][0] = 0;
             nui[0][1][1] = 0;
             */
            
            
            
            
            // Option 2  -- NB: I'm currently suspicious of these
            /*
             etar[0][0][0] = 0.5 * dp_dr_BC * radius_cm_HR[0] / pressure_HR[1];
             etar[0][0][1] = 0.0;
             etar[0][1][0] = 0.5 * radius_cm_HR[0] * ( ( (log(flux_HR[1]) - log(flux_HR[0]))/(radius_cm_HR[1] - radius_cm_HR[0]) )  -  4.0*dlnT_dr_BC );
             etar[0][1][1] = 0.5 * flux_BC / flux_HR[0]; // 0.5;
             
             etai[0][0][0] = 0.0;
             etai[0][0][1] = 0.0;
             etai[0][1][0] = 0.0;
             etai[0][1][1] = 0.0;
             
             mur[0][0][0] = 0.5 * dp_dr_BC * radius_cm_HR[1] / pressure_HR[1];
             mur[0][0][1] = 0.0;
             mur[0][1][0] = 0.5 * radius_cm_HR[1] * ( ( (log(flux_HR[1]) - log(flux_HR[0]))/(radius_cm_HR[1] - radius_cm_HR[0]) )  -  4.0*dlnT_dr_BC );
             mur[0][1][1] = 0.5 * flux_BC / flux_HR[1]; // 0.5;
             
             mui[0][0][0] = 0.0;
             mui[0][0][1] = 0.0;
             mui[0][1][0] = 0.0;
             mui[0][1][1] = 0.0;
             
             nur[0][0][0] = 1.0;
             nur[0][0][1] = 0.0;
             nur[0][1][0] = 0.0;
             nur[0][1][1] = -4.0;
             
             nui[0][0][0] = 0.0;
             nui[0][0][1] = 0.0;
             nui[0][1][0] = 0.0;
             nui[0][1][1] = 0.0;
             */
            
            
            
            
            // Option 2.i
            /*
             etar[0][0][0] = 0.5 * dp_dr_BC * R / pressure_HR[1];
             etar[0][0][1] = 0.0;
             etar[0][1][0] = 0.5 * R * ( ( (log(flux_HR[1]) - log(flux_HR[0]))/(radius_cm_HR[1] - radius_cm_HR[0]) )  -  4.0*dlnT_dr_BC );
             etar[0][1][1] = 0.5 * flux_BC / flux_HR[0]; // 0.5;
             
             etai[0][0][0] = 0.0;
             etai[0][0][1] = 0.0;
             etai[0][1][0] = 0.0;
             etai[0][1][1] = 0.0;
             
             mur[0][0][0] = 0.5 * dp_dr_BC * R / pressure_HR[1];
             mur[0][0][1] = 0.0;
             mur[0][1][0] = 0.5 * R * ( ( (log(flux_HR[1]) - log(flux_HR[0]))/(radius_cm_HR[1] - radius_cm_HR[0]) )  -  4.0*dlnT_dr_BC );
             mur[0][1][1] = 0.5 * flux_BC / flux_HR[1]; // 0.5;
             
             mui[0][0][0] = 0.0;
             mui[0][0][1] = 0.0;
             mui[0][1][0] = 0.0;
             mui[0][1][1] = 0.0;
             
             nur[0][0][0] = 1.0;
             nur[0][0][1] = 0.0;
             nur[0][1][0] = 0.0;
             nur[0][1][1] = -4.0;
             
             nui[0][0][0] = 0.0;
             nui[0][0][1] = 0.0;
             nui[0][1][0] = 0.0;
             nui[0][1][1] = 0.0;
             */
            
            
            
            
            
            
            // Option 2.i, but with the first BC being that the lagrangian derivative of rho = 0 (instead of p)
            /*
             
             etar[0][0][0] = 0.5 * dlnRho_dr_BC * R;
             etar[0][0][1] = 0.0;
             etar[0][1][0] = 0.5 * R * ( ( (log(flux_HR[1]) - log(flux_HR[0]))/(radius_cm_HR[1] - radius_cm_HR[0]) )  -  4.0*dlnT_dr_BC );
             etar[0][1][1] = 0.5 * flux_BC / flux_HR[0]; // 0.5;
             
             etai[0][0][0] = 0.0;
             etai[0][0][1] = 0.0;
             etai[0][1][0] = 0.0;
             etai[0][1][1] = 0.0;
             
             mur[0][0][0] = 0.5 * dlnRho_dr_BC * R;
             mur[0][0][1] = 0.0;
             mur[0][1][0] = 0.5 * R * ( ( (log(flux_HR[1]) - log(flux_HR[0]))/(radius_cm_HR[1] - radius_cm_HR[0]) )  -  4.0*dlnT_dr_BC );
             mur[0][1][1] = 0.5 * flux_BC / flux_HR[1]; // 0.5;
             
             mui[0][0][0] = 0.0;
             mui[0][0][1] = 0.0;
             mui[0][1][0] = 0.0;
             mui[0][1][1] = 0.0;
             
             nur[0][0][0] = 1.0/chiRho_HR[1];
             nur[0][0][1] = - chiT_HR[1] / chiRho_HR[1];
             nur[0][1][0] = 0.0;
             nur[0][1][1] = -4.0;
             
             nui[0][0][0] = 0.0;
             nui[0][0][1] = 0.0;
             nui[0][1][0] = 0.0;
             nui[0][1][1] = 0.0;
             
             */
            
            // option 3 -- JP gave this to me saying that it avoids singularity at the surface, which would be good
            /*
             
             etar[0][0][0] = 0.5 * R * ( (-omega*omega/grav_HR[0]) + (4.0/radius_cm_HR[0]) - ((l*(l+1.0)*grav_HR[0])/(radius_cm_HR[0]*radius_cm_HR[0]*omega*omega)) );
             etar[0][0][1] = 0.0;
             etar[0][1][0] = 0.5 * R * ( ( (log(flux_HR[1]) - log(flux_HR[0]))/(radius_cm_HR[1] - radius_cm_HR[0]) )  -  4.0*dlnT_dr_BC );
             etar[0][1][1] = 0.5 * flux_BC / flux_HR[0]; // 0.5;
             
             etai[0][0][0] = 0.0;
             etai[0][0][1] = 0.0;
             etai[0][1][0] = 0.0;
             etai[0][1][1] = 0.0;
             
             mur[0][0][0] = 0.5 * R * ( (-omega*omega/grav_HR[1]) + (4.0/radius_cm_HR[1]) - ((l*(l+1.0)*grav_HR[1])/(radius_cm_HR[1]*radius_cm_HR[1]*omega*omega)) );
             mur[0][0][1] = 0.0;
             mur[0][1][0] = 0.5 * R * ( ( (log(flux_HR[1]) - log(flux_HR[0]))/(radius_cm_HR[1] - radius_cm_HR[0]) )  -  4.0*dlnT_dr_BC );
             mur[0][1][1] = 0.5 * flux_BC / flux_HR[1]; // 0.5;
             
             mui[0][0][0] = 0.0;
             mui[0][0][1] = 0.0;
             mui[0][1][0] = 0.0;
             mui[0][1][1] = 0.0;
             
             nur[0][0][0] = -(1.0 + ((pressure_HR[1]*l*(l+1.0))/(rho_HR[1]*rmid_cm_HR[1]*rmid_cm_HR[1]*omega*omega)) );
             nur[0][0][1] = 0.0;
             nur[0][1][0] = 0.0;
             nur[0][1][1] = -4.0;
             
             nui[0][0][0] = 0.0;
             nui[0][0][1] = 0.0;
             nui[0][1][0] = 0.0;
             nui[0][1][1] = 0.0;
             
             */
            
            
            
            
            
            // option 4 -- same BC as Pfahl et al -- hydrostatic condition
            /*
             
             etar[0][0][0] = 0.5 * R * ( ((pressure_HR[1] - pressure_HR[0])/(rmid_cm_HR[1] - rmid_cm_HR[0]))  -  (2.0/(grav_HR[1] + grav_HR[0]))*((grav_HR[1] - grav_HR[0])/(rmid_cm_HR[1] - rmid_cm_HR[0]))  +  ((opacity_HR[1] - opacity_HR[0])/(rmid_cm_HR[1] - rmid_cm_HR[0])) )   -   f*R/grav_HR[0];
             etar[0][0][1] = 0.0;
             etar[0][1][0] = 0.5 * R * ( ( (log(flux_HR[1]) - log(flux_HR[0]))/(radius_cm_HR[1] - radius_cm_HR[0]) )  -  4.0*dlnT_dr_BC );
             etar[0][1][1] = 0.5 * flux_BC / flux_HR[0]; // 0.5;
             
             etai[0][0][0] = 0.0;
             etai[0][0][1] = 0.0;
             etai[0][1][0] = 0.0;
             etai[0][1][1] = 0.0;
             
             mur[0][0][0] = 0.5 * R * ( ((pressure_HR[1] - pressure_HR[0])/(rmid_cm_HR[1] - rmid_cm_HR[0]))  -  (2.0/(grav_HR[1] + grav_HR[0]))*((grav_HR[1] - grav_HR[0])/(rmid_cm_HR[1] - rmid_cm_HR[0]))  +  ((opacity_HR[1] - opacity_HR[0])/(rmid_cm_HR[1] - rmid_cm_HR[0])) )   -   f*R/grav_HR[1];
             mur[0][0][1] = 0.0;
             mur[0][1][0] = 0.5 * R * ( ( (log(flux_HR[1]) - log(flux_HR[0]))/(radius_cm_HR[1] - radius_cm_HR[0]) )  -  4.0*dlnT_dr_BC );
             mur[0][1][1] = 0.5 * flux_BC / flux_HR[1]; // 0.5;
             
             mui[0][0][0] = 0.0;
             mui[0][0][1] = 0.0;
             mui[0][1][0] = 0.0;
             mui[0][1][1] = 0.0;
             
             nur[0][0][0] = 1.0  +  dkap_dlnrho_face_HR[1]/(opacity_HR[1]*chiRho_HR[1]);
             nur[0][0][1] = (  dkap_dlnT_face_HR[1]  -  (chiT_HR[1]/chiRho_HR[1])*dkap_dlnrho_face_HR[1]  ) / opacity_HR[1];
             nur[0][1][0] = 0.0;
             nur[0][1][1] = -4.0;
             
             nui[0][0][0] = 0.0;
             nui[0][0][1] = 0.0;
             nui[0][1][0] = 0.0;
             nui[0][1][1] = 0.0;
             */
            
            
            
            

            // Default is 2i
            

            etar[0][0][0] = 0.5 * dp_dr_BC * R / pressure_HR[1];
            etar[0][0][1] = 0.0;
            etar[0][1][0] = 0.5 * R * ( ( (log(flux_HR[1]) - log(flux_HR[0]))/(radius_cm_HR[1] - radius_cm_HR[0]) )  -  4.0*dlnT_dr_BC );
            etar[0][1][1] = 0.5 * flux_BC / flux_HR[0]; // 0.5;
            
            etai[0][0][0] = 0.0;
            etai[0][0][1] = 0.0;
            etai[0][1][0] = 0.0;
            etai[0][1][1] = 0.0;
            
            mur[0][0][0] = 0.5 * dp_dr_BC * R / pressure_HR[1];
            mur[0][0][1] = 0.0;
            mur[0][1][0] = 0.5 * R * ( ( (log(flux_HR[1]) - log(flux_HR[0]))/(radius_cm_HR[1] - radius_cm_HR[0]) )  -  4.0*dlnT_dr_BC );
            mur[0][1][1] = 0.5 * flux_BC / flux_HR[1]; // 0.5;
            
            mui[0][0][0] = 0.0;
            mui[0][0][1] = 0.0;
            mui[0][1][0] = 0.0;
            mui[0][1][1] = 0.0;
            
            nur[0][0][0] = 1.0;
            nur[0][0][1] = 0.0;
            nur[0][1][0] = 0.0;
            nur[0][1][1] = -4.0;
            
            nui[0][0][0] = 0.0;
            nui[0][0][1] = 0.0;
            nui[0][1][0] = 0.0;
            nui[0][1][1] = 0.0;

            
            
            
            
            
            
            
            /*
             xr[0][0][0] = 0.0;
             xr[0][1][0] = 0.0;
             
             xi[0][0][0] = 0.0;
             xi[0][1][0] = 0.0;
             */
            
            xr[0][0][0] = 0.0; // etar[0][0][0]*a[0] + etar[0][0][1]*b[0] + mur[0][0][0]*a[1] + mur[0][0][1]*b[1] + nur[0][0][0]*c[1] + nur[0][0][1]*d[1];  // 0.0;
            xr[0][1][0] = 0.0; // etar[0][1][0]*a[0] + etar[0][1][1]*b[0] + mur[0][1][0]*a[1] + mur[0][1][1]*b[1] + nur[0][1][0]*c[1] + nur[0][1][1]*d[1];  // 0.0;
            
            xi[0][0][0] = 0.0; // etai[0][0][0]*a[0] + etai[0][0][1]*b[0] + mui[0][0][0]*a[1] + mui[0][0][1]*b[1] + nui[0][0][0]*c[1] + nui[0][0][1]*d[1];  // 0.0;
            xi[0][1][0] = 0.0; // etai[0][1][0]*a[0] + etai[0][1][1]*b[0] + mui[0][1][0]*a[1] + mui[0][1][1]*b[1] + nui[0][1][0]*c[1] + nui[0][1][1]*d[1];  // 0.0;
            
            
            
            
            
            
            
            // This rescales the eta - mu matrices
            
            for (i = 0; i < 2; i = i + 1) {
                
                
                // This rescales the T_a bits
                
                CompMult(&etar[0][i][0], &etai[0][i][0], &T_ar, &T_ai, &dummyr, &dummyi);
                
                etar[0][i][0] = dummyr;
                etai[0][i][0] = dummyi;
                
                
                CompMult(&mur[0][i][0], &mui[0][i][0], &T_ar, &T_ai, &dummyr, &dummyi);
                
                mur[0][i][0] = dummyr;
                mui[0][i][0] = dummyi;
                
                
                
                
                // This rescales the T_b bits
                
                CompMult(&etar[0][i][1], &etai[0][i][1], &T_br, &T_bi, &dummyr, &dummyi);
                
                etar[0][i][1] = dummyr;
                etai[0][i][1] = dummyi;
                
                
                CompMult(&mur[0][i][1], &mui[0][i][1], &T_br, &T_bi, &dummyr, &dummyi);
                
                mur[0][i][1] = dummyr;
                mui[0][i][1] = dummyi;
                
                
                
                
                // This rescales the T_c bits
                
                CompMult(&nur[0][i][0], &nui[0][i][0], &T_cr, &T_ci, &dummyr, &dummyi);
                
                nur[0][i][0] = dummyr;
                nui[0][i][0] = dummyi;
                
                
                
                
                
                // This rescales the T_d bits
                
                CompMult(&nur[0][i][1], &nui[0][i][1], &T_dr, &T_di, &dummyr, &dummyi);
                
                nur[0][i][1] = dummyr;
                nui[0][i][1] = dummyi;
                
                
                
                
                
                
                
                
            }
            
            
            
            
            
            
        }
        
        
        
        
        
        
        
        
        /*
         
         This writes the data to the file, in the order:
         1 - radius_cm / R
         2, 3, 4, 5 - A
         6, 7, 8, 9 - C
         10,11,12,13- D
         14,15,16,17- E
         18,19,20,21- F
         22,23,24,25- H
         26, 27     - M
         28, 29     - N
         30 - detA
         31 - detC
         32 - det D
         33 - detE
         34 - detF
         35 - detH
         
         */
        
        matrixfile << radius_cm_HR_output[k]/R << "\t\t\t" << Ar[0][0][0] << "\t\t" << Ar[0][0][1] << "\t\t" << Ar[0][1][0] << "\t\t" << Ar[0][1][1] << "\t\t" << Cr[0][0][0] << "\t\t" << Cr[0][0][1] << "\t\t" << Cr[0][1][0] << "\t\t" << Cr[0][1][1] << "\t\t" << Dr[0][0][0] << "\t\t" << Dr[0][0][1] << "\t\t" << Dr[0][1][0] << "\t\t" << Dr[0][1][1] << "\t\t" << Er[0][0][0] << "\t\t" << Er[0][0][1] << "\t\t" << Er[0][1][0] << "\t\t" << Er[0][1][1] << "\t\t" << Fr[0][0][0] << "\t\t" << Fr[0][0][1] << "\t\t" << Fr[0][1][0] << "\t\t" << Fr[0][1][1] << "\t\t" << Hr[0][0][0] << "\t\t" << Hr[0][0][1] << "\t\t" << Hr[0][1][0] << "\t\t" << Hr[0][1][1] << "\t\t" << Mr[0][0][0] << "\t\t" << Mr[0][1][0] << "\t\t" << Nr[0][0][0] << "\t\t" << Nr[0][1][0] << "\t\t" << Ar[0][0][0]*Ar[0][1][1] - Ar[0][0][1]*Ar[0][1][0] << "\t\t" << Cr[0][0][0]*Cr[0][1][1] - Cr[0][0][1]*Cr[0][1][0] << "\t\t" << Dr[0][0][0]*Dr[0][1][1] - Dr[0][0][1]*Dr[0][1][0] << "\t\t" << Er[0][0][0]*Er[0][1][1] - Er[0][0][1]*Er[0][1][0] << "\t\t" << Fr[0][0][0]*Fr[0][1][1] - Fr[0][0][1]*Fr[0][1][0] << "\t\t" << Hr[0][0][0]*Hr[0][1][1] - Hr[0][0][1]*Hr[0][1][0] << "\t\t" << "\n";
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        // Here the matrix magic happens - they must be multiplied and inverted such that they fit the useful equations above, once again in a sum up to J
        
        
        
        
        
        
        // NOTE - as this calculates alpha_{i+1} from i things, it runs from k=0 to k=J-2, and the final step (with k=J-2) is used just to initialise the matrices for use with the BCs after the for loop closes.
        
        
        // This initialises P
        
        CMatrixInv(Ar,Ai,dumMAr,dumMAi,0,0);
        
        CMatrixMult(Er,Ei,dumMAr,dumMAi,Pr,Pi,0,0,0);
        
        /*        cout << "\nThis is Pr_{" << k << "} \n";
         cout << Pr[k][0][0] << "\t" << Pr[k][0][1] << "\n";
         cout << Pr[k][1][0] << "\t" << Pr[k][1][1] << "\n\n";
         */
        
        // This initialises Q - need to do it for each iteration of k, because it depends on alpha!
        
        
        CMatrixInv(Fr,Fi,dumMAr,dumMAi,0,0);
        
        CMatrixMult(alphar,alphai,dumMAr,dumMAi,Qr,Qi,k,0,0);
        
        /*        cout << "\nThis is Qr_{" << k << "} \n";
         cout << Qr[k][0][0] << "\t" << Qr[k][0][1] << "\n";
         cout << Qr[k][1][0] << "\t" << Qr[k][1][1] << "\n\n";
         */
        
        //This initialises R - need to do it for each iteration of k, because it depends on alpha!
        
        
        CMatrixMult(Qr,Qi,Er,Ei,dumMAr,dumMAi,0,0,0);
        
        /*        cout << "\nThis is (QE)r_{" << k << "} \n";
         cout << dumMAr[0][0][0] << "\t" << dumMAr[0][0][1] << "\n";
         cout << dumMAr[0][1][0] << "\t" << dumMAr[0][1][1] << "\n\n";
         */
        dumMCr[0][0][0] = dumMAr[0][0][0] - 1.0;
        dumMCr[0][0][1] = dumMAr[0][0][1] - 0.0;
        dumMCr[0][1][0] = dumMAr[0][1][0] - 0.0;
        dumMCr[0][1][1] = dumMAr[0][1][1] - 1.0;
        
        dumMCi[0][0][0] = dumMAi[0][0][0] - 0.0;
        dumMCi[0][0][1] = dumMAi[0][0][1] - 0.0;
        dumMCi[0][1][0] = dumMAi[0][1][0] - 0.0;
        dumMCi[0][1][1] = dumMAi[0][1][1] - 0.0;
        
        CMatrixInv(dumMCr,dumMCi,Rr,Ri,0,0);
        
        /*        cout << "\nThis is Rr_{" << k << "} \n";
         cout << Rr[k][0][0] << "\t" << Rr[k][0][1] << "\n";
         cout << Rr[k][1][0] << "\t" << Rr[k][1][1] << "\n\n";
         */
        
        
        /*
         
         // Here the recursion arrays are calculated for this value of i
         
         CMatrixInv(Ar,Ai,dumMAr,dumMAi,0,0);
         
         
         dumMBr[0][0][0] = - dumMAr[0][0][0];
         dumMBr[0][0][1] = - dumMAr[0][0][1];
         dumMBr[0][1][0] = - dumMAr[0][1][0];
         dumMBr[0][1][1] = - dumMAr[0][1][1];
         
         dumMBi[0][0][0] = - dumMAi[0][0][0];
         dumMBi[0][0][1] = - dumMAi[0][0][1];
         dumMBi[0][1][0] = - dumMAi[0][1][0];
         dumMBi[0][1][1] = - dumMAi[0][1][1];
         
         
         //This assigns RECu
         CMatrixMult(dumMBr,dumMBi,Cr,Ci,RECur,RECui,0,0,k);
         
         //This assigns RECv
         CMatrixMult(dumMBr,dumMBi,Dr,Di,RECvr,RECvi,0,0,k);
         
         //This assigns RECc
         CVectorMult(dumMAr,dumMAi,Mr,Mi,RECcr,RECci,0,0,k);
         
         */
        
        
        
        // Case I equations (remember that you also need to change u to v when these are actually implemented)
        
        
        CMatrixMult(Er,Ei,alphar,alphai,dumMAr,dumMAi,0,k,0);
        
        //        cout << "E = " << Er[0][0][0] << "\t" << Er[0][0][1] << "\t" << Er[0][1][0] << "\t" << Er[0][1][1] << "\n";
        
        //        cout << "alpha = " << alphar[k][0][0] << "\t" << alphar[k][0][1] << "\t" << alphar[k][1][0] << "\t" << alphar[k][1][1] << "\n";
        
        //        cout << "dumMA = " << dumMAr[0][0][0] << "\t" << dumMAr[0][1][0] << "\t" << dumMAr[0][1][0] << "\t" << dumMAr[0][1][1] << "\n";
        
        dumMBr[0][0][0] = Fr[0][0][0] - dumMAr[0][0][0];
        dumMBr[0][0][1] = Fr[0][0][1] - dumMAr[0][0][1];
        dumMBr[0][1][0] = Fr[0][1][0] - dumMAr[0][1][0];
        dumMBr[0][1][1] = Fr[0][1][1] - dumMAr[0][1][1];
        
        dumMBi[0][0][0] = Fi[0][0][0] - dumMAi[0][0][0];
        dumMBi[0][0][1] = Fi[0][0][1] - dumMAi[0][0][1];
        dumMBi[0][1][0] = Fi[0][1][0] - dumMAi[0][1][0];
        dumMBi[0][1][1] = Fi[0][1][1] - dumMAi[0][1][1];
        
        //        cout << "dumMB = " << dumMBr[0][0][0] << "\t" << dumMBr[0][1][0] << "\t" << dumMBr[0][1][0] << "\t" << dumMBr[0][1][1] << "\n";
        
        // dumMA = ( F - E alpha_{k})^{-1}
        CMatrixInv(dumMBr,dumMBi,dumMAr,dumMAi,0,0);
        
        //        cout << "dumMA (second use) = " << dumMAr[0][0][0] << "\t" << dumMAr[0][1][0] << "\t" << dumMAr[0][1][0] << "\t" << dumMAr[0][1][1] << "\n";
        
        // This gives RECu = 0
        RECur[k][0][0] = 0.0;
        RECur[k][0][1] = 0.0;
        RECur[k][1][0] = 0.0;
        RECur[k][1][1] = 0.0;
        
        RECui[k][0][0] = 0.0;
        RECui[k][0][1] = 0.0;
        RECui[k][1][0] = 0.0;
        RECui[k][1][1] = 0.0;
        
        dumMBr[0][0][0] = - Hr[0][0][0];
        dumMBr[0][0][1] = - Hr[0][0][1];
        dumMBr[0][1][0] = - Hr[0][1][0];
        dumMBr[0][1][1] = - Hr[0][1][1];
        
        dumMBi[0][0][0] = - Hi[0][0][0];
        dumMBi[0][0][1] = - Hi[0][0][1];
        dumMBi[0][1][0] = - Hi[0][1][0];
        dumMBi[0][1][1] = - Hi[0][1][1];
        
        //        cout << "dumMB (second use) = " << dumMBr[0][0][0] << "\t" << dumMBr[0][0][1] << "\t" << dumMBr[0][1][0] << "\t" << dumMBr[0][1][1] << "\n";
        
        // This gives RECv = - ( F - E alpha_{k})^{-1} H
        CMatrixMult(dumMAr,dumMAi,dumMBr,dumMBi,RECvr,RECvi,0,0,k);
        
        CVectorMult(Er,Ei,gammar,gammai,dumVAr,dumVAi,0,k,0);
        
        dumVBr[0][0][0] = Nr[0][0][0] + dumVAr[0][0][0];
        dumVBr[0][1][0] = Nr[0][1][0] + dumVAr[0][1][0];
        
        dumVBi[0][0][0] = Ni[0][0][0] + dumVAi[0][0][0];
        dumVBi[0][1][0] = Ni[0][1][0] + dumVAi[0][1][0];
        
        // This gives RECc = ( F - E alpha_{k})^{-1} ( N + E gamma_{k} )
        CVectorMult(dumMAr,dumMAi,dumVBr,dumVBi,RECcr,RECci,0,0,k);
        
        
        
        
        
        
        // Now we can iterate up to get alpha and gamma for all positions in the star
        
        //alpha_{i+1}
        
        
        CMatrixMult(Ar,Ai,Rr,Ri,dumMAr,dumMAi,0,0,0); // dumMA wil remain AR until the end of calculating gamma_{i+1}
        
        CMatrixMult(Qr,Qi,Hr,Hi,dumMBr,dumMBi,0,0,0);
        
        CMatrixMult(dumMAr,dumMAi,dumMBr,dumMBi,dumMCr,dumMCi,0,0,0);
        
        dumMDr[0][0][0] = Dr[0][0][0] - dumMCr[0][0][0];
        dumMDr[0][0][1] = Dr[0][0][1] - dumMCr[0][0][1];
        dumMDr[0][1][0] = Dr[0][1][0] - dumMCr[0][1][0];
        dumMDr[0][1][1] = Dr[0][1][1] - dumMCr[0][1][1];
        
        dumMDi[0][0][0] = Di[0][0][0] - dumMCi[0][0][0];
        dumMDi[0][0][1] = Di[0][0][1] - dumMCi[0][0][1];
        dumMDi[0][1][0] = Di[0][1][0] - dumMCi[0][1][0];
        dumMDi[0][1][1] = Di[0][1][1] - dumMCi[0][1][1];
        
        CMatrixInv(Cr,Ci,dumMCr,dumMCi,0,0); // dumMC is C^{-1} and will remain to be so for the calculation of gamma_{i+1} as well
        
        CMatrixMult(dumMCr,dumMCi,dumMDr,dumMDi,alphar,alphai,0,0,k+1);
        
        
        
        //gamma_{i+1}
        
        
        CVectorMult(Qr,Qi,Nr,Ni,dumVAr,dumVAi,0,0,0);
        
        dumVBr[0][0][0] = gammar[k][0][0] + dumVAr[0][0][0];
        dumVBr[0][1][0] = gammar[k][1][0] + dumVAr[0][1][0];
        
        dumVBi[0][0][0] = gammai[k][0][0] + dumVAi[0][0][0];
        dumVBi[0][1][0] = gammai[k][1][0] + dumVAi[0][1][0];
        
        CVectorMult(dumMAr,dumMAi,dumVBr,dumVBi,dumVAr,dumVAi,0,0,0); /* Note that I am reusing dumVA, but as an output once again. */
        
        dumVDr[0][0][0] = dumVAr[0][0][0] - Mr[0][0][0];
        dumVDr[0][1][0] = dumVAr[0][1][0] - Mr[0][1][0];
        
        dumVDi[0][0][0] = dumVAi[0][0][0] - Mi[0][0][0];
        dumVDi[0][1][0] = dumVAi[0][1][0] - Mi[0][1][0];
        
        CVectorMult(dumMCr,dumMCi,dumVDr,dumVDi,gammar,gammai,0,0,k+1);
        
        //cout << " " << k << "\t" << "Done a matrix magic iteration \n";
        
        
        
    } // This closes the for loop
    
    
    
    
    // This closes the output data file for the Matrix data
    matrixfile.close();
    

    
    
    
    
    
    cout << "Flag Exited Matrix Magic \n";
    
    
    cout << "Outermost centre of cell is at r/R = " << rmid_cm_HR[1]/R << "\n";
    
    
    // Here the outer boundary conditions must be used in order to get v_{J-1}
    
    
    
    
    
    
    
    
    
    
    
    /*
     
     
     ********************************
     *                              *
     *                              *
     *                              *
     *                              *
     *        ORIGINAL METHOD       *
     *                              *
     *                              *
     *                              *
     *                              *
     ********************************
     
     
     
     
     The outer boundary conditions can be expressed in the form of a vector equation as:
     
     \eta u_{J-2} + \mu u_{J-1} + \nu v_{J-1} = x
     
     which is valid for the centre of the outermost cell, labelled J-1
     
     If this is used in coordination with the oscillation equations and the u-v relation for J-2, as:
     
     A_{J-2,J-1} u_{J-2} + C_{J-2,J-1} u_{J-1} + D_{J-2,J-1} v_{J-1} = M_{J-2,J-1}
     
     E_{J-2,J-1} u_{J-2} + F_{J-2,J-1} v_{J-2} + H_{J-2,J-1} v_{J-1} = N_{J-2,J-1}
     
     u_{J-2} + alpha_{J-2} v_{J-2} + gamma_{J-2} = 0
     
     we can eliminate everything except v_{J-1}, which gives us the (very long) expression:
     
     v_{J-1} = [ C^{-1} D - mu^{-1} nu + ( mu^{-1} eta - C^{-1} A ) ( alpha_{J-2} F^{-1} E - 1 )^{-1} alpha_{J-2} F^{-1} H ]^{-1}
     [ C^{-1} M - mu^{-1} x + ( mu^{-1} eta - C^{-1} A ) ( alpha_{J-2} F^{-1} E - 1 )^{-1} ( alpha_{J-2} F^{-1} N + gamma_{J-2} ) ]
     
     Introducing the following variables breaks down the expression somewhat:
     
     R = [ alpha_{i} F^{-1} E - 1 ]^{-1} = [ Q E - 1 ]^{-1}  --  NB: this has already been defined and used earlier, this is just a reminder
     
     BCa = ( mu^{-1} eta - C^{-1} A ) R
     
     BCb = mu^{-1}  --  NB this is to reduce re-doing the inversion, which will marginally speed things up
     
     BCc = C^{-1}  --  NB this is just the same reasoning as for BCb
     
     
     which gives
     
     v_{J-1} = [ BCc D - BCb nu + BCa Q H ]^{-1} [ BCc M - BCb x + BCa ( Q N + gamma_{J-2} ) ]
     
     
     */
    
    
    
    // Boundary condition simplification variables calculated here
    
    // BCb defined here
    CMatrixInv(mur,mui,BCbr,BCbi,0,0);
    
    CMatrixInv(Cr,Ci,dumMAr,dumMAi,0,0);
    
    // BCc defined here
    BCcr[0][0][0] = dumMAr[0][0][0];
    BCcr[0][0][1] = dumMAr[0][0][1];
    BCcr[0][1][0] = dumMAr[0][1][0];
    BCcr[0][1][1] = dumMAr[0][1][1];
    
    BCci[0][0][0] = dumMAi[0][0][0];
    BCci[0][0][1] = dumMAi[0][0][1];
    BCci[0][1][0] = dumMAi[0][1][0];
    BCci[0][1][1] = dumMAi[0][1][1];
    
    
    CMatrixMult(BCbr,BCbi,etar,etai,dumMAr,dumMAi,0,0,0);
    
    CMatrixMult(BCcr,BCci,Ar,Ai,dumMBr,dumMBi,0,0,0);
    
    dumMCr[0][0][0] = dumMAr[0][0][0] - dumMBr[0][0][0];
    dumMCr[0][0][1] = dumMAr[0][0][1] - dumMBr[0][0][1];
    dumMCr[0][1][0] = dumMAr[0][1][0] - dumMBr[0][1][0];
    dumMCr[0][1][1] = dumMAr[0][1][1] - dumMBr[0][1][1];
    
    dumMCi[0][0][0] = dumMAi[0][0][0] - dumMBi[0][0][0];
    dumMCi[0][0][1] = dumMAi[0][0][1] - dumMBi[0][0][1];
    dumMCi[0][1][0] = dumMAi[0][1][0] - dumMBi[0][1][0];
    dumMCi[0][1][1] = dumMAi[0][1][1] - dumMBi[0][1][1];
    
    // BCa defined here
    CMatrixMult(dumMCr,dumMCi,Rr,Ri,BCar,BCai,0,0,0);
    
    
    // This initialises P
    
    CMatrixInv(Ar,Ai,dumMAr,dumMAi,0,0);
    
    CMatrixMult(Er,Ei,dumMAr,dumMAi,Pr,Pi,0,0,0);
    
    /*        cout << "\nThis is Pr_{" << k << "} \n";
     cout << Pr[k][0][0] << "\t" << Pr[k][0][1] << "\n";
     cout << Pr[k][1][0] << "\t" << Pr[k][1][1] << "\n\n";
     */
    
    // This initialises Q - need to do it for each iteration of k, because it depends on alpha!
    
    
    CMatrixInv(Fr,Fi,dumMAr,dumMAi,0,0);
    
    CMatrixMult(alphar,alphai,dumMAr,dumMAi,Qr,Qi,k,0,0);
    
    /*        cout << "\nThis is Qr_{" << k << "} \n";
     cout << Qr[k][0][0] << "\t" << Qr[k][0][1] << "\n";
     cout << Qr[k][1][0] << "\t" << Qr[k][1][1] << "\n\n";
     */
    
    //This initialises R - need to do it for each iteration of k, because it depends on alpha!
    
    
    CMatrixMult(Qr,Qi,Er,Ei,dumMAr,dumMAi,0,0,0);
    
    /*        cout << "\nThis is (QE)r_{" << k << "} \n";
     cout << dumMAr[0][0][0] << "\t" << dumMAr[0][0][1] << "\n";
     cout << dumMAr[0][1][0] << "\t" << dumMAr[0][1][1] << "\n\n";
     */
    dumMCr[0][0][0] = dumMAr[0][0][0] - 1.0;
    dumMCr[0][0][1] = dumMAr[0][0][1] - 0.0;
    dumMCr[0][1][0] = dumMAr[0][1][0] - 0.0;
    dumMCr[0][1][1] = dumMAr[0][1][1] - 1.0;
    
    dumMCi[0][0][0] = dumMAi[0][0][0] - 0.0;
    dumMCi[0][0][1] = dumMAi[0][0][1] - 0.0;
    dumMCi[0][1][0] = dumMAi[0][1][0] - 0.0;
    dumMCi[0][1][1] = dumMAi[0][1][1] - 0.0;
    
    CMatrixInv(dumMCr,dumMCi,Rr,Ri,0,0);
    
    /*        cout << "\nThis is Rr_{" << k << "} \n";
     cout << Rr[k][0][0] << "\t" << Rr[k][0][1] << "\n";
     cout << Rr[k][1][0] << "\t" << Rr[k][1][1] << "\n\n";
     */
    
    
    
    
    
    
    
    
    // Here we calculate v_{J-1}
    
    CMatrixMult(BCcr,BCci,Dr,Di,dumMAr,dumMAi,0,0,0);
    
    /*    cout << "(BCc D)r \n";
     cout << dumMAr[0][0][0] << "\t" << dumMAr[0][0][1] << "\n";
     cout << dumMAr[0][1][0] << "\t" << dumMAr[0][1][1] << "\n\n";
     */
    CMatrixMult(BCbr,BCbi,nur,nui,dumMBr,dumMBi,0,0,0);
    
    /*    cout << "(BCb nu)r \n";
     cout << dumMBr[0][0][0] << "\t" << dumMBr[0][0][1] << "\n";
     cout << dumMBr[0][1][0] << "\t" << dumMBr[0][1][1] << "\n\n";
     */
    CMatrixMult(Qr,Qi,Hr,Hi,dumMCr,dumMCi,0,0,0);
    
    CMatrixMult(BCar,BCai,dumMCr,dumMCi,dumMDr,dumMDi,0,0,0);
    
    /*    cout << "(BCa Q H)r \n";
     cout << dumMDr[0][0][0] << "\t" << dumMDr[0][0][1] << "\n";
     cout << dumMDr[0][1][0] << "\t" << dumMDr[0][1][1] << "\n\n";
     */
    dumMEr[0][0][0] = dumMAr[0][0][0] - dumMBr[0][0][0] + dumMDr[0][0][0];
    dumMEr[0][0][1] = dumMAr[0][0][1] - dumMBr[0][0][1] + dumMDr[0][0][1];
    dumMEr[0][1][0] = dumMAr[0][1][0] - dumMBr[0][1][0] + dumMDr[0][1][0];
    dumMEr[0][1][1] = dumMAr[0][1][1] - dumMBr[0][1][1] + dumMDr[0][1][1];
    
    dumMEi[0][0][0] = dumMAi[0][0][0] - dumMBi[0][0][0] + dumMDi[0][0][0];
    dumMEi[0][0][1] = dumMAi[0][0][1] - dumMBi[0][0][1] + dumMDi[0][0][1];
    dumMEi[0][1][0] = dumMAi[0][1][0] - dumMBi[0][1][0] + dumMDi[0][1][0];
    dumMEi[0][1][1] = dumMAi[0][1][1] - dumMBi[0][1][1] + dumMDi[0][1][1];
    
    // This gets us to: dumMA[0] = ( BCc D - BCb nu + BCa Q H )^{-1}
    CMatrixInv(dumMEr,dumMEi,dumMAr,dumMAi,0,0);
    
    /*    cout << "( BCc D - BCb nu + BCa Q H )^{-1} r \n";
     cout << dumMAr[0][0][0] << "\t" << dumMAr[0][0][1] << "\n";
     cout << dumMAr[0][1][0] << "\t" << dumMAr[0][1][1] << "\n\n";
     */
    
    CVectorMult(BCcr,BCci,Mr,Mi,dumVAr,dumVAi,0,0,0);
    
    CVectorMult(BCbr,BCbi,xr,xi,dumVBr,dumVBi,0,0,0);
    
    CVectorMult(Qr,Qi,Nr,Ni,dumVCr,dumVCi,0,0,0);
    
    dumVDr[0][0][0] = dumVCr[0][0][0] + gammar[J-2][0][0];
    dumVDr[0][1][0] = dumVCr[0][1][0] + gammar[J-2][1][0];
    
    dumVDi[0][0][0] = dumVCi[0][0][0] + gammai[J-2][0][0];
    dumVDi[0][1][0] = dumVCi[0][1][0] + gammai[J-2][1][0];
    
    CVectorMult(BCar,BCai,dumVDr,dumVDi,dumVCr,dumVCi,0,0,0);
    
    // This gives us: BCc M - BCb x + BCa ( Q N + gamma_{J-2} )
    dumVDr[0][0][0] = dumVAr[0][0][0] - dumVBr[0][0][0] + dumVCr[0][0][0];
    dumVDr[0][1][0] = dumVAr[0][1][0] - dumVBr[0][1][0] + dumVCr[0][1][0];
    
    dumVDi[0][0][0] = dumVAi[0][0][0] - dumVBi[0][0][0] + dumVCi[0][0][0];
    dumVDi[0][1][0] = dumVAi[0][1][0] - dumVBi[0][1][0] + dumVCi[0][1][0];
    
    
    
    
    
    
    // This step actually gets us v_{J-1}
    CVectorMult(dumMAr,dumMAi,dumVDr,dumVDi,vr,vi,0,0,J-1);
    
    cout << "Original method:\n" << "Got vr_{J-1} \t \n";
    
    
    
    
    
    // We use v_{J-1} to calculate u_{J-1} too
    // alpha_{i} v_{i}
    CVectorMult(alphar,alphai,vr,vi,dumVAr,dumVAi,J-1,J-1,0);
    // Gives u_{J-1}
    ur[J-1][0][0] = - dumVAr[0][0][0] - gammar[J-1][0][0];
    ur[J-1][1][0] = - dumVAr[0][1][0] - gammar[J-1][1][0];
    
    ui[J-1][0][0] = - dumVAi[0][0][0] - gammai[J-1][0][0];
    ui[J-1][1][0] = - dumVAi[0][1][0] - gammai[J-1][1][0];
    
    cout << "Got u_{J-1} ( J-1 = " << J-1 <<" ) \n";
    
    cout << "\nThis is ur_{J-1} \n";
    cout << ur[J-1][0][0] << "\n";
    cout << ur[J-1][1][0] << "\n\n";
    
    cout << "\nThis is vr_{J-1} \n";
    cout << vr[J-1][0][0] << "\n";
    cout << vr[J-1][1][0] << "\n\n";
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*
     
     
     ********************************
     *                              *
     *                              *
     *                              *
     *                              *
     *           METHOD I           *
     *                              *
     *                              *
     *                              *
     *                              *
     ********************************
     
     
     
     
     The outer boundary conditions can be expressed in the form of a vector equation as:
     
     \eta u_{J-2} + \mu u_{J-1} + \nu v_{J-1} = x
     
     which is valid for the centre of the outermost cell, labelled J-1
     
     If this is used in coordination with the oscillation equations and the u-v relation for J-2, as:
     
     A_{J-2,J-1} u_{J-2} + C_{J-2,J-1} u_{J-1} + D_{J-2,J-1} v_{J-1} = M_{J-2,J-1}
     
     E_{J-2,J-1} u_{J-2} + F_{J-2,J-1} v_{J-2} + H_{J-2,J-1} v_{J-1} = N_{J-2,J-1}
     
     u_{J-2} + alpha_{J-2} v_{J-2} + gamma_{J-2} = 0
     
     we can eliminate everything except v_{J-1}, which gives us the (very long) expression:
     
     v_{J-1} = [ H - ( E - F alpha_{J-2}^{-1} ) ( A - C mu^{-1} eta )^{-1} ( D - C mu^{-1} nu ) ]^{-1} [ N_{J-2} + F alpha_{J-2}^{-1} gamma_{J-2} - ( E - F alpha_{J-2}^{-1} ) ( A - C mu^{-1} eta )^{-1} ( M - C mu^{-1} x ) ]
     
     
     Introducing the following variables breaks down the expression somewhat:
     
     BCa = F alpha_{J-2}^{-1}
     
     BCb = C mu^{-1}
     
     BCc = ( E - F alpha_{J-2}^{-1} ) ( A - C mu^{-1} eta )^{-1} = ( E - BCa ) ( A - BCb eta )^{-1}
     
     
     which gives
     
     v_{J-1} = [ H - BCc ( D - BCb nu ) ]^{-1} [ N_{J-2} + BCa gamma_{J-2} - BCc ( M_{J-2} - BCb x ) ]
     
     
     */
    
    
    
    // Boundary condition simplification variables calculated here
    
    
    // BCa defined here
    //cout << "FLAG line 3000\n";
    CMatrixInv(alphar,alphai,dumMAr,dumMAi,J-2,0);
    
    CMatrixMult(Fr,Fi,dumMAr,dumMAi,BCar,BCai,0,0,0);
    
    
    // BCb defined here
    
    CMatrixInv(mur,mui,dumMAr,dumMAi,0,0);
    
    CMatrixMult(Cr,Ci,dumMAr,dumMAi,BCbr,BCbi,0,0,0);
    
    
    // BCc defined here
    
    CMatrixMult(BCbr,BCbi,etar,etai,dumMAr,dumMAi,0,0,0);
    
    dumMCr[0][0][0] = Ar[0][0][0] - dumMAr[0][0][0];
    dumMCr[0][0][1] = Ar[0][0][1] - dumMAr[0][0][1];
    dumMCr[0][1][0] = Ar[0][1][0] - dumMAr[0][1][0];
    dumMCr[0][1][1] = Ar[0][1][1] - dumMAr[0][1][1];
    
    dumMCi[0][0][0] = Ai[0][0][0] - dumMAi[0][0][0];
    dumMCi[0][0][1] = Ai[0][0][1] - dumMAi[0][0][1];
    dumMCi[0][1][0] = Ai[0][1][0] - dumMAi[0][1][0];
    dumMCi[0][1][1] = Ai[0][1][1] - dumMAi[0][1][1];
    
    CMatrixInv(dumMCr,dumMCi,dumMBr,dumMBi,0,0);
    
    dumMAr[0][0][0] = Er[0][0][0] - BCar[0][0][0];
    dumMAr[0][0][1] = Er[0][0][1] - BCar[0][0][1];
    dumMAr[0][1][0] = Er[0][1][0] - BCar[0][1][0];
    dumMAr[0][1][1] = Er[0][1][1] - BCar[0][1][1];
    
    dumMAi[0][0][0] = Ei[0][0][0] - BCai[0][0][0];
    dumMAi[0][0][1] = Ei[0][0][1] - BCai[0][0][1];
    dumMAi[0][1][0] = Ei[0][1][0] - BCai[0][1][0];
    dumMAi[0][1][1] = Ei[0][1][1] - BCai[0][1][1];
    
    CMatrixMult(dumMAr,dumMAi,dumMBr,dumMBi,BCcr,BCci,0,0,0);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    // Here we calculate v_{J-1}
    
    
    CMatrixMult(BCbr,BCbi,nur,nui,dumMAr,dumMAi,0,0,0);
    
    dumMBr[0][0][0] = Dr[0][0][0] - dumMAr[0][0][0];
    dumMBr[0][0][1] = Dr[0][0][1] - dumMAr[0][0][1];
    dumMBr[0][1][0] = Dr[0][1][0] - dumMAr[0][1][0];
    dumMBr[0][1][1] = Dr[0][1][1] - dumMAr[0][1][1];
    
    dumMBi[0][0][0] = Di[0][0][0] - dumMAi[0][0][0];
    dumMBi[0][0][1] = Di[0][0][1] - dumMAi[0][0][1];
    dumMBi[0][1][0] = Di[0][1][0] - dumMAi[0][1][0];
    dumMBi[0][1][1] = Di[0][1][1] - dumMAi[0][1][1];
    
    CMatrixMult(BCcr,BCci,dumMBr,dumMBi,dumMAr,dumMAi,0,0,0);
    
    dumMBr[0][0][0] = Hr[0][0][0] - dumMAr[0][0][0];
    dumMBr[0][0][1] = Hr[0][0][1] - dumMAr[0][0][1];
    dumMBr[0][1][0] = Hr[0][1][0] - dumMAr[0][1][0];
    dumMBr[0][1][1] = Hr[0][1][1] - dumMAr[0][1][1];
    
    dumMBi[0][0][0] = Hi[0][0][0] - dumMAi[0][0][0];
    dumMBi[0][0][1] = Hi[0][0][1] - dumMAi[0][0][1];
    dumMBi[0][1][0] = Hi[0][1][0] - dumMAi[0][1][0];
    dumMBi[0][1][1] = Hi[0][1][1] - dumMAi[0][1][1];
    
    // This gives: dumMA  ==  [ H - ( E - F alpha_{J-2}^{-1} ) ( A - C mu^{-1} eta )^{-1} ( D - C mu^{-1} nu ) ]^{-1}  ==  [ H - BCc ( D - BCb nu ) ]^{-1}
    CMatrixInv(dumMBr,dumMBi,dumMCr,dumMCi,0,0);
    
    
    
    CVectorMult(BCar,BCai,gammar,gammai,dumVAr,dumVAi,0,J-2,0);
    
    CVectorMult(BCbr,BCbi,xr,xi,dumVCr,dumVCi,0,0,0);
    
    dumVBr[0][0][0] = Mr[0][0][0] - dumVCr[0][0][0];
    dumVBr[0][1][0] = Mr[0][1][0] - dumVCr[0][1][0];
    
    dumVBi[0][0][0] = Mi[0][0][0] - dumVCi[0][0][0];
    dumVBi[0][1][0] = Mi[0][1][0] - dumVCi[0][1][0];
    
    CVectorMult(BCcr,BCci,dumVBr,dumVBi,dumVCr,dumVCi,0,0,0);
    
    // This gives: dumVD  ==  N_{J-2} + F alpha_{J-2}^{-1} gamma_{J-2} - ( E - F alpha_{J-2}^{-1} ) ( A - C mu^{-1} eta )^{-1} ( M - C mu^{-1} x )  ==  [ N_{J-2} + BCa gamma_{J-2} - BCc ( M_{J-2} - BCb x ) ]
    dumVDr[0][0][0] = Nr[0][0][0] + dumVAr[0][0][0] - dumVCr[0][0][0];
    dumVDr[0][1][0] = Nr[0][1][0] + dumVAr[0][1][0] - dumVCr[0][1][0];
    
    dumVDi[0][0][0] = Ni[0][0][0] + dumVAi[0][0][0] - dumVCi[0][0][0];
    dumVDi[0][1][0] = Ni[0][1][0] + dumVAi[0][1][0] - dumVCi[0][1][0];
    
    
    
    // This step actually gets us v_{J-1}
    CVectorMult(dumMCr,dumMCi,dumVDr,dumVDi,vr,vi,0,0,J-1);
    
    cout << "Method I:\n" << "Got vr_{J-1} \t \n";
    
    /*    cout << "\nThis is vr_{J-1} \n";
     cout << vr[J-1][0][0] << "\n";
     cout << vr[J-1][1][0] << "\n\n";
     */
    
    
    
    // We use v_{J-1} to calculate u_{J-1} too
    // alpha_{i} v_{i}
    CVectorMult(alphar,alphai,vr,vi,dumVAr,dumVAi,J-1,J-1,0);
    // Gives u_{J-1}
    ur[J-1][0][0] = - dumVAr[0][0][0] - gammar[J-1][0][0];
    ur[J-1][1][0] = - dumVAr[0][1][0] - gammar[J-1][1][0];
    
    ui[J-1][0][0] = - dumVAi[0][0][0] - gammai[J-1][0][0];
    ui[J-1][1][0] = - dumVAi[0][1][0] - gammai[J-1][1][0];
    
    cout << "Got u_{J-1} ( J-1 = " << J-1 <<" ) \n";
    
    cout << "\nThis is ur_{J-1} \n";
    cout << ur[J-1][0][0] << "\n";
    cout << ur[J-1][1][0] << "\n\n";
    
    cout << "\nThis is vr_{J-1} \n";
    cout << vr[J-1][0][0] << "\n";
    cout << vr[J-1][1][0] << "\n\n";
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*
     
     
     ********************************
     *                              *
     *                              *
     *                              *
     *                              *
     *          METHOD II           *
     *                              *
     *                              *
     *                              *
     *                              *
     ********************************
     
     
     
     
     The outer boundary conditions can be expressed in the form of a vector equation as:
     
     \eta u_{J-2} + \mu u_{J-1} + \nu v_{J-1} = x
     
     which is valid for the centre of the outermost cell, labelled J-1
     
     If this is used in coordination with the oscillation equations and the u-v relation for J-2, as:
     
     A_{J-2,J-1} u_{J-2} + C_{J-2,J-1} u_{J-1} + D_{J-2,J-1} v_{J-1} = M_{J-2,J-1}
     
     E_{J-2,J-1} u_{J-2} + F_{J-2,J-1} v_{J-2} + H_{J-2,J-1} v_{J-1} = N_{J-2,J-1}
     
     u_{J-2} + alpha_{J-2} v_{J-2} + gamma_{J-2} = 0
     
     we can eliminate everything except v_{J-1}, which gives us the (very long) expression:
     
     v_{J-1} = [ H - ( E - F alpha_{J-2}^{-1} ) ( eta - mu C^{-1} A )^{-1} ( nu - mu C^{-1} D ) ]^{-1} [ N_{J-2} + F alpha_{J-2}^{-1} gamma_{J-2} - ( E - F alpha_{J-2}^{-1} ) ( eta - mu C^{-1} A )^{-1} ( x - mu C^{-1} M ) ]
     
     
     Introducing the following variables breaks down the expression somewhat:
     
     BCa = F alpha_{J-2}^{-1}
     
     BCb = mu C^{-1}
     
     BCc = ( E - F alpha_{J-2}^{-1} ) ( eta - mu C^{-1} A )^{-1} = ( E - BCa ) ( eta - BCb A )^{-1}
     
     
     which gives
     
     v_{J-1} = [ H - BCc ( nu - BCb D ) ]^{-1} [ N_{J-2} + BCa gamma_{J-2} - BCc ( x - BCb M_{J-2} ) ]
     
     
     */
    
    
    
    // Boundary condition simplification variables calculated here
    
    
    // BCa defined here
    
    CMatrixInv(alphar,alphai,dumMAr,dumMAi,J-2,0);
    cout << "Flag inverse 1\n";
    
    CMatrixMult(Fr,Fi,dumMAr,dumMAi,BCar,BCai,0,0,0);
    
    
    // BCb defined here
    
    CMatrixInv(Cr,Ci,dumMAr,dumMAi,0,0);
    cout << "Flag inverse 2\n";
    
    CMatrixMult(mur,mui,dumMAr,dumMAi,BCbr,BCbi,0,0,0);
    
    
    // BCc defined here
    
    CMatrixMult(BCbr,BCbi,Ar,Ai,dumMAr,dumMAi,0,0,0);
    
    dumMCr[0][0][0] = etar[0][0][0] - dumMAr[0][0][0];
    dumMCr[0][0][1] = etar[0][0][1] - dumMAr[0][0][1];
    dumMCr[0][1][0] = etar[0][1][0] - dumMAr[0][1][0];
    dumMCr[0][1][1] = etar[0][1][1] - dumMAr[0][1][1];
    
    dumMCi[0][0][0] = etai[0][0][0] - dumMAi[0][0][0];
    dumMCi[0][0][1] = etai[0][0][1] - dumMAi[0][0][1];
    dumMCi[0][1][0] = etai[0][1][0] - dumMAi[0][1][0];
    dumMCi[0][1][1] = etai[0][1][1] - dumMAi[0][1][1];
    
    CMatrixInv(dumMCr,dumMCi,dumMBr,dumMBi,0,0);
    cout << "Flag inverse 3\n";
    
    dumMAr[0][0][0] = Er[0][0][0] - BCar[0][0][0];
    dumMAr[0][0][1] = Er[0][0][1] - BCar[0][0][1];
    dumMAr[0][1][0] = Er[0][1][0] - BCar[0][1][0];
    dumMAr[0][1][1] = Er[0][1][1] - BCar[0][1][1];
    
    dumMAi[0][0][0] = Ei[0][0][0] - BCai[0][0][0];
    dumMAi[0][0][1] = Ei[0][0][1] - BCai[0][0][1];
    dumMAi[0][1][0] = Ei[0][1][0] - BCai[0][1][0];
    dumMAi[0][1][1] = Ei[0][1][1] - BCai[0][1][1];
    
    CMatrixMult(dumMAr,dumMAi,dumMBr,dumMBi,BCcr,BCci,0,0,0);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    // Here we calculate v_{J-1}
    
    
    CMatrixMult(BCbr,BCbi,Dr,Di,dumMAr,dumMAi,0,0,0);
    
    dumMBr[0][0][0] = nur[0][0][0] - dumMAr[0][0][0];
    dumMBr[0][0][1] = nur[0][0][1] - dumMAr[0][0][1];
    dumMBr[0][1][0] = nur[0][1][0] - dumMAr[0][1][0];
    dumMBr[0][1][1] = nur[0][1][1] - dumMAr[0][1][1];
    
    dumMBi[0][0][0] = nui[0][0][0] - dumMAi[0][0][0];
    dumMBi[0][0][1] = nui[0][0][1] - dumMAi[0][0][1];
    dumMBi[0][1][0] = nui[0][1][0] - dumMAi[0][1][0];
    dumMBi[0][1][1] = nui[0][1][1] - dumMAi[0][1][1];
    
    CMatrixMult(BCcr,BCci,dumMBr,dumMBi,dumMAr,dumMAi,0,0,0);
    
    dumMBr[0][0][0] = Hr[0][0][0] - dumMAr[0][0][0];
    dumMBr[0][0][1] = Hr[0][0][1] - dumMAr[0][0][1];
    dumMBr[0][1][0] = Hr[0][1][0] - dumMAr[0][1][0];
    dumMBr[0][1][1] = Hr[0][1][1] - dumMAr[0][1][1];
    
    dumMBi[0][0][0] = Hi[0][0][0] - dumMAi[0][0][0];
    dumMBi[0][0][1] = Hi[0][0][1] - dumMAi[0][0][1];
    dumMBi[0][1][0] = Hi[0][1][0] - dumMAi[0][1][0];
    dumMBi[0][1][1] = Hi[0][1][1] - dumMAi[0][1][1];
    
    // This gives: dumMA  ==  [ H - ( E - F alpha_{J-2}^{-1} ) ( eta - mu C^{-1} A )^{-1} ( nu - mu C^{-1} D ) ]^{-1}  ==  [ H - BCc ( nu - BCb D ) ]^{-1}
    CMatrixInv(dumMBr,dumMBi,dumMCr,dumMCi,0,0);
    cout << "Flag inverse 4\n";
    
    cout << "\n\nIs this the problem matrix?\n";
    cout << dumMBr[0][0][0] << "\t" << dumMBr[0][0][1] << "\n";
    cout << dumMBr[0][1][0] << "\t" << dumMBr[0][1][1] << "\n\n";
    
    
    cout << dumMCr[0][0][0] << "\t" << dumMCr[0][0][1] << "\n";
    cout << dumMCr[0][1][0] << "\t" << dumMCr[0][1][1] << "\n\n";
    
    
    CVectorMult(BCar,BCai,gammar,gammai,dumVAr,dumVAi,0,J-2,0);
    
    CVectorMult(BCbr,BCbi,Mr,Mi,dumVCr,dumVCi,0,0,0);
    
    dumVBr[0][0][0] = xr[0][0][0] - dumVCr[0][0][0];
    dumVBr[0][1][0] = xr[0][1][0] - dumVCr[0][1][0];
    
    dumVBi[0][0][0] = xi[0][0][0] - dumVCi[0][0][0];
    dumVBi[0][1][0] = xi[0][1][0] - dumVCi[0][1][0];
    
    CVectorMult(BCcr,BCci,dumVBr,dumVBi,dumVCr,dumVCi,0,0,0);
    
    // This gives: dumVD  ==  N_{J-2} + F alpha_{J-2}^{-1} gamma_{J-2} - ( E - F alpha_{J-2}^{-1} ) ( eta - mu C^{-1} A )^{-1} ( x - mu C^{-1} M )  ==  [ N_{J-2} + BCa gamma_{J-2} - BCc ( x - BCb M_{J-2} ) ]
    dumVDr[0][0][0] = Nr[0][0][0] + dumVAr[0][0][0] - dumVCr[0][0][0];
    dumVDr[0][1][0] = Nr[0][1][0] + dumVAr[0][1][0] - dumVCr[0][1][0];
    
    dumVDi[0][0][0] = Ni[0][0][0] + dumVAi[0][0][0] - dumVCi[0][0][0];
    dumVDi[0][1][0] = Ni[0][1][0] + dumVAi[0][1][0] - dumVCi[0][1][0];
    
    
    
    // This step actually gets us v_{J-1}
    CVectorMult(dumMCr,dumMCi,dumVDr,dumVDi,vr,vi,0,0,J-1);
    
    cout << "Method II:\n" << "Got vr_{J-1} \t \n";
    
    /*    cout << "\nThis is vr_{J-1} \n";
     cout << vr[J-1][0][0] << "\n";
     cout << vr[J-1][1][0] << "\n\n";
     */
    
    
    
    // We use v_{J-1} to calculate u_{J-1} too
    // alpha_{i} v_{i}
    CVectorMult(alphar,alphai,vr,vi,dumVAr,dumVAi,J-1,J-1,0);
    // Gives u_{J-1}
    ur[J-1][0][0] = - dumVAr[0][0][0] - gammar[J-1][0][0];
    ur[J-1][1][0] = - dumVAr[0][1][0] - gammar[J-1][1][0];
    
    ui[J-1][0][0] = - dumVAi[0][0][0] - gammai[J-1][0][0];
    ui[J-1][1][0] = - dumVAi[0][1][0] - gammai[J-1][1][0];
    
    cout << "Got u_{J-1} ( J-1 = " << J-1 <<" ) \n";
    
    cout << "\nThis is ur_{J-1} \n";
    cout << ur[J-1][0][0] << "\n";
    cout << ur[J-1][1][0] << "\n\n";
    
    cout << "\nThis is vr_{J-1} \n";
    cout << vr[J-1][0][0] << "\n";
    cout << vr[J-1][1][0] << "\n\n";
    
    
    
    
    
    
    
    
    
    
    
    
    /*
     
     
     ********************************
     *                              *
     *                              *
     *                              *
     *                              *
     *          METHOD III          *
     *                              *
     *                              *
     *                              *
     *                              *
     ********************************
     
     
     
     
     The outer boundary conditions can be expressed in the form of a vector equation as:
     
     \eta u_{J-2} + \mu u_{J-1} + \nu v_{J-1} = x
     
     which is valid for the centre of the outermost cell, labelled J-1
     
     If this is used in coordination with the oscillation equations and the u-v relation for J-2, as:
     
     A_{J-2,J-1} u_{J-2} + C_{J-2,J-1} u_{J-1} + D_{J-2,J-1} v_{J-1} = M_{J-2,J-1}
     
     E_{J-2,J-1} u_{J-2} + F_{J-2,J-1} v_{J-2} + H_{J-2,J-1} v_{J-1} = N_{J-2,J-1}
     
     u_{J-1} + alpha_{J-1} v_{J-1} + gamma_{J-1} = 0
     
     we can eliminate everything except v_{J-1} by substituting the BC into the A equation, then using the J-1 recurrence relation to eliminate u_{J-2}, which gives us the (very long) expression:
     
     v_{J-1} = [ D - C alpha_{J-1} + A eta^{-1} ( mu alpha_{J-1} - nu ) ]^{-1} [ M + C gamma_{J-1} - A eta^{-1} ( x + mu gamma_{J-1} ) ]
     
     
     Introducing the following variables breaks down the expression somewhat:
     
     BCa = A eta^{-1}
     
     BCb = unnecessary
     
     BCc = unnecessary
     
     
     which gives
     
     v_{J-1} = [ D - C alpha_{J-1} + BCa ( mu alpha_{J-1} - nu ) ]^{-1} [ M + C gamma_{J-1} - BCa ( x + mu gamma_{J-1} ) ]
     
     
     */
    
    
    
    // Boundary condition simplification variables calculated here
    
    
    // BCa defined here
    
    CMatrixInv(etar,etai,dumMAr,dumMAi,0,0);
    cout << "Flag inverse 1\n";
    
    CMatrixMult(Ar,Ai,dumMAr,dumMAi,BCar,BCai,0,0,0);
    
    
    
    
    
    
    
    // Here we calculate v_{J-1}
    
    
    CMatrixMult(mur,mui,alphar,alphai,dumMAr,dumMAi,0,J-1,0);
    
    dumMBr[0][0][0] = - nur[0][0][0] + dumMAr[0][0][0];
    dumMBr[0][0][1] = - nur[0][0][1] + dumMAr[0][0][1];
    dumMBr[0][1][0] = - nur[0][1][0] + dumMAr[0][1][0];
    dumMBr[0][1][1] = - nur[0][1][1] + dumMAr[0][1][1];
    
    dumMBi[0][0][0] = - nui[0][0][0] + dumMAi[0][0][0];
    dumMBi[0][0][1] = - nui[0][0][1] + dumMAi[0][0][1];
    dumMBi[0][1][0] = - nui[0][1][0] + dumMAi[0][1][0];
    dumMBi[0][1][1] = - nui[0][1][1] + dumMAi[0][1][1];
    
    CMatrixMult(BCar,BCai,dumMBr,dumMBi,dumMAr,dumMAi,0,0,0);
    
    CMatrixMult(Cr,Ci,alphar,alphai,dumMBr,dumMBi,0,J-1,0);
    
    dumMCr[0][0][0] = Dr[0][0][0] - dumMBr[0][0][0] + dumMAr[0][0][0];
    dumMCr[0][0][1] = Dr[0][0][1] - dumMBr[0][0][1] + dumMAr[0][0][1];
    dumMCr[0][1][0] = Dr[0][1][0] - dumMBr[0][1][0] + dumMAr[0][1][0];
    dumMCr[0][1][1] = Dr[0][1][1] - dumMBr[0][1][1] + dumMAr[0][1][1];
    
    dumMCi[0][0][0] = Di[0][0][0] - dumMBi[0][0][0] + dumMAi[0][0][0];
    dumMCi[0][0][1] = Di[0][0][1] - dumMBi[0][0][1] + dumMAi[0][0][1];
    dumMCi[0][1][0] = Di[0][1][0] - dumMBi[0][1][0] + dumMAi[0][1][0];
    dumMCi[0][1][1] = Di[0][1][1] - dumMBi[0][1][1] + dumMAi[0][1][1];
    
    // This gives: dumMA  ==  [ D - C alpha_{J-1} + A eta^{-1} ( mu alpha_{J-1} - nu ) ]^{-1}  ==  [ D - C alpha_{J-1} + BCa ( mu alpha_{J-1} - nu ) ]^{-1}
    CMatrixInv(dumMCr,dumMCi,dumMAr,dumMAi,0,0);
    
    
    CVectorMult(mur,mui,gammar,gammai,dumVAr,dumVAi,0,J-1,0);
    
    dumVBr[0][0][0] = xr[0][0][0] + dumVAr[0][0][0];
    dumVBr[0][1][0] = xr[0][1][0] + dumVAr[0][1][0];
    
    dumVBi[0][0][0] = xi[0][0][0] + dumVAi[0][0][0];
    dumVBi[0][1][0] = xi[0][1][0] + dumVAi[0][1][0];
    
    CVectorMult(BCar,BCai,dumVBr,dumVBi,dumVCr,dumVCi,0,0,0);
    
    CVectorMult(Cr,Ci,gammar,gammai,dumVBr,dumVBi,0,J-1,0);
    
    // This gives: dumVA  ==  [ M + C gamma_{J-1} - A eta^{-1} ( x + mu gamma_{J-1} ) ]  ==  [ M + C gamma_{J-1} - BCa ( x + mu gamma_{J-1} ) ]
    dumVAr[0][0][0] = Mr[0][0][0] + dumVBr[0][0][0] - dumVCr[0][0][0];
    dumVAr[0][1][0] = Mr[0][1][0] + dumVBr[0][1][0] - dumVCr[0][1][0];
    
    dumVAi[0][0][0] = Mi[0][0][0] + dumVBi[0][0][0] - dumVCi[0][0][0];
    dumVAi[0][1][0] = Mi[0][1][0] + dumVBi[0][1][0] - dumVCi[0][1][0];
    
    
    
    // This step actually gets us v_{J-1}
    CVectorMult(dumMAr,dumMAi,dumVAr,dumVAi,vr,vi,0,0,J-1);
    
    cout << "Method III:\n" << "Got vr_{J-1} \t \n";
    
    /*    cout << "\nThis is vr_{J-1} \n";
     cout << vr[J-1][0][0] << "\n";
     cout << vr[J-1][1][0] << "\n\n";
     */
    
    
    
    // We use v_{J-1} to calculate u_{J-1} too
    // alpha_{i} v_{i}
    CVectorMult(alphar,alphai,vr,vi,dumVAr,dumVAi,J-1,J-1,0);
    // Gives u_{J-1}
    ur[J-1][0][0] = - dumVAr[0][0][0] - gammar[J-1][0][0];
    ur[J-1][1][0] = - dumVAr[0][1][0] - gammar[J-1][1][0];
    
    ui[J-1][0][0] = - dumVAi[0][0][0] - gammai[J-1][0][0];
    ui[J-1][1][0] = - dumVAi[0][1][0] - gammai[J-1][1][0];
    
    cout << "Got u_{J-1} ( J-1 = " << J-1 <<" ) \n";
    
    cout << "\nThis is ur_{J-1} \n";
    cout << ur[J-1][0][0] << "\n";
    cout << ur[J-1][1][0] << "\n\n";
    
    cout << "\nThis is vr_{J-1} \n";
    cout << vr[J-1][0][0] << "\n";
    cout << vr[J-1][1][0] << "\n\n";
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*
     
     
     ********************************
     *                              *
     *                              *
     *                              *
     *                              *
     *          METHOD IV           *
     *                              *
     *                              *
     *                              *
     *                              *
     ********************************
     
     
     
     
     The outer boundary conditions can be expressed in the form of a vector equation as:
     
     \eta u_{J-2} + \mu u_{J-1} + \nu v_{J-1} = x
     
     which is valid for the centre of the outermost cell, labelled J-1
     
     If this is used in coordination with the oscillation equations and the u-v relation for J-2, as:
     
     A_{J-2,J-1} u_{J-2} + C_{J-2,J-1} u_{J-1} + D_{J-2,J-1} v_{J-1} = M_{J-2,J-1}
     
     E_{J-2,J-1} u_{J-2} + F_{J-2,J-1} v_{J-2} + H_{J-2,J-1} v_{J-1} = N_{J-2,J-1}
     
     u_{J-1} + alpha_{J-1} v_{J-1} + gamma_{J-1} = 0
     
     we can eliminate everything except v_{J-1} by substituting the A equation into the BC, then using the J-1 recurrence relation to eliminate u_{J-2}, which gives us the (very long) expression:
     
     v_{J-1} = [ nu - mu alpha_{J-1} + eta A^{-1} ( C alpha_{J-1} - D ) ]^{-1} [ x + mu gamma_{J-1} - eta A^{-1} ( M + C gamma_{J-1} ) ]
     
     
     Introducing the following variables breaks down the expression somewhat:
     
     BCa = eta A^{-1}
     
     BCb = unnecessary
     
     BCc = unnecessary
     
     
     which gives
     
     v_{J-1} = [ D - C alpha_{J-1} + BCa ( mu alpha_{J-1} - nu ) ]^{-1} [ M + C gamma_{J-1} - BCa ( x + mu gamma_{J-1} ) ]
     
     
     */
    
    
    
    // Boundary condition simplification variables calculated here
    
    
    // BCa defined here
    
    CMatrixInv(Ar,Ai,dumMAr,dumMAi,0,0);
    cout << "Flag inverse 1\n";
    
    CMatrixMult(etar,etai,dumMAr,dumMAi,BCar,BCai,0,0,0);
    
    
    
    
    
    
    
    // Here we calculate v_{J-1}
    
    
    CMatrixMult(Cr,Ci,alphar,alphai,dumMAr,dumMAi,0,J-1,0);
    
    dumMBr[0][0][0] = - Dr[0][0][0] + dumMAr[0][0][0];
    dumMBr[0][0][1] = - Dr[0][0][1] + dumMAr[0][0][1];
    dumMBr[0][1][0] = - Dr[0][1][0] + dumMAr[0][1][0];
    dumMBr[0][1][1] = - Dr[0][1][1] + dumMAr[0][1][1];
    
    dumMBi[0][0][0] = - Di[0][0][0] + dumMAi[0][0][0];
    dumMBi[0][0][1] = - Di[0][0][1] + dumMAi[0][0][1];
    dumMBi[0][1][0] = - Di[0][1][0] + dumMAi[0][1][0];
    dumMBi[0][1][1] = - Di[0][1][1] + dumMAi[0][1][1];
    
    CMatrixMult(BCar,BCai,dumMBr,dumMBi,dumMAr,dumMAi,0,0,0);
    
    CMatrixMult(mur,mui,alphar,alphai,dumMBr,dumMBi,0,J-1,0);
    
    dumMCr[0][0][0] = nur[0][0][0] - dumMBr[0][0][0] + dumMAr[0][0][0];
    dumMCr[0][0][1] = nur[0][0][1] - dumMBr[0][0][1] + dumMAr[0][0][1];
    dumMCr[0][1][0] = nur[0][1][0] - dumMBr[0][1][0] + dumMAr[0][1][0];
    dumMCr[0][1][1] = nur[0][1][1] - dumMBr[0][1][1] + dumMAr[0][1][1];
    
    dumMCi[0][0][0] = nui[0][0][0] - dumMBi[0][0][0] + dumMAi[0][0][0];
    dumMCi[0][0][1] = nui[0][0][1] - dumMBi[0][0][1] + dumMAi[0][0][1];
    dumMCi[0][1][0] = nui[0][1][0] - dumMBi[0][1][0] + dumMAi[0][1][0];
    dumMCi[0][1][1] = nui[0][1][1] - dumMBi[0][1][1] + dumMAi[0][1][1];
    
    // This gives: dumMA  ==  [ nu - mu alpha_{J-1} + eta A^{-1} ( C alpha_{J-1} - D ) ]^{-1}  ==  [ D - C alpha_{J-1} + BCa ( mu alpha_{J-1} - nu ) ]^{-1}
    CMatrixInv(dumMCr,dumMCi,dumMAr,dumMAi,0,0);
    
    
    CVectorMult(Cr,Ci,gammar,gammai,dumVAr,dumVAi,0,J-1,0);
    
    dumVBr[0][0][0] = Mr[0][0][0] + dumVAr[0][0][0];
    dumVBr[0][1][0] = Mr[0][1][0] + dumVAr[0][1][0];
    
    dumVBi[0][0][0] = Mi[0][0][0] + dumVAi[0][0][0];
    dumVBi[0][1][0] = Mi[0][1][0] + dumVAi[0][1][0];
    
    CVectorMult(BCar,BCai,dumVBr,dumVBi,dumVCr,dumVCi,0,0,0);
    
    CVectorMult(mur,mui,gammar,gammai,dumVBr,dumVBi,0,J-1,0);
    
    // This gives: dumVA  ==  [ x + mu gamma_{J-1} - eta A^{-1} ( M + C gamma_{J-1} ) ]  ==  [ M + C gamma_{J-1} - BCa ( x + mu gamma_{J-1} ) ]
    dumVAr[0][0][0] = xr[0][0][0] + dumVBr[0][0][0] - dumVCr[0][0][0];
    dumVAr[0][1][0] = xr[0][1][0] + dumVBr[0][1][0] - dumVCr[0][1][0];
    
    dumVAi[0][0][0] = xi[0][0][0] + dumVBi[0][0][0] - dumVCi[0][0][0];
    dumVAi[0][1][0] = xi[0][1][0] + dumVBi[0][1][0] - dumVCi[0][1][0];
    
    
    
    // This step actually gets us v_{J-1}
    CVectorMult(dumMAr,dumMAi,dumVAr,dumVAi,vr,vi,0,0,J-1);
    
    
    cout << "Method IV:\n" << "Got vr_{J-1} \t \n";
    
    /*    cout << "\nThis is vr_{J-1} \n";
     cout << vr[J-1][0][0] << "\n";
     cout << vr[J-1][1][0] << "\n\n";
     */
    
    
    
    // We use v_{J-1} to calculate u_{J-1} too
    // alpha_{i} v_{i}
    CVectorMult(alphar,alphai,vr,vi,dumVAr,dumVAi,J-1,J-1,0);
    // Gives u_{J-1}
    ur[J-1][0][0] = - dumVAr[0][0][0] - gammar[J-1][0][0];
    ur[J-1][1][0] = - dumVAr[0][1][0] - gammar[J-1][1][0];
    
    ui[J-1][0][0] = - dumVAi[0][0][0] - gammai[J-1][0][0];
    ui[J-1][1][0] = - dumVAi[0][1][0] - gammai[J-1][1][0];
    
    cout << "Got u_{J-1} ( J-1 = " << J-1 <<" ) \n";
    
    cout << "\nThis is ur_{J-1} \n";
    cout << ur[J-1][0][0] << "\n";
    cout << ur[J-1][1][0] << "\n\n";
    
    cout << "\nThis is vr_{J-1} \n";
    cout << vr[J-1][0][0] << "\n";
    cout << vr[J-1][1][0] << "\n\n";
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    // Now that all of the alphas and gammas have been produced, we can evaluate u_{i} from u_{i+1} and v_{i+1}.  This means that we will need to use u_{i} with alpha and gamma to get v_{i} before we move onto the next iteration (as in, it's got to happen within the same loop).
    
    // The equation we use is: u_{i} = RECu u_{i+1} + RECv v_{i+1} + RECc
    
    
    //u_{i}
    for (k=J-2; k >= 0; k = k - 1) {
        
        CVectorMult(RECur,RECui,ur,ui,dumVAr,dumVAi,k,k+1,0);
        
        CVectorMult(RECvr,RECvi,vr,vi,dumVBr,dumVBi,k,k+1,0);
        
        vr[k][0][0] = dumVAr[0][0][0] + dumVBr[0][0][0] + RECcr[k][0][0];
        vr[k][1][0] = dumVAr[0][1][0] + dumVBr[0][1][0] + RECcr[k][1][0];
        
        vi[k][0][0] = dumVAi[0][0][0] + dumVBi[0][0][0] + RECci[k][0][0];
        vi[k][1][0] = dumVAi[0][1][0] + dumVBi[0][1][0] + RECci[k][1][0];
        
        //        cout << "vr[" << k << "][0][0] = " << vr[k][0][0] << "\n";
        //        cout << "vr[" << k << "][1][0] = " << vr[k][1][0] << "\n";
        
        //        cout << "RECu = " << RECur[k][0][0] << "\t" << RECur[k][0][1] << "\t" << RECur[k][1][0] << "\t" << RECur[k][1][1] << "\n";
        
        //        cout << "RECv = " << RECvr[k][0][0] << "\t" << RECvr[k][0][1] << "\t" << RECvr[k][1][0] << "\t" << RECvr[k][1][1] << "\n";
        
        //        cout << "dumVAr = " << dumVAr[0][0][0] << "\n";
        //        cout << "dumVAr = " << dumVAr[0][1][0] << "\n";
        
        //        cout << "dumVBr = " << dumVBr[0][0][0] << "\n";
        //        cout << "dumVBr = " << dumVBr[0][1][0] << "\n\n";
        
        
        CVectorMult(alphar,alphai,vr,vi,dumVAr,dumVAi,k,k,0);
        
        ur[k][0][0] = - dumVAr[0][0][0] - gammar[k][0][0];
        ur[k][1][0] = - dumVAr[0][1][0] - gammar[k][1][0];
        
        ui[k][0][0] = - dumVAi[0][0][0] - gammai[k][0][0];
        ui[k][1][0] = - dumVAi[0][1][0] - gammai[k][1][0];
        
        
        /*
         //Now need to follow up from here and then do use the gammas and alphas to get the v for the u that we've just calculated.
         
         // alpha_{i}^{-1}
         CMatrixInv(alphar,alphai,dumMAr,dumMAi,k,0);
         
         // - ( gamma_{i} + u_{i} )
         dumVCr[0][0][0] = - gammar[k][0][0] - ur[k][0][0];
         dumVCr[0][1][0] = - gammar[k][1][0] - ur[k][1][0];
         
         dumVCi[0][0][0] = - gammai[k][0][0] - ui[k][0][0];
         dumVCi[0][1][0] = - gammai[k][1][0] - ui[k][1][0];
         
         // Gives v_{i}
         CVectorMult(dumMAr,dumMAi,dumVCr,dumVCi,vr,vi,0,0,k);
         */
        
        
        
        
        //        cout << "\t Done the v_{" << k << "} and u_{" << k << "} loop \n";
        
        
        
    }
    
    
    
    
    
    
    // Now the rescaling must be undone
    
    for (k=0; k<J; k=k+1) {
        
        // This rescales a
        CompMult(&T_ar, &T_ai, &ur[k][0][0], &ui[k][0][0], &dummyr, &dummyi);
        ur[k][0][0] = dummyr;
        ui[k][0][0] = dummyi;
        
        // This rescales b
        CompMult(&T_br, &T_bi, &ur[k][1][0], &ui[k][1][0], &dummyr, &dummyi);
        ur[k][1][0] = dummyr;
        ui[k][1][0] = dummyi;
        
        // This rescales c
        CompMult(&T_cr, &T_ci, &vr[k][0][0], &vi[k][0][0], &dummyr, &dummyi);
        vr[k][0][0] = dummyr;
        vi[k][0][0] = dummyi;
        
        // This rescales d
        CompMult(&T_dr, &T_di, &vr[k][1][0], &vi[k][1][0], &dummyr, &dummyi);
        vr[k][1][0] = dummyr;
        vi[k][1][0] = dummyi;
        
        
        
    }
    
    
    
    
    
    
    
    
    
    
    // This bit opens the file to write the data into
    ofstream outfile;
    outfile.open("Output/Henyey_BPT.dat", ios::out);
    
    // This sets the precision at which values are printed at to the named file output
    outfile.precision(10);
    
    // This is for the Pfahl comparison
    ofstream comparison_file;
    comparison_file.open("Output/comparison_file.dat", ios::out);
    
    comparison_file.precision(10);
    
    double xi_h_r, xi_h_pprime_part, xi_r_eq, H_rho, H_p, mod_xi_radial, xi_h_real, xi_h_imaginary, mod_xi_h, delta_P_r, delta_P_i, mod_delta_P, delta_P_r_old, delta_P_i_old, delta_P_r_next, delta_P_i_next;
    double ddelta_P_dr_r, ddelta_P_dr_i, V_div_xi_r_r, V_div_xi_r_i, num_r, num_i, denom_r, denom_i, dgrr_dr, otherVdiv_r, otherVdiv_i, pprime_comp_r, pprime_comp_i;
    double Gvar_r, Gvar_i, Hvar_r, Hvar_i, gradient, pprime_comp_second_r, pprime_comp_second_i, rhoprime_r, rhoprime_i, dp0dr, dplus, dminus, dp0dr_old, dp0dr_avg, dp0dr_sum;
    double delta_P_r_new, delta_P_i_new, xi_h_over_xi_r_real, xi_h_over_xi_r_im, log_mod_xi_r, log_mod_xi_h, d2p_dr2, dp_dr, d2p0_dr2, dp0dr_plus, dp0dr_minus;
    double dDeltaP_dr_b_r, dDeltaP_dr_b_i, dpprime_dr_b_r, dpprime_dr_b_i, dxi_r_dr_b_r, dxi_r_dr_b_i, dp0_dr_b, d2p0_dr2_b, V_div_xi_r_b_r, V_div_xi_r_b_i;
    double xi_r_analytic_r, xi_r_analytic_i, D_analytic, V_analytic_r, V_analytic_i, Fprime_analytic_r, Fprime_analytic_i, dTprime_dT_r, dTprime_dT_i, dT0_dr, dxi_r_dr_r, dxi_r_dr_i, deltaP_hydro_r, deltaP_hydro_i, deltaP_hydro_next_r, deltaP_hydro_next_i, ddeltaP_dr_hydro_r, ddeltaP_dr_hydro_i, xi_r_analytic_hydro_r, xi_r_analytic_hydro_i, V_analytic_hydro_r, V_analytic_hydro_i, xi_r_analytic_hydro_next_r, xi_r_analytic_hydro_next_i, V_analytic_hydro_next_r, V_analytic_hydro_next_i, D_analytic_next;
    double ddeltaP_dr_hydro_Jminus3_r, ddeltaP_dr_hydro_Jminus3_i, ddeltaP_dr_hydro_Jminus4_r, ddeltaP_dr_hydro_Jminus4_i, V_cont_r, V_cont_i, s_prime_r, s_prime_i, delta_s_r, delta_s_i, F_conv_prime_div_F_conv_r, F_conv_prime_div_F_conv_i,V_rho_r,V_rho_i;
    double ds_prime_ds0_r, ds_prime_ds0_i, dr_dT_Fc_K0, x_conv, s_prime_JP, xi_r_JP, V_cont_r_old, V_analytic_hydro_old_r, xi_r_JP_term, s_prime_equation_JP;
    
    int k_start_JP;
    
    s_prime_JP = 0.0;
    k_start_JP = 0;
    
    //
    //
    //
    //
    //
    //
    //
    //
    // This is the loop for output
    //
    //
    //
    //
    //
    //
    //
    //
    //
    //
    
    kold = 0;
    
    for (k=0; k<J; k=k+1) {
        
        xi_h_r = m*omega*((mp + Mstar)/(100.0*mp)) * (1.0/(rmid_cm_HR_output[k] * m * m * omega * omega)) * ( (pressure_HR_output[k] * vr[k][0][0] / rho_HR_output[k]) + (f * rmid_cm_HR_output[k] * rmid_cm_HR_output[k]) );
        
        xi_h_pprime_part = m*omega*((mp + Mstar)/(100.0*mp)) * (1.0/(rmid_cm_HR_output[k] * m * m * omega * omega)) * ( (pressure_HR_output[k] * vr[k][0][0] / rho_HR_output[k]) );
        
        
        // Here the value for kold is calculated
        if ( rmid_cm_HR_output[k] < 0.5*(rmid_cm[0]+rmid_cm[1]) ) {
            // We need to use the 0th, 1st and 2nd old cells, so set kold = 1
            kold = 1;
            
        } else {
            
            if (rmid_cm_HR_output[k] > 0.5*(rmid_cm[Jold-2]+rmid_cm[Jold-1]) ) {
                // We need to use the Jold-3, Jold-2 and Jold-1 cells, so set kold = Jold-2
                kold = Jold - 2;
                
            } else {
                
                if ( rmid_cm_HR_output[k] > 0.5*(rmid_cm[kold-1]+rmid_cm[kold]) && rmid_cm_HR_output[k] < 0.5*(rmid_cm[kold]+rmid_cm[kold+1]) ) {
                    // If our current value for kold works, we keep it
                    kold = kold;
                    
                } else {
                    
                    for (kold=0; rmid_cm_HR_output[k] > 0.5*(rmid_cm[kold]+rmid_cm[kold+1]); kold = kold + 1 ){
                        // Don't need anything in here, as it is all covered by the for loop itself
                    }
                    
                }
                
                
                
            }
            
            
        }
        
        // Checking the value of kold that has been chosen
        
        if ( rmid_cm_HR_output[k] > 0.5*(rmid_cm[kold-1]+rmid_cm[kold]) && rmid_cm_HR_output[k] < 0.5*(rmid_cm[kold]+rmid_cm[kold+1]) ) {
            // Everything is okay here
        } else {
            cout << "Oh dear! Problem with the value of kold for k = " << k << ", at rmid_cm = " << rmid_cm_HR_output[k] << "\n";
            
        }
        
        
        if (k % 10000000000 == 0) {
            cout << 0.5*(rmid_cm[kold-1]+rmid_cm[kold])/R << "\t\t" << rmid_cm_HR_output[k]/R << "\t\t" << 0.5*(rmid_cm[kold]+rmid_cm[kold+1])/R << "\n";
            
            cout << "rmid_cm_HR_output[k] = " << rmid_cm_HR_output[k] << "\n";
            cout << "rmid_cm_[kold-1] = " << rmid_cm[kold-1] << "\n";
            cout << "rmid_cm[kold] = " << rmid_cm[kold] << "\n";
            cout << "rmid_cm[kold+1] = " << rmid_cm[kold+1] << "\n";
            cout << "pressure[kold-1] = " << pressure[kold-1] << "\n";
            cout << "pressure[kold] = " << pressure[kold] << "\n";
            cout << "pressure[kold+1] = " << pressure[kold+1] << "\n";
        }
        
        // Here the first and second derivatives of the background pressure are calculated, I reuse the delta[3] array from earlier
        
        delta[0] = rmid_cm_HR_output[k] - rmid_cm[kold-1];
        delta[1] = rmid_cm_HR_output[k] - rmid_cm[kold];
        delta[2] = rmid_cm[kold+1] - rmid_cm_HR_output[k];
        
        dp_dr = pressure[kold+1]*(delta[1] - delta[0])/( (delta[2]+delta[0])*(delta[1]-delta[2]) )   +   pressure[kold]*(delta[0]-delta[2])/( (delta[1]-delta[2])/(delta[1]+delta[0]) )   -   pressure[kold-1]*(delta[1]+delta[2])/( (delta[2]+delta[0])*(delta[1]+delta[0]) );
        
        d2p_dr2 = -2.0*pressure[kold+1]/( (delta[1]-delta[2])*(delta[0]+delta[2]) )   +   2.0*pressure[kold]/( (delta[1]-delta[2])*(delta[1]+delta[0]) )   +   2.0*pressure[kold-1]/( (delta[2]+delta[0])*(delta[1]+delta[0]) );
        
        
        d2p0_dr2 = ((2.0)/( (delta[2]+delta[0])*delta[2]*delta[0] ))*(  pressure[kold+1]*delta[0]  -  pressure_HR_output[k]*(delta[2]+delta[0]) + pressure[kold-1]*delta[2]  );
        
        
        
        
        
        // These are the variables which I am using for this version of getting dDeltaP_dr_b_r etc, because I want to see if I can get it to be a little bit smoother
        // dDeltaP_dr_b_r, dDeltaP_dr_b_i, dpprime_dr_b_r, dpprime_dr_b_i, dxi_r_dr_r, dxi_r_dr_b_i, dp0_dr_b, d2p0_dr2_b
        
        if (k == 0) {
            
            dpprime_dr_b_r = (pressure_HR_output[1]*vr[1][0][0] - pressure_HR_output[0]*vr[0][0][0])/(rmid_cm_HR_output[1] - rmid_cm_HR_output[0]);
            
            dpprime_dr_b_i = (pressure_HR_output[1]*vi[1][0][0] - pressure_HR_output[0]*vi[0][0][0])/(rmid_cm_HR_output[1] - rmid_cm_HR_output[0]);
            
            dxi_r_dr_b_r = R*(ur[k+1][0][0] - ur[k][0][0])/(radius_cm_HR_output[k+1] - radius_cm_HR_output[k]);
            
            dxi_r_dr_b_i = R*(ui[k+1][0][0] - ui[k][0][0])/(radius_cm_HR_output[k+1] - radius_cm_HR_output[k]);
            
            dp0_dr_b = -0.5*( grav_HR_output[k]*rho_HR_output[k] + grav_HR_output[k+1]*rho_HR_output[k+1] ); //(pressure_HR_output[1] - pressure_HR_output[0])/(rmid_cm_HR_output[1] - rmid_cm_HR_output[0]);
            
            // Just reuse this from earlier, at least for the moment. See if it's good enough and then maybe upgrade it later
            d2p0_dr2_b = d2p0_dr2;
            
            // We don't include xi_r in this one, as it's at the centre, xo xi_r should = 0
            dDeltaP_dr_b_r = dpprime_dr_b_r  +  ( dxi_r_dr_b_r * dp0_dr_b ); // + ( R*ur[k][0][0]*d2p0_dr2_b );
            
            dDeltaP_dr_b_r = dpprime_dr_b_i  +  ( dxi_r_dr_b_i * dp0_dr_b ); // + ( R*ui[k][0][0]*d2p0_dr2_b );
            
        } else {
            
            if (k == J-1) {
                
                dpprime_dr_b_r = (pressure_HR_output[J-1]*vr[J-1][0][0] - pressure_HR_output[J-2]*vr[J-2][0][0])/(rmid_cm_HR_output[J-1] - rmid_cm_HR_output[J-2]);
                
                dpprime_dr_b_i = (pressure_HR_output[J-1]*vi[J-1][0][0] - pressure_HR_output[J-2]*vi[J-2][0][0])/(rmid_cm_HR_output[J-1] - rmid_cm_HR_output[J-2]);
                
                dxi_r_dr_b_r = R*(ur[k][0][0] - ur[k-1][0][0])/(radius_cm_HR_output[k] - radius_cm_HR_output[k-1]);
                
                dxi_r_dr_b_i = R*(ui[k][0][0] - ui[k-1][0][0])/(radius_cm_HR_output[k] - radius_cm_HR_output[k-1]);
                
                dp0_dr_b = -0.5*( grav_HR_output[k-1]*rho_HR_output[k-1] + grav_HR_output[k]*rho_HR_output[k] ); //(pressure_HR_output[k] - pressure_HR_output[k-1])/(rmid_cm_HR_output[k] - rmid_cm_HR_output[k-1]);
                
                // Just reuse this from earlier, at least for the moment. See if it's good enough and then maybe upgrade it later
                d2p0_dr2_b = d2p0_dr2;
                
                
                
            } else {
                
                delta[0] = rmid_cm_HR_output[k] - rmid_cm_HR_output[k-1];
                delta[1] = rmid_cm_HR_output[k] - rmid_cm_HR_output[k];
                delta[2] = rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k];
                
                dpprime_dr_b_r = pressure_HR_output[k+1]*vr[k+1][0][0]*delta[0]/(delta[2]*( delta[0] + delta[2] ))   +   pressure_HR_output[k]*vr[k][0][0]*(delta[2] - delta[0])/(delta[0]*delta[2])   -   pressure_HR_output[k-1]*vr[k-1][0][0]*delta[2]/(delta[0]*( delta[0]+delta[2] ));
                
                dpprime_dr_b_i = pressure_HR_output[k+1]*vi[k+1][0][0]*delta[0]/(delta[2]*( delta[0] + delta[2] ))   +   pressure_HR_output[k]*vi[k][0][0]*(delta[2] - delta[0])/(delta[0]*delta[2])   -   pressure_HR_output[k-1]*vi[k-1][0][0]*delta[2]/(delta[0]*( delta[0]+delta[2] ));
                
                dxi_r_dr_b_r = R*(ur[k][0][0] - ur[k-1][0][0])/(radius_cm_HR_output[k] - radius_cm_HR_output[k-1]);
                
                dxi_r_dr_b_i = R*(ui[k][0][0] - ui[k-1][0][0])/(radius_cm_HR_output[k] - radius_cm_HR_output[k-1]);
                
                dp0_dr_b = -0.5*( grav_HR_output[k]*rho_HR_output[k] + grav_HR_output[k+1]*rho_HR_output[k+1] ); //-0.5*( grav_HR_output[k]*rho_HR_output[k] + grav_HR_output[k+1]*rho_HR_output[k+1] ); //pressure_HR_output[k+1]*delta[0]/(delta[2]*( delta[0] + delta[2] ))   +   pressure_HR_output[k]*(delta[2] - delta[0])/(delta[0]*delta[2])   -   pressure_HR_output[k-1]*delta[2]/(delta[0]*( delta[0]+delta[2] ));
                
                // Just reuse this from earlier, at least for the moment. See if it's good enough and then maybe upgrade it later
                d2p0_dr2_b = d2p0_dr2;
                
            }
            
            
            // Here the cases for k != 0 are all covered
            dDeltaP_dr_b_r = dpprime_dr_b_r  +  ( dxi_r_dr_b_r * dp0_dr_b )  +  ( 0.5*R*(ur[k][0][0] + ur[k-1][0][0])*d2p0_dr2_b );
            
            dDeltaP_dr_b_i = dpprime_dr_b_i  +  ( dxi_r_dr_b_i * dp0_dr_b )  +  ( 0.5*R*(ui[k][0][0] + ui[k-1][0][0])*d2p0_dr2_b );
            
            
            
            
        }
        
        
        
        // The above stuff is used to give V_div_xi_r later on, once the derivatice of g*r^2 has been calculated
        
        
        
        
        
        
        
        // This saves the delta_P values from the previous step with the suffix _old
        delta_P_r_old = delta_P_r;
        
        delta_P_i_old = delta_P_i;
        
        
        
        if (k==J-1) {
            
            // We avoid surface problems here (requires k=J, which isn't a thing)
            
            H_p = H_p;
            
            H_rho = H_rho;
            
            delta_P_r = delta_P_r;
            
            delta_P_i = delta_P_i;
            
            
        } else {
            
            if (k==0) {
                
                // For when k-1 is not an option
                
                // gradient of pressure
                dp0dr = (pressure_HR_output[k+1] - pressure_HR_output[k])/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k]);
                
                dp0dr_old = dp0dr;
                
                dp0dr_avg = dp0dr;
                
                // H_p is the pressure scale height
                H_p = - 0.5 * ( pressure_HR_output[k+1] + pressure_HR_output[k] ) * ( rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k] ) / ( pressure_HR_output[k+1] - pressure_HR_output[k] );
                
                // H_rho is the density scale height
                H_rho = - 0.5 * ( rho_HR_output[k+1] + rho_HR_output[k] ) * ( rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k] ) / ( rho_HR_output[k+1] - rho_HR_output[k] );
                
                // Lagrangian perturbation to pressure (real part)
                delta_P_r = (pressure_HR_output[k]*vr[k][0][0]) + (R*ur[k][0][0])*( (pressure_HR_output[k+1] - pressure_HR_output[k])/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k]) );
                
                // Lagrangian perturbation to pressure (imaginary part)
                delta_P_i = (pressure_HR_output[k]*vi[k][0][0]) + (R*ui[k][0][0])*( (pressure_HR_output[k+1] - pressure_HR_output[k])/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k]) );
                
                
                
            } else {
                
                
                
                
                // For all case which can have k+1 and k-1 and still be all okay
                
                
                
                
                if (k > 21 && k < J-22) {
                    
                    dp0dr_sum = 0.0;
                    
                    
                    for (i=k-20; i<=k+20; i=i+1) {
                        
                        // dplus and dminus are used in caluclating the pressure gradient
                        dplus = rmid_cm_HR_output[i+1] - rmid_cm_HR_output[i];
                        dminus = rmid_cm_HR_output[i] - rmid_cm_HR_output[i-1];
                        
                        // pressure gradient
                        dp0dr = pressure_HR_output[i+1]*( dminus / (dplus*(dplus+dminus)) )  +  pressure_HR_output[i]*( (1.0/dminus) - (1.0/dplus) )  -  pressure_HR_output[i-1] * ( dplus / (dminus*(dplus+dminus)) );
                        
                        dp0dr_sum = dp0dr_sum + dp0dr;
                        
                    }
                    
                    dp0dr_avg = dp0dr_sum / 41.0;
                    
                } else {
                    
                    // dplus and dminus are used in caluclating the pressure gradient
                    dplus = rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k];
                    dminus = rmid_cm_HR_output[k] - rmid_cm_HR_output[k-1];
                    
                    // pressure gradient
                    dp0dr = pressure_HR_output[k+1]*( dminus / (dplus*(dplus+dminus)) )  +  pressure_HR_output[k]*( (1.0/dminus) - (1.0/dplus) )  -  pressure_HR_output[k-1] * ( dplus / (dminus*(dplus+dminus)) );
                    
                    dp0dr_avg = dp0dr;
                }
                
                // dplus and dminus are used in caluclating the pressure gradient
                dplus = rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k];
                dminus = rmid_cm_HR_output[k] - rmid_cm_HR_output[k-1];
                
                // pressure gradient
                dp0dr = pressure_HR_output[k+1]*( dminus / (dplus*(dplus+dminus)) )  +  pressure_HR_output[k]*( (1.0/dminus) - (1.0/dplus) )  -  pressure_HR_output[k-1] * ( dplus / (dminus*(dplus+dminus)) );
                
                dp0dr_old = (pressure_HR_output[k+1] - pressure_HR_output[k-1])/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k-1]);
                
                // H_p is the pressure scale height
                H_p = - 0.5 * ( pressure_HR_output[k+1] + pressure_HR_output[k-1] ) * ( rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k-1] ) / ( pressure_HR_output[k+1] - pressure_HR_output[k-1] );
                
                // H_rho is the density scale height
                H_rho = - 0.5 * ( rho_HR_output[k+1] + rho_HR_output[k-1] ) * ( rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k-1] ) / ( rho_HR_output[k+1] - rho_HR_output[k-1] );
                
                // Lagrangian perturbation to pressure (real part)
                delta_P_r = (pressure_HR_output[k]*vr[k][0][0]) + (R*0.5*(ur[k][0][0] + ur[k-1][0][0]))*( (pressure_HR_output[k+1] - pressure_HR_output[k-1])/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k-1]) );
                
                delta_P_r_new = (pressure_HR_output[k]*vr[k][0][0]) + (R*0.5*(ur[k][0][0] + ur[k-1][0][0]))*dp0dr_avg;
                
                // Lagrangian perturbation to pressure (imaginary part)
                delta_P_i = (pressure_HR_output[k]*vi[k][0][0]) + (R*0.5*(ui[k][0][0] + ui[k-1][0][0]))*( (pressure_HR_output[k+1] - pressure_HR_output[k-1])/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k-1]) );
                
                delta_P_i_new = (pressure_HR_output[k]*vi[k][0][0]) + (R*0.5*(ui[k][0][0] + ui[k-1][0][0]))*dp0dr_avg;
                
                
            }
            
        }
        
        
        
        // This calculates delta_P for k+1, so that we can get the gradient of it
        i = k+1;
        
        if (i==J-1) {
            
            // We avoid surface problems here (requires k=J, which isn't a thing)
            
            
            delta_P_r_next = delta_P_r;
            
            delta_P_i_next = delta_P_i;
            
            
        } else {
            
            
            // For all case which can have i+1 and i-1 and still be all okay
            
            // Lagrangian perturbation to pressure (real part)
            delta_P_r_next = (pressure_HR_output[i]*vr[i][0][0]) + (R*0.5*(ur[i][0][0] + ur[i-1][0][0]))*( (pressure_HR_output[i+1] - pressure_HR_output[i-1])/(rmid_cm_HR_output[i+1] - rmid_cm_HR_output[i-1]) );
            
            // Lagrangian perturbation to pressure (imaginary part)
            delta_P_i_next = (pressure_HR_output[i]*vi[i][0][0]) + (R*0.5*(ui[i][0][0] + ui[i-1][0][0]))*( (pressure_HR_output[i+1] - pressure_HR_output[i-1])/(rmid_cm_HR_output[i+1] - rmid_cm_HR_output[i-1]) );
            
            
        }
        
        
        
        // This calculates the gradient of delta_P and g/r^2 (issues will arise for k=J-1 and probably for k=J-1 as well, though, so bear that in mind when you plot it)
        
        if (k==J-1 | k==J-2) {
            
            // If k==J-1, we leave them alone, as we can't get a gradient
            
            ddelta_P_dr_r = ddelta_P_dr_r;
            
            ddelta_P_dr_i = ddelta_P_dr_i;
            
            dgrr_dr = dgrr_dr;
            
        } else {
            
            ddelta_P_dr_r = (delta_P_r_next - delta_P_r)/(rmid_cm_HR[1] - rmid_cm_HR[0]);
            
            ddelta_P_dr_i = (delta_P_i_next - delta_P_i)/(rmid_cm_HR[1] - rmid_cm_HR[0]);
            
            if (k==0) {
                
                dgrr_dr = ( (grav_HR_output[k+1]/(rmid_cm_HR_output[k+1]*rmid_cm_HR_output[k+1])) - (grav_HR_output[k]/(rmid_cm_HR_output[k]*rmid_cm_HR_output[k])) )/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k]);
                
            } else {
                
                dgrr_dr = ( (grav_HR_output[k+1]/(rmid_cm_HR_output[k+1]*rmid_cm_HR_output[k+1])) - (grav_HR_output[k-1]/(rmid_cm_HR_output[k-1]*rmid_cm_HR_output[k-1])) )/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k-1]);
                
            }
            
        }
        
        
        
        // Here we get the numerator and denomiator in real and imaginary parts for the caluclation of V_div_xi_r
        
        num_r = (grav_HR_output[k]/(rmid_cm_HR_output[k]*m*m*omega*omega)) * (   -ddelta_P_dr_r   +   (rho_HR_output[k]*f*rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*( -(2.0/rmid_cm_HR_output[k]) + (rmid_cm_HR_output[k]*rmid_cm_HR_output[k]/grav_HR_output[k])*dgrr_dr ) )  );
        
        num_i = (grav_HR_output[k]/(rmid_cm_HR_output[k]*m*m*omega*omega)) * (   -ddelta_P_dr_i   );
        
        denom_r = - (  ddelta_P_dr_r + 2.0*rho_HR_output[k]*f*rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*( (1.0/rmid_cm_HR_output[k]) + ((3.0*grav_HR_output[k])/(rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*m*m*omega*omega)) )  );
        
        denom_i = - ddelta_P_dr_i;
        
        V_div_xi_r_r = (  (num_r*denom_r)  +  (num_i*denom_i)  )/(  (denom_r*denom_r)  +  (denom_i*denom_i)  );
        
        V_div_xi_r_i = (  num_i*denom_r  -  num_r*denom_i  )/(  denom_r*denom_r  +  denom_i*denom_i  );
        
        
        
        
        // Here we do it for the other attempt at getting the derivative of DeltaP (that is, attempt b)
        
        num_r = (grav_HR_output[k]/(rmid_cm_HR_output[k]*m*m*omega*omega)) * (   -dDeltaP_dr_b_r   +   (rho_HR_output[k]*f*rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*( -(2.0/rmid_cm_HR_output[k]) + (rmid_cm_HR_output[k]*rmid_cm_HR_output[k]/grav_HR_output[k])*dgrr_dr ) )  );
        
        num_i = (grav_HR_output[k]/(rmid_cm_HR_output[k]*m*m*omega*omega)) * (   -dDeltaP_dr_b_i   );
        
        denom_r = - (  dDeltaP_dr_b_r + 2.0*rho_HR_output[k]*f*rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*( (1.0/rmid_cm_HR_output[k]) + ((3.0*grav_HR_output[k])/(rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*m*m*omega*omega)) )  );
        
        denom_i = - dDeltaP_dr_b_i;
        
        V_div_xi_r_b_r = (  (num_r*denom_r)  +  (num_i*denom_i)  )/(  (denom_r*denom_r)  +  (denom_i*denom_i)  );
        
        V_div_xi_r_b_i = (  num_i*denom_r  -  num_r*denom_i  )/(  denom_r*denom_r  +  denom_i*denom_i  );
        
        
        
        
        // Here the analytic expression for xi_r is caluclated
        
        D_analytic = (1.0  -  l*(l+1.0)*grav_HR_output[k]*grav_HR_output[k]/(rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*m*m*m*m*omega*omega*omega*omega)   -   (rmid_cm_HR_output[k]*rmid_cm_HR_output[k]/(m*m*omega*omega))*dgrr_dr  );
        
        
        xi_r_analytic_r = (-1.0/(m*m*omega*omega*rho_HR_output[k]*D_analytic)) * (  -dDeltaP_dr_b_r - ( l*(l+1.0)*grav_HR_output[k]*delta_P_r/(rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*m*m*omega*omega) )  -  rho_HR_output[k]*( 2.0*f*rmid_cm_HR_output[k] + l*(l+1.0)*grav_HR_output[k]*f/(m*m*omega*omega) )  );
        
        xi_r_analytic_i = (-1.0/(m*m*omega*omega*rho_HR_output[k]*D_analytic)) * (  -dDeltaP_dr_b_i - ( l*(l+1.0)*grav_HR_output[k]*delta_P_i/(rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*m*m*omega*omega) )  );
        
        
        // Here the analytic expression for V is calculated
        
        V_analytic_r = (1.0/(D_analytic*m*m*omega*omega*radius_cm_HR_output[k]*rho_HR_output[k]))  *  (   (grav_HR_output[k]*dDeltaP_dr_b_r)/(m*m*omega*omega)   +   delta_P_r*( 1.0- (radius_cm_HR_output[k]*radius_cm_HR_output[k]*dgrr_dr)/(m*m*omega*omega) )   +   rho_HR_output[k]*( 2.0*f*radius_cm_HR_output[k]*grav_HR_output[k]/(m*m*omega*omega)  + f*radius_cm_HR_output[k]*radius_cm_HR_output[k]*( 1.0 - (radius_cm_HR_output[k]*radius_cm_HR_output[k]*dgrr_dr)/(m*m*omega*omega) ) )   );
        
        V_analytic_i = (1.0/(D_analytic*m*m*omega*omega*radius_cm_HR_output[k]*rho_HR_output[k]))  *  (   (grav_HR_output[k]*dDeltaP_dr_b_i)/(m*m*omega*omega)   +   delta_P_i*( 1.0- (radius_cm_HR_output[k]*radius_cm_HR_output[k]*dgrr_dr)/(m*m*omega*omega) )   );
        
        
        // Here the analytic expression for F' is calculated according to the earlier definitions such that F' = F0 * ( Fprime_forT * T'/T0  +  dT'/dT - dxi/dr  +  Fprime_for_p * p'/p0 )
        
        if ( k == J-1) {
            // dTprime_dT and dT0_dr are unchanged
            dTprime_dT_r = dTprime_dT_r;
            dTprime_dT_i = dTprime_dT_i;
            
            dT0_dr = dT0_dr;
            
            dxi_r_dr_r = dxi_r_dr_r;
            dxi_r_dr_i = dxi_r_dr_i;
            
            delta[0] = (rmid_cm_HR_output[k] - radius_cm_HR_output[k])/(rmid_cm_HR_output[k] - rmid_cm_HR_output[k-1]);
            delta[1] = (radius_cm_HR_output[k] - rmid_cm_HR_output[k-1])/(rmid_cm_HR_output[k] - rmid_cm_HR_output[k-1]);
            
            Fprime_analytic_r = ( Fprime_for_T[k]*(delta[0]*vr[k-1][1][0] + delta[1]*vr[k][1][0])  +  ( flux_HR_output[k] )*(dTprime_dT_r)  +  Fprime_for_p[k]*(delta[0]*vr[k-1][0][0] + delta[1]*vr[k][0][0]) );
            
            Fprime_analytic_i = ( Fprime_for_T[k]*(delta[0]*vi[k-1][1][0] + delta[1]*vi[k][1][0])  +  ( flux_HR_output[k] )*(dTprime_dT_r)  +  Fprime_for_p[k]*(delta[0]*vi[k-1][0][0] + delta[1]*vi[k][0][0]) );
            
        } else {
            
            dTprime_dT_r = ( temperature_HR_output[k+1]*vr[k+1][1][0]  -  temperature_HR_output[k]*vr[k][1][0] )/( temperature_HR_output[k+1]  -  temperature_HR_output[k] );
            dTprime_dT_i = ( temperature_HR_output[k+1]*vi[k+1][1][0]  -  temperature_HR_output[k]*vi[k][1][0] )/( temperature_HR_output[k+1]  -  temperature_HR_output[k] );
            
            dT0_dr = (temperature_HR_output[k+1] - temperature_HR_output[k])/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k]);
            
            dxi_r_dr_r = R*( ur[k+1][0][0] - ur[k][0][0] )/( radius_cm_HR_output[k+1] - radius_cm_HR_output[k] );
            
            dxi_r_dr_i = R*( ui[k+1][0][0] - ui[k][0][0] )/( radius_cm_HR_output[k+1] - radius_cm_HR_output[k] );
            
            delta[0] = (rmid_cm_HR_output[k+1] - radius_cm_HR_output[k])/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k]);
            delta[1] = (radius_cm_HR_output[k] - rmid_cm_HR_output[k])/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k]);
            
            Fprime_analytic_r = ( Fprime_for_T[k]*(delta[0]*vr[k][1][0] + delta[1]*vr[k+1][1][0])  +  ( flux_HR_output[k] )*(dTprime_dT_r)  +  Fprime_for_p[k]*(delta[0]*vr[k][0][0] + delta[1]*vr[k+1][0][0]) );
            
            Fprime_analytic_i = ( Fprime_for_T[k]*(delta[0]*vi[k][1][0] + delta[1]*vi[k+1][1][0])  +  ( flux_HR_output[k] )*(dTprime_dT_r)  +  Fprime_for_p[k]*(delta[0]*vi[k][0][0] + delta[1]*vi[k+1][0][0]) );
            
        }
        
        //cout << "Flag by Fprime, " << k << "\n";
        
        
        
        
        
        
        
        
        
        
        //This gets V_div stuff via the momentum equation
        
        otherVdiv_r = ( R*ur[k][0][0]*(pressure_HR_output[k]*vr[k][0][0]  +  rho_HR_output[k]*f*rmid_cm_HR_output[k]*rmid_cm_HR_output[k])   +   pressure_HR_output[k]*vi[k][0][0]*R*ui[k][0][0] )/( rho_HR_output[k]*m*m*omega*omega*R*R*( ur[k][0][0]*ur[k][0][0] + ui[k][0][0]*ui[k][0][0] ) );
        
        
        otherVdiv_i = ( R*ur[k][0][0]*(pressure_HR_output[k]*vi[k][0][0])   -   (pressure_HR_output[k]*vr[k][0][0]  +  rho_HR_output[k]*f*rmid_cm_HR_output[k]*rmid_cm_HR_output[k])*R*ur[k][0][0] )/( rho_HR_output[k]*m*m*omega*omega*R*R*( ur[k][0][0]*ur[k][0][0] + ui[k][0][0]*ui[k][0][0] ) );
        
        
        // the magnitude of the radial displacement
        mod_xi_radial = R * sqrt( (ur[k][0][0]*ur[k][0][0]) + (ui[k][0][0]*ui[k][0][0]) );
        
        // the magnitude of the lagrangian pressure perturbation
        mod_delta_P = sqrt( (delta_P_r*delta_P_r) + (delta_P_i*delta_P_i) );
        
        
        xi_h_real = (1.0/(rmid_cm_HR_output[k] * m * m * omega * omega)) * ( (pressure_HR_output[k] * vr[k][0][0] / rho_HR_output[k]) + (f * rmid_cm_HR_output[k] * rmid_cm_HR_output[k]) );
        
        xi_h_imaginary = (1.0/(rmid_cm_HR_output[k] * m * m * omega * omega)) * ( (pressure_HR_output[k] * vi[k][0][0] / rho_HR_output[k]) );
        
        mod_xi_h = sqrt(xi_h_real*xi_h_real + xi_h_imaginary*xi_h_imaginary);
        
        
        
        
        
        
        xi_h_over_xi_r_real = (1.0/R)*(xi_h_real*ur[k][0][0] + xi_h_imaginary*ui[k][0][0])/(ur[k][0][0]*ur[k][0][0] + ui[k][0][0]*ui[k][0][0]);
        
        xi_h_over_xi_r_im = (1.0/R)*(xi_h_imaginary*ur[k][0][0] - xi_h_real*ui[k][0][0])/(ur[k][0][0]*ur[k][0][0] + ui[k][0][0]*ui[k][0][0]);
        
        
        
        
        
        log_mod_xi_r = log10(R*sqrt((ur[k][0][0]*ur[k][0][0]) + (ui[k][0][0]*ui[k][0][0])));
        
        log_mod_xi_h = log10(sqrt((xi_h_real*xi_h_real) + (xi_h_imaginary*xi_h_imaginary)));
        
        
        
        
        // Here I get deltaP by using hydrostatic equilibrium to say dp/dr = - rho*grav
        if (k==0) {
            
            deltaP_hydro_r = pressure_HR_output[k]*vr[k][0][0] - R*0.5*(ur[k+1][0][0] + ur[k][0][0])*rho_HR_output[k]*grav_HR_output[k];
            
            deltaP_hydro_i = pressure_HR_output[k]*vi[k][0][0] - R*0.5*(ui[k+1][0][0] + ui[k][0][0])*rho_HR_output[k]*grav_HR_output[k];
            
            
            deltaP_hydro_next_r = pressure_HR_output[k+1]*vr[k+1][0][0] - R*0.5*(ur[k+2][0][0] + ur[k+1][0][0])*rho_HR_output[k+1]*grav_HR_output[k+1];
            
            deltaP_hydro_next_i = pressure_HR_output[k+1]*vi[k+1][0][0] - R*0.5*(ui[k+2][0][0] + ui[k+1][0][0])*rho_HR_output[k+1]*grav_HR_output[k+1];
            
            
            ddeltaP_dr_hydro_r = (deltaP_hydro_next_r - deltaP_hydro_r)/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k]);
            
            ddeltaP_dr_hydro_i = (deltaP_hydro_next_i - deltaP_hydro_i)/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k]);
            
            
        } else {
            
            deltaP_hydro_r = pressure_HR_output[k]*vr[k][0][0] - R*0.5*(ur[k][0][0] + ur[k-1][0][0])*rho_HR_output[k]*grav_HR_output[k];
            
            deltaP_hydro_i = pressure_HR_output[k]*vi[k][0][0] - R*0.5*(ui[k][0][0] + ui[k-1][0][0])*rho_HR_output[k]*grav_HR_output[k];
            
            if (k==J-1 || k==J-2) {
                
                // At the outer edge, we extrapolate its value from the inner values, assuming the gradient is constant
                // I use ddeltaP_dr_hydro from J-4 and J-3, and extrapolate from there
                
                ddeltaP_dr_hydro_r = ddeltaP_dr_hydro_Jminus3_r + (radius_cm_HR_output[k] - radius_cm_HR_output[J-3])*( (ddeltaP_dr_hydro_Jminus3_r - ddeltaP_dr_hydro_Jminus4_r)/(radius_cm_HR_output[J-3] - radius_cm_HR_output[J-4]) );
                
                ddeltaP_dr_hydro_i = ddeltaP_dr_hydro_i + (radius_cm_HR_output[k] - radius_cm_HR_output[J-3])*( (ddeltaP_dr_hydro_Jminus3_i - ddeltaP_dr_hydro_Jminus4_i)/(radius_cm_HR_output[J-3] - radius_cm_HR_output[J-4]) );

                
            } else {
                
                // Otherwise, we take the derivative numerically
                
                deltaP_hydro_next_r = pressure_HR_output[k+1]*vr[k+1][0][0] - R*0.5*(ur[k+1][0][0] + ur[k][0][0])*rho_HR_output[k+1]*grav_HR_output[k+1];
                
                deltaP_hydro_next_i = pressure_HR_output[k+1]*vi[k+1][0][0] - R*0.5*(ui[k+1][0][0] + ui[k][0][0])*rho_HR_output[k+1]*grav_HR_output[k+1];
                
                
                ddeltaP_dr_hydro_r = (deltaP_hydro_next_r - deltaP_hydro_r)/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k]);
                
                ddeltaP_dr_hydro_i = (deltaP_hydro_next_i - deltaP_hydro_i)/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k]);
                
                
            }
            
            
            if (k == J-4) {
                
                ddeltaP_dr_hydro_Jminus4_r = ddeltaP_dr_hydro_r;
                
                ddeltaP_dr_hydro_Jminus4_i = ddeltaP_dr_hydro_i;
                
            }
            
            if (k == J-3) {
                
                ddeltaP_dr_hydro_Jminus3_r = ddeltaP_dr_hydro_r;
                
                ddeltaP_dr_hydro_Jminus3_i = ddeltaP_dr_hydro_i;
                
            }
            
            
        }
        
        
        // This saves the old value so I can do a derivative if needed.
        V_analytic_hydro_old_r = V_analytic_hydro_r;
        
        if (k==0 || k==J-1) {
            
            D_analytic = (1.0  -  l*(l+1.0)*grav_HR_output[k]*grav_HR_output[k]/(rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*m*m*m*m*omega*omega*omega*omega)   -   (rmid_cm_HR_output[k]*rmid_cm_HR_output[k]/(m*m*omega*omega))*dgrr_dr  );
            
            // Here the analytic expression for xi_r is caluclated using the hydrostatic equilibrium assumption
            
            xi_r_analytic_hydro_r = (-1.0/(m*m*omega*omega*rho_HR_output[k]*D_analytic)) * (  -ddeltaP_dr_hydro_r - ( l*(l+1.0)*grav_HR_output[k]*deltaP_hydro_r/(rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*m*m*omega*omega) )  -  rho_HR_output[k]*( 2.0*f*rmid_cm_HR_output[k] + l*(l+1.0)*grav_HR_output[k]*f/(m*m*omega*omega) )  );
            
            xi_r_analytic_hydro_i = (-1.0/(m*m*omega*omega*rho_HR_output[k]*D_analytic)) * (  -ddeltaP_dr_hydro_i - ( l*(l+1.0)*grav_HR_output[k]*deltaP_hydro_i/(rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*m*m*omega*omega) )  );
            
            
            // Here the analytic expression for V is calculated
            
            V_analytic_hydro_r = (1.0/(D_analytic*m*m*omega*omega*radius_cm_HR_output[k]*rho_HR_output[k]))  *  (   (grav_HR_output[k]*ddeltaP_dr_hydro_r)/(m*m*omega*omega)   +   deltaP_hydro_r*( 1.0 - (radius_cm_HR_output[k]*radius_cm_HR_output[k]*dgrr_dr)/(m*m*omega*omega) )   +   rho_HR_output[k]*( 2.0*f*radius_cm_HR_output[k]*grav_HR_output[k]/(m*m*omega*omega)  + f*radius_cm_HR_output[k]*radius_cm_HR_output[k]*( 1.0 - (radius_cm_HR_output[k]*radius_cm_HR_output[k]*dgrr_dr)/(m*m*omega*omega) ) )   );
            
            V_analytic_hydro_i = (1.0/(D_analytic*m*m*omega*omega*radius_cm_HR_output[k]*rho_HR_output[k]))  *  (   (grav_HR_output[k]*ddeltaP_dr_hydro_i)/(m*m*omega*omega)   +   deltaP_hydro_i*( 1.0 - (radius_cm_HR_output[k]*radius_cm_HR_output[k]*dgrr_dr)/(m*m*omega*omega) )   );
            
        } else {
            
            
            // Here I set up delta to be used for this case (because xi_r and xi_h are both evaluated at the cell face)
            
            delta[0] = ( rmid_cm_HR_output[k+1] - radius_cm_HR_output[k] )/( rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k] );
            
            delta[1] = ( radius_cm_HR_output[k] - rmid_cm_HR_output[k] )/( rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k] );
            
            
            D_analytic = (1.0  -  l*(l+1.0)*grav_HR_output[k]*grav_HR_output[k]/(rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*m*m*m*m*omega*omega*omega*omega)   -   (rmid_cm_HR_output[k]*rmid_cm_HR_output[k]/(m*m*omega*omega))*dgrr_dr  );
            
            D_analytic_next = (1.0  -  l*(l+1.0)*grav_HR_output[k+1]*grav_HR_output[k+1]/(rmid_cm_HR_output[k+1]*rmid_cm_HR_output[k+1]*m*m*m*m*omega*omega*omega*omega)   -   (rmid_cm_HR_output[k+1]*rmid_cm_HR_output[k+1]/(m*m*omega*omega))*dgrr_dr  );
            
            // Here the analytic expression for xi_r is caluclated using the hydrostatic equilibrium assumption
            
            xi_r_analytic_hydro_r = (-1.0/(m*m*omega*omega*rho_HR_output[k]*D_analytic)) * (  -ddeltaP_dr_hydro_r - ( l*(l+1.0)*grav_HR_output[k]*deltaP_hydro_r/(rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*m*m*omega*omega) )  -  rho_HR_output[k]*( 2.0*f*rmid_cm_HR_output[k] + l*(l+1.0)*grav_HR_output[k]*f/(m*m*omega*omega) )  );
            
            xi_r_analytic_hydro_i = (-1.0/(m*m*omega*omega*rho_HR_output[k]*D_analytic)) * (  -ddeltaP_dr_hydro_i - ( l*(l+1.0)*grav_HR_output[k]*deltaP_hydro_i/(rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*m*m*omega*omega) )  );
            
            xi_r_analytic_hydro_next_r = (-1.0/(m*m*omega*omega*rho_HR_output[k+1]*D_analytic_next)) * (  -ddeltaP_dr_hydro_r - ( l*(l+1.0)*grav_HR_output[k+1]*deltaP_hydro_next_r/(rmid_cm_HR_output[k+1]*rmid_cm_HR_output[k+1]*m*m*omega*omega) )  -  rho_HR_output[k+1]*( 2.0*f*rmid_cm_HR_output[k+1] + l*(l+1.0)*grav_HR_output[k+1]*f/(m*m*omega*omega) )  );
            
            xi_r_analytic_hydro_next_i = (-1.0/(m*m*omega*omega*rho_HR_output[k+1]*D_analytic_next)) * (  -ddeltaP_dr_hydro_i - ( l*(l+1.0)*grav_HR_output[k+1]*deltaP_hydro_next_i/(rmid_cm_HR_output[k+1]*rmid_cm_HR_output[k+1]*m*m*omega*omega) )  );
            
            // Here the xi_r values are combined as a weighted average
            
            xi_r_analytic_hydro_r = delta[0]*xi_r_analytic_hydro_r + delta[1]*xi_r_analytic_hydro_next_r;
            
            xi_r_analytic_hydro_i = delta[0]*xi_r_analytic_hydro_i + delta[1]*xi_r_analytic_hydro_next_i;
            
            
            // Here the analytic expression for V is calculated
            
            V_analytic_hydro_r = (1.0/(D_analytic*m*m*omega*omega*radius_cm_HR_output[k]*rho_HR_output[k]))  *  (   (grav_HR_output[k]*ddeltaP_dr_hydro_r)/(m*m*omega*omega)   +   deltaP_hydro_r*( 1.0 - (radius_cm_HR_output[k]*radius_cm_HR_output[k]*dgrr_dr)/(m*m*omega*omega) )   +   rho_HR_output[k]*( 2.0*f*radius_cm_HR_output[k]*grav_HR_output[k]/(m*m*omega*omega)  + f*radius_cm_HR_output[k]*radius_cm_HR_output[k]*( 1.0 - (radius_cm_HR_output[k]*radius_cm_HR_output[k]*dgrr_dr)/(m*m*omega*omega) ) )   );
            
            V_analytic_hydro_i = (1.0/(D_analytic*m*m*omega*omega*radius_cm_HR_output[k]*rho_HR_output[k]))  *  (   (grav_HR_output[k]*ddeltaP_dr_hydro_i)/(m*m*omega*omega)   +   deltaP_hydro_i*( 1.0 - (radius_cm_HR_output[k]*radius_cm_HR_output[k]*dgrr_dr)/(m*m*omega*omega) )   );
            
            V_analytic_hydro_next_r = (1.0/(D_analytic_next*m*m*omega*omega*radius_cm_HR_output[k+1]*rho_HR_output[k+1]))  *  (   (grav_HR_output[k+1]*ddeltaP_dr_hydro_r)/(m*m*omega*omega)   +   deltaP_hydro_next_r*( 1.0 - (radius_cm_HR_output[k+1]*radius_cm_HR_output[k+1]*dgrr_dr)/(m*m*omega*omega) )   +   rho_HR_output[k+1]*( 2.0*f*radius_cm_HR_output[k+1]*grav_HR_output[k+1]/(m*m*omega*omega)  + f*radius_cm_HR_output[k+1]*radius_cm_HR_output[k+1]*( 1.0 - (radius_cm_HR_output[k+1]*radius_cm_HR_output[k+1]*dgrr_dr)/(m*m*omega*omega) ) )   );
            
            V_analytic_hydro_next_i = (1.0/(D_analytic_next*m*m*omega*omega*radius_cm_HR_output[k+1]*rho_HR_output[k+1]))  *  (   (grav_HR_output[k+1]*ddeltaP_dr_hydro_i)/(m*m*omega*omega)   +   deltaP_hydro_next_i*( 1.0 - (radius_cm_HR_output[k+1]*radius_cm_HR_output[k+1]*dgrr_dr)/(m*m*omega*omega) )   );
            
            // Here the V values are combined as a weighted average
            
            V_analytic_hydro_r = delta[0]*V_analytic_hydro_r + delta[1]*V_analytic_hydro_next_r;
            
            V_analytic_hydro_i = delta[0]*V_analytic_hydro_i + delta[1]*V_analytic_hydro_next_i;
            
        }
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        // This is for the comparison between p' and the first-order equation
        
        pprime_comp_r = ( m*m*omega*omega*xi_h_real - f*radius_cm_HR_output[k]*radius_cm_HR_output[k] ) *rho_HR_output[k];
        
        pprime_comp_i = ( m*m*omega*omega*xi_h_imaginary ) *rho_HR_output[k];
        
        
        
        // Here I produce the functions G and H, which are used in second order corrections for P'
        
        if (k==J-1) {
            
            // Because we can't get it without a cell outside, we have issues at the surface
            Gvar_r = Gvar_r;
            Gvar_i = Gvar_i;
            Hvar_r = Hvar_r;
            Hvar_i = Hvar_i;
            
        } else {
            
            gradient = ( (rho_face_HR_output[k+1] - rho_face_HR_output[k])/(radius_cm_HR_output[k+1] - radius_cm_HR_output[k]) ); // (  ( (radius_cm_HR_output[k+1]*radius_cm_HR_output[k+1]*rho_face_HR_output[k+1]*ur[k+1][0][0])  -  (radius_cm_HR_output[k]*radius_cm_HR_output[k]*rho_face_HR_output[k]*ur[k][0][0]) ) / (  radius_cm_HR_output[k+1] - radius_cm_HR_output[k]  )  );
            
            //cout << radius_cm_HR_output[k]*radius_cm_HR_output[k]*rho_face_HR_output[k]*ur[k][0][0] << "\n";
            
            Gvar_r = ((m*m*omega*omega)/(l*(l+1.0))) * (R/rmid_cm_HR_output[k]) * gradient;
            
            Gvar_i = ((m*m*omega*omega)/(l*(l+1.0))) * (R/rmid_cm_HR_output[k]) * (  ( (radius_cm_HR_output[k+1]*radius_cm_HR_output[k+1]*rho_face_HR_output[k+1]*ui[k+1][0][0])  -  (radius_cm_HR_output[k]*radius_cm_HR_output[k]*rho_face_HR_output[k]*ui[k][0][0]) ) / (  radius_cm_HR_output[k+1] - radius_cm_HR_output[k]  )  );
            
            Hvar_r = (  (m*m*omega*omega*radius_cm_HR_output[k]*radius_cm_HR_output[k]*radius_cm_HR_output[k]*radius_cm_HR_output[k]*f) / (grav_HR_output[k]*l*(l+1.0)) )  * ( (rho_face_HR_output[k+1] - rho_face_HR_output[k])/(radius_cm_HR_output[k+1] - radius_cm_HR_output[k]) ) ;
            
            Hvar_i = 0.0;
            
            
        }
        
        // Now the second order expressions for p' are evaluated
        
        pprime_comp_second_r = (  - f*radius_cm_HR_output[k]*radius_cm_HR_output[k] ) *rho_HR_output[k]  +  Gvar_r  +  Hvar_r;
        
        pprime_comp_second_i =   Gvar_i  +  Hvar_i;
        
        
        
        
        // Equilibrium tide
        xi_r_eq = (-f*radius_cm_HR_output[k]*radius_cm_HR_output[k]/grav_HR_output[k]);
        
        // This saves the previous value of V_cont_r, so I can do gradients of it later
        
        V_cont_r_old = V_cont_r;
        
        
        if (k == J-1 || k == 0) {
            
            if (k==0) {
                
                V_cont_r = (radius_cm_HR_output[k]/(l*(l+1.0))) * (  0.5*( (1.0/chiRho_HR_output[k])*( vr[k][0][0] - chiT_HR_output[k]*vr[k][1][0] ) + (1.0/chiRho_HR_output[k+1])*( vr[k+1][0][0] - chiT_HR_output[k+1]*vr[k][1][0] ) )    +    R*ur[k][0][0]*(log(rho_HR_output[k+1]) - log(rho_HR_output[k]))/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k])   +   (R/(radius_cm_HR_output[k]*radius_cm_HR_output[k]))*( (radius_cm_HR_output[k+1]*radius_cm_HR_output[k+1]*ur[k+1][0][0] - radius_cm_HR_output[k]*radius_cm_HR_output[k]*ur[k][0][0])/(radius_cm_HR_output[k+1] - radius_cm_HR_output[k]) )  );
                
                V_cont_i = (radius_cm_HR_output[k]/(l*(l+1.0))) * (  0.5*( (1.0/chiRho_HR_output[k])*( vi[k][0][0] - chiT_HR_output[k]*vi[k][1][0] ) + (1.0/chiRho_HR_output[k+1])*( vi[k+1][0][0] - chiT_HR_output[k+1]*vi[k][1][0] ) )    +    R*ui[k][0][0]*(log(rho_HR_output[k+1]) - log(rho_HR_output[k]))/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k])   +   (R/(radius_cm_HR_output[k]*radius_cm_HR_output[k]))*( (radius_cm_HR_output[k+1]*radius_cm_HR_output[k+1]*ui[k+1][0][0] - radius_cm_HR_output[k]*radius_cm_HR_output[k]*ui[k][0][0])/(radius_cm_HR_output[k+1] - radius_cm_HR_output[k]) )  );
                
                
            } else {
                
                V_cont_r = (radius_cm_HR_output[k]/(l*(l+1.0))) * (  (1.0/chiRho_HR_output[k])*( vr[k][0][0] - chiT_HR_output[k]*vr[k][1][0] )    +    R*ur[k][0][0]*(log(rho_HR_output[k]) - log(rho_HR_output[k-1]))/(rmid_cm_HR_output[k] - rmid_cm_HR_output[k-1])   +   (R/(radius_cm_HR_output[k]*radius_cm_HR_output[k]))*( (radius_cm_HR_output[k]*radius_cm_HR_output[k]*ur[k][0][0] - radius_cm_HR_output[k-1]*radius_cm_HR_output[k-1]*ur[k-1][0][0])/(radius_cm_HR_output[k] - radius_cm_HR_output[k-1]) )  );
                
                V_cont_i = (radius_cm_HR_output[k]/(l*(l+1.0))) * (  (1.0/chiRho_HR_output[k])*( vi[k][0][0] - chiT_HR_output[k]*vi[k][1][0] )    +    R*ui[k][0][0]*(log(rho_HR_output[k]) - log(rho_HR_output[k-1]))/(rmid_cm_HR_output[k] - rmid_cm_HR_output[k-1])   +   (R/(radius_cm_HR_output[k]*radius_cm_HR_output[k]))*( (radius_cm_HR_output[k]*radius_cm_HR_output[k]*ui[k][0][0] - radius_cm_HR_output[k-1]*radius_cm_HR_output[k-1]*ui[k-1][0][0])/(radius_cm_HR_output[k] - radius_cm_HR_output[k-1]) )  );
                
            }
            
            
        } else {
            
            delta[0] = ( rmid_cm_HR_output[k+1] - radius_cm_HR_output[k] )/( rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k] );
            
            delta[1] = ( radius_cm_HR_output[k] - rmid_cm_HR_output[k] )/( rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k] );
            
            V_cont_r = (radius_cm_HR_output[k]/(l*(l+1.0))) * (  ( (delta[0]/chiRho_HR_output[k])*( vr[k][0][0] - chiT_HR_output[k]*vr[k][1][0] ) + (delta[1]/chiRho_HR_output[k+1])*( vr[k+1][0][0] - chiT_HR_output[k+1]*vr[k][1][0] ) )    +    R*ur[k][0][0]*(log(rho_HR_output[k+1]) - log(rho_HR_output[k]))/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k])   +   (R/(radius_cm_HR_output[k]*radius_cm_HR_output[k]))*( (radius_cm_HR_output[k+1]*radius_cm_HR_output[k+1]*ur[k+1][0][0] - radius_cm_HR_output[k-1]*radius_cm_HR_output[k-1]*ur[k-1][0][0])/(radius_cm_HR_output[k+1] - radius_cm_HR_output[k-1]) )  );
            
            V_cont_i = (radius_cm_HR_output[k]/(l*(l+1.0))) * (  ( (delta[0]/chiRho_HR_output[k])*( vi[k][0][0] - chiT_HR_output[k]*vi[k][1][0] ) + (delta[1]/chiRho_HR_output[k+1])*( vi[k+1][0][0] - chiT_HR_output[k+1]*vi[k][1][0] ) )    +    R*ui[k][0][0]*(log(rho_HR_output[k+1]) - log(rho_HR_output[k]))/(rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k])   +   (R/(radius_cm_HR_output[k]*radius_cm_HR_output[k]))*( (radius_cm_HR_output[k+1]*radius_cm_HR_output[k+1]*ui[k+1][0][0] - radius_cm_HR_output[k-1]*radius_cm_HR_output[k-1]*ui[k-1][0][0])/(radius_cm_HR_output[k+1] - radius_cm_HR_output[k-1]) )  );
            
        }
        
        

        
        // Looking at the entropy perturbation
        
        s_prime_r = cp_HR_output[k] * ( vr[k][1][0] - grada_HR_output[k]*vr[k][0][0] );
        s_prime_i = cp_HR_output[k] * ( vi[k][1][0] - grada_HR_output[k]*vi[k][0][0] );
        
        delta_s_r = s_prime_r + ( cp_HR_output[k] * (brunt_A_HR_output[k]/rmid_cm_HR_output[k]) * (chiRho_HR_output[k]/chiT_HR_output[k]) * R * ur[k][0][0] );
        delta_s_i = s_prime_i + ( cp_HR_output[k] * (brunt_A_HR_output[k]/rmid_cm_HR_output[k]) * (chiRho_HR_output[k]/chiT_HR_output[k]) * R * ui[k][0][0] );
        
        if (brunt_A_HR_output[k] <= 0.0) {
            
            if (k < J-1) {
                
                F_conv_prime_div_F_conv_r = ( 1.0/( cp_HR_output[k] * (brunt_A_HR_output[k]/rmid_cm_HR_output[k]) * (chiRho_HR_output[k]/chiT_HR_output[k]) ) ) * ( cp_HR_output[k+1] * ( vr[k+1][1][0] - grada_HR_output[k+1]*vr[k+1][0][0] )  -  cp_HR_output[k] * ( vr[k][1][0] - grada_HR_output[k]*vr[k][0][0] ) ) / ( rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k] );
                
                F_conv_prime_div_F_conv_i = ( 1.0/( cp_HR_output[k] * (brunt_A_HR_output[k]/rmid_cm_HR_output[k]) * (chiRho_HR_output[k]/chiT_HR_output[k]) ) ) * ( cp_HR_output[k+1] * ( vi[k+1][1][0] - grada_HR_output[k+1]*vi[k+1][0][0] )  -  cp_HR_output[k] * ( vi[k][1][0] - grada_HR_output[k]*vi[k][0][0] ) ) / ( rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k] );
                
                prefactor = ( flux_tot_HR_output[k] * conv_L_div_L_HR_output[k] * ( (1.0/(K_HR_output[k] * cp_HR_output[k] * (brunt_A_HR_output[k]/rmid_cm_HR_output[k]) * (chiRho_HR_output[k]/chiT_HR_output[k]) )) ) / ( temperature_HR_output[k+1] - temperature_HR_output[k] ) );
                
                
            } else {
                
                F_conv_prime_div_F_conv_r = ( 1.0/( cp_HR_output[k] * (brunt_A_HR_output[k]/rmid_cm_HR_output[k]) * (chiRho_HR_output[k]/chiT_HR_output[k]) ) ) * ( cp_HR_output[k] * ( vr[k][1][0] - grada_HR_output[k]*vr[k][0][0] )  -  cp_HR_output[k-1] * ( vr[k-1][1][0] - grada_HR_output[k-1]*vr[k-1][0][0] ) ) / ( rmid_cm_HR_output[k] - rmid_cm_HR_output[k-1] );
                
                F_conv_prime_div_F_conv_i = ( 1.0/( cp_HR_output[k] * (brunt_A_HR_output[k]/rmid_cm_HR_output[k]) * (chiRho_HR_output[k]/chiT_HR_output[k]) ) ) * ( cp_HR_output[k] * ( vi[k][1][0] - grada_HR_output[k]*vi[k][0][0] )  -  cp_HR_output[k-1] * ( vi[k-1][1][0] - grada_HR_output[k-1]*vi[k-1][0][0] ) ) / ( rmid_cm_HR_output[k] - rmid_cm_HR_output[k-1] );
                
                // We just leave this one
                prefactor = prefactor;
                
            }
            
            
        } else {
            
            F_conv_prime_div_F_conv_r = 0.0;
            
            F_conv_prime_div_F_conv_i = 0.0;
            
        }
        
        
        
        
        if (k<J-1) {
            
            V_rho_r = ( radius_cm_HR_output[k+1]*radius_cm_HR_output[k+1]*rho_face_HR_output[k+1]*R*ur[k+1][0][0]  -  radius_cm_HR_output[k]*radius_cm_HR_output[k]*rho_face_HR_output[k]*R*ur[k][0][0] )/( l*(l+1.0)*rmid_cm_HR_output[k+1]*rho_HR_output[k+1]*( radius_cm_HR_output[k+1] - radius_cm_HR_output[k] ) );
            
            V_rho_i = ( radius_cm_HR_output[k+1]*radius_cm_HR_output[k+1]*rho_face_HR_output[k+1]*R*ui[k+1][0][0]  -  radius_cm_HR_output[k]*radius_cm_HR_output[k]*rho_face_HR_output[k]*R*ui[k][0][0] )/( l*(l+1.0)*rmid_cm_HR_output[k+1]*rho_HR_output[k+1]*( radius_cm_HR_output[k+1] - radius_cm_HR_output[k] ) );
            
            ds_prime_ds0_r = ( ( cp_HR_output[k+1] * ( vr[k+1][1][0] - grada_HR_output[k+1]*vr[k+1][0][0] ) - cp_HR_output[k] * ( vr[k][1][0] - grada_HR_output[k]*vr[k][0][0] ) )  /  (  rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k]  ) )   /   (   cp_HR_output[k] * (brunt_A_HR_output[k] / rmid_cm_HR_output[k]) * (chiRho_HR_output[k] / chiT_HR_output[k])   );
            
            ds_prime_ds0_i = ( ( cp_HR_output[k+1] * ( vi[k+1][1][0] - grada_HR_output[k+1]*vi[k+1][0][0] ) - cp_HR_output[k] * ( vi[k][1][0] - grada_HR_output[k]*vi[k][0][0] ) )  /  (  rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k]  ) )   /   (   cp_HR_output[k] * (brunt_A_HR_output[k] / rmid_cm_HR_output[k]) * (chiRho_HR_output[k] / chiT_HR_output[k])   );
            
            dr_dT_Fc_K0 = ( (rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k])/(temperature_HR_output[k+1] - temperature_HR_output[k]) ) * ( flux_tot_HR_output[k] * conv_L_div_L_HR_output[k] / K_HR_output[k] );
            
        } else {
            
            V_rho_r = V_rho_r;
            
            V_rho_i = V_rho_i;
            
            ds_prime_ds0_r = ( ( cp_HR_output[k] * ( vr[k][1][0] - grada_HR_output[k]*vr[k][0][0] ) - cp_HR_output[k-1] * ( vr[k-1][1][0] - grada_HR_output[k-1]*vr[k-1][0][0] ) )  /  (  rmid_cm_HR_output[k] - rmid_cm_HR_output[k-1]  ) )   /   (   cp_HR_output[k] * (brunt_A_HR_output[k] / rmid_cm_HR_output[k]) * (chiRho_HR_output[k] / chiT_HR_output[k])   );
            
            ds_prime_ds0_i = ( ( cp_HR_output[k] * ( vi[k][1][0] - grada_HR_output[k]*vi[k][0][0] ) - cp_HR_output[k-1] * ( vi[k-1][1][0] - grada_HR_output[k-1]*vi[k-1][0][0] ) )  /  (  rmid_cm_HR_output[k] - rmid_cm_HR_output[k-1]  ) )   /   (   cp_HR_output[k] * (brunt_A_HR_output[k] / rmid_cm_HR_output[k]) * (chiRho_HR_output[k] / chiT_HR_output[k])   );
            
            dr_dT_Fc_K0 = ( (rmid_cm_HR_output[k] - rmid_cm_HR_output[k-1])/(temperature_HR_output[k] - temperature_HR_output[k-1]) ) * ( flux_tot_HR_output[k] * conv_L_div_L_HR_output[k] / K_HR_output[k] );
            
        }
        
        
        
        x_conv = 24*5.67e-5*temperature_HR_output[k]*temperature_HR_output[k]*temperature_HR_output[k] / ( rho_HR_output[k]*rho_HR_output[k] * cp_HR_output[k] * 1e8 * 1e5 * opacity_HR_output[k] );
        
        if (radius_cm_HR_output[k]/radius_cm_HR_output[J-1] > 0.995) {
            
            if (k_start_JP == 0) {
                
                k_start_JP = k;
                
            }
            
            s_prime_JP = s_prime_JP + ur[k_start_JP][1][0]*cp_HR_output[k]*(brunt_A_HR_output[k]/rmid_cm_HR_output[k])*(chiRho_HR_output[k]/chiT_HR_output[k])*(radius_cm_HR_output[k]-radius_cm_HR_output[k-1]);
            
        }
        
        
        if (k==0) {
            
            xi_r_JP = 0.0;
            
            xi_r_JP_term = 0.0;
            
            s_prime_equation_JP = 0.0;
            
        } else {
            
            if (k > 50 && k+50 < J-1) {
                
                xi_r_JP =   (1.0/gamma1_HR_output[k]) * (    xi_h_r*( 1.0 + gamma1_HR_output[k]*brunt_A_HR_output[k] )  +  radius_cm_HR_output[k] * (1.0/(rmid_cm_HR_output[k] * m * m * omega * omega)) * ( ( (( (pressure_HR_output[k+50] * vr[k+50][0][0] / rho_HR_output[k+50]) + (f * rmid_cm_HR_output[k+50] * rmid_cm_HR_output[k+50]) )) - ( ( (pressure_HR_output[k-50] * vr[k-50][0][0] / rho_HR_output[k-50]) + (f * rmid_cm_HR_output[k-50] * rmid_cm_HR_output[k-50]) ) ) )/(radius_cm_HR_output[k+50]-radius_cm_HR_output[k-50]) ) + s_prime_JP*( grav_HR_output[k] / (m*m*omega*omega) )*( gamma1_HR_output[k]/cp_HR_output[k] )*( chiT_HR_output[k]/chiRho_HR_output[k] )   );
                
                xi_r_JP_term =   (1.0/gamma1_HR_output[k]) * (    radius_cm_HR_output[k] * (1.0/(rmid_cm_HR_output[k] * m * m * omega * omega)) * ( ( (( (pressure_HR_output[k+50] * vr[k+50][0][0] / rho_HR_output[k+50]) + (f * rmid_cm_HR_output[k+50] * rmid_cm_HR_output[k+50]) )) - ( ( (pressure_HR_output[k-50] * vr[k-50][0][0] / rho_HR_output[k-50]) + (f * rmid_cm_HR_output[k-50] * rmid_cm_HR_output[k-50]) ) ) )/(radius_cm_HR_output[k+50]-radius_cm_HR_output[k-50]) )    );
                
                s_prime_equation_JP = ( cp_HR_output[k] / grav_HR_output[k] ) * ( chiRho_HR_output[k] / chiT_HR_output[k] ) * (   f*rmid_cm_HR_output[k]*brunt_A_HR_output[k]   +   m*m*omega*omega*(  R*ur[k][0][0]  -  brunt_A_HR_output[k]*xi_h_r  )  -  (1.0/gamma1_HR_output[k])*(  ( (pressure_HR_output[k+50]*vr[k+50][0][0]/rho_HR_output[k+50]  +  f*radius_cm_HR_output[k+50]*radius_cm_HR_output[k+50])  -  (pressure_HR_output[k-50]*vr[k-50][0][0]/rho_HR_output[k-50]  +  f*radius_cm_HR_output[k-50]*radius_cm_HR_output[k-50])  )/(  radius_cm_HR_output[k+50]  -  radius_cm_HR_output[k-50]  )  )     );
                
            } else {
                
                xi_r_JP =   (1.0/gamma1_HR_output[k]) * (    xi_h_r*( 1.0 + gamma1_HR_output[k]*brunt_A_HR_output[k] )  +  radius_cm_HR_output[k] * (1.0/(rmid_cm_HR_output[k] * m * m * omega * omega)) * ( ( (( (pressure_HR_output[k] * vr[k][0][0] / rho_HR_output[k]) + (f * rmid_cm_HR_output[k] * rmid_cm_HR_output[k]) )) - ( ( (pressure_HR_output[k-1] * vr[k-1][0][0] / rho_HR_output[k-1]) + (f * rmid_cm_HR_output[k-1] * rmid_cm_HR_output[k-1]) ) ) )/(radius_cm_HR_output[k]-radius_cm_HR_output[k-1]) ) + s_prime_JP*( grav_HR_output[k] / (m*m*omega*omega) )*( gamma1_HR_output[k]/cp_HR_output[k] )*( chiT_HR_output[k]/chiRho_HR_output[k] )   );
                
                xi_r_JP_term =   (1.0/gamma1_HR_output[k]) * (    radius_cm_HR_output[k] * (1.0/(rmid_cm_HR_output[k] * m * m * omega * omega)) * ( ( (( (pressure_HR_output[k] * vr[k][0][0] / rho_HR_output[k]) + (f * rmid_cm_HR_output[k] * rmid_cm_HR_output[k]) )) - ( ( (pressure_HR_output[k-1] * vr[k-1][0][0] / rho_HR_output[k-1]) + (f * rmid_cm_HR_output[k-1] * rmid_cm_HR_output[k-1]) ) ) )/(radius_cm_HR_output[k]-radius_cm_HR_output[k-1]) )    );
                
                s_prime_equation_JP = ( cp_HR_output[k] / grav_HR_output[k] ) * ( chiRho_HR_output[k] / chiT_HR_output[k] ) * (   f*rmid_cm_HR_output[k]*brunt_A_HR_output[k]   +   m*m*omega*omega*(  R*ur[k][0][0]  -  brunt_A_HR_output[k]*xi_h_r  )  -  (1.0/gamma1_HR_output[k])*(  ( (pressure_HR_output[k]*vr[k][0][0]/rho_HR_output[k]  +  f*radius_cm_HR_output[k+50]*radius_cm_HR_output[k+50])  -  (pressure_HR_output[k-50]*vr[k-50][0][0]/rho_HR_output[k-50]  +  f*radius_cm_HR_output[k-1]*radius_cm_HR_output[k-1])  )/(  radius_cm_HR_output[k]  -  radius_cm_HR_output[k-1]  )  )     );
                
            }
            
            
            // xi_h_real = (1.0/(rmid_cm_HR_output[k] * m * m * omega * omega)) * ( (pressure_HR_output[k] * vr[k][0][0] / rho_HR_output[k]) + (f * rmid_cm_HR_output[k] * rmid_cm_HR_output[k]) );
            
        }
        
        
        
        /*
         cout << radius_cm[k] << "\tReal parts\n";
         cout << ur[k][0][0] << "\t\t" << vr[k][0][0] << "\t\t\t\t" << alphar[k][0][0] << "\t\t" << alphar[k][0][1] << "\t\t\t\t" << gammar[k][0][0] << "\n" ;
         cout << ur[k][1][0] << "\t\t" << vr[k][1][0] << "\t\t\t\t" << alphar[k][1][0] << "\t\t" << alphar[k][1][1] << "\t\t\t\t" << gammar[k][1][0] << "\n\n" ;
         
         cout << radius_cm[k] << "\tImaginary parts\n";
         cout << ui[k][0][0] << "\t\t" << vi[k][0][0] << "\t\t\t\t" << alphai[k][0][0] << "\t\t" << alphai[k][0][1] << "\t\t\t\t" << gammai[k][0][0] << "\n" ;
         cout << ui[k][1][0] << "\t\t" << vi[k][1][0] << "\t\t\t\t" << alphai[k][1][0] << "\t\t" << alphai[k][1][1] << "\t\t\t\t" << gammai[k][1][0] << "\n\n" ;
         */
        
        
        /* This writes the data to the file, in the order:
         1 - radius_cm
         2 - a
         3 - b
         4 - c
         5 - d
         6 - alphar[0][0]
         7 - alphar[0][1]
         8 - alphar[1][0]
         9 - alphar[1][1]
         10- alphai[0][0]
         11- alphai[0][1]
         12- alphai[1][0]
         13- alphai[1][1]
         14- gammar[0][0]
         15- gammar[1][0]
         16- gammai[0][0]
         17- gammai[1][0]
         18- /xi_{r} m omega in units of mp/(mp + Mstar)
         19- test
         20- radius_cm_HR
         21- flux_HR
         22- pressure_HR
         23- temperature_HR
         24- rmid_cm/R
         25- xi  (that is, a * R)
         26- F' (b * flux_BC)
         27- p' (c * pressure)
         28- T' (d * temperature)
         29- m * omega * xi_r in units given in Terquem, 1998 (figure 1): ( mp / (mp + Mstar) ) m/s
         30- a (imaginary part)
         31- b (imaginary part)
         32- c (imaginary part)
         33- d (imaginary part)
         34- mod(a) - currently all positive, as it's just a straight up mod, which makes the graph look a bit funky
         35- phase(a)
         36- xi_theta_r * tan(theta) (Terquem units)
         37- xi_theta_i * tan(theta) (Terquem units)
         38- xi_phi_r (Terquem units)
         39- xi_phi_i (Terquem units)
         40- xi_h_r (defined in Terquem, 98 and in Terquem units)
         41- xi_h_r (defined in Terquem, 98 and in Terquem units)
         42- rho_HR_output
         43- - f * r * r * rho    -->   (p'_{eq})
         44- f * r * r * rho / (dp/dr)   -->   (xi_{r, eq})
         45- p'_im
         46- xi_h_r (in Terquem units)
         47- xi_h_pprime_part (in Terquem units)
         48- eta_Terquem_output - measure of non-adiabaticity
         49- (r/h)*xi_r,real
         50- (r/h)*xi_r,im
         51- (r/h)*mod(xi_r)
         52- xi_h_real
         53- xi_h_imaginary
         54- mod(xi_h)
         55- H_p - pressure scale height
         56- H_rho - density scale height
         57- delta_P_r (real part of Lagrnagian pressure perturbation)
         58- delta_P_i (imaginaruy part of lagrangian pressure perturbation)
         59- mod(delta_P) (magnitude of lagrangian pressure perturbation)
         60- V_div_xi_r_r
         61- V_div_xi_r_i
         62- otherVdiv_r
         63- otherVdiv_i
         64- pprime_comp_r
         65- pprime_comp_i
         66- Gvar_r
         67- Gvar_i
         68- gradient
         69- rho_face_HR_output
         70- Hvar_r
         71- pprime_comp_second_r (second order expression for p' to be compared against the modelled value)
         72- pprime_comp_second_i (second order expression for p')
         73- num_r
         74- num_i
         75- denom_r
         76- denom_i
         77- ddelta_P_dr_r
         78- dgrr_dr * r*r/g
         79- grav/(r*m*m*omega*omega)
         80- from eq 53, real part of (xi_h / (V/xi_r))
         81- from eq 53, imaginary part of (xi_h / (V/xi_r))
         82- dp0dr (pressure gradient)
         83- dp0dr_old (simpler version of pressure gradient, for comparison)
         84- dp0dr_avg
         85- delta_P_r_new  -- BEWARE as this breaks when the cells change size, and it breaks in quite a big way
         86- delta_P_i_new
         87- xi_r_eq = - f * r * r / g
         88- (xi_h/xi_r)_real
         89- (xi_h/xi_r)_imaginary
         90- log10(mod(xi_r))
         91- log10(mod(xi_h))
         92- mod(xi_h/xi_r)
         93- mod(V/xi_r)
         94- dp_dr (second order interpolation for gradient of background pressure)
         95- d2p_dr2 (second order interpolation for secondderivative of background pressure)
         96- d2p0_dr2 (essentially first order interpolation kind of thing)
         97- dDeltaP_dr_b_r
         98- dDeltaP_dr_b_i
         99- dpprime_dr_b_r
         100- dpprime_dr_b_i
         101- dxi_r_dr_b_r
         102- dxi_r_dr_b_i
         103- dp0_dr_b
         104- d2p0_dr2_b
         105- V_div_xi_r_b_r
         106- V_div_xi_r_b_i
         107- xi_r_analytic_r
         108- xi_r_analytic_i
         109- V_analytic_r
         110- V_analytic_i
         111- Fprime_analytic_r
         112- Fprime_analytic_i
         113- dTprime_dT_r
         114- dTprime_dT_i
         115- Fprime_analytic_r T' part
         116- Fprime_analytic_i T' part
         117- Fprime_analytic_r p' part
         118- Fprime_analytic_i p' part
         119- Fprime_analytic_r gradT part
         120- Fprime_analytic_i gradT part
         121- flux_recovered  ->  radiative flux directly computed
         122- -K*dT0_dr
         123- rho*g
         124- deltaP_hydro_r  -> this uses hydrostatic equilibrium to use dp0_dr = -rho*g
         125- deltaP_hydro_i
         126- (m omega)^2 rho xi_h real
         127- (m omega)^2 rho xi_h imaginary
         128- other side of that equation ^ (real)
         129- other side of that equation ^ (imaginary)
         130- log10( mod(xi_radial) )
         131- log10( mod(F') )
         132- log10( mod(p') )
         133- log10( mod(T') )
         134- xi_r_r / xi_r_eq
         135- xi_r_i / xi_r_eq
         136- xi_h_real / xi_r_eq
         137- xi_h_imaginary / xi_r_eq
         138- xi_r_eq
         139- xi_r_i
         140- ddeltaP_dr_hydro_r
         141- ddeltaP_dr_hydro_i
         142- xi_r_analytic_hydro_r
         143- xi_r_analytic_hydro_i
         144- V_analytic_hydro_r
         145- V_analytic_hydro_i
         146- -l(l+1) K_0 T'
         147- manual derivative of r^2 F'
         148- K0
         149- DeltaF
         150- V_cont_r
         151- V_cont_i
         152- rho_prime_r / rho0
         153- rho_prime_i / rho0
         154- (1/r^2) d(r^2 xi_r)/dr
         155- (1/r^2) d(r^2 xi_i)/dr
         156- (xi_r_r / rho_0) * d(rho_0)/dr
         157- (xi_r_i / rho_0) * d(rho_0)/dr
         158- -l(l+1) * xi_{h,r} / r
         159- -l(l+1) * xi_{h,i} / r
         160- delta_rho_r
         161- delta_rho_i
         162- grav_HR_output
         163- H/R = L/(4 pi r^2 P omega R)   --->   thickness of non-adiabatic layer
         164- 163, but using total flux (approximately valid outside of energy generating regions)
         165- 163, but using the convective flux (approximately valid outside of energy generating regions)
         166- F_{rad} / F_{tot} (but F_{tot} is actually a fudge which is only valid once outside the sphere whithin which fusion occurs)
         167- s_prime_r
         168- s_prime_i
         169- delta_s_r
         170- delta_s_i
         171- conv_L_div_L
         172- F_tot
         173- F_conv
         174- F_conv_prime_div_F_conv_r
         175- F_conv_prime_div_F_conv_i
         176- prefactor (the prefactor of F'_conv)
         177- V_rho_r
         178- V_rho_i
         179- ds_prime_ds0_r
         180- ds_prime_ds0_i
         181- dr_dT_Fc_K0 = (dr/dT) * (F_conv / K0)
         182- x_conv
         183- chi_T
         184- chi_rho
         185- dkap_dlnrho_face_HR_output
         186- dkap_dlnT_face_HR_output
         187- opacity (kappa)
         188- s_prime_JP (from integrating)
         189- xi_r_JP
         190- xi_r_JP_term (to look at individual terms and see what is going on inside)
         191- F'_{conv} / F_{conv,0} (assuming v'_{c} = 0)
         192- s' from JP's equation (no integration)
         193- kappa' / kappa_{0}
         194- cp_HR_output
         195- F' von Ziepel's theorem
         196- drho_dr / rho
        
         
         */
        
        //	Fprime_analytic_r = dT0_dr*( Fprime_for_T[k]*vr[k][1][0]  +  ( - 4.0*7.5657e-15*2.99792458e10*temperature_HR_output[k]*temperature_HR_output[k]*temperature_HR_output[k]/( 3.0*opacity_HR_output[k]*rho_HR_output[k] ) )*(dTprime_dT_r  -  dxi_r_dr_r)  +  Fprime_for_p[k]*vr[k][0][0] );
        
        //	Fprime_analytic_i = dT0_dr*( Fprime_for_T[k]*vi[k][1][0]  +  ( - 4.0*7.5657e-15*2.99792458e10*temperature_HR_output[k]*temperature_HR_output[k]*temperature_HR_output[k]/( 3.0*opacity_HR_output[k]*rho_HR_output[k] ) )*(dTprime_dT_i  -  dxi_r_dr_i)  +  Fprime_for_p[k]*vi[k][0][0] );
        
        // 1 to 10
        outfile << radius_cm_HR_output[k]/R << "\t\t\t" << ur[k][0][0] << "\t\t" << ur[k][1][0] << "\t\t" << vr[k][0][0] << "\t\t" << vr[k][1][0] << "\t\t" << alphar[k][0][0] << "\t" << alphar[k][0][1] << "\t" << alphar[k][1][0] << "\t" << alphar[k][1][1] << "\t" << alphai[k][0][0] << "\t";
        
        // 11 to 20
        outfile << alphai[k][0][1] << "\t" << alphai[k][1][0] << "\t" << alphai[k][1][1] << "\t" << gammar[k][0][0] << "\t" << gammar[k][1][0] << "\t" << gammai[k][0][0] << "\t" << gammai[k][1][0] <<  "\t" << ur[k][0][0]*radius_cm_HR_output[k]*(1.0/100.0)*m*omega*((mp + Mstar)/(mp)) << "\t" << test[k] << "\t" << radius_cm_HR_output[k] << "\t";
        
        // 21 to 30
        outfile << flux_HR_output[k] << "\t" << pressure_HR_output[k] << "\t" << temperature_HR_output[k] << "\t" << rmid_cm_HR_output[k]/R << "\t\t\t" << R*ur[k][0][0] << "\t\t" << flux_BC*ur[k][1][0] << "\t\t" << pressure_HR_output[k]*vr[k][0][0] << "\t\t" << temperature_HR_output[k]*vr[k][1][0] << "\t\t\t" << R*ur[k][0][0]*m*omega*(mp + Mstar)/(100.0*mp) << "\t\t\t" << ui[k][0][0] << "\t\t";
        
        // 31 to 40
        outfile << ui[k][1][0] << "\t\t" << vi[k][0][0] << "\t\t" << vi[k][1][0] << "\t\t" << sqrt((ur[k][0][0]*ur[k][0][0]) + (ui[k][0][0]*ui[k][0][0])) << "\t\t" << atan(-ui[k][0][0] / ur[k][0][0]) << "\t\t" << m*omega*((mp + Mstar)/(100.0*mp))*((3.0)/(radius_cm_HR_output[k]*m*m*omega*omega))*( (pressure_HR_output[k]*vr[k][0][0]/rho_HR_output[k]) + f*radius_cm_HR_output[k]*radius_cm_HR_output[k] ) << "\t\t" << m*omega*((mp + Mstar)/(100.0*mp))*((3.0)/(radius_cm_HR_output[k]*m*m*omega*omega))*( (pressure_HR_output[k]*vi[k][0][0]/rho_HR_output[k]) ) << "\t\t" << m*omega*((mp + Mstar)/(100.0*mp))*((2.0)/(radius_cm_HR_output[k]*m*m*omega*omega))*( (pressure_HR_output[k]*vi[k][0][0]/rho_HR_output[k])) << "\t\t" << m*omega*((mp + Mstar)/(100.0*mp))*((2.0)/(radius_cm_HR_output[k]*m*m*omega*omega))*(-1.0)*( (pressure_HR_output[k]*vr[k][0][0]/rho_HR_output[k]) + f*radius_cm_HR_output[k]*radius_cm_HR_output[k] ) <<  "\t\t" << m*omega*((mp + Mstar)/(100.0*mp))*((1.0)/(radius_cm_HR_output[k]*m*m*omega*omega))*( (pressure_HR_output[k]*vr[k][0][0]/rho_HR_output[k]) + f*radius_cm_HR_output[k]*radius_cm_HR_output[k] ) << "\t\t";
        
        // 41 to 50
        outfile << m*omega*((mp + Mstar)/(100.0*mp))*((1.0)/(radius_cm_HR_output[k]*m*m*omega*omega))*( (pressure_HR_output[k]*vi[k][0][0]/rho_HR_output[k]) ) << "\t\t" << rho_HR_output[k] << "\t\t" << - f * rmid_cm_HR_output[k]*rmid_cm_HR_output[k]*rho_HR_output[k] << "\t\t" << f*radius_cm_HR_output[k]*radius_cm_HR_output[k]*rho_HR_output[k]*( (rmid_cm_HR_output[k+1] - rmid_cm_HR_output[k])/(pressure_HR_output[k+1] - pressure_HR_output[k]) ) << "\t\t" << pressure_HR_output[k]*vi[k][0][0] << "\t\t" << xi_h_r << "\t\t" << xi_h_pprime_part << "\t\t" << eta_Terquem_output[k] << "\t\t" << (radius_cm_HR_output[k]/H_rho)*R*ur[k][0][0] << "\t\t" << (radius_cm_HR_output[k]/H_rho)*R*ui[k][0][0] << "\t\t";
        
        // 51 to 60
        outfile << (radius_cm_HR_output[k]/H_rho)*mod_xi_radial << "\t\t" << xi_h_real << "\t\t" << xi_h_imaginary << "\t\t" << mod_xi_h << "\t\t" << H_p << "\t\t" << H_rho << "\t\t" << delta_P_r << "\t\t" << delta_P_i << "\t\t" << mod_delta_P << "\t\t" << V_div_xi_r_r << "\t\t";
        
        // 61 to 70
        outfile << V_div_xi_r_i << "\t\t" << otherVdiv_r << "\t\t" << otherVdiv_i << "\t\t" << pprime_comp_r << "\t\t" << pprime_comp_i << "\t\t" << Gvar_r << "\t\t" << Gvar_i << "\t\t" << gradient << "\t\t" << rho_face_HR_output[k] << "\t\t" << Hvar_r << "\t\t";
        
        // 71 to 80
        outfile << pprime_comp_second_r << "\t\t" << pprime_comp_second_i << "\t\t" << num_r << "\t\t" << num_i << "\t\t" << denom_r << "\t\t" << denom_i << "\t\t" << ddelta_P_dr_r << "\t\t" << dgrr_dr*radius_cm_HR_output[k]*radius_cm_HR_output[k]/grav_HR_output[k] << "\t\t" << grav_HR_output[k]/(m*m*omega*omega*radius_cm_HR_output[k]) << "\t\t" << ( (xi_h_real*V_div_xi_r_r) + (xi_h_imaginary*V_div_xi_r_i) )/( (V_div_xi_r_r*V_div_xi_r_r) + (V_div_xi_r_i*V_div_xi_r_i) ) << "\t\t";
        
        // 81 to 90
        outfile << ( (xi_h_imaginary*V_div_xi_r_r) - (xi_h_real*V_div_xi_r_i) )/( (V_div_xi_r_r*V_div_xi_r_r) + (V_div_xi_r_i*V_div_xi_r_i) ) << "\t\t" << dp0dr << "\t\t" << dp0dr_old << "\t\t" << dp0dr_avg << "\t\t" << delta_P_r_new << "\t\t" << delta_P_i_new << "\t\t" << -f*radius_cm_HR_output[k]*radius_cm_HR_output[k]/grav_HR_output[k] << "\t\t" << xi_h_over_xi_r_real << "\t\t" << xi_h_over_xi_r_im << "\t\t" << log_mod_xi_r << "\t\t";
        
        // 91 to 100
        outfile << log_mod_xi_h << "\t\t" << sqrt((xi_h_over_xi_r_real*xi_h_over_xi_r_real) + (xi_h_over_xi_r_im*xi_h_over_xi_r_im)) << "\t\t" << sqrt((V_div_xi_r_r*V_div_xi_r_r) + (V_div_xi_r_i*V_div_xi_r_i)) << "\t\t" << dp_dr << "\t\t" << d2p_dr2 << "\t\t" << d2p0_dr2 << "\t\t" << dDeltaP_dr_b_r << "\t\t" << dDeltaP_dr_b_i << "\t\t" << dpprime_dr_b_r << "\t\t" << dpprime_dr_b_i << "\t\t";
        
        // 101 to 110
        outfile << dxi_r_dr_b_r << "\t\t" << dxi_r_dr_b_i << "\t\t" << dp0_dr_b << "\t\t" << d2p0_dr2_b << "\t\t" << V_div_xi_r_b_r << "\t\t" << V_div_xi_r_b_i << "\t\t" << xi_r_analytic_r << "\t\t" << xi_r_analytic_i << "\t\t" << V_analytic_r << "\t\t" << V_analytic_i << "\t\t";
        
        // 111 to 120
        outfile << Fprime_analytic_r << "\t\t" << Fprime_analytic_i << "\t\t" << dTprime_dT_r << "\t\t" << dTprime_dT_i << "\t\t" << ( Fprime_for_T[k]*vr[k][1][0]) << "\t\t" << ( Fprime_for_T[k]*vi[k][1][0]) << "\t\t" << Fprime_for_p[k]*vr[k][0][0] << "\t\t" << Fprime_for_p[k]*vi[k][0][0] << "\t\t" << (flux_HR_output[k])*(dTprime_dT_r  -  0.0*dxi_r_dr_r) << "\t\t" << (flux_HR_output[k])*(dTprime_dT_i  -  dxi_r_dr_i) << "\t\t";
        
        // 121 to 130
        outfile << dT0_dr*( - 4.0*7.5657e-15*2.99792458e10*temperature_HR_output[k]*temperature_HR_output[k]*temperature_HR_output[k]/( 3.0*opacity_HR_output[k]*rho_HR_output[k] ) ) << "\t\t" << -K_HR_output[k]*dT0_dr << "\t\t" << -rho_HR_output[k]*grav_HR_output[k] << "\t\t" << deltaP_hydro_r << "\t\t" << deltaP_hydro_i << "\t\t" << m*m*omega*omega*rho_HR_output[k]*xi_h_real << "\t\t" << m*m*omega*omega*rho_HR_output[k]*xi_h_imaginary << "\t\t" << (deltaP_hydro_r/rmid_cm_HR_output[k]) + (grav_HR_output[k]*rho_HR_output[k]*R*0.5*(ur[k][0][0]+ur[k-1][0][0])/rmid_cm_HR_output[k]) + (rho_HR_output[k]*f*rmid_cm_HR_output[k]) << "\t\t" << (deltaP_hydro_i/rmid_cm_HR_output[k]) + (grav_HR_output[k]*rho_HR_output[k]*R*0.5*(ui[k][0][0]+ui[k-1][0][0])/rmid_cm_HR_output[k]) << "\t\t" << 0.5*log10(R*R*( ur[k][0][0]*ur[k][0][0] + ui[k][0][0]*ui[k][0][0]  )) << "\t\t";
        
        // 131 to 140
        outfile << 0.5*log10(flux_BC*flux_BC*( ur[k][1][0]*ur[k][1][0] + ui[k][1][0]*ui[k][1][0]  )) << "\t\t" << 0.5*log10( pressure_HR_output[k]*pressure_HR_output[k]*(vr[k][0][0]*vr[k][0][0] + vi[k][0][0]*vi[k][0][0]  )) << "\t\t" << 0.5*log10( temperature_HR_output[k]*temperature_HR_output[k]*(vr[k][1][0]*vr[k][1][0] + vi[k][1][0]*vi[k][1][0]  )) << "\t\t" << (R*ur[k][0][0])/xi_r_eq << "\t\t" << (R*ui[k][0][0])/xi_r_eq << "\t\t" << xi_h_real/xi_r_eq << "\t\t" << xi_h_imaginary/xi_r_eq << "\t\t" << xi_r_eq << "\t\t" << ui[k][0][0]*R << "\t\t" << ddeltaP_dr_hydro_r << "\t\t";
        
        // 141 to 150
        outfile << ddeltaP_dr_hydro_i << "\t\t" << xi_r_analytic_hydro_r << "\t\t" << xi_r_analytic_hydro_i << "\t\t" << V_analytic_hydro_r << "\t\t" << V_analytic_hydro_i << "\t\t" << -l*(l+1.0)*K_HR_output[k]*temperature_HR_output[k]*vr[k][1][0] << "\t\t" << (radius_cm_HR_output[k+1]*radius_cm_HR_output[k+1]*flux_BC*ur[k+1][1][0]  -  radius_cm_HR_output[k]*radius_cm_HR_output[k]*flux_BC*ur[k][1][0])/(radius_cm_HR_output[k+1]  -  radius_cm_HR_output[k]) << "\t\t" << K_HR_output[k] << "\t\t" << flux_BC*ur[k][1][0] + (R*ur[k][0][0])*((flux_tot_HR_output[k+1] - flux_tot_HR_output[k-1])/(radius_cm_HR_output[k+1]-radius_cm_HR_output[k-1])) << "\t\t" << V_cont_r << "\t\t";
        
        // 151 to 160
        outfile << V_cont_i << "\t\t" << (1.0/chiRho_HR_output[k])*( vr[k][0][0] - chiT_HR_output[k]*vr[k][1][0] ) << "\t\t" << (1.0/chiRho_HR_output[k])*( vi[k][0][0] - chiT_HR_output[k]*vi[k][1][0] ) << "\t\t" << (R/(radius_cm_HR_output[k]*radius_cm_HR_output[k]))*( (radius_cm_HR_output[k+1]*radius_cm_HR_output[k+1]*ur[k+1][0][0] - radius_cm_HR_output[k-1]*radius_cm_HR_output[k-1]*ur[k-1][0][0])/(radius_cm_HR_output[k+1] - radius_cm_HR_output[k-1]) ) << "\t\t" << (R/(radius_cm_HR_output[k]*radius_cm_HR_output[k]))*( (radius_cm_HR_output[k+1]*radius_cm_HR_output[k+1]*ui[k+1][0][0] - radius_cm_HR_output[k-1]*radius_cm_HR_output[k-1]*ui[k-1][0][0])/(radius_cm_HR_output[k+1] - radius_cm_HR_output[k-1]) ) << "\t\t" << (R*ur[k][0][0]/rho_HR_output[k])*( (rho_HR_output[k+1]-rho_HR_output[k])/(rmid_cm_HR_output[k+1]-rmid_cm_HR_output[k]) ) << "\t\t" << (R*ui[k][0][0]/rho_HR_output[k])*( (rho_HR_output[k+1]-rho_HR_output[k])/(rmid_cm_HR_output[k+1]-rmid_cm_HR_output[k]) ) << "\t\t" << -l*(l+1.0)*xi_h_real/radius_cm_HR_output[k] << "\t\t" << -l*(l+1.0)*xi_h_imaginary/radius_cm_HR_output[k] << "\t\t" << (rho_HR_output[k]/chiRho_HR_output[k])*( vr[k][0][0] - chiT_HR_output[k]*vr[k][1][0] ) + R*ur[k][0][0]*( (rho_HR_output[k+1]-rho_HR_output[k])/(rmid_cm_HR_output[k+1]-rmid_cm_HR_output[k]) ) << "\t\t";
        
        // 161 to 170
        outfile << (rho_HR_output[k]/chiRho_HR_output[k])*( vi[k][0][0] - chiT_HR_output[k]*vi[k][1][0] ) + R*ui[k][0][0]*( (rho_HR_output[k+1]-rho_HR_output[k])/(rmid_cm_HR_output[k+1]-rmid_cm_HR_output[k]) ) << "\t\t" << grav_HR_output[k] << "\t\t" << (flux_HR_output[k]/(omega*pressure_HR_output[k]))/R << "\t\t" << (flux_BC/(omega*pressure_HR_output[k]))/R << "\t\t" << ((flux_BC - flux_HR_output[k])/(omega*pressure_HR_output[k]))/R << "\t\t" << (flux_HR_output[k]/flux_BC) * ( (radius_cm_HR_output[k]*radius_cm_HR_output[k])/(R*R) ) << "\t\t" << s_prime_r << "\t\t" << s_prime_i << "\t\t" << delta_s_r << "\t\t" << delta_s_i << "\t\t";
        
        // 171 to 180
        outfile << conv_L_div_L_HR_output[k] << "\t\t" << flux_HR_output[k]/(1.0 - conv_L_div_L_HR_output[k]) << "\t\t" << flux_HR_output[k]*conv_L_div_L_HR_output[k]/(1.0 - conv_L_div_L_HR_output[k]) << "\t\t" << F_conv_prime_div_F_conv_r << "\t\t" << F_conv_prime_div_F_conv_i << "\t\t" << prefactor << "\t\t" << V_rho_r << "\t\t" << V_rho_i << "\t\t" << ds_prime_ds0_r << "\t\t" << ds_prime_ds0_i << "\t\t";
        
        // 181 to 190
        outfile << dr_dT_Fc_K0 << "\t\t" << x_conv << "\t\t" << chiT_HR_output[k] << "\t\t" << chiRho_HR_output[k] << "\t\t" << dkap_dlnrho_face_HR_output[k] << "\t\t" << dkap_dlnT_face_HR_output[k] << "\t\t" << opacity_HR_output[k] << "\t\t" << s_prime_JP << "\t\t" << xi_r_JP << "\t\t" << xi_r_JP_term << "\t\t";
        
        // 191 to 200
        outfile << ds_prime_ds0_r + vr[k][1][0]*( (1.0 - 2.0*x_conv)/(1.0 + x_conv) ) + (1.0/chiRho_HR_output[k])*( vr[k][0][0] - chiT_HR_output[k]*vr[k][1][0] )*( (1.0 + 3.0*x_conv)/(1.0 + x_conv) )  +  (   (dkap_dlnrho_face_HR_output[k]*(1.0/chiRho_HR_output[k])*( vr[k][0][0] - chiT_HR_output[k]*vr[k][1][0] )/opacity_HR_output[k])  +  (dkap_dlnrho_face_HR_output[k]*vr[k][1][0]/opacity_HR_output[k])   )*( (x_conv)/(1.0 + x_conv) ) << "\t\t" << s_prime_equation_JP << "\t\t" << (   ( (1.0/chiRho_HR_output[k])*( vr[k][0][0] - chiT_HR_output[k]*vr[k][1][0] ) * ( dkap_dlnrho_face_HR_output[k] / opacity_HR_output[k] )  +  vr[k][1][0]*( dkap_dlnT_face_HR_output[k] / opacity_HR_output[k] ) )   ) << "\t\t" << cp_HR_output[k] << "\t\t" << 2.0*f*radius_cm_HR_output[k]/grav_HR_output[k] << "\t\t" << (1.0/rho_HR_output[k])*( (rho_HR_output[k+1]-rho_HR_output[k])/(rmid_cm_HR_output[k+1]-rmid_cm_HR_output[k]) )  << "\n";
        
        
        
        
        xi_r_eq = - f * radius_cm_HR_output[k] * radius_cm_HR_output[k] * radius_cm_HR_output[k] * radius_cm_HR_output[k] / (G * Mstar) ; // mp*R*R*R*R/(4.0*Mstar*D*D*D);
        
        /*
         1 - log(pressure / ( G M^2 R^-4))
         2 - xi_r_real / R (scaled so that 1 at R = equilibrium tide)
         3 - xi_r_im / R (scaled so that 1 at R = equilibrium tide)
         4 - test
         5 - mod(xi_r)
         6 - mod(xi_h)
         7 - log10(mod(xi_r))
         8 - log10(mod(xi_h))
         9 - r / R
         10- log(mod(F_r'))
         11- log(mod(p'))
         12- log(mod(T'))
         13- mod(xi_r/R)
         14- mod(F'/flux_BC)
         15- log(mod(xi_r/xi_eq))
         16- log(mod(xi_h/xi_eq))
         17- mod(xi_h_cont)
         18- br
         19- bi
         
         */
        
        comparison_file << log10(pressure_HR_output[k]/((G*Mstar*Mstar/(R*R*R*R)))) << "\t\t" << log10( ur[k][0][0]*R*R / (xi_r_eq * radius_cm_HR_output[k]) ) << "\t\t" << log10( ui[k][0][0]*R*R / (xi_r_eq * radius_cm_HR_output[k]) ) << "\t\t" << test[k] << "\t\t" << sqrt( R*R*( ur[k][0][0]*ur[k][0][0] + ui[k][0][0]*ui[k][0][0] ) ) << "\t\t" << sqrt( xi_h_real*xi_h_real + xi_h_imaginary*xi_h_imaginary ) << "\t\t" << log10(sqrt( R*R*( ur[k][0][0]*ur[k][0][0] + ui[k][0][0]*ui[k][0][0] ) + 0.000000000001 )) << "\t\t" << log10(sqrt( xi_h_real*xi_h_real + xi_h_imaginary*xi_h_imaginary )) << "\t\t" << radius_cm_HR_output[k]/R << "\t\t" << log10(flux_BC*sqrt( ur[k][1][0]*ur[k][1][0]  +  ui[k][1][0]*ui[k][1][0] + 0.00000000000000000001)) << "\t\t" << log10( pressure_HR_output[k]*sqrt( vr[k][0][0]*vr[k][0][0]  +  vi[k][0][0]*vi[k][0][0]  ) ) << "\t\t" << log10( temperature_HR_output[k]*sqrt( vr[k][1][0]*vr[k][1][0]  +  vi[k][1][0]*vi[k][1][0]  )  ) << "\t\t" << sqrt( ur[k][0][0]*ur[k][0][0]  +  ui[k][0][0]*ui[k][0][0] ) << "\t\t" <<  sqrt( ur[k][1][0]*ur[k][1][0]  +  ui[k][1][0]*ui[k][1][0])  << "\t\t" << log10(sqrt(ur[k][0][0]*ur[k][0][0]  +  ui[k][0][0]*ui[k][0][0] + 0.000000000000000000001 )*R/xi_r_eq) << "\t\t" << log10(sqrt(xi_h_real*xi_h_real  +  xi_h_imaginary*xi_h_imaginary)/xi_r_eq) << "\t\t" << sqrt(V_cont_r*V_cont_r + V_cont_i*V_cont_i) << "\t\t" << ur[k][1][0] << "\t\t" << ui[k][1][0] << "\n";
        
        
        
        
        
    }
    
    // This closes the output data file
    outfile.close();
    
    // This closes the comparison file
    comparison_file.close();
    
    
    cout << "ur[J-1][0][0] =\t" << ur[J-1][0][0] << "\n";
    
    dumMAr[0][0][0] = 1.0;
    dumMAr[0][0][1] = 2.0;
    dumMAr[0][1][0] = 3.0;
    dumMAr[0][1][1] = 4.0;
    
    dumMAi[0][0][0] = 0.0;
    dumMAi[0][0][1] = 0.0;
    dumMAi[0][1][0] = 0.0;
    dumMAi[0][1][1] = 0.0;
    
    
    CMatrixDiagInv(dumMAr,dumMAi,dumMBr,dumMBi,0,0);
    
    
    
    
    
    
    sum = count[Jold-1] - Jold;
    
    //cout << "\n\n   total number of new cells needed = " << sum << "\n\n";
    
    cout << "\n\n   total number of cells = " << J << "\n\n";
    
    sum = 3.2 / 1.01;
    
    cout << "\n\n 3.2 / 1.01 = " << sum << "\n";
    
    //cout << "\n\ndx_max_surface_tracker = " << dx_max_surface_tracker << "\n\n";
    
    
    
    
    // This outputs the necessary data for the photometry file
    
    // This bit opens the file to write the data into
    ofstream photodata;
    photodata.open("Output/photo_data_Ubuntu.dat", ios::out);
    
    photodata.precision(15);
    
    /*
     This writes, in the following order:
     
     1 - f0 (equilibrium flux at the surface)
     2 - f_r (real part of perturbed flux)
     3 - f_i (imaginary part of perturbed flux)
     4 - df_dr (gradient of equilibrium flux with radius)
     5 - R (equilibrium surface radius)
     6 - xi_r (real displacement at the surface)
     7 - xi_i (imaginary displacement at the surface)
     
     */
    
    photodata <<  flux_HR_output[J-1] << "\n" << flux_BC*ur[J-1][1][0] << "\n" << flux_BC*ui[J-1][1][0] << "\n" << (flux_HR_output[J-1] - flux_HR_output[J-2])/(radius_cm_HR_output[J-1] - radius_cm_HR_output[J-2]) << "\n" << radius_cm_HR_output[J-1] << "\n" << R*ur[J-1][0][0] << "\n" << R*ui[J-1][0][0] << "\n";
    
    photodata.close();
    
    
    
    // This is written to create the input for the observables code
    // This depends upon the fact that we have just finished the loop writing everything to file, which runs from the centre to the surface, so the last time these variables were calculated they were at the surface. (Might want to check exactly this, just to be safe, though.)
    
    // This opens the appropriately named file
    
    string observables_filename;
    
    observables_filename = "Output/observables_input.dat";
    
    ofstream observables_file;
    observables_file.open("Output/observables_input.dat", ios::out);
    
    observables_file.precision(10);
    
    /*
     The outputs are:
     
     1 - f0 (equilibrium flux at the surface)
     2 - f_r (real part of perturbed flux)
     3 - f_i (imaginary part of perturbed flux)
     4 - df_dr (gradient of equilibrium flux with radius)
     5 - R (equilibrium surface radius)
     6 - xi_r (real displacement at the surface)
     7 - xi_i (imaginary displacement at the surface)
     8 - xi_h_re (V real)
     9 - xi_h_im (V imaginary)
     10- xi_r_eq (for comparison)
     11- omega
     
     */
    
    observables_file <<  flux_HR_output[J-1] << "\n" << flux_BC*ur[J-1][1][0] << "\n" << flux_BC*ui[J-1][1][0] << "\n" << (flux_HR_output[J-1] - flux_HR_output[J-2])/(radius_cm_HR_output[J-1] - radius_cm_HR_output[J-2]) << "\n" << radius_cm_HR_output[J-1] << "\n" << R*ur[J-1][0][0] << "\n" << R*ui[J-1][0][0] << "\n" << xi_h_real << "\n" << xi_h_imaginary << "\n" << xi_r_eq << "\n" << omega << "\n";
    
    
    
    cout << "The last value of lambda used = " << parameter_nonad << "\n";
    cout << "This was at x = " << location_nonad << "\n";
    
    
    
    
    
}
// main's  curly brackets have just been closed









/*
 Here are the definitions of the functions, in the order:
 CompMult - does complex multiplication
 CompDiv - does complex division
 MatrixMult - matrix multiplication
 MatrixInv - inverts a matrix
 VectorMult - multiplies a vector by a matrix
 CMatrixMult - does complex matrix multiplication
 CMatrixInv - does complex matrix inversion
 CVectorMult - does complex vector multiplication
 
 
 */

int CompMult(double* xr, double* xi, double* yr, double* yi, double* zr, double* zi)
{
    // x*y = z, each of which is split into the real and imaginary parts, as denoted by the second letter (r or i)
    
    (*zr) = (*xr)*(*yr) - (*xi)*(*yi);
    (*zi) = (*xi)*(*yr) + (*xr)*(*yi);
    
    return (0);
    
}

int CompDiv(double* xr, double* xi, double* yr, double* yi, double* zr, double* zi)
{
    // (x/y) = z, each of which is split into the real and imaginary parts, as denoted by the second letter (r or i)
    (*zr) = ((*xr)*(*yr) + (*xi)*(*yi))/(((*yr)*(*yr))+((*yi)*(*yi)));
    (*zi) = ((*xi)*(*yr) - (*xr)*(*yi))/(((*yr)*(*yr))+((*yi)*(*yi)));
    
    return (0);
    
}

int MatrixMult(double X[][2][2], double Y[][2][2], double Z[][2][2], int i, int j, int k)
{
    // i is the index of the matrices which are to be multipled - eg A_{i,i+1} B_{i,i+1}
    // i, j and k starts at 0, and go to J-1 in general, although this is not the case with the boundary condition matrices
    // i is the index of the first matrix, j of the second, and k of the output matrix
    // X and Y are the matrices to be multiplied, as XY = Z
    // Z is the output matrix
    
    
    Z[k][0][0] = X[i][0][0]*Y[j][0][0] + X[i][0][1]*Y[j][1][0];
    Z[k][0][1] = X[i][0][0]*Y[j][0][1] + X[i][0][1]*Y[j][1][1];
    Z[k][1][0] = X[i][1][0]*Y[j][0][0] + X[i][1][1]*Y[j][1][0];
    Z[k][1][1] = X[i][1][0]*Y[j][0][1] + X[i][1][1]*Y[j][1][1];
    
    return (0);
    
}




int MatrixInv(double X[][2][2], double Z[][2][2], int i)
{
    // i is the index of the matrices which are to be multipled - eg A_{i,i+1} B_{i,i+1}
    // i starts at 0, and goes to J-1
    // X is the matrix to be inverted, such that Z = X^(-1)
    // Z is the output matrix
    
    double det;
    
    det = X[i][0][0]*X[i][1][1] - X[i][0][1]*X[i][1][0] ;
    
    Z[i][0][0] = X[i][1][1]/det;
    Z[i][0][1] = -X[i][0][1]/det;
    Z[i][1][0] = -X[i][1][0]/det;
    Z[i][1][1] = X[i][0][0]/det;
    
    if (det == 0) {
        cout << "Determinant = 0 \n";
    }
    
    return (0);
    
}



int VectorMult(double X[][2][2], double Y[][2][1], double Z[][2][1], int i, int j, int k)
{
    // i is the index of the matrix which are to be multipled - eg A_{i,i+1}
    // j is the index of the vector to be multiplied - eg u_{j}
    // k is the index of the output vector - eg x_{k}
    // i, j and k generally start at 0, and go to J-1, but this is not the case when using boundary condition vectors
    // X and Y are the matrix and vector to be multiplied, as X_{i}Y_{j} = Z_{k}
    // Z is the output vector
    
    
    Z[k][0][0] = X[i][0][0]*Y[j][0][0] + X[i][0][1]*Y[j][1][0];
    Z[k][1][0] = X[i][1][0]*Y[j][0][0] + X[i][1][1]*Y[j][1][0];
    
    
    return (0);
    
}



int CMatrixMult(double Xr[][2][2], double Xi[][2][2] , double Yr[][2][2], double Yi[][2][2] , double Zr[][2][2], double Zi[][2][2], int i, int j, int k)
{
    // i is the index of the first matrix which is to be multipled - eg X_{i,i+1}
    // j is the index of the second matrix which is to be multipled - eg Y_{j,j+1}
    // k is the index of the output matrix - eg Z_{k,k+1}
    // i, j & k start at 0, and go to J-1
    // Z = XY
    // Z is the output matrix
    // The second letters on the matrix names refer to whether they contain the real or imaginary part
    
    double ar, ai, br, bi;
    
    //    cout << "The first matrix is:\n\n";
    //    cout << Xr[i][0][0] << " + " << Xi[i][0][0] << " i \t\t\t" << Xr[i][0][1] << " + " << Xi[i][0][1] << " i \n";
    //    cout << Xr[i][1][0] << " + " << Xi[i][1][0] << " i \t\t\t" << Xr[i][1][1] << " + " << Xi[i][1][1] << " i \n\n";
    
    //    cout << "The second matrix is:\n\n";
    //    cout << Yr[j][0][0] << " + " << Yi[j][0][0] << " i \t\t\t" << Yr[j][0][1] << " + " << Yi[j][0][1] << " i \n";
    //    cout << Yr[j][1][0] << " + " << Yi[j][1][0] << " i \t\t\t" << Yr[j][1][1] << " + " << Yi[j][1][1] << " i \n\n\n";
    
    CompMult(&Xr[i][0][0],&Xi[i][0][0],&Yr[j][0][0],&Yi[j][0][0],&ar,&ai);
    
    CompMult(&Xr[i][0][1],&Xi[i][0][1],&Yr[j][1][0],&Yi[j][1][0],&br,&bi);
    
    // This gives the 00 component
    Zr[k][0][0] = ar + br;
    Zi[k][0][0] = ai + bi;
    
    
    CompMult(&Xr[i][0][0],&Xi[i][0][0],&Yr[j][0][1],&Yi[j][0][1],&ar,&ai);
    
    CompMult(&Xr[i][0][1],&Xi[i][0][1],&Yr[j][1][1],&Yi[j][1][1],&br,&bi);
    
    // This gives the 01 component
    Zr[k][0][1] = ar + br;
    Zi[k][0][1] = ai + bi;
    
    
    CompMult(&Xr[i][1][0],&Xi[i][1][0],&Yr[j][0][0],&Yi[j][0][0],&ar,&ai);
    
    CompMult(&Xr[i][1][1],&Xi[i][1][1],&Yr[j][1][0],&Yi[j][1][0],&br,&bi);
    
    // This gives the 10 component
    Zr[k][1][0] = ar + br;
    Zi[k][1][0] = ai + bi;
    
    
    CompMult(&Xr[i][1][0],&Xi[i][1][0],&Yr[j][0][1],&Yi[j][0][1],&ar,&ai);
    
    CompMult(&Xr[i][1][1],&Xi[i][1][1],&Yr[j][1][1],&Yi[j][1][1],&br,&bi);
    
    // This gives the 11 component
    Zr[k][1][1] = ar + br;
    Zi[k][1][1] = ai + bi;
    
    
    
    return (0);
    
}



int CMatrixInv(double Xr[][2][2], double Xi[][2][2] , double Zr[][2][2], double Zi[][2][2], int i, int j)
{
    // i is the index of the matrix which is to be inverted - eg A_{i,i+1}
    // j is the index of the output matrix
    // i & j generally start at 0, and goes to J-1
    // X is the matrix to be inverted, such that Z = X^(-1)
    // Z is the output matrix
    // The second letters on the matrix names refer to whether they contain the real or imaginary part
    
    double detr, deti, ar, ai, br, bi, condnum, norm, norminv, absA, absB, absC, absD;
    
    //    cout << "The matrix to be inverted is:\n\n";
    //    cout << Xr[i][0][0] << " + " << Xi[i][0][0] << " i \t\t\t" << Xr[i][0][1] << " + " << Xi[i][0][1] << " i \n";
    //    cout << Xr[i][1][0] << " + " << Xi[i][1][0] << " i \t\t\t" << Xr[i][1][1] << " + " << Xi[i][1][1] << " i \n";
    
    CompMult(&Xr[i][0][0],&Xi[i][0][0],&Xr[i][1][1],&Xi[i][1][1],&ar,&ai);
    
    //    cout << "(" << Xr[i][0][0] << " + " << Xi[i][0][0] << " i) * (" << Xr[i][1][1] << " + " << Xi[i][1][1] << " i) = "<< ar << " + " << ai << " i \n";
    
    CompMult(&Xr[i][0][1],&Xi[i][0][1],&Xr[i][1][0],&Xi[i][1][0],&br,&bi);
    
    //    cout << "(" << Xr[i][0][1] << " + " << Xi[i][0][1] << " i) * (" << Xr[i][1][0] << " + " << Xi[i][1][0] << " i) = "<< br << " + " << bi << " i \n";
    
    detr = ar - br;
    deti = ai - bi;
    
    // This gives the 00 component
    CompDiv(&Xr[i][1][1],&Xi[i][1][1],&detr,&deti,&Zr[j][0][0],&Zi[j][0][0]);
    
    // This gives the 11 component
    CompDiv(&Xr[i][0][0],&Xi[i][0][0],&detr,&deti,&Zr[j][1][1],&Zi[j][1][1]);
    
    CompDiv(&Xr[i][0][1],&Xi[i][0][1],&detr,&deti,&Zr[j][0][1],&Zi[j][0][1]);
    
    CompDiv(&Xr[i][1][0],&Xi[i][1][0],&detr,&deti,&Zr[j][1][0],&Zi[j][1][0]);
    
    Zr[j][0][1] = (-1.0)*(Zr[j][0][1]);
    Zi[j][0][1] = (-1.0)*(Zi[j][0][1]);
    
    Zr[j][1][0] = (-1.0)*(Zr[j][1][0]);
    Zi[j][1][0] = (-1.0)*(Zi[j][1][0]);
    
    if (detr == 0 && deti == 0) {
        cout << "Determinant = 0 \n" << "Input matrix was:\n";
        cout << Xr[i][0][0] << " + " << Xi[i][0][0] << " i \t\t\t" << Xr[i][0][1] << " + " << Xi[i][0][1] << " i \n";
        cout << Xr[i][1][0] << " + " << Xi[i][1][0] << " i \t\t\t" << Xr[i][1][1] << " + " << Xi[i][1][1] << " i \n";
    }
    
    
    // This calculates the condition number
    
    absA = sqrt( (Xr[i][0][0]+Xi[i][0][0])*(Xr[i][0][0]-Xi[i][0][0]) );
    absB = sqrt( (Xr[i][0][1]+Xi[i][0][1])*(Xr[i][0][1]-Xi[i][0][1]) );
    absC = sqrt( (Xr[i][1][0]+Xi[i][1][0])*(Xr[i][1][0]-Xi[i][1][0]) );
    absD = sqrt( (Xr[i][1][1]+Xi[i][1][1])*(Xr[i][1][1]-Xi[i][1][1]) );
    
    
    norm = max( absA + absC , absB + absD );
    
    
    absA = sqrt( (Zr[j][0][0]+Zi[j][0][0])*(Zr[j][0][0]-Zi[j][0][0]) );
    absB = sqrt( (Zr[j][0][1]+Zi[j][0][1])*(Zr[j][0][1]-Zi[j][0][1]) );
    absC = sqrt( (Zr[j][1][0]+Zi[j][1][0])*(Zr[j][1][0]-Zi[j][1][0]) );
    absD = sqrt( (Zr[j][1][1]+Zi[j][1][1])*(Zr[j][1][1]-Zi[j][1][1]) );
    
    
    norminv = max( absA + absC , absB + absD );
    
    condnum = norm*norminv;
    
    if ( log10(condnum) > 20 && detr < 0.0000000000001) {
        
        cout << "\n\nLarge condition number:\t\t"  << condnum << "\n" ;
        cout << "Small determinant:\t\t"  << detr << "\n" ;
        // CMatrixDiagInv(Xr,Xi,Zr,Zi,i,j);
        
        
        cout << "\nThe matrix to be inverted is:\n";
        cout << Xr[i][0][0] << " + " << Xi[i][0][0] << " i \t\t\t" << Xr[i][0][1] << " + " << Xi[i][0][1] << " i \n";
        cout << Xr[i][1][0] << " + " << Xi[i][1][0] << " i \t\t\t" << Xr[i][1][1] << " + " << Xi[i][1][1] << " i \n";
        
        cout << "\nThe inverted matrix is:\n";
        cout << Zr[j][0][0] << " + " << Zi[j][0][0] << " i \t\t\t" << Zr[j][0][1] << " + " << Zi[j][0][1] << " i \n";
        cout << Zr[j][1][0] << " + " << Zi[j][1][0] << " i \t\t\t" << Zr[j][1][1] << " + " << Zi[j][1][1] << " i \n";
        
        
        
        
        
        //cout << "The inverted matrix using the diagonalisation is:\n\n";
        //cout << Zr[i][0][0] << " + " << Zi[i][0][0] << " i \t\t\t" << Zr[i][0][1] << " + " << Zi[i][0][1] << " i \n";
        //cout << Zr[i][1][0] << " + " << Zi[i][1][0] << " i \t\t\t" << Zr[i][1][1] << " + " << Zi[i][1][1] << " i \n";
        
        
        
    }
    
    
    
    
    
    return (0);
    
}



int CVectorMult(double Xr[][2][2], double Xi[][2][2] , double Yr[][2][1], double Yi[][2][1] , double Zr[][2][1], double Zi[][2][1], int i, int j, int k)
{
    // i is the index of the matrix which is to be multipled - eg X_{i,i+1}
    // j is the index of the vector which is to be multipled - eg Y_{j,j+1}
    // k is the index of the output vector - eg Z_{k,k+1}
    // Z = XY
    // Z is the output vector
    // The second letters on the matrix names refer to whether they contain the real or imaginary part
    
    double ar, ai, br, bi;
    
    //    cout << "The matrix is:\n\n";
    //    cout << Xr[i][0][0] << " + " << Xi[i][0][0] << " i \t\t\t" << Xr[i][0][1] << " + " << Xi[i][0][1] << " i \n";
    //    cout << Xr[i][1][0] << " + " << Xi[i][1][0] << " i \t\t\t" << Xr[i][1][1] << " + " << Xi[i][1][1] << " i \n\n";
    
    //    cout << "The vector is:\n\n";
    //    cout << Yr[j][0][0] << " + " << Yi[j][0][0] << " i \n";
    //    cout << Yr[j][1][0] << " + " << Yi[j][1][0] << " i \n\n\n";
    
    CompMult(&Xr[i][0][0],&Xi[i][0][0],&Yr[j][0][0],&Yi[j][0][0],&ar,&ai);
    
    CompMult(&Xr[i][0][1],&Xi[i][0][1],&Yr[j][1][0],&Yi[j][1][0],&br,&bi);
    
    // This gives the 00 component
    Zr[k][0][0] = ar + br;
    Zi[k][0][0] = ai + bi;
    
    
    
    CompMult(&Xr[i][1][0],&Xi[i][1][0],&Yr[j][0][0],&Yi[j][0][0],&ar,&ai);
    
    CompMult(&Xr[i][1][1],&Xi[i][1][1],&Yr[j][1][0],&Yi[j][1][0],&br,&bi);
    
    // This gives the 10 component
    Zr[k][1][0] = ar + br;
    Zi[k][1][0] = ai + bi;
    
    
    
    
    return (0);
    
}






int CMatrixDiagInv(double Xr[][2][2], double Xi[][2][2] , double Zr[][2][2], double Zi[][2][2], int i, int j)
{
    // i is the index of the matrix which is to be inverted - eg A_{i,i+1}
    // j is the index of the output matrix
    // i & j generally start at 0, and goes to J-1
    // X is the matrix to be inverted, such that Z = X^(-1)
    // Z is the output matrix
    // The second letters on the matrix names refer to whether they contain the real or imaginary part
    // This method diagonalises the matrix, inverts it, then transforms it back
    
    double detr, deti, ar, ai, br, bi, condnum, norm, norminv, absA, absB, absC, absD;
    double lambdaPr,lambdaPi,lambdaMr,lambdaMi; // for lambda+ and lambda-, the eigenvalues
    double ur, ui, vr, vi, xr, xi, yr, yi, zr, zi,AP,AM; // These are dummy variables
    double MRr[1][2][2], MRi[1][2][2], MLr[1][2][2], MLi[1][2][2], MRInvr[1][2][2], MRInvi[1][2][2], MLInvr[1][2][2], MLInvi[1][2][2]; // These are the right and left matrices and their inverses
    double Dr[1][2][2], Di[1][2][2]; // This is the diagonal matrix, which is then inverted
    double dumMAr[1][2][2], dumMAi[1][2][2];  // This is a dummy matrix
    
    
    cout << "\nThe matrix to be inverted is:\n";
    cout << Xr[i][0][0] << " + " << Xi[i][0][0] << " i \t\t\t" << Xr[i][0][1] << " + " << Xi[i][0][1] << " i \n";
    cout << Xr[i][1][0] << " + " << Xi[i][1][0] << " i \t\t\t" << Xr[i][1][1] << " + " << Xi[i][1][1] << " i \n\n";
    
    CompMult(&Xr[i][0][0],&Xi[i][0][0],&Xr[i][1][1],&Xi[i][1][1],&ar,&ai);
    
    //    cout << "(" << Xr[i][0][0] << " + " << Xi[i][0][0] << " i) * (" << Xr[i][1][1] << " + " << Xi[i][1][1] << " i) = "<< ar << " + " << ai << " i \n";
    
    CompMult(&Xr[i][0][1],&Xi[i][0][1],&Xr[i][1][0],&Xi[i][1][0],&br,&bi);
    
    //    cout << "(" << Xr[i][0][1] << " + " << Xi[i][0][1] << " i) * (" << Xr[i][1][0] << " + " << Xi[i][1][0] << " i) = "<< br << " + " << bi << " i \n";
    
    detr = ar - br;
    deti = ai - bi;
    
    // lambda = (a+d)/2.0 +- sqrt( ((a+d)/2.0)^2 - det )
    
    ur = 0.5*( Xr[i][0][0]+Xr[i][1][1] );
    ui = 0.5*( Xi[i][0][0]+Xi[i][1][1] );
    
    CompMult(&ur,&ui,&ur,&ui,&vr,&vi);
    
    xr = vr - detr;
    xi = vi - deti;
    
    CompSQRT(&xr,&xi,&yr,&yi);
    
    lambdaPr = ur + yr;
    lambdaPi = ui + yi;
    
    lambdaMr = ur - yr;
    lambdaMi = ui - yi;
    
    cout << "lambdaP = \t" << lambdaPr << "\t+\t" << lambdaPi << " i\n";
    
    cout << "lambdaM = \t" << lambdaMr << "\t+\t" << lambdaMi << " i\n";
    
    //Now we calculate the right matrix, MR
    ur = lambdaPr - Xr[i][0][0];
    ui = lambdaPi - Xi[i][0][0];
    
    AP = sqrt( ( Xr[i][0][1]*Xr[i][0][1] + Xi[i][0][1]*Xi[i][0][1] )/( ( Xr[i][0][1]*Xr[i][0][1] + Xi[i][0][1]*Xi[i][0][1] ) + ( ur*ur + ui*ui ) ) );
    
    cout << "AP = \t" << AP << "\n";
    
    ur = lambdaMr - Xr[i][0][0];
    ui = lambdaMi - Xi[i][0][0];
    
    AM = sqrt( ( Xr[i][0][1]*Xr[i][0][1] + Xi[i][0][1]*Xi[i][0][1] )/( ( Xr[i][0][1]*Xr[i][0][1] + Xi[i][0][1]*Xi[i][0][1] ) + ( ur*ur + ui*ui ) ) );
    
    cout << "AM = \t" << AM << "\n";
    
    
    ur = lambdaPr - Xr[i][0][0];
    ui = lambdaPi - Xi[i][0][0];
    
    CompDiv(&ur,&ui,&Xr[i][0][1],&Xi[i][0][1],&vr,&vi);
    
    MRr[0][0][0] = AP;
    
    MRi[0][0][0] = 0.0;
    
    MRr[0][1][0] = AP*vr;
    
    MRi[0][1][0] = AP*vi;
    
    
    ur = lambdaMr - Xr[i][0][0];
    ui = lambdaMi - Xi[i][0][0];
    
    CompDiv(&ur,&ui,&Xr[i][0][1],&Xi[i][0][1],&vr,&vi);
    
    MRr[0][0][1] = AM;
    
    MRi[0][0][1] = 0.0;
    
    MRr[0][1][1] = AM*vr;
    
    MRi[0][1][1] = AM*vi;
    
    cout << "The right matrix is:\n";
    
    cout << MRr[0][0][0] << " + " << MRi[0][0][0] << " i \t" << MRr[0][0][1] << " + " << MRi[0][0][1] << " i \n";
    cout << MRr[0][1][0] << " + " << MRi[0][1][0] << " i \t" << MRr[0][1][1] << " + " << MRi[0][1][1] << " i \n";
    
    CMatrixInv(MRr,MRi,MRInvr,MRInvi,0,0);
    
    
    
    
    //Now we calculate the left matrix, ML
    ur = lambdaPr - Xr[i][0][0];
    ui = lambdaPi - Xi[i][0][0];
    
    AP = sqrt( ( Xr[i][1][0]*Xr[i][1][0] + Xi[i][0][1]*Xi[i][0][1] )/( ( Xr[i][1][0]*Xr[i][1][0] + Xi[i][1][0]*Xi[i][1][0] ) + ( ur*ur + ui*ui ) ) );
    
    cout << "AP = \t" << AP << "\n";
    
    ur = lambdaMr - Xr[i][0][0];
    ui = lambdaMi - Xi[i][0][0];
    
    AM = sqrt( ( Xr[i][1][0]*Xr[i][1][0] + Xi[i][0][1]*Xi[i][0][1] )/( ( Xr[i][1][0]*Xr[i][1][0] + Xi[i][1][0]*Xi[i][1][0] ) + ( ur*ur + ui*ui ) ) );
    
    cout << "AM = \t" << AM << "\n";
    
    
    ur = lambdaPr - Xr[i][0][0];
    ui = lambdaPi - Xi[i][0][0];
    
    CompDiv(&ur,&ui,&Xr[i][1][0],&Xi[i][1][0],&vr,&vi);
    
    
    MLr[0][0][0] = AP;
    
    MLi[0][0][0] = 0.0;
    
    MLr[0][0][1] = AP*vr;
    
    MLi[0][0][1] = AP*vi;
    
    
    ur = lambdaMr - Xr[i][0][0];
    ui = lambdaMi - Xi[i][0][0];
    
    CompDiv(&ur,&ui,&Xr[i][1][0],&Xi[i][1][0],&vr,&vi);
    
    MLr[0][1][0] = AM;
    
    MLi[0][1][0] = 0.0;
    
    MLr[0][1][1] = AM*vr;
    
    MLi[0][1][1] = AM*vi;
    
    cout << "The left matrix is:\n";
    
    cout << MLr[0][0][0] << " + " << MLi[0][0][0] << " i \t" << MLr[0][0][1] << " + " << MLi[0][0][1] << " i \n";
    cout << MLr[0][1][0] << " + " << MLi[0][1][0] << " i \t" << MLr[0][1][1] << " + " << MLi[0][1][1] << " i \n";
    
    CMatrixInv(MLr,MLi,MLInvr,MLInvi,0,0);
    
    
    
    
    Dr[0][0][0] = lambdaPr;
    Dr[0][0][1] = 0.0;
    Dr[0][1][0] = 0.0;
    Dr[0][1][1] = lambdaMr;
    
    Di[0][0][0] = lambdaPi;
    Di[0][0][1] = 0.0;
    Di[0][1][0] = 0.0;
    Di[0][1][1] = lambdaMi;
    
    vr = 1.0;
    vi = 0.0;
    
    CompDiv(&vr,&vi,&lambdaPr,&lambdaPi,&Dr[0][0][0],&Di[0][0][0]);
    
    CompDiv(&vr,&vi,&lambdaMr,&lambdaMi,&Dr[0][1][1],&Di[0][1][1]);
    
    CMatrixMult(Dr,Di,MLr,MLi,dumMAr,dumMAi,0,0,0);
    
    CMatrixMult(MLInvr,MLInvi,dumMAr,dumMAi,Zr,Zi,0,0,j);
    
    
    cout << "\n\nThe left-inverted matrix is:\n";
    cout << Zr[j][0][0] << " + " << Zi[j][0][0] << " i \t" << Zr[j][0][1] << " + " << Zi[j][0][1] << " i \n";
    cout << Zr[j][1][0] << " + " << Zi[j][1][0] << " i \t" << Zr[j][1][1] << " + " << Zi[j][1][1] << " i \n";
    
    
    
    
    
    CMatrixMult(Dr,Di,MRInvr,MRInvi,dumMAr,dumMAi,0,0,0);
    
    CMatrixMult(MRr,MRi,dumMAr,dumMAi,Zr,Zi,0,0,j);
    
    
    cout << "\n\nThe right-inverted matrix is:\n";
    cout << Zr[j][0][0] << " + " << Zi[j][0][0] << " i \t" << Zr[j][0][1] << " + " << Zi[j][0][1] << " i \n";
    cout << Zr[j][1][0] << " + " << Zi[j][1][0] << " i \t" << Zr[j][1][1] << " + " << Zi[j][1][1] << " i \n";
    
    
    return (0);
    
}







int CompSQRT(double* ar, double* ai , double* br, double* bi)
{
    // This outputs the positive square root of the complex input, a
    
    double r, theta, rroot;
    
    r = sqrt( ((*ar)*(*ar)) + ((*ai)*(*ai)) );
    
    
    rroot = sqrt(r);
    
    
    theta = atan( (*ai)/(*ar) );
    
    if ((*ar) < 0.0) {
        theta = theta + 3.141592653589793238462643383279502884;
    }
    
    
    
    *br = rroot * cos( 0.5*theta );
    
    *bi = rroot * sin( 0.5*theta );
    
    
    
    return (0);
    
}










int FifthOrderExtrap(double var[], double r[], double C)
{
    /*
     This extrapolates using the following setup:
     
     dq = r_{q} - r_{x}
     where q is a stand-in for a to e, and r_{q} is the radius evaluated at the position of variable q, and r_{x} is the radius at the position of x, which we seek to extrapolate to get.
     
     var is the variable which we are seeking to extrapolate.
     
     To extrapolate for the centremost cell, we will use the known values of var[1] to var[5].
     
     r is to differentiate between the variables evaluated at either radius_cm or rmid_cm.
     
     C is a measure of how much closer the extrapolated cell will be to the centre: r[0] = r[1] / C
     
     */
    
    
    double a, b, c, d, e, f, g, h, da, db, dc, dd, de, df, dg, dh;
    
    r[0] = r[1]/C;
    
    a = var[1];
    b = var[2];
    c = var[3];
    d = var[4];
    e = var[5];
    f = var[6];
    g = var[7];
    h = var[8];
    
    
    da = r[1] - r[0];
    db = r[2] - r[0];
    dc = r[3] - r[0];
    dd = r[4] - r[0];
    de = r[5] - r[0];
    df = r[6] - r[0];
    dg = r[7] - r[0];
    dh = r[8] - r[0];
    
    
    // This does the fourth order approximation
    
    var[0] = ( ( a*db*dc*dd*de )/( (da - db)*(da - dc)*(da - dd)*(da - de) ) )  +  ( ( da*b*dc*dd*de )/( (db - da)*(db - dc)*(db - dd)*(db - de) ) )  +  ( ( da*db*c*dd*de )/( (dc - da)*(dc - db)*(dc - dd)*(dc - de) ) )  +  ( ( da*db*dc*d*de )/( (dd - da)*(dd - db)*(dd - dc)*(dd - de) ) )  +  ( ( da*db*dc*dd*e )/( (de - da)*(de - db)*(de - dc)*(de - dd) ) );
    
    
    // This does the second order approximation
    
    var[0] = (  (a*db*dc)/((da - db)*(da - dc))  )  +  (  (da*b*dc)/((db - da)*(db - dc))  )  +  (  (da*db*c)/((dc - da)*(dc - db))  );
    
    // This does the seventh order approximation
    
    // var[0] = ( ( a*db*dc*dd*de*df*dg*dh )/( (da - db)*(da - dc)*(da - dd)*(da - de)*(da - df)*(da - dg)*(da - dh) ) )  +  ( ( da*b*dc*dd*de*df*dg*dh )/( (db - da)*(db - dc)*(db - dd)*(db - de)*(db - df)*(db - dg)*(db - dh) ) )  +  ( ( da*db*c*dd*de*df*dg*dh )/( (dc - da)*(dc - db)*(dc - dd)*(dc - de)*(dc - df)*(dc - dg)*(dc - dh) ) )  +  ( ( da*db*dc*d*de*df*dg*dh )/( (dd - da)*(dd - db)*(dd - dc)*(dd - de)*(dd - df)*(dd - dg)*(dd - dh) ) )  +  ( ( da*db*dc*dd*e*df*dg*dh )/( (de - da)*(de - db)*(de - dc)*(de - dd)*(de - df)*(de - dg)*(de - dh) ) )  +  ( ( da*db*dc*dd*de*f*dg*dh )/( (df - da)*(df - db)*(df - dc)*(df - dd)*(df - de)*(df - dg)*(df - dh) ) )  +  ( ( da*db*dc*dd*de*df*g*dh )/( (dg - da)*(dg - db)*(dg - dc)*(dg - dd)*(dg - de)*(dg - df)*(dg - dh) ) )  +  ( ( da*db*dc*dd*de*df*dg*h )/( (dh - da)*(dh - db)*(dh - dc)*(dh - dd)*(dh - de)*(dh - df)*(dh - dg) ) );
    
    
    // This does a first order approximation (that is, linear)
    
    var[0] = (a*db)/(db - da)  -  (da*b)/(db - da);
    
    
    
    return(0);
    
    
    
}




int FunctionF(double* xadd, double* f)
{
    // At the value of x, this outputs the function defined within to variable f
    
    double x;
    
    x = *xadd;
    
    *f = 0.00006 + 0.00013*x*x; //0.00006 + 0.00013*x*x; //0.01+0.03*x*x; //0.00005 + 0.00015*x*x;
    
    return (0);
}


int FunctionG(double* xadd, double* g)
{
    // At the value of x, this outputs the function defined within to variable g
    
    double x;
    
    x = *xadd;
    
    *g = 0.0000035; //(0.0001 + 0.00015*x*x)*( 1.0 - exp(-(x-0.4999)/0.3) ); // 0.000000027; //0.0000035;
    
    return (0);
}


int FunctionY(double* xadd, double* y)
{
    // At the value of x, this outputs the function defined within to variable y
    
    double f, g, beta, x0, x, near;
    
    x = *xadd;
    
    // This defines where the transition from f to g takes place
    x0 = 0.985; //0.985; //0.98001; //0.9995; //0.98;
    
    // This defines the scale over which the change takes place
    //BEWARE!! If beta becomes too small, then we get issues with the exponentials
    // BUT!! That has now been sorted by manually avoiding the exponentials once far from the transition point
    beta = 0.0005; //0.0005; //0.0000005; //0.0005;
    
    near = std::abs((x-x0)/beta);
    
    FunctionF(&x,&f);
    
    FunctionG(&x,&g);
    
    // We only want to include the exponentials if we are near the transition point,
    // otherwise we just assert what the value is, as either f or g, as appropriate
    if (near < 30 ) {
        
        *y = ( (f*exp(-(x-x0)/beta))/(1.0 + exp(-(x-x0)/beta)) ) + ( (g*exp((x-x0)/beta))/(1.0 + exp((x-x0)/beta)) );
        
    } else {
        
        if ( x < x0 ) {
            
            *y = f;
            
        } else {
            
            *y = g;
            
        }
        
    }
    
    return (0);
}




int MeasureGrid(int* J_new_add, double* ratio_max_add)
{
    // grid_old[] is the array containing the locations of all of the previous cell outer edges
    // J_old is the previous number of cells
    // J_new is the new total number of cells
    
    // y is the guide for cell width
    // f and g are the function which make up y
    
    int alpha, j, k, too_big, too_small;
    
    double x, y, f, g, dx, dx_old, ratio_max;
    
    ratio_max = *ratio_max_add;
    
    
    
    
    
    // file with all of the sample locations
    ofstream outfile;
    
    outfile.open("Output/Interpolation.dat", ios::out);
    
    outfile.precision(10);
    
    // This tells us that we work from the outside edge back in
    x=1.0;
    
    outfile << "x" << "\t\t\t\t" << "f" << "\t\t\t\t" << "g" << "\t\t\t\t" << "y" << "\t\t\t\t" << "dx" << "\n";
    
    // These keep track of, respectively:
    // how many cells have been created so far
    // how many times the next cell width wanted to exceed ratio_max
    // how many times the next cell width wanted to go smaller than 1/ratio_max
    j = 0;
    too_big = 0;
    too_small = 0;
    
    
    while (x > 0.0) {
        
        FunctionF(&x, &f);
        
        FunctionG(&x, &g);
        
        FunctionY(&x, &y);
        
        dx = y;
        
        // This initialises dx_old for the first time through
        if( j == 0 ) {
            dx_old = dx;
        }
        
        if( dx > (ratio_max*dx_old) ) {
            dx = ratio_max*dx_old;
            too_big = too_big + 1;
        }
        
        if( dx < (dx_old/ratio_max) ) {
            dx = dx_old/ratio_max;
            too_small = too_small + 1;
        }
        
        outfile << x << "\t\t\t" << f << "\t\t\t" << g << "\t\t\t" << y << "\t\t\t" << dx << "\n";
        
        x = x - dx;
        j = j + 1;
        dx_old = dx;
        
    }
    
    outfile.close();
    
    
    cout << "j = " << j << "\n";
    cout << "too_big = " << too_big << "\n";
    cout << "too_small = " << too_small << "\n";
    
    *J_new_add = j;
    
    
    
    
    
    return (0);
    
}




int MakeGrid(int J_new, double ratio_max)
{
    // grid_old[] is the array containing the locations of all of the previous cell outer edges
    // J_old is the previous number of cells
    // The locations are written to a file, and then inverted so they can be read from the centre outwards
    // J_new is the new total number of cells
    
    // y is the guide for cell width
    // f and g are the function which make up y
    
    int j, k, too_big, too_small;
    
    double x, x0, y, dx, dx_old;
    
    //cout << "flag a \n";
    
    // This tells us that we work from the outside edge back in
    x=1.0;
    
    
    // These keep track of, respectively:
    // how many times the next cell width wanted to exceed ratio_max
    // how many times the next cell width wanted to go smaller than 1/ratio_max
    too_big = 0;
    too_small = 0;
    
    j = 0;
    
    k = J_new-1;
    
    std::ofstream gridfile;
    gridfile.open("Memory/grid_new_surface.txt");
    
    gridfile.precision(12);
    
    
    //cout << "flag b \n\n";
    
    while (x > 0.0) {
        
        
        gridfile << x << "\n";
        
        FunctionY(&x, &y);
        
        dx = y;
        
        // This initialises dx_old for the first time through
        if( j == 0 ) {
            dx_old = dx;
        }
        
        if( dx > (ratio_max*dx_old) ) {
            dx = ratio_max*dx_old;
            too_big = too_big + 1;
        }
        
        if( dx < (dx_old/ratio_max) ) {
            dx = dx_old/ratio_max;
            too_small = too_small + 1;
        }
        
        
        //cout << grid_new[k] << ", " << y << ", " << dx << ", " << too_small << "; ";
        
        x = x - dx;
        j = j + 1;
        k = k - 1;
        dx_old = dx;
        
        if (j%10000 == 0) {
            
            gridfile.close();
            
            gridfile.open("Memory/grid_new_surface.txt", std::ofstream::app);
            
            gridfile.precision(12);
            
        }
        
        
    }
    
    cout << "At the end of using MakeGrid, j = " << j << "\n";
    cout << "and k = " << k << "\n";
    cout << "and the last x value used was " << x+dx_old << "\n";
    
    
    //cout << "\nflag c \n";
    
    gridfile.close();
    
    // This inverts the file, so that it goes from the centre outwards, instead of from the surface inwards
    system("awk '{a[i++]=$0} END {for (j=i-1; j>=0;) print a[j--] > \"Memory/grid_new_centre.txt\" }' Memory/grid_new_surface.txt");
    
    //system("tail -r \"Memory/grid_new_surface.txt\" > \"Memory/grid_new_centre.txt\"");
    
    /*
     string instruction1, instruction2, full_instruction, J_new_string;
     
     instruction1 = "tail -r -n ";
     
     instruction2 = " \"Memory/grid_new_surface.txt\" > \"Memory/grid_new_centre.txt\"";
     
     J_new_string = to_string(J_new);
     
     full_instruction = instruction1+J_new_string+instruction2;
     
     system(full_instruction.c_str());
     
     */
    
    
    
    
    
    return (0);
    
}










