#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


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


int main()
{

    int J,k,z;
    double C;


    // J must be defined and given a value before any arrays which need it are defined, or else you'll get a segmentation fault because the arrays won't know how big they are.
    // So J is given a temporary value here until some MESA data is actually read in and use to define the size of the vectors instead.

    J=17500; // At the moment, it seems that this produces a seg fault when J >= 17500, which seems to be a result of filling up a memory limit, potentially the RAM?  But that would seem unlikely... It was helped by changing the dummy matrices to be [1][2][2] instead of [J][2][2], so it is a total memory issue rather than the memory used by any given array.
    z=J-1;


    // This section is to do with reading input from a file

    ifstream infile;
    infile.open("Input/profiles_Henyey_sun2.txt");

    k=0;

    double input;
    int no_of_lines;
    string line;

    no_of_lines = 0;

    ifstream linefile;
    linefile.open("Input/profiles_Henyey_sun2.txt");

    getline(linefile,line);
    while (linefile) {
        if (line != "") {
            no_of_lines = no_of_lines + 1;
        }

        getline(linefile,line);

    }

    cout << "The number of non-empty lines in this file is: " << no_of_lines << "\n" ;
    
    

    
    J = no_of_lines;
    
    
    
    
    z = J-1;


    // Here the input arrays are defined, and J cannot be changed again

    double zone[J], lnT[J], lnRho[J], grav[J], radius_cm[J], rmid_cm[J];
    double temperature[J], rho[J], pressure[J], grada[J], cp[J], chiRho[J], chiT[J];
    double opacity[J], dkap_dlnrho_face[J], dkap_dlnT_face[J], flux[J], brunt_A[J];
    double K[J], rho_face[J];

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

            // lnP[z] = input*2.30258509299404568401799; // This has been commented out as it is not being used, but is taking up quite a lot of memory

            infile >> input;

            grav[z] = input;

            infile >> input;

            radius_cm[z] = input;

            infile >> input;

            rmid_cm[z] = input;

            infile >> input;

            //dr[z] = input; // This has been commented out as it is not being used, but is taking up quite a lot of memory

            infile >> input;

            temperature[z] = input;

            infile >> input;

            rho[z] = input;

            infile >> input;

            //entropy[z] = input; // This has been commented out as it is not being used, but is taking up quite a lot of memory

            infile >> input;

            pressure[z] = input;

            infile >> input;

            grada[z] = input;

            infile >> input;

            cp[z] = input;

            infile >> input;

            chiRho[z] = input;

            infile >> input;

            chiT[z] = input;

            infile >> input;

        opacity[z] = input*10.0;

            infile >> input;

            dkap_dlnrho_face[z] = input;

            infile >> input;

            dkap_dlnT_face[z] = input;

            infile >> input;

            flux[z] = input;

            infile >> input;

            //brunt_N2[z] = input; // This has been commented out as it is not being used, but is taking up quite a lot of memory

            infile >> input;

            brunt_A[z] = input;

            infile >> input;

            //brunt_N[z] = input; // This has been commented out as it is not being used, but is taking up quite a lot of memory

            infile >> input;


        K[z] = 4.0 * 7.565767e-15 * 2.99792458e10 * temperature[z]*temperature[z]*temperature[z] / (3.0 * opacity[z] * rho[z]);


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

    prop = 0.000;

    G = 6.67428e-8;

    mp = 1.0; // Planetary mass in terms of Jupiter masses
    mp = 1.8986e30 * mp; // converted into g
    
    Mstar = 1.0; // Stellar mass in solar masses
    Mstar = Mstar * 1.9892e33; // Converted into g
    

    
    
    D = 0.0512; // Orbital radius of planet in AU, 1.58740105197 = 4^(1/3), which acts to double the period of the orbit
        
    D = 1.495978921e13 * D; // Converted into cm

    f = -(G * mp) / (4.0 * D * D * D);



    omega = sqrt(G * ((Mstar + mp)/(D*D*D)));

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



    // count[kold] = tyhe number of cells used to get to the outer edge of cell kold (including that edge).  Therefore J = count[Jold-1], and the extra number of cells = J - Jold
    int N[J], count[J], Jold, k_start, k_end;
    
    int U, E, L, Nmax, n;
    
    double dx_max;
    
    Jold = J;
    
    count[0] = 1;
    
    
    
    
    // This sets N[0] and N[Jold-1] =0, as the first and last cells will be left alone
    N[0] = 0;
    N[Jold-1] = 0;
    
    
    // This sums over the extra cells needed to achieve at least the minimum resolution, given by the maximum cell size: dx_max
    // It only goes up to k < Jold - 2 because we don't want to split the final cell, that is the cell from radius_cm[Jold-2] to radius_cm[Jold-1]
    
    for (k=1; k < Jold-1; k = k + 1) {
        
        // by defining dx_max inside this loop, we can make it a function of radius_cm, for instance

        dx_max = 0.000015 + 0.0004*(radius_cm[k]/R)*(radius_cm[k]/R);
        
        N[k] = ((radius_cm[k] - radius_cm[k-1])/(R*dx_max));
        
        count[k] = count[k-1] + 1 + N[k];
        
        //cout << "At k = " << k << ", N = " << N << "\n";
      
        
    }
    
    count[Jold-1] = count[Jold-2] + 1;
    
    
    cout << "\n\ncount[Jold-1] = " << count[Jold-1] << "\n\n";

    // This gives the new value of J, by taking into account all of the extra cells needed
    J = count[Jold-1];
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    






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

    double count_doub,N_doub,k_doub;

    // These are defined to keep only two variables at any one time - k and k+1, which will be re-written for each new zone
    double lnT_HR[2], lnRho_HR[2], grav_HR[2], radius_cm_HR[2], rmid_cm_HR[2];
    double temperature_HR[2], rho_HR[2], pressure_HR[2], grada_HR[2], cp_HR[2], chiRho_HR[2], chiT_HR[2];
    double opacity_HR[2], dkap_dlnrho_face_HR[2], dkap_dlnT_face_HR[2], flux_HR[2], brunt_A_HR[2];
    double K_HR[2], rho_face_HR[2];

    // This is defined for the k and k+1 values of delta which help with linear interpolation: F_{at rmid[1]} = delta[0]*flux[0] + delta[1]*flux[1]
    // The third one is for interpolating cell-central variables, because you need to involve k+2 variables for that
    double delta[3];

    // This is an exception, as this needs to be output at the end for plotting purposes
    double radius_cm_HR_output[J],rmid_cm_HR_output[J],flux_HR_output[J],pressure_HR_output[J],temperature_HR_output[J],test[J];
    
    
    
    
    
    
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
    
    
    
    
    
    // This is the start of the big for loop in which all of the matrix stuff is setup and first used
    
    for (k=0; k < J-1; k=k+1) {
        
        kold = 0;
        
        if (k == 0) {
            kold = 0;
        } else {
            
            sum = 0;

            
            // This finds what kold is for this loop, by going to the point where radius_cm[kold] has first exceeded radius_cm_HR_output[k-1]
            for (kold = 0; k > sum; kold = kold + 1) {

                // It's count[kold+1] because of the fact that kold is incremented at the end of the loop, and the new value of kold is actually the one we want to test against sum
                sum = count[kold+1] - 1;
                
            }
            
        }
        

        /*
        // This sets the number of extra cells being added to each original cell, with both the innermost and outermost cells being left alone
        if (kold == 0 || kold == Jold - 1 ) {
            N = 0;
        } else {
            N = ( radius_cm[kold] - radius_cm[kold-1] ) / ( R * dx_max );
        }
        */
        
        N_doub = (double) (N[kold]);
        // Note that count_doub is NOT for kold, but for kold-1
        count_doub = (double) (count[kold-1]);
        k_doub = (double) (k);
     
        
        
        
        
        
            // These give the variables at the surface of the cell as: q[k] = delta[0]*q0[kold-1] + delta[1]*q0[kold]
            if (N[kold] == 0) {
                delta[0] = 0.0;
                delta[1] = 1.0;
                
            } else {
                
                // We now need to calculate where we are in the cell, using the last interpolation point and the known value of N
                
                // We use the fact that the fractional distance through the cell of the last known point is: ( ( radius_cm_HR_output[k-1] - radius_cm[kold - 1] )/( radius_cm[kold] radius_cm[kold-1] ) )
                // and the fractional step size within the cell is 1.0/(N_doub+1.0)
                // so we just increase the distance through the cell by that step
                
                
                delta[0] = 1.0 - (k_doub - count_doub + 1.0)/(N_doub + 1.0);
                delta[1] = (k_doub - count_doub + 1.0)/(N_doub + 1.0);
            }
            
            // This prevents issues with calling the -1st element of an array when kold == 0
            if (kold == 0) {
                
                radius_cm_HR_output[k] = delta[1]*radius_cm[kold];
                
                radius_cm_HR[0] = delta[1]*radius_cm[kold];
                flux_HR[0] = delta[1]*flux[kold];
                dkap_dlnrho_face_HR[0] = delta[1]*dkap_dlnrho_face[kold];
                dkap_dlnT_face_HR[0] = delta[1]*dkap_dlnT_face[kold];
                opacity_HR[0] = delta[1]*opacity[kold];
                rho_face_HR[0] = delta[1]*rho_face[kold];

                
            } else {
                
                radius_cm_HR_output[k] = delta[0]*radius_cm[kold-1] + delta[1]*radius_cm[kold];
                
                radius_cm_HR[0] = delta[0]*radius_cm[kold-1] + delta[1]*radius_cm[kold];
                flux_HR[0] = delta[0]*flux[kold-1] + delta[1]*flux[kold];
                dkap_dlnrho_face_HR[0] = delta[0]*dkap_dlnrho_face[kold-1] + delta[1]*dkap_dlnrho_face[kold];
                dkap_dlnT_face_HR[0] = delta[0]*dkap_dlnT_face[kold-1] + delta[1]*dkap_dlnT_face[kold];
                opacity_HR[0] = delta[0]*opacity[kold-1] + delta[1]*opacity[kold];
                rho_face_HR[0] = delta[0]*rho_face[kold-1] + delta[1]*rho_face[kold];

            }
            
            
            
            // This gets rmid_cm_HR_output, which will be used to determine exactly how to interpolate for the variables at the cell centres
            if (k==0) {
                rmid_cm_HR_output[k] = rmid_cm[0];
            } else {
                rmid_cm_HR_output[k] = 0.5*(radius_cm_HR_output[k] + radius_cm_HR_output[k-1]);
            }
            
            
            
            
            // These give the variables at the centre of the cell as: p[k] = delta[0]*p0[kold-1] + delta[1]*p0[kold] + delta[2]*p0[kold+1]
            if (N[kold] == 0) {
                
                // Here the cell-central variables are defined for the case where no further cell-division is needed (including the innermost and outermost cells as exceptions to the rule (set as N = 0))
                
                rmid_cm_HR[0] = rmid_cm[kold];
                lnRho_HR[0] = lnRho[kold];
                rho_HR[0] = rho[kold];
                cp_HR[0] = cp[kold];
                temperature_HR[0] = temperature[kold];
                lnT_HR[0] = lnT[kold];
                pressure_HR[0] = pressure[kold];
                grav_HR[0] = grav[kold];
                K_HR[0] = K[kold];
                chiRho_HR[0] = chiRho[kold];
                chiT_HR[0] = chiT[kold];
                grada_HR[0] = grada[kold];
                brunt_A_HR[0] = brunt_A[kold];
                
                
                
            } else {
                if (rmid_cm_HR_output[k] <= rmid_cm[kold]) {
                    
                    delta[0] = 1.0 - ( ( rmid_cm_HR_output[k] - rmid_cm[kold-1] )/( rmid_cm[kold] - rmid_cm[kold-1] ) );
                    delta[1] = ( ( rmid_cm_HR_output[k] - rmid_cm[kold-1] )/( rmid_cm[kold] - rmid_cm[kold-1] ) );
                    delta[2] = 0.0;
                    
                } else {
                    
                    delta[0] = 0.0;
                    delta[1] = 1.0 - ( ( rmid_cm_HR_output[k] - rmid_cm[kold] )/( rmid_cm[kold+1] - rmid_cm[kold] ) );
                    delta[2] = ( ( rmid_cm_HR_output[k] - rmid_cm[kold] )/( rmid_cm[kold+1] - rmid_cm[kold] ) );
                    
                }
                
                // Here the cell-central variables are actually defined
                
                rmid_cm_HR[0] = delta[0]*rmid_cm[kold-1] + delta[1]*rmid_cm[kold] + delta[2]*rmid_cm[kold+1];
                lnRho_HR[0] = delta[0]*lnRho[kold-1] + delta[1]*lnRho[kold] + delta[2]*lnRho[kold+1];
                rho_HR[0] = delta[0]*rho[kold-1] + delta[1]*rho[kold] + delta[2]*rho[kold+1];
                cp_HR[0] = delta[0]*cp[kold-1] + delta[1]*cp[kold] + delta[2]*cp[kold+1];
                temperature_HR[0] = delta[0]*temperature[kold-1] + delta[1]*temperature[kold] + delta[2]*temperature[kold+1];
                lnT_HR[0] = delta[0]*lnT[kold-1] + delta[1]*lnT[kold] + delta[2]*lnT[kold+1];
                pressure_HR[0] = delta[0]*pressure[kold-1] + delta[1]*pressure[kold] + delta[2]*pressure[kold+1];
                grav_HR[0] = delta[0]*grav[kold-1] + delta[1]*grav[kold] + delta[2]*grav[kold+1];
                K_HR[0] = delta[0]*K[kold-1] + delta[1]*K[kold] + delta[2]*K[kold+1];
                chiRho_HR[0] = delta[0]*chiRho[kold-1] + delta[1]*chiRho[kold] + delta[2]*chiRho[kold+1];
                chiT_HR[0] = delta[0]*chiT[kold-1] + delta[1]*chiT[kold] + delta[2]*chiT[kold+1];
                grada_HR[0] = delta[0]*grada[kold-1] + delta[1]*grada[kold] + delta[2]*grada[kold+1];
                brunt_A_HR[0] = delta[0]*brunt_A[kold-1] + delta[1]*brunt_A[kold] + delta[2]*brunt_A[kold+1];
                
  
            
        }
        
        
        
        
        //cout << "radius_cm_HR_output[k]/R = " << radius_cm_HR_output[k]/R << "     rmid_cm_HR_output[k]/R = " << rmid_cm_HR_output[k]/R << "     kold = " << kold << "    k = " << k << "\n";
        
        
        
        
        
        
        //
        //
        //
        // Here the variables are defined for k+1, for which we will use i=k+1
        //
        //
        //
        
        i=k+1;
        
        
        
        if (i == 0) {
            kold = 0;
        } else {
            if ( i == J-1) {
                cout << "Flag i == J-1\n";
                kold = Jold - 1;
            } else {
                
                sum = 0;
                
                // This finds what kold is for this loop
                // sum = the total number of cells required to describe everything up to kold
                for (kold = 0; i > sum; kold = kold + 1) {
                    
                    // It's count[kold+1] because of the fact that kold is incremented at the end of the loop, and the new value of kold is actually the one we want to test against sum
                    sum = count[kold+1] - 1;
                    
                }
                
            }
            
        }
        
        
        /*
        // This sets the number of extra cells being added to each original cell, with both the innermost and outermost cells being left alone
        if (kold == 0 || kold == Jold - 1 ) {
            N = 0;
        } else {
            N = ( radius_cm[kold] - radius_cm[kold-1] ) / ( R * dx_max );
        }
        */
        
        
        N_doub = (double) (N[kold]);
        // Note that count_doub is NOT for kold, but for kold-1
        count_doub = (double) (count[kold-1]);
        k_doub = (double) (i);
        
        
        
        
        
        
        // These give the variables at the surface of the cell as: q[i] = delta[0]*q0[kold] + delta[1]*q0[kold+1]
        if (N[kold] == 0) {
            delta[0] = 0.0;
            delta[1] = 1.0;
            
        } else {
            
            // We now need to calculate where we are in the cell, using the last interpolation point and the known value of N
            
            // We use the fact that the fractional distance through the cell of the last known point is: ( ( radius_cm_HR_output[k-1] - radius_cm[kold - 1] )/( radius_cm[kold] radius_cm[kold-1] ) )
            // and the fractional step size within the cell is 1.0/(N_doub+1.0)
            // so we just increase the distance through the cell by that step
            
            
            delta[0] = 1.0 - (k_doub - count_doub + 1.0)/(N_doub + 1.0);
            delta[1] = (k_doub - count_doub + 1.0)/(N_doub + 1.0);
        }
        
        // This prevents issues with calling the -1st element of an array when kold == 0
        if (kold == 0 || kold == J-1) {
            
            radius_cm_HR_output[i] = delta[1]*radius_cm[kold];
            
            radius_cm_HR[1] = delta[1]*radius_cm[kold];
            flux_HR[1] = delta[1]*flux[kold];
            dkap_dlnrho_face_HR[1] = delta[1]*dkap_dlnrho_face[kold];
            dkap_dlnT_face_HR[1] = delta[1]*dkap_dlnT_face[kold];
            opacity_HR[1] = delta[1]*opacity[kold];
            rho_face_HR[1] = delta[1]*rho_face[kold];
            
            
        } else {
            
            radius_cm_HR_output[i] = delta[0]*radius_cm[kold-1] + delta[1]*radius_cm[kold];
            
            radius_cm_HR[1] = delta[0]*radius_cm[kold-1] + delta[1]*radius_cm[kold];
            flux_HR[1] = delta[0]*flux[kold-1] + delta[1]*flux[kold];
            dkap_dlnrho_face_HR[1] = delta[0]*dkap_dlnrho_face[kold-1] + delta[1]*dkap_dlnrho_face[kold];
            dkap_dlnT_face_HR[1] = delta[0]*dkap_dlnT_face[kold-1] + delta[1]*dkap_dlnT_face[kold];
            opacity_HR[1] = delta[0]*opacity[kold-1] + delta[1]*opacity[kold];
            rho_face_HR[1] = delta[0]*rho_face[kold-1] + delta[1]*rho_face[kold];
            
        }
        
        
        
        // This gets rmid_cm_HR_output, which will be used to determine exactly how to interpolate for the variables at the cell centres
        if (i==0) {
            rmid_cm_HR_output[i] = rmid_cm[0];
        } else {
            rmid_cm_HR_output[i] = 0.5*(radius_cm_HR_output[i] + radius_cm_HR_output[i-1]);
        }
        
        
        
        
        
        // These give the variables at the centre of the cell as: p[i] = delta[0]*p0[kold-1] + delta[1]*p0[kold] + delta[2]*p0[kold+1]
        if (N[kold] == 0 || kold == J-1) {
            
            // Here the cell-central variables are defined for the case where no further cell-division is needed (including the innermost and outermost cells as exceptions to the rule (set as N = 0))
            
            //cout << "Flag N == 0 case\n";
            
            //cout << "kold = " << kold << "   rmid_cm[kold]/R = " << rmid_cm[kold]/R << "\n";
            
            rmid_cm_HR[1] = rmid_cm[kold];
            lnRho_HR[1] = lnRho[kold];
            rho_HR[1] = rho[kold];
            cp_HR[1] = cp[kold];
            temperature_HR[1] = temperature[kold];
            lnT_HR[1] = lnT[kold];
            pressure_HR[1] = pressure[kold];
            grav_HR[1] = grav[kold];
            K_HR[1] = K[kold];
            chiRho_HR[1] = chiRho[kold];
            chiT_HR[1] = chiT[kold];
            grada_HR[1] = grada[kold];
            brunt_A_HR[1] = brunt_A[kold];
            
            
            
        } else {
            
            //cout << "Flag N =/= 0 case\n";
            
            //cout << "Flag with information: Jold=" << Jold << "   J=" << J << "   k=" << k << "   i=" << i << "   radius_cm_HR_output[k]/R=" << radius_cm_HR_output[k]/R << "   kold=" << kold << "\n";
            
            if (rmid_cm_HR_output[i] <= rmid_cm[kold]) {
                
                //cout << "Flag if case, rmid_cm[kold]/R = " << rmid_cm[kold]/R << "\n";
                
                delta[0] = 1.0 - ( ( rmid_cm_HR_output[i] - rmid_cm[kold-1] )/( rmid_cm[kold] - rmid_cm[kold-1] ) );
                delta[1] = ( ( rmid_cm_HR_output[i] - rmid_cm[kold-1] )/( rmid_cm[kold] - rmid_cm[kold-1] ) );
                delta[2] = 0.0;
                
            } else {
                
                //cout << "Flag else case, rmid_cm[kold]/R = " << rmid_cm[kold]/R << "\n";
                
                delta[0] = 0.0;
                delta[1] = 1.0 - ( ( rmid_cm_HR_output[i] - rmid_cm[kold] )/( rmid_cm[kold+1] - rmid_cm[kold] ) );
                delta[2] = ( ( rmid_cm_HR_output[i] - rmid_cm[kold] )/( rmid_cm[kold+1] - rmid_cm[kold] ) );
                
            }
            
            // Here the cell-central variables are actually defined
            
            //cout << "Flag pre-definition\n";
            
            //cout << "kold = " << kold << "   rmid_cm[kold]/R = " << rmid_cm[kold]/R << "\n";
            
            rmid_cm_HR[1] = delta[0]*rmid_cm[kold-1] + delta[1]*rmid_cm[kold] + delta[2]*rmid_cm[kold+1];
            lnRho_HR[1] = delta[0]*lnRho[kold-1] + delta[1]*lnRho[kold] + delta[2]*lnRho[kold+1];
            rho_HR[1] = delta[0]*rho[kold-1] + delta[1]*rho[kold] + delta[2]*rho[kold+1];
            cp_HR[1] = delta[0]*cp[kold-1] + delta[1]*cp[kold] + delta[2]*cp[kold+1];
            temperature_HR[1] = delta[0]*temperature[kold-1] + delta[1]*temperature[kold] + delta[2]*temperature[kold+1];
            lnT_HR[1] = delta[0]*lnT[kold-1] + delta[1]*lnT[kold] + delta[2]*lnT[kold+1];
            pressure_HR[1] = delta[0]*pressure[kold-1] + delta[1]*pressure[kold] + delta[2]*pressure[kold+1];
            grav_HR[1] = delta[0]*grav[kold-1] + delta[1]*grav[kold] + delta[2]*grav[kold+1];
            K_HR[1] = delta[0]*K[kold-1] + delta[1]*K[kold] + delta[2]*K[kold+1];
            chiRho_HR[1] = delta[0]*chiRho[kold-1] + delta[1]*chiRho[kold] + delta[2]*chiRho[kold+1];
            chiT_HR[1] = delta[0]*chiT[kold-1] + delta[1]*chiT[kold] + delta[2]*chiT[kold+1];
            grada_HR[1] = delta[0]*grada[kold-1] + delta[1]*grada[kold] + delta[2]*grada[kold+1];
            brunt_A_HR[1] = delta[0]*brunt_A[kold-1] + delta[1]*brunt_A[kold] + delta[2]*brunt_A[kold+1];
            
            //cout << "Flag post-definition\n";
            
            
            
        }
        
        
        
        
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





        



        // This checks the interpolation of any given variable array
        test[k] = (radius_cm_HR[1] - radius_cm_HR[0])/R;
        test[k+1] = test[k];

        
        
        
        flux_HR_output[k] = flux_HR[0];
        flux_HR_output[k+1] = flux_HR[1];
        
        pressure_HR_output[k] = pressure_HR[0];
        pressure_HR_output[k+1] = pressure_HR[1];
        
        temperature_HR_output[k] = temperature_HR[0];
        temperature_HR_output[k+1] = temperature_HR[1];

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
        Ai[0][1][0] = ( brunt_A_HR[1]/( R ) ) * ( chiRho_HR[1] / chiT_HR[1] ) * 0.5 * radius_cm_HR[0]; // ( brunt_A_HR[1]/( R*R ) ) * ( chiRho_HR[1] / chiT_HR[1] ) * 0.5 * radius_cm_HR[0] * radius_cm_HR[0];
        Ai[0][1][1] = 0.0;
        
        Cr[0][0][0] = radius_cm_HR[1]*radius_cm_HR[1]*rho_face_HR[1] / ( rho_HR[1]*R*( radius_cm_HR[1] - radius_cm_HR[0] ) ); // radius_cm_HR[1]*radius_cm_HR[1]*radius_cm_HR[1]*rho_face_HR[1] / ( rho_HR[1]*R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) );
        Cr[0][0][1] = 0.0;
        Cr[0][1][0] = 0.0;
        Cr[0][1][1] = ( radius_cm_HR[1]*radius_cm_HR[1]*flux_BC )/( m*omega*rho_HR[1]*temperature_HR[1]*cp_HR[1]*R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) ); // ( radius_cm_HR[1]*radius_cm_HR[1]*flux_HR[1] )/( m*omega*rho_HR[1]*temperature_HR[1]*cp_HR[1]*R*R*( radius_cm_HR[1] - radius_cm_HR[0] ) );
        
        Ci[0][0][0] = 0.0;
        Ci[0][0][1] = 0.0;
        Ci[0][1][0] = ( brunt_A_HR[1]/( R ) ) * ( chiRho_HR[1] / chiT_HR[1] ) * 0.5 * radius_cm_HR[1]; // ( brunt_A_HR[1]/( R*R ) ) * ( chiRho_HR[1] / chiT_HR[1] ) * 0.5 * radius_cm_HR[1] * radius_cm_HR[1];
        Ci[0][1][1] = 0.0;
        
        Dr[0][0][0] = ( rmid_cm_HR[1]*rmid_cm_HR[1] / ( chiRho_HR[1]*R*R ) )  -  ( (l*(l+1.0)*pressure_HR[1])/(m*m*omega*omega*R*R*rho_HR[1]) );
        Dr[0][0][1] = - (chiT_HR[1] / chiRho_HR[1]) * ( rmid_cm_HR[1]*rmid_cm_HR[1] / (R*R) );
        Dr[0][1][0] = 0.0;
        Dr[0][1][1] = ( l * ( l + 1.0 ) * K_HR[1] )/( m * omega * rho_HR[1] * cp_HR[1] * R * R );
        
        Di[0][0][0] = 0.0;
        Di[0][0][1] = 0.0;
        Di[0][1][0] = - grada_HR[1] * rmid_cm_HR[1] * rmid_cm_HR[1] / ( R*R );
        Di[0][1][1] = rmid_cm_HR[1] * rmid_cm_HR[1] / ( R * R );
        
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
        // This sets the equations to be adiabatic for x > 0.7 (that is, once the convective envelope starts)
        if ( 0.75 < radius_cm_HR[0]/R  && radius_cm_HR[0]/R < 1.1) {
            
            prop = 0.0;
            
            if ( 0.75 < radius_cm_HR[0]/R  && radius_cm_HR[0]/R < 0.8 ) {
                
                prop = 1.0 - ( ( (radius_cm_HR[0]/R) - 0.7) / 0.05);
                
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
        
        
        
        
        
        Mr[0][0][0] = ( ( ( l*(l+1.0)*rmid_cm_HR[1]*rmid_cm_HR[1]*f )/( m*m*omega*omega*R*R ) ) )*( 1.0 + prop*sin(radius_cm_HR[0]*rmid_cm_HR[1]*rho_HR[0]*pressure_HR[1])); // Ar[0][0][0]*a[0] + Ar[0][0][1]*b[0] + Cr[0][0][0]*a[1] + Cr[0][0][1]*b[1] + Dr[0][0][0]*c[1] + Dr[0][0][1]*d[1];  // ( ( l*(l+1.0)*rmid_cm_HR[1]*rmid_cm_HR[1]*f )/( m*m*omega*omega*R*R ) );
        Mr[0][1][0] = 0.0; // Ar[0][1][0]*a[0] + Ar[0][1][1]*b[0] + Cr[0][1][0]*a[1] + Cr[0][1][1]*b[1] + Dr[0][1][0]*c[1] + Dr[0][1][1]*d[1];  // 0.0;
        
        Mi[0][0][0] = 0.0; // Ai[0][0][0]*a[0] + Ai[0][0][1]*b[0] + Ci[0][0][0]*a[1] + Ci[0][0][1]*b[1] + Di[0][0][0]*c[1] + Di[0][0][1]*d[1];  // 0.0;
        Mi[0][1][0] = 0.0; // Ai[0][1][0]*a[0] + Ai[0][1][1]*b[0] + Ci[0][1][0]*a[1] + Ci[0][1][1]*b[1] + Di[0][1][0]*c[1] + Di[0][1][1]*d[1];  // 0.0;
        
        
        Nr[0][0][0] = 0.0; // Er[0][0][0]*a[0] + Er[0][0][1]*b[0] + Fr[0][0][0]*c[0] + Fr[0][0][1]*d[0] + Hr[0][0][0]*c[1] + Hr[0][0][1]*d[1];  // 0.0;
        Nr[0][1][0] = ( - 2.0 * f * radius_cm_HR[0] / ( m*m*omega*omega*R ) ) * ( 1.0 + prop*sin(temperature_HR[1]*rho_HR[0]*radius_cm_HR[0]*temperature_HR[0]*rho_HR[1])); // Er[0][1][0]*a[0] + Er[0][1][1]*b[0] + Fr[0][1][0]*c[0] + Fr[0][1][1]*d[0] + Hr[0][1][0]*c[1] + Hr[0][1][1]*d[1];  // - 2.0 * f * radius_cm_HR[0] / ( m*m*omega*omega*R );
        
        Ni[0][0][0] = 0.0; // Ei[0][0][0]*a[0] + Ei[0][0][1]*b[0] + Fi[0][0][0]*c[0] + Fi[0][0][1]*d[0] + Hi[0][0][0]*c[1] + Hi[0][0][1]*d[1];  // 0.0;
        Ni[0][1][0] = 0.0; // Ei[0][1][0]*a[0] + Ei[0][1][1]*b[0] + Fi[0][1][0]*c[0] + Fi[0][1][1]*d[0] + Hi[0][1][0]*c[1] + Hi[0][1][1]*d[1];  // 0.0;
        
        
        
        
        
        
        
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
    outfile.open("Output/Henyey_dx_max_perturb.dat", ios::out);

    // This sets the precision at which values are printed at to the named file output
    outfile.precision(10);


    for (k=0; k<J; k=k+1) {

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
         */
        
        
        
        outfile << radius_cm_HR_output[k]/R << "\t\t\t" << ur[k][0][0] << "\t\t" << ur[k][1][0] << "\t\t" << vr[k][0][0] << "\t\t" << vr[k][1][0] << "\t\t" << alphar[k][0][0] << "\t" << alphar[k][0][1] << "\t" << alphar[k][1][0] << "\t" << alphar[k][1][1] << "\t" << alphai[k][0][0] << "\t" << alphai[k][0][1] << "\t" << alphai[k][1][0] << "\t" << alphai[k][1][1] << "\t" << gammar[k][0][0] << "\t" << gammar[k][1][0] << "\t" << gammai[k][0][0] << "\t" << gammai[k][1][0] <<  "\t" << ur[k][0][0]*radius_cm_HR_output[k]*(1.0/100.0)*m*omega*((mp + Mstar)/(mp)) << "\t" << test[k] << "\t" << radius_cm_HR_output[k] << "\t" << flux_HR_output[k] << "\t" << pressure_HR_output[k] << "\t" << temperature_HR_output[k] << "\t" << rmid_cm_HR_output[k]/R << "\t\t\t" << R*ur[k][0][0] << "\t\t" << flux_BC*ur[k][1][0] << "\t\t" << pressure_HR_output[k]*vr[k][0][0] << "\t\t" << temperature_HR_output[k]*vr[k][1][0] << "\t\t\t" << R*ur[k][0][0]*m*omega*(mp + Mstar)/(100.0*mp) << "\t\t\t" << ui[k][0][0] << "\t\t" << ui[k][1][0] << "\t\t" << vi[k][0][0] << "\t\t" << vi[k][1][0] << "\t\t" << sqrt((ur[k][0][0]*ur[k][0][0]) + (ui[k][0][0]*ui[k][0][0])) << "\t\t" << atan(-ui[k][0][0] / ur[k][0][0]) << "\n";


    }

    // This closes the output data file
    outfile.close();


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
    
    cout << "\n\n   total number of new cells needed = " << sum << "\n\n";
    
    cout << "\n\n   total number of cells = " << J << "\n\n";
    
    sum = 3.2 / 1.01;
    
    cout << "\n\n 3.2 / 1.01 = " << sum << "\n";
    







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





