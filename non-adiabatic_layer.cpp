#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>


using namespace std;




int main()
{
    
    int k, z, J;
    double C;
    
    

    
    ifstream infile;
    infile.open("Input/non-adiabatic_layer_1_0.txt");
    
    k=0;
    
    double input;
    int no_of_lines;
    string line;
    
    no_of_lines = 0;
    
    ifstream linefile;
    linefile.open("Input/non-adiabatic_layer_1_0.txt");
    
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
    
    double radius_cm[J], rmid_cm[J];
    double temperature[J], rho[J];
    double flux[J], dr[J], luminosity[J], mu[J];
    double cv[J], gamma1[J], gamma3[J], conv_L_div_L[J];
    
    
    
    cout << "FLAG - about to take the input\n";
    
    
    // Get the first line's input
    infile >> input;
    
    // Then process the first line's input, then try to get another line.
    // This works because infile==true if the last thing it tried to read successfully gave it something.
    // Therefore, the while loop will only end when the condition is tested after the read has failed.
    // Therefore we need to keep the reading of line m and the processing of line m separated by a test of the condition
    while (infile) {
        
        radius_cm[z] = input;
        
        infile >> input;
        
        rmid_cm[z] = input;
        
        infile >> input;
        
        temperature[z] = input;
        
        infile >> input;
        
        rho[z] = input;
        
        infile >> input;
        
        flux[z] = input;
        
        infile >> input;
        
        cv[z] = input;
        
        infile >> input;
        
        gamma1[z] = input;
        
        infile >> input;
        
        gamma3[z] = input;
        
        infile >> input;
        
        conv_L_div_L[z] = input;
        
        infile >> input;
        
        dr[z] = input;
        
        infile >> input;
        
        mu[z] = input;
        
        infile >> input;
        
        luminosity[z] = input;

        
        infile >> input; // This is to get the input for the zone for the next loop (and therefore to decide whether or not to do the next loop, too)
        
        

        
        
        z = z - 1;
    }
    
    
    
    double omega, pi, m_H, N_A;
    
    pi = 3.14159265358979;
    
    omega = 2 * pi / ( 4.23 * 24 * 3600 );
    
    m_H = 1.6605402e-24;
    
    N_A = 6.0221367e23;
    
    
    
    
    
    double dE, E_internal, dE_scaled, E_internal_scaled;
    
    dE = 0.0;
    
    E_internal = 0.0;
    
    dE_scaled = 0.0;
    
    E_internal_scaled = 0.0;
    
    
    ofstream write_file;
    
    write_file.open("Output/non-adiabatic_layer_1_0.txt", ios::out);
    
    write_file.precision(10);
    
    
    
    for (k=J-1; k>=0; k=k-1) {
        
        dE = ( 4 * pi / (m_H * N_A) ) * (  ( cv[k] * rho[k] * temperature[k] * rmid_cm[k] * rmid_cm[k] ) / ( mu[k] )  ) * dr[k];
        
        E_internal = E_internal + dE;
        
        
        dE_scaled = dE * 2.0 * omega / luminosity[k];
        
        E_internal_scaled = E_internal_scaled + dE_scaled;
        
        
        // 1 - radius_cm
        // 2 - E_internal
        // 3 - dE
        // 4 - E_internal_scaled (that is, E_internal * omega / L )
        // 5 - dE_scaled (that is, dE * omega / L )
        
        write_file << radius_cm[k] << "\t\t" << E_internal << "\t\t" << dE << "\t\t" << E_internal_scaled << "\t\t" << dE_scaled << "\n";
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}
