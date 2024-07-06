/* 
This is a script to test out generation of generalized ChIMES model that does not require explicit parmeters for different element types
R. K. Lindsey (2024)

This only works for 2-body models.

This script evaluates energy for a given configuration, based on an un-formatted parameter file.

This is a horribly written script and should really be integrated into the test_chemfit utility...


// Run with: 
// g++ <this file>
// ./a.out Oe Od params.txt elemtypeidx_1 elemtypeidx_2

*/


#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<vector>
#include<cmath>
#include<cstdlib>
#include<algorithm>
#include<sstream>

using namespace std;

// Helper utilities/objects


int split_line(string line, vector<string> & items)
{
    // Break a line up into tokens based on space separators.
    // Returns the number of tokens parsed.
    
    string       contents;
    stringstream sstream;

    // Strip comments beginining with ! or ## and terminal new line

    int pos = line.find('!');
      
    if ( pos != string::npos ) 
        line.erase(pos, line.length() - pos);

    pos = line.find("##");
    if ( pos != string::npos ) 
        line.erase(pos, line.length()-pos);

    pos = line.find('\n');
    if ( pos != string::npos ) 
        line.erase(pos, 1);

    sstream.str(line);
     
    items.clear();

    while ( sstream >> contents ) 
        items.push_back(contents);

    return items.size();
}


bool get_next_line(istream& str, string & line)
{
    // Read a line and return it, with error checking.
    
        getline(str, line);

        if(!str)
            return false;
    
    return true;
}

struct xyz
{
    double x;
    double y;
    double z;
};


// Hyperparameters class and member functions

class hyperparams
{
    public:
        
        vector<string> element;
        vector<double> vdw;
        vector<double> eion;
        
        int O2B_element;
        int O2B_distance;
        int max_order;
        int nparams;
        
        
        hyperparams(int Oe, int Od);
        ~hyperparams();
        
        int assign_elemtype(string element); 
        void universalize_dist(double rvdwi, double rvdwj, double & morse_lambda, double & rcin, double & rout);
        void universalize_elem();
        void transform_dist(double rij, double morse_lambda, double rcin, double rout, double & sij, double & dsij_drij);
        double transform_eion(double eioni, double eionj);
                       
};

hyperparams::hyperparams(int Oe, int Od)
{
    element.resize(7);
    vdw    .resize(7);
    eion   .resize(7);
    
    
    // vdw rad is in Ang, eion is in eV, converted to kcal/mol
    element[0] = "Si "; vdw[0] = 2.0*2.10; eion[0] = 8.1520 * 23.0609; // Other properties: unpaired valence = 4; electoneg = 1.90 (Pauling Scale); eletronaffin = 1.385 (eV)
    element[1] = "C"  ; vdw[1] = 2.0*1.70; eion[1] = 11.260 * 23.0609; // Other properties: unpaired valence = 4; electoneg = 2.55 (Pauling Scale); eletronaffin = 1.263 (eV)
    element[2] = "N"  ; vdw[2] = 2.0*1.55; eion[2] = 14.534 * 23.0609; // Other properties: unpaired valence = 3; electoneg = 3.04 (Pauling Scale); eletronaffin = 0.000 (eV)
    element[3] = "O"  ; vdw[3] = 2.0*1.52; eion[3] = 13.618 * 23.0609; // Other properties: unpaired valence = 2; electoneg = 3.44 (Pauling Scale); eletronaffin = 1.461 (eV)
    element[4] = "H"  ; vdw[4] = 2.0*1.20; eion[4] = 13.598 * 23.0609; // Other properties: unpaired valence = 1; electoneg = 2.20 (Pauling Scale); eletronaffin = 0.754 (eV)
    element[5] = "He" ; vdw[5] = 2.0*1.40; eion[5] = 24.587 * 23.0609; // Other properties: unpaired valence = 0; electoneg = 0.00 (Pauling Scale); eletronaffin = 0.000 (eV)
    element[6] = "Ar" ; vdw[6] = 2.0*1.88; eion[6] = 15.760 * 23.0609; // Other properties: unpaired valence = 0; electoneg = 0.00 (Pauling Scale); eletronaffin = 0.000 (eV)
     
    O2B_element  = Oe;
    O2B_distance = Od;
    
    nparams = Oe * Od;
        
    max_order = O2B_element;
    if (O2B_distance > max_order)
        max_order = O2B_distance;    
}

hyperparams::~hyperparams(){}

int hyperparams::assign_elemtype(string element_in)
{
    auto it = find(element.begin(), element.end(), element_in); 
    
    if (it != element.end()) // If element was found 
        return it - element.begin();   // calculating the index of element_in 
    else 
    { 
         cout << "Could not find a hyperparameter definition for element type " << element_in << endl;
         exit(0);
    } 
}


void hyperparams::universalize_dist(double vdwi, double vdwj, double & morse_lambda, double & rcin, double & rout)
{    
    double sigma_i      = vdwi;
    double sigma_j      = vdwj;
    double sigma_ij     = 0.5 * (sigma_i + sigma_j);
    
    morse_lambda = pow(2.0,1.0/6.0)*sigma_ij;
    rout  = 3.5  * sigma_ij;
    rcin  = 0.65*morse_lambda; // 0.45 * morse_lambda;
} 


void hyperparams::transform_dist(double rij, double morse_lambda, double rcin, double rout, double & sij, double & dsij_drij)
{
    double x     = exp(-1 * rij  / morse_lambda);   
    double xin   = exp(-1 * rout / morse_lambda);  
    double xout  = exp(-1 * rcin / morse_lambda);  
    double xavg  = 0.5 * (xout + xin);
    double xdiff = 0.5 * (xout - xin);
    
    sij       = (x - xavg)/xdiff;                   
    dsij_drij = -1 * x / morse_lambda / xdiff;
}

double hyperparams::transform_eion(double eioni, double eionj)
{
    double eionij    = sqrt(eioni*eionj);
    double eion_min  = 3.894 * 23.0609;    // Ionization energy for Cs, converted from eV to kcal/mol
    double eion_max  = 24.587 * 23.0609;   // Ionization energy for He, converted from eV to kcal/mol
    double eion_del  = eion_max - eion_min;
    double eion_avg  = 0.5*(eion_max + eion_min);
    double eion_diff = 0.5*(eion_max - eion_min);

    //return pow( (eionij - eion_min) / eion_del, 3.0);
    //return pow( (eionij - eion_avg) / eion_diff, 3.0);
    return (eionij - eion_avg) / eion_diff;
}



class params : public hyperparams
{
    public:
    
        vector<int> O2BD;
        vector<int> O2BE;
        vector<double> coeff;
        
        params(int Oe, int Od) : hyperparams(Oe, Od){};
        ~params(){};
        void set_hyperparams(int O2B_element_in, int O2B_distance_in, int max_order_in, int nparams_in);
        void read_params(string paramfile);
};

void params::set_hyperparams(int O2B_element_in, int O2B_distance_in, int max_order_in, int nparams_in)
{
    O2B_element  = O2B_element_in;
    O2B_distance = O2B_distance_in;
    max_order    = max_order_in;
    nparams      = nparams_in;
    
}
    
void params::read_params(string paramfile)
{
    cout << "# Reading parameter file: " << paramfile << endl;
    
    ifstream parfile;
    parfile.open(paramfile);
    
    string line;
    
    int nparams = 0;

    for ( int Or=0; Or<O2B_distance; Or++ ) 
    {
        for (int Oe=0; Oe<O2B_element; Oe++)
        {
            get_next_line(parfile, line);
            
            O2BD .push_back(Or);
            O2BE .push_back(Oe);
            coeff.push_back(stod(line));

            nparams++; 
        }
    }
    
    parfile.close();

}

// Frame class and member functions

class frame
{
    public:
        
        // Box properties
        
        vector<double> avec;
        vector<double> bvec;
        vector<double> cvec;
        
        double vol;
        
        vector<double> stress;   //  sxx syy szz sxy sxz syz
        
        double energy;           
        
        // Coordinates
        
        int natoms;
        vector<string> element;
        vector<int>    elemtype; // Type index for element. Assign via the hyperparametrs class
        vector<xyz>    coords;
        vector<xyz>    force;
        
        // Member functions
        
        frame();
        ~frame();
        
        void read_frame(ifstream & infile, hyperparams & hyperparams);
        double get_vol();
        void get_MIC_dist(int i, int j, xyz & rvec, double & rlen);
        
};

frame::frame()
{
    avec  .resize(3);
    bvec  .resize(3);
    cvec  .resize(3);
    stress.resize(9);
}

frame::~frame(){}

void frame::read_frame(ifstream & infile, hyperparams & hyperparams)
{
    // Line parsing vars
    
    string line;
    vector<string> tokens;
    int ntokens;
    int ntokens_read = 0;
    
    // Read number of atoms
    
    get_next_line(infile, line);
    natoms = stod(line);
    
    // Read box length line, determine format
    
    get_next_line(infile, line);
    ntokens = split_line(line,tokens);

    
    if (tokens[0] == "NON_ORTHO") // Read boxdims
    {
        avec[0] = stod(tokens[1]); avec[1] = stod(tokens[2]); avec[2] = stod(tokens[3]);
        bvec[0] = stod(tokens[4]); bvec[1] = stod(tokens[5]); bvec[2] = stod(tokens[6]);
        cvec[0] = stod(tokens[7]); cvec[1] = stod(tokens[8]); cvec[2] = stod(tokens[9]);
        ntokens_read += 10;
    }
    else
    {        
        avec[0] = stod(tokens[0]); avec[1] = 0.0            ; avec[2] = 0.0            ;
        bvec[0] = 0.0            ; bvec[1] = stod(tokens[1]); bvec[2] = 0.0            ;
        cvec[0] = 0.0            ; cvec[1] = 0.0            ; cvec[2] = stod(tokens[2]);   
        ntokens_read += 3;    
    }

    
    // Read the atom coordinate lintes
    
    element .resize(natoms);
    coords  .resize(natoms);
    force   .resize(natoms);
    elemtype.resize(natoms);

    for (int i=0; i<natoms; i++)
    {
        get_next_line(infile, line);
        split_line(line,tokens);
        
        element[i]   =      tokens[0];
        coords [i].x = stod(tokens[1]);
        coords [i].y = stod(tokens[2]);
        coords [i].z = stod(tokens[3]);
        force  [i].x = stod(tokens[4]) * 627.50961 * 1.889725989; // Read in H/B, convert to kcal/mol/Ang
        force  [i].y = stod(tokens[5]) * 627.50961 * 1.889725989;
        force  [i].z = stod(tokens[6]) * 627.50961 * 1.889725989;    
        
        elemtype[i] = hyperparams.assign_elemtype(element[i]);
        
        //cout << "Assigned element: " << element[i] << " type index: " << elemtype[i] << endl; 
                       
    }

}


void frame::get_MIC_dist(int i, int j, xyz & rvec, double & rlen) // Assumes an orthorhombic box for now
{
    double dx = coords[j].x - coords[i].x;
    double dy = coords[j].y - coords[i].y;
    double dz = coords[j].z - coords[i].z;
    
    rvec.x = dx - rint(dx / avec[0]) * avec[0];
    rvec.y = dy - rint(dy / bvec[1]) * bvec[1];
    rvec.z = dz - rint(dz / cvec[2]) * cvec[2];
    
    rlen = sqrt(rvec.x*rvec.x + rvec.y*rvec.y + rvec.z*rvec.z);
}


// Remining functions


void smooth(double rij, double rout, double & fs, double & fsderiv)
{
    fs      =  pow(1.0 - rij/rout,3.0); 
    fsderiv = -3.0*pow(rout - rij,2.0)  / pow(rout,3.0); 
}

void calc_2b(params & parameters, int etypei, int etypej)
{

	static double *Tr, *Trd;    // Polynomials for distance
    static double *Te;//, *Ted;     // Polynomials for ionization energy
	static bool called_before = false;

	double fcut; 
	double fcutderiv; 				
	double deriv;
	double tmp_doub; 	

    int parameter_idx;
    
    int stress_count;

	if ( ! called_before ) 
	{
		called_before = true;
		int dim = 0;

		dim  = parameters.max_order + 1;
		Te   = new double [dim]; 		Tr   = new double [dim];
		/*Ted  = new double [dim];*/ 		Trd  = new double [dim];
	}
    
    // Variables that get updated each loop:
    
    double morse_lambda, rcin, rout, sij, dsij_drij, sigmaij, fs, fsderiv;
    
    // Doing this the inefficient way - with a double loop
    


    int Oe, Or;
    double coeff;
    double tmp_ener = 0;

           
    parameters.universalize_dist(parameters.vdw[etypei], parameters.vdw[etypej],  morse_lambda, rcin, rout);
    
    double step = (rout - rcin)/100;
    double rlen;
    
for (int i=0; i<100; i++)
{
        rlen = rcin + step*i;
        
        //cout << rlen << " " << rcin << " " << rout << endl;
    
    if ( (rlen > rcin) && (rlen < rout) ) // No special handling for r < rcin, since not allowed. Should be included later
    {
        
        ////// Do distance transformation
        
        parameters.transform_dist(rlen, morse_lambda, rcin, rout, sij, dsij_drij);
        
   
        ////// Do element transformation (new)
        
        //cout << config.elemtype[i] << " --> " << parameters.eion[config.elemtype[i]] << endl;
            
        sigmaij = parameters.transform_eion(parameters.eion[etypei], parameters.eion[etypej]);
        
        //cout << "Calculated eion sigmaij: " << sigmaij << endl;
        
        
        ////// Setup the Chebyshev polynomials    

        Tr[0] = 1.0; Te[0] = 1.0;        // First two 1st-kind Chebys:
        Tr[1] = sij; Te[1] = sigmaij;

        Trd[0] = 1.0;       //Ted[0] = 1.0;       // Start the derivative setup. Set the first two 1st-kind Cheby's equal to the first two of the 2nd-kind
        Trd[1] = 2.0 * sij; //Ted[1] = 2.0 * sigmaij; 
	
        // Use recursion to set up the higher n-value Tn and Tnd's

        for ( int k = 2; k <= parameters.max_order; k++ ) 
        {
        	Tr[k]  = 2.0 * sij *  Tr[k-1] -  Tr[k-2];
        	Trd[k] = 2.0 * sij * Trd[k-1] - Trd[k-2];
            
        	Te[k]  = 2.0 * sigmaij *  Te[k-1] -  Te[k-2];
        	//Ted[k] = 2.0 * sigmaij * Ted[k-1] - Ted[i-2];
        }
        
        // For Tr polynomials, need to include rij transformation chain rule contribution

        for ( int k = parameters.max_order; k >= 1; k-- ) 
        	Trd[k] = k * dsij_drij * Trd[k-1];

        Trd[0] = 0.0;


        ////// Compute the smoothing function and its derivative to make life easy, just use the cubic form for now
        
        smooth(rlen, rout, fs, fsderiv);
        
        
        ////// Put it all together
        

        
        for (int p=0; p<parameters.coeff.size(); p++)
        {
            coeff = parameters.coeff[p];
            Oe    = parameters.O2BE[p];
            Or    = parameters.O2BD[p];

            tmp_doub = Te[Oe+1] * ( fs*Trd[Or+1] + fsderiv*Tr[Or+1]);   
            
            tmp_ener += fs * coeff * Te[Oe+1] * Tr[Or+1];
                    
            //To exclude eion dimension, comment out the line above and uncomment line below. Then, set Oe to 1 in main:
            //tmp_doub = ( fs*Trd[Or+1] + fsderiv*Tr[Or+1]); 
        }
    }

    cout << rlen<< " " << tmp_ener << endl;
}
    return;
}




int main(int argc, char* argv[])
{
    ifstream infile;

    params parameters(stoi(argv[1]), stoi(argv[2])); // Arguments: eion order, rij order

    parameters.read_params(argv[3]); // Parameter file - params.txt
    
    calc_2b(parameters, stoi(argv[4]), stoi(argv[5])); // Last two arguments: element type indices. 5 = He, 6 = Ar

    infile.close();

    return 0;
    
}