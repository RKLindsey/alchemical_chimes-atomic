/* 
This is a script to test out generation of generalized ChIMES model that does not require explicit parmeters for different element types
R. K. Lindsey (2024)

Got ionization energies and vdw radii from: https://pubchem.ncbi.nlm.nih.gov/element/89

Note: Setting Oe to 1 should return behavior of the normal ChIMES code, when used with a cubic smoothing function. 
Cutoffs will need to be adjusted to match what this code uses as well.

Run this code with:
g++ <this script>
./a.out <.xyzf file> <number of frames> <Oe> <Or>
Produces files A.txt and b.txt

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

///////////////////////////////////
// Helper utilities/objects
///////////////////////////////////

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


///////////////////////////////////
// Hyperparameters class and member functions
///////////////////////////////////

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


// This constructor defines the element descriptors, e.g., particles size (2*vdw radius) and energy descriptor, 
// which is currently taken as the ionization energy. Alternative energy descriptors are provided in comments
//  Note that these values are only hard-coded for a small number of atom types at the moment (see the function definition)
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

// THIS FUNCTION HAS HARD-CODED CHOICES THAT SHOULD BE EVALUATED AND MADE INTO OPTIONS
void hyperparams::universalize_dist(double vdwi, double vdwj, double & morse_lambda, double & rcin, double & rout)
{    
    double sigma_i      = vdwi;
    double sigma_j      = vdwj;
    double sigma_ij     = 0.5 * (sigma_i + sigma_j);
    
    morse_lambda = pow(2.0,1.0/6.0)*sigma_ij;
    rout  = 3.5  * sigma_ij;                            // <-- 3.5  = HARDCODED CHOICE
    rcin  = 0.65*morse_lambda; //0.45 * morse_lambda;   // <-- 0.65 = HARDCODED CHOICE
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

///////////////////////////////////
// Frame class and member functions
///////////////////////////////////

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

    // Calculate box volume....
    
    vol = get_vol();
    
    if (ntokens-ntokens_read >= 6) // Read stresses. Expects as:xx yy zz xy xz yz // old: xx, xy, xz, yy, yx, yz, zx, zy, zz
    {
        for(int i=0; i<6; i++)
            stress[0] = stod(tokens[ntokens_read+i])/6.9476955; // Convert from GPa to kcal/mol/Ang^3
        
        ntokens_read += 9;
    }
    else if (ntokens-ntokens_read >= 3)
    {
        cout << "OVER HERE" << endl;
        
        for(int i=0; i<3; i++)
            stress[0] = stod(tokens[ntokens_read+i])/6.9476955; // Convert from GPa to kcal/mol/Ang^3
        
        ntokens_read += 3;
    }

    if (ntokens-ntokens_read > 0) // Read in kcal/mol
        energy = stod(tokens[ntokens_read]);
    
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

double frame::get_vol() // Assume an orthorhombic box for now
{
    return avec[0]*bvec[1]*cvec[2];
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

///////////////////////////////////
// Amat class and member functions
///////////////////////////////////

class design_matrix
{
    public:
        
        vector<vector<double > > amat;        
        vector<double> bvec;
        
        double rows;
        double cols;
        
        design_matrix();
        ~design_matrix();
        void add_frame(vector<vector<double> > & tmp_avec, vector<double> & tmp_bvec);
        void print_amat();
        void print_bvec();
};

design_matrix::design_matrix(){}

design_matrix::~design_matrix(){}

void design_matrix::add_frame(vector<vector<double> > & tmp_avec, vector<double> & tmp_bvec)
{
    amat.insert(amat.end(), tmp_avec.begin(), tmp_avec.end());
    bvec.insert(bvec.end(),tmp_bvec.begin(), tmp_bvec.end());
}

void design_matrix::print_amat()
{
    ofstream outfile; 
    
    rows = amat.size();
    cols = amat[0].size();
        
    outfile.open("A.txt");
    
    for(int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
            outfile << amat[i][j] << " ";
        outfile << endl;
    }
    
    outfile.close();
}

void design_matrix::print_bvec()
{
    rows = bvec.size();
    
    ofstream outfile; 
    
    outfile.open("b.txt");
    
    for(int i=0; i<rows; i++)
            outfile << bvec[i] << endl;
        
    outfile.close();
}


///////////////////////////////////
// Remining functions
///////////////////////////////////

void smooth(double rij, double rout, double & fs, double & fsderiv)
{
    fs      =  pow(1.0 - rij/rout,3.0); 
    fsderiv = -3.0*pow(rout - rij,2.0)  / pow(rout,3.0); 
}

void deriv_2b(design_matrix  & design, hyperparams & hyper, frame & config, int fit_stress, bool fit_energy)
{
    /* Note: Don't need the derivative of the Te polynomial. Only the following derivatives matter here:
    
    dE/dc where c are the fitting coefficients - this is used in the A-matrix for energy
    
    dE/dr where r is the distance - this is the definition of negative force
    
    d/dc (dE/dr) - this is used in the A-matrix for force
    
    Nowhere in these expressions does d/de factor in, where e is the energy dimension descriptor (e.g., ionization energy)
    */
    
	xyz    rvec; 
	double rlen;
    
	static double *Tr, *Trd;    // Polynomials for distance
    static double *Te;//, *Ted;     
	static bool called_before = false;

	double fcut; 
	double fcutderiv; 				
	double deriv;
	double tmp_doub; 	

	double inv_vol = 1.0 / config.vol;
    
    int parameter_idx;
    
    int stress_count;
    int force_count = 3*config.natoms;
    

	if ( ! called_before ) 
	{
		called_before = true;
		int dim = 0;

		dim  = hyper.max_order + 1;
		Te   = new double [dim]; 		Tr   = new double [dim];
		/*Ted  = new double [dim];*/ 	Trd  = new double [dim];
	}
    
    // Variables that get updated each loop:
    
    double morse_lambda, rcin, rout, sij, dsij_drij, sigmaij, fs, fsderiv;
    
    // Doing this the inefficient way - with a double loop
    
    
    // Generate temporary amat/bmat
    
    double new_rows = config.natoms*3.0;
    
    if (fit_energy)
        new_rows += 1;
    if (fit_stress == 1)
        new_rows += 3;
    if (fit_stress == 2)
        new_rows += 6;
    
    vector<vector<double> > tmp_avec(new_rows, vector<double>(hyper.nparams,0));
    vector<double>          tmp_bvec(new_rows);
   
    // Populate the amat

    for (int i=0; i<config.natoms; i++)
    {    
        tmp_bvec[3*i+0] = config.force[i].x;
        tmp_bvec[3*i+1] = config.force[i].y;
        tmp_bvec[3*i+2] = config.force[i].z;       
        
        for (int j=i+1; j<config.natoms; j++)
        {
            config.get_MIC_dist(i, j, rvec, rlen); 
                        
            hyper.universalize_dist(hyper.vdw[config.elemtype[i]], hyper.vdw[config.elemtype[j]],  morse_lambda, rcin, rout);
            
            if ( (rlen > rcin) && (rlen < rout) ) // No special handling for r < rcin, since not allowed. Should be included later
            {
                parameter_idx = -1;
                force_count = 0;
                stress_count = 0;
                
                ////// Do distance transformation
                
                hyper.transform_dist(rlen, morse_lambda, rcin, rout, sij, dsij_drij);
                
   
                ////// Do element transformation (new)
                             
                sigmaij = hyper.transform_eion(hyper.eion[config.elemtype[i]], hyper.eion[config.elemtype[j]]);
                                
                ////// Setup the Chebyshev polynomials    

                Tr[0] = 1.0; Te[0] = 1.0;        // First two 1st-kind Chebys:
                Tr[1] = sij; Te[1] = sigmaij;

                Trd[0] = 1.0;       //Ted[0] = 1.0;       // Start the derivative setup. Set the first two 1st-kind Cheby's equal to the first two of the 2nd-kind
                Trd[1] = 2.0 * sij; //Ted[1] = 2.0 * sigmaij; 
	
                // Use recursion to set up the higher n-value Tn and Tnd's

                for ( int k = 2; k <= hyper.max_order; k++ ) 
                {
                	Tr[k]  = 2.0 * sij *  Tr[k-1] -  Tr[k-2];
                	Trd[k] = 2.0 * sij * Trd[k-1] - Trd[k-2];
                    
                	Te[k]  = 2.0 * sigmaij *  Te[k-1] -  Te[k-2];
                	//Ted[k] = 2.0 * sigmaij * Ted[k-1] - Ted[i-2];
                }
                
                // For Tr polynomials, need to include rij transformation chain rule contribution

                for ( int k = hyper.max_order; k >= 1; k-- ) 
                	Trd[k] = k * dsij_drij * Trd[k-1];

                Trd[0] = 0.0;


                ////// Compute the smoothing function and its derivative to make life easy, just use the cubic form for now
                
                smooth(rlen, rout, fs, fsderiv);
                
                
                ////// Put it all together
                
                
                for ( int Or=0; Or<hyper.O2B_distance; Or++ ) 
                {
                    for (int Oe=0; Oe<hyper.O2B_element; Oe++)
                    {
                        parameter_idx++;
                        
                        // tmp_doub is d/dc (d/dr), i.e. the derivative of the force magnitude with repsect to coefficients, which is 
                        // the negative derivative of energy with respect to distance and coefficients
                        
                        tmp_doub = -1.0 * Te[Oe+1] * ( fs*Trd[Or+1] + fsderiv*Tr[Or+1]);   
                        
    					// Finally, account for the x, y, and z unit vectors

    					deriv = tmp_doub * rvec.x / rlen;
    					tmp_avec[3*i+0][parameter_idx] += deriv;    // fx on atom i
    					tmp_avec[3*j+0][parameter_idx] -= deriv;    // fx on atom j

    					deriv = tmp_doub * rvec.y / rlen; 
    					tmp_avec[3*i+1][parameter_idx] += deriv;
    					tmp_avec[3*j+1][parameter_idx] -= deriv;

    					deriv = tmp_doub * rvec.z / rlen;
    					tmp_avec[3*i+2][parameter_idx] += deriv;
    					tmp_avec[3*j+2][parameter_idx] -= deriv;

    					if (fit_stress == 1)    //  STRESS INCLUSION UNTESTED!
    					{
    						tmp_avec[force_count+0][parameter_idx] -= tmp_doub * rvec.x * rvec.x / rlen;
    						tmp_avec[force_count+2][parameter_idx] -= tmp_doub * rvec.y * rvec.y / rlen;
    						tmp_avec[force_count+3][parameter_idx] -= tmp_doub * rvec.z * rvec.z / rlen;    
                            stress_count = 3;
    					}
 
    					else if (fit_stress == 2)    //  STRESS INCLUSION UNTESTED!
    					{
    						tmp_avec[force_count+0][parameter_idx] -= tmp_doub * rvec.x * rvec.x / rlen;   // xx
    						tmp_avec[force_count+2][parameter_idx] -= tmp_doub * rvec.x * rvec.y / rlen;   // xy
    						tmp_avec[force_count+3][parameter_idx] -= tmp_doub * rvec.x * rvec.z / rlen;   // xz

    						tmp_avec[force_count+0][parameter_idx] -= tmp_doub * rvec.y * rvec.y / rlen;   // yy
    						tmp_avec[force_count+2][parameter_idx] -= tmp_doub * rvec.y * rvec.z / rlen;   // yz
    						tmp_avec[force_count+3][parameter_idx] -= tmp_doub * rvec.z * rvec.z / rlen;   // zz   
                            stress_count = 6;
    					} 
 
    					if(fit_energy)    //  ENERGY INCLUSION UNTESTED!
    						tmp_avec[force_count+stress_count][parameter_idx] += fs * Te[Oe+1] * Tr[Or+1];
 
                     }
                }
                
                if (i==1)
                {
                    if ( (parameter_idx+1) != (hyper.nparams) )
                    {
                        cout << "ERROR: Parameter numbermismatch";
                        cout << "Expected " << hyper.nparams << " parameters." << endl;
                        cout << "Counted  " << parameter_idx+1 << " parameters." << endl;
                        exit(0);
                    }
                }
            }
        }
    }
    

    
    if (fit_stress == 1)
    {	 
    	 for ( int i = 0; i < hyper.nparams; i++ ) 
    	 {
    		 tmp_avec[force_count+0][parameter_idx] *= inv_vol;
    		 tmp_avec[force_count+2][parameter_idx] *= inv_vol;
    		 tmp_avec[force_count+3][parameter_idx] *= inv_vol;	
    	 }
    }
    else if (fit_stress == 2)
    {		
    	for ( int i = 0; i < hyper.nparams; i++ ) 
    	{
    		tmp_avec[force_count+0][parameter_idx] *= inv_vol;
    		tmp_avec[force_count+2][parameter_idx] *= inv_vol;
    		tmp_avec[force_count+3][parameter_idx] *= inv_vol;    
    		
    		tmp_avec[force_count+0][parameter_idx] *= inv_vol;
    		tmp_avec[force_count+2][parameter_idx] *= inv_vol;
    		tmp_avec[force_count+3][parameter_idx] *= inv_vol;
    	}
    }

    design.add_frame(tmp_avec, tmp_bvec);

    return;
}




int main(int argc, char* argv[])
{
    ifstream infile;    
    infile.open(argv[1]);//"test-Ar+He_50F.xyzf"
    int  nframes    = stoi(argv[2]); // 50;
    int  fit_stress = 0;
    bool fit_energy = 0;

    hyperparams hyperparam(stoi(argv[3]) , stoi(argv[4])); // Arguments: eion order, rij order
    
    design_matrix design;
    
    for(int i=0; i<nframes; i++)
    {
        cout << "Working on frame " << i << endl;
        
        frame config;
        config.read_frame(infile, hyperparam);
        
        cout << "   ...reading frame complete" << endl; 
        
        deriv_2b(design, hyperparam, config, fit_stress, fit_energy);
        
        cout << "   ...frame deriv calculation complete" << endl; 

    }
    
    // Print out design matrix
    
    design.print_amat();
    design.print_bvec();
        
    infile.close();

    return 0;
}