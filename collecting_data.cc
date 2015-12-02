#include "nr3.h"

using namespace std;

//test space for fluid to move through
struct Flow {
    //number of verticies per side for grid of test space
    static const int dim = 65;
    //free flow velocity of fluid
    double v0;
    //relaxation factor
    double w;
    //will hold norm of all residuals for psi values
    double psi_resid_norm;
    //distance between verticies
    static const double delta = 1.0/double(dim-1);
    //viscosity value of fluid
    double visc;
    
    //location of rectangular obstruction in test space (given in x and y coordinates instead of j and l)
    static const double plate_upstream = 0.25;
    static const double plate_downstream = 0.375;
    static const double plate_top = 0.25;
    
    //matricies that will all psi, xi, and respective residual values
    MatDoub psi;
    MatDoub xi;
    MatDoub residuals_psi;
    MatDoub residuals_xi;
    
    //Constructor
    Flow(double w_in, double v0_in, double visc_in) {
        //set all parameters
        w = w_in;
        v0 = v0_in;
        visc= visc_in;
        //makes all interior points 0 unless otherwise corrected later
        xi.assign(dim, dim, 0.0);
        psi.assign(dim, dim, 0.0);
        residuals_psi.assign(dim, dim, 0.0);
        residuals_xi.assign(dim, dim, 0.0);
        

        //Initialize Psi to free flow state. (j and l correspond to grid)
        for (int ell=0; ell<dim; ell++) {
            double y = delta*double(ell);
            
            double v0_y = v0 * y;
            
            for (int j=0; j<dim; j++) {
                
                double x = delta*double(j);
                
                //psi also 0 along and over plate surface (boundary condition)
                if(x >= plate_upstream && x <= plate_downstream && y <= plate_top){
                    v0_y = 0.0;
                }
                //if not in/ on plate, make free flow
                else{
                    v0_y = v0*y;
                }
                
                psi[j][ell] = v0_y;
            }
        }
    }
    
    //world's most thrilling deconstructor
    ~Flow() {}
    
    
    /* Whenever running code withn new value of w/v/visc, means new test space (so
    things need to be reset, which is the only time this is called).
    Seemed less computationally expensive than making a new Flow object each time
    (but this is literally the code from constructor) */
    void reset_everything(){
        //makes all interior points 0 unless otherwise corrected later
        xi.assign(dim, dim, 0.0);
        psi.assign(dim, dim, 0.0);
        residuals_psi.assign(dim, dim, 0.0);
        residuals_xi.assign(dim, dim, 0.0);
        
        
        //Initialize Psi to free flow state. (j and l correspond to grid)
        for (int ell=0; ell<dim; ell++) {
            double y = delta*double(ell);
            
            double v0_y = v0 * y;
            
            for (int j=0; j<dim; j++) {
                
                double x = delta*double(j);
                
                //psi also 0 along and over plate surface (boundary condition)
                if(x >= plate_upstream && x <= plate_downstream && y <= plate_top){
                    v0_y = 0.0;
                }
                //if not in/ on plate, make free flow
                else{
                    v0_y = v0*y;
                }
                
                psi[j][ell] = v0_y;
            }
        }
    }
    
    //for updating any of the parameters we're testing
    void setW(double w_in){
        w = w_in;
        reset_everything();
    }
    
    void setVel(double v0_in){
        v0 = v0_in;
        reset_everything();
    }
    
    void setVisc(double visc_in){
        visc = visc_in;
        reset_everything();
    }
    
    
    //prints out all psi values
    void printPsi() {
        cout << "\n\nPrinting Psi . . . \n";
        
        for (int j=dim-1; j>-1; j--) {
            for (int ell=0; ell<dim; ell++) {
                cout << setw(8)<< setprecision(2) << psi[ell][j];
            }
            cout << endl;
        }
    }
    
    //prints out all xi values
    void printXi() {
        cout << "\n\nPrinting Xi . . .\n";
        
        for (int j=dim-1; j>-1; j--) {
            for (int ell=0; ell<dim; ell++) {
                cout << setw(8)<< setprecision(2) << xi[ell][j];
            }
            cout << endl;
        }
    }
    
    //prints all local residual values for psi and xi
    void printResiduals() {
        cout<< "\n\nResiduals for Psi: \n";
        for (int j=dim-1; j>-1; j--) {
            for (int ell=0; ell<dim; ell++) {
                cout << setw(8)<< setprecision(2) << residuals_psi[ell][j];
            }
            cout << endl;
        }
        
        cout << "\n\nResiduals for Xi: \n";
        for (int j=dim-1; j>-1; j--) {
            for (int ell=0; ell<dim; ell++) {
                cout << setw(8)<< setprecision(2) << residuals_xi[ell][j];
            }
            cout << endl;
        }
    }
    
    //uses boundary conditions and SOR to update all local psi values
    void psi_update() {
        for (int j=1; j<dim-1; j++) {
            double x = delta*j;
            
            for (int ell=1; ell<dim; ell++) {
                double y = delta * ell;

                //Make sure obstruction is set to
                //	zero for all points on it.
                if(x >= plate_upstream && x <= plate_downstream && y <= plate_top){
                    psi[j][ell] = 0.0;
                }
                
                //Implementing upper boundary
                //	condition.
                else if (ell == dim-1) {
                    psi[j][ell] = v0 * y;
                }
                

                //Update interior points...
                else {
                    psi[j][ell] = w * .25 * ( psi[j+1][ell] + psi[j-1][ell]
                                             + psi[j][ell+1] + psi[j][ell-1] - delta * delta * xi[j][ell] )
                    + (1.0 - w ) * psi[j][ell];
                    
                    //Upstream bdry: - dpsi/dx = 0.
                    //can't update this point independently so need to set it at same time
                    if (j==1) {
                        psi[0][ell] = psi[1][ell];
                    }
                    
                    //Downstream bdry: dpsi/dx = 0.
                    //can't update this point independently so need to set it at same time
                    else if (j == dim-2) {
                        psi[j+1][ell] = psi[j][ell];
                        
                    }
                }
            }
        }
    }
    
    //uses boundary conditions and SOR to update all local xi values
    void xi_update() {
        
        //Top of obstruction bdry condition.
        int el = int(plate_top/delta);
        for (int j=int(plate_upstream/delta); j < int(plate_downstream/delta); j++){
            xi[j][el] = psi[j][el+1] * 2.0 / (delta * delta);
        }
        
        //Left and right obstruction bdry condition.
        int j_left = int(plate_upstream / delta);
        int j_right = int(plate_downstream / delta);
        for (int ell = 0; ell <= int(plate_top / delta); ell++) {
            xi[j_left][ell] = psi[j_left - 1][ell] * (2.0 / (delta*delta) );
            xi[j_right][ell] = psi[j_right+1][ell] * (2.0 / (delta*delta) );
        }
        
        //Update interior points
        for (int j=1; j<dim-1; j++) {
            double x = delta*j;
            
            for (int ell=1; ell<dim-1; ell++) {
                double y = delta * ell;
                
                //If point is ON obstruction boundary (sides and top),
                //	continue.
                if ( ( (x == plate_upstream || x == plate_downstream ) && y <= plate_top)) {
                    continue;
                }
                
                if (x >= plate_upstream && x <= plate_downstream && y == plate_top) {
                    continue;
                }
                
                //using centered finite differencing to find derivatives
                double dpsi_dx = (psi[j+1][ell]-psi[j-1][ell]) / (2.0 * delta);
                double dpsi_dy = (psi[j][ell+1]-psi[j][ell-1]) / (2.0 * delta);
                
                double dxi_dx = (xi[j+1][ell] - xi[j-1][ell]) / (2.0 * delta);
                double dxi_dy = (xi[j][ell+1] - xi[j][ell-1]) / (2.0 * delta);
                
                double source = (1.0 / visc) * (dpsi_dy * dxi_dx - dpsi_dx * dxi_dy);
                
                //xi is 0 inside obstruction always
                if(x > plate_upstream && x < plate_downstream && y < plate_top){
                    xi[j][ell] = 0.0;
                }
                //using SOR
                else {
                    xi[j][ell] = w * .25 * ( xi[j+1][ell] + xi[j-1][ell]
                                            + xi[j][ell+1] + xi[j][ell-1] - delta * delta * source )
                                            + (1.0 - w ) * xi[j][ell];
                }
            }
        }
        
        //Downstream bdry condition:
        //	dxi / dx = 0.
        for (int ell=0; ell<dim; ell++) {
            xi[dim-1][ell] = xi[dim-2][ell];
        }
        
        
    }
    
    
    /*calculates residual values for psi and xi by measuring
    difference between lhs and rhs of equations given in write up*/
    void residuals() {
        for (int j=0; j<dim; j++) {
            for (int ell=0; ell<dim; ell++) {
                double x = delta*double(j);
                double y = delta*double(ell);
  
                //
                //Centerline boundaries:
                //
                if (ell==0){
                    residuals_psi[j][ell] = psi[j][ell];
                    residuals_xi[j][ell] = xi[j][ell];
                }
                //
                //Downstream obstruction bdry.
                //
                else if (j == int( plate_downstream / delta)
                         && ell < int(plate_top / delta) ) {
                    
                    residuals_psi[j][ell] = psi[j][ell];
                    residuals_xi[j][ell] = xi[j][ell] - (2.0 / (delta * delta) ) *
                    psi[j+1][ell];
                }
                //
                //Upstream obstruction bdry.
                //
                else if (j == int( plate_upstream / delta)
                         && ell < int(plate_top / delta) ) {
                    
                    residuals_psi[j][ell] = psi[j][ell];
                    residuals_xi[j][ell] = xi[j][ell] - (2.0 / (delta * delta) ) *
                    psi[j-1][ell];
                }
                //
                //Top obstruction bdry.
                //
                else if (ell == int(plate_top / delta) && j >= int(plate_upstream / delta)
                         && j <= int(plate_downstream / delta) ) {
                    residuals_psi[j][ell] = psi[j][ell];
                    residuals_xi[j][ell] = xi[j][ell] - (2.0 / (delta*delta) ) *
                    psi[j][ell+1];
                }
                //
                //Upstream Bdry.
                //
                else if (j == 0) {
                    residuals_psi[j][ell] = - (psi[j+1][ell] - psi[j][ell]) / ( delta);
                    residuals_xi[j][ell] = xi[j][ell];
                }
                //
                //Downstream Bdry.
                //
                else if (j == dim-1) {
                    residuals_psi[j][ell] = - (psi[j-1][ell] - psi[j][ell]) / (delta);
                    residuals_xi[j][ell] = xi[j][ell];
                }
                //
                //Top bdry.
                //
                else if (ell == dim-1) {
                    residuals_psi[j][ell] = (psi[j][ell] - psi[j][ell-1]) / (delta)
                    - v0;
                    residuals_xi[j][ell] = xi[j][ell];
                }
                //
                //Inside of obstruction.
                //
                else if (x > plate_upstream && x < plate_downstream && y < plate_top) {
                    residuals_xi[j][ell] = xi[j][ell];
                    residuals_psi[j][ell]=psi[j][ell];
                }
                //
                //Interior points
                //
                else {
                    //
                    //Compute second derivatives...
                    //
                    double ddpsi_ddx = (psi[j+1][ell] - 2.0*psi[j][ell] + psi[j-1][ell])
                    / (delta*delta);
                    double ddpsi_ddy = (psi[j][ell+1] - 2.0*psi[j][ell] + psi[j][ell-1])
                    / (delta*delta);
                    double ddxi_ddx = (xi[j+1][ell] - 2.0*xi[j][ell] + xi[j-1][ell])
                    / (delta*delta);
                    double ddxi_ddy = (xi[j][ell+1] - 2.0*xi[j][ell] + xi[j][ell-1])
                    / (delta*delta);
                    
                    //
                    //...and first derivatives.
                    //
                    double dpsi_dx = (psi[j+1][ell]-psi[j-1][ell]) / (2.0 * delta);
                    double dpsi_dy = (psi[j][ell+1]-psi[j][ell-1]) / (2.0 * delta);
                    
                    double dxi_dx = (xi[j+1][ell] - xi[j-1][ell]) / (2.0 * delta);
                    double dxi_dy = (xi[j][ell+1] - xi[j][ell-1]) / (2.0 * delta);
                    
                    double source = (1.0 / visc) * (dpsi_dy * dxi_dx - dpsi_dx * dxi_dy);
                    
                    //and now finally appy that to getting residuals
                    //(see write up for where these equations come from)
                    double laplace_psi = ddpsi_ddx + ddpsi_ddy;
                    residuals_psi[j][ell] = - xi[j][ell] + laplace_psi;
                    
                    double laplace_xi = ddxi_ddx + ddxi_ddy;
                    residuals_xi[j][ell] = laplace_xi - source;
                }
            }
        }
    }
    
    //finds norm of residuals of psi
    void residualNormPsi() {
        //Calculate integral using Root mean square.
        double sum=0.0;
        for (int j=0; j<dim; j++) {
            for (int ell=0; ell<dim; ell++) {
                sum += (residuals_psi[j][ell] * residuals_psi[j][ell])*delta*delta;
            }
        }
        sum = sqrt(sum);
        psi_resid_norm = sum;
    }
    
    //runs a single sweep
    void one_sweep() {
        psi_update();
        xi_update();
        
        residuals();
        residualNormPsi();
    }
    
    //does the sweeps here so that way it's only one function call from main
    double allTheSweeps() {
        int repeats = 0;
        double lastNorm;
        
        for(int sweep = 0; sweep < 10000; sweep++){
            lastNorm = psi_resid_norm;
            one_sweep();
            
            //checks for converging
            if(psi_resid_norm == lastNorm){
                repeats++;
                
                //if residual has converged (avoiding extra computation time)
                //stop doing sweeps and return value
                if(repeats > 99){
                    
                    cout << "converged early at " << sweep << " iteration\n";
                    return psi_resid_norm;
                }
            }
            //reset amount of times value has repeated if residual has updated
            else{
                repeats = 0;
            }
        }
        
        //otherwise current value is returned after 10,000 sweeps
        return psi_resid_norm;
    }
    
    
};

int main() {
    
    //will hold current residual so it can be written to file
    double resid;
    //width for each collumn in data file
    int width = 16;
    
    //values from example
    static const double w_eg = 1.5;
    static const double vel_eg = 1.0;
    static const double visc_eg = 0.1;
    
    
    //creates original flow object with values from example
    Flow flow(w_eg,  vel_eg, visc_eg);
    
    
    //evaluates just w values
    ofstream outfile;
    outfile.open("w_values.dat");
    outfile.setf(ios::left);
    outfile << "#================================"
    << "================================================================\n";
    outfile << "#" << setw(width) << "w" << setw(width) << "residual norm\n";
    
    for(double w = 1.0; w <= 2.0; w+= 0.1){
        cout << "\n\n";
        flow.setW(w);
        resid = flow.allTheSweeps();
        outfile << setw(width) << w << setw(width) << resid << endl;
    }
    outfile.close();
    

    
    //return back to example value
    Flow flow2(w_eg,  vel_eg, visc_eg);
    
    
    //evaluates just velocity values
    ofstream out;
    out.open("v_values.dat");
    out.setf(ios::left);
    out << "#================================"
    << "================================================================\n";
    out << "#" << setw(width) << "velocity" << setw(width) << "residual norm\n";
    
    for(double v = 0.5; v <= 10.0; v+= 0.5){
        flow2.setVel(v);
        resid = flow2.allTheSweeps();
        out << setw(width) << v << setw(width) << resid << endl;
    }
    out.close();

    
    
    //return back to example value
    Flow flow3(w_eg,  vel_eg, visc_eg);
    
    
    
    //evaluates just viscosity values
    ofstream outf;
    outf.open("visc_values.dat");
    outf.setf(ios::left);
    outf << "#================================"
    << "================================================================\n";
    outf << "#" << setw(width) << "viscosity" << setw(width) << "residual norm\n";
    
    for(double visc = 0.5; visc <= 10.0; visc+= 0.5){
        flow3.setVisc(visc);
        resid = flow3.allTheSweeps();
        outf << setw(width) << visc << setw(width) << resid << endl;
    }
    outf.close();
    
    //return back to example value
    Flow flow4(w_eg,  vel_eg, visc_eg);
    
    
    
    //now evaluates how all variables affect each other
    //by looping through all possible combinations
    ofstream otf;
    otf.open("all_the_values.dat");
    otf << "#================================"
    << "================================================================\n";
    otf << "#" << setw(width) << "w" << setw(width) << "velocity" << setw(width) << "viscosity" << setw(width) << "residual norm\n";
    
    for(double w = 1.0; w <= 2.0; w+= 0.1){
        flow4.setW(w);
        
        for(double v = 0.5; v <= 10.0; v+= 0.5){
            flow4.setVel(v);
            
            for(double visc = 0.5; visc <= 10.0; visc+= 0.5){
                flow4.setVisc(visc);
                resid = flow4.allTheSweeps();
                otf << setw(width) << w << setw(width) << v << setw(width) << visc << setw(width) << resid << endl;
            }
        }
    }
    
    otf.close();
    
    return 0; 
}
