#include "nr3.h"

using namespace std;

struct Flow {
	const int dim = 17;
	const double v0 = 1.0;
	double w;
	const double delta = 1.0/double(dim-1);
	const double visc = 0.1;

	const double plate_lower_x = 0.25;
	const double plate_upper_x = 0.375;
	const double plate_top = 0.25;

	MatDoub psi;
	MatDoub xi;

	Flow(double w_in) {
		w = w_in;
		xi.assign(dim, dim, 0.0); 
		psi.assign(dim, dim, 0.0);

		//
		//Initialize Psi to free flow state. 
		//
		for (int ell=0; ell<dim; ell++) {
			double y = delta*double(ell);

			double v0_y = v0 * y;

			for (int j=0; j<dim; j++) {

				double x = delta*double(j);


				//psi also 0 along and over plate surface
				if(x >= plate_lower_x && x <= plate_upper_x && y <= plate_top){
					v0_y = 0.0;

				}else{
					v0_y = v0*y;
				}


				psi[j][ell] = v0_y;
			}
		}



	}

	~Flow() {}

	void printPsi() {
		for (int j=dim-1; j>-1; j--) {
			for (int ell=0; ell<dim; ell++) {
				cout << setw(8)<< setprecision(2) << psi[ell][j]; 
			}
			cout << endl;
		}
	}

	void printXi() {
		for (int j=dim-1; j>-1; j--) {
			for (int ell=0; ell<dim; ell++) {
				cout << setw(8)<< setprecision(2) << xi[ell][j]; 
			}
			cout << endl;
		}
	}

	void psi_update() {
		for (int j=1; j<dim-1; j++) {
			double x = delta*j;

			for (int ell=1; ell<dim; ell++) {	
				double y = delta * ell; 

				//
				//Make sure surface is set to 
				//	zero for all points. 
				// 
				if(x >= plate_lower_x && x <= plate_upper_x && y <= plate_top){
					psi[j][ell] = 0.0;
				}		
				//
				//Implementing upper boundary 
				//	condition. 
				//
				else if (ell == dim-1) {
					psi[j][ell] = v0 * y;
				}	
				//
				//Update...
				//
				else {
					psi[j][ell] = w * .25 * ( psi[j+1][ell] + psi[j-1][ell]
						+ psi[j][ell+1] + psi[j][ell-1] - delta * delta * xi[j][ell] )
					+ (1.0 - w ) * psi[j][ell];
				}


			}
		}
	}

	void xi_update() {
		for (int j=1; j<dim-1; j++) {
			double x = delta*j;

			for (int ell=1; ell<dim-1; ell++) {	
				double y = delta * ell; 

				double dpsi_dx = (psi[j+1][ell]-psi[j-1][ell]) / (2.0 * delta);
				double dpsi_dy = (psi[j][ell+1]-psi[j][ell-1]) / (2.0 * delta);

				double dxi_dx = (xi[j+1][ell] - xi[j-1][ell]) / (2.0 * delta);
				double dxi_dy = (xi[j][ell+1] - xi[j][ell-1]) / (2.0 * delta);

				double source = (1.0 / visc) * (dpsi_dy * dxi_dx - dpsi_dx * dxi_dy);

				if(x >= plate_lower_x && x <= plate_upper_x && y <= plate_top){
					xi[j][ell] = 0.0;
				}			
				else {
					xi[j][ell] = w * .25 * ( xi[j+1][ell] + xi[j-1][ell]
						+ xi[j][ell+1] + xi[j][ell-1] - delta * delta * xi[j][ell] )
					+ (1.0 - w ) * delta * delta * source;
				}

			}
		}

		for (int ell=0; ell<dim; ell++) {
			xi[dim-1][ell] = 0.0;
		}


	}

	void sweep() {
		psi_update(); 
		xi_update(); 

	}


};

int main() {
	Flow flow(1.5);
	flow.printPsi();
	cout << endl;
	flow.printXi();

	cout<< endl;

	flow.sweep();

	flow.printPsi(); 
	cout<< endl;
	flow.printXi();


	return 0; 
}
