#include "nr3.h"

using namespace std;

struct Flow {
	const int dim = 65;
	const double v0 = 1.0;
	double w;
	const double delta = 1.0/double(dim-1);
	const double visc = 0.1;
	ofstream outfile; 

	const double plate_upstream = 0.25;
	const double plate_downstream = 0.375;
	const double plate_top = 0.25;

	MatDoub psi;
	MatDoub xi;
	MatDoub residuals_psi;
	MatDoub residuals_xi;

	Flow(double w_in) {
		w = w_in;
		xi.assign(dim, dim, 0.0); 
		psi.assign(dim, dim, 0.0);
		residuals_psi.assign(dim, dim, 0.0);
		residuals_xi.assign(dim, dim, 0.0);

		//
		//Initialize Psi to free flow state. 
		//
		for (int ell=0; ell<dim; ell++) {
			double y = delta*double(ell);

			double v0_y = v0 * y;

			for (int j=0; j<dim; j++) {

				double x = delta*double(j);


				//psi also 0 along and over plate surface
				if(x >= plate_upstream && x <= plate_downstream && y <= plate_top){
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
		cout << "\n\nPrinting Psi . . . \n";

		for (int j=dim-1; j>-1; j--) {
			for (int ell=0; ell<dim; ell++) {
				cout << setw(8)<< setprecision(2) << psi[ell][j]; 
			}
			cout << endl;
		}
	}

	void printXi() {
		cout << "\n\nPrinting Xi . . .\n";

		for (int j=dim-1; j>-1; j--) {
			for (int ell=0; ell<dim; ell++) {
				cout << setw(8)<< setprecision(2) << xi[ell][j]; 
			}
			cout << endl;
		}
	}

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

	void psi_update() {
		for (int j=1; j<dim-1; j++) {
			double x = delta*j;

			for (int ell=1; ell<dim; ell++) {	
				double y = delta * ell; 

				//
				//Make sure surface is set to 
				//	zero for all points. 
				// 
				if(x >= plate_upstream && x <= plate_downstream && y <= plate_top){
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
				//Update interior points...
				//
				else {
					psi[j][ell] = w * .25 * ( psi[j+1][ell] + psi[j-1][ell]
						+ psi[j][ell+1] + psi[j][ell-1] - delta * delta * xi[j][ell] )
					+ (1.0 - w ) * psi[j][ell];

					//
					//Upstream bdry: - dpsi/dx = 0. 
					//
					if (j==1) {
						psi[0][ell] = psi[1][ell];
					}
					
					//
					//Downstream bdry: dpsi/dx = 0.
					//
					else if (j == dim-2) {
						psi[j+1][ell] = psi[j][ell];

					}
				} 				
			}
		}
	}

	void xi_update() {

		//
		//Top of obstruction bdry condition. 
		//
		int el = int(plate_top/delta);
		for (int j=int(plate_upstream/delta); j < int(plate_downstream/delta); j++){
			xi[j][el] = psi[j][el+1] * 2.0 / (delta * delta);
		}

		//
		//Left and right obstruction bdry condition. 
		//
		int j_left = int(plate_upstream / delta);
		int j_right = int(plate_downstream / delta); 
		for (int ell = 0; ell <= int(plate_top / delta); ell++) {	
			xi[j_left][ell] = psi[j_left - 1][ell] * (2.0 / (delta*delta) );
			xi[j_right][ell] = psi[j_right+1][ell] * (2.0 / (delta*delta) );
		}

		//
		//Update interior points
		//
		for (int j=1; j<dim-1; j++) {
			double x = delta*j;

			for (int ell=1; ell<dim-1; ell++) {	
				double y = delta * ell; 

				//
				//If point is ON obstruction boundary, 
				//	continue. 
				//
				if ( ( (x == plate_upstream || x == plate_downstream ) && y <= plate_top)) {
					continue;
				}

				if (x >= plate_upstream && x <= plate_downstream && y == plate_top) {
					continue;
				}

				double dpsi_dx = (psi[j+1][ell]-psi[j-1][ell]) / (2.0 * delta);
				double dpsi_dy = (psi[j][ell+1]-psi[j][ell-1]) / (2.0 * delta);

				double dxi_dx = (xi[j+1][ell] - xi[j-1][ell]) / (2.0 * delta);
				double dxi_dy = (xi[j][ell+1] - xi[j][ell-1]) / (2.0 * delta);

				double source = (1.0 / visc) * (dpsi_dy * dxi_dx - dpsi_dx * dxi_dy);

				if(x > plate_upstream && x < plate_downstream && y < plate_top){
					xi[j][ell] = 0.0;
				}			
				else {
					xi[j][ell] = w * .25 * ( xi[j+1][ell] + xi[j-1][ell]
						+ xi[j][ell+1] + xi[j][ell-1] - delta * delta * source )
					+ (1.0 - w ) * xi[j][ell];
				}
			}
		}

		//
		//Downstream bdry condition: 
		//	dxi / dx = 0.
		//
		for (int ell=0; ell<dim; ell++) {
			xi[dim-1][ell] = xi[dim-2][ell];
		}


	}

	void residuals() {
		double residual_norm;
		int counter = 0;
		for (int j=0; j<dim; j++) {
			for (int ell=0; ell<dim; ell++) {
				double x = delta*double(j);
				double y = delta*double(ell);				
				//cout << "inner iterated: ell = " << ell << endl;
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
		//Should we be dividing by
		//	2.0 to get the derivative?!
		//
				else if (j == 0) {
					//cout << "on upstream. iteration no.: " << counter << "\n";
					counter++;
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

					double laplace_psi = ddpsi_ddx + ddpsi_ddy;
					residuals_psi[j][ell] = - xi[j][ell] + laplace_psi;

					double laplace_xi = ddxi_ddx + ddxi_ddy;
					residuals_xi[j][ell] = laplace_xi - source;
				}
			}
		}
	}

	double residualNorm() {
		double norm_psi_reimann=0.0;
		double norm_xi_reimann=0.0;

		//
		//Calculate integral using reimann sum. 
		//
		double integral=0.0;
		for (int j=0; j<dim-1; j++) {
			for (int ell=0; ell<dim-1; ell++) {
	

				double rave = 0.25*(residuals_psi[j][ell] + residuals_psi[j+1][ell] 
								+ residuals_psi[j][ell+1] + residuals_psi[j+1][ell+1]);

				integral += delta*delta*rave*rave;
			}
		}


		norm_psi_reimann = sqrt(integral); 
		cout << "Using reimann ... " << norm_psi_reimann << endl;

		//
		//Calculate integral using Root mean square. 
		//
		double sum=0.0;
		for (int j=0; j<dim; j++) {
			for (int ell=0; ell<dim; ell++) {
				sum += (residuals_psi[j][ell] * residuals_psi[j][ell])*delta*delta;
			}
		}
		sum = sqrt(sum);

		cout << "Using averaged... " << sum << endl; 

		return sum;

	}

	double sweep() {
		psi_update(); 
		xi_update(); 



		cout << "finished updates...\n";

		residuals();
		return residualNorm();
		 
	}

	void data_out() {
		outfile.open("data.dat");
		int width = 16;
		outfile.setf(ios::left);
		outfile << "#================================"
						<< "================================"
                        << "================================\n";
        outfile << "#" << setw(width) << "x" << setw(width) << "y"
        				<< setw(width) << "Psi" << setw(width) << "Res_Psi"
        				<< setw(width) << "Xi" << "Res_Xi\n";
		outfile << "#================================"
						<< "================================"
                        << "================================\n";

        for (int j = 0; j<dim; j++) {
        	double x = delta*double(j);

        	for (int ell=0; ell<dim; ell++) {
        		double y = delta*double(ell);

        		outfile << setw(width) << x << setw(width) << y
        				<< setw(width) << psi[j][ell] << setw(width) << residuals_psi[j][ell]
        				<< setw(width) << xi[j][ell] << setw(width) << residuals_xi [j][ell]
        				<< endl;
        	}
        	outfile << endl;
        }
        outfile.close();        			

	}

};

int main() {
	Flow flow(1.6);
	//flow.printPsi();
	cout << endl;
	//flow.printXi();

	cout<< endl;

	for(int patrick = 0; patrick < 10000; patrick++){
		cout << flow.sweep() << endl << endl;
		//flow.printPsi();
		cout << endl << endl;
		//flow.printXi();
		cout << "\n\n";
		//flow.printResiduals();
	}

	//flow.printPsi();
	//flow.printXi();
	flow.data_out();

	//flow.printPsi(); 
	cout<< endl;
	//flow.printXi();

	flow.printResiduals();


	return 0; 
}