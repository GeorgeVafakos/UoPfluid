#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
	// Define variables
	int i,SIZE=5, D, Nx[SIZE], Ny[SIZE];
	float r, ro, r1, Fs, GCI_12, GCI_23, GCI_est, r_ast;
	float f[SIZE], f_RE[SIZE], p[SIZE], GCI_fine[SIZE], GCI_coarse[SIZE], ratio[SIZE];
	
	// Grid parameters
	D=2;													// Dimensions of problem
	SIZE=3;													// Number of test cases
	Fs=1.25;													// Factor of safety
	
	Nx[0]=168;												// Number of cells in x_direction in the 1st grid
	Ny[0]=126;												// Number of cells in y_direction in the 1st grid
	Nx[1]=202;												// Number of cells in x_direction in the 2nd grid
	Ny[1]=152;												// Number of cells in y_direction in the 2nd grid
	Nx[2]=242;												// Number of cells in x_direction in the 3rd grid
	Ny[2]=182;												// Number of cells in y_direction in the 3rd grid
	
	
	
	f[0]=1.004229;											// The quantity value of the 1st grid
	f[1]=1.004972;											// The quantity value of the 2nd grid
	f[2]=1.006468;											// The quantity value of the 3rd grid

		
	// Calculate effective refinement ratio
	r=pow(((1.0*Nx[1]*Ny[1])/(Nx[0]*Ny[0])),(1.0/D));
	
	for (i=2; i<SIZE; i++)
	{
		// Calculate order of grid convergencer
		p[i]=log((f[i-2]-f[i-1])/(f[i-1]-f[i]))/log(r);
		
		// p-th order Richardson Extrapolation
		f_RE[i]=f[i] + (f[i]-f[i-1])/(pow(r,p[i])-1);
		
		// Grid Convergence Index on fine grid (GCI-fine)
		GCI_fine[i]=100*Fs*fabs((f[i-1]-f[i])/f[i])/(pow(r,p[i])-1);
		
		// Grid Convergence Index on coarse grid (GCI-coarse)
		GCI_coarse[i]=100*Fs*fabs((f[i-1]-f[i-2])/f[i-1])/(pow(r,p[i])-1);
		
		// Asymptotic ratio of grid convergence
		ratio[i]=GCI_coarse[i]/(pow(r,p[i])*GCI_fine[i]);
	}
	
	// Print results
	printf("--- GRID CONVERGENCE STUDY ---\n\n");
	printf("Number of data sets read \t SIZE=%d\n",SIZE);
	printf("Effective Refinement Ratio \t r=%.2f\n",r);
	printf("Factor of Safety \t \t Fs=%.2f\n",Fs);
	printf("\n");
	printf("\t Grid Size \t\t Quantity \t Order p \t Rich. Extr. \t GCI-fine(%%) \t CGI-coarse(%%) \t Conv. Ratio \t Change(%%)\n");
	for (i=0 ; i<SIZE ; i++)
	{
		if (i==0 || i==1)
			printf("%d.\t %dx%d \t \t %.4f \t \t \t \t \t \t \t \t \t \t \t %.2f\n",i+1,Nx[i],Ny[i],f[i],100*fabs(f[SIZE-1]-f[i])/fabs(f[SIZE-1]));
		else
			printf("%d.\t %dx%d \t \t %.4f \t %.2f \t \t %.4f \t %.2f \t \t %.2f \t \t %.4f \t %.2f\n",i+1,Nx[i],Ny[i],f[i],p[i],f_RE[i],GCI_fine[i],GCI_coarse[i],ratio[i],100*fabs(f[SIZE-1]-f[i])/fabs(f[SIZE-1]));
	}
	printf("\n\n");
	printf("Parameter Explanation\n\n");
	printf("Grid Size   : \tThe number of points or cells in x,y direction.\n");
	printf("Order p     : \tThe order of the frid convergence calculated on the coarse, medium and fine grid.\n");
	printf("Rich. Extr. : \tThe value of the quantity for zero grid spacing, as calculated from the Generalized p-th order Richardson\n");
	printf("\t \tExtrapolation. Asymptorically, as the grid gets thicher, the value of the quantity must be equivalent to this value.\n");
	printf("GCI-fine    : \tThe Crid Convergence Index calculated from the finer and medium grid.\n");
	printf("GCI-coarse  : \tThe Crid Convergence Index calculated from the medium and coarse grid.\n");
	printf("Conv. Ratio : \tThe asymptotic range of convergenve. A ratio of 1.0 indicates asymptotic range.\n");
	printf("Change(%%)   : \tA parmeter that indicates the percentage change of the quantity on each grid compared to the value of the\n");
	printf("\t \tquantity on the denser grid.\n");
	printf("\n");
	printf("\n");
	printf("\n--- END OF STUDY ---\n\n");

}
