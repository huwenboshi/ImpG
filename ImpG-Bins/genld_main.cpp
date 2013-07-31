#include <stdio.h>
#include <cmath>

#include "genld.h"

int main(int argc, char **argv) {
	// get command line
	char *IN_TYPED_SNP_FILE = NULL;
	char *IN_GEN_FILE = NULL, *OUT_LD_FILE = NULL;
	get_gen_ld_cmd_line(argc, argv, &IN_TYPED_SNP_FILE, 
		&IN_GEN_FILE, &OUT_LD_FILE);

	// load typed snps
	vector<typed_snp> typed_snps;
	size_t num_typed_snps = load_typed_snps(IN_TYPED_SNP_FILE, typed_snps);

	// load the genotype file
	vector<string> gens;
	size_t num_indv = load_genotypes(IN_GEN_FILE, num_typed_snps, gens);
	
	// compute the correlation between typed snps
	FILE *fout;
	fout = safe_fopen(OUT_LD_FILE, "w");
	fprintf(fout, "LD matrix for %s\n", IN_GEN_FILE);
	
	// iterate through snp i
	for(size_t i = 0; i < num_typed_snps; i++) {
		// iterate through snp j
		for(size_t j = i+1; j < num_typed_snps; j++) {
			double x11 = 0.0, x12 = 0.0, x21 = 0.0, x22 = 0.0;
			// iterate through individuals
			for(size_t k = 0; k < num_indv; k++) {
				char gi = gens[k][i];
				char gj = gens[k][j];
				if(gi == '0' && gj == '0') {
					x11 += 2.0;
				}
				else if (gi == '1' && gj == '1') {
					x11 += 1.0;
					x22 += 1.0;
				}
				else if (gi == '2' && gj == '2') {
					x22 += 2.0;
				}
				else if (gi == '0' && gj == '1') {
					x11 += 1.0;
					x12 += 1.0;
				}
				else if (gi == '1' && gj == '0') {
					x11 += 1.0;
					x21 += 1.0;
				}
				else if (gi == '0' && gj == '2') {
					x12 += 2.0;
				}
				else if (gi == '2' && gj == '0') {
					x21 += 2.0;
				}
				else if (gi == '1' && gj == '2') {
					x12 += 1.0;
					x22 += 1.0;
				}
				else if (gi == '2' && gj == '1') {
					x21 += 1.0;
					x22 += 1.0;
				}
				else {
					fprintf(stderr, "error! exiting...\n");
					exit(1);
				}
			} // close k
			x11 /= ((double)(num_indv)*2.0);
			x12 /= ((double)(num_indv)*2.0);
			x21 /= ((double)(num_indv)*2.0);
			x22 /= ((double)(num_indv)*2.0);
			double p1 = x11 + x12;
			double p2 = x21 + x22;
			double q1 = x11 + x21;
			double q2 = x12 + x22;
			double D = x11*x22 - x12*x21;
			double r = D/(sqrt(p1*p2*q1*q2));
			fprintf(fout, "%.8lf ", r);
		} // close j
		fprintf(fout, "\n");
	} // close i
	fclose(fout);
	
	// clean up
	free(IN_TYPED_SNP_FILE);
	free(IN_GEN_FILE);
	free(OUT_LD_FILE);
}
