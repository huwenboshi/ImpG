#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>


#include "impg.h"
#include "linsubs.h"

using namespace std;

int main(int argc, char **argv) {
	double maf_th = 0.01;
    
    string filename, prefix;

	// get command line input
	char *IN_HAP_FILE = NULL, *IN_ALL_SNP_FILE = NULL;
	char *IN_TYPED_SNP_FILE = NULL, *OUT_FILE_PREFIX = NULL;
    char *MAF_TH = NULL;
	get_gen_beta_cmd_line(argc, argv, &IN_HAP_FILE, &IN_ALL_SNP_FILE,
		&IN_TYPED_SNP_FILE, &OUT_FILE_PREFIX, &MAF_TH);
	if(MAF_TH != NULL) {
        maf_th = atof(MAF_TH);
    }
    if(verbose) {
        printf("Info: Using SNPs with MAF >= %.8lf for constructing sigma matrix...\n", maf_th);
    }

    prefix = string(OUT_FILE_PREFIX);


	// load typed snps
	vector<typed_snp> typed_snps;
	size_t num_typed_snps = load_typed_snps(IN_TYPED_SNP_FILE, typed_snps);
	
	// load all snps
	vector<ref_snp> all_snps;
	size_t num_total_snps = load_all_snps(IN_ALL_SNP_FILE, all_snps);

	// mark snps whose z-scores need to be converted in the imputation step
	vector<char> convert_flags;
	convert_flags.resize(num_typed_snps, 0);
	vector<char> impute_flags;
	impute_flags.resize(num_total_snps, 1);
	mark_snps(typed_snps, all_snps, convert_flags, impute_flags);

	// load haplotypes
	vector<string> haps;
	load_haplotypes(IN_HAP_FILE, haps, num_total_snps);
	
	// compute allele frequencies for all snps
	vector<double> freqs;
	freqs.resize(num_total_snps, 0.0);
	get_all_freqs(haps, freqs);
	
	// estimate sigma_t matrix using maf filtered typed snps
	size_t sigma_t_tmp_size = num_typed_snps*num_typed_snps;
	double *sigma_t_tmp = (double*)safe_calloc(sigma_t_tmp_size,
		sizeof(double));
	
	for(size_t i = 0; i < num_typed_snps; i++) {
		size_t idxi = typed_snps[i].idx;
		double pi = freqs[idxi];
		for(size_t j = i+1; j < num_typed_snps; j++) {
			size_t idxj = typed_snps[j].idx;
			double pj = freqs[idxj];
			double pij = get_h_freq(haps, idxi, idxj);
			double r = (pij-pi*pj)/sqrt(pi*(1.0-pi)*pj*(1.0-pj));  
			sigma_t_tmp[i*num_typed_snps+j] = r;
			sigma_t_tmp[j*num_typed_snps+i] = sigma_t_tmp[i*num_typed_snps+j];
		}
		sigma_t_tmp[i*num_typed_snps+i] = 1.0;
	}
	
	// output the betas and vars
	get_beta_var(prefix, haps, all_snps, typed_snps,
		freqs, impute_flags, sigma_t_tmp, LAMBDA, maf_th);

	// clean up
	free(sigma_t_tmp);	
	free(IN_HAP_FILE);
	free(IN_ALL_SNP_FILE);
	free(IN_TYPED_SNP_FILE);
	free(OUT_FILE_PREFIX);
    if(MAF_TH != NULL) {
        free(MAF_TH);
    }
}
