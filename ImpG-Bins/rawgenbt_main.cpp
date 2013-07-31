#include <iostream>
#include <cmath>
#include <stdio.h>

#include "impg.h"
#include "linsubs.h"

using namespace std;

int main(int argc, char **argv) {
	FILE *fout;
	string filename, prefix;

	// get command line input
	char *IN_HAP_FILE = NULL, *IN_ALL_GEN_FILE = NULL;
	char *IN_TYPED_GEN_FILE = NULL, *OUT_FILE_PREFIX = NULL;
	get_gen_beta_cmd_line(argc, argv, &IN_HAP_FILE, &IN_ALL_GEN_FILE,
		&IN_TYPED_GEN_FILE, &OUT_FILE_PREFIX);
	prefix = string(OUT_FILE_PREFIX);

	// load genotype encodings for all SNPs
	vector<genotype> all_genos;
	size_t num_total_snps = load_genotypes(IN_ALL_GEN_FILE, all_genos);
	
	// load genotype encodings for typed SNPs
	vector<genotype> typed_genos;
	size_t num_typed_snps = load_genotypes(IN_TYPED_GEN_FILE, typed_genos);
	
	vector<char> convert_flags;
	convert_flags.resize(num_typed_snps, 0);
	vector<char> impute_flags;
	impute_flags.resize(num_total_snps, 1);
	
	// mark for conversion and check if need imputation
	for(size_t i = 0; i < num_typed_snps; i++) {
		size_t pos = typed_genos[i].snp_pos;
		size_t idx;
		bool found = search_by_pos(all_genos, pos, idx);
		if(found == false) {
			printf("Error: SNP not found!\n");
			exit(1);
		}
		impute_flags[idx] = 0;
		typed_genos[i].idx = idx;
		if(all_genos[idx].zero == typed_genos[i].two &&
			all_genos[idx].two == typed_genos[i].zero) {
			convert_flags[i] = 1;
			if(verbose) {
				printf("Info: Marked %s for genotype conversion...\n",
					typed_genos[i].snp_name.c_str());
			}
		}
	}
	
	// print convert flags to file
	filename = prefix+".genotype_convert_flags";
	fout = safe_fopen(filename.c_str(), "w");
	fprintf(fout, "Y/N\n");
	for(size_t i = 0; i < num_typed_snps; i++) {
		fprintf(fout, "%d\n", (int)convert_flags[i]);
	}
	
	// print impute flags to file
	filename = prefix+".genotype_impute_flags";
	fout = safe_fopen(filename.c_str(), "w");
	fprintf(fout, "Y/N\n");
	for(size_t i = 0; i < num_total_snps; i++) {
		fprintf(fout, "%d\n", (int)impute_flags[i]);
	}
	
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
		size_t idxi = typed_genos[i].idx;
		double pi = freqs[idxi];
		for(size_t j = i+1; j < num_typed_snps; j++) {
			size_t idxj = typed_genos[j].idx;
			double pj = freqs[idxj];
			double pij = get_h_freq(haps, idxi, idxj);
			double r = (pij-pi*pj)/sqrt(pi*(1.0-pi)*pj*(1.0-pj));  
			sigma_t_tmp[i*num_typed_snps+j] = r;
			sigma_t_tmp[j*num_typed_snps+i] = sigma_t_tmp[i*num_typed_snps+j];
		}
		sigma_t_tmp[i*num_typed_snps+i] = get_var(haps, freqs, i);
	}
	
	// output the betas and vars
	get_beta_var_raw(prefix, haps, all_genos, 
		typed_genos, freqs, sigma_t_tmp, 0.0);

	// clean up
	free(sigma_t_tmp);	
	free(IN_HAP_FILE);
	free(IN_ALL_GEN_FILE);
	free(IN_TYPED_GEN_FILE);
	free(OUT_FILE_PREFIX);
}
