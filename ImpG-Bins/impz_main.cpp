#include <iostream>
#include <cmath>
#include <stdio.h>

#include "impg.h"

int main(int argc, char **argv) {

	// get command line input
	char *PREFIX = NULL;
	char *IN_ALL_SNP_FILE = NULL;
	char *IN_TYPED_SNP_FILE = NULL;
	char *OUT_FILE = NULL;
	get_imp_cmd_line(argc, argv, &PREFIX, &IN_ALL_SNP_FILE, 
		&IN_TYPED_SNP_FILE, &OUT_FILE);

	string filename, prefix;
	prefix = string(PREFIX);
	
	FILE *fin, *fout;

	// load typed snps
	vector<typed_snp> typed_snps;
	size_t num_typed_snps = load_typed_snps(IN_TYPED_SNP_FILE, typed_snps);
	
	// load all the snps
	vector<ref_snp> all_snps;
	size_t num_total_snps = load_all_snps(IN_ALL_SNP_FILE, all_snps);

	// mark snps whose z-scores need to be converted in the imputation step
	vector<char> convert_flags;
	convert_flags.resize(num_typed_snps, 0);
	vector<char> impute_flags;
	impute_flags.resize(num_total_snps, 1);
	mark_snps(typed_snps, all_snps, convert_flags, impute_flags);
    if(verbose) {
        size_t convert_cnt = 0, impute_cnt = 0;
        for(size_t i = 0; i < num_typed_snps; i++) {
            if((int)convert_flags[i] == 1) {
                convert_cnt++;
            }
        }
        printf("Info: Marked %u SNPs for Z-score conversion\n",
                (unsigned int) convert_cnt);
        for(size_t i = 0; i < num_total_snps; i++) {
            if((int)impute_flags[i] == 1) {
                impute_cnt++;
            }
        }
        printf("Info: %u out of %u SNPs need Z-score imputation\n",
                (unsigned int) impute_cnt, (unsigned int) num_total_snps);
    }

    // convert z-score
	for(size_t i = 0; i < num_typed_snps; i++) {
		if(convert_flags[i] == 1) {
			typed_snps[i].zscore *= -1.0;
		}
	}
	
	// load z-score typed snps
	vector<zscore_typed_snp> ztyped_snps;
	filename = prefix + ".snp";
	size_t num_ztyped_snps = load_zscore_typed_snps(filename.c_str(),
		typed_snps, ztyped_snps);
	
	// impute for all snps
	double exp_zscore_mean = 0.0;
	
	filename = prefix+".beta";
	fin = safe_fopen(filename.c_str(), "r");
	skip_first_line(fin);
	
	fout = safe_fopen(OUT_FILE, "w");
	fprintf(fout, "SNP_name SNP_pos Ref_Allele Alt_Allele Z-Score r2pred\n");
	
	size_t typed_snp_idx = 0;
	
	for(size_t idx = 0; idx < num_total_snps; idx++) {
	
		// get rid of snp name and snp pos in the first two columns
		char snp_name[1024];
		pos_t snp_pos;
		int rc = fscanf(fin,"%s %lld", snp_name, &snp_pos);
		if(rc <= 0) {
			fprintf(stderr, "error in zscore file! exiting...\n");
			exit(1);
		}
		
		// get beta and impute z-score
		double imp_zscore = 0.0;
		double beta;
		for(size_t i = 0; i < num_ztyped_snps; i++){
			rc = fscanf(fin, "%lf", &beta);
			if(rc <= 0) {
				fprintf(stderr, "error in zscore file! exiting...\n");
				exit(1);
			}
			imp_zscore += beta*(ztyped_snps[i].zscore - exp_zscore_mean);
		}
		
		// get variance
		double var;
		rc = fscanf(fin, "%lf", &var);
		if(rc <= 0) {
			fprintf(stderr, "error in zscore file! exiting...\n");
			exit(1);
		}
		
		// use original z-score for snps that don't need imputation
		if(impute_flags[idx] == 0) {
			imp_zscore = typed_snps[typed_snp_idx].zscore;
			var = 1.0;
			typed_snp_idx++;
		}
		
		// print to file
		fprintf(fout, "%s %lld %c %c %.6f %.6f\n",
			all_snps[idx].snp_name.c_str(), snp_pos, all_snps[idx].ref_allele,
			all_snps[idx].alt_allele, imp_zscore, var);
	}
	
	fclose(fin);
	fclose(fout);
	
	// clean up
	free(PREFIX);
	free(IN_ALL_SNP_FILE);
	free(IN_TYPED_SNP_FILE);
	free(OUT_FILE);
}
