#include <string.h>
#include <cmath>

#include "impg.h"

using namespace std;

// load z-score typed snps file
size_t load_zscore_typed_snps(const char *filename,
const vector<typed_snp>& typed_snps, vector<zscore_typed_snp> &ztyped_snps) {
	FILE *fin = safe_fopen(filename, "r");
	
	// skip first line
	skip_first_line(fin);
	
	// read in the rest
	size_t snp_idx = 0;
	char snp_name[1024];
	pos_t snp_pos;
	while(fscanf(fin, "%s %lld", snp_name, &snp_pos)>0) {
		zscore_typed_snp tmp;
		tmp.idx = snp_idx;
		tmp.snp_name = string(snp_name);
		tmp.snp_pos = snp_pos;
		size_t idx;
		bool found = search_by_pos(typed_snps, snp_pos, idx);
		if(!found) {
			fprintf(stderr, "Error: %s is not found!\n", snp_name);
			exit(1);
		}
		tmp.zscore = typed_snps[idx].zscore;
		ztyped_snps.push_back(tmp);
		snp_idx++;
	}
	fclose(fin);
	
	size_t num_ztyped_snps = ztyped_snps.size();
	if(verbose) {
		printf("Info: Loaded %u Z-score typed SNPs from %s...\n",
			(unsigned int)num_ztyped_snps, filename);
	}
	
	return num_ztyped_snps;
}

// load genotypes
size_t load_genotypes(const char *filename, 
	vector<genotype> &genotypes) {
	FILE *fin = safe_fopen(filename, "r");
	
	// skip first line
	skip_first_line(fin);
	
	// read in the rest
	size_t snp_idx = 0;
	char snp_name[1024];
	char zero[8];
	char one[8];
	char two[8];
	pos_t snp_pos;
	while(fscanf(fin, "%s %lld %s %s %s", snp_name, &snp_pos,
			zero, one, two)>0) {
		genotype tmp;
		tmp.idx = snp_idx;
		tmp.snp_name = string(snp_name);
		tmp.snp_pos = snp_pos;
		tmp.zero = string(zero);
		tmp.one = string(one);
		tmp.two = string(two);
		genotypes.push_back(tmp);
		snp_idx++;
	}
	fclose(fin);
	
	size_t num_genotypes = genotypes.size();
	if(verbose) {
		printf("Info: Loaded %u genotypes from %s...\n",
			(unsigned int)num_genotypes, filename);
	}
	
	return num_genotypes;
}

// load the ld matrix from file to memory
void load_ld_mat(const char *filename, size_t num_typed_snps, double *ld_mat) {
	FILE *fin = safe_fopen(filename, "r");
	
	// skip first line
	skip_first_line(fin);
	
	for(size_t i = 0; i < num_typed_snps; i++) {
		for(size_t j = i+1; j < num_typed_snps; j++) {
			double corr;
			if(fscanf(fin, "%lf", &corr)) {
				ld_mat[i*num_typed_snps+j] = corr;
				ld_mat[j*num_typed_snps+i] = corr;
			}
			else {
				fprintf(stderr, "Error: Wrong LD matrix...\n");
				exit(1);
			}
		}
		ld_mat[i*num_typed_snps+i] = 1.0;
	}
	
	if(verbose) {
		printf("Info: Loaded LD matrix from %s...\n", filename);
	}
}

// compute allele frequencies for all the snps
void get_all_freqs(const vector<string>& haps, vector<double> &freqs) {
	size_t num_total_snps = freqs.size();
	for(size_t i = 0; i < num_total_snps; i++) {
		freqs[i] = get_freq(haps, i);
	}
}

// compute the allele frequency of a snp, specified by snp_idx
double get_freq(const vector<string>& haps, size_t snp_idx) {
	double freq = 0.0;
	size_t nhaps = haps[0].size();
	size_t i;
	for(i = 0; i < nhaps; i++) {
		if(haps[snp_idx][i] == '1') {
			freq += 1.0;
		}
	}
	freq = freq/((double)nhaps);
	return freq;
}

// compute the h frequency of two snps, specified by snp_idx1, snp_idx2
double get_h_freq(const vector<string>& haps, 
			size_t snp_idx1, size_t snp_idx2) {
	double freq = 0.0;
	size_t nhaps = haps[0].size();
	size_t i;
	for(i = 0 ; i < nhaps; i++) {
		if(haps[snp_idx1][i] == '1' && haps[snp_idx2][i] == '1') {
			freq += 1.0;
		}
	}
	freq = freq/((double)nhaps);
	return freq;
}

double get_var(const vector<string>& haps,
		const vector<double>& freqs, int idx) {
	double mn = freqs[idx];
	double vr = 0.0;
	size_t nhaps = haps.size();
	for(size_t i = 0; i < nhaps; i++) {
		if(haps[i][idx] == '1') {
			vr += (1.0-mn)*(1.0-mn);
		}
		else {
			vr += (0.0-mn)*(0.0-mn);
		}
	}
	return (vr/((double)nhaps-1.0));
}

// output betas and vars
void get_beta_var(const string& prefix, const vector<string>& haps,
	const vector<ref_snp>& all_snps, const vector<typed_snp>& typed_snps,
	const vector<double>& freqs, const vector<char>& impute_flags,
	double *sigma_t_tmp, double lambda) {
	
	string filename;
	FILE *fout;
	
	size_t num_typed_snps = typed_snps.size();
	size_t num_total_snps = all_snps.size();
	
	// for filtering
	vector<char> filter_flags;
	filter_flags.resize(num_typed_snps, 0);

	// filter perfectly correlated snps
	/*
	for(size_t i = 0; i < num_typed_snps; i++) {
		for(size_t j = i+1; j < num_typed_snps; j++) {
			if(sigma_t_tmp[i*num_typed_snps+j] > 0.99) {
				filter_flags[j] = 1;
			}
		}
	}
	*/
	
	// filter maf snps
	for(size_t i = 0; i < num_typed_snps; i++) {
		double freq = freqs[typed_snps[i].idx];
		if(freq < 0.01 || freq > 0.99) {
			filter_flags[i] = 1;
		}
	}
	
	// count the number of unfiltered snps
	size_t num_snps = 0;
	for(size_t i = 0; i < num_typed_snps; i++) {
		if(filter_flags[i] != 1) {
			num_snps++;
		}
	}
	
	// for debugging
	if(verbose) {
		printf("Info: Using %u SNPs for sigma_tt matrix...\n",
			(unsigned int)num_snps);
	}
	
	size_t sigma_t_size = num_snps*num_snps;
	double *sigma_t = (double*)safe_calloc(sigma_t_size, sizeof(double));
	double *sigma_t_inv = (double*)safe_calloc(sigma_t_size, sizeof(double));
	double *sigma_it = (double*)safe_calloc(num_snps, sizeof(double));
	double *beta = (double*)safe_calloc(num_snps, sizeof(double));

	// construct sigma for typed snps
	size_t ii = 0;
	for(size_t i = 0; i < num_typed_snps; i++) {
		for(size_t j = 0; j < num_typed_snps; j++) {
			if(filter_flags[i] != 1 && filter_flags[j] != 1) {
				sigma_t[ii] = sigma_t_tmp[i*num_typed_snps+j];
				if(i == j) {
					sigma_t[ii] += lambda;
				}
				ii++;
			}
		}
	}

	// compute the inverse
	pdinv(sigma_t_inv, sigma_t, num_snps);
			
	// print out typed snps used for estimating sigma
	filename = prefix+".snp";
	fout = safe_fopen(filename.c_str(),"w");
	fprintf(fout, "SNP_name SNP_pos\n");
	for(size_t i = 0; i < num_typed_snps; i++) {
		if(filter_flags[i] != 1) {
			fprintf(fout,"%s %lld\n", typed_snps[i].snp_name.c_str(),
				(long long int)typed_snps[i].snp_pos);
		}
	}
	fclose(fout);
	
	// print out beta and variance
	filename = prefix+".beta";
	fout = safe_fopen(filename.c_str(),"w");
	fprintf(fout, "SNP_name SNP_pos ");
	for(size_t i = 0; i < num_typed_snps; i++) {
		if(filter_flags[i] != 1) {
			fprintf(fout,"%s ", typed_snps[i].snp_name.c_str());
		}
	}
	fprintf(fout, "r2pred\n");
	
	for(size_t idx = 0; idx < num_total_snps; idx++) {
		// for SNPs that need imputation
		if(impute_flags[idx] == 1) {
			// get sigma_it
			ii = 0;
			double pi = freqs[idx];
			for(size_t i = 0; i < num_typed_snps; i++) {
				if(filter_flags[i] != 1) {
					size_t sidx = typed_snps[i].idx;
					double pj = freqs[sidx];
					double r = 0.0;
					if(pi < 0.01 || pi > 0.99) {
						// do nothing
					}
					else {
						double pij = get_h_freq(haps, idx, sidx);
						r = (pij-pi*pj)/sqrt(pi*(1.0-pi)*pj*(1.0-pj));
					}
					sigma_it[ii] = r;
					ii++;
				}
			}
		}
		// for SNPs that don't need imputation
		else {
			// do nothing
		}
		
		// print beta
		fprintf(fout,"%s %lld", all_snps[idx].snp_name.c_str(),
			all_snps[idx].snp_pos);
		
		// for SNPs that need imputation
		if(impute_flags[idx] == 1) {
			for(size_t i = 0; i < num_snps; i++) {
				beta[i] = 0.0;
				for(size_t j = 0; j < num_snps; j++) {
					// valid because sigma_t_inv is symmetric
					beta[i] += sigma_it[j]*sigma_t_inv[num_snps*i+j];
				}
				fprintf(fout," %.8f",beta[i]);
			}
		}
		// for SNPs that don't need imputation
		else {
			for(size_t i = 0; i < num_snps; i++) {
				if(typed_snps[i].idx == idx) {
					fprintf(fout," %.8f", 1.0);
				}
				else {
					fprintf(fout," %.8f", 0.0);
				}
			}
		}
		
		// print variance
		// for SNPs that need imputation
		if(impute_flags[idx] == 1) {
			double var = 0.0;
			for(size_t i = 0; i < num_snps; i++) {
				for(size_t j = 0; j < num_snps; j++) {
					// valid because sigma_t is symmetric			
					var += beta[i]*beta[j]*sigma_t[num_snps*i+j];
				}
			}
			fprintf(fout," %.8f\n", var);
		}
		// for SNPs that don't need imputation
		else {
			fprintf(fout," %.8f\n", 1.0);
		}
	}
	fclose(fout);
	
	free(sigma_t);
	free(sigma_t_inv);
	free(sigma_it);
	free(beta);
}

// output betas and vars
void get_beta_var_raw(const string& prefix, const vector<string>& haps,
	const vector<genotype>& all_snps, const vector<genotype>& typed_snps,
	const vector<double>& freqs, double *sigma_t_tmp, double lambda) {
	
	string filename;
	FILE *fout;
	
	size_t num_typed_snps = typed_snps.size();
	size_t num_total_snps = all_snps.size();
	
	// for filtering
	vector<char> filter_flags;
	filter_flags.resize(num_typed_snps, 0);
	
	// filter perfectly correlated snps
	for(size_t i = 0; i < num_typed_snps; i++) {
		for(size_t j = i+1; j < num_typed_snps; j++) {
			if(sigma_t_tmp[i*num_typed_snps+j] > 0.99) {
				filter_flags[j] = 1;
			}
		}
	}
	
	// filter maf snps
	for(size_t i = 0; i < num_typed_snps; i++) {
		double freq = freqs[typed_snps[i].idx];
		if(freq < 0.01 || freq > 0.99) {
			filter_flags[i] = 1;
		}
	}
	
	// count the number of unfiltered snps
	size_t num_snps = 0;
	for(size_t i = 0; i < num_typed_snps; i++) {
		if(filter_flags[i] != 1) {
			num_snps++;
		}
	}
	
	// for debugging
	if(verbose) {
		printf("Info: Using %u SNPs for sigma_tt matrix...\n",
			(unsigned int)num_snps);
	}
	
	size_t sigma_t_size = num_snps*num_snps;
	double *sigma_t = (double*)safe_calloc(sigma_t_size, sizeof(double));
	double *sigma_t_inv = (double*)safe_calloc(sigma_t_size, sizeof(double));
	double *sigma_it = (double*)safe_calloc(num_snps, sizeof(double));
	double *beta = (double*)safe_calloc(num_snps, sizeof(double));

	// construct sigma for typed snps
	size_t ii = 0;
	for(size_t i = 0; i < num_typed_snps; i++) {
		for(size_t j = 0; j < num_typed_snps; j++) {
			if(filter_flags[i] != 1 && filter_flags[j] != 1) {
				sigma_t[ii] = sigma_t_tmp[i*num_typed_snps+j];
				if(i == j) {
					sigma_t[ii] += lambda;
				}
				ii++;
			}
		}
	}

	// compute the inverse
	pinv_jacobi(sigma_t_inv, sigma_t, num_snps);

	// print out typed snps used for estimating sigma
	filename = prefix+".snp";
	fout = safe_fopen(filename.c_str(),"w");
	fprintf(fout, "SNP_name SNP_pos\n");
	for(size_t i = 0; i < num_typed_snps; i++) {
		if(filter_flags[i] != 1) {
			fprintf(fout,"%s %lld\n", typed_snps[i].snp_name.c_str(),
				(long long int)typed_snps[i].snp_pos);
		}
	}
	fclose(fout);
	
	// print out beta and variance
	filename = prefix+".beta";
	fout = safe_fopen(filename.c_str(),"w");
	fprintf(fout, "SNP_name SNP_pos ");
	for(size_t i = 0; i < num_typed_snps; i++) {
		if(filter_flags[i] != 1) {
			fprintf(fout,"%s ", typed_snps[i].snp_name.c_str());
		}
	}
	fprintf(fout, "Var\n");
	
	for(size_t idx = 0; idx < num_total_snps; idx++) {
		// get sigma_it
		ii = 0;
		double pi = freqs[idx];
		for(size_t i = 0; i < num_typed_snps; i++) {
			if(filter_flags[i] != 1) {
				size_t sidx = typed_snps[i].idx;
				double pj = freqs[sidx];
				double r = 0.0;
				if(pi < 0.01 || pi > 0.99) {
					// do nothing
				}
				else {
					double pij = get_h_freq(haps, idx, sidx);
					r = (pij-pi*pj)/sqrt(pi*(1.0-pi)*pj*(1.0-pj));
				}
				sigma_it[ii] = r;
				ii++;
			}
		}
		
		// print beta
		fprintf(fout,"%s %lld", all_snps[idx].snp_name.c_str(),
			all_snps[idx].snp_pos);
		for(size_t i = 0; i < num_snps; i++) {
			beta[i] = 0.0;
			for(size_t j = 0; j < num_snps; j++) {
				beta[i] += sigma_it[j]*sigma_t_inv[num_snps*j+i];
			}
			fprintf(fout," %.8f",beta[i]);
		}
		
		// print variance
		double var = 0.0;
		for(size_t i = 0; i < num_snps; i++) {
			for(size_t j = 0; j < num_snps; j++) {
				var += beta[i]*beta[j]*sigma_t[num_snps*j+i];
			}
		}
		fprintf(fout," %.8f\n", var);
	}
	fclose(fout);
	
	free(sigma_t);
	free(sigma_t_inv);
	free(sigma_it);
	free(beta);
}

// get the command line input for gen_beta program
int get_gen_beta_cmd_line(int argc, char **argv, char **IN_HAP_FILE,
		char **IN_ALL_SNP_FILE, char **IN_TYPED_SNP_FILE,
		char **OUT_FILE_PREFIX) {
	int nflags = 0;
	char opt;
	while((opt = getopt(argc,argv,"h:m:t:p:")) != -1) {
		switch(opt) {
			case 'h': // haplotype file
				*IN_HAP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 'm': // all snps file
				*IN_ALL_SNP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 't': // typed snps file
				*IN_TYPED_SNP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 'p': // prefix for output files
				*OUT_FILE_PREFIX = (char *) strdup(optarg);
				nflags++;
				break;
		}
	}
	if(!(*IN_HAP_FILE) || !(*IN_ALL_SNP_FILE) || !(*IN_TYPED_SNP_FILE) ||
       !(*OUT_FILE_PREFIX)) {
		fprintf(stderr, "Usage:\n");
        fprintf(stderr, "\t-h (required) specify haplotype file\n");
		fprintf(stderr, "\t-m (required) specify SNP mapping file\n");
		fprintf(stderr, "\t-t (required) specify typed SNP file\n");
		fprintf(stderr, "\t-p (required) specify output file prefix\n");
		exit(1);
	}
	
	return nflags;
}

// get command line input for impute z-scores
int get_imp_cmd_line(int argc, char **argv, char **PREFIX,
		char **IN_ALL_SNP_FILE, char **IN_TYPED_SNP_FILE, char **OUT_FILE) {
	int nflags = 0;
	char opt;
	while((opt = getopt(argc,argv,"p:m:t:o:")) != -1) {
		switch(opt) {
			case 'p': // prefix used in gen_beta
				*PREFIX = (char *) strdup(optarg);
				nflags++;
				break;
			case 'm': // snp mapping file
				*IN_ALL_SNP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 't': // typed snp file
				*IN_TYPED_SNP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 'o': // output file name
				*OUT_FILE = (char *) strdup(optarg);
				nflags++;
				break;
		}
	}
	if(!(*PREFIX) || !(*IN_ALL_SNP_FILE) || !(*IN_TYPED_SNP_FILE) ||
       !(*OUT_FILE)) {
		fprintf(stderr, "Usage:\n");
        fprintf(stderr, "\t-p (required) specify input file prefix\n");
		fprintf(stderr, "\t-m (required) specify SNP mapping file\n");
		fprintf(stderr, "\t-t (required) specify typed SNP file\n");
		fprintf(stderr, "\t-o (required) specify output file name\n");
		exit(1);
	}
	
	return nflags;
}

// get the command line input for gen_beta with ld program
int get_gen_beta_ld_cmd_line(int argc, char **argv, char **IN_HAP_FILE,
		char **IN_ALL_SNP_FILE, char **IN_TYPED_SNP_FILE, char **IN_LD_FILE,
		char **OUT_FILE_PREFIX) {
	int nflags = 0;
	char opt;
	while((opt = getopt(argc,argv,"h:m:t:p:l:")) != -1) {
		switch(opt) {
			case 'h': // haplotype file
				*IN_HAP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 'm': // snp mapping file
				*IN_ALL_SNP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 't': // typed snps file
				*IN_TYPED_SNP_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 'l': // name for ld
				*IN_LD_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 'p': // prefix for output files
				*OUT_FILE_PREFIX = (char *) strdup(optarg);
				nflags++;
				break;
		}
	}
	if(!(*IN_HAP_FILE) || !(*IN_ALL_SNP_FILE) || !(*IN_TYPED_SNP_FILE) ||
       !(*IN_LD_FILE) || !(*OUT_FILE_PREFIX)) {
		fprintf(stderr, "Usage:\n");
        fprintf(stderr, "\t-h (required) specify haplotype file\n");
		fprintf(stderr, "\t-m (required) specify SNP mapping file\n");
		fprintf(stderr, "\t-t (required) specify typed SNP file\n");
		fprintf(stderr, "\t-l (required) specify LD statistics file\n");
		fprintf(stderr, "\t-p (required) specify output file prefix -p\n");
		exit(1);
	}
	
	return nflags;
}
