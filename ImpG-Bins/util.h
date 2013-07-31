#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string>
#include <vector>

using namespace std;

#ifndef UTILS
#define UTILS

const int verbose = 1;

typedef long long int pos_t;

// provided by user
typedef struct {
	size_t idx;
	string snp_name;
	pos_t snp_pos;
	char ref_allele;
	char alt_allele;
	double zscore;
} typed_snp;

// from internal file
typedef struct {
	size_t idx;
	string snp_name;
	pos_t snp_pos;
	char ref_allele;
	char alt_allele;
} ref_snp;

// from internal file
typedef struct {
	size_t idx;
	string snp_name;
	pos_t snp_pos;
	double zscore;
} zscore_typed_snp;


// to store a genotype
typedef struct {
	size_t idx;
	string snp_name;
	pos_t snp_pos;
	string zero;
	string one;
	string two;
} genotype;


FILE* safe_fopen(const char* filename, const char *op);

void* safe_calloc(size_t nelements, size_t sizeof_element);

void skip_first_line(FILE *fin);

// load typed snps file
size_t load_typed_snps(const char *filename, 
					vector<typed_snp> &typed_snps);
					
// load all snps file
size_t load_all_snps(const char *filename, vector<ref_snp> &snp_ref_alt);

// search for the index of the snp whose pos is specified by pos
// assume the snps in typed_snps are sorted by snp positions
bool search_by_pos(const vector<typed_snp> &typed_snps,
			pos_t pos, size_t &idx);

// search for the index of the snp whose pos is specified by pos
// assume the snps in all_snps are sorted by snp positions
bool search_by_pos(const vector<ref_snp> &all_snps, pos_t pos, size_t &idx);

// search for the index of the snp whose pos is specified by pos
// assume the snps in all_snps are sorted by snp positions
bool search_by_pos(const vector<genotype> &all_genos, pos_t pos, size_t &idx);

// mark snps that are typed in snp_flags to 1
// mark snps that need to convert z-scores in convert_flags to 1
size_t mark_snps(vector<typed_snp> &typed_snps, 
			const vector<ref_snp> &snp_ref_alt,vector<char>& convert_flags,
			vector<char>& impute_flags);

// load haplotypes from the 1000 genome project reference panel file
size_t load_haplotypes(const char *filename,
					vector<string> &haps, size_t num_total_snps);

#endif
