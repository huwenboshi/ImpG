#include <vector>
#include "util.h"

using namespace std;

// load genotype file, return the number of individuals
size_t load_genotypes(const char* filename,
			size_t num_typed_snps, vector<string> &gens);

// get the command line input for gen_ld
int get_gen_ld_cmd_line(int argc, char **argv, char **IN_TYPED_SNP,
			char **IN_GEN_FILE, char **OUT_LD_FILE);
