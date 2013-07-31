#include <string.h>

#include "genld.h"

// load genotype file, return the number of snps
size_t load_genotypes(const char* filename,
			size_t num_typed_snps, vector<string> &gens) {
	FILE *fin = safe_fopen(filename, "r");
	
	// skip first line
	skip_first_line(fin);
	
	// read in the rest
	size_t nbytes = 1024;
	char *buf = (char*) malloc(nbytes*sizeof(char));
	ssize_t num_read;
	size_t line = 2;
	while((num_read = getline(&buf, &nbytes, fin)) > 0) {
		buf[num_read-1] = '\0';
		string gen = string(buf);
		gens.push_back(gen);
		if(num_read != (ssize_t)num_typed_snps+1) {
			fprintf(stderr, "Error: Expecting %u of SNPs at line %u, ",
				(unsigned int)num_typed_snps, (unsigned int)line);
			fprintf(stderr, "instead of %u SNPs...\n",
				(unsigned int)(num_read-1));
			exit(1);
		}
		line++;
	}
	free(buf);
	fclose(fin);	
	
	size_t ninds = gens.size();
	if(verbose) {
		printf("Info: Loaded genotypes of %u individuals from %s...\n", 
			(unsigned int)ninds, filename);
	}
	
	return ninds;
}

// get the command line input for gen_ld
int get_gen_ld_cmd_line(int argc, char **argv, char **IN_TYPED_SNP,
			char **IN_GEN_FILE, char **OUT_LD_FILE) {
	int nflags = 0;
	char opt;
	while((opt = getopt(argc,argv,"a:t:g:o:")) != -1) {
		switch(opt) {
			case 't': // typed snp file
				*IN_TYPED_SNP = (char*) strdup(optarg);
				nflags++;
				break;
			case 'g': // genotype file
				*IN_GEN_FILE = (char *) strdup(optarg);
				nflags++;
				break;
			case 'o': // output file
				*OUT_LD_FILE = (char *) strdup(optarg);
				nflags++;
				break;
		}
	}
	if(!(*IN_TYPED_SNP) || !(*IN_GEN_FILE) || !(*OUT_LD_FILE)) {
		fprintf(stderr, "Usage:\n");
        fprintf(stderr, "\t-t (required) specify typed SNPs file\n");
		fprintf(stderr, "\t-g (required) specify genotype file\n");
		fprintf(stderr, "\t-o (required) specify output file name\n");
		exit(1);
	}
	
	return nflags;
}
