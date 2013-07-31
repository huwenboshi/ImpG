#include <string.h>
#include "util.h"

FILE* safe_fopen(const char* filename, const char *op) {
	FILE* f = fopen(filename, op);
	if(!f) {
		fprintf(stderr,"Error: Could not open file %s!\n",filename);
		exit(1);
	}
	return f;
}

void* safe_calloc(size_t nelements, size_t sizeof_element) {
	void* ptr = calloc(nelements, sizeof_element);
	if(!ptr) {
		size_t nbytes = nelements*sizeof_element;
		fprintf(stderr, "Error: Could not allocate ");
		fprintf(stderr, "%u bytes of memory!\n", (unsigned int)nbytes);
		exit(1);
	}
	return ptr;
}

void skip_first_line(FILE *fin) {
	size_t nbytes = 1024;
	char *buf = (char*) malloc(nbytes*sizeof(char));
	ssize_t num_read = getline(&buf, &nbytes, fin);
	if(num_read <= 0) {
		fprintf(stderr, "Error: Empty file!\n");
		exit(1);
	}
	free(buf);
}

// load typed snps
size_t load_typed_snps(const char *filename, 
					vector<typed_snp> &typed_snps) {
	FILE *fin = safe_fopen(filename, "r");
	
	// skip first line
	skip_first_line(fin);
	
	// read in the rest
	char snp_name[1024];
	pos_t snp_pos;
	char ref_allele;
	char alt_allele;
	double zscore;
	while(fscanf(fin, "%s %lld %c %c %lf",
	snp_name, &snp_pos, &ref_allele, &alt_allele, &zscore)>0) {
		typed_snp tmp;
		tmp.snp_name = string(snp_name);
		tmp.snp_pos = snp_pos;
		tmp.ref_allele = ref_allele;
		tmp.alt_allele = alt_allele;
		tmp.zscore = zscore;
		typed_snps.push_back(tmp);
	}
	
	fclose(fin);
	
	size_t num_typed_snps = typed_snps.size();
	
	if(verbose) {
		printf("Info: Loaded %u typed SNPs from %s...\n",
			(unsigned int)num_typed_snps, filename);
	}
	
	return num_typed_snps;
}

// load encoding file
size_t load_all_snps(const char *filename, vector<ref_snp> &snp_ref_alt) {
	FILE *fin = safe_fopen(filename, "r");
	
	// skip first line
	skip_first_line(fin);
	
	// read in the rest
	size_t snp_idx = 0;
	char snp_name[1024];
	pos_t snp_pos;
	char ref_allele;
	char alt_allele;
	while(fscanf(fin, "%s %lld %c %c",
	snp_name, &snp_pos, &ref_allele, &alt_allele) > 0) {
		ref_snp rs;
		rs.idx = snp_idx;
		rs.snp_name = string(snp_name);
		rs.snp_pos = snp_pos;
		rs.ref_allele = ref_allele;
		rs.alt_allele = alt_allele;
		snp_ref_alt.push_back(rs);
		snp_idx++;
	}
	fclose(fin);
	
	size_t num_total_snps = snp_ref_alt.size();
	if(verbose) {
		printf("Info: Loaded %u SNPs from %s...\n", 
			(unsigned int)num_total_snps, filename);
	}
	
	return num_total_snps;
}

// search for the index of the snp whose pos is specified by pos
// assume the snps in typed_snps are sorted by snp positions
bool search_by_pos(const vector<typed_snp> &typed_snps,
			pos_t pos, size_t &idx) {
	size_t num_typed_snps = typed_snps.size();
	
	size_t left = 0;
	size_t right = num_typed_snps - 1;
	bool found = false;
	
	while(left <= right) {
		idx = left+(right-left)/2;
		if(idx > right) {
			break;
		}
		if(typed_snps[idx].snp_pos == pos) {
			found = true;
			break;
		}
		else if(typed_snps[idx].snp_pos < pos) {
			left = idx + 1;
		}
		else {
			right = idx - 1;
		}
	}
	
	return found;
}

// search for the index of the snp whose pos is specified by pos
// assume the snps in all_snps are sorted by snp positions
bool search_by_pos(const vector<ref_snp> &all_snps, pos_t pos, size_t &idx) {
	size_t num_total_snps = all_snps.size();
	
	size_t left = 0;
	size_t right = num_total_snps - 1;
	bool found = false;
	
	while(left <= right) {
		// to avoid overflow
		idx = left+(right-left)/2;
		
		if(idx > right) {
			break;
		}
		// found it
		if(all_snps[idx].snp_pos == pos) {
			found = true;
			break;
		}
		// too small
		else if(all_snps[idx].snp_pos < pos) {
			left = idx + 1;
		}
		// too large
		else {
			right = idx - 1;
		}
	}
	
	return found;
}

// search for the index of the snp whose pos is specified by pos
// assume the snps in all_snps are sorted by snp positions
bool search_by_pos(const vector<genotype> &all_genos, pos_t pos, size_t &idx) {
	size_t num_total_snps = all_genos.size();
	
	size_t left = 0;
	size_t right = num_total_snps - 1;
	bool found = false;
	
	while(left <= right) {
		// to avoid overflow
		idx = left+(right-left)/2;
		
		if(idx > right) {
			break;
		}
		// found it
		if(all_genos[idx].snp_pos == pos) {
			found = true;
			break;
		}
		// too small
		else if(all_genos[idx].snp_pos < pos) {
			left = idx + 1;
		}
		// too large
		else {
			right = idx - 1;
		}
	}
	
	return found;
}

// mark snps that are typed in snp_flags to 1
// mark snps that need to convert z-scores in convert_flags to 1
size_t mark_snps(vector<typed_snp> &typed_snps, 
			const vector<ref_snp> &snp_ref_alt,vector<char>& convert_flags,
			vector<char>& impute_flags) {
	size_t num_total_snps = impute_flags.size();
	size_t nsnps = typed_snps.size();
	size_t count = 0;
	
	for(size_t i = 0; i < nsnps; i++) {
		pos_t snp_pos = typed_snps[i].snp_pos;
		
		// search for the snp
		size_t j;
		bool found = search_by_pos(snp_ref_alt, snp_pos, j);
		
		if(!found) {
			fprintf(stderr, "Error: %s not found in encoding file!\n",
				typed_snps[i].snp_name.c_str());
			exit(1);
		}
		
		size_t idx = snp_ref_alt[j].idx;
		if(idx < num_total_snps) {
			impute_flags[idx] = 0;
		}
		else {
			fprintf(stderr, "Error! Flag index out of bound!\n");
			exit(1);
		}
		
		// remember the typed snp idx
		typed_snps[i].idx = idx;
		
		char ref = snp_ref_alt[j].ref_allele;
		char alt = snp_ref_alt[j].alt_allele;
		
		// matched
		if(ref == typed_snps[i].ref_allele 
			&& alt == typed_snps[i].alt_allele) {
			// do nothing
		}
		// switched
		else if(alt == typed_snps[i].ref_allele 
			&& ref == typed_snps[i].alt_allele) {
			convert_flags[i] = 1;
			count++;
		}
		// error
		else {
			fprintf(stderr, "Error: Undefined encoding for SNP %s!\n", 
				typed_snps[i].snp_name.c_str());
			exit(1);
		}
	}
	
	return count;
}

// load haplotypes from the 1000 genome project reference panel file
size_t load_haplotypes(const char *filename,
					vector<string> &haps, size_t num_total_snps) {
	FILE *fin = safe_fopen(filename, "r");
	
	// skip first line
	skip_first_line(fin);	
	
	// read in the rest
	size_t nbytes = 1024;
	char *buf = (char*) malloc(nbytes*sizeof(char));
	ssize_t num_read;
	while((num_read = getline(&buf, &nbytes, fin)) > 0) {
		buf[num_read-1] = '\0';
		string line = string(buf);
		size_t blank_pos = line.find(" ");
		string hap = line.substr(blank_pos+1, line.size()-blank_pos);
		haps.push_back(hap);
	}
	free(buf);	
	fclose(fin);
	
	// check if there is any error
	size_t num_snps = haps.size();
	if(num_snps != num_total_snps) {
	    fprintf(stderr, "Error: Number of SNPs in haplotype file ");
	    fprintf(stderr, "does not match the number of SNPs in SNP file\n");
	    exit(1);
	}
	
	size_t num_haps = haps[0].size();	
	if(verbose) {
		printf("Info: Loaded %u haplotypes from %s...\n",
			(unsigned int)num_haps, filename);
	}
	
	return num_haps;
}

