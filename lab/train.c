#include <stdio.h>
#include "genomikon.h"

char *gkmer(int, int);
int gidx(const char*, int);

int main(int argc, char **argv) {
	/*
	char *seq = argv[1];
	kn_pipe io = = gnk_pipe_open(file, "r");
	*/
	
	char *file = argv[1];
	gkn_pipe fh = gkn_pipe_open(file, "r");
	gkn_fasta ff = gkn_fasta_read(fh);
	
	char *seq = ff->seq;
	
	int n = 4;
	int num_kmer = pow(4,n);
	
	//initialize emission count
	double *emp = malloc(num_kmer*sizeof(double));
	for (int i = 0; i < num_kmer; i++) emp[i] = 0;
	
	//count
	int counter = 0;
	for (int i = 0; i < strlen(seq)-(n-1); i++) {
		char* lkmer = malloc(n*sizeof(char));
		strncpy(lkmer, seq+i, n);
		emp[gidx(lkmer,n)]++;
		counter++;
	}
	
	for (int i = 0; i < num_kmer/4; i++) {
		double local_sum = emp[i*4] + emp[i*4+1] + emp[i*4+2] + emp[i*4+3];
		for(int j = 0; j < 4; j++) emp[i*4+j] = emp[i*4+j]/local_sum;
	}
	
	//display emission probabilities
	for (int i = 0; i < num_kmer; i++) {
		printf("%s	%f\n", gkmer(i,n), emp[i]);
		if (i%4 == 3) printf("\n");
	}
	
	
	
	
}

char *gkmer(int idx, int k) {
	char *kmer = malloc(k);
	for (int i = 0; i < k; i++) {
		switch((idx/(int)pow(4,(k-i-1)))%4) {
			case 0: kmer[i] = 'A'; break;
			case 1: kmer[i] = 'C'; break;
			case 2: kmer[i] = 'G'; break;
			case 3: kmer[i] = 'T'; break;
			default: kmer[i] = 'n';
		}
	}
	return kmer;
}

int gidx(const char *kmer, int k) {
	int idx = 0;
	for (int i = 0; i < k; i++) {
		switch(kmer[i]) {
			case 'A': idx += 0*pow(4,k-i-1); break;
			case 'C': idx += 1*pow(4,k-i-1); break;
			case 'G': idx += 2*pow(4,k-i-1); break;
			case 'T': idx += 3*pow(4,k-i-1); break;
			default: idx = -1;
		}
	}
	return idx;
}
