#include <stdio.h>
#include <string.h>

#include "genomikon.h"

static char *usage = "\
sw - align sequences\n\n\
usage: sw <file1> <file2> [options]\n\
options:\n\
  -m <int>	match score [1]\n\
  -n <int>	mismatch score [-1]\n\
  -g <int>	gap score [-1]\n\
  -v		verbose\n";
int x = -1;
int y = -1;

/*
Output structure:
	- alignment start
	- alignment end	
	- alignment length
	- alignment score

*/



struct DPMatrix {
	int    len1;
	int    len2;
	char  *seq1;
	char  *seq2;
	int  **score;
	char **trace;
};

typedef struct DPMatrix * dpmatrix;

dpmatrix new_dpmatrix(const gkn_fasta ffa, const gkn_fasta ffb) {
	dpmatrix mat = malloc(sizeof(struct DPMatrix));
	mat->seq1 = NULL;
	mat->seq2 = NULL;
	mat->len1 = ffa->length + 1;
	mat->len2 = ffb->length + 1;
	
	mat->seq1 = malloc(strlen(ffa->seq)+1);
	strcpy(mat->seq1,ffa->seq);
	mat->seq2 = malloc(strlen(ffb->seq)+1);
	strcpy(mat->seq2,ffb->seq);
	
	mat->score = malloc(sizeof(int*)*mat->len1);
	mat->trace = malloc(sizeof(char*)*mat->len1);
	for (int i = 0; i < mat->len1; i++) {
		mat->score[i] = malloc(sizeof(int)*mat->len2);
		mat->trace[i] = malloc(sizeof(char*)*mat->len2);
	}
	return mat;
}

dpmatrix new_blank(void) {
	dpmatrix mat = malloc(sizeof(struct DPMatrix));
	mat->len1 = 0;
	mat->len2 = 0;
	mat->seq1 = NULL;
	mat->seq2 = NULL;
	mat->score = NULL;
	mat->trace = NULL;
	return mat;
}

void free_dpmatrix(dpmatrix matrix) {
	if (matrix->seq1) free(matrix->seq1);
	if (matrix->seq2) free(matrix->seq2);
	
	if (matrix->score) {
		for(int i = 0; i < matrix->len1; i++) {
			free(matrix->score[i]);
			free(matrix->trace[i]);
		}
		free(matrix->score);
		free(matrix->trace);
	}
}

void align(dpmatrix mat, int M, int N, int G) {
	
	//pre-initialization
	for (int i = 0; i < mat->len1; i++) for (int j = 0; j < mat->len2; j++) {
		mat->score[i][j] = 0;
		mat->trace[i][j] = 'N';
	}

	//fill
	int mm; // match/mismatch
	int lgap; // left gap
	int tgap; // top gap
	int max_score = 0;
	for (int i = 1; i < mat->len1; i++) {
		for (int j = 1; j < mat->len2; j++) {
			//determine match or mismatch
			char c1 = mat->seq1[i-1];
			char c2 = mat->seq2[j-1];
			
			if (c1 == c2) mm = mat->score[i-1][j-1] + M; else mm = mat->score[i-1][j-1] + N;
			lgap = mat->score[i][j-1] + G;
			tgap = mat->score[i-1][j] + G;
			
			if (lgap > 0 || tgap > 0 || mm > 0) {
				if (lgap > tgap && lgap > mm) {
					mat->score[i][j] = lgap;
					mat->trace[i][j] = 'L';
					if (lgap>max_score) max_score = lgap;
				} else if (tgap > mm) {
					mat->score[i][j] = tgap;
					mat->trace[i][j] = 'T';
					if (tgap>max_score) max_score = tgap;
				} else {
					mat->score[i][j] = mm;
					mat->trace[i][j] = 'M';
					if (mm>max_score) max_score = mm;
				}
			}	
		}
	}
	
	//display matrix
	/*
	for (int i = 0; i < num_row; i++) {
		for (int j = 0; j < num_col; j++) {
			printf("%d	", score[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	for (int i = 0; i < num_row; i++) {
		for (int j = 0; j < num_col; j++) {
			printf("%c	", trace[i][j]);
		}
		printf("\n");
	}
	*/
	
	printf("\nmax score = %d\n\n", max_score);
	
	//traceback
	for (int i = 1; i < mat->len1; i++) {
		for (int j = 1; j < mat->len2; j++) {
			x = i;
			y = j;
			if (mat->score[i][j] == max_score) {
				char *s1_aligned = malloc(mat->len1 + mat->len2);
				char *s2_aligned = malloc(mat->len1 + mat->len2);
				char *alignment = malloc(mat->len1 + mat->len2);
				int counter = 0;
				//printf("(x,y): (%d,%d)", x, y);
				while (mat->score[x][y]!=0) {
					if (mat->trace[x][y]=='L') {
						s2_aligned[counter] = mat->seq2[y-1];
						s1_aligned[counter] = '-';
						alignment[counter] = ' ';
						counter++;
						y--;
					}
					else if (mat->trace[x][y]=='T') {
						s1_aligned[counter] = mat->seq1[x-1];
						s2_aligned[counter] = '-';
						alignment[counter] = ' ';
						counter++;
						x--;
					}
					else {
						s1_aligned[counter] = mat->seq1[x-1];
						s2_aligned[counter] = mat->seq2[y-1];
						alignment[counter] = (s1_aligned[counter] == s2_aligned[counter])
							? '|'
							: ':';
						counter++;
						x--;
						y--;
					}	
				} 
				//Print alignments
				for (int k = counter-1; k >= 0; k--) printf("%c ", s1_aligned[k]);
				printf("\n");
				for (int k = counter-1; k >= 0; k--) printf("%c ", alignment[k]);
				printf("\n");
				for (int k = counter-1; k >= 0; k--) printf("%c ", s2_aligned[k]);
				printf("\n\n");

				
				//free memory
				free(s1_aligned);
				free(s2_aligned);
				free(alignment);
			}
		}
	}
}


	
int main(int argc, char **argv) {
	int M = 1;
	int N = -1;
	int G = -1;
	int verbose = 0;
	// Command Line Interface
	gkn_set_program_name(argv[0]);
	gkn_register_option("-m", 1);
	gkn_register_option("-n", 1);
	gkn_register_option("-g", 1);
	gkn_register_option("-v", 0);
	gkn_parse_options(&argc, argv);

	if (argc == 1) gkn_exit("%s", usage);

	char *file1 = argv[1];
	char *file2 = argv[2];
	if (gkn_option("-m")) M = atoi(gkn_option("-m"));
	if (gkn_option("-n")) N = atoi(gkn_option("-n"));
	if (gkn_option("-g")) G = atoi(gkn_option("-g"));
	if (gkn_option("-v")) verbose = 1;
	
	//open fasta files
	gkn_pipe fha = gkn_pipe_open(file1, "r"); //fha = file handler a
	gkn_pipe fhb = gkn_pipe_open(file2, "r"); // r = read
	
	gkn_fasta ffa = gkn_fasta_read(fha); //ffa = fasta file a
	gkn_fasta ffb = gkn_fasta_read(fhb);
	
	printf("%s %s %s %s\n", ffa->def, ffa->seq, ffb->def, ffb->seq);
	printf("match = %d, mismatch = %d gap = %d verbose = %d\n", M, N, G, verbose);
	
	
	
	dpmatrix mat = new_dpmatrix(ffa, ffb);
	align(mat,M,N,G);
	free_dpmatrix(mat);
	
	dpmatrix blank_mat = new_blank();
	free(blank_mat);
}


