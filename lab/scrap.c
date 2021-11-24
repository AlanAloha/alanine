#include <stdio.h>
#include <string.h>
#include "genomikon.h"

char *reverse_string(char *str) {
	int len = strlen(str);
	char *rev = malloc(len);
	
	for(int i = 0; i < len; i++) {
		rev[i] = str[len-i-1];
	}
	
	return rev;
}

struct Sequence {
	int len;
	char *seq;
};

typedef struct Sequence * sequence;

sequence new_sequence(const char *seq) {
	sequence my_seq = malloc(sizeof(struct Sequence));
	my_seq->seq = malloc(strlen(seq)+1);
	strcpy(my_seq->seq, seq);
	my_seq->len = strlen(seq);
	return my_seq;
	
}

int main(int argc, char **argvm) {
	gkn_tvec test = gkn_tvec_new();
	gkn_tvec_push(test, "CSD");
	gkn_tvec_push(test, "A");
	printf("%s, %d\n", test->elem[0],test->size);
	
	sequence s1 = new_sequence("AGGTCA");
	sequence s2 = new_sequence("TGCCATA");
	printf("%s %d\n", s1->seq, s1->len);
	printf("%s %d\n", s2->seq, s2->len);
	
	gkn_vec *seqs = gkn_vec_new();
	gkn_vec_push(test,seqs);
	gkn_vec_push(test,seqs);

}
