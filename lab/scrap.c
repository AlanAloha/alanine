#include <stdio.h>
#include <string.h>
#include "genomikon.h"

int main(int argc, char **argv) {
	gkn_tvec test = gkn_tvec_new();
	gkn_tvec_push(test, "CSD");
	gkn_tvec_push(test, "A");
	printf("%s, %d\n", test->elem[0],test->size);
	
	char* seq = malloc(20);
	seq = "ACCGT";
	printf("%s \n", seq);
	
	revseq = strrev(seq);
	printf("%s \n", revseq);
}
