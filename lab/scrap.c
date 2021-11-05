#include <stdio.h>
#include "genomikon.h"

int main(int argc, char **argv) {
	gkn_tvec test = gkn_tvec_new();
	gkn_tvec_push(test, "CSD");
	gkn_tvec_push(test, "A");
	printf("%s, %d\n", test->elem[0],test->size);
}
