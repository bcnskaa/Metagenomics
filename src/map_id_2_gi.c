#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "ctype.h"



char** map_gi2ctx;


/**
 *
 */
int count_line(fn)
 char *fn;
{
	int l_n;
	FILE *fp;
	char ch;

	l_n = 0;
	fp = fopen(fn, "r");

	while(!feof(fp))
	{
		ch = fgetc(fp);
		if(ch == '\n')
		{
			l_n++;
		}
	}
	fclose(fp); fp = 0;
	return l_n;
}





/**
 *
 */
int main(argc, argv)
  int argc;
  char **argv;
{
	  int l_n;

	  printf("Reading from %s\n", argv[1]);

	  l_n = count_line(argv[1]);

	  printf("Number of line: %d\n", l_n);

	  return 0;
}

