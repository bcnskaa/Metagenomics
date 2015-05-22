#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "ctype.h"



char** map_gi2ctx;

typedef struct {
	char **gi;
	int *line;
	size_t gi_n;
} gi2line_map;

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
char* get_col_val(str, tok_pos, del)
	char *str;		// String to be processed
	int tok_pos; 	// Token position to be extracted
	char *del; 		// Delimiting character
{
	char *ptr;
	int tok_n;

	tok_n = 0
	ptr = str;
	while(ptr != 0)
	{
		if(ptr[0] == del)
		{
			tok_n++;
		}
	}

}

//char* get_col_val(str, tok_pos, del)
//	char *str;		// String to be processed
//	int tok_pos; 	// Token position to be extracted
//	char *del; 		// Delimiting character
//{
//	char *val;
//	char *tok = strtok(str, del);
//	int cur_idx;
//
//	val = 0;
//	cur_idx = 1;
//
//
//	while(tok)
//	{
//		printf("%d: %s\n", cur_idx, tok);
//		if(tok_pos == cur_idx)
//		{
//			val = strdup(tok);
//			printf("===%d:%s\n", tok_pos, val);
//			break;
//		}
//		tok = strtok(0, del);
//		cur_idx++;
//	}
//
//	return val;
//}


int separate_ids(str, vals, del)
	char *str;
	char ***vals;
	char *del;
{
	int val_count;
	char *tok = strtok(str, del);
	val_count = 0;

	while(tok)
	{
		*vals = (char**)realloc(*vals, (val_count + 1) * sizeof(char*));
		(*vals)[val_count] = strdup(tok);

		val_count++;

		tok = strtok(0, del);
	}

	return val_count;

}


/**
 *
 */
char** import_dat(fn, idx_col, dat_col, line_count)
	char *fn;
	int idx_col;
	int dat_col;
	int line_count;
{
	char** dat, **val;
	FILE *fp;
	char ch;
	char *buf, *tok, *tok_val;
	size_t read_char_n, buf_size, tok_idx, cur_line, val_count;
	int i;

	dat = (char**)calloc(line_count, sizeof(char*));

	buf = 0;
	read_char_n = 0;
	cur_line = 0;
	buf_size = 0;
	fp = fopen(fn, "r");

	// Read in the file
	while((read_char_n = getline(&buf, &buf_size, fp)) != -1)
	{
		//printf("Line %d: %s\n", cur_line, buf);
		tok_val = get_col_val(buf, dat_col, "\t");

		val = 0;

		//printf("Line %d: %s\n", cur_line, buf);
		val_count = separate_ids(tok_val, &val, ";");

		printf("Line %d, ID_n=%d (%s)\n", cur_line, val_count, val[0]);

		dat[cur_line] = strdup(val[0]);

		for(i = 0; i < val_count; i++)
			free(val[i]);
		free(val); val = 0;

		free(tok_val); tok_val = 0;

		cur_line++;
	}
	fclose(fp); fp = 0;
	free(buf); buf = 0;

	return dat;
}


/**
 *
 */
int main(argc, argv)
  int argc;
  char **argv;
{
	  char *data_fn;
	  char **dat;
	  int l_n, i;
	  int dat_col, idx_col;

	  if(argc != 4)
		  return 0;

	  data_fn = argv[1];
	  idx_col = atoi(argv[2]);
	  dat_col = atoi(argv[3]);


	  printf("Reading from %s (idx_col=%d, dat_col=%d)\n", data_fn, idx_col, dat_col);

	  l_n = count_line(data_fn);
	  printf("Number of line: %d\n", l_n);

	  printf("Importing data...\n");
	  dat = import_dat(data_fn, idx_col, dat_col, l_n);

	  for(i = 0; i < l_n;i++)
		  free(dat[i]);
	  free(dat); dat = 0;

	  return 0;
}

