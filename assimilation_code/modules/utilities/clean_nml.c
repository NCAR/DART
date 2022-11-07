/* 
 * DART software - Copyright UCAR. This open source software is provided
 * by UCAR, "as is", without charge, subject to all terms of use at
 * http://www.image.ucar.edu/DAReS/DART/DART_download
 *
 * DART $Id$
 */
 

/*
 * read in a fortran namelist on stdin, reformat and sort, then
 * output on stdout.  make spacing, commas, etc consistent.
 * the goal is to make it easier to compare 2 different namelist
 * files which have diverged in order of namelists inside the file,
 * formatting, etc.
 *
 * nsc 18nov2011 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* 
 * an F90 namelist starts with &name and ends with a /
 * the contents are formatted as 'name = value' pairs
 * it is common for lines to end with comma but doesn't seem
 * to be required.  strings are enclosed in single or double quotes.
 * lines outside the &name, / delimiters are ignored and often
 * used as comments.
 *
 * separators seem to be commas, in which case multiple values
 * can occur on the same line.  arrays of values can seem to
 * occur on multiple lines.   bother.
 *
 * this program needs to read in multiple namelists (so it can
 * sort them by name) and the contents, and the lines outside the
 * namelists, and output them in a consistent spacing, order, etc
 * so we can diff them easily.
 *
 * usually there are a few more prominent namelist names that
 * should appear at the top of the file.  first pass, just do
 * them all alphabetically, but later on it would be nice to
 * be able to indicate which ones should be at the top of the
 * output file.
 *
 * it figures out the longest name and lines up the ='s
 *
 * comment lines (outside any namelist) are collected and all
 * output at the end - kinda a pain since it removes them from
 * the namelist they were about.  but it does preserve them.
 * it also calls tolower() so all text is lower case.
 * i don't remember if i put in support yet for:
 *  &namelist name=value / 
 * (start, val, stop all on the same line) which is technically legal.
 */

#define MAXNMLS      200
#define MAXENTRIES  1000
#define MAXVALUES   2000
#define MAXCOMMENTS 1000

/* data structs */
struct nv_pairs {
   char *name;
   int nvalues;
   char **value;
};

struct nml {
   char *name;
   int nitems;
   struct nv_pairs nvp[MAXENTRIES];
};

struct nmllist {
  struct nml nmll[MAXNMLS];
  int nmllcount;
  char *comment[MAXCOMMENTS];
  int commentcount;
  int sort_list[MAXNMLS];
};

struct nmllist l;

/* make these long */
#define MAXLINE 1024
#define MAXTOKEN 1024

char linebuf[MAXLINE];
char nbuf[MAXTOKEN];
char vbuf[MAXTOKEN];

void setup(void);
void takedown(void);
void readin(void);
void do_sort(void);
void writeout(int sortme);
void printnml(struct nml *nl);
int nmlstart(char *line, int linelen, char **name);
int nmlend(char *line, int linelen);
int onlyslash(char *line, int linelen);
int emptyline(char *line, int linelen);
int splitme(char *line, int linelen, char **name, char **value);
char *haschar(char *line, int linelen, char target);
int longestname(struct nml *nl);
int nextname(char *line, int linelen, int offset, int *start, int *end);
int nextvalue(char *line, int linelen, int offset, int *start, int *end);

int main(int argc, char **argv)
{
    if (argc > 1) {
       fprintf(stderr, "usage: %s < stdin > stdout\n", argv[0]);
       fprintf(stderr, "    takes no arguments\n");
       exit (-1);
    }

    setup();

    /* read stdin, add new namelists for each & encountered.
     * read contents, adding a new item for each name=value pair.
     * sort
     * output
     */

    readin();

    do_sort();
 
    /* set arg to 0 to avoid alphabetical sort */
    writeout(1);

    takedown();

    exit(0);
}  


void setup(void)
{
   /* allocate here? */ 
}

void takedown(void)
{
   /* deallocate here */ 
}

void readin(void)
{
    char *name, *val;
    int in_nml, linelen, action;
    struct nv_pairs *nvp;
    struct nml *n;

    n = NULL;
    in_nml = 0;
    while (fgets(linebuf, sizeof(linebuf), stdin) != NULL) {
        linelen = strlen(linebuf);
/*printf("before line: '%s' \n", linebuf); */
/*printf("before linelen = %d\n", linelen); */
/*printf("before char[n] = '%c' \n", linebuf[linelen-1]); */
        if (linebuf[linelen-1] == '\n') {
            linebuf[linelen-1] = '\0';
            linelen--;
        }
/*printf("after line: '%s' \n", linebuf); */
/*printf("after linelen = %d\n", linelen); */
/*printf("after char[n] = '%c' \n", linebuf[linelen-1]); */
        if (linelen <= 0) continue;
        action = 0;

        if (!in_nml) {
            /* not currently in namelist definition */
            if (nmlstart(linebuf, linelen, &name) > 0) {
                l.nmllcount++;
                n = &(l.nmll[l.nmllcount-1]);
                n->name = name;
                in_nml = 1;
                action = 1;
            } else {
                /* comment line outside nmls - keep or toss? */
                /* keep and output at end of file. */
                if (!emptyline(linebuf, linelen)) {
                    l.commentcount++;
                    l.comment[l.commentcount-1] = malloc(linelen+1);
                    strncpy(l.comment[l.commentcount-1], linebuf, linelen);
                }
            }
        }

        if (in_nml) {
            /* in a namelist */
            if (n == NULL) {
                fprintf(stderr, "internal inconsistency: in_nml true, n null\n");
                exit (-1);
            }
            /* definition of item */
            if (splitme(linebuf, linelen, &name, &val) > 0) {
                n->nitems++;
                nvp = &(l.nmll[l.nmllcount-1].nvp[n->nitems-1]);
                nvp->name = name;
                nvp->value = malloc(MAXVALUES * sizeof(char *));
                nvp->value[0] = val;
                nvp->nvalues = 1;
                action = 1;
            }
            /* only a slash on a single line? */
            if (onlyslash(linebuf, linelen) > 0) {
                n = NULL;
                nvp = NULL;
                in_nml = 0;
                action = 1;
            }
            /* continuation line - append to current name */
            if (action == 0 && nvp != NULL) {
                if (justvalue(linebuf, linelen, &val)) {
                    nvp->value[nvp->nvalues] = val;
                    nvp->nvalues++;
                }
            }
            /* trailing slash on same line as value? */
            if (nmlend(linebuf, linelen) > 0) {
                n = NULL;
                nvp = NULL;
                in_nml = 0;
                action = 1;
            }
        }
    }
}

void do_sort()
{
    int i, j, tmp;

    for (i=0; i<l.nmllcount; i++) 
        l.sort_list[i] = i;

    for (i=0; i<l.nmllcount; i++) {
        for (j=0; j<l.nmllcount-1; j++) {
            if (strcmp(l.nmll[l.sort_list[j]].name, l.nmll[l.sort_list[j+1]].name) > 0) {
                tmp = l.sort_list[j];
                l.sort_list[j] = l.sort_list[j+1];
                l.sort_list[j+1] = tmp;
            }
        } 
    }
}


void writeout(int sortme)
{
    int i;

    for (i=0; i<l.nmllcount; i++) {
        if (sortme)
            printnml(&l.nmll[l.sort_list[i]]);
        else
            printnml(&l.nmll[i]);
    }
    for (i=0; i<l.commentcount; i++) 
        printf("%s\n", l.comment[i]);

    printf("\n\n");
}

/* lcase the left name, lcase .true. and .false.? */
void printnml(struct nml *nl)
{
    int i, j, len;
    char formatE[32], formatEc[32], formatS[32], formatSc[32];

    printf("&%s\n", nl->name);
    len = longestname(nl);
    sprintf(formatE,  "  %%-%ds = %%s\n",  len);
    sprintf(formatEc, "  %%-%ds = %%s,\n", len);
    sprintf(formatS,  "  %%-%ds   %%s\n",  len);
    sprintf(formatSc, "  %%-%ds   %%s,\n", len);

    /* call longestname() here and set name format len */
    for (i=0; i<nl->nitems; i++) {
        if (nl->nvp[i].nvalues > 1) 
            printf(formatEc, nl->nvp[i].name, nl->nvp[i].value[0]);
        else
            printf(formatE,  nl->nvp[i].name, nl->nvp[i].value[0]);

        if (nl->nvp[i].nvalues > 1) {
            for (j=1; j<nl->nvp[i].nvalues-1; j++) 
                printf(formatSc, "", nl->nvp[i].value[j]);
            printf(formatS, "", nl->nvp[i].value[j]);
        }
    }
    printf("/\n");
    printf("\n\n");
}

/* make sure the & isn't in quotes; stop the name at the next whitespace.
 */
int nmlstart(char *line, int linelen, char **name) 
{
    int i, len;
    char *e, c;

    e = haschar(line, linelen, '&');
    len = linelen - (e-line) - 1;
    if (e != NULL) {
        for (i=(e-line)+1; i<linelen; i++) {
            c = line[i];
/*printf("nmlstart: i %d, c '%c', (e-line) %ld, linelen %d\n", i, c, (e-line), linelen);*/
/*printf("line: '%s'\n", line);*/
            if (isspace(c)) {
                len = i - (e-line) - 1;
                break;
            }
        }
/*printf("len now %d\n", len);*/

        *name = malloc(len + 1);
        strncpy(*name, e+1, len);
        /* lowercase name */
        for (i=0; i<len; i++)
            (*name)[i] = (char)tolower((int)(*name)[i]);
        return 1;
    } else 
        return 0;
}

/* ok, these need to get smarter.  if there are quotes, either single or
 * double, slashes don't count.
 */
int nmlend(char *line, int linelen)
{
    int i, len;
    char *e;

    e = haschar(line, linelen, '/');
    if (e != NULL)
        return 1;
    else
        return 0;
}

int onlyslash(char *line, int linelen)
{
    int i, len;
    char c;

    for (i=0; i<linelen; i++) {
        c = line[i];
        if (isspace(c)) continue;

        if (c == '/') continue;
      
        return 0;
    }
        
    return 1;
}

int emptyline(char *line, int linelen)
{
    int i, len;
    char c;

    for (i=0; i<linelen; i++) {
        c = line[i];
        if (isspace(c)) continue;
      
        return 0;
    }
        
    return 1;
}

/* stop at commas (outside of quotes), stop name at whitespace */
int splitme(char *line, int linelen, char **name, char **value)
{
    int i, len, startc, endc;
    char *e;

    e = haschar(line, linelen, '=');
    if (e == NULL) {
        *name = NULL;
        *value = NULL;
        return 0;
    } 

    /* FIXME: this should start at =, work back, skip any
     * initial whitespace, then count chars until you get
     * whitespace again.  &nmlname nam=val /
     * is unusual but technically a valid namelist and
     * this code doesn't handle it right (yet).
     */
    len = (e - line) - 1;
    /* change 0 to len, search backwards */
    if (nextname(line, linelen, 0, &startc, &endc)) {
        len = endc - startc + 1;
        *name = malloc(len + 1);
        strncpy(*name, line+startc, len);
        /* lowercase name */
        for (i=0; i<len; i++)
            (*name)[i] = (char)tolower((int)(*name)[i]);
    } else {
        *name = NULL;
    }

    len = (e - line) + 1;
   
    if (nextvalue(line, linelen, len, &startc, &endc)) {
        len = endc - startc + 1;
        *value = malloc(len + 1);
        strncpy(*value, line+startc, len);
/*printf("nextvalue returns '%s'\n", *value); */
    } else {
        *value = NULL;
    }
    return 1;
}
   
int nextname(char *line, int linelen, int offset, 
             int *start, int *end)
{
    int i, j;
    char c;

    if (offset >= linelen) 
       return 0;

    for (i=offset; i<linelen; i++) {
        c = line[i];
        if (isspace(c)) continue;

        *start = i;
        for (j=i+1; j<linelen; j++) {
            c = line[j];
            if (isspace(c) || c == '=') {
                *end = j-1;
                return 1;
            }
        }
    }
    
    return 0;
}

int nextvalue(char *line, int linelen, int offset, 
              int *start, int *end)
{
    int i, j;
    int in_squote, in_dquote, in_value;
    char c;

    if (offset >= linelen) 
       return 0;

    for (i=linelen-1; i>=offset; --i) {
        c = line[i];
        if (isspace(c)) continue;
        if (c == ',') continue;
        if (c == '/') continue;
        break;
    }

    *start = -1;
    *end = -1;
    in_squote = 0;
    in_dquote = 0;
    in_value = 0;
    /* leave i alone and start there */
    for (   ; i>=offset; --i) {
/* printf("nv: i, c = %d, '%c'\n", i, line[i]);  */
        if (in_squote) {
            if (line[i] != '\'') continue;
            in_squote = 0;
            if (in_value) {
                *start = i;
                in_value = 0;
            }
            continue;
        }
        if (line[i] == '\'') {
            in_squote = 1; 
            if (!in_value) {
                if (*end < 0)
                    *end = i;
                in_value = 1;
            }
            continue;
        }
        if (in_dquote) {
            if (line[i] != '"') continue;
            in_dquote = 0;
            if (in_value) {
                *start = i;
                in_value = 0;
            }
            continue;
        }
        if (line[i] == '"') {
            in_dquote = 1; 
            if (!in_value) {
                if (*end < 0)
                    *end = i;
                in_value = 1;
            }
            continue;
        }

        c = line[i];
        if (!in_value) {
            if (isspace(c)) continue;
/*printf("in value\n"); */
            in_value = 1;
            if (*end < 0)
                *end = i;
        } else {
            if (!isspace(c)) continue;
/*printf("done value\n"); */
            in_value = 0;
            *start = i+1;
        }
    }
    if (in_value) {
        *start = offset;
        return 1;
    }
    
    if (*start > 0 && *end > 0) 
        return 1;
    else
        return 0;
}

/* return a pointer to the first occurrence of char -- outside of
 * single or double quotes.
 */
char *haschar(char *line, int linelen, char target)
{
    int i;
    int in_squote, in_dquote;

    in_squote = 0;
    in_dquote = 0;
    for (i=0; i<=linelen; i++) {
        if (in_squote) {
            if (line[i] != '\'') continue;
            in_squote = 0;
            continue;
        }
        if (line[i] == '\'') {
            in_squote = 1; 
            continue;
        }
        if (in_dquote) {
            if (line[i] != '"') continue;
            in_dquote = 0;
            continue;
        }
        if (line[i] == '"') {
            in_dquote = 1; 
            continue;
        }
        if (line[i] == target) {
            return line+i;
        }
    }

    return NULL; 
}

int justvalue(char *line, int linelen, char **value)
{
    int i, len, startc, endc;
    char *e;

    if (nextvalue(line, linelen, 0, &startc, &endc)) {
        len = endc - startc + 1;
        *value = malloc(len + 1);
        strncpy(*value, line+startc, len);
/* printf("justvalue returns '%s'\n", *value); */
    } else {
        *value = NULL;
    }
    return 1;
}

/* return the length of the longest name in an nml */
int longestname(struct nml *nl)
{
    int i, j;
    int longest = 0;

    for (i=0; i<nl->nitems; i++) {
        if (strlen(nl->nvp[i].name) > longest)
            longest = strlen(nl->nvp[i].name);
    }

    return longest;
}

/* <next few lines under version control, do not edit>
 * $URL$
 * $Revision$
 * $Date$
 */
