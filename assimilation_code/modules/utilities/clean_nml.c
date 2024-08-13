/* 
 * DART software - Copyright UCAR. This open source software is provided
 * by UCAR, "as is", without charge, subject to all terms of use at
 * http://www.image.ucar.edu/DAReS/DART/DART_download
 */
 

/*
 * read in a fortran namelist on stdin, reformat and sort, then
 * output on stdout.  make spacing, commas, etc consistent.
 * the goal is to make it easier to compare 2 different namelist
 * files which have diverged in order of namelists inside the file,
 * formatting, etc.
 *
 * usage: [ -no_sort_nmls ] [ -no_sort_entries ] [ -no_remove_comments]  < stdin > stdout
 * 
 * the default is to sort all namelists by name, then sort
 * the contents of each list, and remove all comments lines
 * (those that appear outside any namelists).
 *
 * each of the arguments changes those defaults.
 * 
 * nsc 18nov2011 
 *     13feb2017
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
 *
 * i have not yet put in support yet for:
 *  &namelist name=value / 
 * (start, val, stop all on the same line) which is technically legal
 * but in practice is hardly ever used..
 *
 * other things this could do:
 *
 * - find and handle duplicate namelists; either complain, compare
 * and output a single copy of duplicates if they have identical
 * contents, or output them all.
 *
 * - give an option to sort the values inside a namelist array item, 
 * for use in comparing two different namelists to see if the values
 * are the same or not.  must be careful because this breaks things
 * if the order of the values in the item are significant.
 *
 * - an option to preserve comments near the namelist they were
 * read in from.  one wrinkle - it's possible to have comments
 * before or after a namelist; there's no unambiguous way to
 * know which namelist they should stay associated with.
 * 
 * - rearrange the code so the array limits reallocate themselves
 * if they're too small.  right now they are a fixed size
 * allocated at the start of run time.
 *
 */

#define MAXNMLS      400
#define MAXENTRIES  4000
#define MAXVALUES   9000
#define MAXCOMMENTS 5000

/* data structs */

struct nv_pairs {    /* name-value pairs */
   char *name;
   int nvalues;
   char **value;
};

struct nml {         /* namelist */
   char *name;
   int nitems;
   struct nv_pairs nvp[MAXENTRIES];
    int sort_list[MAXNMLS];
};

struct nmllist {    /* list of namelists */
  struct nml nmll[MAXNMLS];
  int nmllcount;
  char *comment[MAXCOMMENTS];
  int commentcount;
  int sort_list[MAXNMLS];
};

struct nmllist l;

/* make these long */
#define MAXLINE 2048
#define MAXTOKEN 1024

char linebuf[MAXLINE];
char nbuf[MAXTOKEN];
char vbuf[MAXTOKEN];

/* subroutine declarations */
void setup(void);
void takedown(void);
void readin(void);
void do_nml_sort(void);
void writeout(int sortnml, int sortitems, int suppresscomments);
void printnml(struct nml *nl, int sorted);
int nmlstart(char *line, int linelen, char **name);
int nmlend(char *line, int linelen);
int onlyslash(char *line, int linelen);
int emptyline(char *line, int linelen);
void sortmyitems(struct nml *nl);
int samecontents(struct nml *nl1, struct nml *nl2);
int splitme(char *line, int linelen, char **name, char **value);
char *haschar(char *line, int linelen, char target, int nolead);
int longestname(struct nml *nl);
int nextname(char *line, int linelen, int offset, int *start, int *end);
int nextvalue(char *line, int linelen, int offset, int *start, int *end);
int justvalue(char *line, int linelen, char **value);

int main(int argc, char **argv)
{
    int sort_nmls = 1;
    int sort_entries = 1;
    int remove_comments = 1;

 /* usage: [ -no_sort_nmls ] [ -no_sort_entries ] [ -no_remove_comments] */

    while (argc > 1) {
      if (!strcmp(argv[1], "-no_sort_nmls")) 
        sort_nmls = 0;
      else if (!strcmp(argv[1], "-no_sort_entries")) 
        sort_entries = 0;
      else if (!strcmp(argv[1], "-no_remove_comments")) 
        remove_comments = 0;
      else {
        fprintf(stderr, 
          "usage: %s [ -no_sort_nmls ] [ -no_sort_entries ] [ -no_remove_comments] < stdin > stdout\n", argv[0]);
        fprintf(stderr, 
          "  defaults are to sort all namelists, then sort items inside each namelist, and to delete\n");
        fprintf(stderr, 
          "  all lines outside of a namelist.  use the arguments to change these defaults.\n");
       exit (-1);
    }
      --argc; argv++;
    }

    setup();

    /* read stdin, add new namelists for each & encountered.
     * read contents, adding a new item for each name=value pair.
     * sort
     * output
     */

    readin();

    do_nml_sort();
 
    /* set first arg to 0 to avoid alphabetical sort, set second arg to 0
     * to avoid sorting the contents of each namelist, last arg to 0 to
     * remove comments instead of appending them to the end of the output.
     */
    writeout(sort_nmls, sort_entries, remove_comments);

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
/* printf("before line: '%s' \n", linebuf); */
/* printf("before linelen = %d\n", linelen); */
/* printf("before char[n] = '%c' \n", linebuf[linelen-1]); */
        if (linebuf[linelen-1] == '\n') {
            linebuf[linelen-1] = '\0';
            linelen--;
        }
/* printf("after line: '%s' \n", linebuf); */
/* printf("after linelen = %d\n", linelen); */
/* printf("after char[n] = '%c' \n", linebuf[linelen-1]); */
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
                    (l.comment[l.commentcount-1])[linelen] = '\0';
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

void do_nml_sort()
{
    int i, j, tmp;

    for (i=0; i<l.nmllcount; i++) {
        l.sort_list[i] = i;
        sortmyitems(l.nmll+i);
    }

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


void writeout(int sortnml, int sortitems, int suppresscomments)
{
    int i, next, nextp1;

    /* all nmls first */
    for (i=0; i<l.nmllcount; i++) {
        
        if (!sortnml)
            /* you would have to search to find dups if you don't sort,
             * so for now don't even try if you're not sorting.
             */
            printnml(l.nmll+i, sortitems);

        else {

            next = l.sort_list[i];

            if (i < l.nmllcount - 1) {
                nextp1 = l.sort_list[i+1];

                if (strcmp(l.nmll[next].name, l.nmll[nextp1].name) == 0) {
                    /* compare to be sure contents are identical */
                    if (samecontents(l.nmll+next, l.nmll+nextp1)) {
                        continue;
                    } else 
                        fprintf(stderr, "warning: duplicate of namelist %s found but contents are\n", l.nmll[next].name);
                        fprintf(stderr, " NOT identical.  both namelists are being written to output.\n");
    }
            }

            printnml(l.nmll+next, sortitems);
        }
    }
  
    if (!suppresscomments) {
        /* all comments at end of file. */
        /* should be option to attach them to the nml
         * they are closest to (hard because the comments
         * can before or after the nml.)
         */
    for (i=0; i<l.commentcount; i++) 
        printf("%s\n", l.comment[i]);
    }

    printf("\n\n");
}

/* lcase the left name, lcase .true. and .false.? */
void printnml(struct nml *nl, int sorted)
{
    int i, j, next, len;
    char formatE[32], formatEc[32], formatS[32], formatSc[32];

    printf("&%s\n", nl->name);
    len = longestname(nl);
    sprintf(formatE,  "  %%-%ds = %%s\n",  len);
    sprintf(formatEc, "  %%-%ds = %%s,\n", len);
    sprintf(formatS,  "  %%-%ds   %%s\n",  len);
    sprintf(formatSc, "  %%-%ds   %%s,\n", len);

    /* call longestname() here and set name format len */
    for (i=0; i<nl->nitems; i++) {
        /* sort contents or not? */
        if (sorted) 
            next = nl->sort_list[i];
        else
            next = i;

        if (nl->nvp[next].nvalues > 1) 
            printf(formatEc, nl->nvp[next].name, nl->nvp[next].value[0]);
        else
            printf(formatE,  nl->nvp[next].name, nl->nvp[next].value[0]);

        if (nl->nvp[next].nvalues > 1) {
            for (j=1; j<nl->nvp[next].nvalues-1; j++) 
                printf(formatSc, "", nl->nvp[next].value[j]);
            printf(formatS, "", nl->nvp[next].value[j]);
        }
    }
    printf("/\n");
    printf("\n\n");
}

/* make sure the & isn't in quotes; stop the name at the next whitespace.
 * there shouldn't be anything but whitespace before the &, and there
 * must be a char immediately after the &.   e.g. this was being flagged
 * as the start of a namelist but it clearly shouldn't have been:
 * # grid interpolation routines will wrap over the north & south poles.
 */
int nmlstart(char *line, int linelen, char **name) 
{
    int i, len;
    char *e, c;

    /* find the location of the (first) & in the line */
    e = haschar(line, linelen, '&', 1);
    /* is there an & in the line?  if so, let's parse further */
    if (e != NULL) {
        /* len is going to eventually be the length of the 
         * namelist name.   the longest it can be is the entire
         * rest of the line, minus the trailing null (or newline?)
         */
        len = linelen - (e-line) - 1;
        for (i=(e-line)+1; i<linelen; i++) {
            c = line[i];
/* printf("nmlstart: i %d, c '%c', (e-line) %ld, linelen %d\n", i, c, (e-line), linelen);*/
/* printf("line: '%s'\n", line);*/
            if (isspace(c)) {
                len = i - (e-line) - 1;
                break;
            }
        }
/* printf("len now %d\n", len);*/

        *name = malloc(len + 1);
        strncpy(*name, e+1, len);
        /* lowercase name */
        for (i=0; i<len; i++)
            (*name)[i] = (char)tolower((int)(*name)[i]);
        (*name)[len] = '\0';
        return 1;
    } else 
        /* no &, return that this is not the start of a namelist */
        return 0;
}

/* ok, these need to get smarter.  if there are quotes, either single or
 * double, slashes don't count.
 */
int nmlend(char *line, int linelen)
{
    int i, len;
    char *e;

    e = haschar(line, linelen, '/', 0);
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

/* compare two namelists and report if their contents are the same. */
/* assumes the contents have already been sorted */
int samecontents(struct nml *nl1, struct nml *nl2)
{
    int i, j, next1, next2;

    /* start with the fast compares here */
    if (nl1->nitems != nl2->nitems) return 0;
    if (strcmp(nl1->name, nl2->name) != 0) return 0;

    for (i=0; i<nl1->nitems; i++) {
        next1 = nl1->sort_list[i];
        next2 = nl2->sort_list[i];
        if (strcmp(nl1->nvp[next1].name, nl2->nvp[next2].name) != 0) return 0;
        if (nl1->nvp[next1].nvalues != nl2->nvp[next2].nvalues) return 0;
        for (j=0; j<nl1->nvp[next1].nvalues; j++) 
            if (strcmp(nl1->nvp[next1].value[j], nl2->nvp[next2].value[j]) != 0) return 0;
    }
    return 1;
}

/* sort the interior contents of a namelist so all items are listed
 * in alphabetical order.  could be useful to compare two namelists
 * which have drifted apart.
 */
void sortmyitems(struct nml *nl)
{
    int i, j, tmp;

    for (i=0; i<nl->nitems; i++) 
        nl->sort_list[i] = i;

    /* here's where we'd print out/flag duplicate items in the same list */
    for (i=0; i<nl->nitems; i++) {
        for (j=0; j<nl->nitems-1; j++) {
            if (strcmp(nl->nvp[nl->sort_list[j]].name, nl->nvp[nl->sort_list[j+1]].name) > 0) {
                tmp = nl->sort_list[j];
                nl->sort_list[j] = nl->sort_list[j+1];
                nl->sort_list[j+1] = tmp;
            }
        } 
    }
}


/* stop at commas (outside of quotes), stop name at whitespace */
int splitme(char *line, int linelen, char **name, char **value)
{
    int i, len, startc, endc;
    char *e;

    e = haschar(line, linelen, '=', 0);
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
        (*name)[len] = '\0';
    } else {
        *name = NULL;
    }

    len = (e - line) + 1;
   
    if (nextvalue(line, linelen, len, &startc, &endc)) {
        len = endc - startc + 1;
        *value = malloc(len + 1);
        strncpy(*value, line+startc, len);
        (*value)[len] = '\0';
/* printf("nextvalue returns '%s'\n", *value); */
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
/* printf("in value\n"); */
            in_value = 1;
            if (*end < 0)
                *end = i;
        } else {
            if (!isspace(c)) continue;
/* printf("done value\n"); */
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
 * single or double quotes.  if nolead is 1, there can't be any
 * leading non-whitespace chars before the target char.  if nolead
 * is 0, it's ok to have intervening non-whitespace chars, but
 * the rules about outside of quotes still holds.
 */
char *haschar(char *line, int linelen, char target, int nolead)
{
    int i;
    int in_squote, in_dquote, leading_white;

    in_squote = 0;
    in_dquote = 0;
    leading_white = 1;

    for (i=0; i<=linelen; i++) {
        if (in_squote) {
            if (line[i] != '\'') continue;
            in_squote = 0;
            continue;
        }
        if (line[i] == '\'') {
            in_squote = 1; 
            leading_white = 0;
            continue;
        }
        if (in_dquote) {
            if (line[i] != '"') continue;
            in_dquote = 0;
            continue;
        }
        if (line[i] == '"') {
            in_dquote = 1; 
            leading_white = 0;
            continue;
        }
        if (line[i] == target) {
            /* we found the target char, but if there were
             * intervening non-whitespace chars between the start
             * of the line and this char, and the user set the 'nolead'
             * flag, then don't say we found the char.
             */
            if (nolead && !leading_white)
                return NULL;

            return line+i;
        }
        if (!isspace(line[i])) 
            leading_white = 0;
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
        (*value)[len] = '\0';
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
