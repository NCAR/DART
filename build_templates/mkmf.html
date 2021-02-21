<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>mkmf user's guide</title>
</head>
<body>
<!-- title using title stylespec -->
<h1>mkmf</h1>
<h2>Introduction</h2>
<p><code>mkmf</code> is a tool written in perl version 5 that
constructs a makefile from distributed source. <code>mkmf</code>
typically produces a makefile that can compile a single executable
program. But it is extensible to create a makefile for any purpose
at all.</p>
<h3>Features of mkmf</h3>
<ul>
<li>It understands dependencies in f90 (<code>module</code>s and
<code>use</code>), the fortran <code>include</code> statement, and
the cpp <code>#include</code> statement in any type of source.</li>
<li>There are no restrictions on filenames, module names, etc.</li>
<li>It supports the concept of overlays (where source is maintained
in layers of directories with a defined precedence).</li>
<li>It can keep track of changes to <code>cpp</code> flags, and
knows when to recompile affected source (i.e, files containing
<code>#ifdef</code>s that have been changed since the last
invocation).</li>
<li>It will run on any unix platform that has perl version 5
installed.</li>
<li>It is free, and released under GPL. GFDL users can copy (or,
better still, directly invoke) the file
<code>/net/vb/public/bin/mkmf</code>.</li>
</ul>
<p>It can be downloaded via <a href=
"https://github.com/NOAA-GFDL/mkmf">GitHub</a>. <code>mkmf</code>
is pronounced <i>make-make-file</i> or <i>make-m-f</i> or even
<i>McMuff</i> (Paul Kushner's suggestion).</p>
<h2>Syntax</h2>
<p>The calling syntax is:</p>
<code>mkmf [-a abspath] [-c cppdefs] [-d] [-f] [-m makefile] [-p
program] [-t template] [-v] [-w] [-x] [args]</code>
<ol>
<li><code>-a abspath</code> attaches the <code>abspath</code> at
the <i>front</i> of all <i>relative</i> paths to sourcefiles.</li>
<li><code>cppdefs</code> is a list of <code>cpp</code>
<code>#define</code>s to be passed to the source files: affected
object files will be selectively removed if there has been a change
in this state.</li>
<li><code>-d</code> is a debug flag to <code>mkmf</code> (much more
verbose than <code>-v</code>, but probably of use only if you are
modifying <code>mkmf</code> itself).</li>
<li><code>-f</code> is a formatting flag to restrict lines in the
makefile to 256 characters. This was introduced in response to a
customer who wanted to edit his makefiles using <code>vi</code>).
Lines longer than that will use continuation lines as needed.</li>
<li><code>makefile</code> is the name of the makefile written
(default <code>Makefile</code>).</li>
<li><code>template</code> is a file containing a list of make
macros or commands written to the beginning of the makefile.</li>
<li><code>program</code> is the name of the final target (default
<code>a.out</code>)</li>
<li><code>-v</code> is a verbosity flag to <code>mkmf</code></li>
<li><code>-w</code> generates compile rules which use the `wrapper'
commands MPIFC and MPILD instead of FC and LD. These can then be
defined as the mpif90 compile scripts to ease changing between an
MPI and non-MPI version.</li>
<li><code>-x</code> executes the makefile immediately.</li>
<li><code>args</code> are a list of directories and files to be
searched for targets and dependencies.</li>
</ol>
<h2>Makefile structure</h2>
<p>A <i>sourcefile</i> is any file with a source file suffix
(currently <code>.F, .F90, .c, .f. .f90</code>). An
<i>includefile</i> is any file with an include file suffix
(currently <code>.H, .fh, .h, .inc</code>). A valid sourcefile can
also be an includefile.</p>
<p>Each sourcefile in the list is presumed to produce an object
file with the same basename and a <code>.o</code> extension in the
current working directory. If more than one sourcefile in the list
would produce identically-named object files, only the first is
used and the rest are discarded. This permits the use of overlays:
if <code>dir3</code> contained the basic source code,
<code>dir2</code> contained bugfixes, and <code>dir1</code>
contained mods for a particular run, <code>mkmf dir1 dir2
dir3</code> would create a makefile for correct compilation. Please
note that precedence <i>descends</i> from left to right. This is
the conventional order used by compilers when searching for
libraries, includes, etc: left to right along the command line,
with the first match invalidating all subsequent ones. See the
<a href="#examples">Examples</a> section for a closer look at
precedence rules.</p>
<p>The makefile currently runs <code>$(FC)</code> on fortran files
and <code>$(CC)</code> on C files (unless the <code>-w</code> flag
is specified). Flags to the compiler can be set in
<code>$(FFLAGS)</code> or <code>$(CFLAGS)</code>. The final loader
step executes <code>$(LD)</code>. Flags to the loader can be set in
<code>$(LDFLAGS)</code>. Preprocessor flags are used by
<code>.F</code>, <code>.F90</code> and <code>.c</code> files, and
can be set in <code>$(CPPFLAGS)</code>. These macros have a default
meaning on most systems, and can be modified in the template file.
The predefined macros can be discovered by running <code>make
-p</code>.</p>
<p>In addition, the macro <code>$(CPPDEFS)</code> is applied to the
preprocessor. This can contain the <code>cpp #define</code>s which
may change from run to run. <code>cpp</code> options that do not
change between compilations should be placed in
<code>$(CPPFLAGS)</code>.</p>
<p>If the <code>-w</code> flag is given the commands run are
<code>$(MPIFC)</code> on fortran files, <code>$(MPICC)</code> on C
files, and <code>$(MPILD)</code> for the loader step. The flags
retain their same values with or without the <code>-w</code> flag.
(This is a local addition.)</p>
<p>Includefiles are recursively searched for embedded includes.</p>
<p>For <code>emacs</code> users, the make target <code>TAGS</code>
is always provided. This creates a TAGS file in the current working
directory with a cross-reference table linking all the sourcefiles.
If you don't know about emacs tags, please consult the emacs help
files! It is an incredibly useful feature.</p>
<p>The default action for non-existent files is to
<code>touch</code> them (i.e create null files of that name) in the
current working directory.</p>
<p>All the object files are linked to a single executable. It is
therefore desirable that there be a single main program source
among the arguments to <code>mkmf</code>, otherwise, the loader is
likely to complain.</p>
<h2>Treatment of [args]</h2>
<p>The argument list <code>args</code> is treated sequentially from
left to right. Arguments can be of three kinds:</p>
<ul>
<li>If an argument is a sourcefile, it is added to the list of
sourcefiles.</li>
<li>If an argument is a directory, all the sourcefiles in that
directory are added to the list of sourcefiles.</li>
<li>If an argument is a regular file, it is presumed to contain a
list of sourcefiles. Any line not containing a sourcefile is
discarded. If the line contains more than one word, the last word
on the line should be the sourcefile name, and the rest of the line
is a file-specific compilation command. This may be used, for
instance, to provide compiler flags specific to a single file in
the sourcefile list.</li>
</ul>
<pre>
<code>
a.f90
b.f90
f90 -Oaggress c.f90
</code>
</pre>
<p>This will add <code>a.f90, b.f90</code> and <code>c.f90</code>
to the sourcefile list. The first two files will be compiled using
the generic command <code>$(FC) $(FFLAGS)</code>. But when the make
requires <code>c.f90</code> to be compiled, it will be compiled
with <code>f90 -Oaggress</code>.</p>
<p>The current working directory is always the first (and
top-precedence) argument, even if <code>args</code> is not
supplied.</p>
<h2>Treatment of [-c cppdefs]</h2>
<p>The argument <code>cppdefs</code> is treated as follows.
<code>cppdefs</code> should contain a comprehensive list of the
<code>cpp</code> <code>#define</code>s to be preprocessed. This
list is compared against the current "state", maintained in the
file <code>.cppdefs</code> in the current working directory. If
there are any changes to this state, <code>mkmf</code> will remove
all object files affected by this change, so that the subsequent
<code>make</code> will recompile those files. Previous versions of
<code>mkmf</code> attempted to <code>touch</code> the relevant
source, an operation that was only possible with the right
permissions. The current version works even with read-only
source.</p>
<p>The file <code>.cppdefs</code> is created if it does not exist.
If you wish to edit it by hand (don't!) it merely contains a list
of the <code>cpp</code> flags separated by blanks, in a single
record, with no newline at the end.</p>
<p><code>cppdefs</code> also sets the <code>make</code> macro
<code>CPPDEFS</code>. If this was set in a template file and also
in the <code>-c</code> flag to <code>mkmf</code>, the value in
<code>-c</code> takes precedence. Typically, you should set only
<code>CPPFLAGS</code> in the template file, and
<code>CPPDEFS</code> via <code>mkmf -c</code>.</p>
<h2>Treatment of includefiles</h2>
<p>Include files are often specified without an explicit path,
e.g:</p>
<pre>
<code>
#include "config.h"
</code>
</pre>
<p><code>mkmf</code> first attempts to locate the includefile in
the same directory as the source file. If it is not found there, it
looks in the directories listed as arguments, maintaining the same
left-to-right precedence as described above.</p>
<p>This follows the behaviour of most f90 compilers: includefiles
inherit the path to the source, or else follow the order of include
directories specified from left to right on the <code>f90</code>
command line, with the <code>-I</code> flags <i>descending</i> in
precedence from left to right.</p>
<p>If you have includefiles in a directory <code>dir</code> other
than those listed above, you can specify it yourself by including
<code>-Idir</code> in <code>$(FFLAGS)</code> in your template file.
Includepaths in the template file take precedence over those
generated by <code>mkmf</code>. (I suggest using
<code>FFLAGS</code> for this rather than <code>CPPFLAGS</code>
because fortran <code>include</code>s can occur even in source
requiring no preprocessing).</p>
<h2 id="examples">Examples</h2>
<p>The template file for the SGI MIPSpro compiler contains:</p>
<pre>
<code>
FC = f90
LD = f90
CPPFLAGS = -macro_expand
FFLAGS = -d8 -64 -i4 -r8 -mips4 -O3
LDFLAGS = -64 -mips4 $(LIBS)
LIST = -listing
</code>
</pre>
<p>The meaning of the various flags may be divined by reading the
manual. A line defining the <code>make</code> macro LIBS, e.g:</p>
<pre>
<code>
LIBS = -lmpi
</code>
</pre>
<p>may be added anywhere in the template to have it added to the
link command line.</p>
<p>Sample template files for different OSs and compilers are
available in the directory <code>/net/vb/public/bin</code>.</p>
<p>This example illustrates the effective use of
<code>mkmf</code>'s precedence rules. Let the current working
directory contain a file named <code>path_names</code> containing
the lines:</p>
<pre>
<code>
updates/a.f90
updates/b.f90
</code>
</pre>
<p>The directory <code>/home/src/base</code> contains the
files:</p>
<pre>
<code>
a.f90
b.f90
c.f90
</code>
</pre>
<p>Typing <code>mkmf path_names /home/src/base</code> produces the
following <code>Makefile</code>:</p>
<pre>
<code>
# Makefile created by mkmf

.DEFAULT:
        -touch $@
all: a.out
c.o: /home/src/base/c.f90
        $(FC) $(FFLAGS) -c      /home/src/base/c.f90
a.o: updates/a.f90
        $(FC) $(FFLAGS) -c      updates/a.f90
b.o: updates/b.f90
        $(FC) $(FFLAGS) -c      updates/b.f90
./c.f90: /home/src/base/c.f90
        cp /home/src/base/c.f90 .
./a.f90: updates/a.f90
        cp updates/a.f90 .
./b.f90: updates/b.f90
        cp updates/b.f90 .
SRC = /home/src/base/c.f90 updates/a.f90 updates/b.f90
OBJ = c.o a.o b.o
OFF = /home/src/base/c.f90 updates/a.f90 updates/b.f90
clean: neat
        -rm -f .cppdefs $(OBJ) a.out
neat:
        -rm -f $(TMPFILES)
localize: $(OFF)
        cp $(OFF) .
TAGS: $(SRC)
        etags $(SRC)
tags: $(SRC)
        ctags $(SRC)
a.out: $(OBJ)
        $(LD) $(OBJ) -o a.out $(LDFLAGS)
</code>
</pre>
<p>Note that when files of the same name recur in the target list,
the files in the <code>updates</code> directory (specified in
<code>path_names</code>) are used rather than those in the base
source repository <code>/home/src/base</code>.</p>
<p>Assume that now you want to test some changes to
<code>c.f90</code>. You don't want to make changes to the base
source repository itself prior to testing; so you make yourself a
local copy.</p>
<pre>
<code>
$ make ./c.f90
</code>
</pre>
<p>You didn't even need to know where <code>c.f90</code> originally
was.</p>
<p>Now you can make changes to your local copy
<code>./c.f90</code>. To compile using your changed copy, type:</p>
<pre>
<code>
$ mkmf path_names /home/src/base
$ make
</code>
</pre>
<p>The new Makefile looks like this:</p>
<pre>
<code>
# Makefile created by mkmf

.DEFAULT:
        -touch $@
all: a.out
c.o: c.f90
        $(FC) $(FFLAGS) -c      c.f90
a.o: updates/a.f90
        $(FC) $(FFLAGS) -c      updates/a.f90
b.o: updates/b.f90
        $(FC) $(FFLAGS) -c      updates/b.f90
./a.f90: updates/a.f90
        cp updates/a.f90 .
./b.f90: updates/b.f90
        cp updates/b.f90 .
SRC = c.f90 updates/a.f90 updates/b.f90
OBJ = c.o a.o b.o
OFF = updates/a.f90 updates/b.f90
clean: neat
        -rm -f .cppdefs $(OBJ) a.out
neat:
        -rm -f $(TMPFILES)
localize: $(OFF)
        cp $(OFF) .
TAGS: $(SRC)
        etags $(SRC)
tags: $(SRC)
        ctags $(SRC)
a.out: $(OBJ)
        $(LD) $(OBJ) -o a.out $(LDFLAGS)
</code>
</pre>
<p>Note that you are now using your local copy of
<code>c.f90</code> for the compile, since the files in the current
working directory always take precedence. To revert to using the
base copy, just remove the local copy and run <code>mkmf</code>
again.</p>
<p>This illustrates the use of <code>mkmf -c</code>:</p>
<pre>
<code>
$ mkmf -c "-Dcppflag -Dcppflag2=2 -Dflag3=string ..."
</code>
</pre>
<p>will set <code>CPPDEFS</code> to this value, and also save this
state in the file <code>.cppdefs</code>. If the argument to
<code>-c</code> is changed in a subsequent call:</p>
<pre>
<code>
$ mkmf -c "-Dcppflag -Dcppflag2=3 -Dflag3=string ..."
</code>
</pre>
<p><code>mkmf</code> will scan the source list for sourcefiles that
make references to <code>cppflag2</code>, and the corresponding
object files will be removed.</p>
<h2>Caveats</h2>
<p>In F90, the module name must occur on the same source line as
the <code>module</code> or <code>use</code> keyword. That is to
say, if your code contained:</p>
<pre>
<code>
use &
this_module
</code>
</pre>
<p>it would confuse <code>mkmf</code>. Similarly, a fortran
<code>include</code> statement must not be split across lines.</p>
<p>Two <code>use</code> statements on the same line is not
currently recognized, that is:</p>
<pre>
<code>
use module1; use module2
</code>
</pre>
<p>is to be avoided.</p>
<p><code>mkmf</code> provides a default action for files listed as
dependencies but not found. In this case, <code>mkmf</code> will
<code>touch</code> the file, creating a null file of that name in
the current directory. It is the least annoying way to take care of
a situation when cpp <code>#include</code>s buried within obsolete
<code>ifdef</code>s ask for files that don't exist:</p>
<pre>
<code>
#ifdef obsolete
#include "nonexistent.h"
#endif
</code>
</pre>
<p>If the formatting flag <code>-f</code> is used, long lines will
be broken up at intervals of 256 characters. This can lead to
problems if individual paths are longer than 256 characters.</p>
</body>
</html>
