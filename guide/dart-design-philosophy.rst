DART’s design philosophy
========================

In this section we cover DART’s design philosophy. Understanding this philosophy
will make it easier to get started with DART, as you will quickly be able to
predict how and where to find a particular feature of DART.

The main design goals of DART are to:

1. Create a system that is **coherent** and **easy to understand**. DART is
   carefully engineered to have self-contained programs that each do one job and
   do it well. Likewise, DART *just does DA, and does it well.*
2. Release source code that is **as compatible as possible** with the widest
   possible number of systems. The code is written in Fortran 90, which is one
   of the lowest possible common denominators available on virtually all
   systems. See the section `Why Fortran?`_ if this seems like a
   questionable decision to you in this modern world of Matlab, C++, Java,
   Python, Go, etc.
3. Strive to **limit library dependencies**. There is only one required
   dependency of DART: netCDF. Many modern systems have 10s or 100s of
   dependencies, each of which introduces complexity and the potential for bugs,
   lack of support, broken backwards compatibility, etc. If you’ve ever been
   frustrated struggling to debug relationships to packages you’ve never even
   heard of, you are likely to appreciate this DART design goal. Of course,
   there is nothing to stop you from using whatever dependencies you require,
   for example, to collect observations for the ``obs_seq.out`` in an OSE case,
   but DART by design will remain separate from that dependency for you and all
   other users.
4. Only **compile the code you need**. If you are only using a single model for
   your experiments, there is no reason to compile or even touch code for
   another model you never plan to use. Likewise, if you are not using a
   particular observation operator in your experiment, there is also no need to
   compile it or let it cause you headaches. DART recognizes this fact, and
   through the use of the *mkmf* utility and the *preprocess* program, only what
   you need will ever be compiled.
5. Use **explicit interfaces** to enforce contract programming. In practice this
   means that it is easy to add new models, observations operators, data
   assimilation algorithms, etc. as long as they can implement the required
   interface. This approach *allows all of the benefits of object-oriented
   programming without the added complexity for the end user.*
6. Provide results that are **reliable** and meaningful. The DART algorithms are
   carefully tested and maintained in order to be quickly published along with
   appropriate analysis. In a world of chaos, being able to quantify and shrink
   forecast uncertainty via data assimilation in a reliable way is a valuable
   tool for research and operations and everything in between.

In short, DART is designed at each step to make it as easy as possible for users
to get up and running with their models, observations, and possibly even data
assimilation algorithm advances.

Why Fortran?
------------

Many users new to scientific computing such as graduate students raise their
eyebrows when they first hear that a program uses Fortran for active
development. Fortran is considered by many outside (and some inside) of the
scientific computing community to be a dinosaur, old and decrepit, and not
worthy of serious attention. However, this view is short-sighted. There is a
Chinese idiom 喜新厭舊, which means “to love the new and loathe the old,”
indicating that just because something is old does not automatically make it
bad.

While Fortran does have some outdated features that are far removed from the
mainstream of software engineering (such as implicit typing by first initial of
the variable), these can all be disabled, and the stylistic rules for
easy-to-read, modern Fortran are always followed by DART. On the other hand,
Fortran has many other attractive features that make it a top choice for modern
scientific computing. In particular, Fortran offers vectorization of matrices
that make it possible to operate on entire elements of an array at once or
perform linear algebra operations on multi-dimensional arrays. With or without
the use of the *colon operator* (:), Fortran multi-dimensional array support
makes mathematical algorithms easier to read than the equivalent code written in
many other languages. This highly intuitive Fortran syntax was adopted by
Matlab, NumPy, and other languages. Furthermore, for parallel programs using
distributed memory in *MPI*, Fortran remains a top choice along with C and C++
when considering performance. Python code, for example, remains difficult to
parallelize via MPI, not to mention the difficulties in supporting Python 2,
Python 3, pip, anaconda, virtualenv, …

Altogether, for large mathematically-oriented programs that need to be parallel,
Fortran remains a top choice, especially considering the needs of DART:

1. DART does data assimilation, which is primarily mathematically-oriented
   operations on large data sets.
2. DART needs to be parallel with MPI to run on modern supercomputers.
3. Many users of DART are not software development professionals and appreciate
   straightforward and easily understandable code.
4. DART source distributions should be easy to compile and run reliably on many
   different systems. In practice this means avoiding software features that
   might not be supported on all compilers or systems.

With these considerations in mind, the choice of Fortran for DART development is
clear. DART remains highly successful by keeping things simple and *not fixing
what is not broken* even if it isn’t shiny and new.
