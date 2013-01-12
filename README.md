Parana 2 (Reviewer Pre-release) 
===============================

This software accompanies ISMB 2013 submission 217 "Predicting protein
interactions via parsimonious ancestral network reconstruction" and is for
reviewer purposes only.  It should not currently be redistributed.

Installation
------------

The [binary release](http://www.cs.cmu.edu/~ckingsf/software/parana2/Parana2-0.9.0-Darwin.tar.bz2) 
should be ready to use (if you're on OSX).

The [source version](http://www.cs.cmu.edu/~ckingsf/software/parana2/Parana2-0.9.0-src.tar.bz2) uses 
the [CMake](http://www.cmake.org) build system.  It depends on
a C++11 compatible (e.g. G++-4.7) compiler and the following external libraries:

[Bio++](http://biopp.univ-montp2.fr/)  
[Boost](http://www.boost.org/)  
[GMP](http://gmplib.org/)  
[MPFR](http://www.mpfr.org/)  
[pugixml](http://pugixml.org/)  

From the top level directory, one should execute the following commands 
(> designates the prompt):


    > mkdir build
    > cd build
    > cmake ..
    > make


Any errors are likely the result of a missing library.

Running the program
-------------------

The current _pre-release_ version of the program exposes a number of options that 
are not fully implemented or which may be removed in the final version.  However, all
of the functionality necessary to reproduce the results of the paper are present and
working.  A typical invocation of the program looks something like this:

    > ./parana2 pars -u -t [target network] -d [dup. hist] -o [output file]

### The target network

The target network should be in the NetworkX [adjacency list format](http://networkx.lanl.gov/reference/readwrite.adjlist.html).  The comment character "#" is respected, and the file should list nodes with no incident edges on
lines by themselves.

### The duplication history

The duplication history should be in [PhyloXML](http://www.phyloxml.org/) format (the NHX format will be supported shortly).
Networks currently in the NHX format (specifically those that result from Notung) can be easily converted
to NHX format using the [phyloXML converter](https://sites.google.com/site/cmzmasek/home/software/forester/phyloxml-converter).  For example, if one currently has a file "tree.ntg" (a NHX format tree that is the result of
a gene-species tree reconciliation using Notung), she can obtain the appropriate phyloXML format file
(assuming forester.jar is in the current directory) by running the command:

    java -cp forester.jar org.forester.application.phyloxml_converter -f=nn tree.ntg tree.xml

This will output an appropriately formated phyloXML file, "tree.xml".

### The output file

The output file is a list of posterior scores for all potential edges (for those with scores > 0).
The format is very simple, each line is given as follows:

	p1	p2	et	s

Where p1 and p2 are the names of two proteins, et is the edge-type (currently, this is just 'b'), and
s is the score assigned to this edge.

### Additionally

Additional arguments are supported which affect the scores of different ancestral histories (e.g.
the branch-length penalty and the ratio of cost of edge creation to deletion).  These parameters 
are briefly documented with the program's help option (which can be viewd by running ./parana2 pars -h).
The following parameters all currently work:

    -r [--ratio]
    -b [--beta] (this is called gamma in the paper --- will be fixed)
    -k [--numOpt]
    -p [--timePenalty]

Finally, to save the effort of rebuilding the hypergraph for every experiment when performing
a cross-validation experiment (only the costs of the leaf hyperverties are affected in such
experients), a cross validation flag (`-c [--cross]`) is available.  This command modifies the
meaning of the `-o [--output]` option, so that argument passed to this option will now become
a directory, under which the results of all separate cross validation experiments will be written.