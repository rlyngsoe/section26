/***************************************************************************

    kwarg.c
  
    Implementation of front end for a greedy heuristic for finding plausible
    evolutionary histories with a low number of recombinations under the
    infinite sites assumption for an SNP data set.

    Rune Lyngsø (lyngsoe@stats.ox.ac.uk), March 2007

****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "kwarg.h"
#include "gene.h"
#include "bounds.h"
#include "exact.h"
#include "common.h"
#include "backtrack.h"
#include "expression.h"

static void _print_usage(FILE *f, char *name)
{
  fprintf(f, "Usage: %s [options] [input]\n", name);
  pretty_print(f, "The program reads data from the input file specified (or from standard input if no input is specified) and greedily constructs a history with a low number of recombinations under the infinite sites assumption. The history is constructed by stepping backwards in time using mutation, coalescence, and recombination events. At each point in this process, all possible next events (strictly speaking, only a subset of all possible next events is used - this subset has as least one event that can lead to a history with a minimum number of recombinations from the current point, though) are considered, and the resulting ancestral states are scored. The scores are used to choose an event either at random or to proceed to an ancestral state with minimum score (see options -r and -T). This process is NOT guaranteed to lead to a history with a minimum number of recombinations.", 80, 0);
  fprintf(f, "Legal options are:\n");
  print_option(f, "-#[name]", "Print a # for every recombination as it is used to file name (stdout if no file is specified). Only the last destination specified will actually receive the #s, though any previous destinations will be initialised (i.e. any previous information stored in files will be deleted).", 80, -1);
  print_option(f, "-s", "Suppress output of number of recombinations used", 80, -1);
  print_option(f, "-r", "Select next ancestral state back in time randomly from all considered states, with probability proportional to its score (or if a positive temperature has been specified - see the -T option - with probability proportional to the exponential of minus its score divided by the temperature). Default is to just choose an ancestral state with minimum score uniformly at random. Negative scores are truncated to zero. If no scores are positive, a configuration with maximal score will be chosen uniformly at random.", 80, -1);
  print_option(f, "-TF", "Treat scores as pseudo-energies, where the next configuration is selected according to the corresponding Boltzmann distribution, i.e. the probability of choosing a configuration is proportional to the exponential of minus its score divided by the temperature F. A low temperature will cause a strong bias towards minimum scoring configurations, while a high temperature will make the selection of the next configuration close to uniformly random. Invoking this option with a positive temperature automatically triggers random selection. If it is invoked with temperature equal to zero, selection reverts to selecting among the configurations with minimum score. A negative temperature cancels the pseudo-energy treatment (though random selection will not be cancelled if in effect, even if it was triggered by specifying a positive temperature).", 80, -1);
  print_option(f, "-Eexp", "Score configurations according to expression exp. An expression should consist of zero or more definitions followed by a single statement, all separated by semicolons. The score is the value of the final statement, while the definitions allow variables and functions to be defined and used. A variable is defined using the syntax `name = statement' where name is the name of the variable that is set to the value of statement. A function is defined using the syntax `f(x, y, z) = statement' where f is the function name, x, y, and z the parameters of the function (at least one parameter is required, but there is no upper bound on the number of parameters), and statement the body of the function. Variable, function and parameter names are required to be a sequence of letters and digits, beginning with a letter. Statements can be\n- real or integer numbers\n- variables\n- function calls, with syntax `f(statement1, statement2, statement3)' where the values of arguments statement1, statement2, and statement3 are determined when the function is called; the number of arguments in a function call is required to match the number of parameters in the function definition\n- arithmetic operations: +, -, *, /, and ^ (for exponentiation - this operation is also available using the synonym **) that require numerical arguments and yield numerical results\n- comparisons: <, >, = (with synonym ==), <= (with synonym =<), >= (with synonym =>), and != (with synonyms <> and ><) that require numerical arguments and yield Boolean results\n- Boolean operations: & (and, with synonym &&), | (or, with synonym ||), ! (not, with synonym ¬), and \\ (exclusive or, with synonym \\\\) that require Boolean arguments and yield Boolean values\n- conditionals with the syntax `test ? consequence : alternative' where test, consequence and alternative are statements; the value of test, which is required to be a Boolean, is first determined, and only the relevant branch - consequence if test is True, alternative if test is False - is then evaluated to determine the value of the conditional.\nWhite space - space, tabs and newlines - is ignored. The statement defining a variable cannot refer to the variable itself, only previously defined variables and functions. Conversely, a function is allowed to be recursive, so e.g. factorials can be computed as `fac(n) = n < 2 ? 1 : n * fac(n - 1)'. Reusing a variable name does not change the value of the original definition but introduces a new variable with the same name and the new value. So the expression `a = 1; f(x) = a; a = 2; f(0)' will evaluate to 1. Observe that we need to use a dummy parameter x for f as functions are required to have at least one parameter. Arguments to functions are allowed to have any type as long as the corresponding parameter is used in concordance with the argument type. Functions can even take other functions as arguments, e.g. `fixpoint(f, a) = (f(a) = a ? a : fixpoint(f, f(a)))' finds a fixpoint for f by reapplying f starting from a until the value no longer changes.\nSome built-in functions and variables are available for convenience and to allow probing of the configuration being scored. The names of these can of course be be reused to define new variables and functions. This will hide the built-in functions, making them inaccessible for the remainder of the expression. The convenience functions are log (base 10 logarithm), ln (natural logarithm) and sqrt (square root) and the trigonometric functions sin, cos, tan, asin, acos, and atan. Convenience constants are e and pi. Constants for probing the configuration are\n- r: number of recombinations in step leading to configuration, which is always 0, 1, or 2 - this constant also has synonym recombinations\n- am: total amount of ancestral material in the sequences of the configuration, with synonym ancestralmaterial\n- maxam: maximum value of am over all configurations considered at this point, with synonym maximumancestralmaterial\n- minam: minimum value of am, with synonym minimumancestralmaterial\n- seq: number of sequences in configuration, with synonym numberofsequences\n- maxseq: maximum value of seq over all possible configurations considered at this point, with synonym maximumnumberofsequences\n- minseq: minimum value of seq, with synonym minimumnumberofsequences\n- len: number of sites remaining in configuration - uninformative sites are always removed, so this may be smaller than the original number of sites; it has synonym numberofsites\n- maxlen: maximum value of len over all possible configurations considered at this point, with synonym maximumnumberofsites\n- minlen:  minimum value of len, with synonym minimumnumberofsites\n- rmin: minimum number of recombinations required for configuration - note that for all but small data sets this will be extremely slow to compute; it has synonyms Rmin, beagle and Beagle\n- hk: Hudson-Kaplan bound on number of recombinations required for configuration, with synonyms HK, hudsonkaplan and HudsonKaplan\nFunctions for probing configurations are\n- hb: haplotype bound on number of recombinations required for configuration, with synonyms HB, haplotypebound and HaplotypeBound. It takes three arguments specifying the maximum size of subsets of sites considered, the maximum length of an interval the sites in a subset can be spread over, and the minimal increment of number of sequence types required for a subset of sites to be further expanded, in this order; for large data sets haplotype bound computation can be time consuming, and it is suggested to experiment with the haplotypebound program to find suitable values for the arguments\n- eagl: lower bounds based on Rmin values for subintervals combined with Hudson-Kaplan recombination requirements between pairs of sites using the composite method, with synonym Eagl. It takes one argument specifying the maximum length of intervals for which Rmin is computed; again, this can be time consuming for large data sets.\nBuilt-in constants and function calls that are time consuming to compute are only computed when needed and then stored for future reference.\nNumerical values are represented using the C double data type. Expressions where intermediate values are very large or very small should be avoided, as these may incur significant rounding errors. For example defining the combinatorial choose function as `choose(m, n) = fac(m) / (fac(n) * fac(m - n))' using the factorial function fac defined above will not be as robust as defining it as `choose(m, n) = n = 0 ? 1 : (m / n) * choose(m - 1, n - 1)' or even better as `min(a, b) = a < b ? a : b; choose2(m, n) = n = 0 ? 1 : (m / n) * choose2(m - 1, n - 1); choose(m, n) = choose2(m, min(n, m - n))'. One should also beware of the approximate nature of computers' real numbers. The data type (C doubles) used for numbers can behave unexpectedly, as witnessed by the value of `7 * 2.15 = 15.05 ? 1 : 0' on at least some machines (where 7 * 2.15 does not equal 15.05), and if possible tests for equality or inequality should be replaced by less than or greater than comparisons.\nThe expression used as default scoring scheme is `" KWARG_DEFAULTSCORE "'. It computes a proxy for the minimum number of recombinations required for a configuration and among the configurations minimising the proxy chooses one with a minimum amount of ancestral material left. The proxy is computed as twice the Hudson-Kaplan bound, except for very small data sets where the true minimum number of recombinations required is used. To account for the step leading to a configuration, the number of recombinations in this step is added to the proxy.", 80, -1);
  print_option(f, "-Fname", "Read expression to score configurations by from file name. The format of an expression is as described for the -E option.", 80, -1);
  print_option(f, "-m", "Normally the decision of the next event is taken as the resulting configurations are enumerated. This is very memory efficient, but if any of the built in variables referring to minimum or maximum value of an entity are used in the scoring scheme it will result in an increase of running time up to a factor of two. Use this option to first enumerate all next events and resulting configurations before choosing among them.", 80, -1);
  print_option(f, "-b[name]", "Output a minimum recombination history to file name (stdout if no file is specified).", 80, -1);
  print_option(f, "-d[name]", "Output ancestral recombination graph of minimum recombination history in dot format to file name (stdout if no file is specified).", 80, -1);
  print_option(f, "-g[name]", "Output ancestral recombination graph of minimum recombination history in GDL format to file name (stdout if no file is specified).", 80, -1);
  print_option(f, "-j[name]", "Output ancestral recombination graph of minimum recombination history in GML format to file name (stdout if no file is specified).", 80, -1);
  print_option(f, "-t[name]", "Output list of marginal phylogenies for each site in Newick's 8:45 format to file name (stdout if no file is specified).", 80, -1);
  print_option(f, "-D[name]", "Output list of marginal phylogenies for each site in dot format to file name (stdout if no file is specified).", 80, -1);
  print_option(f, "-G[name]", "Output list of marginal phylogenies for each site in GDL format to file name (stdout if no file is specified).", 80, -1);
  print_option(f, "-J[name]", "Output list of marginal phylogenies for each site in GML format to file name (stdout if no file is specified).", 80, -1);
  print_option(f, "-I", "Marginal trees are only output one for each of the intervals between two recombination points, instead of one for each site.", 80, -1);
  print_option(f, "-vnodelabel", "Use nodelabel convention for labelling nodes in ancestral recombination graphs. The possible conventions are\nnone - do not label nodes\nid - only label nodes representing sampled sequences, using their sequence ids from the data file, and nodes representing recombinations, indicating the recombination point\nsequence - label nodes with the inferred sequences; these sequences will be in binary format, with 0 representing wild type and 1 representing mutant type, even if the original data is not in binary format\nboth - label nodes with both id and inferred sequence\none - use only one label for a node, sequence id or recombination point if avaialable and otherwise the inferred sequence\nDefault convention is id. The colour coding scheme used for the nodes is red for sequences in the input data, blue for recombination nodes, green for standard coalescent nodes, and yellow for the final coalescence into the most recent common ancestor.", 80, -1);
  print_option(f, "-i", "Sequences not having a sequence id in the data file are assigned their index in the data file as id, e.g. the first sequence in the data file would be assigned `1' as id.", 80, -1);
  print_option(f, "-e", "Label edges in ancestral recombination graphs with the sites undergoing mutation along the edge.", 80, -1);
  print_option(f, "-k", "Assume that the common ancestral sequence is known, i.e. that we know which is the wild type and which is the mutant in each site. In binary format the all-0 sequence is assumed to be the common ancestral sequence, while the common ancestral sequence has to be specified directly for amino acid and nucleic sequences (see options -a and -n)", 80, -1);
  print_option(f, "-o", "Assume input data is in own format. Default is to first try to parse data in own format, and if that fails to try to parse it in fasta format. Specifying this option, no attempt will be made to try to parse the data in fasta format.", 80, -1);
  print_option(f, "-f", "Assume input data is in fasta format. No attempt will be made to try to parse the data in own format. Note that the -o and the -f options override each other, so only the last one occurring in the command line will have an effect.", 80, -1);
  print_option(f, "-a", "Assume input consists of amino acid (protein) sequences using the one letter amino acid alphabet; anything not in the amino acid one letter alphabet is treated as an unresolved site (default is to assume sequences in binary, i.e. 0/1, format where anything but a 0 or a 1 is considered an unresolved site). If the most recent common ancestor is assumed known, the first sequence in the input data is considered to specify the wild type of each site and is not included as part of the sample set.", 80, -1);
  print_option(f, "-n", "Assume input consists of nucleic sequences; anything but a/c/g/t/u is considered an unresolved site (default is to assume binary, i.e. 0/1, format where anything but a 0 or a 1 is considered an unresolved site). If the most recent common ancestor is assumed known, the first sequence in the input data is considered to specify the wild type of each site and is not included as part of the sample set.", 80, -1);
  print_option(f, "-h, -H -?", "Print this information and stop.", 80, -1);
}

/* Use the random function to draw an integer between 0 and n - 1 */
static int _unbiased_random(int n)
{
  long int l = XRAND_MAX / n;
  long int i;

  do{
    i = xrandom() / l;
  } while (i >= n); /* i ought to always be at most n, but just to make sure */

  return i;
}

/* Determine whether to select a value according to a random minimal
 * value scheme. The minimum value seen so far is maintained in _ms_v
 * and the number of times this has been encountered is maintained in
 * _ms_w. Return value indicates whether the value should be selected.
 */
static int _ms_w = 0;
static double _ms_v = DBL_MAX;
static int _minimum_select(double a)
{
  if (a == _ms_v)
    /* The two values are equal - choose one at random and increment
     * number of times we have seen this value.
     */
    return (_unbiased_random(++_ms_w) == 0);
  else if (a < _ms_v){
    /* The value of a is new minimum - reset count and report this */
    _ms_v = a;
    _ms_w = 1;
    return 1;
  }

  /* The old value is still the minimum */
  return 0;
}

/* Determine whether to select a value according to a scheme where
 * values are selected with probability proportional to their
 * value. The sum of values seen so far is maintained in _rs_w. Values
 * less than zero or truncated to zero (as long as all values seen so
 * far are non-positive a least negative will be the
 * selection). Return value indicates whether the value should be
 * selected.
 */
static double _rs_w = DBL_MIN;
static int _rs_n = 0;
static int _random_select(double a)
{
  /* Check for non-positive values */
  if (a <= 0){
    if (_rs_w < a){
      _rs_w = a;
      _rs_n = 1;
      return 1;
    }
    else if (_rs_w == a)
      return (_unbiased_random(++_rs_n) == 0);
    else
      return 0;
  }

  /* Update _rs_w and perform random selection */
  if (_rs_w < 0)
    _rs_w = a;
  else
    _rs_w += a;

  return (a * XRAND_MAX > _rs_w * xrandom());
}

/* Determine whether to select a value according to a scheme where
 * values are interpreted as energies and selected according to the
 * corresponding Boltzmann distribution at temperature _prs_kT / k.
 * The partition function so far is maintained in _srs_w and values
 * are shifted such that this is maintained to be close to 1. The
 * temperature is required to be positive. Return value indicates
 * whether the value should be selected.
 */
static double _prs_kT = 1;
static double _prs_Z = 0;
static double _prs_offset = 0;
static int _pseudoenergy_random_select(double a)
{
  if (_prs_Z == 0){
    /* First value seen */
    /* Choose offset such that partition function initially is 1 */
    _prs_offset = a;
    _prs_Z = 1;
    return 1;
  }
  else{
    /* Include contribution of a in partition function and update
     * offset if necessary.
     */
    if (a < _prs_offset){
      /* It's a good idea to change offset before proceeding */
      _prs_Z = exp((a - _prs_offset) / _prs_kT) * _prs_Z + 1;
      _prs_offset = a;
    }
    else
      /* Update partition function */
      _prs_Z += exp((_prs_offset - a) / _prs_kT);
    /* Check whether we should change the offset after the fact */
    if (_prs_Z > 2){
      /* By the precheck, inclusion of a can at most increase the
       * value of the partition function by 1 so reducing by a factor
       * of e should be sufficient.
       */
      _prs_offset -= _prs_kT;
      _prs_Z = _prs_Z / M_E;
    }
  }

  /* Determine whether to choose a */
  if (exp((_prs_offset - a) / _prs_kT) * XRAND_MAX < xrandom() * _prs_Z)
    return 0;
  else
    return 1;
}

/* Reset variables for the various selection functions */
static void _reset_selections()
{
  _ms_w = 0;
  _ms_v = DBL_MAX;
  _rs_w = DBL_MIN;
  _rs_n = 0;
  _prs_Z = 0;
  _prs_offset = 0;
}

/* Parse a floating point option argument and store it in value.
 * Return value states whether the argument could be parsed in full.
 */
static int _parse_double(char *s, double *value)
{
  int i;

  if (sscanf(s, "%lf%n", value, &i) != 1)
    return 0;

  return s[i] == '\0';
}

/* Read entire content of file name into a string */
static char *_read_file(char *name)
{
  FILE *f = stdin;
  char *s = (char *)xmalloc(8 * sizeof(char));
  int size = 0, capacity = 8;

  /* Open file */
  if ((name != NULL) && ((f = fopen(name, "r")) == NULL))
    return NULL;

  /* Read file character by character */
  while ((s[size++] = fgetc(f)) != EOF)
    if (size == capacity - 1){
      /* Ran out of buffer capacity, double its size */
      capacity *= 2;
      s = (char *)xrealloc(s, capacity * sizeof(char));
    }

  /* Zero-terminate string, shrink buffer to fit it, and close file */
  s[size - 1] = '\0';
  s = xrealloc(s, size);
  if (f != stdin)
    fclose(f);

  return s;
}

int main(int argc, char **argv)
{
  Genes *g;
  AnnotatedGenes *a;
  int i, n,
    silent = 0,
    intervals = 0,
    score_from_file = 0;
  Gene_Format format = GENE_ANY;
  Gene_SeqType seqtype = GENE_BINARY;
  FILE * print_progress = NULL;
  char *score = NULL;
  int (*select)(double) = _minimum_select;
  FILE *fp;
  LList *history_files = MakeLList(),
    *dot_files = MakeLList(),
    *gml_files = MakeLList(),
    *gdl_files = MakeLList(),
    *tree_files = MakeLList(),
    *dottree_files = MakeLList(),
    *gmltree_files = MakeLList(),
    *gdltree_files = MakeLList();
  ARG *arg = NULL;
  ARGLabels nodelabel = ARGLABEL;
  int edgelabel = 0;
  int generate_id = 0;
  int ontheflyselection = 1;

  /* Initialise random number generator */
  initialise_xrandom();

#ifdef ENABLE_VERBOSE
  set_verbose(1);
#endif

  /* Analyse command line options */
#define KWARG_OPTIONS "#::srT:E:F::mb::d::g::j::t::D::G::J::Iv:iekofanhH?"

  /* Parse command line options */
  while ((i = getopt(argc, argv, KWARG_OPTIONS)) >= 0){
    switch(i){
    case '#':
      if ((print_progress != NULL) && (print_progress != stdout))
	/* Close existing progress file */
	fclose(print_progress);
      /* Was a file name specified? */
      if (optarg != 0){
	if ((print_progress = fopen(optarg, "w")) == NULL){
	  _print_usage(stderr, argv[0]);
	  fprintf(stderr, "Could not open file %s for output\n", optarg);
	  exit(1);
	}
      }
      else
	print_progress = stdout;
      break;
    case 's':
      silent = 1;
      break;
    case 'r':
      /* Allow random choices of non-minimal steps */
      if (select != _pseudoenergy_random_select)
	select = _random_select;
      break;
    case 'T':
      /* Set new temperature */
      if (!_parse_double(optarg, &_prs_kT)){
	_print_usage(stderr, argv[0]);
	fprintf(stderr, "Could not parse %s as temperature\n", optarg);
	exit(1);
      }
      if (_prs_kT > 0)
	select = _pseudoenergy_random_select;
      else if (_prs_kT == 0)
	select = _minimum_select;
      else
	select = _random_select;
      break;
    case 'E':
      if (score_from_file){
	free(score);
	score_from_file = 0;
      }
      score = optarg;
      break;
    case 'F':
      if (score_from_file)
	free(score);
      if ((score = _read_file(optarg == 0 ? NULL : optarg)) == NULL){
	_print_usage(stderr, argv[0]);
	fprintf(stderr, "\nCould not open file %s for reading\n", optarg);
	exit(1);
      }
      score_from_file = 1;
      break;
    case 'm':
      ontheflyselection = 0;
      break;
    case 'b':
      /* Backtrack history leading to minimum number of recombinations */
      /* Was a file name specified? */
      if (optarg != 0){
	/* Check whether file can be written before initiating compuation */
	if ((fp = fopen(optarg, "w")) == NULL){
	  fprintf(stderr, "Could not open file %s for output\n", optarg);
	  _print_usage(stderr, argv[0]);
	  exit(1);
	}
	fclose(fp);
	Enqueue(history_files, (void *)optarg);
      }
      else
	Enqueue(history_files, stdout);
      break;
    case 'd':
      /* Output ancestral recombination graph of history leading to
       * minimum number of recombinations in dot format.
       */
      /* Was a file name specified? */
      if (optarg != 0){
	/* Check whether file can be written before initiating computation */
	if ((fp = fopen(optarg, "w")) == NULL){
	  fprintf(stderr, "Could not open file %s for output\n", optarg);
	  _print_usage(stderr, argv[0]);
	  exit(1);
	}
	fclose(fp);
	Enqueue(dot_files, (void *)optarg);
      }
      else
	Enqueue(dot_files, stdout);
      break;
    case 'g':
      /* Output ancestral recombination graph of history leading to
       * minimum number of recombinations in gdl format.
       */
      /* Was a file name specified? */
      if (optarg != 0){
	/* Check whether file can be written before initiating computation */
	if ((fp = fopen(optarg, "w")) == NULL){
	  fprintf(stderr, "Could not open file %s for output\n", optarg);
	  _print_usage(stderr, argv[0]);
	  exit(1);
	}
	fclose(fp);
	Enqueue(gdl_files, (void *)optarg);
      }
      else
	Enqueue(gdl_files, stdout);
      break;
    case 'j':
      /* Output ancestral recombination graph of history leading to
       * minimum number of recombinations in gml format.
       */
      /* Was a file name specified? */
      if (optarg != 0){
	/* Check whether file can be written before initiating computation */
	if ((fp = fopen(optarg, "w")) == NULL){
	  fprintf(stderr, "Could not open file %s for output\n", optarg);
	  _print_usage(stderr, argv[0]);
	  exit(1);
	}
	fclose(fp);
	Enqueue(gml_files, (void *)optarg);
      }
      else
	Enqueue(gml_files, stdout);
      break;
    case 't':
      /* Output marginal trees in ancestral recombination graph of
       * history leading to minimum number of recombinations in
       * Newick's 8:45 format.
       */
      /* Was a file name specified? */
      if (optarg != 0){
	/* Check whether file can be written before initiating computation */
	if ((fp = fopen(optarg, "w")) == NULL){
	  fprintf(stderr, "Could not open file %s for output\n", optarg);
	  _print_usage(stderr, argv[0]);
	  exit(1);
	}
	fclose(fp);
	Enqueue(tree_files, (void *)optarg);
      }
      else
	Enqueue(tree_files, stdout);
      break;
    case 'D':
      /* Output marginal trees in ancestral recombination graph of
       * history leading to minimum number of recombinations in
       * dot format.
       */
      /* Was a file name specified? */
      if (optarg != 0){
	/* Check whether file can be written before initiating computation */
	if ((fp = fopen(optarg, "w")) == NULL){
	  fprintf(stderr, "Could not open file %s for output\n", optarg);
	  _print_usage(stderr, argv[0]);
	  exit(1);
	}
	fclose(fp);
	Enqueue(dottree_files, (void *)optarg);
      }
      else
	Enqueue(dottree_files, stdout);
      break;
    case 'G':
      /* Output marginal trees in ancestral recombination graph of
       * history leading to minimum number of recombinations in
       * GDL format.
       */
      /* Was a file name specified? */
      if (optarg != 0){
	/* Check whether file can be written before initiating computation */
	if ((fp = fopen(optarg, "w")) == NULL){
	  fprintf(stderr, "Could not open file %s for output\n", optarg);
	  _print_usage(stderr, argv[0]);
	  exit(1);
	}
	fclose(fp);
	Enqueue(gdltree_files, (void *)optarg);
      }
      else
	Enqueue(gdltree_files, stdout);
      break;
    case 'J':
      /* Output marginal trees in ancestral recombination graph of
       * history leading to minimum number of recombinations in
       * GML format.
       */
      /* Was a file name specified? */
      if (optarg != 0){
	/* Check whether file can be written before initiating computation */
	if ((fp = fopen(optarg, "w")) == NULL){
	  fprintf(stderr, "Could not open file %s for output\n", optarg);
	  _print_usage(stderr, argv[0]);
	  exit(1);
	}
	fclose(fp);
	Enqueue(gmltree_files, (void *)optarg);
      }
      else
	Enqueue(gmltree_files, stdout);
      break;
    case 'I':
      intervals = 1;
      break;
    case 'v':
      if (!strcmp(optarg, "none"))
	nodelabel = ARGNONE;
      else if (!strcmp(optarg, "id"))
	nodelabel = ARGLABEL;
      else if (!strcmp(optarg, "sequence"))
	nodelabel = ARGSEQUENCE;
      else if (!strcmp(optarg, "both"))
	nodelabel = ARGBOTH;
      else if (!strcmp(optarg, "one"))
	nodelabel = ARGLABELFIRST;
      else{
	fprintf(stderr, "Unrecognised nodelabel convention (%s) for -v option\n", optarg);
	_print_usage(stderr, argv[0]);
	exit(1);
      }
      break;
    case 'i':
      generate_id = 1;
      break;
    case 'e':
      edgelabel = 1;
      break;
    case 'k':
      gene_knownancestor = 1;
      break;
    case 'o':
      format = GENE_BEAGLE;
      break;
    case 'f':
      format = GENE_FASTA;
      break;
    case 'a':
      seqtype = GENE_AMINO;
      break;
    case 'n':
      seqtype = GENE_NUCLEIC;
      break;
    case 'h':
    case 'H':
    case '?':
      _print_usage(stdout, argv[0]);
      /* Clean up */
      DestroyLList(dot_files);
      DestroyLList(gml_files);
      DestroyLList(gdl_files);
      DestroyLList(tree_files);
      DestroyLList(dottree_files);
      DestroyLList(gmltree_files);
      DestroyLList(gdltree_files);
      DestroyLList(history_files);
      if (score_from_file)
	free(score);
      exit(0);
    case ':':
      _print_usage(stderr, argv[0]);
      exit(1);
    }
  }

  /* Set default scoring scheme if none has been specified */
  if (score == NULL){
    score = KWARG_DEFAULTSCORE;
    if (select == _random_select){
      /* Default scoring scheme is always treated as pseudoenergy for
       * random selection.
       */
      select = _pseudoenergy_random_select;
      if (_prs_kT <= 0)
	_prs_kT = 1;
    }
  }

  /* Read data */
  if (argc > optind){
    if ((a = read_genes(argv[optind], format, seqtype)) == NULL){
      fprintf(stderr, "Could not parse %s as valid data\n", argv[optind]);
      exit(1);
    }
  }
  else{
    if ((a = read_genes(NULL, format, seqtype)) == NULL){
      fprintf(stderr, "Could not parse input as valid data\n");
      exit(1);
    }
  }
  if ((gene_knownancestor) && (seqtype != GENE_BINARY))
    /* First sequnce only included to specify known common ancestor */
    remove_annotatedgene(a, 0);
  g = a->g;

#ifdef ENABLE_VERBOSE
  output_annotatedgenes(a, NULL, NULL);
#endif

  /* Set up structures for computation */
  if ((Length(history_files) > 0) || (Length(dot_files) > 0)
      || (Length(gml_files) > 0) || (Length(gdl_files) > 0)
      || (Length(tree_files) > 0) || (Length(dottree_files) > 0)
      || (Length(gmltree_files) > 0) || (Length(gdltree_files) > 0))
    eventlist = MakeLList();

  /* Find plausible evolutionary history */
  n = greedy(a->g, print_progress, score, select, _reset_selections,
	     ontheflyselection);

  /* Report on computation */
  /* Start with number of recombinations, if not suppressed */
  if (!silent)
    printf("Found history with %d recombination%s\n", n, (n != 1 ? "s" : ""));

  /* Output inferred ARG */
  if ((Length(history_files) > 0) || (Length(dot_files) > 0)
	   || (Length(gml_files) > 0) || (Length(gdl_files) > 0)
	   || (Length(tree_files) > 0) || (Length(dottree_files) > 0)
      || (Length(gmltree_files) > 0) || (Length(gdltree_files) > 0)){
    while ((fp = (FILE *)Pop(history_files)) != NULL){
      if (fp != stdout)
	/* Open named file for output */
	if ((fp = fopen((char *)fp, "w")) == NULL){
	  fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
	  continue;
	}
      /* Only remember last ARG constructed (they should all be the same) */
      if (arg != NULL)
	arg_destroy(arg);
      arg = eventlist2history(a, fp);
    }
    if (arg == NULL)
      arg = eventlist2history(a, NULL);
    if (arg != NULL){
      /* Output ARG in dot format */
      while ((fp = (FILE *)Pop(dot_files)) != NULL){
	if (fp != stdout)
	  /* Open named file for output */
	  if ((fp = fopen((char *)fp, "w")) == NULL)
	    fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
	arg_output(arg, a, fp, ARGDOT, nodelabel, edgelabel, generate_id);
	if (fp != stdout)
	  fclose(fp);
      }
      /* Output ARG in GML format */
      while ((fp = (FILE *)Pop(gml_files)) != NULL){
	if (fp != stdout)
	  /* Open named file for output */
	  if ((fp = fopen((char *)fp, "w")) == NULL)
	    fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
	arg_output(arg, a, fp, ARGGML, nodelabel, edgelabel, generate_id);
	if (fp != stdout)
	  fclose(fp);
      }
      /* Output ARG in dot format */
      while ((fp = (FILE *)Pop(gdl_files)) != NULL){
	if (fp != stdout)
	  /* Open named file for output */
	  if ((fp = fopen((char *)fp, "w")) == NULL)
	    fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
	arg_output(arg, a, fp, ARGGDL, nodelabel, edgelabel, generate_id);
	if (fp != stdout)
	  fclose(fp);
      }
      /* Output marginal trees of ARG in Newick's 8:45 format */
      while ((fp = (FILE *)Pop(tree_files)) != NULL){
	if (fp != stdout)
	  /* Open named file for output */
	  if ((fp = fopen((char *)fp, "w")) == NULL)
	    fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
	arg_output(arg, a, fp, TREENEWICK, nodelabel, edgelabel, generate_id,
		   intervals);
	if (fp != stdout)
	  fclose(fp);
      }
      /* Output marginal trees of ARG in dot format */
      while ((fp = (FILE *)Pop(dottree_files)) != NULL){
	if (fp != stdout)
	  /* Open named file for output */
	  if ((fp = fopen((char *)fp, "w")) == NULL)
	    fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
	arg_output(arg, a, fp, TREEDOT, nodelabel, edgelabel, generate_id,
		   intervals);
	if (fp != stdout)
	  fclose(fp);
      }
      /* Output marginal trees of ARG in GML format */
      while ((fp = (FILE *)Pop(gmltree_files)) != NULL){
	if (fp != stdout)
	  /* Open named file for output */
	  if ((fp = fopen((char *)fp, "w")) == NULL)
	    fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
	arg_output(arg, a, fp, TREEGML, nodelabel, edgelabel, generate_id,
		   intervals);
	if (fp != stdout)
	  fclose(fp);
      }
      /* Output marginal trees of ARG in GDL format */
      while ((fp = (FILE *)Pop(gdltree_files)) != NULL){
	if (fp != stdout)
	  /* Open named file for output */
	  if ((fp = fopen((char *)fp, "w")) == NULL)
	    fprintf(stderr, "Could not open file %s for output\n", (char *)fp);
	arg_output(arg, a, fp, TREEGDL, nodelabel, edgelabel, generate_id,
		   intervals);
	if (fp != stdout)
	  fclose(fp);
      }
      arg_destroy(arg);
    }
    if (eventlist != NULL){
      while (Length(eventlist) > 0)
	free(Pop(eventlist));
      DestroyLList(eventlist);
    }
  }

  /* Clean up */
  DestroyLList(dot_files);
  DestroyLList(gml_files);
  DestroyLList(gdl_files);
  DestroyLList(tree_files);
  DestroyLList(dottree_files);
  DestroyLList(gmltree_files);
  DestroyLList(gdltree_files);
  DestroyLList(history_files);
  free_annotatedgenes(a);
  if (score_from_file)
    free(score);

  return 0;
}
