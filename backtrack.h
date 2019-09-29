/*******************************************************************

    backtrack.h

    Description of function and data types for backtracking a minimum
    recombination history

    Rune Lyngsø (lyngsoe@stats.ox.ac.uk), August 2005

********************************************************************/

#ifndef BACKTRACK_H
#define BACKTRACK_H

#include "arg.h"
#include "gene.h"

typedef enum {SUBSTITUTION, COALESCENCE, RECOMBINATION, REMOVE,
	      COLLAPSE, SWAP, LOOKUP} EventType;

typedef struct _Event {
  EventType type;
  union {
    struct {
      int seq;
      int site;
    } s;
    struct {
      int s1;
      int s2;
    } c;
    struct {
      int seq;
      int pos;
    } r;
    struct {
      int s1;
      int s2;
    } swap;
    int collapse;
    int remove;
    int lookup;
  } event;
} Event;

typedef struct _HistoryFragment {
  Genes *g;           /* End configuration */
  LList *event;       /* List of events leading from start
		       * configuration to end configuration.
		       */
  int recombinations; /* Number of recombination events */
} HistoryFragment;

#ifdef DEBUG
extern HashTable *ancestral_state_trace;
void output_fragment(HistoryFragment *h);
#endif
ARG *eventlist2history(AnnotatedGenes *a, FILE *output);

#endif
