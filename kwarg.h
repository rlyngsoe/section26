/*******************************************************************

    kwarg.h
  
    Defintion of default expression for kwarg
		
    Rune Lyngsoe (lyngsoe@stats.ox.ac.uk), July 2007

********************************************************************/

#ifndef _KWARG_H
#define _KWARG_H

/* Default expression for scoring configurations */
#define KWARG_DEFAULTSCORE "proxy = (maxam < 50 ? r + rmin : .2 * (r + hk)); floor(proxy * maxam) + am"

#endif
