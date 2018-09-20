\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{mixstock News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 0.9.9 (1/9/2014)}{
  \itemize{
    \item don't use \code{uml.grad} either (this means
    that CML/UML estimation will be a little bit slower
    and more likely to fail; please contact the maintainer
    if this is an important issue)
  }
}
\section{CHANGES IN VERSION 0.9.8 (1/9/2014)}{
  \itemize{
    \item don't use \code{cml.grad} any more (subtle bug in C code that
    leads to errors under valgrind/segfaults under some CRAN tests)
    \item vignette cleanup
  }
}

\section{CHANGES IN VERSION 0.9.7 (9/3/2014)}{
  \itemize{
    \item cleanup for CRAN
  }

}
\section{CHANGES IN VERSION 0.9.6 (3/3/2014)}{
  \itemize{
    \item update Depends/Suggests/Imports
    \item change vignette to knitr, update/clean up a bit
  }
}
\section{CHANGES IN VERSION 0.9.5.3 (23/10/2012)}{
 \itemize{
   \item fix bugs in many-to-many confint() conversion, naming
   \item add as.data.frame for many-to-many mixstock.est objects (add plyr dependency)
 }
}

\section{CHANGES IN VERSION 0.9.5.2 (17/07/2012)}{
 \itemize{
    \item clean out spurious RCS/auto files; version bump
 }
}

\section{CHANGES IN VERSION 0.9.5.1 (17/07/2012)}{
 \itemize{
    \item clean out spurious .RData file; version bump
 }
}

\section{CHANGES IN VERSION 0.9.5 (8/2/2012)}{
  \itemize{
    \item add NAMESPACE
    \item tweak source code to avoid warnings; comment out unused BUGSout, etc.
    \item allow r*(m+1) informative prior in many-to-many (suggestion of Julia Reisser)
  }
}

\section{CHANGES IN VERSION 0.9.3 (4/2011)}{
  \itemize{
  \item document/boost/test for hardcoded limits in C code (thanks to Jennifer Mandel)
  \item show.data arg in summary.mixstock.est
  \item added resample.index to internal storage
  \item improved xyplot
  \item added hooks for coeftab/coefplot2 package
  \item changed to NEWS.Rd format
}
}

\section{CHANGES IN VERSION 0.9.2 2008/10/02}{
  \itemize{
    \item minor tweaks to Rd files
    \item fix manymany.txt and TO code for double definition
    \item changed working directory for recent R2WinBUGS changes
    \item JAGS now works (but not with "BB" code)
    }}

\section{CHANGES IN VERSION 0.9.1 2008/08/18}{
\itemize{
\item minor vignette improvements
\item use layout() for legends, check for text truncation
\item added R2jags support
\item tweaked manymany.txt for R2jags/WinBUGS compatibility
\item changed "manymany.bug" call to "manymany.txt" call to work around R2WinBUGS attempt to write to system dir
\item added silent=TRUE to several try() clauses
}}


\section{CHANGES IN VERSION 0.9 (2008/05/14)}{
\itemize{
  \item version bump for move to r-forge
  }}

  \section{CHANGES IN VERSION 0.8.5}{
  \itemize{
    \item another plotting bug
    \item remove WINE/useWINE/WINEPATH from pm.wbugs, mm.wbugs,
  sourcesize.wbugs  (can be passed to bugs via ... ?)
  \item	tweaked plot.mixstock.data (horizontal plotting now default, more legend options)
  \item tweaked plot.mixstock.est (more options for many-to-many plots)
  \item added prior settings for mixed stocks in many-to-many analysis
}}

\section{CHANGES IN VERSION 0.8.4 }{
  \itemize{
    \item identified (but didn't fix) long-file-path bug for R2WinBUGS
    \item messed with internal structure of resample, confint, etc. for
  many-to-many results. $resample is the sims.matrix object (flat
  matrix, parameters x chain); $resamplist is the sims.list object
  (chain x parameter matrix)
  \item slight changes to summary()
  \item added nsource(), nmix(), nmark() to try to abstract some of this info
  \item record timing info for R2WinBUGS runs
  \item tiny tweaks to C code to satisfy '-Wall -pedantic' flags
  }}



  \section{CHANGES IN VERSION 0.8.4.2}{
    \itemize{
      \item minor changes for documentation compliance with
    checkDocStyle, codetools:::checkUsagePackage("mixstock")
}}

\section{CHANGES IN VERSION 0.8.4.1}{
  \itemize{
    \item fixed plotting bug (!)
}}

\section{CHANGES IN VERSION 0.8.3}{
  \itemize{
    \item small bugs in as.mixstock.est for many-to-many (had M and R reversed in several places)
    \item started to add support for sourcesize in mixstock.data objects
    \item moved bugs code locators inside functions (LazyLoad problem)
}}

\section{CHANGES IN VERSION 0.8.2}{
  \itemize{
    \item fixed mixed stock labeling stuff
    \item tweaking many-to-many code (TO-style BUGS code)
    \item tiny tweak in utils.c/h to avoid sign() conflict with R
    \item simmixstock2 allows random start
    \item added stub for UML start for mm.wbugs
}}

\section{CHANGES IN VERSION 0.8.1 }{
  \itemize{
    \item added a bit more support for many-to-many
}}

\section{CHANGES IN VERSION 0.8}{
  \itemize{
    \item first version called "mixstock"
    \item change names to mixstock, de-turtle-ize
    \item adding multiple-mixed-stock support
    \item changed n.iter/n.burn/n.thin to match R2WinBUGS default
}}

\section{CHANGES IN VERSION 5.2.03}{
  \itemize{
    \item fixed bug in plot.turtle.est (changed cex to cex.axis)
}}

\section{CHANGES IN VERSION 0.3.0}{
  \itemize{
    \item add preliminary R2WinBUGS functionality
}}

\section{CHANGES IN VERSION 0.2.0}{
  \itemize{
    \item more general cleanup, remove bbmisc dependency (add plotrix
    dependency), change plot formats, fix up vignette ...
    }}

    \section{CHANGES IN VERSION 0.1.3.3}{
      \itemize{
	\item version bump for R 2.0.x  should fix all the stuff!
	}}
	\section{CHANGES IN VERSION 0.1.4}{
	  \itemize{
    \item general cleanup
}}










	
