\name{pie1}
\alias{pie1}
\title{
  Pie Charts
}
\description{
  Draw a pie chart with percentage
}
\usage{
pie1(x, labels = names(x), edges = 200, 
                  radius = 0.8, clockwise = FALSE, 
                  init.angle = if (clockwise) 90 else 0, 
                  density = NULL, angle = 45, 
                  col = NULL, border = NULL, lty = NULL, 
                  main = NULL, percentage=TRUE, rawNumber=FALSE, 
                  digits=3, cutoff=0.01, 
                  legend=FALSE, legendpos="topright", legendcol=2, ...)
}
\arguments{
  \item{x}{
  a vector of non-negative numerical quantities. 
  The values in x are displayed as the areas of pie slices.
}
  \item{labels}{	
  one or more expressions or character strings giving names for the slices. 
  Other objects are coerced by as.graphicsAnnot. For empty or NA (after 
  coercion to character) labels, no label nor pointing line is drawn.
}
  \item{edges}{	
  the circular outline of the pie is approximated by a polygon with this 
  many edges.
}
  \item{radius}{	
  the pie is drawn centered in a square box whose sides range from -1 to 1. 
  If the character strings labeling the slices are long it may be necessary 
  to use a smaller radius.
}
  \item{clockwise}{	
  logical indicating if slices are drawn clockwise or counter clockwise (i.e.,
  mathematically positive direction), the latter is default.
}
  \item{init.angle}{	
  number specifying the starting angle (in degrees) for the slices. Defaults 
  to 0 (i.e., "3 o'clock") unless clockwise is true where init.angle defaults 
  to 90 (degrees), (i.e., "12 o'clock").
}
  \item{density}{
  the density of shading lines, in lines per inch. The default value of NULL 
  means that no shading lines are drawn. Non-positive values of density also 
  inhibit the drawing of shading lines.
}
  \item{angle}{
  the slope of shading lines, given as an angle in degrees (counter-clockwise).
}
  \item{col}{
  a vector of colors to be used in filling or shading the slices. If missing a 
  set of 6 pastel colours is used, unless density is specified when par("fg") 
  is used.
}
  \item{border, lty}{
  (possibly vectors) arguments passed to polygon which draws each slice.
}
  \item{main}{	
  an overall title for the plot.
} 
  \item{percentage}{
    logical. Add percentage in the figure or not. default TRUE.
  }
  \item{rawNumber}{
    logical. Instead percentage, add raw number in the figure or not. default 
    FALSE.
  }
  \item{digits}{
    When set percentage as TRUE, how many significant digits are to be used for
    percentage. see \link[base]{format}. default 3.
  }
  \item{cutoff}{
    When percentage is TRUE, if the percentage is lower than cutoff, it will
    NOT be shown. default 0.01.
  }
  \item{legend}{
    logical. Instead of lable, draw legend for the pie. default, FALSE.
  }
  \item{legendpos, legendcol}{
    legend position and legend columns. see \link[graphics]{legend}
  }
  \item{...}{	
  graphical parameters can be given as arguments to pie. They will affect the main title and labels only.
}
}
\author{
 Jianhong Ou
}
\seealso{
\code{\link[graphics]{pie}}
}
\examples{
pie1(1:5)
}

\keyword{ misc }
