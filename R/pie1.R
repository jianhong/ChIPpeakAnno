#' Pie Charts
#' 
#' Draw a pie chart with percentage
#' 
#' 
#' @param x a vector of non-negative numerical quantities.  The values in x are
#' displayed as the areas of pie slices.
#' @param labels one or more expressions or character strings giving names for
#' the slices.  Other objects are coerced by as.graphicsAnnot. For empty or NA
#' (after coercion to character) labels, no label nor pointing line is drawn.
#' @param edges the circular outline of the pie is approximated by a polygon
#' with this many edges.
#' @param radius the pie is drawn centered in a square box whose sides range
#' from -1 to 1.  If the character strings labeling the slices are long it may
#' be necessary to use a smaller radius.
#' @param clockwise logical indicating if slices are drawn clockwise or counter
#' clockwise (i.e., mathematically positive direction), the latter is default.
#' @param init.angle number specifying the starting angle (in degrees) for the
#' slices. Defaults to 0 (i.e., "3 o'clock") unless clockwise is true where
#' init.angle defaults to 90 (degrees), (i.e., "12 o'clock").
#' @param density the density of shading lines, in lines per inch. The default
#' value of NULL means that no shading lines are drawn. Non-positive values of
#' density also inhibit the drawing of shading lines.
#' @param angle the slope of shading lines, given as an angle in degrees
#' (counter-clockwise).
#' @param col a vector of colors to be used in filling or shading the slices.
#' If missing a set of 6 pastel colours is used, unless density is specified
#' when par("fg") is used.
#' @param border,lty (possibly vectors) arguments passed to polygon which draws
#' each slice.
#' @param main an overall title for the plot.
#' @param percentage logical. Add percentage in the figure or not. default
#' TRUE.
#' @param rawNumber logical. Instead percentage, add raw number in the figure
#' or not. default FALSE.
#' @param digits When set percentage as TRUE, how many significant digits are
#' to be used for percentage. see \link[base]{format}. default 3.
#' @param cutoff When percentage is TRUE, if the percentage is lower than
#' cutoff, it will NOT be shown. default 0.01.
#' @param legend logical. Instead of lable, draw legend for the pie. default,
#' FALSE.
#' @param legendpos,legendcol legend position and legend columns. see
#' \link[graphics]{legend}
#' @param radius.innerlabel position of percentage or raw number label relative
#' to the circle.
#' @param ...  graphical parameters can be given as arguments to pie. They will
#' affect the main title and labels only.
#' @author Jianhong Ou
#' @seealso \code{\link[graphics]{pie}}
#' @keywords misc
#' @export
#' @importFrom grDevices as.graphicsAnnot dev.hold dev.flush 
#' @importFrom graphics plot.new plot.window polygon lines text title
#' @examples
#' 
#' pie1(1:5)
#' 
pie1 <- function (x, labels = names(x), edges = 200, 
                  radius = 0.8, clockwise = FALSE, 
                  init.angle = if (clockwise) 90 else 0, 
                  density = NULL, angle = 45, 
                  col = NULL, border = NULL, lty = NULL, 
                  main = NULL, percentage=TRUE, rawNumber=FALSE, 
                  digits=3, cutoff=0.01, 
                  legend=FALSE, legendpos="topright", legendcol=2, 
                  radius.innerlabel = radius, ...) 
{
    if (!is.numeric(x) || any(is.na(x) | x < 0)) 
        stop("'x' values must be positive.")
    if (is.null(labels)) 
        labels <- as.character(seq_along(x))
    else labels <- as.graphicsAnnot(labels)
    rawX <- x
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    plot.new()
    pin <- par("pin")
    xlim <- ylim <- c(-1, 1)
    if (pin[1L] > pin[2L]) 
        xlim <- (pin[1L]/pin[2L]) * xlim
    else ylim <- (pin[2L]/pin[1L]) * ylim
    dev.hold()
    on.exit(dev.flush())
    plot.window(xlim, ylim, "", asp = 1)
    if (is.null(col)) 
        col <- if (is.null(density)) 
            c("white", "lightblue", "mistyrose", "lightcyan", 
              "lavender", "cornsilk", "pink")
    else par("fg")
    if (!is.null(col)) 
        col <- rep_len(col, nx)
    if (!is.null(border)) 
        border <- rep_len(border, nx)
    if (!is.null(lty)) 
        lty <- rep_len(lty, nx)
    angle <- rep(angle, nx)
    if (!is.null(density)) 
        density <- rep_len(density, nx)
    twopi <- if (clockwise) 
        -2 * pi
    else 2 * pi
    t2xy <- function(t) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = radius * cos(t2p), y = radius * sin(t2p))
    }
    for (i in 1L:nx) {
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
        polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
                border = border[i], col = col[i], lty = lty[i])
        if(!legend){
            P <- t2xy(mean(x[i + 0:1]))
            lab <- as.character(labels[i])
            if (!is.na(lab) && nzchar(lab)) {
                lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
                text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE, 
                     adj = ifelse(P$x < 0, 1, 0), ...)
            }
        }
    }
    if (percentage) {
        for (i in 1L:nx){
            if(dx[i]>cutoff){
                P <- t2xy(mean(x[i + 0:1]))
                text(radius.innerlabel * P$x, radius.innerlabel * P$y, 
                     paste(formatC(dx[i]*100, digits=digits), "%", sep=""), 
                     xpd = TRUE, 
                     adj = .5, ...)
            }
        }
    }else{
        if(rawNumber){
            for (i in 1L:nx){
                if(dx[i]>cutoff){
                    P <- t2xy(mean(x[i + 0:1]))
                    text(radius.innerlabel * P$x, radius.innerlabel * P$y, 
                         rawX[i], xpd = TRUE, 
                         adj = .5, ...)
                }
            }
        }
    }
    if(legend) legend(legendpos, legend=labels, fill=col, 
                      border="black", bty="n", ncol = legendcol)
    title(main = main, ...)
    invisible(NULL)
}
