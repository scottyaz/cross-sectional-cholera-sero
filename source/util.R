## some generic helper functions
## most of these are taken from random places online
## they have been passed along within my util functions for a
## while with no attribution so sorry in advance.

wrapper for sending output from function to device
##' (note: this is used by to pdf to create pdfs)
##' @param expr
##' @param dev
##' @param filename
##' @param ...
##' @param verbose
##' @return
##' @author http://nicercode.github.io/blog/2013-07-09-figure-functions/
to.dev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  dev(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

## to make a pdf
to.pdf <- function(expr, filename, ...)
    to.dev(expr, pdf, filename, ...)

## to make a png
to.png <- function(expr, filename, ...)
    to.dev(expr, png, filename,...)


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


##' Adds alpha to a set of colors
##' @title
##' @param COLORS
##' @param ALPHA
##' @return
AddAlpha <- function(COLORS, ALPHA){
    if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
    RGB <- col2rgb(COLORS, alpha=TRUE)
    RGB[4,] <- round(RGB[4,]*ALPHA)
    NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
    return(NEW.COLORS)
}
