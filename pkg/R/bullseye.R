##
## PLOT A FAN TREE, WITH BULLSEYE LEGEND AND AXIS, AND OPTIONAL COLORS
## FOR TIPS
##
## Author: Thibaut Jombart, May 2013.
## t.jombart@imperial.ac.uk
##

############
## bullseye
############
bullseye <- function(phy, traits=NULL, type="fan", col.tips.by=NULL, col.pal=seasun,
                     n.circ=6, bg.circ=transp("royalblue",.1), circ.unit=NULL,
                     legend=TRUE, posi.leg="bottomleft", leg.title="",
                     ...){
    ## CHECKS ##
    if(inherits(phy, c("phylo4","phylo4d"))) phy <- as(phy, "phylo")
    if(!is.list(col.pal)) col.pal <- c(col.pal)

    ## REORDER DATA BY TIP LABEL ##
    ## traits
    if(!is.null(traits) && !is.null(row.names(traits))){
        traits <- traits[phy$tip.label,,drop=FALSE]
    }
    ## col.tips.by
    if(!is.null(col.tips.by) && !is.null(names(col.tips.by))){
        col.tips.by <- col.tips.by[phy$tip.label]
    }


    ## PLOT THE PHYLOGENY
    ## handle color info
    leg.txt <- NULL
    if(!is.null(col.tips.by)){
        tip.col.info <- any2col(col.tips.by, col.pal=col.pal[[1]])
        plot(phy, type="fan", tip.col=tip.col.info$col, ...)
        leg.col <- tip.col.info$leg.col
        leg.txt <- tip.col.info$leg.txt
    } else{
        plot(phy, type="fan", ...)
    }

    ## HANDLE THE 'BULLSEYE' ##
    ## window setting
    oxpd <- par("xpd")
    par(xpd=TRUE)
    on.exit(par(oxpd))

    ## annot info
    if(is.null(circ.unit)){
        annot.max <- 0.5*diff(par("usr")[1:2])
        annot.dist <- seq(from=0, to=annot.max, length=n.circ)
    } else {
        annot.dist <- seq(from=0, by=circ.unit, length=n.circ)
        annot.max <- max(annot.dist)
    }
    
    ## trace the disks
    symbols(rep(0,n.circ), rep(0,n.circ), circ=annot.dist, inches=FALSE,
        bg=bg.circ, fg=NA, add=TRUE)

    ## axis annotation
    segments(-annot.dist[2],0,-annot.dist[3],0)
    text(-mean(annot.dist[2:3]),-annot.dist[2]/5,
         label=format(annot.dist[2], scientific=TRUE, digits=3),cex=.7)

    ## legend info
    if(!is.null(legend) && !is.null(leg.txt)){
        legend(x=posi.leg, legend=leg.txt, fill=leg.col, title=leg.title)
    }

    leg.info <- list(posi=posi.leg, col=leg.col, txt=leg.txt)
    return(invisible(leg.info))
} # end bullseye
