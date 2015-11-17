##################################################################################################################
# Copied from RQTL package wihout any change
#----------------------------------------------------------------------#
#  R/qtl															   #
#  copyright (c) 2001-2014, Karl W Broman							   #
#  http://www.rqtl.org												   #
#																	   #		
#  Authors: Karl W Broman <kbroman@biostat.wisc.edu> and Hao Wu, with  #
#         ideas from Gary Churchill and Saunak Sen and contributions   #
#         from Danny Arends, Ritsert Jansen, Pjotr Prins, Laura        #
#         Shannon, and Brian Yandell								   #
#														               #
#  The R/qtl package is free software; you can redistribute it and/or  #
#  modify it under the terms of the GNU General Public License,        #
#  version 3, as published by the Free Software Foundation.		       #
#----------------------------------------------------------------------#

probbc <- function (cross,chrom_num, step = 0, off.end = 0, error.prob = 1e-04, map.function = c("haldane", 
    "kosambi", "c-f", "morgan"), stepwidth = c("fixed", "variable", 
    "max")) 
{
    if (!any(class(cross) == "cross")) 
        stop("Input should have class \"cross\".")
    map.function <- match.arg(map.function)
    if (map.function == "kosambi") 
        mf <- mf.k
    else if (map.function == "c-f") 
        mf <- mf.cf
    else if (map.function == "morgan") 
        mf <- mf.m
    else mf <- mf.h
    stepwidth <- match.arg(stepwidth)
    if (error.prob < 1e-50) 
        error.prob <- 1e-50
    if (error.prob > 1) {
        error.prob <- 1 - 1e-50
        warning("error.prob shouldn't be > 1!")
    }
    n.ind <- nind(cross)
    n.chr <- nchr(cross)
    n.mar <- nmar(cross)
    type <- class(cross)[1]
    #for (i in 1:n.chr) {
	for (i in 1:1) {
		i <- chrom_num;
        if (n.mar[i] == 1) 
            temp.offend <- max(c(off.end, 5))
        else temp.offend <- off.end
        chrtype <- class(cross$geno[[i]])
        if (chrtype == "X") 
            xchr <- TRUE
        else xchr <- FALSE
        if (type == "f2") {
            one.map <- TRUE
            if (!xchr) {
                cfunc <- "calc_genoprob_f2"
                n.gen <- 3
                gen.names <- getgenonames("f2", "A", cross.attr = attributes(cross))
            }
            else {
                cfunc <- "calc_genoprob_bc"
                n.gen <- 2
                gen.names <- c("g1", "g2")
            }
        }
        else if (type == "bc") {
            cfunc <- "calc_genoprob_bc"
            n.gen <- 2
            if (!xchr) 
                gen.names <- getgenonames("bc", "A", cross.attr = attributes(cross))
            else gen.names <- c("g1", "g2")
            one.map <- TRUE
        }
        else if (type == "riself" || type == "risib" || type == 
            "dh" || type == "haploid") {
            cfunc <- "calc_genoprob_bc"
            n.gen <- 2
            gen.names <- getgenonames(type, "A", cross.attr = attributes(cross))
            one.map <- TRUE
        }
        else if (type == "4way") {
            cfunc <- "calc_genoprob_4way"
            n.gen <- 4
            one.map <- FALSE
            gen.names <- getgenonames(type, "A", cross.attr = attributes(cross))
        }
        else if (type == "ri8sib" || type == "ri4sib" || type == 
            "ri8self" || type == "ri4self" || type == "bgmagic16") {
            cfunc <- paste("calc_genoprob_", type, sep = "")
            if (type == "bgmagic16") 
                n.gen <- 16
            else n.gen <- as.numeric(substr(type, 3, 3))
            one.map <- TRUE
            gen.names <- LETTERS[1:n.gen]
            if (xchr) 
                warning("calc.genoprob not working properly for the X chromosome for 4- or 8-way RIL.")
        }
        else if (type == "bcsft") {
            one.map <- TRUE
            cfunc <- "calc_genoprob_bcsft"
            cross.scheme <- attr(cross, "scheme")
            if (!xchr) {
                gen.names <- getgenonames("bcsft", "A", cross.attr = attributes(cross))
                n.gen <- 2 + (cross.scheme[2] > 0)
            }
            else {
                cross.scheme[1] <- cross.scheme[1] + cross.scheme[2] - 
                  (cross.scheme[1] == 0)
                cross.scheme[2] <- 0
                gen.names <- c("g1", "g2")
                n.gen <- 2
            }
        }
        else stop("calc.genoprob not available for cross type ", 
            type, ".")
        gen <- cross$geno[[i]]$data
        gen[is.na(gen)] <- 0
        if (one.map) {
            map <- create.map(cross$geno[[i]]$map, step, temp.offend, 
                stepwidth)
            rf <- mf(diff(map))
            if (type == "risib" || type == "riself") 
                rf <- adjust.rf.ri(rf, substr(type, 3, nchar(type)), 
                  chrtype)
            rf[rf < 1e-14] <- 1e-14
            newgen <- matrix(ncol = length(map), nrow = nrow(gen))
            dimnames(newgen) <- list(NULL, names(map))
            newgen[, colnames(gen)] <- gen
            newgen[is.na(newgen)] <- 0
            n.pos <- ncol(newgen)
            marnames <- names(map)
        }
        else {
            map <- create.map(cross$geno[[i]]$map, step, temp.offend, 
                stepwidth)
            rf <- mf(diff(map[1, ]))
            rf[rf < 1e-14] <- 1e-14
            rf2 <- mf(diff(map[2, ]))
            rf2[rf2 < 1e-14] <- 1e-14
            newgen <- matrix(ncol = ncol(map), nrow = nrow(gen))
            dimnames(newgen) <- list(NULL, dimnames(map)[[2]])
            newgen[, colnames(gen)] <- gen
            newgen[is.na(newgen)] <- 0
            n.pos <- ncol(newgen)
            marnames <- colnames(map)
        }
        if (one.map) {
            temp <- as.double(rep(0, n.gen * n.ind * n.pos))
            if (type == "bcsft") 
                temp[1:2] <- cross.scheme
            z <- .C(cfunc, as.integer(n.ind), as.integer(n.pos), 
                as.integer(newgen), as.double(rf), as.double(error.prob), 
                genoprob = as.double(temp), PACKAGE = "qtl")
        }
        else {
            z <- .C(cfunc, as.integer(n.ind), as.integer(n.pos), 
                as.integer(newgen), as.double(rf), as.double(rf2), 
                as.double(error.prob), genoprob = as.double(rep(0, 
                  n.gen * n.ind * n.pos)), PACKAGE = "qtl")
        }
        cross$geno[[i]]$prob <- array(z$genoprob, dim = c(n.ind, 
            n.pos, n.gen))
        dimnames(cross$geno[[i]]$prob) <- list(NULL, marnames, 
            gen.names)
        attr(cross$geno[[i]]$prob, "map") <- map
        attr(cross$geno[[i]]$prob, "error.prob") <- error.prob
        attr(cross$geno[[i]]$prob, "step") <- step
        attr(cross$geno[[i]]$prob, "off.end") <- temp.offend
        attr(cross$geno[[i]]$prob, "map.function") <- map.function
        attr(cross$geno[[i]]$prob, "stepwidth") <- stepwidth
    }
    if (type == "ri4self" || type == "ri4sib" || type == "ri8self" || 
        type == "ri8sib" || type == "bgmagic16") 
        cross <- reorgRIgenoprob(cross)
    #cross
	cross$geno[[chrom_num]][[3]]
}
##################################################################################################################