MYappend2tiplabel<-function (phy, tips, offset, align = FALSE, grid = FALSE, col = "red", 
          bg="black", text = NULL, pch = NULL, cex, ...) 
{
 #finding last tree plotted
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  #checking that it is rightwards or leftwards
  if (any(c("upwards", "downwards") %in% lastPP$direction)) 
    stop("direction '", lastPP$direction, "' is not supported by 'append2tiplabel'. ", 
         "\nYou might try and contact the package maintainer", 
         sep = "")
  #setting default values if parameters are't specified in the command
  if (missing(tips)) 
    tips <- 1:lastPP$Ntip
  if (missing(offset)) 
    offset <- lastPP$label.offset
  if (missing(cex)) 
    cex <- lastPP$cex
  #setting one spacing value to indicate where the current tip labels end
  space1 <- lastPP$label.offset + strwidth(gsub("_", " ", phy$tip.label[tips]), 
                                           "user", cex, lastPP$font)
  #setting another spacing value: dependent on the point size or text length
  space2 <- rep(0, length(tips))
  if (!is.null(pch)) 
    space2 <- strwidth("O", cex = cex)
  if (!is.null(text)) 
    space2 <- strwidth(text, cex = cex)
  #setting a third spacing value to put points an even distance after current tip labels?
  space3 <- strwidth("-", "user", cex, lastPP$font)
 #setting some variables to determine... the ??order of tip labels? not sure
  if (any(c("rightwards", "upwards") %in% lastPP$direction)) {
    dir <- 1
    adj1 <- 1
    minmax <- 2
    gord <- 1:2
  }
  else {
    dir <- -1
    adj1 <- 0
    minmax <- 1
    gord <- 2:1
  }
  #setting some more variables...
  if (any(c("rightwards", "leftwards") %in% lastPP$direction)) {
    #changed ord from 1:2, to 1:3, to include y coords as well
    ord <- 1:3
    srt <- 0
    ratio <- 1
    #if you don't have any variable "new.xx" in the current tree,
    #call "id" "xx" - otherwise call it 'new.xx" (??)
    id <- ifelse(is.null(lastPP$new.xx), "xx", "new.xx")
    new.list.element <- "new.xx"
  }
  else {
    ord <- 2:1
    srt <- 90
    ratio <- (diff(lastPP$y.lim) + space1 + space2 + offset)/diff(lastPP$x.lim)
    #if you don't have any variable "new.yy" in the current tree,
    #call "id" "yy" - otherwise call it 'new.yy" (??)
    id <- ifelse(is.null(lastPP$new.yy), "yy", "new.yy")
    new.list.element <- "new.yy"
  }
  #if you don't have any variables called "new"-something in the variable "id"
  #but you are still showing the tip labels in the current tree, 
  #make fac1 equal to 1: otherwise 0
  fac1 <- ifelse(length(grep("new", id)) == 0 & lastPP$show.tip.label, 
                 1, 0)
#in a rightwards tree, this then selects "xx" in lastPP (the x-coords of the tips)
  xy <- lastPP[[id]]
  #then finds the x-coordinates of the tips you've specified
  #and adds direction, spaces etc.
  #making fac1 0 would position the symbols right next to the tips?
  xy[tips] <- xy[tips] + dir * ratio * (fac1 * space1 + (space2 + 
    offset))
  #making a new variable for this x-coordinate (adjusted)
  cds <- xy[tips]
  #if you've specified the "align" option in the original command,
  #line up all the x-coordinates
  if (align) {
    #finding the max. x-coordinate
    extcds <- range(cds)[minmax]
    #making this x-coordinate the length of cds
    extcds <- rep(extcds, length(cds))
    line.cds <- cbind(cds - dir * (space2 + offset - space3), 
                      extcds - dir * (space2 + space3))
    cds <- extcds
  }
  #otherwise don't align them
  else {
    line.cds <- cbind(cds - dir * (space2 + offset), cds - 
      dir * (space2 + space3))
  }
  #now cbind with the tip numbers 
  #(for rightwards tree, ordered as coord, tip)
  #this originally only included the x locations:
  #cds <- cbind(cds, tips)[, ord]
  
  #not sure if this is right, but including the y locations too.
  cds<-cbind(lastPP$xx[tips]+(0.5*offset),lastPP$yy[tips], tips)[,ord]
  #if you've specified the "grid" option in the original command,
  #will set it up like that...
  if (grid) {
    gid <- which(apply(line.cds[, gord], 1, diff) > 0)
    for (i in gid) lines(line.cds[i, ], rep(tips[i], 2), 
                         lty = 3, cex = cex/2)
  }
  #if you specified the point type, draw points as follows
  if (!is.null(pch)) {
    #temporarily changing x-coord
    #also changeded this bit for the case where there's only 1 point per category
    if (is.vector(cds)==TRUE) {
    points(cds[1], cds[2], pch = pch, col = col, bg = bg, 
         cex = cex, ...)
    }
    else {
      points(cds[, 1], cds[, 2], pch = pch, col = col, bg = bg, 
             cex = cex, ...)
    }
  }
  #if you specified the text type, draw points as follows
  if (!is.null(text)) {
    text(cds[, 1], cds[, 2], text, col = col, cex = cex, 
         adj = c(adj1, 0.5), srt = srt, ...)
  }
  if (max(xy) > max(lastPP$x.lim)) 
    warning("The appended text exceeds the limit of the plotting region.", 
            "\n  Consider setting 'x.lim' in 'plot.phylo' to at least ", 
            round(max(xy), 5), ".")
#This code assigns the points plotted as 'lastPP'
  #but I don't really want that...
  #lastPP[[new.list.element]] <- xy
  #assign("last_plot.phylo", lastPP, envir = .PlotPhyloEnv)
  #invisible(lastPP)
}