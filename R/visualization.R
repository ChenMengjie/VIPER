PlotWeights <- function(Imputed, KnownLabel = NULL){

  suppressMessages(require(gplots))
  weights.mat <- apply(Imputed$sample_weights, 2, function(x){
    y <- x
    y[y>0.01] <- 1
    y[y<=0.01] <- 0
    y
  })

  palette.gr.marray <- colorRampPalette(c("white", "red"))(4)

  if(!is.null(KnownLabel)){

    colorlist <- c("turquoise4", "cyan", "lavender",  "gold", "slateblue1", "violet", "skyblue1", "khaki", "pink", "salmon", "limegreen", "chocolate", "maroon", "purple", "blue", "yellow", "red",  "brown", '#FFAAD4', '#00CC00', '#66B2FF', '#B266FF', '#FFD4AA', '#AAFFD4', '#CC0000', "#B266CC")
    col.color <- colorlist[as.numeric(as.factor(KnownLabel))]
    heatmap.2(as.matrix(weights.mat), trace = "none", col = palette.gr.marray, symbreaks = F,
              labCol = NA, dendrogram = "none", ColSideColors = col.color, labRow = NA, key = F, Colv = FALSE, Rowv = FALSE)

  } else {
    heatmap.2(as.matrix(weights.mat), trace = "none", col = palette.gr.marray, symbreaks = F,
              labCol = NA, dendrogram = "none", labRow = NA, key = F, Colv = FALSE, Rowv = FALSE)


  }

}

PlotGeneMatrix <- function(Imputed, GeneExpression = NULL, KnownLabel = NULL, GeneNum = 2000){

  suppressMessages(require(grid))

  require(gplots)
  require(gridExtra)

  Imputed_log <- Imputed$imputed_log

  p <- nrow(Imputed$imputed_log)
  if(p > GeneNum){
    selected <- round(runif(2500)*p)
    fll <- c(1:p)[selected]
  } else {
    fll <- c(1:p)
  }
  palette.gr.marray <- colorRampPalette(c("blue", "white", "red"))(56)

  if(!is.null(GeneExpression)){
    if(!is.null(KnownLabel)){

      colorlist <- c("turquoise4", "cyan", "lavender",  "gold", "slateblue1", "violet", "skyblue1", "khaki", "pink", "salmon", "limegreen", "chocolate", "maroon", "purple", "blue", "yellow", "red",  "brown", '#FFAAD4', '#00CC00', '#66B2FF', '#B266FF', '#FFD4AA', '#AAFFD4', '#CC0000', "#B266CC")
      col.color <- colorlist[as.numeric(as.factor(KnownLabel))]
      aa <- heatmap.2(as.matrix(GeneExpression[fll, ]), trace = "none", col = palette.gr.marray,
                  symbreaks = T, labCol = NA, dendrogram = "none", ColSideColors = col.color,
                  labRow = NA, key = T, Colv = FALSE, Rowv = TRUE)
      bb <- rev(aa$rowInd)
      cc <- heatmap.2(as.matrix(Imputed_log[fll, ][bb, ]), trace = "none", col = palette.gr.marray,
                symbreaks = T, labCol = NA, dendrogram = "none", ColSideColors = col.color,
                labRow = NA, key = T, Colv = FALSE, Rowv = FALSE)

    } else {
      aa <- heatmap.2(as.matrix(GeneExpression[fll, ]), trace = "none", col = palette.gr.marray,
                      symbreaks = T, labCol = NA, dendrogram = "none", labRow = NA, key = T, Colv = FALSE, Rowv = TRUE)
      bb <- rev(aa$rowInd)
      cc <- heatmap.2(as.matrix(Imputed_log[fll, ][bb, ]), trace = "none", col = palette.gr.marray,
                symbreaks = T, labCol = NA, dendrogram = "none", labRow = NA, key = T, Colv = FALSE, Rowv = FALSE)
      gl <- list(aa, cc)

      grab_grob <- function(){
        grid.grab()
      }

    }
  } else {

    if(!is.null(KnownLabel)){

      colorlist <- c("turquoise4", "cyan", "lavender",  "gold", "slateblue1", "violet", "skyblue1", "khaki", "pink", "salmon", "limegreen", "chocolate", "maroon", "purple", "blue", "yellow", "red",  "brown", '#FFAAD4', '#00CC00', '#66B2FF', '#B266FF', '#FFD4AA', '#AAFFD4', '#CC0000', "#B266CC")
      col.color <- colorlist[as.numeric(as.factor(KnownLabel))]
      eatmap.2(as.matrix(Imputed_log[fll, ]), trace = "none", col = palette.gr.marray,
                symbreaks = T, labCol = NA, dendrogram = "none", ColSideColors = col.color,
                labRow = NA, key = T, Colv = FALSE, Rowv = TRUE)
    } else {
      heatmap.2(as.matrix(Imputed_log[fll, ]), trace = "none", col = palette.gr.marray,
                symbreaks = T, labCol = NA, dendrogram = "none", labRow = NA, key = T, Colv = FALSE, Rowv = TRUE)

    }

  }

}

PlotCV <- function(Imputed, GeneExpression, KnownLabel = NULL, GeneNum = 2000){

  require(gplots)

  logxx <- apply(gene.expression, 2, function(y){log(y + 1)})
  logxx[gene.expression==0] <- 0

  Imputed_log <- Imputed$imputed_log

  CV.all <- function(data){
    CV <- function(mean, sd){
      abs(sd/mean)
    }
    CV.per.gene <- apply(data, 1, function(x){
      CV(mean(x), sd(x))
    })
    CV.per.gene
  }

  CV.nonzero <- function(data){
    CV <- function(mean, sd){
      abs(sd/mean)
    }
    CV.per.gene <- apply(data, 1, function(x){
      CV(mean(x[x!=0]), sd(x[x!=0]))
    })
    CV.per.gene
  }

  types <- unique(KnownLabel)
  len <- length(types)
  show.some <- sample(1:nrow(GeneExpression), GeneNum)

  for(i in 1:len){

    flag <- which(KnownLabel%in%types[i])
    raw.data <- logxx[, flag]
    zero.num <- apply(raw.data, 1, function(x){
      length(x[x==0])
    })
    zero.rate <- round(zero.num/ncol(logxx), 2)*100
    dropout.rate <- apply(raw.data, 1, function(x){
      round(mean(x[x!=0]), 2)
    })

    SIMPLE.CV <- CV.all(Imputed_log[, flag])[show.some]
    raw.cv <- CV.nonzero(logxx[, flag])[show.some]
    zero.rate.selected <- zero.rate[show.some]
    dropout.rate.selected <- dropout.rate[show.some]

    df <- data.frame(value = SIMPLE.CV, compr = rep("Without-imputation", GeneNum), raw = raw.cv, zero = zero.rate.selected,
                     dropout = dropout.rate.selected)

    gg1 <- ggplot2.scatterplot(data = df, xName = 'raw', yName = 'value', size = 0.2,
              backgroundColor = "white", xtitle="CV (Before Imputation)", ytitle="CV (After Imputation)",
               mainTitle = types[i], removePanelGrid=TRUE, removePanelBorder=FALSE, showLegend=TRUE,
              legendTitle = "Percentage \n of Zero", legendTitleFont = c(15, "bold", "black"),
              legendTextFont = c(15, "bold", "black"), mainTitleFont = c(10, "bold", "black"),
              xtitleFont = c(15, "bold", "black"),  ytitleFont = c(15, "bold", "black"),
              xTickLabelFont = c(15, "bold", "white"), yTickLabelFont = c(15, "bold", "black")) + geom_point(aes(colour = zero), size = 0.2) + ylim(0, 5) + xlim(0, 5) + theme(strip.text.x = element_text(size = 15,
              colour = "black", face = "bold")) + geom_abline(col = "brown", linetype = "dashed", size=1)

    gg2 <- ggplot2.scatterplot(data = df, xName = 'raw', yName = 'value', size = 0.2,
           backgroundColor = "white", xtitle="CV (Before Imputation)", ytitle="CV (After Imputation)",
           mainTitle = types[i], removePanelGrid=TRUE, removePanelBorder=FALSE, showLegend=TRUE,
           legendTitle = "Mean of nonzero \n values", legendTitleFont = c(15, "bold", "black"),
           legendTextFont = c(15, "bold", "black"), mainTitleFont = c(10, "bold", "black"),
           xtitleFont = c(15, "bold", "black"),  ytitleFont = c(15, "bold", "black"),
            xTickLabelFont = c(15, "bold", "white"), yTickLabelFont = c(15, "bold", "black")) + geom_point(aes(colour = dropout), size = 0.2) + ylim(0, 5) + xlim(0, 5)  + theme(strip.text.x = element_text(size = 15,colour = "black", face = "bold")) + geom_abline(col = "brown", linetype = "dashed", size=1)
    multiplot(gg1, gg2, cols=2)

  }


}




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
  library(grid)

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
