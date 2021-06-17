#adapted from R-package "visStatistics"

openGraph = function(width = 7,
                     height = 7/2 ,
                     mag = 1 , family="Times New Roman",
                     ...) {
  if (.Platform$OS.type != "windows") {
    # Mac OS, Linux
    # X11( width<-width*mag , height<-height*mag , type="cairo" , ... ) #Unix style
    quartz(width = width * mag , height = height * mag  , ...)
  } else {
    # Windows OS
    windows(width = width * mag , height = height * mag , ...)
  }
}

#short function to save graphics displayed either in quartz(for mac) or windows
saveGraph = function(file = "saveGraphOutput" , type = "pdf" , ...) {
  #Set default values
  if (missing(type))
  {
    type = "pdf"
  }


  if (.Platform$OS.type != "windows") {
    # Mac OS, Linux
    if (any(type == c("png", "jpeg", "jpg", "tiff", "bmp"))) {
      sptype = type
      if (type == "jpg") {
        sptype = "jpeg"
      }
      # savePlot(file = paste(file, ".", type, sep = "") ,
      #         type = sptype ,
      #          ...)
      #dev.copy(sptype,file = paste(file, ".", sptype, sep = ""))
      #functions only on windows/macos with quartz
      quartz.save(file = paste(file, ".", type, sep = ""),
                  type = sptype)

    }
    if (type == "pdf") {
      dev.copy2pdf(file = paste(file, ".", type, sep = "") , ...)
    }
    if (type == "eps") {
      dev.copy2eps(file = paste(file, ".", type, sep = "") , ...)
    }
    if (type == "svg") {
      dev.copy(svg, file = paste(file, ".", type, sep = ""))
      dev.off()
    }

  } else {
    # Windows OS
    file = paste(file, ".", type, sep = "") # force explicit extension
    savePlot(file = file , type = type , ...)
  }
}


### helper for Kruskal Wallis and post-hoc Wilcox-----
sig_diffs_nongauss <- function(samples, fact)
{
  # function to produce a table similar to that produced for TukeyHSD
  # but for non-normally distributed data

  # calculate p values for each data classification

  ufactor = levels(fact)
  pwt = pairwise.wilcox.test(samples, fact)
  factormeans = matrix(0, length(ufactor), 1)
  for (ii in 1:length(ufactor)) {
    pos = which(fact == ufactor[ii])

    factormeans[ii] = mean(samples[pos])

  }

  # make a matrix with a row for every possible combination of
  # 2 data classifications and populate it with the calculated
  # p values

  xcomb = combn(length(ufactor), 2)
  tukeylike = matrix(0, ncol(xcomb), 4)
  colnames(tukeylike) <- c("diff", "lwr", "upr", "p adj")
  tukeynames = vector("list", ncol(xcomb))
  for (ii in 1:ncol(xcomb)) {
    tukeynames[ii] =
      paste(ufactor[xcomb[2, ii]], "-", ufactor[xcomb[1, ii]], sep = "")


    p_value = pwt$p.value[xcomb[2, ii] - 1, xcomb[1, ii]]

    if (is.na(p_value)) {
      p_value = 1
    }
    tukeylike[ii, 4] = p_value
    tukeylike[ii, 1] = 0

    tukeylike[ii, 2] = 0

    tukeylike[ii, 3] = 0

  }
  rownames(tukeylike) = tukeynames

  # re-format the table slightly so it is the same as that produced
  # by TukeyHSD and output

  tukeylike2 = list(tukeylike)
  #print(tukeylike2)
  return(tukeylike2)
}




###### Visualize Kruskal_Wallis ###############################
## performs Kruskal Wallis and post-hoc Wilcoxon:

vis_Kruskal_Wallis_clusters = function(samples,
                                       fact,
                                       alpha = 0.05,
                                       xlab = "",
                                       ylab = "samples",
                                       cex = cexsize,
                                       notch = F,
                                       samplename = "",
                                       lwd=lwdsize) {
  if (missing(lwd))
  {
    lwd=1
  }

  if (missing(alpha))
  {
    alpha = 0.05
  }

  if (missing(xlab))
  {
    xlab = ""
  }
  if (missing(ylab))
  {
    ylab = "samples"
  }

  if (missing(cex))
  {
    cex = 1
  }

  #remove rows with NAs in samples
  samples3 = na.omit(samples)
  fact <- subset(fact,!is.na(samples))
  samples = samples3


  n_classes = length(unique(fact))

  mc=c("red","blue","darkgreen")
  s = tapply(samples, fact, sd)
  m = tapply(samples, fact, mean)

  samples_per_class = c()
  for (i in 1:n_classes) {
    samples_per_class[i] = sum(fact == unique(fact)[i])
  }

  kk = kruskal.test(samples ~ fact)

  maximum = max(samples, na.rm = T)

  minimum = min(samples, na.rm = T)

  sp = maximum - minimum

  mi = minimum - 0.1 * sp
  ma = maximum + 0.1 * sp



  openGraph()
  if (notch == TRUE) {
    b = boxplot(
      samples ~ fact,
      notch = TRUE,
      col = mc,
      las = 1,
      # xlim = c(0, n_classes + 0.04),
      ylim = c(mi, ma),
      # xlab = xlab,
      #ylab = ylab,
      boxwex=.4,

      #changes group names size
      cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex,
      xaxt="n" #no label on x axis

    )
  }
  else
  {
    b = boxplot(
      samples ~ fact,
      notch = FALSE,
      col = mc,
      las = 1,
      lwd=cexsize,
      #xlim = c(0, n_classes +0.05),
      ylim = c(mi, ma),
      xlab = "",
      ylab = ylab,
      #boxwex = 0.4,
      xaxt="n"
    )
  }


  stripchart(
    samples ~ fact,
    vertical = TRUE,
    jitter=0.2,  #decides the amount of jitter
    method="jitter",
    col = rep(mc, n_classes),
    ylab = ylab,
    xlab = "",
    las = 1,
    #horizontal legend,
    add = TRUE,
    xaxt="n"
  )

  #mtext(c("N = ", b$n), at = c(0.7, seq(1, n_classes)), las = 1) #nmber of cases in each group
  tuk = sig_diffs_nongauss(samples, fact)
  library(multcompView)
  s = multcompLetters(tuk[[1]][, 4], threshold = alpha)

  ord = c()

  v = attributes(s$Letters)$names
  f_levels = sort(unique(fact))
  for (i in 1:n_classes) {
    ord[i] = which(v == f_levels[i])
  }
  (ma)
  text(
    seq(1:n_classes + 1),
    mi,
    s$Letters[ord],
    col = "black",
    cex = cexsize,
    lwd = 2
  )

  #title(paste("Kruskal Wallis: P = ", signif(kk$p.value, digits = 4)), outer = TRUE)
  my_list <-
    list("Kruskal_wallis" = kk,
         "adjusted_p_values_wilcoxon" = tuk)
  return(my_list)


}

## Define function to project into First Quadrant ------
projectAngleToFirstQuadrant=function (angle)
{
  #project second and third quadrant to first quadrant
  angle[angle>90 &angle<=270] <-abs(180-angle[angle>90&angle<=270])

  #project fourth quadrant to first quadrant
  angle[angle>270 &angle<=360] <-360-angle[angle>270 &angle<=360]

  #project negative angle
  return(angle)
}
