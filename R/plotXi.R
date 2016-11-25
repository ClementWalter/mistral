plotXi <- function(krig, q, xi = getXi(krig), xplot = c(-80:80)/10, yplot = c(-80:80)/10, lsf, plot.lab = c('x', 'y'), plot.learndb = TRUE) {

  # Fix NOTE issue with R CMD check
  x <- y <- z <- ..level.. <- crit <- NULL
  
  df_plot = data.frame(expand.grid(x=xplot, y=yplot), z = NA)
  if(!missing(lsf)){
    df_plot$z = lsf(t(df_plot[,1:2]))
  }
  learn_db <- krig@X
  colnames(learn_db) <- c('x', 'y')
  p <- ggplot2::ggplot(data = df_plot, aes(x,y)) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlim(-8, 8) + ggplot2::ylim(-8, 8) +
    ggplot2::xlab(plot.lab[1]) + ggplot2::ylab(plot.lab[2])
  if(plot.learndb==TRUE) {
    p <- p + ggplot2::geom_point(data = data.frame(learn_db, z = krig@y), aes(color=z))
  }
  if(!missing(lsf)){
    p <- p + ggplot2::geom_contour(aes(z=z), color="black", alpha = 0.5, breaks = q)
  }
  z_meta = xi(t(df_plot[,1:2]))
  df_plot_meta <- data.frame(df_plot[,1:2], z = z_meta$mean, crit = abs(q - z_meta$mean)/z_meta$sd)
  p_meta <- p + ggplot2::geom_contour(data = df_plot_meta, aes(z=z, color=..level.., alpha = 0.5), breaks = q) +
          ggplot2::geom_contour(data = df_plot_meta, aes(z=crit), breaks = 2, linetype = 'dotted')
  return(p_meta)
}
