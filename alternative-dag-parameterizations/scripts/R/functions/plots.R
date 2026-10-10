imap_plot_diamond <- function() {
  old <- par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2.5, 0.5))
  on.exit(par(old))
  positions <- list(rbind(c(.5,.9),c(.15,.5),c(.85,.5),c(.5,.1)),
                    rbind(c(.1,.85),c(.45,.5),c(.45,.1),c(.85,.85)))
  edges <- list(rbind(c(1,2),c(1,3),c(2,4),c(3,4)),
                rbind(c(1,4),c(1,2),c(4,2),c(1,3),c(4,3),c(2,3)))
  titles <- c('Original causal DAG', 'Ordered I-MAP')
  for (k in 1:2) {
    xy <- positions[[k]]; plot.new(); plot.window(c(0,1), c(0,1), asp = 1)
    title(titles[k], cex.main = 1)
    for (r in seq_len(nrow(edges[[k]]))) {
      u <- edges[[k]][r,1]; v <- edges[[k]][r,2]
      direction <- xy[v,]-xy[u,]; direction <- direction/sqrt(sum(direction^2))
      start <- xy[u,]+.065*direction; end <- xy[v,]-.065*direction
      color <- if (k == 2 && u == 1 && v == 4) 'steelblue4' else 'gray35'
      arrows(start[1],start[2],end[1],end[2],length=.09,col=color,lwd=1.5)
    }
    points(xy,pch=21,bg='white',cex=3.5,lwd=1.3)
    text(xy,labels=1:4,cex=.9)
    if (k == 2) text(.47,.92,expression(tau[14]),col='steelblue4')
  }
}
