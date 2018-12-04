### function for counting number of cycles of all possible sizes
count_cycles <- function(net, nmax = length(V(net))){ # default nmax = size of the network
  #net<-net1
  ### create adjacency matrix
  lm <- laplacian_matrix(net) # matrix representation
  A <- as.matrix(lm)
  A[which(A >= 0)] <- 0
  A[which(A < 0)] <- 1
  ### create table for cycle counts 
  ncyc <- matrix(0, nmax/2, 2)
  ncyc[, 1] <- seq(2, nmax, by=2)
  i <- 1
  k <- 1
  while (i <= nmax){  # while we haven't exceeded number of nodes
    if (i %% 2 == 0){     # for all even numbers
      Ai <- A
      for (j in 2:i){ # calc A^i
        Ai <- Ai %*% A
      }
      t <- tr(Ai)       # calc trace(A^i)
      sum <- 0
      for (j in 1:k){
        if (i %% (2 * j) == 0){    # find all multiple cycles of smaller length
          sum <- sum + ncyc[j, 2] * 2 * j
        }
      }
      ncyc[k, 2] <- (t - sum) / i   # subtract number of diagonal elements from all smaller multiple cycles 
      if (ncyc[k, 2] > 0){      # decrease nmax if cycle is found (each node could be a part of only one cycle) 
        nmax <- nmax - ncyc[k, 2] * i
      }
      k <- k + 1
    }
    i <- i + 1
  }
  colnames(ncyc) <- c("length", "number")
  return(ncyc)
}

### fuction for plotting single network
plot_network <- function(links, nodes, layout = layout.bipartite, vertex.size = 20){ 
  ### create net
  net <- graph_from_data_frame(d = links, vertices = nodes, directed = T) 
  net <- simplify(net, remove.multiple = F, remove.loops = F) 
  
  ### set different colors for Carbon and Nitrogen
  nodes <- cbind(nodes,lapply(nodes, function(x) strsplit(x, split = "")[[1]][1]))
  V(net)$type <- nodes[,2] %in% "N" 
  col <- c("red", "royal blue")
  
  ### try to display concentrations and abundances...
  #sizes <- round((log10(as.numeric(nodes[,2]))-min(log10(as.numeric(nodes[,2])))+1), digits = 1)*10
  #size <- as.numeric(nodes[,2])/max(as.numeric(nodes[,2]))
  #V(net)$size <- size
  #E(net)$weight <- as.numeric(links$V3)
  #V(net)$color <- alpha(V(net)$color, V(net)$size/max(V(net)$size))
  #E(net)$color = alpha("black", E(net)$weight/max(E(net)$weight))
  
  V(net)$color <- col[as.numeric(V(net)$type) + 1]
  E(net)$color <- "black"
  ### plot network
  plot_net <- plot(net, vertex.label.cex = 1, edge.curved = curve_multiple(net), 
                 vertex.label.color = "black", layout = layout, 
                 vertex.size = vertex.size, 
                 #vertex.size = V(net)$size*2,
                 edge.color = E(net)$color,
                 edge.width = 1,#E(net)$width,
                 edge.arrow.size = 1,
                 vertex.color = V(net)$color)
  return(plot_net)
}

### function for comparison of two networks
plot_joint_network<-function(merged, nodes, st, layout = layout.bipartite, vertex.size = 20, start=0.4){
  
  ### create networks
  #nodes<-unique(as.vector(as.matrix(merged[,1:2])))
  net <- graph_from_data_frame(d=merged, vertices=nodes, directed=T) 
  
  ### set different colors for Carbon and Nitrogen
  nodes<-cbind(nodes,lapply(nodes, function(x) strsplit(x, split = "")[[1]][1]))
  V(net)$type <- nodes[,2] %in% "N" 
  vcol <- c("red", "royal blue") # gold, etc
  
  ### set defferent colors for bacteria in different states
  ecolors <- c("black","firebrick","olivedrab3","gold","magenta")
  ecol <- rep(ecolors[1], ecount(net))
  for (j in 2:length(unique(E(net)$state))){
    ecol[which(E(net)$state==j)] <- ecolors[j]
  }
  
  plot(net, edge.arrow.size=1, vertex.label.cex=1, 
       #edge.curved=seq(-0.5, 0.5, length = ecount(net)), 
       edge.curved=autocurve.edges2(net, start=start),
       vertex.label.color="black", #layout=layout.bipartite, 
       vertex.size=vertex.size, edge.color=ecol, 
       edge.width=2, edge.arrow.size=1,
       asp = 1,
       vertex.color = vcol[as.numeric(V(net)$type) + 1])
  legend("topleft", title = "Volume of states",legend = t(st), col = ecolors[1:length(t(st))], 
         lty= 1, lwd = 2, bty="n")
}

### fixed autocurve function for curving multiple edges on graph
autocurve.edges2 <-function (graph, start = 0.5)
{
  cm <- count.multiple(graph)
  mut <-is.mutual(graph)  #are connections mutual?
  el <- apply(get.edgelist(graph, names = FALSE), 1, paste,
              collapse = ":")
  ord <- order(el)
  res <- numeric(length(ord))
  p <- 1
  while (p <= length(res)) {
    m <- cm[ord[p]]
    mut.obs <-mut[ord[p]] #are the connections mutual for this point?
    idx <- p:(p + m - 1)
    if (m == 1 & mut.obs==FALSE) { #no mutual conn = no curve
      r <- 0
    }
    else {
      r <- seq(-start, start, length = m)
    }
    res[ord[idx]] <- r
    p <- p + m
  }
  res
}

### function for nice legend
legend.col <- function(col, lev){
  opar <- par
  n <- length(col)
  bx <- par("usr")
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  
  xx <- rep(box.cx, each = 2)
  
  par(xpd = TRUE)
  for(i in 1:n){
    
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
    
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par <- opar
}