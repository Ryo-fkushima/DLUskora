dluskora1 <- function(fac1, Q, syR, c_ave, D_0, R, T_1, T_2, K_d, Mr, Mt, fg, ft, garsize){
  #-------------------------------------------------------------------------------
  A <- D_0 / fac1
  gar_g <- function(t){ #growth speed of gar (linear)
    return(A * t^(fg))
  }
  time <- (garsize / A)^(1 / fg) #year
  temp <- function(x){
    return((T_2 - T_1) / (time)^ft * x^ft + T_1 + 273.15)
  }
  D <- function(x){
    return(D_0 * exp(-Q / R / temp(x)))
  }
  diftype <- 3 #(1)eq/non-dif (2)eq/fixed-c/flowing (3)diseq/no-flow(appropriate)
  #-------------------------------------------------------------------------------
  delta_r <- syR / Mr
  delta_t <- time / Mt
  rdots <- seq(delta_r / 2, (Mr -delta_r / 2), by = delta_r)
  cter <- numeric(Mt)

  K <- function(x, y){
    return(D(y * delta_t) * delta_t / 2 / x / delta_r / delta_r)
  }

  #-------------------------------------------------------------------------------
  c <- matrix(ncol = Mr, nrow = Mt, 0) #c matrix created

  #initial c defined
  for(i in 1:Mr){
    c[1, i] <- c_ave
  }


  #-------------------------------------------------------------------------------
  mstool <- matrix(ncol = Mr, nrow = Mr, 0)
  avemas <- numeric(Mr)
  for(i in 1:Mr){
    mstool[i,i] <- i * i
    avemas[i] <- i * i
  }
  avemasv <- c_ave * sum(avemas)
  eps <- avemasv * 0.000001
  escape <- 0
  #-------------------------------------------------------------------------------
  for(j in 1:(Mt - 1)){
    cter[j + 1] <- length(rdots[rdots <= gar_g(j * delta_t)])
    cwork <- c[j,]
    cworkm <- matrix(cwork, ncol = 1)

    #Bnew defined
    B <- matrix(ncol = (Mr + 2), nrow = Mr, 0)
    for(i in 1:Mr){
      for(l in 2:(Mr + 1)){
        if(l == (i + 1)){
          B[i,(l - 1)] <- -K(i, j) * (l - 2)
          B[i,l] <- 1 + 2 * K(i, j) * (l - 1)
          B[i,(l + 1)] <- -K(i, j) * l
        }
      }
    }

    Bnew <- B[,2:(Mr + 1)]
    Bnew[Mr, Mr] <- 1 - K(Mr, j) + K(Mr, j) * Mr


    #Enew defined
    E <- matrix(ncol = (Mr + 2), nrow = Mr, 0)
    for(i in 1:Mr){
      for(l in 2:(Mr + 1)){
        if(l == (i + 1)){
          E[i,(l - 1)] <- K(i, j) * (l - 2)
          E[i,l] <- 1 - 2 * K(i, j) * (l - 1)
          E[i,(l + 1)] <- K(i, j) * l
        }
      }
    }

    Enew <- E[,2:(Mr + 1)]
    Enew[Mr, Mr] <- 1 + K(Mr, j) - K(Mr, j) * Mr


    G <- matrix(ncol = 1, nrow = Mr, 0)




    if(cter[j + 1] >= Mr - 1) break

    if(cter[j + 1] == 0){#when garnet didn't grow(beginning of the story)
      afterc <- (solve(Bnew)) %*% cworkm
      for(i in 1:Mr){
        c[(j + 1), i] <- afterc[i,1]
      }
    }

    if(cter[j] == 0 && cter[j + 1] == 1){#when A garnet grow at first
      Bnew[1, 1] <- 1
      Bnew[1, 2] <- 0
      Enew[1, 1] <- 1
      Enew[1, 2] <- 0
      Bnew[2, 1] <- 0
      Bnew[2, 2] <- 1 + 5 * K(2, j)
      Enew[2, 1] <- 0
      Enew[2, 2] <- 1 - 5 * K(2, j)

      c_bdry <- seq(0, c[j, 2], c[j, 2])
      c_kari <- mean(c_bdry)
      G[2, 1] <- 4 * K(2, j) * c_kari
      afterc <- (solve(Bnew)) %*% ((Enew %*% cworkm) + G)
      afterc[1, 1] <- K_d * c_kari
      mbal <- mstool %*% afterc
      mbalv <- sum(mbal)

      while((mbalv - avemasv)^2 > eps^2 && escape < 1000){
        if(mbalv > avemasv){
          c_bdry[2] <- c_kari
        }else{
          c_bdry[1] <- c_kari
        }
        c_kari <- mean(c_bdry)
        G[2, 1] <- 4 * K(2, j) * c_kari
        afterc <- (solve(Bnew)) %*% ((Enew %*% cworkm) + G)
        afterc[1, 1] <- K_d * c_kari
        mbal <- mstool %*% afterc
        mbalv <- sum(mbal)
        escape <- escape + 1
      }
      c_rg <- c_kari

      G[2, 1] <- 4 * K(2, j) * c_rg
      afterc <- (solve(Bnew)) %*% ((Enew %*% cworkm) + G)
      for(i in 1:Mr){
        c[(j + 1), i] <- afterc[i,1]
      }
      c[(j + 1), 1] <- K_d * c_rg #new garnet's c
    }



    if(cter[j] == 0 && cter[j + 1] > 1){#when garnet"s" grow at first
      Bnew[1, 1] <- 1
      Bnew[1, 2] <- 0
      Enew[1, 1] <- 1
      Enew[1, 2] <- 0
      for(i in 2:cter[j + 1]){
        Bnew[i, (i - 1)] <- 0
        Bnew[i, i] <- 1
        Bnew[i, (i + 1)] <- 0
        Enew[i, (i - 1)] <- 0
        Enew[i, i] <- 1
        Enew[i, (i + 1)] <- 0
      }
      Bnew[(cter[j + 1] + 1), cter[j + 1]] <- 0
      Bnew[(cter[j + 1] + 1), (cter[j + 1] + 1)] <- 1 + 3 * K((cter[j + 1] + 1), j) * (cter[j + 1] + 1) - K((cter[j + 1] + 1), j)
      Enew[(cter[j + 1] + 1), cter[j + 1]] <- 0
      Enew[(cter[j + 1] + 1), (cter[j + 1] + 1)] <- 1 - 3 * K((cter[j + 1] + 1), j) * (cter[j + 1] + 1) + K((cter[j + 1] + 1), j)

      c_bdry <- seq(0, c[j, (cter[j + 1] + 1)], c[j, (cter[j + 1] + 1)])
      c_kari <- mean(c_bdry)
      G[(cter[j + 1] + 1), 1] <- 4 * K((cter[j + 1] + 1), j) * c_kari * cter[j + 1]
      afterc <- (solve(Bnew)) %*% ((Enew %*% cworkm) + G)
      for(i in 1:cter[j + 1]){
        afterc[i, 1] <- K_d * c_kari
      }
      mbal <- mstool %*% afterc
      mbalv <- sum(mbal)

      while((mbalv - avemasv)^2 > eps^2 && escape < 1000){
        if(mbalv > avemasv){
          c_bdry[2] <- c_kari
        }else{
          c_bdry[1] <- c_kari
        }
        c_kari <- mean(c_bdry)
        G[(cter[j + 1] + 1), 1] <- 4 * K((cter[j + 1] + 1), j) * c_kari * cter[j + 1]
        afterc <- (solve(Bnew)) %*% ((Enew %*% cworkm) + G)
        for(i in 1:cter[j + 1]){
          afterc[i, 1] <- K_d * c_kari
        }
        mbal <- mstool %*% afterc
        mbalv <- sum(mbal)
        escape <- escape + 1
      }
      c_rg <- c_kari

      G[(cter[j + 1] + 1), 1] <- 4 * K((cter[j + 1] + 1), j) * c_rg * cter[j + 1] # matrix edit
      afterc <- (solve(Bnew)) %*% ((Enew %*% cworkm) + G) #calculation
      for(i in 1:Mr){
        c[(j + 1), i] <- afterc[i,1]
      }
      for(i in (cter[j] + 1):cter[j + 1]){
        c[(j + 1), i] <- K_d * c_rg #new garnet's c
      }
    }


    if(cter[j] > 0 && cter[j + 1] - cter[j] > 0){#when garnet grow
      #c_rg <- c_ave * 1#rim's c
      Bnew[1, 1] <- 1
      Bnew[1, 2] <- 0
      Enew[1, 1] <- 1
      Enew[1, 2] <- 0
      for(i in 2:cter[j + 1]){
        Bnew[i, (i - 1)] <- 0
        Bnew[i, i] <- 1
        Bnew[i, (i + 1)] <- 0
        Enew[i, (i - 1)] <- 0
        Enew[i, i] <- 1
        Enew[i, (i + 1)] <- 0
      }
      Bnew[(cter[j + 1] + 1), cter[j + 1]] <- 0
      Bnew[(cter[j + 1] + 1), (cter[j + 1] + 1)] <- 1 + 3 * K((cter[j + 1] + 1), j) * (cter[j + 1] + 1) - K((cter[j + 1] + 1), j)
      Enew[(cter[j + 1] + 1), cter[j + 1]] <- 0
      Enew[(cter[j + 1] + 1), (cter[j + 1] + 1)] <- 1 - 3 * K((cter[j + 1] + 1), j) * (cter[j + 1] + 1) + K((cter[j + 1] + 1), j)
      G <- matrix(ncol = 1, nrow = Mr, 0)

      c_bdry <- seq(0, c[j, (cter[j + 1] + 1)], c[j, (cter[j + 1] + 1)])
      c_kari <- mean(c_bdry)
      G[(cter[j + 1] + 1), 1] <- 4 * K((cter[j + 1] + 1), j) * c_kari * cter[j + 1]
      afterc <- (solve(Bnew)) %*% ((Enew %*% cworkm) + G)
      for(i in 1:cter[j]){
        afterc[i, 1] <- c[j, i]
      }
      for(i in (cter[j] + 1):cter[j + 1]){
        afterc[i, 1] <- K_d * c_kari
      }
      mbal <- mstool %*% afterc
      mbalv <- sum(mbal)
      escape <- 0

      while((mbalv - avemasv)^2 > eps^2 && escape < 1000){
        if(mbalv > avemasv){
          c_bdry[2] <- c_kari
        }else{
          c_bdry[1] <- c_kari
        }
        c_kari <- mean(c_bdry)
        G[(cter[j + 1] + 1), 1] <- 4 * K((cter[j + 1] + 1), j) * c_kari * cter[j + 1]
        afterc <- (solve(Bnew)) %*% ((Enew %*% cworkm) + G)
        for(i in 1:cter[j]){
          afterc[i, 1] <- c[j, i]
        }
        for(i in (cter[j] + 1):cter[j + 1]){
          afterc[i, 1] <- K_d * c_kari
        }
        mbal <- mstool %*% afterc
        mbalv <- sum(mbal)
        escape <- escape + 1
      }
      c_rg <- c_kari


      G[(cter[j + 1] + 1), 1] <- 4 * K((cter[j + 1] + 1), j) * c_rg * cter[j + 1] # matrix edit
      afterc <- (solve(Bnew)) %*% ((Enew %*% cworkm) + G) #calculation
      for(i in 1:Mr){
        c[(j + 1), i] <- afterc[i,1]
      }
      for(i in (cter[j] + 1):cter[j + 1]){
        c[(j + 1), i] <- K_d * c_rg #new garnet's c
      }
    }

    if(cter[j] == 1 && cter[j + 1] - cter[j] == 0){#when garnet didn't grow(during the story)
      Bnew[1, 1] <- 1
      Bnew[1, 2] <- 0
      Enew[1, 1] <- 1
      Enew[1, 2] <- 0
      Bnew[2, 1] <- 0
      Bnew[2, 2] <- 1 + 5 * K(2, j)
      Enew[2, 1] <- 0
      Enew[2, 2] <- 1 - 5 * K(2, j)


      if(diftype == 1){
        for(i in 1:Mr){
          c[(j + 1), i] <- c[j, i]
        }#if the matrix diff never occurs during garnet's non-growth steps
      }
      if(diftype == 2){
        afterc <- (solve(Bnew)) %*% ((Enew %*% cworkm) + G) #calculation
        for(i in 1:Mr){
          c[(j + 1), i] <- afterc[i,1]
        }#if the matrix diff occurs during garnet's non-growth steps (fixed c)
      }
      if(diftype == 3){
        Bnew[(cter[j + 1] + 1), (cter[j + 1] + 1)] <- 1 + K((cter[j + 1] + 1), j) * (cter[j + 1] + 1) + K((cter[j + 1] + 1), j)
        Enew[(cter[j + 1] + 1), (cter[j + 1] + 1)] <- 1 - K((cter[j + 1] + 1), j) * (cter[j + 1] + 1) - K((cter[j + 1] + 1), j)
        ##G[(cter[j + 1] + 1), 1] <- 0
        afterc <- (solve(Bnew)) %*% ((Enew %*% cworkm) + G) #calculation
        for(i in 1:Mr){
          c[(j + 1), i] <- afterc[i,1]
        }#if the matrix diff occurs non-eq & no flow
      }
    }







    if(cter[j] > 1 && cter[j + 1] - cter[j] == 0){#when garnet didn't grow(during the story)
      Bnew[1, 1] <- 1
      Bnew[1, 2] <- 0
      Enew[1, 1] <- 1
      Enew[1, 2] <- 0
      for(i in 2:cter[j + 1]){
        Bnew[i, (i - 1)] <- 0
        Bnew[i, i] <- 1
        Bnew[i, (i + 1)] <- 0
        Enew[i, (i - 1)] <- 0
        Enew[i, i] <- 1
        Enew[i, (i + 1)] <- 0
      }
      Bnew[(cter[j + 1] + 1), cter[j + 1]] <- 0
      Bnew[(cter[j + 1] + 1), (cter[j + 1] + 1)] <- 1 + 3 * K((cter[j + 1] + 1), j) * (cter[j + 1] + 1) - K((cter[j + 1] + 1), j)
      Enew[(cter[j + 1] + 1), cter[j + 1]] <- 0
      Enew[(cter[j + 1] + 1), (cter[j + 1] + 1)] <- 1 - 3 * K((cter[j + 1] + 1), j) * (cter[j + 1] + 1) + K((cter[j + 1] + 1), j)


      if(diftype == 1){
        for(i in 1:Mr){
          c[(j + 1), i] <- c[j, i]
        }#if the matrix diff never occurs during garnet's non-growth steps
      }
      if(diftype == 2){
        Bnew[1, 1] <- 1
        Bnew[1, 2] <- 0
        Enew[1, 1] <- 1
        Enew[1, 2] <- 0
        for(i in 2:cter[j + 1]){
          Bnew[i, (i - 1)] <- 0
          Bnew[i, i] <- 1
          Bnew[i, (i + 1)] <- 0
          Enew[i, (i - 1)] <- 0
          Enew[i, i] <- 1
          Enew[i, (i + 1)] <- 0
        }
        Bnew[(cter[j + 1] + 1), cter[j + 1]] <- 0
        Bnew[(cter[j + 1] + 1), (cter[j + 1] + 1)] <- 1 + 3 * K((cter[j + 1] + 1), j) * (cter[j + 1] + 1) - K((cter[j + 1] + 1), j)
        Enew[(cter[j + 1] + 1), cter[j + 1]] <- 0
        Enew[(cter[j + 1] + 1), (cter[j + 1] + 1)] <- 1 - 3 * K((cter[j + 1] + 1), j) * (cter[j + 1] + 1) + K((cter[j + 1] + 1), j)
        G[(cter[j + 1] + 1), 1] <- 4 * K((cter[j + 1] + 1), j) * c_rg * cter[j + 1] # matrix edit

        afterc <- (solve(Bnew)) %*% ((Enew %*% cworkm) + G) #calculation
        for(i in 1:Mr){
          c[(j + 1), i] <- afterc[i,1]
        }#if the matrix diff occurs during garnet's non-growth steps (fixed c)
      }
      if(diftype == 3){
        Bnew[(cter[j + 1] + 1), (cter[j + 1] + 1)] <- 1 + K((cter[j + 1] + 1), j) * (cter[j + 1] + 1) + K((cter[j + 1] + 1), j)
        Enew[(cter[j + 1] + 1), (cter[j + 1] + 1)] <- 1 - K((cter[j + 1] + 1), j) * (cter[j + 1] + 1) - K((cter[j + 1] + 1), j)
        ##G[(cter[j + 1] + 1), 1] <- 0
        afterc <- (solve(Bnew)) %*% ((Enew %*% cworkm) + G) #calculation
        for(i in 1:Mr){
          c[(j + 1), i] <- afterc[i,1]
        }#if the matrix diff occurs non-eq & no flow
      }
    }
    if(j == 1){
      cat("time mesh =", Mt, ">")
    }
    if((j == Mt / 2) || (j == (Mt + 1) / 2)){
      cat(">")
    }else{
      cat(".")
    }
  }


  #-------------------------------------------------------------------------------

  x <- numeric(Mr)#x axis

  for(i in 1:Mr){
    x[i] <-  (i-1) * delta_r + delta_r / 2
  }


  p <- seq(0, (1.1 * garsize), by = (1.1 * garsize))
  q <- seq(0, (c_ave * K_d * 1.1), by = (c_ave * K_d * 1.1))
  plot(p, q, ylab = ("REE [ppm]"), xlab = "Radius [cm]", xaxs = "i", yaxs = "i", axes = T, type="n")

  #lines(x, c[(Mt - 1),])#graph making
  lines(x, c[Mt,])
  text(garsize * 0.8, (c_ave * K_d), paste(format(round(time), big.mark=","), "years"))

  mas <- numeric(Mt)
  for(i in 1:Mt){
    mas[i] <- sum(c[i,] %*% mstool)
  }
  (mas[j] / mas[1] - 1) * 100 #gain of the mass(%)
  mas#mass profile
  cter#garnet points number
  return(list(mass_value = mas, garnet_mesh = cter, mass_gain_percent = (mas[j] / mas[1] - 1) * 100))


}
