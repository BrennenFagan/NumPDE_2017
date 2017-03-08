D <- 0.04
nu <- 0.053
beta <- 0.43
gamma <- 0.05
w <- 0.34
hp <- 1
hm <- 1

h <- 0.02
N <- 50
x <- (1:(N-1))*h

#Fix tau
tau <- (h-h*11/12)^2/(2*D)

T <- 10

M <- 10/tau

#M1 is effectively a smoothing value
M1 <- 400
M2 <- round(M/M1)

t <- (0:M1)*M2*tau

wp <- matrix(0, nrow=N-1, ncol=M1+1)
wz <- matrix(0, nrow=N-1, ncol=M1+1)

Pn <- 0.035 + ifelse(x<1/4, (sin(4*pi*x))/16, 0)
Zn <- rep(0.046, length(x))

#fix assignment of Pn, Zn
wp[, 1] <- Pn
wz[, 1] <- Zn


for (j in 1:M1) {
  Pn <- wp[, j]
  Zn <- wz[, j]
  for (j2 in 1:M2) {
    
    P <- Pn
    Pp <- c(P[2:(N-1)], P[N-1])
    Pm <- c(P[1], P[1:(N-2)])
    
    Z <- Zn
    Zp <- c(Z[2:(N-1)], Z[N-1])
    Zm <- c(Z[1], Z[1:(N-2)])
    
    Pn <- P + tau * (
      D/h^2 * (Pp - 2*P + Pm) + 
        hm/(2*h^2) * ((Pp+P)*(Zp-Z)-(P+Pm)*(Z-Zm)) +
        beta*P*(1-P)-Z*P^2/(P^2+nu^2)
    )
    
    Zn <- Z + tau * (
      D/h^2 * (Zp - 2*Z + Zm) - 
        hm/(2*h^2) * ((Zp+Z)*(Pp-P)-(Z+Zm)*(P-Pm)) +
        gamma*Z*(P^2/(P^2+nu^2) - w)
    )
  }
  wp[, j+1] <- Pn
  wz[, j+1] <- Zn
}

library(plot3Drgl)
persp3D(x, t, wz,
        xlab="x", ylab="t", zlab="Z",
        ticktype="detailed", nticks=4)

plotrgl(lighting = TRUE)

persp3D(x, t, wp,
        xlab="x", ylab="t", zlab="P",
        ticktype="detailed", nticks=4)

plotrgl(lighting = TRUE)
