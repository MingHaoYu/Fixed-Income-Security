#title : 405-Computational Methods HW8
#author: Ming-Hao Yu
#date  : 2018-03-03

Bfunction <- function(t, T, kappa) {
    return( (1-exp(-kappa*(T-t)))/kappa )
}

Afunction <- function(t, T, kappa, rbar, sigma) {
    B <- (1-exp(-kappa*(T-t)))/kappa
    return( exp((rbar-sigma^2/(2*kappa^2))*(B-(T-t)) - sigma^2/(4*kappa)*B^2) )
}

Vfunction <- function(t, T, a, b, sigma, eta, rho) {
    V <- sigma^2/a^2*(T-t+2/a*exp(-a*(T-t))-0.5/a*exp(-2*a*(T-t))-1.5/a)+
      eta^2/b^2*(T-t+2/b*exp(-b*(T-t))-0.5/b*exp(-2*b*(T-t))-1.5/b)+
      2*rho*sigma*eta/(a*b)*(T-t+(exp(-a*(T-t))-1)/a+(exp(-b*(T-t))-1)/b-(exp(-(a+b)*(T-t))-1)/(a+b))
    return(V)
}

CalculateRStar <- function(r0, T, k, kappa, rbar, sigma, Coupon, CouponT) {
    stop <- FALSE
    r <- r0
    #FaceValue <- 1000
    #Coupon <- c(rep(30, 7), FaceValue+30)
    #CouponT <- seq(0.5, length.out=length(Coupon), by=0.5)
    #T <- 3/12
    #k <- 980
  
    while(!stop) {
        f <- 0
        f2 <- 0
        for(i in 1:length(CouponT)) {
            f <- f + Coupon[i]*Afunction(t=T, T=CouponT[i], kappa=kappa, rbar=rbar, sigma=sigma)*exp(-Bfunction(t=T, T=CouponT[i], kappa=kappa)*r)
            f2 <- f2 + Coupon[i]*Afunction(t=T, T=CouponT[i], kappa=kappa, rbar=rbar, sigma=sigma)*exp(-Bfunction(t=T, T=CouponT[i], kappa=kappa)*(r+0.001))
        }
        if((f-k)==0) {
            stop <- TRUE
        }else{
            r <- r - (f-k)/((f2-f)/0.001)
        }
    }
    return(r)
}

CIRExplicitPureDiscountBond <- function(rt, t, T, sigma, kappa, rbar) {
    h1 <- sqrt(kappa^2+2*sigma^2)
    h2 <- 0.5*(kappa+h1)
    h3 <- 2*kappa*rbar/(sigma^2)
    
    A <- ( h1*exp(h2*(T-t))/(h2*(exp(h1*(T-t))-1)+h1) )^(h3)
    B <- (exp(h1*(T-t))-1)/(h2*(exp(h1*(T-t))-1)+h1 )
    return(A*exp(-B*rt))
}

G2ppExplicitPireDiscountBond <- function(t, T, phit, xt, yt, a, b, sigma, eta, rho) {
    P <- exp(-phit*(T-t)-xt*(1-exp(-a*(T-t)))/a-yt*(1-exp(-b*(T-t)))/b+0.5*Vfunction(t, T, a, b, sigma, eta, rho))
    return(P)
}

PureDiscountBond <- function(r0, t, T, sigma, kappa, rbar, simulation=1000, type="Vasicek") {
    dt <- 1/252
    steps <- (T-t)/dt
    dW <- sqrt(dt)*rnorm(simulation*steps)
    dWTable <- t(matrix(ncol=simulation, nrow=steps, dW))
    rTable <- matrix(nrow=simulation, ncol=steps+1, r0)
    FaceValue <- 1
    
    if(toupper(type)=="VASICEK") {
        for(i in 1:steps) {
            rTable[, i+1] <- rTable[, i] + kappa*(rbar-rTable[, i])*dt + sigma*dWTable[, i]
        }
    }else if (toupper(type)=="CIR") {
        for(i in 1:steps) {
            rTable[, i+1] <- rTable[, i] + kappa*(rbar-rTable[, i])*dt + sigma*sqrt(rTable[, i])*dWTable[, i]
        }
    }else {cat("No", type, "type of model for interest rate.")}
    
    price <- FaceValue*mean(exp(-(rowSums(rTable)*dt)))
    return(price)
}

CouponPayingBond <- function(r0, t, T, Coupon, simulation=10000) {
    if(length(T) != length(Coupon)) {
        cat("Coupon amount and payment date do not match.")
        return(0)
    }
  
    sigma <- 0.1
    kappa <- 0.82
    rbar <- 0.05
    dt <- 1/252
    steps <- max(T-t)/dt
    
    dW <- sqrt(dt)*rnorm(simulation*steps)
    dWTable <- t(matrix(ncol=simulation, nrow=steps, dW))
    rTable <- matrix(nrow=simulation, ncol=steps+1, r0)
    RTable <- matrix(nrow=simulation, ncol=length(T))
    for(i in 1:steps) {
        rTable[, i+1] <- rTable[, i] + kappa*(rbar-rTable[, i])*dt + sigma*dWTable[, i]
    }
    for(i in 1:length(T)) {
        RTable[, i] <- rowSums(rTable[, 1:((T[i]-t)/dt)])*dt
    }
    
    return(mean(exp(-RTable)%*%Coupon))
}

#Problem 1
Problem1 <- function(){
    r0 <- 0.05
    sigma <- 0.1
    kappa <- 0.82
    rbar <- 0.05
    dt <- 1/252
    
    #(a)
    FaceValue <- 1000
    T <- 0.5
    
    pricea <- FaceValue*PureDiscountBond(r0=r0, t=0, T=0.5, sigma=sigma, kappa=kappa, rbar=rbar)
    
    #(b)
    Coupon <- c(rep(30, 7), FaceValue+30)
    T <- seq(0.5, length.out=length(Coupon), by=0.5)
    
    priceb <- CouponPayingBond(r0=r0, t=0, T=T, Coupon=Coupon, simulation=10000)
    
    #(c)
    T <- 3/12
    S <- 0.5
    k <- 980
    simulation <- 100000
    steps <- T/dt
    dW <- sqrt(dt)*rnorm(simulation*steps)
    dWTable <- t(matrix(ncol=simulation, nrow=steps, dW))
    rTable <- matrix(nrow=simulation, ncol=steps+1, r0)
    for(i in 1:steps) {
        rTable[, i+1] <- rTable[, i] + kappa*(rbar-rTable[, i])*dt + sigma*dWTable[, i]
    }
    R <- exp(-rowSums(rTable)*dt)
    P <- Afunction(t=T, T=S, kappa=kappa, rbar=rbar, sigma=sigma)*exp(-Bfunction(t=T, T=S, kappa=kappa)*rTable[, steps+1])
    payoffc <- pmax(FaceValue*P-k, 0)
    pricec <- mean(R*payoffc)
    
    #(d)
    Coupon <- c(rep(30, 7), FaceValue+30)
    CouponT <- seq(0.5, length.out=length(Coupon), by=0.5)
    T <- 3/12
    k <- 980
    dt <- 1/252
    simulation <- 10000
    steps <- T/dt
    dW <- sqrt(dt)*rnorm(simulation*steps)
    dWTable <- t(matrix(ncol=simulation, nrow=steps, dW))
    rTable <- matrix(nrow=simulation, ncol=steps+1, r0)
    for(i in 1:steps) {
        rTable[, i+1] <- rTable[, i] + kappa*(rbar-rTable[, i])*dt + sigma*dWTable[, i]
    }
    R <- exp(-rowSums(rTable[,1:63])*dt)
    payoffd <- vector(length=simulation)
    for(i in 1:simulation) {
        payoffd[i] <- max(CouponPayingBond(r0=rTable[i, steps+1], t=T, T=CouponT, Coupon=Coupon, simulation=100)-k, 0)
    }
    priced <- mean(R*payoffd)
    
    #(e)
    #t = 0, T = T, Ti = CouponT, for all Ti > T
    ##P <- matrix(nrow=simulation, ncol=length(Coupon))
    ##for(i in 1:length(Coupon)) {
    ##    P[,i] <- Afunction(t=T, T=CouponT[i], kappa=kappa, rbar=rbar, sigma=sigma)*exp(-Bfunction(t=T, T=CouponT[i], kappa=kappa)*rTable[, steps+1])
    ##}
    ##pricee <- mean(pmax(P%*%Coupon-k, 0)*R)
    
    #(e)
    ##  CBOP explicit formula
    ## t = 0, T = T, Ti = CouponT, for all Ti > T
    Coupon <- c(rep(30, 7), FaceValue+30)
    CouponT <- seq(0.5, length.out=length(Coupon), by=0.5)
    T <- 3/12
    k <- 980
    ki <- vector(length=length(CouponT))
    pricee <- 0
    rstar <- CalculateRStar(r0=r0, T=T, k=k, kappa=kappa, rbar=rbar, sigma=sigma, Coupon=Coupon, CouponT=CouponT)
    for(i in 1:length(ki)) {
        ki[i] <- Afunction(t=T, T=CouponT[i], kappa=kappa, rbar=rbar, sigma=sigma)*exp(-Bfunction(t=T, T=CouponT[i], kappa=kappa)*rstar)
    }
    for(i in 1:length(Coupon)) {
        #t = 0, T = T, Ti = CouponT, for all Ti > T
        PTi <- Afunction(t=0, T=CouponT[i], kappa=kappa, rbar=rbar, sigma=sigma)*exp(-Bfunction(t=0, T=CouponT[i], kappa=kappa)*r0)
        PT <- Afunction(t=0, T=T, kappa=kappa, rbar=rbar, sigma=sigma)*exp(-Bfunction(t=0, T=T, kappa=kappa)*r0)
        sigmap <- sigma/kappa*(1-exp(-kappa*(CouponT[i]-T)))*sqrt(0.5*(1-exp(-2*kappa*(T-0)))/kappa)
        d1 <- 1/sigmap*log(PTi/(ki[i]*PT))+sigmap/2
        d2 <- 1/sigmap*log(PTi/(ki[i]*PT))-sigmap/2
        pricee <- pricee + Coupon[i]*(PTi*pnorm(d1) - ki[i]*PT*pnorm(d2))
    }
    
    ans <- c(pricea, priceb, pricec, priced, pricee)
    names(ans) <- c("Problem1_a", "Problem1_b", "Problem1_c", "Problem1_d", "Problem1_e")
    return(ans)
}

Problem2 <- function(){
    r0 <- 0.05
    sigma <- 0.12
    kappa <- 0.92
    rbar <- 0.055
    k <- 980
    dt <- 1/252
    
    #(a)
    FaceValue <- 1000
    T <- 0.5
    S <- 1
    simulation <- 10000
    steps <- T/dt
    dW <- sqrt(dt)*rnorm(simulation*steps)
    dWTable <- t(matrix(ncol=simulation, nrow=steps, dW))
    rTable <- matrix(nrow=simulation, ncol=steps+1, r0)
    for(i in 1:steps) {
        rTable[, i+1] <- rTable[, i] + kappa*(rbar-rTable[, i])*dt + sigma*sqrt(rTable[, i])*dWTable[, i]
    }
    R <- exp(-rowSums(rTable)*dt)
    payoff <- pmax(FaceValue*mapply(FUN=PureDiscountBond, 
                                    rTable[, steps+1], rep(T, simulation), rep(S, simulation), 
                                    rep(sigma, simulation), rep(kappa, simulation), 
                                    rep(rbar, simulation), rep(1000, simulation), rep("CIR", simulation)) -k, 0)
    price1 <- mean(R*payoff)
    
    #(b)
    M <- T/dt
    N <- 80
    dr <- 0.4/N
    j <- seq(N-1,1)
    pu <- -dt*(0.5*sigma^2*j*dr/(dr^2)+kappa*(rbar-j*dr)/(2*dr))
    pm <- 1+dt*(sigma^2*j*dr/(dr^2) + j*dr)
    pd <- -dt*(0.5*sigma^2*j*dr/(dr^2)-kappa*(rbar-j*dr)/(2*dr))
    A <- matrix(nrow=N+1, ncol=N+1, 0)
    F <- matrix(nrow=N+1, ncol=M+1)
    rList <- c(N, j, 0)*dr
    B <- matrix(nrow=N+1, ncol=1, 0)
    B[1,1] <- dr
    P <- FaceValue*mapply(FUN=CIRExplicitPureDiscountBond, 
                                         rList, rep(T, length(rList)), rep(S, length(rList)), 
                                         rep(sigma, length(rList)), rep(kappa, length(rList)), 
                                         rep(rbar, length(rList)))
    F[1:(N+1), M+1] <- pmax(P[1:(N+1)]-k, 0)
    #A matrix
    for(i in 2:N) {
        A[i, (i-1):(i+1)] <- c(pu[i-1], pm[i-1], pd[i-1])
    }
    A[1, 1:3] <- c(1, -1, 0)
    A[N+1, (N-1):(N+1)] <- c(0, 1, -1)
    invA <- solve(A)
    
    for(i in M:1) {
        B[2:N] <- F[2:N, i+1]
        F[, i] <- invA%*%B
    }
    #linear interpolation
    index <- length(rList)
    for(i in 2:length(rList)) {
        if(r0 > rList[i]) { 
            index <- i
            break
        }
    }
    price2 <- F[index] + (F[index-1, 1]-F[index, 1])*(r0 - rList[index]) / (rList[index-1] - rList[index])
    
    #(c)
    payoff3 <- pmax(FaceValue*CIRExplicitPureDiscountBond(rt=rTable[,steps+1], t=T, T=S, sigma=sigma, kappa=kappa, rbar=rbar)-k, 0)
    price3 <- mean(R*payoff3)
    
    #t=0, T=T, S=S
    rt <- r0
    theta <- sqrt(kappa^2+2*sigma^2)
    phi <- 2*theta/(sigma^2*(exp(theta*(T-0))-1))
    psi <- (kappa+theta)/(sigma^2)
    h1 <- sqrt(kappa^2+2*sigma^2)
    h2 <- 0.5*(kappa+h1)
    h3 <- 2*kappa*rbar/(sigma^2)
    
    A <- ( h1*exp(h2*(S-T))/(h2*(exp(h1*(S-T))-1)+h1) )^(h3)
    B <- (exp(h1*(S-T))-1)/(h2*(exp(h1*(S-T))-1)+h1 )
    rstar <- log(FaceValue*A/k)/B
    P_tS <- CIRExplicitPureDiscountBond(t=0, T=S, kappa=kappa, rbar=rbar, sigma=sigma, rt=r0)
    P_tT <- CIRExplicitPureDiscountBond(t=0, T=T, kappa=kappa, rbar=rbar, sigma=sigma, rt=r0)
    pricec <- FaceValue*P_tS*pchisq( q=2*rstar*(phi+psi+B), df=4*kappa*rbar/sigma^2, ncp=(2*phi^2*rt*exp(theta*(T-0)))/(phi+psi+B) )-
      k*P_tT*pchisq( q=2*rstar*(phi+psi), df=4*kappa*rbar/sigma^2, ncp=(2*phi^2*rt*exp(theta*(T-0)))/(phi+psi))
    
    ans <- c(price1, price2, pricec)
    names(ans) <- c("Problem2_a", "Problem2_b", "Problem2_c")
    
    return(ans)
    
}

Problem3 <- function(){
    t <- 0
    T <- 0.5
    S <- 1
    dt <- 1/252
    k <- 950
    FaceValue <- 1000
    r0 <- 0.03
    rho <- 0.7
    a <- 0.1 
    b <- 0.3
    sigma <- 0.03
    eta <- 0.08
    phi <- 0.03
    simulation <- 10000
    steps <- T/dt
    payoff <- vector(length=simulation)
    XTable <- YTable <- matrix(nrow=simulation, ncol=steps+1, 0)
    rTable <- matrix(nrow=simulation, ncol=steps+1, r0)
    dw1 <- sqrt(dt)*rnorm(simulation*steps)
    dw2 <- rho*dw1+sqrt(1-rho^2)*sqrt(dt)*rnorm(simulation*steps)
    dWTable1 <- t(matrix(ncol=simulation, nrow=steps, dw1))
    dWTable2 <- t(matrix(ncol=simulation, nrow=steps, dw2))
    
    for(i in 1:steps) { 
        XTable[, i+1] <- XTable[, i] - a*XTable[, i]*dt + sigma*dWTable1[, i] 
        YTable[, i+1] <- YTable[, i] - b*YTable[, i]*dt + eta*dWTable2[, i]
        rTable[, i+1] <- XTable[, i+1]+YTable[, i+1]+phi
    }
    #rTable <- XTable+YTable+phi
    R <- exp(-rowSums(rTable)*dt)
    
    for(i in 1:simulation) {
        simulation2 <- 300
        steps2 <- (S-T)/dt
        XTable2 <- matrix(nrow=simulation2, ncol=steps2+1, XTable[i, steps+1])
        YTable2 <- matrix(nrow=simulation2, ncol=steps2+1, YTable[i, steps+1])
        rTable2 <- matrix(nrow=simulation2, ncol=steps2+1, rTable[i, steps+1])
        dw3 <- sqrt(dt)*rnorm(simulation2*steps2)
        dw4 <- rho*dw3+sqrt(1-rho^2)*sqrt(dt)*rnorm(simulation2*steps2)
        dWTable3 <- t(matrix(ncol=simulation2, nrow=steps2, dw3))
        dWTable4 <- t(matrix(ncol=simulation2, nrow=steps2, dw4))
        
        for(j in 1:steps2) { 
            XTable2[, j+1] <- XTable2[, j] - a*XTable2[, j]*dt + sigma*dWTable3[, j] 
            YTable2[, j+1] <- YTable2[, j] - b*YTable2[, j]*dt + eta*dWTable4[, j]
            rTable2[, j+1] <- XTable2[, j+1]+YTable2[, j+1]+phi
        }
        #rTable2 <- XTable2+YTable2+phi
        R2 <- exp(-rowSums(rTable2)*dt)
        payoff[i] <- max(k-mean(FaceValue*R2), 0)
    }
    price1 <- mean(payoff*R)
    
    #explicit
    sigmaSquare <- sigma^2/(2*a^3)*(1-exp(-a*(S-T)))^2*(1-exp(-2*a*(T-t)))+eta^2/(2*b^3)*(1-exp(-b*(T-t)))^2*(1-exp(-2*b*(T-t)))+2*rho*sigma*eta/(a*b*(a+b))*(1-exp(-a*(S-T)))*(1-exp(-b*(S-T)))*(1-exp(-(a+b)*(T-t)))
    P_tS <- G2ppExplicitPireDiscountBond(t, S, phi, 0, 0, a, b, sigma, eta, rho)
    P_tT <- G2ppExplicitPireDiscountBond(t, T, phi, 0, 0, a, b, sigma, eta, rho)
    CapSigma <- sqrt(sigmaSquare)
    price <- -FaceValue*P_tS*pnorm(log(k*P_tT/(FaceValue*P_tS))/CapSigma-0.5*CapSigma) + P_tT*k*pnorm(log(k*P_tT/(FaceValue*P_tS))/CapSigma+0.5*CapSigma)
    
    ans <- c(price1, price)
    names(ans) <- c("price_MonteCarlo", "price_explicit")
    return(ans)
}

Problem1()
Problem2()
Problem3()
