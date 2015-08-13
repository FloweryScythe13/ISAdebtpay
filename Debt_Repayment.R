#########################################################################
## Debt Repayment 'Clock': Experimenting with Selling 'Personal-Stock' ##
#########################################################################


## Disclaimer ##
# It is the opinion of the writer of this code that utilizing so-called alternative
# debt repayment plans, such as selling 'personal-stock,' is high-risk, dangerous, and
# that such plans are POTENTIALLY FRADULANT. The writer strongly recommends avoiding
# such contracts. However, if an individual is determined to pursue such options, or is
# in a RARE situation where such action is financially preferable, the following program
# should give an ESTIMATE of future consequences (good or bad).


## Purpose ##
# The model below takes in parameters from the 'Constants'
# section and finds the optimal debt repayment plan (using linear programming)
# for the case where one either wishes to avoid selling 'personal-stocks' or
# is comfortable selling 'personal-stocks'


## Assumptions ##
# Interest on Debt starts Now
# The Interest Rate on Debt is in 'compounded annually' form
# If any 'Personal-Stocks' are sold, that occurs in the first year (technically immediatly).
# Any income left over at the end of the year (after expenses, debt payment, stock payment, etc.) will be invested at the given rate.
# Each year, the user will pay at minimum the interest on the debt owed
# The user starts with no saved money
# Income payed to holders of 'Personal-Stocks' is only a percentage of paycheck and does not include income from other sources such as investment return
#not
#Note to Eric: assumes that income from interest on savings is not paid out to 
#investors, only income from paychecks

## Constants ##
Debt <- 40000 # Dollars
Debt.Interest <- 0.06 # percent in decimal notation, default at 0.06 = 6%
Max.Years.Repay.Debt <- 25 # Maximum number of years in which to repay debt

Income <- 18000 # annual income in dollars
Rate.Income.Increase <- 0.05 # predicted percentage increase in annual income
Expenses <- 0.60 # annual fixed expenses (food, rent, travel, etc.) expressed as a percentage of annual income

Interest.Saved.Income <- 0.03 # annual rate of return on any 'saved' income

Allow.Personal <- FALSE # 'FALSE' if user wants to prevent algorithm from selling 'personal-stocks', 'TRUE' if user wants to allow algoritm to sell 'personal-stocks.' Overriden by 'Force.Personal'

Price.Personal <- 1000 # per-share price of a 'personal-stock.' In other words, the amount of money you receive for each stock you sell.
Rate.Personal <- 0.01 # per-stock percentage of future income which must be payed to owner of 'personal-stock.' So, if the per-stock rate is .1%, then a person who sells 10 stocks will pay 1% of income spread accross all stock holders.
Personal.Start <- 11 # first year in which one must pay percentage of income to owner of stock
Personal.End <- 20 # last year in which one must pay percentage of income to owner 

Force.Personal <- FALSE # 'FALSE' if algorithm should choose whether to sell 'personal-stocks,' 'TRUE' if user wants to mandate that stocks be sold. Overrides 'Allow.Personal'
Force.Number <- 100 # Number of stocks to be sold if we force the algorithm to sell stocks


## Function ##
# The user should avoid making any changes here unless he/she wants to change the 
# functionality of the optimization algorithm
Debt.Clock <- function(D,DI,MYRP,I,RII,E,ISI,AP=F,PP=0,RP=0,PS=0,PE=0,FP=F,FN=0) {
  Num.Yrs <- max(MYRP,PE)
  Ann.Inc <- I * (1+RII)^(1:Num.Yrs)
  Ann.Ext <- Ann.Inc * (1-E)
  
  Personal.Pay <- Ann.Ext[PS:PE] * RP
  
  # vars: Funds 1:Num.Yrs, Debt.Pay 1:MYRP, Personal.Stocks, Debt.Remaining 0:MYRP
  DP <- Num.Yrs+1; PStk <- DP+MYRP; DR <- PStk+1
  
  objective <- c(rep(0,Num.Yrs-1),1,rep(0,MYRP+MYRP+2))
  
  Mat <- matrix(0,ncol=(Num.Yrs+MYRP+MYRP+2),nrow=1+Num.Yrs+1+MYRP+MYRP+Num.Yrs)
  const.dir <- c()
  const.rhs <- c()
  
  # Intialize values
  Mat[1,DR] <- 1; const.dir[1] <- "="; const.rhs[1] <- D
  
  # Add Income
  const.dir <- c(const.dir,rep("=",Num.Yrs)); const.rhs <- c(const.rhs,Ann.Ext)
  Mat[2,1] <- 1
  Mat[2,PStk] <- -PP
  for (i in 2:Num.Yrs) {
    Mat[i+1,i] <- 1
    Mat[i+1,i-1] <- -(1+ISI)
  }
  
  # Add Costs
  for (i in 1:MYRP) {
    Mat[i+1,DP+i-1] <- 1
  }
  for (i in PS:PE) {
    Mat[i+1,PStk] <- Personal.Pay[i+1-PS]
  }
  
  # Drive Debt to 0
  Mat[Num.Yrs+2,DR+MYRP] <- 1
  const.dir <- c(const.dir,"="); const.rhs <- c(const.rhs,0)
  
  # Make Minimum Payments
  for (i in 1:MYRP) {
    Mat[Num.Yrs+2+i,DP+i-1] <- 1
    Mat[Num.Yrs+2+i,DR+i-1] <- -DI
  }
  const.dir <- c(const.dir,rep(">=",MYRP)); const.rhs <- c(const.rhs,rep(0,MYRP))
  
  # Interest on Debt and Reductions in Debt
  for (i in 1:MYRP) {
    Mat[Num.Yrs+2+MYRP+i,DR+i] <- 1
    Mat[Num.Yrs+2+MYRP+i,DR+i-1] <- -(DI+1)
    Mat[Num.Yrs+2+MYRP+i,DP+i-1] <- (DI+1)
  }
  const.dir <- c(const.dir,rep("=",MYRP)); const.rhs <- c(const.rhs,rep(0,MYRP))
  
  # Keep Funds non-negative
  for (i in 1:Num.Yrs) {
    Mat[Num.Yrs+2+MYRP+MYRP+i,i] <- 1
  }
  const.dir <- c(const.dir,rep(">=",Num.Yrs)); const.rhs <- c(const.rhs,rep(0,Num.Yrs))
  
  # Allow or Disallow or Force 'Personal-Stocks'
  if (!AP) {
    New.Mat <- matrix(0,nrow=1,ncol=Num.Yrs+MYRP+MYRP+2)
    const.dir <- c(const.dir,"=")
    New.Mat[1,PStk] <- 1
    Mat <- rbind(Mat,New.Mat)
    if (FP) {
      const.rhs <- c(const.rhs,FN)
    }
    else {
      const.rhs <- c(const.rhs,0)
    }
  }
  if (AP) {
    if (FP) {
      New.Mat <- matrix(0,nrow=1,ncol=Num.Yrs+MYRP+MYRP+2)
      const.dir <- c(const.dir,"=")
      New.Mat[1,PStk] <- 1
      Mat <- rbind(Mat,New.Mat)
      const.rhs <- c(const.rhs,FN)
    }
  }
  
  # Linear Program
  soln <- lp(direction="max",objective.in=objective,const.mat=Mat,
             const.dir=const.dir,const.rhs=const.rhs)
  
  # Assign Solution Values to Relevant Vectors
  Fin.Funds <- c(0,soln$solution[1:Num.Yrs])
  Fin.DebtPayments <- c(0,soln$solution[(Num.Yrs+1):(Num.Yrs+MYRP)])
  Fin.DebtPayments <- c(Fin.DebtPayments,rep(0,Num.Yrs-MYRP))
  Fin.PS <- soln$solution[(Num.Yrs+MYRP+1)]
  Fin.RemainingDebt <- soln$solution[(Num.Yrs+MYRP+2):(Num.Yrs+MYRP+2+MYRP)]
  Fin.RemainingDebt <- c(Fin.RemainingDebt,rep(0,Num.Yrs-MYRP))
  
  # Return Value
  return(list(Fin.Funds,Fin.DebtPayments,Fin.PS,Fin.RemainingDebt,Num.Yrs))
}


## Executible ##
library(linprog); library(lpSolve)

# No 'Personal-Stocks'
mod.1 <- Debt.Clock(Debt,Debt.Interest,Max.Years.Repay.Debt,Income,Rate.Income.Increase,
                    Expenses,Interest.Saved.Income,Allow.Personal,Price.Personal,
                    Rate.Personal,Personal.Start,Personal.End,Force.Personal,Force.Number)
mod.1.Funds <- mod.1[[1]]
mod.1.Debt.Payments <- mod.1[[2]]
mod.1.Personal.Stocks <- mod.1[[3]]; mod.1.Personal.Stocks
mod.1.Remaining.Debt <- mod.1[[4]]
mod.1.max <- max(c(mod.1.Funds,mod.1.Debt.Payments,mod.1.Remaining.Debt))
plot(0:mod.1[[5]],mod.1.Funds/1000,type="l",lwd=3,ylim=c(0,1.2*mod.1.max/1000),
     xlab="Time (Years)",ylab="Dollars (Thousands)",main="Optimal Plan w/o 'PS'")
lines(0:mod.1[[5]],mod.1.Debt.Payments/1000,lwd=3,col=2)
lines(0:mod.1[[5]],mod.1.Remaining.Debt/1000,lwd=3,col=3)
legend(0,1.1*mod.1.max/1000,col=c(1,2,3),lwd=3,lty=1,legend=c("Funds in Account",
                                                              "Annual Debt Payments",
                                                              "Remaining Debt"))

# With 'Personal-Stocks'
Allow.Personal <- TRUE
mod.2 <- Debt.Clock(Debt,Debt.Interest,Max.Years.Repay.Debt,Income,Rate.Income.Increase,
                    Expenses,Interest.Saved.Income,Allow.Personal,Price.Personal,
                    Rate.Personal,Personal.Start,Personal.End,Force.Personal,Force.Number)
mod.2.Funds <- mod.2[[1]]
mod.2.Debt.Payments <- mod.2[[2]]
mod.2.Personal.Stocks <- mod.2[[3]]; mod.2.Personal.Stocks
mod.2.Remaining.Debt <- mod.2[[4]]
mod.2.max <- max(c(mod.2.Funds,mod.2.Debt.Payments,mod.2.Remaining.Debt))
plot(0:mod.2[[5]],mod.2.Funds/1000,type="l",lwd=3,ylim=c(0,1.2*mod.2.max/1000),
     xlab="Time (Years)",ylab="Dollars (Thousands)",main="Optimal Plan w/o 'PS'")
lines(0:mod.2[[5]],mod.2.Debt.Payments/1000,lwd=3,col=2)
lines(0:mod.2[[5]],mod.2.Remaining.Debt/1000,lwd=3,col=3)
legend(0,1.1*mod.2.max/1000,col=c(1,2,3),lwd=3,lty=1,legend=c("Funds in Account",
                                                              "Annual Debt Payments",
                                                              "Remaining Debt"))

# Forcing 'Personal-Stocks'
Force.Personal <- TRUE
mod.3 <- Debt.Clock(Debt,Debt.Interest,Max.Years.Repay.Debt,Income,Rate.Income.Increase,
                    Expenses,Interest.Saved.Income,Allow.Personal,Price.Personal,
                    Rate.Personal,Personal.Start,Personal.End,Force.Personal,Force.Number)
mod.3.Funds <- mod.3[[1]]
mod.3.Debt.Payments <- mod.3[[2]]
mod.3.Personal.Stocks <- mod.3[[3]]; mod.3.Personal.Stocks
mod.3.Remaining.Debt <- mod.3[[4]]
mod.3.max <- max(c(mod.3.Funds,mod.3.Debt.Payments,mod.3.Remaining.Debt))
plot(0:mod.3[[5]],mod.3.Funds/1000,type="l",lwd=3,ylim=c(0,1.2*mod.3.max/1000),
     xlab="Time (Years)",ylab="Dollars (Thousands)",main="Optimal Plan w/o 'PS'")
lines(0:mod.3[[5]],mod.3.Debt.Payments/1000,lwd=3,col=2)
lines(0:mod.3[[5]],mod.3.Remaining.Debt/1000,lwd=3,col=3)
legend(0,1.1*mod.3.max/1000,col=c(1,2,3),lwd=3,lty=1,legend=c("Funds in Account",
                                                              "Annual Debt Payments",
                                                              "Remaining Debt"))

# Comparing Final Funds in Each Case
comp.max <- max(mod.1.Funds,mod.2.Funds,mod.3.Funds)
plot(0:mod.1[[5]],mod.1.Funds/1000,type="l",lwd=3,ylim=c(0,1.2*comp.max/1000),
     xlab="Time (Years)",ylab="Dollars (Thousands)",main="Final Funds Comparison")
lines(0:mod.1[[5]],mod.2.Funds/1000,lwd=3,col=2,lty=2)
lines(0:mod.1[[5]],mod.3.Funds/1000,lwd=3,col=3,lty=3)
legend(0,1.1*comp.max/1000,col=c(1,2,3),lwd=3,lty=c(1,2,3),legend=c("No 'Personal-Stocks",
                                                              "Optimal 'Personal-Stocks'",
                                                              "Forced 'Personal-Stocks'"))
