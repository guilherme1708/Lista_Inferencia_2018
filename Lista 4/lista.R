library(MASS)

x <- c(1.19, 1.33, 1.29, 0.97, 0.57, 0.26, 1.46, 0.73, 0.45, 0.85,
1.67, 0.56, 0.45, 0.35, 0.52, 1.32, 1.22, 1.09, 0.27, 0.34,
0.59, 0.78, 0.55, 1.29, 1.11, 1.04, 1.21, 0.38, 0.61, 1.12,
0.72, 0.55, 0.90, 0.26, 0.90, 0.54, 0.99, 0.67, 1.36, 0.18,
0.58, 0.22, 1.38, 1.36, 0.35, 1.43, 0.04, 0.26, 0.86, 1.06,
1.47, 0.42, 0.62, 0.58, 0.65, 0.54, 0.76, 0.93, 1.15, 0.92,
1.95, 1.29, 0.64, 0.13, 1.70, 1.00, 0.75, 1.09, 1.40, 1.26,
0.87, 0.80, 0.67, 0.47, 0.66, 0.33, 0.56, 1.01, 1.54, 0.46,
1.39, 1.30, 1.17, 1.60, 1.16, 0.93, 1.27, 0.20, 1.17, 0.42,
1.53, 0.31, 1.31, 1.20, 0.75, 0.72, 1.97, 1.26, 0.48, 0.27)

MV = function(x,theta0=c(1,1),epsilon=10^(-5), int=1000){
      n <- length(x)
      
      U <- function(t1,t2){
        U1 <- n/t1 - sum(x^t2)
        U2 <- n/t2 + sum(log(x)) - t1*sum(x^t2)*log(sum(x)) + sum(log(x))
        return (c(U1,U2))
      }
      
      H <- function(t1,t2){
        H11 <- -n/t1^2
        H12 <- sum(x^t2)*log(sum(x))
        H22 <- -n/t2^2 -t1*sum(x^t2)*log(sum(x))^2
        M <- matrix(c(H11,H12,H12,H22),2,2)
        return (M)
      }
      
      erro <- 10
      j <- 0
      t1 <- numeric ()
      t2 <- numeric()
      t1[1]<-theta0[1]
      t2[1]<-theta0[2]
      
      while(erro > epsilon & j < int){
        j <- j+1
        Aux <- c(t1[j],t2[j]) - solve(H(t1[j],t2[j]))%*%U(t1[j],t2[j])
        t1[j+1] <- Aux[1]
        t2[j+1] <- Aux[2]
        erro <- max(abs(c(t1[j+1]-t1[j],t2[j+1]-t2[j])))
        print(erro)
      }
      
      S <- list()
      S$erro <- erro
      S$Iteracoes <- j
      S$theta <- c(t1[j+1],t2[j+1])
      S$H <- H(t1[j+1],t2[j+1])
      S$U <- U(t1[j+1],t2[j+1])
      return(S)
}

#####################################################################

MV1 = function(x,theta0=c(1,1),epsilon=10^(-5), int=1000, prop=0.01, imprimir = FALSE){
  n <-length(x)
  
  U <- function(t1,t2){
    U1 <- n/t1 - sum(x^t2)
    U2 <- n/t2 + sum(log(x)) - t1*sum(x^t2)*log(sum(x)) + sum(log(x))
    return (c(U1,U2))
  }
  
  #Função escore considerando a reparametrização
  U1 = function(aa,bb){
    w1 <- exp(aa)*U(exp(aa),exp(bb))[1]
    w2 <- exp(bb)*U(exp(aa),exp(bb))[2]
    return(c(w1,w2))
  }
  
  H <- function(t1,t2){
    H11 <- -n/t1^2
    H12 <- sum(x^t2)*log(sum(x))
    H22 <- -n/t2^2 -t1*sum(x^t2)*log(sum(x))^2
    M <- matrix(c(H11,H12,H12,H22),2,2)
    return (M)
  }
  
  #Matriz Hessiana considerando a reparametrização
  H1 <- function(aa,bb){
    s11 <- U1(aa,bb)[1] + exp(2*aa)*H(exp(aa),exp(bb))[1,1]
    s12 <- exp(aa+bb)*H(exp(aa),exp(bb))[1,2]
    s22 <- U1(aa,bb)[2] + exp(2*bb)*H(exp(aa),exp(bb))[2,2]            
    return(matrix(c(s11,s12,s12,s22),2,2))
  }
  
  err = 10
  j=0
  aa = numeric()
  bb= numeric()
  aa[1]<-theta0[1]
  bb[1]<-theta0[2]
  while(err>epsilon & j < int){
    j<-j+1
    Aux <- c(aa[j],bb[j]) - prop*solve(H1(aa[j],bb[j]),tol=1e-100)%*%U1(aa[j],bb[j])
    aa[j+1] <- Aux[1]
    bb[j+1] <- Aux[2]
    err = max(abs(c(aa[j+1]-aa[j],bb[j+1]-bb[j])))
    if(imprimir==TRUE){ print(err)}	}
  
  S <- list()
  S$err<- err
  S$Iteracoes <- j
  
  S$theta<- c(aa[j+1],bb[j+1])
  S$H<- H1(aa[j+1],bb[j+1])
  S$U<- U1(aa[j+1],bb[j+1])
  return(S)
}



#Encontrando as estimativas de MV para o exercício

fit = MV(x)
fit1 = MV1(x)

#Estimativas:
fit$theta
fit1$theta
exp(fit1$theta)

#Erro do processo iterativo erro = max{|t1{j+1} - t1{j}|,|t2{j+1} - t2{j}|\}
fit$erro
fit1$erro
#Matriz Hessiana calculada na estimativa de MV
fit$H

#Verificando se os autovalores são negativos
eigen(fit$H)
#Como todos os autovalores são estritamente negativos, temos um ponto de máximo!!


#Vetor Escore calculado na estimativa de MV
fit$U

t1.hat <- k$par[1]
t2.hat <- k$par[2]


f <- function (y) {
  t1.hat * t2.hat * y^(t2.hat-1) * exp(-t1.hat * y^(t2.hat))
}

hist(x, probability = T, breaks = 10, ylab = 'Densidade de Frequência', xlab = 'Dados observados', main = 'Histograma da amostra')
curve(f, add = T)

#######################################################################################


library(alabama) # Pacote de otimização não linear
library(numDeriv)

vero <- function(x,par){
  v <- n*log(par[1])+n*log(par[2])+(par[2]-1)*sum(log(x))-par[1]*sum(x^par[2])
  return(-v)
}


k <- optim(par=c(1,1),fn=vero,
      method="BFGS",x=x
)
