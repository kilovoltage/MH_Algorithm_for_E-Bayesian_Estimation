#alpha BS#
options(digits=5)
lambda0<-4.36534
a<-0.6119
b<-0.15234
h<-1.5
c<-1
beta<-1
lambda=1
h<-1.5
t=0.07
ns<-c(50,80,120)
rs<-c(40,60,90)
ks<-c(30,40,60)
Ts<-c(0.2,0.4)
lambdahat1<-matrix(numeric(48),6,8)
kk=1

for(ii in 1:3){
  n<-ns[ii]
  r<-rs[ii]
  k<-ks[ii]
  
  for(jj in 1:2){
    T<-Ts[jj]
    
    A<-matrix(numeric(80000),10000,8)
    for(j in 1:10000){
      x<-numeric(n)
      U<-runif(n,0,1)
      x<-(log(1-log(1-U)/lambda0))^beta
      x<-x[order(x)]
      D<-sum(x<T)
      
      D<-sum(x<T)
      if(D<k){
        l=k
      }
      if(D>=k&D<r){
        l=D
      }
      if(D>=r){
        l=r
      }
      
      #这里换成咱们昨天写的M的表达#
      S1<-0
      if(D<k){
        for(i in 1:k){S1<-S1-(1-exp(x[i]^beta))} 
        M<-((n-k)*(-1+exp(x[k]^beta))+S1)
      }
      if(D>=k&D<r){
        for(i in 1:D){S1<-S1-(1-exp(x[i]^beta))} 
        M<-((n-D)*(-1+exp(T^beta))+S1)
      }
      if(D>=r){
        for(i in 1:r){S1<-S1-(1-exp(x[i]^beta))} 
        M<-((n-r)*(-1+exp(x[r]^beta))+S1)
      }
      
      J<-function(b){log((b+M)/(h+b+M))}
      JJ<-function(b){(log((b+M)/(h+b+M)))^2}
      J1<-function(b){log((b+M)/(h+b+M))/(b+M)}
      J2<-function(b){1/(b+M)^2}
      E<-function(a){(a+l+1)*(a+l)}
      E1<-function(a){(a+l)^2}
      
      A[j,1]<-(a+l)/(b+M)
      A[j,2]<-(a+l)/(b+M)^2
      A[j,3]<-((2*l+1)/(2*c))*log((M+c)/M)
      A[j,4]<-(2*l+1)/(2*M*(M+c))
      A[j,5]<--(a+l)/h*log((b+M)/(h+b+M))
      A[j,6]<-(a+l+1)*(a+l)/(b+M)^2+2*(a+l)*(a+l)/(h*(b+M))*log((b+M)/(h+b+M))+(a+l)^2/h^2*(log((b+M)/(h+b+M)))^2
      A[j,7]<--(2*l+1)/(2*h*c)*integrate(J,0,c)$value
      A[j,8]<-(integrate(E,0,1)$value*integrate(J2,0,c)$value+2/h*integrate(E1,0,1)$value*integrate(J1,0,c)$value+1/h^2*integrate(E1,0,1)$value*integrate(JJ,0,c)$value)/c
    }
    
    lambdahat1[kk,1]<-mean(A[,1])
    lambdahat1[kk,2]<-mean(A[,2])
    lambdahat1[kk,3]<-mean(A[,3])
    lambdahat1[kk,4]<-mean(A[,4])
    lambdahat1[kk,5]<-mean(A[,5])
    lambdahat1[kk,6]<-mean(A[,6])
    lambdahat1[kk,7]<-mean(A[,7])
    lambdahat1[kk,8]<-mean(A[,8])
    
    kk=kk+1
  }
}

lab<-matrix(numeric(84),7,12)
lab2<-matrix(numeric(84),7,12)
lab3<-matrix(numeric(84),7,12)
lab4<-matrix(numeric(84),7,12)

kk<-1

for(ii in 1:3){
  n<-ns[ii]
  r<-rs[ii]
  k<-ks[ii]
  x<-numeric(n)
  U<-seq(1/(n+1),1-1/(n+1),1/(n+1)) 
  x<-(log(1-log(1-U)/lambda0))^beta
  x<-x[order(x)]
  
  for(jj in 1:2){
    T<-Ts[jj]
    D<-sum(x<T)
    
    
    if(D<k){
      g<-function(t){S<-1
      for(i in 1:k){S<-S*((1+x[i]/lambda)^(-t-1))} 
      p<-(t^(k+a-1))*exp(-b*t)*(((1+x[k]/lambda)^(-t))^(n-k))*S
      return(p)
      }
    }
    if(D>=k&D<r){
      g<-function(t){S<-1
      for(i in 1:D){S<-S*((1+x[i]/lambda)^(-t-1))} 
      p<-(t^(D+a-1))*exp(-b*t)*(((1+T/lambda)^(-t))^(n-k))*S
      return(p)
      }
    }
    if(D>=r){
      g<-function(t){S<-1
      for(i in 1:r){S<-S*((1+x[i]/lambda)^(-t-1))} 
      p<-(t^(r+a-1))*exp(-b*t)*(((1+x[r]/lambda)^(-t))^(n-k))*S
      return(p)
      }
    }
    
    if(D<k){
      l=k
    }
    if(D>=k&D<r){
      l=D
    }
    if(D>=r){
      l=r
    }
    
    S1<-0
    if(D<k){
      for(i in 1:k){S1<-S1+log(1+x[i]/lambda)} 
      M<-((n-k)*log(1+x[k]/lambda)+S1)
    }
    if(D>=k&D<r){
      for(i in 1:D){S1<-S1+log(1+x[i]/lambda)} 
      M<-((n-D)*log(1+T/lambda)+S1)
    }
    if(D>=r){
      for(i in 1:r){S1<-S1+log(1+x[i]/lambda)}
      M<-((n-r)*log(1+x[r]/lambda)+S1)
    }
    
    alphamle<- lambdahat1[kk,1]
    MSE<-lambdahat1[kk,2]
    alphamle1<-lambdahat1[kk,3]
    MSE1<-lambdahat1[kk,4]
    
    al<-numeric(11001)
    al[1]<-alphamle
    als<-abs(rnorm(11001,alphamle,MSE))
    u<-runif(11001,0,1)
    r1<-numeric(11001)
    for(j in 2:11001){
      r1[j]<-min(1,g(als[j])*dnorm(al[j-1],alphamle,MSE)/(g(al[j-1]*dnorm(als[j],alphamle,MSE))))
      if(u[j]-r1[j]<0){alpha<-als[j]
      al[j]<-alpha}
      else{al[j]<-al[j-1]}
    }
    alphahat<-numeric(10000)
    for(i in 1:10000){alphahat[i]<-al[i+1001]}
    
    
    al<-numeric(11001)
    al[1]<-alphamle1
    als<-abs(rnorm(11001,alphamle1,MSE1))
    u<-runif(11001,0,1)
    r1<-numeric(11001)
    for(j in 2:11001){
      r1[j]<-min(1,g(als[j])*dnorm(al[j-1],alphamle1,MSE1)/(g(al[j-1]*dnorm(als[j],alphamle1,MSE1))))
      if(u[j]-r1[j]<0){alpha<-als[j]
      al[j]<-alpha}
      else{al[j]<-al[j-1]}
    }
    
    ealphahat<-numeric(10000)
    for(i in 1:10000){ealphahat[i]<-al[i+1001]}
    
    
    lab[1,kk]<-n
    lab[2,kk]<-T
    lab[3,kk]<-mean(alphahat)
    lab[4,kk]<-var(alphahat)
    lab[5,kk]<-min(alphahat)
    lab[6,kk]<-max(alphahat)
    lab[7,kk]<-max(alphahat)-min(alphahat)
    
    lab[1,kk+6]<-n
    lab[2,kk+6]<-T
    lab[3,kk+6]<-mean(ealphahat)
    lab[4,kk+6]<-var(ealphahat)
    lab[5,kk+6]<-min(ealphahat)
    lab[6,kk+6]<-max(ealphahat)
    lab[7,kk+6]<-max(ealphahat)-min(ealphahat)
    
    RtBS<-numeric(10000)
    RtEBS<-numeric(10000)
    for(i in 1:10000){
      RtBS[i]<-(1+t/lambda)^(-alphahat[i])
      RtEBS[i]<-(1+t/lambda)^(-ealphahat[i])
    }
    
    lab2[1,kk]<-n
    lab2[2,kk]<-T
    lab2[3,kk]<-mean(RtBS)
    lab2[4,kk]<-var(RtBS)
    lab2[5,kk]<-min(RtBS)
    lab2[6,kk]<-max(RtBS)
    lab2[7,kk]<-max(RtBS)-min(RtBS)
    
    lab2[1,kk+6]<-n
    lab2[2,kk+6]<-T
    lab2[3,kk+6]<-mean(RtEBS)
    lab2[4,kk+6]<-var(RtEBS)
    lab2[5,kk+6]<-min(RtEBS)
    lab2[6,kk+6]<-max(RtEBS)
    lab2[7,kk+6]<-max(RtEBS)-min(RtEBS)
    
    kk<-kk+1
    
  }
}

kk<-1

for(ii in 1:3){
  n<-ns[ii]
  r<-rs[ii]
  k<-ks[ii]
  x<-numeric(n)
  U<-seq(1/(n+1),1-1/(n+1),1/(n+1)) 
  x<-(log(1-log(1-U)/lambda0))^beta
  x<-x[order(x)]
  
  for(jj in 1:2){
    T<-Ts[jj]
    D<-sum(x<T)
    
    
    if(D<k){
      g<-function(t){S<-1
      for(i in 1:k){S<-S*((1+x[i]/lambda)^(-t-1))} 
      p<-(t^(k+a-1))*exp(-b*t)*(((1+x[k]/lambda)^(-t))^(n-k))*S
      return(p)
      }
    }
    if(D>=k&D<r){
      g<-function(t){S<-1
      for(i in 1:D){S<-S*((1+x[i]/lambda)^(-t-1))} 
      p<-(t^(D+a-1))*exp(-b*t)*(((1+T/lambda)^(-t))^(n-k))*S
      return(p)
      }
    }
    if(D>=r){
      g<-function(t){S<-1
      for(i in 1:r){S<-S*((1+x[i]/lambda)^(-t-1))} 
      p<-(t^(r+a-1))*exp(-b*t)*(((1+x[r]/lambda)^(-t))^(n-k))*S
      return(p)
      }
    }
    
    if(D<k){
      l=k
    }
    if(D>=k&D<r){
      l=D
    }
    if(D>=r){
      l=r
    }
    
    S1<-0
    if(D<k){
      for(i in 1:k){S1<-S1+log(1+x[i]/lambda)} 
      M<-((n-k)*log(1+x[k]/lambda)+S1)
    }
    if(D>=k&D<r){
      for(i in 1:D){S1<-S1+log(1+x[i]/lambda)} 
      M<-((n-D)*log(1+T/lambda)+S1)
    }
    if(D>=r){
      for(i in 1:r){S1<-S1+log(1+x[i]/lambda)}
      M<-((n-r)*log(1+x[r]/lambda)+S1)
    }
    
    alphamle<- lambdahat1[kk,5]
    MSE<-lambdahat1[kk,6]
    alphamle1<-lambdahat1[kk,7]
    MSE1<-lambdahat1[kk,8]
    
    al<-numeric(11001)
    al[1]<-alphamle
    als<-abs(rnorm(11001,alphamle,MSE))
    u<-runif(11001,0,1)
    r1<-numeric(11001)
    for(j in 2:11001){
      r1[j]<-min(1,g(als[j])*dnorm(al[j-1],alphamle,MSE)/(g(al[j-1]*dnorm(als[j],alphamle,MSE))))
      if(u[j]-r1[j]<0){alpha<-als[j]
      al[j]<-alpha}
      else{al[j]<-al[j-1]}
    }
    alphahat<-numeric(10000)
    for(i in 1:10000){alphahat[i]<-al[i+1001]}
    
    
    al<-numeric(11001)
    al[1]<-alphamle1
    als<-abs(rnorm(11001,alphamle1,MSE1))
    u<-runif(11001,0,1)
    r1<-numeric(11001)
    for(j in 2:11001){
      r1[j]<-min(1,g(als[j])*dnorm(al[j-1],alphamle1,MSE1)/(g(al[j-1]*dnorm(als[j],alphamle1,MSE1))))
      if(u[j]-r1[j]<0){alpha<-als[j]
      al[j]<-alpha}
      else{al[j]<-al[j-1]}
    }
    
    ealphahat<-numeric(10000)
    for(i in 1:10000){ealphahat[i]<-al[i+1001]}
    
    
    lab3[1,kk]<-n
    lab3[2,kk]<-T
    lab3[3,kk]<-mean(alphahat)
    lab3[4,kk]<-var(alphahat)
    lab3[5,kk]<-min(alphahat)
    lab3[6,kk]<-max(alphahat)
    lab3[7,kk]<-max(alphahat)-min(alphahat)
    
    lab3[1,kk+6]<-n
    lab3[2,kk+6]<-T
    lab3[3,kk+6]<-mean(ealphahat)
    lab3[4,kk+6]<-var(ealphahat)
    lab3[5,kk+6]<-min(ealphahat)
    lab3[6,kk+6]<-max(ealphahat)
    lab3[7,kk+6]<-max(ealphahat)-min(ealphahat)
    
    RtBS<-numeric(10000)
    RtEBS<-numeric(10000)
    for(i in 1:10000){
      RtBS[i]<-(1+t/lambda)^(-alphahat[i])
      RtEBS[i]<-(1+t/lambda)^(-ealphahat[i])
    }
    
    lab4[1,kk]<-n
    lab4[2,kk]<-T
    lab4[3,kk]<-mean(RtBS)
    lab4[4,kk]<-var(RtBS)
    lab4[5,kk]<-min(RtBS)
    lab4[6,kk]<-max(RtBS)
    lab4[7,kk]<-max(RtBS)-min(RtBS)
    
    lab4[1,kk+6]<-n
    lab4[2,kk+6]<-T
    lab4[3,kk+6]<-mean(RtEBS)
    lab4[4,kk+6]<-var(RtEBS)
    lab4[5,kk+6]<-min(RtEBS)
    lab4[6,kk+6]<-max(RtEBS)
    lab4[7,kk+6]<-max(RtEBS)-min(RtEBS)
    
    kk<-kk+1
    
  }
}

data<-rbind(lab,lab2,lab3,lab4)

write.csv(data,file = "D:\\R and CRIs.csv",row.names = F)


