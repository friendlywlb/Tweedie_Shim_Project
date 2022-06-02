Q_r<-function(y,mu,power,lamb_r,lamb_b,r,b){
    den<-tweedie.dev(y=y,mu=mu,power=power)
    Q<-sum(den,na.rm=T)+sum(abs(r))*lamb_r+sum(abs(b[-1]))*lamb_b
    return(Q)
}

beta_comb<-function(beta){
  n=length(beta)
  if (n>2){
    index=combn(c(1:n),2)
    count=1
    comb_result=c()
    while (count<=ncol(index)){
      comb_result=c(comb_result,beta[index[1,count]]*beta[index[2,count]])
      count=count+1
    }
    return(c(beta,comb_result))}
  else{
    index=combn(c(1:n),2)
    comb_result=c()
    comb_result=c(comb_result,beta[index[1,1]]*beta[index[2,1]])
    return(c(beta,comb_result))
  }
}

beta_update <- function(beta,r,power,y,X,lambda){
    beta_x<-beta[-1]
    n_beta = length(beta_x)
    beta_update_result<-rep(0,n_beta)
    X=as.data.frame(X)
    non_zero_beta<-which(beta_x!=0)
    for (i in non_zero_beta){
        X_i<-X[,-i]
        X_i_m <- model.matrix(~.^2-1, data=X_i)
        beta_i<-beta_x[-i]
        index_r<-which(apply(combn(c(1:n_beta),m=2), 2,function(x) i %in% x))            
        coefficient_r<-r[index_r]*beta_i
        index_r1<-which(apply(combn(c(1:n_beta),m=2), 2,function(x) i %!in% x))
        coefff<-r_beta2(beta_i,r[index_r1])  
        tmp_offset<-X_i_m%*%coefff
        X_matrix_0<-as.matrix(X_i*X[,i])
        b_matrix_0<-X[,i]+X_matrix_0%*%coefficient_r
        b_matrix_1<- as.matrix(cbind(0, b_matrix_0))
        b_fit<-glmnet(x=b_matrix_1,y=y,family=tweedie(link.power=0,var.power=power)
                      ,offset=tmp_offset,intercept=FALSE,lambda=lambda,control=list(maxit=500))
        b_update_tmp<-coef(b_fit)[3,1]
        beta_update_result[i]<-b_update_tmp
        beta_x[i]<-b_update_tmp
            }
    b_0_x<-model.matrix(~.^2-1,data=as.data.frame(X))    
    coeffff<-r_beta2(beta_x,r)
    b_0_offset<-b_0_x%*%coeffff
    b0_fit<-glm(y~1,offset=b_0_offset,family=tweedie(link.power=0,var.power=power),control=list(maxit=500))
    b0<-b0_fit$coefficients
    beta_update_result<-c(b0,beta_update_result)                                
    return(beta_update_result)
}
                             
r_update<-function(beta,power,y,X,lambda){
    offset=as.matrix(cbind(1,X))%*%beta
    beta_x<-beta[-1]
    n<-length(beta_x)
    design_matrix<-c()
    index<-combn(c(1:n),m=2)
    for (i in 1:ncol(index)){
        beta_tmp<-beta_x[index[1,i]]*beta_x[index[2,i]]
        design_matrix<-cbind(design_matrix,beta_tmp*X[,index[1,i]]*X[,index[2,i]]  )
    }
    design_matrix<-as.matrix(design_matrix)
    r_fit<-glmnet(x=design_matrix,y=y,family=tweedie(link.power=0,var.power=power),offset=offset,intercept=FALSE,lambda=lambda,control=list(maxit=500))
    return(coef(r_fit)[-1,1])
}                   
 
                             
                             
shim_lasso_update_2<-function(x,y,power,lamb_r,lamb_b,initial_r,initial_b){
    initial_per<-999
    tmp_r<-initial_r
    tmp_b<-initial_b
    n<-dim(x)[2]
    index<-combn(c(1:n),m=2)
    n_c=0
    while(initial_per>=0.01){
        total_effect=0
        for (i in 1:ncol(index)){
            tmp_effect<-tmp_r[i]*tmp_b[index[1,i]+1]*tmp_b[index[2,i]+1]*x[,index[1,i]]*x[,index[2,i]]
            total_effect=total_effect+tmp_effect
        }
        tmp_fitted<-exp(as.matrix(x)%*%tmp_b[-1]+total_effect+tmp_b[1])
        Q_1<-Q_r(y,tmp_fitted,power,lamb_r,lamb_b,tmp_r,tmp_b)
        design_matrix<-c()
        for (i in 1:ncol(index)){
            beta_tmp<-tmp_b[index[1,i]+1]*tmp_b[index[2,i]+1]
            design_matrix<-cbind(design_matrix,beta_tmp*x[,index[1,i]]*x[,index[2,i]]  )
        }
        if (all(design_matrix==0)){
            result=list('beta'=tmp_b,'r'=tmp_r,'n'=rep(n_c,length(tmp_b)))
            return(result)
            break
        }
        else{
            tmp_r<-r_update(tmp_b,power=power,y=y,X=x,lambda=lamb_r)
            if (all(tmp_r==0)){
                tmp_fit<-glmnet(x=as.matrix(x),y=y,family=tweedie(link.power=0,var.power=power),control=list(maxit=500),lambda=lamb_b)
                tmp_b<-coef(tmp_fit)[,1]
                result=list('beta'=tmp_b,'r'=tmp_r,'n'=rep(n_c,length(tmp_b)))
            return(result)
            break}
            else{
                tmp_b<-beta_update(tmp_b,tmp_r,power=power,y=y,X=x,lambda=lamb_b) 
                total_effect=0
                for (i in 1:ncol(index)){
                    tmp_effect<-tmp_r[i]*tmp_b[index[1,i]+1]*tmp_b[index[2,i]+1]*x[,index[1,i]]*x[,index[2,i]]
                    total_effect=total_effect+tmp_effect
                }
                tmp_fitted<-exp(as.matrix(x)%*%tmp_b[-1]+total_effect+tmp_b[1])
                Q_2<-Q_r(y,tmp_fitted,power,lamb_r,lamb_b,tmp_r,tmp_b)
                initial_per<-abs(Q_1-Q_2)/abs(Q_1)
                n_c=n_c+1}
    }
    result=list('beta'=tmp_b,'r'=tmp_r,'n'=rep(n_c,length(tmp_b)))
    return(result)}
    }
         
                         
ts_fitted<-function(x,beta,r){
    n<-length(beta[-1])
    total_effect=0
    index<-combn(c(1:n),m=2)
    for (i in 1:ncol(index)){
        tmp_effect<-r[i]*beta[index[1,i]+1]*beta[index[2,i]+1]*x[,index[1,i]]*x[,index[2,i]]
        total_effect=total_effect+tmp_effect
    }
    tmp_fitted<-exp(as.matrix(x)%*%beta[-1]+total_effect+beta[1])
    return(tmp_fitted)
    
}
                             
                             
r_beta<-function(beta,r){
    n<-length(beta[-1])
    total_effect=c()
    index<-combn(c(1:n),m=2)
    for (i in 1:ncol(index)){
        tmp_effect<-r[i]*beta[index[1,i]+1]*beta[index[2,i]+1]
        total_effect=c(total_effect,tmp_effect)
    }
    final_beta<-c(beta,total_effect)
    return(final_beta)
    
}
                
                             
                             
        
c_x<- function(x,index,beta){
    index<-as.matrix(index)
    result=0
    for (i in 1:ncol(index)){
        result=result+beta[index[1,i]]*beta[index[2,i]]*x[,index[1,i]]*x[,index[2,i]]
    }
    return(result)
}
                                  
mysd <- function(z) sqrt(sum((z-mean(z))^2)/(length(z)-1))
                                  
`%!in%` <- Negate(`%in%`)
                                  
                                                               
r_beta2<-function(beta,r){
    n<-length(beta)
    total_effect=c()
    index<-combn(c(1:n),m=2)
    for (i in 1:ncol(index)){
        tmp_effect<-r[i]*beta[index[1,i]]*beta[index[2,i]]
        total_effect=c(total_effect,tmp_effect)
    }
    final_beta<-c(beta,total_effect)
    return(final_beta)
    
}
    