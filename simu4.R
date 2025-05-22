## simu 4 ##
## basis function
# rm(list=ls())
base_val = 1
phi_k <- function(t=NULL, params=NULL, Ti=1){
    base_val + params[1] + params[2]*sqrt(2)*sin(2*pi*t)
}

source("simuFun_simple.R")
# ntrain  = 500
#p = 20 
# mc.cores = 2
#sig = 0.25
b = 2
blocksize = 4
n = ntrain

## covariates needed for the hawkes function ##
# ncov = 2
# beta = matrix(c(0,0),ncol=ncov,nrow=p)

## covariance matrix

cor_v =  0.45
mat =  matrix(0, b*blocksize, b*blocksize)

for(i in 1:blocksize){
    for(j in i:blocksize){
        if(i==j){
            mat[1:b + (i-1)*b ,1:b + (j-1)*b] = diag(1,b,b) 
        }else if(j-i==1){
            mat[1:b + (i-1)*b ,1:b + (j-1)*b] = diag(cor_v,b,b)*rep(c(1,-1),len=b) 
        }else{}
    }
}
mat = (mat + t(mat))
diag(mat) = 1



#mat = (sig_vec%*%t(sig_vec) ) * mat

cov = lapply(1:(p/blocksize), function(i){
    mat
})
cov = as.matrix(bdiag(cov))
graph_true = get_norm(cov)
graph_true[graph_true > 0] = 1
#cov = sig^2 * solve(cov)
#cov = cov2cor(cov)
S_true_val = get_norm(cov2cor(solve(cov)))

#mu = rep(c(1,rep(0,b-1)), p)
mu = rep(0,b*p)
cov_mat = cov2cor(solve(cov))
sig_vec = sqrt(rep(c(2^2,2^2),p))
cov_mat =  (sig_vec%*%t(sig_vec) ) * cov_mat
param_mat = mvrnorm(n, mu=mu, Sigma=cov_mat)

#base_seq = seq(1,b*p , by=b)
#param_mat[,base_seq] = 2

data.list = mclapply(1:n, function(i){
    #print(i)
    #Xi = runif(ncov,-1,1)
    #xb = as.numeric(beta%*%Xi)
    obstime = 1
    TT = c(0, obstime)
    
    eventlist = lapply(1:p, function(pp){
        params = param_mat[i, (1:b) + (pp-1)*b]
        model = mpp(data=NULL, gif = fexp_gif, mark = list(NULL, NULL),
                    params= params, TT = TT, 
                    gmap=expression(params), mmap=NULL)
        
        res = NULL
        max.rate = optimize(bhz, interval=c(0,obstime),
                            params=params, maximum=TRUE)$objective
        max.rate = 1.1 *max.rate
        
        while(any(class(res)=="try-error") | any(is.null(res)) ){
            res = try({simulate(model,  max.rate = max.rate)},silent = TRUE)
            max.rate = max.rate*2
            #cat("increased max rate!\n")
            #cat(xb , "\n")
        }
        #print(max.rate)
        df = NULL
        if(!is.null(res$data)){
            df=data.frame(feature_id=rep(pp,length(res$data$time)),
                          time=res$data$time)
        }
        df
    })
    
    df = do.call(rbind, eventlist)
    #list(eventlist=eventlist,X=Xi,obstime=obstime)
    df = data.table(subject_num=i, df)
    df
}, mc.cores = ncores)

data_df = do.call(rbind, data.list)

