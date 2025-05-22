## simu functions ##
require(PtProcess)
require(parallel)
require(MASS)
require(Matrix)
require(data.table)

## calculate true graph
get_norm <- function(x=NULL){
    res = matrix(0,p,p)
    for(i in 1:p){
        for(j in 1:p){
            res[i,j] = sqrt(sum( x[1:b+ (i-1)*b, 1:b+ (j-1)*b]^2 ))
        }
    }
    res
}

## intensity function
# bhz <- function(t=NULL, params=NULL, Ti=1){
#     B = lapply(seq_along(params), function(k){
#         params[k] * phi_k(t=t, k=k, Ti=Ti)
#     })
#     tmp = exp(Reduce("+", B))
#     tmp[tmp>500] = 500
#     tmp
# }

bhz <- function(t=NULL, params=NULL, Ti=1){
    # B = lapply(seq_along(params), function(k){
    #     params[k] * phi_k(t=t, k=k, Ti=Ti)
    # })
    tmp = exp( phi_k(t=t, params=params, Ti=Ti) )
    tmp[tmp>500] = 500
    tmp
}


fexp_gif <- function(data=NULL, evalpts=NULL, params=NULL, TT = NA, 
                     tplus = FALSE){
    
    if(any(is.na(TT))) {
        if (is.vector(evalpts)) Time <- evalpts
        else Time <- evalpts[, "time"]
        ci = bhz(t=Time, params=params)
    }
    else {
        
        ci = integrate(bhz, TT[1], TT[2], params=params)
        if (ci$message == "OK")
            ci <- ci$value
        else
            stop(paste("Problems with Numerical Integration: ", ci$message))
    }
    
    ci <- as.vector(ci)
    return(ci)
}

attr(fexp_gif, "rate") = "bounded"


gmap <- expression(params)
mmap <- expression(params)
