#' Combine independent posterior samples from multiple subsets
#' 
#' @param Samples A list containing posterior samples from each subset. Samples[[i]] contains a matrix of posterior samples for i-th set. Each line of the matrix is a sample and each column corresponds to a parameter.
#' @param num.sets The number of subsets. If not specified, the length of Samples will be used.
#' @param para.dim The number of parameters. If not specified, The number of columns of Samples[[1]] will be used. 
#' @param method Methods for combining the subset samples. Default is "reject" for rejection sampling. Alternative could be "importance" for importance sampling.
#' @param kernel The kernel function used for rejection sampling and importance sampling. The default is "uniform". Alternative could be "gaussian". Importance sampling always uses "gaussian" no matter what value is specified for kernel.
#' @param kernel.corr Indicator of whether a correlated kernel should be used for each subset. If TRUE, the algorithm will use a correlated kernels where the correlation equal to the sample correlation of posteriors on each subset. If FALSE, independent kernels are used. Default is TRUE.
#' @param accept The acceptance rate for rejection sampling and importance resampling. For more information please refer to the details. Default value is 0.1.
#' @param weight Indicator of whether subsets should be weighted. If TRUE, inverse variance weighting will be used. Depending on kernel.corr, the weights could be either matrices or scalars. Default is TRUE.
#' @param average Indicator of how samples are combined for initial processing. If TRUE, the subset samples will be (weighted) averaged before entering the rejection stage. If FALSE, mixture will be used instead, i.e., randomly (or weighted) select one sample among all subsets as the candidate for rejection sampling. The parameter is always set to TRUE for importance sampling. For more information please refer to the details. Default is TRUE.
#' @param matching Number of samples regenerated at each iteration (and final output). If TRUE, the algorithm will match the number of samples generated at each iteration to be the same as the number of input posterior samples (on one subset). If a number is given, the algorithm will use that number for each iteration. If FALSE, the algorithm will adjust its acceptance rate at each iteration to ensure the final number of combined samples equal to the input size times acceptance rate. For more information please refer to the details. Default is TRUE.
#' @param resampling Indicate whether resampling should be used for importance sampling. Default is TRUE.
#' @param resampling.corr If resampling is set to be TRUE, resampling.corr determines whether a joint kernel should be used for resampling, otherwise independent kernels are used. Default is FALSE.
#' @param display Indicate whether processing information should be printed. Default is TRUE.
#' 
#'@return A matrix contain combined posterior samples. Each row is a sample and each column corresponds to a parameter.
#' 
#' @details Weierstrass is the implementation of the Weierstrass rejection sampler proposed in the paper. Weierstrass sampler is a 'divide-conquer-combine' type parallel sampler, including three distinct samplers. Different from the refining sampler and the sequential rejection sampler, Weierstrass rejection sampler is a pure post-processing algorithm in the sense that it directly works on posterior samples obtained from multiple subsets, while the other two also relies on the subset sampling procedure and the likelihood. The algorithm applies rejection sampling and importance sampling on the formula described in Theorem 3 in the paper.
#'
#' To combine the subset posterior samples, the algorithm adopts the 'pairwise-combining' strategy, i.e, it will first combine the subset pairwisely to obtain half numbers of new subsets and then repeat the procedure until obtaining the final one. Such procedure, though, adding a logarithm factor to the original run time, achieves substantially better result. Based on this process, the parameter "matching" is introduced as for determining how many samples should be preserved (or regenerated) at each iteration. If "matching" is set to be FALSE, the algorithm will increase the acceptance rate at each iteration in order to obtain enough numbers of combined samples at the final output. If "matching" is set to be TRUE or any number, the algorithm will use the same acceptance rate at each iteration but regenerate samples to match the number required. For the final output, the number of combined samples will be equal to "matching" (if a number is assigned to "matching") or the number of input samples (If "matching" is equal to TRUE).
#'
#' Subset samples are required to be initially combined before entering the rejection sampling phase. Two different ways for initial combining are introduced. The first is (weighted) average(if "average" is TRUE, default) and the second is mixture (otherwise). Typically, weighted average achieves better accuracy, which is thereby recommended for most cases. However, there are exceptions. First, when the support of posterior samples are discontinuous. For example, integers or confined regions, in which case, any averaging procedure will fail. Mixture is a good alternative for these examples as it will always be in the correct region. Second, when the posterior is multimodal and each parameter is likely to stuck at one mode. For example, Gaussian mixture model without label switching. In this case, using mixture is likely to recover the multimodel shape of posteriors at the cost of losing (a certain level) accuracy. It is worth noting that when the parameter "accept" is set to be 1, the algorithm will combine the subset posteriors just by (weighted) averaging or mixture, which might serve for other interest.
#' 
#' "accept" is the only parameter that needs to be manually determined by the user. It determines the acceptance rate for rejection sampling (and importance sampling). The bandwidth will also be tuned automatically by the algorithm based on the value of acceptance rate. It is crucial taht "accept" have impacts on both the accuracy and the complexity of the algorithm. When "matching" is set to be TRUE, the run time of the algorithm will be proportional to 1/accept. When the "matching" is set to be FALSE, the number of actual obtained combined posterior samples is equal to the number of input posterior samples times "accept". In both cases, the value of "accept" cannot be set too small. In addition, even in the scenario when "matching" is TRUE, "accept" will also subtly affect the effective sample size. Therefoe, we recommend the value of "accept" to be chosen within (0.01, 0.1) or just set equal to 2/(number of subsets). The default is 0.1.
#' @seealso \code{\link{logitTest}}, \code{\link{BinTest}}
#' @examples
#' \dontrun{Two examples are provided as separate functions. \code{\link{logitTest}} is an example for logistic regression which depends on the package "BayesLogit" and \code{\link{BinTest}} is an example for binomial distributions.}
#' \dontrun{A simple example of linear regression is given as follows.}
#' library('MASS')
#' p = 5
#' n = 10000
#' m = 20
#' V = n/m
#' beta = matrix(c(1,-1,2,4,0),p,1)
#' X = matrix(rnorm(n*p), n, p)
#' err = matrix(rnorm(n),n,1)
#' Y = X%*%beta + err
#' piror.mu = matrix(rep(0,p),p,1)
#' prior.var = diag(rep(10,p),p,p)
#'
#' Samples = list()
#' for(i in 1:m){
#' Samples[[i]] = mvrnorm(n = 10000, mu = solve(t(X[(V*(i-1)+1):(V*i),])%*%X[(V*(i-1)+1):(V*i),] + solve(prior.var)/m)%*%t(X[(V*(i-1)+1):(V*i),])%*%Y[(V*(i-1)+1):(V*i),], Sigma = solve(t(X[(V*(i-1)+1):(V*i),])%*%X[(V*(i-1)+1):(V*i),] + solve(prior.var)/m))
#' }
#' true.posterior = mvrnorm(n = 10000, mu = solve(t(X)%*%X + solve(prior.var)/m)%*%t(X)%*%Y, Sigma = solve(t(X)%*%X + solve(prior.var)/m))
#' CombSample = weierstrass(Samples)
#'
#' par(mfrow = c(2,3))
#' for(i in 1:5){
#'    plot(density(true.posterior[,i]),col = 'red',main = paste('posterior for beta',toString(i)))
#'    lines(density(CombSample[,i]),col = 'blue')
#'    legend('topleft', c('True posterior', 'Weierstrass'), col = c('red','blue'), lty = c(1,1))
#' }
#' 

weierstrass <- function(Samples, num.sets = NULL, para.dim = NULL, method = 'reject', kernel = 'uniform', kernel.corr = TRUE, accept = 0.1, weight = TRUE, average = TRUE, matching = TRUE, resampling = TRUE, resampling.corr = FALSE, display = TRUE, level = 1, sets = 1){
    m = length(Samples)
    if(!is.null(num.sets) && m!=num.sets){
        cat('Number of sets not consistent')
        break
    }
    if(is.null(dim(Samples[[1]])))
        p = 1
    else
        p = dim(Samples[[1]])[2]
    
    if(!is.null(para.dim) && p!=para.dim){
        cat('Number of parameters not consistent')
        break
    }
    if(p == 1)
        kernel.corr = FALSE
    if(m == 1)
        return(Samples[[1]])
    else if(m == 2){
        X1 = matrix(Samples[[1]],ncol = p)
        X2 = matrix(Samples[[2]],ncol = p)
        N = min(dim(X1)[1], dim(X2)[1])

        X1 = matrix(X1[1:N,],N,p)
        X2 = matrix(X2[1:N,],N,p)
        if(matching == TRUE) matching = N
        
        output = matrix(0, N, p)
        X = NULL
        w1 = w = matrix(1/2, p, 2)
        sigmah = matrix(0, p, 1)
        for(j in 1:p){
            w1[j,] = c(1/var(X1[,j]), 1/var(X2[,j]))
            sigmah[j,1] =1/(1/var(X1[,j])+1/var(X2[,j]))
            w1[j,] = w1[j,]/sum(w1[j,])
        }
        if(weight == TRUE) w = w1
        else sigmah = matrix(rep(mean(sigmah),p),p,1)
        
        if(kernel.corr == TRUE){
            #print(dim(X1))
            invC1 = solve(cov(X1))
            invC2 = solve(cov(X2))
            newC = solve(invC1 + invC2)
            if(weight == TRUE){
                W1 = invC1%*%newC
                W2 = invC2%*%newC
            }else{
                W1 = W2 = diag(rep(1/2,p), p,p)
            }
        }
        
        if(method == 'reject'){
            if(matching != FALSE)
                iter = ceiling(matching/N/accept)
            else
                iter = 1
            for(k in 1:iter){
            permut = sample(1:N, N)
            X2 = matrix(X2[permut,],N,p)
            if(average == TRUE){
                if(kernel.corr == TRUE){
                    #print(dim(X1), dim(X2), dim(W1), dim(W2))
                    output = X1%*%W1 + X2%*%W2
                }
                else
                    output = X1*repmat(t(w[,1]),N,1) + X2*repmat(t(w[,2]),N,1)
            }
            else{
                judge = matrix(runif(N*p), N, p)<repmat(t(w[,1]),N,1)
                output = X2
                output[judge] = X1[judge]
            }

            if(kernel == 'uniform'){
                if(kernel.corr == TRUE)
                    diff = abs(X1 - X2)%*%(invC1+invC2)
                else
                    diff = abs(X1 - X2)/repmat(t(sigmah),N,1)
                max.diff = apply(diff, 1, max)
                cutoff = quantile(max.diff, prob = accept)
                keep = max.diff<=cutoff
                #print(length(keep))
            }

            if(kernel == 'gaussian'){
                if(kernel.corr == TRUE){
                    max.diff = (X1%*%W1%*%(invC1+invC2)*X1 + X2%*%W2%*%(invC1+invC2)*X2
                                - output%*%(invC1+invC2)*output)%*%matrix(1,p,1)
                }
                else{
                    diff = (X1^2*repmat(t(w[,1]),N,1) + X2^2*repmat(t(w[,2]),N,1)
                            - output^2)/repmat(t(sigmah^2),N,1)
                    max.diff = diff%*%matrix(1,p,1)
                }
                #find h
                minh = 0
                maxh = 1
                while(mean(exp(-max.diff/maxh))<accept){
                    minh = maxh
                    maxh = maxh*2
                }
                while(abs(maxh-minh)>0.01){
                    h = (minh + maxh)/2
                    if(mean(exp(-max.diff/h))<accept)
                        minh = h
                    else
                        maxh = h
                }
                h = (minh + maxh)/2
                #print(sum(exp(-max.diff/h)))
                keep = (runif(N) - exp(-max.diff/h)) <=0
            }

            X = rbind(X, matrix(output[keep,], sum(keep), p))
        }
        }

        if(method == 'importance'){
            if(matching != FALSE){
                iter = ceiling(matching/N/accept)
                accept0 = accept
            }
            else{
                iter = 1
                accept0 = 1
            }
            for(k in 1:iter){
            permut = sample(1:N, N)
            X2 = matrix(X2[permut,],N,p)
            if(kernel.corr == TRUE){
                output = X1%*%W1 + X2%*%W2
                max.diff = (X1%*%W1%*%(invC1+invC2)*X1 + X2%*%W2%*%(invC1+invC2)*X2
                            - output%*%(invC1+invC2)*output)%*%matrix(1,p,1)
            }else{
                output = X1*repmat(t(w[,1]),N,1) + X2*repmat(t(w[,2]),N,1)
                output2 = X1^2*repmat(t(w[,1]),N,1) + X2^2*repmat(t(w[,2]),N,1)
                diff = (output2 - output^2)/repmat(t(sigmah), N, 1)
                max.diff = diff%*%matrix(1,p,1)
            }
            
            if(k==1){
                #find h
                minh = 0
                maxh = 1
                
                while(mean(exp(-max.diff/maxh))<accept){
                    minh = maxh
                    maxh = maxh*2
                }
                while(abs(maxh-minh)>0.01){
                    h = (minh + maxh)/2
                    if(mean(exp(-max.diff/h))<accept)
                        minh = h
                    else
                        maxh = h
                }
            }
            w2 = exp(-max.diff/h)
            w2 = w2/sum(w2)
            keep = sample(1:N, size = ceiling(N*accept0), replace = TRUE, prob = w2)
            temp = matrix(output[keep,], length(keep), p)
            if(resampling == TRUE){
                require(MASS)
                if(resampling.corr == TRUE){
                    for(it in 1:length(keep)){
                        temp[it,] = mvrnorm(n=1,mu = temp[it,],Sigma = newC*sqrt(h)/8)
                    }
                }
                else{
                    temp = matrix(rnorm(p*length(keep), temp, repmat(t(sigmah)*sqrt(h)/8)), length(keep), p)
                }
            }
            X = rbind(X, temp)
        }
        }            
        return(X)
    }
    else if(m==3){
        Comb1 = weierstrass(Samples[1:2], num.sets = num.sets, para.dim = para.dim, method = method, kernel = kernel, kernel.corr = kernel.corr, accept = accept, weight = weight, average = average, matching = TRUE, resampling = resampling, resampling.corr = resampling.corr, display = display, level = level + 1)
        
        twoSamples = list()
        twoSamples[[1]] = Comb1
        twoSamples[[2]] = Samples[[3]]

        FinalComb = weierstrass(twoSamples, num.sets = num.sets, para.dim = para.dim, method = method, kernel = kernel, kernel.corr = kernel.corr, accept = accept, weight = weight, average = average, matching = matching, resampling = resampling, resampling.corr = resampling.corr, display = display, level = level)
        return(FinalComb)
    }
    else{
        if(matching == FALSE)
            accept1 = sqrt(accept)
        else
            accept1 = accept
        if(display == TRUE)
            cat(paste('starting Level',toString(level),', Set',toString(sets),'-',toString(sets-1+ceiling(m/2)),'\n'))
        Comb1 = weierstrass(Samples[1:ceiling(m/2)], num.sets = num.sets, para.dim = para.dim, method = method, kernel = kernel, kernel.corr = kernel.corr, accept = accept1, weight = weight, average = average, matching = matching, resampling = resampling, resampling.corr = resampling.corr, display = display, level = level + 1, sets = sets)
        
        if(display == TRUE)
            cat(paste('Level',toString(level),', Set',toString(sets),'-',toString(sets-1+ceiling(m/2)),'completed. Starting set', toString(sets-1+ceiling(m/2)+1), '-', toString(sets-1+m),'\n'))
        Comb2 = weierstrass(Samples[(ceiling(m/2)+1):m], num.sets = num.sets, para.dim = para.dim, method = method, kernel = kernel, kernel.corr = kernel.corr, accept = accept1, weight = weight, average = average, matching = matching, resampling = resampling, resampling.corr= resampling.corr, display = display, level = level + 1, sets = sets+ceiling(m/2))
        
        if(display == TRUE)
            cat(paste('Level',toString(level),', Set',toString(sets-1+ceiling(m/2)+1),'-',toString(sets-1+m),'completed. Combining two final sets.\n'))
        twoSamples = list()
        twoSamples[[1]] = Comb1
        twoSamples[[2]] = Comb2

        FinalComb = weierstrass(twoSamples, num.sets = num.sets, para.dim = para.dim, method = method, kernel = kernel, kernel.corr = kernel.corr, accept = accept1, weight = weight, average = average, matching = matching, resampling = resampling, resampling.corr = resampling.corr, display = display, level = level)
        if(display == TRUE)
            cat(paste('Level',toString(level),'completed, retrieving upper level.\n'))
        return(FinalComb)
    }
}


repmat <- function(mat, nrow = 1, ncol = 1){
    if(is.null(dim(mat))){
        cat('it is not a matrix')
        break
    }

    r = dim(mat)[1]
    c = dim(mat)[2]
    return(matrix(rep(mat, nrow*ncol), nrow = nrow*r, ncol = ncol*c))
}
