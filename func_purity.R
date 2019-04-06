

## NN <- N; MM_B <- M_B; nn <- n; SS <- S
fn.preprocessing <- function(NN, MM_B, nn, SS)
{
    if(sd(MM_B) > 0)
    {
        tmp <- lm(NN~MM_B)
        coef <- summary(tmp)$coeff[,1]
        
        N_hat <- coef[1] + MM_B*coef[2]  ## estimated N from the linear regression
        u.b <- summary(tmp)$sigma*3 + N_hat
        l.b <- -summary(tmp)$sigma*3 + N_hat
        
    }else{
        u.b <- sd(NN)*3 + mean(NN)
        l.b <- -sd(NN)*3 + mean(NN)
    }
    
    ind_remove <- (((NN > u.b)|(NN < l.b)))
    snv_set <- (1:SS)[ind_remove]  ## set of snv's removed
    
    if(length(snv_set) > 0)
    {
        NN <- as.matrix(NN[-snv_set,])
        MM_B <- as.matrix(MM_B[-snv_set,])
        nn <- as.matrix(nn[-snv_set,])
    }
    
    #print(paste("# of removed loci", length(snv_set)))
    
    return(list(N=NN, M_B=MM_B, n=nn, S=nrow(nn),removed_ind=((1:SS)%in%snv_set)))
}



##  NN <- NN[,i_t]; nn <- nn[,i_t]; #MM_B <- MM_B[,i_t]
fn.get.wstar <- function(NN, nn, MM_B)
{
    ##########################   1.  Prepare Data Set   #####################
    ###                                                                 ###
    ### Take the SNVs between the 40% percentile and 60% percentile of N's, and treat these SNVs as having copy number neutral. Then cluster n/N for these SNVs using DP.
    set.seed(919919)
  
    #plot_data(NN, nn, MM_B)
    #browser()
  
    nn2 <- jitter(nn)
    NN2 <- jitter(NN)
    nN2 <- nn2/NN2

    keep.id <- (nN2 < 1 & NN2 > quantile(NN2, 0.01))

    nN2 <- nN2[keep.id]
    MM_B2 <- MM_B[keep.id]
    
    if(sum(MM_B2==2) > 0)
    {
        selnN <- nN2[MM_B2==2]
    }else{
        N40 <- quantile(NN2, 0.40)
        N60 <- quantile(NN2, 0.60)
        selnN <- nN2[(NN2 > N40)& (NN2 < N60)]
    }
    ########################################################################
    
    ##########################   2.  Set up Dirichlet process mixture #####
    ###                                                                 ###
    nburn <- 500
    nsave <- 1000
    nskip <- 5
    ndisplay <- 500
    mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)
    
    priordp <- list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1),
    psiinv2=solve(diag(0.5,1)),
    nu1=4,nu2=4,tau1=1,tau2=100)  ### Used the prior4 setting in the DPdensity example, which gives noninformative priors
    
    fit <- DPdensity(y=selnN,prior=priordp,mcmc=mcmc,state=state,status=TRUE)  ### Fit an initial DP MCMC iteration
    
    #########    3. Trick DP to save clustering membership every 10 iterations for 100 samples  #######
    ### Note: DPpackage does not automatically save the MCMC samples for clustering memberships    ###
    
    state.save <- NULL
    input.state <- fit$state$ss
    
    #### Change start. March 5, 2016 YJ ########
    #n.sim <- 100
    n.sim <- 100 ### change to less to run quicker @subhajit
    mcmc1 <- list(nburn=100, nsave=10, nskip=1, ndisplay=10)  ##
    
    for(i in 1:n.sim){
        cat(i, " ")
        fit=DPdensity(y=selnN, state=input.state,prior=priordp, mcmc=mcmc1, status=TRUE)  ### Change end. March 5, 2016 YJ ###
        input.state=fit$state$ss
        state.save <- cbind(state.save, input.state)        
    } #for(i in 1:n.sim){
    
   
    list_data = NULL
    list_data$selnN = selnN
    list_data$state.save = state.save
    #save(list_data,file ="save.fit.4.RData")
    #### 4. For each MCMC sample, use the right cluster and middle but <0.5 cluster to estimate w-star. 
    ### The right cluster corresponds to homozygous somatic mutations, and middle but <0.5 cluster corresponds 
    ### to heterozygous somatic mutations.
    
    right.clust <- NULL
    mid.clust <- NULL
    right.size <- NULL
    mid.size <- NULL
    
    var=NULL
    var$selnN = selnN
    var$state.save = state.save
    save(var,file="./sim1_state_save.RData")
    
    stop("done")
    
    pur_arr = NULL
    
    for(i in 1:n.sim){
        curr.state <- state.save[,i]

         ### chaged subhajit
         #max.c <- max(curr.state)
         
        curr.mean = NULL
        curr.size = NULL
        curr.sd = NULL
        #for(j in 1:max.c){
        #    curr.mean <- c(curr.mean, mean(selnN[curr.state==j]))
        #}
        
        max.c <- max(curr.state)
        #U1 = unique(curr.state)
        for(j in 1:max.c){
            curr.mean <- c(curr.mean, mean(selnN[curr.state==j]))
            curr.size = c(curr.size,sum(curr.state==j))
            curr.sd = c(curr.sd,sd(selnN[curr.state==j]))
        }

        print(curr.mean)
        print(curr.size)
        ##### remove smaller clusters based on their size; if they are not at least 1%
        #DF = as.data.frame(table(curr.state))
        #browser()
        #DF1 = DF[order(DF[,2], decreasing = TRUE), ]
        ##I1=which(DF1[,2]>0.01*length(selnN))
        #loc_max = length(DF1[I1,2])
        #val_cluster = DF1[I1,1]
        #curr.mean <- NULL
        
        #for(j in 1:loc_max){
        #    curr.mean <- c(curr.mean, mean(selnN[curr.state==val_cluster[j]]))
        #}
       
        lower.range = 0.5
        upper.range = 1.0
        iter = 1
        
        perc1 = 0.1
        perc2 = 0.01
        
        while(lower.range > 1/16){
          
            print(paste0("iter ### ",iter))
            
            id1 = which(curr.mean > lower.range & curr.mean < upper.range)
            mean.arr = curr.mean[id1]
            max.mean1 = max(mean.arr)  ### location for homozygous if it is > 1% in size
            h = which(curr.mean==max.mean1)
            if(length(h) == 0)
                break
            if(iter == 1){ ### for the case where homozygous cluster > 0.5
                if(curr.size[h] >= perc2*length(selnN)){
                    m1 = max.mean1 ## homozygous mean
                    
                    l_limit_1 = max((m1/2) - (sd_clust/4),lower.range/2)
                    u_limit_1 = min((m1/2) + (sd_clust/4),upper.range/2)
                    
                    cluster.id = NULL
                    w = NULL
                    m1_by_2_arr = NULL
                    cnt1 = 0
                    for(j in 1:max.c){
                        if((curr.mean[j] > l_limit_1) & (curr.mean[j] < u_limit_1) &(j != h)){
                            cnt1 = cnt1+1
                            w[cnt1] = curr.size[j]
                            m1_by_2_arr[cnt1] = 2*curr.mean[j]
                        }
                    }
                    if(cnt1 == 0 & sum(w) < curr.size[h])
                        break
                    
                    cnt1 = cnt1+1
                    w[cnt1] = curr.size[h]
                    m1_by_2_arr[cnt1] = m1
                    
                    w1 = w/sum(w)
                     
                    pur_tmp = t(w)%*%m1_by_2_arr
                    pur_arr = c(pur_arr, pur_tmp)
                    
                    break
                }else
                    break
            }
            
            if(iter > 1){
                if(curr.size[h] >= perc1*length(selnN)){
                    pur_arr = c(pur_arr,2*max.mean1)    
                    break
                }        
            
                if(curr.size[h] >= perc2*length(selnN) & curr.size[h] < perc2*length(selnN)){
                
                    m1 = max.mean1
                    
                    l_limit_1 = max((m1/2) - (sd_clust/4),lower.range/2)
                    u_limit_1 = min((m1/2) + (sd_clust/4),upper.range/2)
                    
                    cluster.id = NULL
                    w = NULL
                    m1_by_2_arr = NULL
                    cnt1 = 0
                    for(j in 1:max.c){
                        if((curr.mean[j] > l_limit_1) & (curr.mean[j] < u_limit_1) &(j != h)){
                            cnt1 = cnt1+1
                            w[cnt1] = curr.size[j]
                            m1_by_2_arr[cnt1] = 2*curr.mean[j]
                        }
                    }
                    
                    if(cnt1 == 0 & sum(w) < curr.size[h])
                        break
                    
                    cnt1 = cnt1+1
                    w[cnt1] = curr.size[h]
                    m1_by_2_arr[cnt1] = m1
                    
                    w1 = w/sum(w)
                    
                    pur_tmp = t(w)%*%m1_by_2_arr
                    pur_arr = c(pur_arr, pur_tmp)
                     
                    break
                }else
                    break
            }
            
            if(curr.size[h] < perc2*length(selnN)){ ### cluster size is less than 1% of all the SNVs
                lower.range = lower.range/2
                upper.range = upper.range/2
                iter = iter+1
            }
        }
    }
    pur = median(pur_arr)
    return(pur)

#         n_cnt <- 2
#         mid_loc <- NA
#         
#         #while((n_cnt <= 8)&(is.na(mid_loc)))  ###  at least (0.25, 0.125)
#         while((n_cnt <= 16)&(is.na(mid_loc)))  ###  at least (0.125, 0.0625)
#         {
#             if(max.mean1 > (1.0/n_cnt))
#             {
#                 mid_loc <- 1.0/n_cnt
#             }else{
#                 n_cnt <- n_cnt*2.0
#             }
#         }
#         
#         if((!is.na(max.mean1))&(!is.na(mid_loc)))
#         {
#             right.ix <- which(curr.mean == max.mean1)
#             ## the size of the right cluster -- homozygous mutations
#             right.size <- c(right.size, length(selnN[curr.state %in% right.ix])) 
#             ## the mean of the right cluster as an estimate for 1-wstar
#             right.clust <- c(right.clust, max.mean1)  
#             
#             ind.tmp <- ((curr.mean - mid_loc)<=0)
#             
#             if(sum(ind.tmp) > 0)  ## nothing is below mid_loc
#             {
#                 max.mean2 <- max(curr.mean[ind.tmp])   ### location for heterozygous
#                 
#                 mid.ix <- which(curr.mean== max.mean2)
#                 mid.size <- c(mid.size, length(selnN[curr.state %in% mid.ix])) ## the size of the middle cluster -- heterozygous mutations
#                 mid.clust <- c(mid.clust, max.mean2) ## the mean of the middle cluster as an estimate for wstar/2
#             }else{  ### if(sum(ind.tmp) > 0)
#                 
#                 right.size[i] <- NA
#                 right.clust[i] <- NA
#                 mid.size <- c(mid.size, NA)
#                 mid.clust <- c(mid.clust, NA)
#             }  ###   if(sum(ind.tmp) > 0)
#         }else{
#             
#             right.size <- c(right.size, NA)
#             right.clust <- c(right.clust, NA)
#             mid.size <- c(mid.size, NA)
#             mid.clust <- c(mid.clust, NA)
#         } ##if((!is.na(max.mean1))&(!is.na(mid_loc)))
#         
#     } ##for(i in 1:n.sim){
#     
#     w <- mid.size/(mid.size+right.size)  ## w-star is estiamted as a weighted average of the means of the middle and right clusters, with weight proportional to the cluster size
#  
#     w.star <- w*mid.clust*2 + (1-w)*right.clust
#     print(paste0("Estimated Purity is ", median(w.star, na.rm = TRUE)))
#     
#     return(1-median(w.star, na.rm = TRUE)) ## return the posterior median
}

fn.get.wstar.new <- function(NN, nn, MM_B,a0,MIN_CNT)
{
    ##########################   1.  Prepare Data Set   #####################
    ###                                                                 ###
    ### Take the SNVs between the 40% percentile and 60% percentile of N's, and treat these SNVs as having copy number neutral. Then cluster n/N for these SNVs using DP.
    set.seed(919919)
    
    #plot_data(NN, nn, MM_B)
    
    nn2 <- jitter(nn)
    NN2 <- jitter(NN)
    nN2 <- nn2/NN2
    
    keep.id <- (nN2 < 1 & NN2 > quantile(NN2, 0.01))
    
    nN2 <- nN2[keep.id]
    MM_B2 <- MM_B[keep.id]
    
    if(sum(MM_B2==2) > 0)
    {
        selnN <- nN2[MM_B2==2]
    }else{
        N40 <- quantile(NN2, 0.40)
        N60 <- quantile(NN2, 0.60)
        selnN <- nN2[(NN2 > N40)& (NN2 < N60)]
    }
    ########################################################################
    
    ##########################   2.  Set up Dirichlet process mixture #####
    ###                                                                 ###
    nburn <- 500
    nsave <- 100
    nskip <- 5
    ndisplay <- 500
    mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)
    
    #priordp <- list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1),
    #                psiinv2=solve(diag(0.5,1)),
    #                nu1=4,nu2=4,tau1=1,tau2=100)  ### Used the prior4 setting in the DPdensity example, which gives noninformative priors
    
    priordp <- list(a0=a0,b0=1,m2=rep(0,1),s2=diag(100000,1),
                    psiinv2=solve(diag(0.5,1)),
                    nu1=4,nu2=4,tau1=1,tau2=100)  ### Used the prior4 setting in the DPdensity example, which gives noninformative priors
    
    
    fit <- DPdensity(y=selnN,prior=priordp,mcmc=mcmc,state=state,status=TRUE)  ### Fit an initial DP MCMC iteration
    
    #########    3. Trick DP to save clustering membership every 10 iterations for 100 samples  #######
    ### Note: DPpackage does not automatically save the MCMC samples for clustering memberships    ###
    
    state.save <- NULL
    input.state <- fit$state$ss
    
    #### Change start. March 5, 2016 YJ ########
    #n.sim <- 100
    n.sim <- 100 ### change to less to run quicker @subhajit
    mcmc1 <- list(nburn=100, nsave=10, nskip=1, ndisplay=10)  ##
    
    for(i in 1:n.sim){
        cat(i, " ")
        fit=DPdensity(y=selnN, state=input.state,prior=priordp, mcmc=mcmc1, status=TRUE)  ### Change end. March 5, 2016 YJ ###
        input.state=fit$state$ss
        state.save <- cbind(state.save, input.state)        
    } #for(i in 1:n.sim){
    
    
    list_data = NULL
    list_data$selnN = selnN
    list_data$state.save = state.save
    #save(list_data,file ="save.fit.4.RData")
    #### 4. For each MCMC sample, use the right cluster and middle but <0.5 cluster to estimate w-star. 
    ### The right cluster corresponds to homozygous somatic mutations, and middle but <0.5 cluster corresponds 
    ### to heterozygous somatic mutations.
    
    right.clust <- NULL
    mid.clust <- NULL
    right.size <- NULL
    mid.size <- NULL
    
    var=NULL
    #var$selnN = selnN
    #var$state.save = state.save
    #save(var,file="./sim14_state_save.RData")
    
    #load("./sim14_state_save.RData")
    #selnN = var$selnN
    #state.save = var$state.save
    #browser()
    #n.sim = 100
    pur_arr = NULL
    
    for(i in 1:n.sim){
        curr.state <- state.save[,i]
        
        ### chaged subhajit
        #max.c <- max(curr.state)
        
        curr.mean = NULL
        curr.size = NULL
        curr.sd = NULL
        
        max.c <- max(curr.state)
        #U1 = unique(curr.state)
        for(j in 1:max.c){
            curr.mean <- c(curr.mean, mean(selnN[curr.state==j]))
            curr.size = c(curr.size,sum(curr.state==j))
            curr.sd = c(curr.sd,sd(selnN[curr.state==j]))
        }
       
		### het cluster should be the largest
		htc_id =  which.max(curr.size)
		m_het = curr.mean[htc_id]
 
        print(curr.mean)
        print(curr.size)
        #a1 = which(curr.mean>0.5)
        #m1 = max(curr.mean[a1])
        #h1=which(curr.mean==m1)
        #print((curr.size[h1])/sum(curr.size) )
        #print(curr.size[h1])
        ##### remove smaller clusters based on their size; if they are not at least 1%
        #DF = as.data.frame(table(curr.state))
        #browser()
        #DF1 = DF[order(DF[,2], decreasing = TRUE), ]
        ##I1=which(DF1[,2]>0.01*length(selnN))
        #loc_max = length(DF1[I1,2])
        #val_cluster = DF1[I1,1]
        #curr.mean <- NULL
        
        #for(j in 1:loc_max){
        #    curr.mean <- c(curr.mean, mean(selnN[curr.state==val_cluster[j]]))
        #}
        
        lower.range = 0.5
        upper.range = 1.0
        iter = 1
        MIN_CUTOFF = MIN_CNT
        perc1 = 0.05
        perc2 = 0.005
        
        while(lower.range > 1/32){
            
            print(paste0("iter ### ",iter))
            
            id1 = which(curr.mean > lower.range & curr.mean < upper.range)
            #print(paste0("length id1 = ",length(id1)))
            if(length(id1) == 0){
                lower.range = lower.range/2
                upper.range = upper.range/2
                iter = iter+1
                #print("### inside if ###")
                next
            }
            
            #cat("###length id1 = ",(id1))    
            mean.arr = curr.mean[id1]
            #browser()
            max.mean1 = max(mean.arr)  ### location for homozygous if it is > 1% in size
            h = which(curr.mean==max.mean1)
            
            if(iter == 1){ ### for the case where homozygous cluster > 0.5
                if((curr.size[h] >= perc2*length(selnN)) | (curr.size[h] > MIN_CUTOFF) ){
                    m1 = max.mean1 ## homozygous mean
                    m1_by_2 = m1/2  ### heterozygous mean
					cat("###m1 ",m1)
                    sd_clust = curr.sd[h]
                    l_limit_1 = max((m1/2) - (sd_clust/4),lower.range/2)
                    u_limit_1 = min((m1/2) + (sd_clust/4),upper.range/2)
                    
                    cluster.id = NULL
                    w = NULL
                    m1_by_2_arr = NULL
                    cnt1 = 0
                    for(j in 1:max.c){
                        if((curr.mean[j] > l_limit_1) & (curr.mean[j] < u_limit_1) &(j != h)){
                            cnt1 = cnt1+1
                            w[cnt1] = curr.size[j]
                            m1_by_2_arr[cnt1] = 2*curr.mean[j]
                        }
                    }
                    #if(cnt1 == 0 & sum(w) < curr.size[h]){
                    if((cnt1 == 0) | (sum(w) < curr.size[h])){
                        lower.range = lower.range/2
                        upper.range = upper.range/2
                        iter = iter+1
                        next
                    }
                    
                    
                    cnt1 = cnt1+1
                    w[cnt1] = curr.size[h]
                    m1_by_2_arr[cnt1] = m1
                    
                    w1 = w/sum(w)
                    print(paste0("w1 = ",w1))
                    print(paste0("m1_by_2_arr = ",m1_by_2_arr))
                    pur_tmp = t(w1)%*%m1_by_2_arr
                    print(paste0("1)pur_tmp = ",pur_tmp))
                    pur_arr = c(pur_arr, pur_tmp)
                    break
                    
                }else{
                    lower.range = lower.range/2
                    upper.range = upper.range/2
                    iter = iter+1
                    next
                }
            }
            
            if(iter > 1){
                if(curr.size[h] >= perc1*length(selnN)){
                    pur_arr = c(pur_arr,2*max.mean1)    
                    break
                }        
                
                if((curr.size[h] >= perc2*length(selnN) & curr.size[h] < perc2*length(selnN)) | (curr.size[h] > MIN_CUTOFF) ){
                    #if((curr.size[h] >= perc2*length(selnN) & curr.size[h] < perc2*length(selnN)) ){    
                    m1 = max.mean1
                    sd_clust = curr.sd[h]
                    
                    l_limit_1 = max((m1/2) - (sd_clust/4),lower.range/2)
                    u_limit_1 = min((m1/2) + (sd_clust/4),upper.range/2)
                    
                    cluster.id = NULL
                    w = NULL
                    m1_by_2_arr = NULL
                    cnt1 = 0
                    for(j in 1:max.c){
                        if((curr.mean[j] > l_limit_1) & (curr.mean[j] < u_limit_1) &(j != h)){
                            cnt1 = cnt1+1
                            w[cnt1] = curr.size[j]
                            m1_by_2_arr[cnt1] = 2*curr.mean[j]
                        }
                    }
                    
                    if((cnt1 == 0) | (sum(w) < curr.size[h])){
                        lower.range = lower.range/2
                        upper.range = upper.range/2
                        iter = iter+1
                        next
                    }
                    
                    cnt1 = cnt1+1
                    w[cnt1] = curr.size[h]
                    m1_by_2_arr[cnt1] = m1
                    
                    w1 = w/sum(w)
                    
                    pur_tmp = t(w1)%*%m1_by_2_arr
                    pur_arr = c(pur_arr, pur_tmp)
                    
                    break
                }else{
                    lower.range = lower.range/2
                    upper.range = upper.range/2
                    iter = iter+1
                    next
                }
            }
            
            if(curr.size[h] < perc2*length(selnN)){ ### cluster size is less than 1% of all the SNVs
                lower.range = lower.range/2
                upper.range = upper.range/2
                iter = iter+1
                next
            }
        }
    }
    
    print(pur_arr)
    pur = median(pur_arr)
    cat("pur = ",pur,"\n")
    return(pur)

}

fn_get_purity <- function(N, n, M_B){
    ##########################   1.  Prepare Data Set   #####################
    ###                                                                 ###
    ### Take the SNVs between the 40% percentile and 60% percentile of N's, and treat these SNVs as having copy number neutral. Then cluster n/N for these SNVs using DP.
    set.seed(939919)
    
    n2 = jitter(n)
    N2 = jitter(N)
    VAF = n2/N2
    
    keep.id = (VAF < 1 & N2 > quantile(N2, 0.01))
    
    VAF = VAF[keep.id]
    M_B2 = M_B[keep.id]
    
    
    ## region of VAF: need to revise
    if(sum(M_B2 == 2) > 0){
        VAF = VAF[M_B2==2]
    }else{
        N40 = quantile(N2, 0.40)
        N60 = quantile(N2, 0.60)
        VAF = VAF[(N2 > N40)& (N2 < N60)]
    }
    
    ########################################################################
    ##########################   2.  Set up Dirichlet process mixture #####
    ###                                                                 ###
    
    nburn = 500
    nsave = 1000
    nskip = 10
    ndisplay = 500
    mcmc = list(nburn = nburn, nsave = nsave, nskip = nskip, ndisplay = ndisplay)
    
    #priordp <- list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1),
    #                psiinv2=solve(diag(0.5,1)),
    #                nu1=4,nu2=4,tau1=1,tau2=100)  ### Used the prior4 setting in the DPdensity example, which gives noninformative priors
    
    priordp = list(a0 = 2, b0 = 1, 
                    m2=rep(0,1),s2=diag(100000,1),
                    psiinv2=solve(diag(0.5,1)),
                    nu1=4,nu2=4,tau1=1,tau2=100)  ### Used the prior4 setting in the DPdensity example, which gives noninformative priors
    

    fit = DPdensity(y = VAF, prior = priordp, mcmc = mcmc, state = NULL, status = TRUE)  ### Fit an initial DP MCMC iteration
    
    #########    3. Trick DP to save clustering membership every 10 iterations for 100 samples  #######
    ### Note: DPpackage does not automatically save the MCMC samples for clustering memberships    ###
    
    
    # each row records the samples for iteration i. # of rows: nsave
    # # of columns: 2 * # of data points + 3
    n_data = length(VAF)
    
    param_matrix = fit$save.state$randsave[ , 1:(n_data*2)]

    
    mu_matrix = param_matrix[ , (1:n_data)*2 - 1]
    sigma_matrix = param_matrix[ , (1:n_data)*2]
    
    purity = rep(-1, nsave)
    
    for(i in 1:nsave){
        mu_table = c(table(mu_matrix[i, ]))
        mu_unique = as.numeric(names(mu_table))
        # names(mu_table) = NULL
        size_cluster = mu_table
        
        if(length(mu_unique[mu_unique < 0.5]) == 0){
            left_max = 0
            left_size = 0
        }else{
            left_max = max(mu_unique[mu_unique < 0.5])
            left_size = size_cluster[toString(left_max)]
        }
        
        
        if(length(mu_unique[mu_unique >= 0.5 & mu_unique < 1]) == 0){
            right_max = 0
            right_size = 0
        }else{
            right_max = max(mu_unique[mu_unique >= 0.5 & mu_unique < 1])
            right_size = size_cluster[toString(right_max)]
        }
        
        purity[i] = (left_size * 2 * left_max + right_size * left_max)/(left_size + right_size)        
    }
    
    
    
    return(median(purity, na.rm = TRUE)) ## return the posterior median
}



fn.elicit.w.star.B <- function(MM_B, NN, nn, SS)
{
    loci_set <- (1:SS)[(MM_B > 1.8)&(MM_B < 2.2)]
    #length(loci_set) > 0.3*SS
    p_set <- (nn/NN)[loci_set,]
    
    p_set1 <- p_set[(p_set < 0.5)]
    p_set1_cutoff <- quantile(p_set1, 0.8)
    w_star_1 <- 1 - 2*mean(p_set1[p_set1_cutoff < p_set1])
    
    
    p_set2 <- p_set[(p_set > 0.5)]
    p_set2_cutoff <- quantile(p_set2, 0.8)
    w_star_2 <- 1 - mean(p_set2[p_set2_cutoff < p_set2])
    
    w_star <- (length(p_set1)*w_star_1 + length(p_set2)*w_star_2)/(length(p_set1) + length(p_set2))
    
    return(w_star)
    
}



fn.elicit.w.star <- function(NN, nn, SS, TT)
{
    w_star <- rep(NA, TT)
    
    for(i_t in 1:TT)
    {
        N_t <- NN[,i_t]
        n_t <- nn[,i_t]
        
        m_N <- median(N_t)
        loci_set <- (1:SS)[(N_t > m_N*0.9)&(N_t < m_N*1.1)]
        #length(loci_set) > 0.3*SS
        p_set <- (n_t/N_t)[loci_set]
        
        
        p_set1 <- p_set[(p_set < 0.5)]
        p_set1_cutoff <- quantile(p_set1, 0.8)
        w_star_1 <- 1 - 2*mean(p_set1[p_set1_cutoff < p_set1])
        
        p_set2 <- p_set[(p_set > 0.5)]
        p_set2_cutoff <- quantile(p_set2, 0.8)
        w_star_2 <- 1 - mean(p_set2[p_set2_cutoff < p_set2])
        
        w_star[i_t] <- (length(p_set1)*w_star_1 + length(p_set2)*w_star_2)/(length(p_set1) + length(p_set2))
    }
    return(w_star)
}

fn_calc_purity <- function(NN, nn, SS, TT, MM_B, ind_LM, ind_DP,a0,MIN_CNT)
{
    if(sum(MM_B==2) > 0)
    {
      phi_tmp <- median(NN[MM_B==2])
    }else{
      phi_tmp <- median(NN)
    }
    
    for(i_t in 1:TT)
    {
        w_star_tmp <- rep(NA, 2)
        
        if(ind_LM == 1)
        {
            tmp <- summary(lm(NN[,i_t]~MM_B[,i_t]))$coeff[,1]
            w_star_tmp1 <- tmp[1]/phi_tmp
            w_star_tmp2 <- 1-2*tmp[2]/phi_tmp
            w_star_tmp[1] <- (w_star_tmp1 + w_star_tmp2)/2
            
            print(paste("w_star estimate from LM   ", w_star_tmp[1]))
        }
        
        if(ind_DP == 1)
        {
            ##  USE DP to estimate w_star
            #w_star_tmp[2] <- fn.get.wstar(NN[,i_t], nn[,i_t], MM_B[,i_t])
            w_star_tmp[2] <- 1-(fn.get.wstar.new(NN[,i_t], nn[,i_t], MM_B[,i_t],a0,MIN_CNT))
            
            print(paste0("(1-w_star) estimate from DP  ", 1-w_star_tmp[2]))
        }
        
        w_star_tmp <- mean(w_star_tmp, na.rm=TRUE)
        #print(c(i_t, w_star_tmp))
        
    }
    
    return(w_star_tmp)

}

plot_data <-function(N,n,M_B){
  plot(n/N,N)
}

