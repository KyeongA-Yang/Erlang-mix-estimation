
  library(pracma)
  library(gtools)

  #finding initial value
  Erlang_initial = function(x, breaks=8, shape=NULL){
    x <- sort(x)
    n <- length(x)
    theta <- (x[ceiling(0.9*n)] - x[floor(0.1*n)])/breaks
    
    if(!is.null(shape) && length(shape) > 0){  
      M = length(shape)
      weight <- rep(NA, M)
      for (i in 1:M){
        weight[i] = length(x[x <= i*theta & x > (i-1)*theta]) / length(x)
      }
      weight[weight == 0] = 0.1  
      weight = weight / sum(weight)  
    
    } else {  
      M = ceiling(x[n] / theta)
      shape = 1:M
      weight <- rep(NA, M)
      for (i in 1:M){
        weight[i] = length(x[x <= i*theta & x > (i-1)*theta]) / length(x)
      }
      shape = shape[weight > 0] 
      weight = weight[weight > 0] 
    }
    
    return(list(shape = shape, weight = weight, theta = theta))
  }

  
  #calculate log-like
  Erlang_lik=function(x, support, weight, theta){
    mix=Erlang_mix(x,support, weight, theta)
    log_lik = sum(log(mix))
    
    if(is.infinite(log_lik)){
      log_lik= Erlang_mix_loglik(x, support, weight, theta)
    }
    return(log_lik)
  }
  
  Erlang_mix_loglik <- function(x, support, weight, theta) {
    n <- length(x)
    p <- length(support)
    log_mix_lik <- 0
    
    for (i in 1:n) {
      log_probs <- rep(0, p)
      
      for (j in 1:p) {
        log_probs[j] <- log(weight[j]) + dgamma(x[i], shape = support[j], scale = theta, log = TRUE)
      }
      
      max_log_prob <- max(log_probs)
      log_sum_exp <- max_log_prob + log(sum(exp(log_probs - max_log_prob)))
      log_mix_lik <- log_mix_lik + log_sum_exp
    }
    return(log_mix_lik)
  }
  
  
  
  #algorithm for updating weight and scale
  update_Erlang=function(x,support,weight,theta,update_theta=TRUE){
    n=length(x)
    p=length(support)
    lik=Erlang_lik(x, support, weight, theta)
    
    max_x = max(x)
    new_weight = weight
    new_theta = theta
    for (itr in 1:200){
      if(p==1){
        weight=1
        if(update_theta){
          theta = mean(x)/sum(support*new_weight)
        }
        break
      }else{
        S=compute_S(x,support,weight,theta)
        
        new_weight=tryCatch(lsqlincon(S,matrix(2,n,1),NULL,NULL,matrix(1,1,p),1,lb=0*t(weight)+1.0e-14,ub=t(support/support)-1.0e-14), error=function(err){Erlang_EM(x,support,weight,theta,itr_max = 5,update_theta=FALSE)$weight}, finally = NULL)
        if(any(new_weight<0)){new_weight[new_weight<0]=1e-14}
        new_weight= new_weight/sum(new_weight)
        if(update_theta){
          new_theta = mean(x)/(sum(support*new_weight))
          new_theta =max(0.05, new_theta)
        }
        
        new_lik=Erlang_lik(x,support, new_weight, new_theta)    
        
        if(new_lik<lik){
          em_result = Erlang_EM(x,support, weight,theta,itr_max = 10,update_theta)
          
          new_weight=em_result$weight
          new_theta=em_result$theta
          new_theta =max(0.05, new_theta)
          
        }
        
        if(itr>=2 && (new_lik-lik)<1.0e-6){
          break
        }
        weight=new_weight
        theta=new_theta
        lik=new_lik
      }
    }
    return(list(support=support,weight=weight, theta=theta,lik=lik))
  }
  
  ## Main function that returns a set of initial values with fixed theta
  find_initial=function(x, r=1, weight=1, theta, comp, r_range=1:100,max_num_shape=round(1.5*comp), ...){
    small_comp = FALSE
    old_lik=Erlang_lik(x, r,weight,theta)
    target_var=var(x)
    target_mean = mean(x)
    
    update_theta=TRUE
    update_decision=0
    scale_itr=1
    for (itr in 1:101){
      y=0*r_range
      for (i in r_range){
        y[i]=Erlang_gradient(x,r_range[i],r,weight,theta)
      }
      if(max(y)<0.1 && update_decision==0){
        
        if(!update_theta){
          final_update = update_Erlang(x,r,weight,theta,update_theta = TRUE)
          weight= final_update$weight
          theta= final_update$theta
          
          if (min(weight)<1.1e-14){
            r=r[-which(weight<1.1e-14)]
            weight=weight[-which(weight<1.1e-14)]
            weight=weight/sum(weight)
            theta=mean(x)/sum(r*weight)
            update_result$lik=Erlang_lik(x,r,weight,theta)
          }
        }
        
        if(length(r)>round(1.2*comp) || theta <= 0.05){
          update_decision=1
        }else{
          if(scale_itr<8){
            theta=max(theta*(0.8^scale_itr), 0.05)
            scale_itr=scale_itr+1
            weight=update_Erlang(x,r,weight,theta,update_theta = FALSE)$weight
            
            y=0*r_range
            for (i in r_range){
              y[i]=Erlang_gradient(x,r_range[i],r,weight,theta)
            }
          }else{
            update_decision=1
          }
        }
      }
      
      ## if maximum gradient is less than 0.1, return values 
      if(update_decision==1 || itr==100){
        if (max(y)<.1 || itr==100){
          if (length(r)<=comp){
            cat('Caution: NPMLE is achieved with smaller number of components or exact number. Maximum value of gradeint function is ', max(y), '\n')
            small_comp=TRUE
            
            ini_final=find_best_shape_scale_comb(x, r, weight, theta, r_range=r_range)
            
            return(list(ini_r=ini_final$r,ini_weight=ini_final$weight, ini_theta=ini_final$theta,ini_lik=ini_final$lik, max_y=max(y),small_comp=small_comp))
          }
          
          r_set=weighted_kmeans(r,weight,comp)
          w_set=r_set/r_set/comp;theta_set = NULL
          likelihood_set=NULL;
          for (i in 1:length(r_set[,1])){
            update_result=update_Erlang(x,r_set[i,],w_set[i,], theta,update_theta=TRUE)
            w_set[i,]=update_result$weight
            theta_set[i]=update_result$theta
            likelihood_set[i]=update_result$lik
          }
          initial_r=initial_weight=initial_theta=initial_lik=NULL
          for (j in 1:dim(r_set)[1]){
            if (min(w_set[j,])>1.1e-14){
              initial_r=rbind(initial_r,r_set[j,])
              initial_weight=rbind(initial_weight,w_set[j,])
              initial_theta=c(initial_theta, theta_set[j])
              initial_lik=c(initial_lik,likelihood_set[j])
            }
          }
          if (is.null(initial_r)){
            r_set=t(all_comb);w_set=r_set/r_set/comp
            for (i in 1:length(r_set[,1])){
              update_result=update_Erlang(x,r_set[i,],w_set[i,],theta,update_theta=TRUE)
              w_set[i,]=update_result$weight
              theta_set[i]=update_result$theta
              likelihood_set[i]=update_result$lik
            }
            initial_r=initial_weight=initial_theta=initial_lik=NULL
            for (j in 1:dim(r_set)[1]){
              if (min(w_set[j,])>1.1e-14){
                initial_r=rbind(initial_r,r_set[j,])
                initial_weight=rbind(initial_weight,w_set[j,])
                initial_theta=c(initial_theta, theta_set[j])
                initial_lik=c(initial_lik,likelihood_set[j])
              }
            }
            max_lik=NULL
            if(is.null(initial_r)){
              max_lik = which.max(likelihood_set)
              initial_r = r_set[max_lik,]
              initial_weight=w_set[max_lik,]
              initial_theta=theta_set[max_lik]
            }
          }
          
          if(length(initial_lik)>1){
            max_indx=which.max(initial_lik)
            initial_r=initial_r[max_indx,]
            initial_weight=initial_weight[max_indx,]
            initial_theta=initial_theta[max_indx]
          }
          
          ini_final=find_best_shape_scale_comb(x, initial_r, initial_weight, initial_theta, r_range=r_range)
          
          
          return(list(ini_r=ini_final$r,ini_weight=ini_final$weight, ini_theta=ini_final$theta, ini_lik=ini_final$lik, max_y=max(y),small_comp=small_comp))
        }
      }
      
      ## find all local maximizers
      if (y[1]>y[2] & y[1]>0.1){
        new_r=1
      } else {
        new_r=NULL
      }
      dif_y=diff(y)
      for (i in 1:(length(dif_y)-1)){
        if (dif_y[i]>0 & dif_y[i+1]<0){
          if (y[i]>y[i+1] & y[i]>0.1) {
            new_r=c(new_r,i)
          }
          if (y[i]<y[i+1] & y[i+1]>0.1) {
            new_r=c(new_r,i+1)
          }
        }
      }
      ## update weight
      r=sort(c(r,new_r))
      r=unique(r)
      weight=rep(1,length(r))/length(r)
      update_result=update_Erlang(x,r,weight,theta,update_theta)
      weight=update_result$weight
      theta=update_result$theta
      
      if (min(weight)<1.1e-14){
        r=r[-which(weight<1.1e-14)]
        weight=weight[-which(weight<1.1e-14)]
        weight=weight/sum(weight)
        
        if(update_theta){
          theta=mean(x)/sum(r*weight)
          update_result$lik=Erlang_lik(x,r,weight,theta)
        }
      }
      
      # print(old_lik - update_result$lik)
      if(update_theta==TRUE){
        if(length(r)>max_num_shape){
          update_theta=FALSE
        }
      }
      
      old_lik=update_result$lik
      
    }
    
  }
  
  
  #EM
  Erlang_EM = function(x,support,weight,theta,itr_max=1000, epsil=1e-6, update_theta=TRUE){
    n=length(x)
    p=length(support)
    new_weight=weight
    new_theta=theta
    old_lik=Erlang_lik(x,support,weight,theta)
    for (itr in 1:itr_max){
      for(j in 1:p){
        new_weight[j]=new_weight[j]*(sum(Erlang_gradient(x,support[j],support,weight,theta))+n)/n
      }
      if(update_theta){
        new_theta=mean(x)/sum(support*new_weight)
      }
      
      new_lik=Erlang_lik(x,support,new_weight,new_theta)
      if (new_lik-old_lik>epsil){
        weight=new_weight
        theta=new_theta
        old_lik=new_lik
      }else{
        break
      }
    }
    return(list(r=support,weight=weight,theta=theta,lik=old_lik))
  }
  
  #weightd K-means
  weighted_kmeans <- function(r, weight, k, max_iter = 100) {
    all_combinations <- combn(r, k)  
    all_results <- list() 
    unique_centers <- matrix(nrow = 0, ncol = k)  
    
    for (i in 1:ncol(all_combinations)) {
      centers <- all_combinations[, i]  
      clusters <- rep(0, length(r))
      
      for (iter in 1:max_iter) {
        distances <- outer(r, centers, FUN = function(x, y) (x - y)^2)
        clusters <- apply(distances, 1, which.min)
        
        new_centers <- sapply(1:k, function(j) {
          if (sum(weight[clusters == j]) > 0) {
            sum(r[clusters == j] * weight[clusters == j]) / sum(weight[clusters == j])  
          } else {
            centers[j]
          }
        })
        
        if (all(abs(new_centers - centers) < 1e-6)){break}
        centers <- new_centers
      }
      
      sorted_centers <- sort(centers)
      all_results[[i]] <- round(sorted_centers)
    }
    
    unique_centers <- unique(do.call(rbind, all_results))
    
    return(unique_centers)
  }
  
  
  find_neighbor=function(x, can_r, weight, theta) {
    max_lik = -Inf 
    best_r = NULL
    best_weight = NULL
    best_theta = NULL
    
    for (j in 1:nrow(can_r)) {
      r_candidate = can_r[j, ]
      update_result = update_Erlang(x, r_candidate, weight, theta)
      lik_value = Erlang_lik(x, r_candidate, update_result$weight, update_result$theta)
      
      if (lik_value > max_lik) {
        max_lik = lik_value
        best_r = r_candidate
        best_weight = update_result$weight
        best_theta = update_result$theta
      }
    }
    
    return(list(r = best_r, weight = best_weight, theta = best_theta, lik = max_lik))
  }
  
  
  escape=function(x,support,weight,theta,r_range=1:100, max_itr=1,add=0, ...){
    for(itr in 1:max_itr){
      r=support
      y=0*r_range
      for (i in 1:length(r_range)){
        y[i]=Erlang_gradient(x,r_range[i],r,weight,theta)
      }
      
      ## Add global maximizer
      if(add==0){
        if (max(y)>1e-6 && !which.max(y) %in% r){
          new_r1=c(r,which.max(y))
        } else {
          lik=Erlang_lik(x,r,weight,theta)
          return(list(r=r,weight=weight,theta=theta,lik=lik,max_y=max(y)))
        }
      }else{
        if (max(y)>1e-6 && !which.max(y) %in% r){
          new_r1=c(r,which.max(y))
        } else {
          while(max(y)<=1e-6){
            theta=theta*0.8
            for (i in 1:length(r_range)){
              y[i]=Erlang_gradient(x,r_range[i],r,weight,theta)
            }
          }
          new_r1=c(r,which.max(y))
        }
      }

      new_r=new_r1
      can_r=t(combn(new_r,length(r)+add))
      can_r=t(apply(can_r,1,sort))
      
      #print(can_r)
      n_ind=dim(can_r)[1]
      weight=rep(1/(length(r)+add),length(r)+add)
      lik=Erlang_lik(x,r,weight,theta)
      
      for (j in 1:n_ind){
        update_result=update_Erlang(x,can_r[j,],weight,theta)
        if (update_result$lik>lik){
          r=update_result$support
          weight=update_result$weight
          theta=update_result$theta
          lik=update_result$lik
        }
      }
      ord=order(r)
      r=r[ord]
      weight=weight[ord]
    }
    return(list(r=r,weight=weight,theta=theta,lik=lik,max_y=max(y)))
  }
  
  
  gen_comb=function(vec,comp,extra=NULL){
    ind=permutations(n = length(vec), r = comp, v=vec, repeats.allowed=TRUE)
    if (!is.null(extra) && length(extra) > 0) {
      for (j in 1:length(extra)){
        ind2=ind
        ind2=ind2[-which(ind2[,extra[j]]==0),]
        ind2[,extra[j]]=-1
        ind=rbind(ind,ind2)
      }
    }
    
    return(ind)
  }  
  
  #Full function
  est_erlang=function(x,support=NULL,weight=NULL, theta=NULL,initial=TRUE,setseed=1,comp=NULL,r_range=1:100,...){
    
    if(initial==TRUE){
      if(is.null(comp)){
        cat('ERROR: Please fill in the "comp"! ', '\n')
      }
      if(length(support)==0){
        set.seed(setseed)
        support=sort(sample(1:50,comp,replace=FALSE))
        init <- Erlang_initial(x=x, breaks=8, shape=support)
        support <- init$shape
        weight <- init$weight
        theta <- init$theta
        
      }else if(length(support)>0){
        init <- Erlang_initial(x=x, breaks=8,shape=support)
        support <- init$shape
        weight <- init$weight
        theta <- init$thet
      }
    }
    
    ini_support=support
    ini_weight=weight
    ini_theta=theta
    comp=length(ini_support)
    max_shape=length(r_range)
    
    #initial using weighted K means
    ini=find_initial(x=x, r=ini_support, weight=ini_weight, theta=ini_theta, r_range=r_range,comp=comp)
    
    r=ini$ini_r;weight=ini$ini_weight
    theta=ini$ini_theta
    old_lik=ini$ini_lik
    
    for (itr in 1:200){
      add=0
      if(length(ini_support)!=length(r)){
        add=1
      }
      tmp=escape(x,r,weight,theta,add=add,r_range = r_range)
      r=tmp$r;weight=tmp$weight;theta=tmp$theta;old_lik=tmp$lik
      
      Er_gra=Er_gra2=0*r
      for (j in 1:length(r)){
        Er_gra[j]=Erlang_dgradient(x,r[j],r,weight,theta)
        Er_gra2[j]=Erlang_d2gradient(x,r[j],r,weight,theta)
      }
      direc=sign(Er_gra)
      extra=NULL
      if(sum(Er_gra2>0)>0){extra=which(Er_gra2>0 | abs(weight*Er_gra/sqrt(n))<1/n)}
      # ind=permutations(n = 2, r = comp, v=c(0,1), repeats.allowed=TRUE)
      ind=gen_comb(c(0,1),length(r),extra)
      n_ind=dim(ind)[1]
      can_r=0*ind
      can_weight=can_r
      lik_check=rep(0,n_ind)
      
      for (j in 1:n_ind){
        new_direc=direc*ind[j,]
        if (length(unique(r+new_direc))==length(r)){
          can_r[j,]=r+new_direc
        } else {
          can_r[j,]=r
        }
      }
      
      #delete shape = 0 or over max_shape
      rows0=which(apply(can_r, 1, function(row) 0 %in% row))
      rowsmax=which(apply(can_r, 1, function(row) (max_shape+1) %in% row))
      
      if(length(rows0)>0){
        can_r = can_r[-rows0,]
      }
      if(length(rowsmax)>0){
        can_r = can_r[-rowsmax,]
      }
      
      sorted_rows <- t(apply(can_r, 1, sort))
      if(dim(sorted_rows)[2]!=length(r)){sorted_rows=t(sorted_rows)}
      can_r <- unique(sorted_rows)
      
      findd=find_neighbor(x,can_r, weight, theta)
      new_r=findd$r
      new_weight=findd$weight
      new_theta = findd$theta
      new_lik = findd$lik
      
      if (min(new_weight)<1.1e-14){
        new_r=new_r[-which(new_weight<1.1e-14)]
        new_weight=new_weight[-which(new_weight<1.1e-14)]
        new_weight=new_weight/sum(new_weight)
        new_theta=mean(x)/sum(new_r*new_weight)
        findd$lik=Erlang_lik(x,new_r,new_weight,new_theta)
      }
      
      #if all shapes are not change, stop
      if(length(new_r)==comp){
        if (sum(new_r==r)==length(r) && abs(new_lik - old_lik)<1e-2){
          
          check_for_scale=find_best_shape_scale_comb(x, new_r, new_weight, new_theta,r_range=r_range)
          
          if(sum(check_for_scale$r==new_r)==length(r)){
            r=new_r;weight=new_weight;theta=new_theta
            old_lik = new_lik
            break
          }
          
          new_r= check_for_scale$r
          new_weight=check_for_scale$weight
          new_theta=check_for_scale$theta
          new_lik=check_for_scale$lik
        }
      }
      
      if(new_lik>=old_lik){
        r=new_r;weight=new_weight;theta=new_theta
        old_lik = new_lik
      }
      
    }
    return(list(setseed=setseed,ini_support=ini_support, ini_weight=ini_weight, ini_theta=ini_theta, r=r,weight=weight,theta=theta,lik=old_lik))
  }
  
  
  #second moment matching
  find_shape_scale_comb <- function(x, shape, weight, theta, r_range = 1:100, scale_grid=0.01) {
    num_components <- length(shape)
    target_means=shape*theta
    target_var=var(x)
    updated_weight=weight
    udpated_theta=theta
    scale_range=seq(max(0.05,theta*(1-0.5)), theta*(1+0.5), scale_grid)
    min_shape=r_range[1]
    max_shape=r_range[length(r_range)]
    
    all_valid_combinations <- lapply(scale_range, function(scale) {
      # Calculate the shape values by dividing the target means by the current scale value
      shape_vec <- round(target_means / scale)
      
      if (all(shape_vec >= min(r_range) & shape_vec <= max(r_range)) && length(unique(shape_vec))==num_components) {
        updated_result=update_Erlang(x,shape_vec ,weight = weight,theta=scale,update_theta)
        updated_weight=updated_result$weight
        updated_scale=updated_result$theta
        var <- Erlang_mean_var(shape_vec, updated_weight,updated_scale)$mixture_variance
        
        if( var/target_var<=1.01 & var/target_var>=0.99){
          return(list(r = shape_vec, weight=updated_weight, theta = updated_scale, lik=updated_result$lik))
        }else{
          return(NULL)
        }
      } else {
        return(NULL)
      }
    })
    # Remove NULL elements from the list
    all_valid_combinations <- Filter(Negate(is.null), all_valid_combinations)
    
    return(all_valid_combinations)
  }
  
  
  
  select_best_scale_shape <- function(x, valid_combinations) {
    max_likelihood <- -Inf
    best_combination <- NULL
    
    # Loop through each valid combination to calculate the likelihood
    for (comb in 1:length(valid_combinations)) {
      result<- update_Erlang(x, valid_combinations[[comb]]$r,  valid_combinations[[comb]]$weight, valid_combinations[[comb]]$theta,update_theta=TRUE)
      likelihood <- result$lik
      
      # Update the best combination if the current one has a higher likelihood
      if (likelihood > max_likelihood && min(weight)>1.1e-14) {
        max_likelihood <- likelihood
        best_combination <- result
      }
    }
    
    return(best_combination)
  }
  
  
  find_best_shape_scale_comb <- function(x, shape, weight, theta, r_range = 1:100,scale_prod=0.8, scale_grid=NULL, scale_length=50,var_epsil=0.05){
    num_components <- length(shape)
    target_means <- shape * theta
    target_var <- var(x)
    if(!is.null(scale_grid)){
      scale_length = (theta*(1+scale_prod)-max(0.05,theta*(1-scale_prod)))/scale_grid +1
    }
    scale_range=seq(max(0.05,theta*(1-scale_prod)), theta*(1+scale_prod), length.out=scale_length)
    min_shape <- r_range[1]
    max_shape <- r_range[length(r_range)]
    
    max_likelihood <- Erlang_lik(x,shape, weight, theta)
    best_combination <-  list(r = shape, weight = weight, theta = theta, lik = max_likelihood)
    
    
    # Iterate through scale range
    for (scale in scale_range) {
      # Calculate the shape values
      shape_vec <- round(target_means / scale)
      
      # Validate shape values
      if (!all(shape_vec >= min_shape & shape_vec <= max_shape) || length(unique(shape_vec)) != num_components) {
        next
      }
      
      # Update weight and scale
      updated_result <- update_Erlang(x, shape_vec, weight = weight, theta = scale, update_theta = TRUE)
      updated_weight <- updated_result$weight
      if(any(updated_weight < 1.1e-14)){
        next
      }
      
      updated_scale <- updated_result$theta
      likelihood <- updated_result$lik
      
      # Validate variance
      var <- Erlang_mean_var(shape_vec, updated_weight, updated_scale)$mixture_variance
      if (var / target_var > 1+var_epsil || var / target_var < 1-var_epsil) {
        next
      }
      
      # Update best combination if likelihood improves
      if (likelihood > max_likelihood) {
        max_likelihood <- likelihood
        best_combination <- list(r = shape_vec, weight = updated_weight, theta = updated_scale, lik = likelihood)
      }
    }
    
    return(best_combination)
  }
  
  
  #Calculate the mean and variance of Erlang mixtures
  Erlang_mean_var = function(shape, weight, theta){
    component_means <- shape * theta
    component_variances <- shape * theta^2
    
    mixture_mean <- sum(weight * component_means)
    
    mixture_variance <- sum(weight * (component_variances + (component_means - mixture_mean)^2))
    
    return(list(
      component_means = component_means,
      component_variances = component_variances,
      mixture_mean = mixture_mean,
      mixture_variance = mixture_variance
    ))
  }
  
  #bic function
  BIC <- function(loglike, M, n){
    BIC <- -2 * loglike + M * log(n)
    return(BIC)
  }

  