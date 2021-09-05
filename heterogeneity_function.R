## gt taken from https://www.medrxiv.org/content/10.1101/2020.07.23.20160317v4.full.pdf
## 10.56, 1.85
install.packages("SparseM")

library(zoo)
library(dplyr)
library(metR)
library(ggplot2)
library(magrittr)
library(grid)
library(gridExtra)
library(directlabels)
library(emdbook)
library(egg)
library(SparseM)
library(EpiEstim)
library(stringr)

gT_generator <- function(distribution_type, distribution_param_1, distribution_param_2=NULL){
  day_sequence <- seq(0,100,0.01)
  
  if(distribution_type == "gamma") gT_dist <- dgamma(day_sequence, shape=distribution_param_1, rate=distribution_param_2)/100
  if(distribution_type == "exponential") gT_dist <- dexp(day_sequence, rate=distribution_param_1)/sum(dexp(seq(1,100,1), rate=distribution_param_1))/100
  if(distribution_type == "step") gT_dist <- dunif(day_sequence, min=distribution_param_1-0.005, max=distribution_param_2+0.005)/100
  if(distribution_type == "delta") gT_dist <- dunif(day_sequence, min=distribution_param_1-0.005, max=distribution_param_1+0.005)/100
  if(distribution_type == "normal") gT_dist <- dnorm(day_sequence, mean=distribution_param_1, sd=distribution_param_2)/100
  
  #print(length(gT_dist))
  
  list <- list()
  list[["distribution_type"]] = distribution_type
  list[["day_sequence"]] = day_sequence
  list[["gT_dist"]] = gT_dist
  
  return(list)
}

gT_proposal <- function(mean_difference, distribution_type, distribution_param_1, distribution_param_2, plot=FALSE){
  if(distribution_type == "gamma"){
    gT_dist <- ((dgamma(seq(0,100,0.01),shape=distribution_param_1, rate=distribution_param_2) * (1/mean_difference)^(distribution_param_1-1) * exp(-(1/mean_difference-1)*distribution_param_2*seq(0,100,0.01)))/100)/mean_difference
    
    #print(sum(gT_dist))
  } 
  if(distribution_type == "step"){
    previous_mean <- (distribution_param_1+distribution_param_2)/2
    previous_variance <- distribution_param_2-distribution_param_1
    new_mean <- previous_mean*mean_difference
    
    new_min <- new_mean - previous_variance/2
    new_max <- new_mean + previous_variance/2
    
    print(c(new_min = new_min, new_max=new_max))
    
    gT_dist <- dunif(seq(0,100,0.01),min=new_min, max=new_max)/100
    #print(sum(gT_dist))
    
  }
  if(distribution_type == "delta"){
    previous_mean <- distribution_param_1
    new_mean <- distribution_param_1 * mean_difference
    
    print(c(new_mean = new_mean))
    
    gT_dist <- dunif(seq(0,100,0.01), min=new_mean-0.005, max=new_mean+0.005)/100
    #print(sum(gT_dist))
  }
  if(distribution_type == "normal"){
    previous_mean <- distribution_param_1
    new_mean <- distribution_param_1 * mean_difference
    
    print(c(new_mean = new_mean))
    
    gT_dist <- dnorm(seq(0,100,0.01), new_mean, distribution_param_2)/100
    #print(sum(gT_dist))
  }
  old_dist <- gT_generator(distribution_type, distribution_param_1, distribution_param_2)$gT_dist
  
  if(plot==TRUE){
    plot(gT_dist, type="l", xlim=c(0,1000), ylim=c(0,max(gT_dist,old_dist)), col="red")
    lines(old_dist) 
  }
  
  return(gT_dist)  
}

gT_proposal(mean_difference=0.5, distribution_type="gamma", distribution_param_1=10.56, distribution_param_2=1.85, plot=TRUE)

gT_proposal_trunc <- function(truncation_point, distribution_type, distribution_param_1, distribution_param_2, plot=FALSE, normalise=TRUE){
  old <- gT_generator(distribution_type, distribution_param_1, distribution_param_2)
  old_dist <- old$gT_dist
  old_cumulative_dist <- cumsum(old_dist)
  
  new_dist <- old_dist
  new_dist[c(which(old_cumulative_dist>=truncation_point)[1]:length(new_dist))] <- 0
  new_dist <- new_dist
  
  if(normalise==TRUE) new_dist <- new_dist/sum(new_dist)
  
  if(plot == TRUE){
    plot(old_dist, type="l", xlim=c(0,1000), ylim=c(0,max(new_dist,old_dist)))
    lines(new_dist, col="red") 
  }
  return(new_dist)
} 

gT_weighted_sum <- function(mean_difference, risk_trans, proposal, distribution_type, distribution_param_1, distribution_param_2, plot=FALSE){
  gT_1 <- gT_generator(distribution_type = "gamma", distribution_param_1, distribution_param_2)$gT_dist
  if(proposal == "relative") gT_2 <- gT_proposal(mean_difference = mean_difference, distribution_type = "gamma", distribution_param_1, distribution_param_2, plot)
  else if(proposal == "trunc") gT_2 <- gT_proposal_trunc(truncation_point = mean_difference, distribution_type = "gamma", distribution_param_1, distribution_param_2, plot, normalise = TRUE)
  
  gT_overall = cbind(gT_1, gT_2)
  gT_weighted = rowSums(t(t(gT_overall)*risk_trans))
  
  gT_normalised = gT_weighted/sum(gT_weighted)
  
  if(plot==TRUE) plot(gT_normalised, type="l")
  
  return(gT_normalised)
  
}

gT_weighted_sum(0.2, risk_trans=c(1, 0.5), proposal = "trunc", "gamma", 10.56, 1, plot=TRUE)

infectiousness_trunc_plotter <- function(R, truncation_vec, distribution_type, distribution_param_1, distribution_param_2, lower_letter, upper_letter){
  gT_init <- gT_generator(distribution_type, distribution_param_1, distribution_param_2)$gT_dist
  gT_series <- matrix(NA, ncol=length(truncation_vec)+1, nrow=length(gT_init))
  gT_series[,1] <- gT_init
  
  for(j in 1:(ncol(gT_series)-1)) gT_series[,j+1] = gT_proposal_trunc(truncation_vec[j], distribution_type = "gamma", distribution_param_1 = distribution_param_1, distribution_param_2 = distribution_param_2, normalise=FALSE)
  
  colnames(gT_series) <- c(1,truncation_vec)
  gT_series_plot <- as.data.frame(gT_series) %>%
    mutate(day = seq(0,100,0.01)) %>%
    tidyr::gather(truncation_point, "value", 1:(ncol(gT_series)))
  
  dat_text <- data.frame(
    label = LETTERS[lower_letter:upper_letter],
    truncation_point  = as.factor(c(truncation_vec,1)),
    x     = c(rep(min(gT_series_plot$day),upper_letter-lower_letter+1)),
    y     = c(rep(max(R*gT_series_plot$value)*100,upper_letter-lower_letter+1))
  )
  
  p <- ggplot(data=gT_series_plot, aes(x=as.numeric(day), y=R*as.numeric(value)*100, color=truncation_point)) +
    scale_color_grey(start = 0.5, end = 0.1) +
    geom_line() +
    xlim(c(0,25)) +
    theme_bw() +
    facet_wrap(~truncation_point, nrow=1) +
    labs(x = expression(tau), y = expression(paste(beta,"(",t,",",tau,")"))) +
    geom_text(data = dat_text, mapping = aes(x = x, y = y, label = label, fontface="bold"), size=5, color="black") + 
    theme(strip.text.x = element_blank(),
          axis.text.y = element_blank())
  
  ggsave(plot=p, filename=paste("graphs/trunc_inf_dist.png"), width=10, height=2.5)
  
  return(p)
}

trunc_plots <- infectiousness_trunc_plotter(R=3,truncation_vec=c(0.25,0.5,0.75), distribution_type="gamma", distribution_param_1=10.56, distribution_param_2=1.85, lower_letter=1, upper_letter=4)
trunc_plots

symp_asymp_gt_plotter <- function(mean_diff_vec, distribution_type, distribution_param_1, distribution_param_2, lower_letter, upper_letter){
  gT_init <- gT_generator(distribution_type, distribution_param_1, distribution_param_2)$gT_dist
  gT_series1  <- matrix(NA, ncol=length(mean_diff_vec), nrow=length(gT_init))
  gT_series2 <- matrix(NA, ncol=length(mean_diff_vec), nrow=length(gT_init))
  
  for(j in 1:(ncol(gT_series1))) gT_series1[,j] = gT_proposal(mean_diff_vec[j], distribution_type = "gamma", distribution_param_1, distribution_param_2) 
  for(k in 1:(ncol(gT_series2))) gT_series2[,k] = gT_proposal(1, distribution_type = "gamma", distribution_param_1, distribution_param_2)
  
  
  colnames(gT_series1) <- c(mean_diff_vec)
  gT_series_plot1 <- as.data.frame(gT_series1) %>%
    mutate(day = seq(0,100,0.01),
           group="asymptomatic") %>%
    tidyr::gather(mean_diff, "value", 1:(ncol(gT_series1))) 
  
  colnames(gT_series2) <- c(mean_diff_vec)
  gT_series_plot2 <- as.data.frame(gT_series2) %>%
    mutate(day = seq(0,100,0.01), 
           group = "symptomatic") %>%
    tidyr::gather(mean_diff,"value", 1:(ncol(gT_series2)))
  
  gT_series_plot <- rbind(gT_series_plot1,gT_series_plot2) 
  
  dat_text <- data.frame(
    label = LETTERS[lower_letter:upper_letter],
    mean_diff   = mean_diff_vec,
    x     = c(rep(min(gT_series_plot$day),upper_letter-lower_letter+1)),
    y     = c(rep(max(gT_series_plot$value)*100,upper_letter-lower_letter+1))
  )
  
  p <- ggplot(data=gT_series_plot, aes(x=as.numeric(day), y=as.numeric(value)*100)) +
    geom_line(aes(color = group, linetype=group)) +
    scale_color_manual(values = c("dark grey", "black")) +
    scale_linetype_manual(values = c("longdash","solid")) +
    xlim(c(0,25)) +
    theme_bw() +
    facet_wrap(~mean_diff, labeller = label_both, nrow=1) +
    labs(x = expression(tau), y = expression(paste(omega,"(",tau,")"))) + 
    theme(strip.text.x = element_blank(),
          axis.text.y = element_blank()) +
    geom_text(data = dat_text, mapping = aes(x = x, y = y, label = label, fontface="bold"), size=5)
  
  ggsave(plot=p, filename=paste("graphs/symp_gT_dist.png"), width=17, height=5)
  
  return(p)
}

symp_asymp_gt <- symp_asymp_gt_plotter(mean_diff_vec=c(0.5,2), distribution_type="gamma", distribution_param_1=10.56, distribution_param_2=1.85, lower_letter=1, upper_letter=2)
symp_asymp_gt

vaccine_impact_plotter <- function(R, vaccine_impact_vec, distribution_type, distribution_param_1, distribution_param_2, lower_letter, upper_letter){
  gT_init <- gT_generator(distribution_type, distribution_param_1, distribution_param_2)$gT_dist
  gT_series1  <- matrix(NA, ncol=length(vaccine_impact_vec), nrow=length(gT_init))
  gT_series2 <- matrix(NA, ncol=length(vaccine_impact_vec), nrow=length(gT_init))
  
  for(j in 1:(ncol(gT_series1))) gT_series1[,j] = gT_proposal(vaccine_impact_vec[j], distribution_type = "gamma", distribution_param_1 =10.56, distribution_param_2 = 1) * vaccine_impact_vec[j]^2
  for(k in 1:(ncol(gT_series2))) gT_series2[,k] = gT_proposal(1, distribution_type = "gamma", distribution_param_1 =10.56, distribution_param_2 = 1)
  
  
  colnames(gT_series1) <- c(vaccine_impact_vec)
  gT_series_plot1 <- as.data.frame(gT_series1) %>%
    mutate(day = seq(0,100,0.01)) %>%
    tidyr::gather(vaccine_impact, "value", 1:(ncol(gT_series1))) %>%
    mutate(reduction_factor = as.numeric(vaccine_impact))
  
  colnames(gT_series2) <- c(vaccine_impact_vec)
  gT_series_plot2 <- as.data.frame(gT_series2) %>%
    mutate(day = seq(0,100,0.01), 
           vaccine_impact = 1) %>%
    tidyr::gather(reduction_factor, "value", 1:(ncol(gT_series2)))
  
  
  gT_series_plot <- rbind(gT_series_plot1,gT_series_plot2) %>% mutate(group = ifelse(vaccine_impact == 1, "unvaccinated", "vaccinated"))
  
  dat_text <- data.frame(
    label = LETTERS[lower_letter:upper_letter],
    reduction_factor  = vaccine_impact_vec,
    x     = c(rep(min(gT_series_plot$day),upper_letter-lower_letter+1)),
    y     = c(rep(max(R*gT_series_plot$value)*100,upper_letter-lower_letter+1))
  )
  
  p <- ggplot(data=gT_series_plot, aes(x=as.numeric(day), y=R*as.numeric(value)*100)) +
    geom_line(aes(linetype=group, color=group)) +
    scale_color_manual(values = c("black", "dark grey")) +
    scale_linetype_manual(values = c("solid", "longdash")) +
    xlim(c(0,20)) +
    theme_bw() +
    facet_wrap(~reduction_factor, labeller = label_both, nrow=1) +
    labs(x = expression(tau), y = expression(paste(beta,"(",t,",",tau,")"))) +
    theme(strip.text.x = element_blank(),
          axis.text.y = element_blank()) +
    geom_text(data = dat_text, mapping = aes(x = x, y = y, label = label, fontface="bold"), size=5)
  
  
  ggsave(plot=p, filename=paste("graphs/vaccine_inf_dist.png"), width=17, height=5)
  
  return(p)
}

vaccine_impact <- vaccine_impact_plotter(R=3, vaccine_impact_vec=seq(0.25,0.75,0.25), distribution_type="gamma", distribution_param_1=10.56, distribution_param_2=1.85, lower_letter = 1, upper_letter=3)
vaccine_impact


mult_format <- function() {
  function(x) format(x/(1-x),digits = 2) 
}

laplace_transform <- function(w, generation_int, day_granularity) generation_int %*% exp(-w*seq(0,((length(generation_int)-1))*day_granularity, by=day_granularity)) 

#laplace_transform(w=0.1, generation_int = gt_gamma1, day_granularity = 0.01)

laplace_transform_matrix <- function(w, gT_list, day_granularity){
  laplace_matrix <- matrix(NA, ncol=ncol(gT_list), nrow=ncol(gT_list))
  for(j in 1:ncol(gT_list)){
    laplace_matrix[,j] <- laplace_transform(w, gT_list[,j], day_granularity)
  }
  return(laplace_matrix)
}

contact_matrix_generator <- function(population_2_size, associativity_vec, i, log=F){
  if(log==T){
    if(population_2_size < 0.5)       contact_matrix=matrix(c(1-population_2_size *(1-population_2_size)^(associativity_vec[i]-1),   (1-population_2_size)^associativity_vec[i], 
                                                              population_2_size *(1-population_2_size)^(associativity_vec[i]-1), 1-(1-population_2_size)^associativity_vec[i]), nrow=2, byrow=T)
    
    else if(population_2_size >= 0.5) contact_matrix=matrix(c(1-population_2_size^associativity_vec[i],   population_2_size^(associativity_vec[i]-1)*(1-population_2_size), 
                                                              population_2_size^associativity_vec[i], 1-population_2_size^(associativity_vec[i]-1)*(1-population_2_size)), nrow=2, byrow=T)
  }
  
  if(log==F){
    if(associativity_vec[i] >= 0){
      if(population_2_size <= 0.5)        contact_matrix=matrix(c(1-population_2_size + associativity_vec[i]*population_2_size/(length(associativity_vec[associativity_vec >= 0])-1), 1-population_2_size - associativity_vec[i]*(1-population_2_size)/(length(associativity_vec[associativity_vec >= 0])-1),
                                                                  population_2_size - associativity_vec[i]*population_2_size/(length(associativity_vec[associativity_vec >= 0])-1),   population_2_size + associativity_vec[i]*(1-population_2_size)/(length(associativity_vec[associativity_vec >= 0])-1)), nrow=2, byrow=T)
      
      else if(population_2_size > 0.5){
        population_2_size = 1-population_2_size
        
        contact_matrix=matrix(c(  population_2_size + associativity_vec[i]*(1-population_2_size)/(length(associativity_vec[associativity_vec >= 0])-1), population_2_size - associativity_vec[i]*population_2_size/(length(associativity_vec[associativity_vec >= 0])-1),
                                  1-population_2_size - associativity_vec[i]*(1-population_2_size)/(length(associativity_vec[associativity_vec >= 0])-1), 1-population_2_size + associativity_vec[i]*population_2_size/(length(associativity_vec[associativity_vec >= 0])-1)), nrow=2, byrow=T)
      }
    }     
    else if(associativity_vec[i] < 0){
      if(population_2_size <= 0.5)         contact_matrix=matrix(c(1-population_2_size + associativity_vec[i]*((population_2_size/(1-population_2_size)-population_2_size)/(length(associativity_vec[associativity_vec <= 0])-1)), (1-population_2_size) - associativity_vec[i]*population_2_size/(length(associativity_vec[associativity_vec <= 0])-1),
                                                                   population_2_size - associativity_vec[i]*((population_2_size/(1-population_2_size)-population_2_size)/(length(associativity_vec[associativity_vec <= 0])-1)),    population_2_size  + associativity_vec[i]*population_2_size/(length(associativity_vec[associativity_vec <= 0])-1) ), nrow=2, byrow=T)
      
      else if(population_2_size > 0.5){
        population_2_size = 1 - population_2_size
        
        contact_matrix=matrix(c(population_2_size  + associativity_vec[i]*population_2_size/(length(associativity_vec[associativity_vec <= 0])-1), population_2_size - associativity_vec[i]*((population_2_size/(1-population_2_size)-population_2_size)/(length(associativity_vec[associativity_vec <= 0])-1)),    
                                (1-population_2_size) - associativity_vec[i]*population_2_size/(length(associativity_vec[associativity_vec <= 0])-1), 1-population_2_size + associativity_vec[i]*((population_2_size/(1-population_2_size)-population_2_size)/(length(associativity_vec[associativity_vec <= 0])-1))), nrow=2, byrow=T)
      }
    }
  }
  
  return(contact_matrix)
}

contact_matrix_generator(population_2_size=0.2, associativity_vec=c(-5,-4,-3,-2,-1,0,1,2,3,4,5), i=11, log=F)
contact_matrix_generator(population_2_size=0.8, associativity_vec=c(-5,-4,-3,-2,-1,0,1,2,3,4,5), i=11, log=F)

contact_matrix_generator(population_2_size=0.5, associativity_vec=c(-5,-4,-3,-2,-1,0,1,2,3,4,5), i=2, log=F)
contact_matrix_generator(population_2_size=0.5, associativity_vec=c(-5,-4,-3,-2,-1,0,1,2,3,4,5), i=10, log=F)

contact_matrix_generator(population_2_size=0.2, associativity_vec=c(-5,-4,-3,-2,-1,0,1,2,3,4,5), i=1, log=F)
contact_matrix_generator(population_2_size=0.8, associativity_vec=c(-5,-4,-3,-2,-1,0,1,2,3,4,5), i=1, log=F)


R_solver_nd <- function(w, gT_list, risk_matrix, day_granularity){
  laplace_matrix <- matrix(NA, nrow=ncol(gT_list), ncol(gT_list))
  for(j in 1:ncol(laplace_matrix)) laplace_matrix[,j] <- laplace_transform(w, gT_list[,j], day_granularity)
  
  eigenvalue <- Re(eigen(risk_matrix*laplace_matrix)$values[1])
  overall_R <- 1/eigenvalue
  
  next_generation_matrix <- overall_R * risk_matrix
  
  R_pop_1 <- colSums(next_generation_matrix)[1]
  R_pop_2 <- colSums(next_generation_matrix)[2]
  
  eigen_vectors <- Re(eigen(next_generation_matrix * laplace_matrix)$vectors)
  
  dominant_eigen <- Re(eigen_vectors[,1])
  #print(dominant_eigen[1])
  
  output <- c(overall_R = overall_R, R1 = R_pop_1, R2 = R_pop_2, k1 = abs(dominant_eigen[1]), k2 = abs(dominant_eigen[2]))
  
  return(output)
}

R_solver_nd(w = 0.2, gT_list = cbind(gt_no_isol, gt_isol, gt_asymp), risk_matrix = risk_matrix_UK, day_granularity = 0.01)

#R_solver_nd(w=0.15, gT_list = cbind(gt1, gt2), risk_matrix = relative_risk_mat_generator(risk_trans=c(1,0.25), contact_matrix = contact_matrix_homo),  day_granularity = 0.01)


w_solver_nd <- function(R, gT_list, risk_matrix, day_granularity){
  k <- 0.1
  w_trials <- seq(-1,1,k)
  
  output_trials <- vector(length=length(w_trials))
  laplace_matrix <- matrix(NA, nrow=ncol(gT_list), ncol=ncol(gT_list))
  
  while(k > 0.0000001){
    for(i in 1:length(output_trials)){
      for(j in 1:ncol(laplace_matrix)) laplace_matrix[,j] <- laplace_transform(w_trials[i], gT_list[,j], day_granularity)
      overall_matrix <- R*risk_matrix*laplace_matrix-diag(ncol(gT_list))
      output_trials[i] <- det(overall_matrix)
      #print(output_trials)
    }
    
    #print(w_trials)
    #print(output_trials)
    if (all(output_trials >0 ) & k > 1e-5) w_trials = seq(w_trials[1], w_trials[length(w_trials)], k/10)
    else if (all(output_trials >0 ) & k <= 1e-5){ 
      w_trials = seq(w_trials[1], w_trials[length(w_trials)], k/10)
      break
    }
    else if(output_trials[1] < 0)  w_trials = seq(w_trials[which(output_trials == max(output_trials[output_trials<0]))], w_trials[which(output_trials == max(output_trials[output_trials<0]))]+k, k/10)
    else if(output_trials[1] > 0) w_trials = seq(w_trials[max(which(output_trials<0))], w_trials[max(which(output_trials<0))]+k, k/10)
    
    output_trials <- vector(length=length(w_trials))
    k <- k/10
  }
  best_w <- w_trials[1]
  return(best_w)
}

relative_risk_mat_generator <- function(risk_trans, contact_matrix, susceptibility_vec=NA){
  relative_risk_matrix <- t(risk_trans * t(contact_matrix)) * ifelse(is.na(susceptibility_vec)==TRUE, 1, susceptibility_vec)
  #print(relative_risk_matrix)
  relative_risk_matrix_normalised <- relative_risk_matrix / Re(eigen(relative_risk_matrix)$values[1])
  #print(relative_risk_matrix_normalised)
  return(relative_risk_matrix_normalised)
}

relative_risk_mat_generator(risk_trans = c(1,0.5), contact_matrix = matrix(c(2,3,4,5),nrow=2,byrow=T), susceptibility_vec = c(1,1))

numerical_outbreak_simulation <- function(R, gT_list, risk_matrix){
  infecteds <- matrix(NA, ncol=ncol(gT_list), nrow=30000)
  infecteds[(1:nrow(gT_list)),] <- 1
  #risk_matrix <- relative_risk_mat_generator(risk_trans, contact_matrix)
  for(i in (nrow(gT_list)+1):nrow(infecteds)){  
    for(j in 1:ncol(infecteds)){
      colSums(infecteds[(i-nrow(gT_list)):(i-1),] * gT_list[nrow(gT_list):1,])
      infecteds[i,j] <- R*colSums(infecteds[(i-nrow(gT_list)):(i-1),] * gT_list[nrow(gT_list):1,]) %*% contact_matrix[j,]
    } 
  }
  total_infecteds <- rowSums(infecteds)
  plot(log(total_infecteds), type="l")
  #print(infecteds)
  w <- ((log(total_infecteds[length(total_infecteds)-1])-log(total_infecteds[2*nrow(gT_list)]))/(length(total_infecteds)-1-2*nrow(gT_list)))
  print(w)
  
  list <- list()
  list[["infecteds"]] <- infecteds
  list[["growth rate"]] <- w
  return(list)
}


R_by_gT_A_w <- function(distribution_type="gamma", distribution_param_1=10.56, distribution_param_2=1.85, NPI="vaccine",
                        w_vec=seq(-0.35,0.35,0.01), associativity_vec = seq(0), population_sizes_vec = seq(0.2, 0.8, 0.2),
                        gT_diff = "relative", mean_diff_vec = seq(0.5, 2, 0.02), relative_inf_vec = seq(0.01, 0.99, 0.02),
                        day_granularity=0.01, log_associativity=F){
  
  #browser()
  gT_group1 <- gT_generator(distribution_type="gamma", distribution_param_1 =10.56, distribution_param_2 = 1)$gT_dist
  
  output_storage <- array(NA, dim=c(length(population_sizes_vec), length(associativity_vec), length(mean_diff_vec), length(relative_inf_vec), length(w_vec), 5), dimnames=list(population_sizes_vec,associativity_vec,mean_diff_vec, relative_inf_vec,w_vec, c("overall_R", "R1", "R2", "k1", "k2")))
  
  for(i in 1:length(population_sizes_vec)){
    print(c("population_size = ", population_sizes_vec[i]), quote=F)
    for(j in 1:length(associativity_vec)){
      print(c("associativity = ", associativity_vec[j]), quote=F)
      for(k in 1:length(mean_diff_vec)){
        print(c("mean_diff_vec =", mean_diff_vec[k]), quote=FALSE)
        
        contact_matrix = contact_matrix_generator(population_2_size=population_sizes_vec[i], associativity_vec, j, log=log_associativity)
        
        for(l in 1:length(relative_inf_vec)){
          print(c("relative_inf_vec =", relative_inf_vec[l]), quote=F)
          if(gT_diff == "relative"){
            risk_trans = c(1-relative_inf_vec[l],relative_inf_vec[l])
            if(NPI=="vaccine") risk_trans = c(1,relative_inf_vec[l])
            risk_matrix <- relative_risk_mat_generator(risk_trans = risk_trans, contact_matrix)
            
            gT_group2 <- gT_proposal(mean_difference = mean_diff_vec[k], distribution_type, distribution_param_1, distribution_param_2)
            
          }
          else if(gT_diff == "trunc"){
            risk_matrix <- relative_risk_mat_generator(risk_trans = c(1, relative_inf_vec[l]),
                                                       contact_matrix=matrix(c(1-population_sizes_vec[i]*(1-population_sizes_vec[i])^(associativity_vec[j]-1),   (1-population_sizes_vec[i])^associativity_vec[j], 
                                                                               population_sizes_vec[i]*(1-population_sizes_vec[i])^(associativity_vec[j]-1), 1-(1-population_sizes_vec[i])^associativity_vec[j]), nrow=2, byrow=T))
            
            gT_group2 <- gT_proposal_trunc(truncation_point = relative_inf_vec[l], distribution_type, distribution_param_1, distribution_param_2, normalise = TRUE)
            
          }
          
          gT_set <- cbind(gT_group1, gT_group2)
          for(m in 1:length(w_vec)){
            #print(w_vec[m])
            output_storage[i,j,k,l,m,] <- R_solver_nd(w_vec[m], gT_set, risk_matrix, day_granularity)
            
          }
        }
      }
    }
  }
  print(output_storage)
  
  output <- as.data.frame.table(output_storage) %>% set_colnames(c("population_size", "associativity", "mean_diff","relative_inf", "w", "measure", "value")) %>% as.data.frame() %>%
    mutate(population_size = as.numeric(as.character(population_size)), 
           associativity = as.numeric(as.character(associativity)),
           mean_diff = as.numeric(as.character(mean_diff)),
           relative_inf = as.numeric(as.character(relative_inf)),
           w = as.numeric(as.character(w)),
           value = as.numeric(as.character(value))) %>%
    tidyr::spread(measure, value) %>% rowwise() %>%
    mutate(R_simple = R_solver_nd(w, gT_list = cbind(gT_group1), risk_matrix = matrix(c(1)), day_granularity)[1], 
           R_relative = overall_R/R_simple) %>%
    mutate(growth_rate = as.factor(as.character(w)), population_size = as.factor(population_size))
  
  return(output)
}


R_by_gT_A_w_simple <- function(distribution_type="gamma", distribution_param_1=10.56, distribution_param_2=1.85, 
                               w_vec=seq(-0.35,0.35,0.01), associativity_vec = seq(0), population_sizes_vec = seq(0.2, 0.8, 0.2),
                               gT_diff = "relative", mean_diff_vec = seq(0.5, 2, 0.02), relative_inf_vec = seq(0.01, 0.99, 0.02),
                               day_granularity=0.01, log_associativity=F, susceptibility_vec = NULL){
  
  #browser()
  gT_group1 <- gT_generator(distribution_type="gamma", distribution_param_1 = distribution_param_1, distribution_param_2 = distribution_param_2)$gT_dist
  
  output_storage <- array(NA, dim=c(length(population_sizes_vec), length(associativity_vec), length(mean_diff_vec), length(relative_inf_vec), length(w_vec), 5), dimnames=list(population_sizes_vec,associativity_vec,mean_diff_vec, relative_inf_vec,w_vec, c("overall_R", "R1", "R2", "k1", "k2")))
  
  output <- as.data.frame.table(output_storage) %>% set_colnames(c("population_size", "associativity", "mean_diff","relative_inf", "w", "measure", "value")) %>% as.data.frame() %>%
    mutate(population_size = as.numeric(as.character(population_size)), 
           associativity = as.numeric(as.character(associativity)),
           mean_diff = as.numeric(as.character(mean_diff)),
           relative_inf = as.numeric(as.character(relative_inf)),
           w = as.numeric(as.character(w)),
           value = as.numeric(as.character(value)), 
           gT_diff = as.character(gT_diff),
           distribution_type = as.character(distribution_type),
           param_1 = as.numeric(distribution_param_1), 
           param_2 = as.numeric(distribution_param_2)) %>%
    tidyr::spread(measure, value) 
  
  if(length(associativity_vec) !=1) output <- output %>% filter(associativity != min(associativity), associativity != max(associativity))
  for(i in 1:nrow(output)){
    
    contact_matrix = contact_matrix_generator(population_2_size=output[i, "population_size"], associativity_vec, i=which(output[i,"associativity"]==associativity_vec), log=log_associativity)
    
    if(gT_diff == "relative"){
      risk_trans = c(1-output[i,"relative_inf"],output[i,"relative_inf"])
      gT_group2 = gT_proposal(mean_difference = output[i,"mean_diff"], distribution_type = output[i,"distribution_type"], distribution_param_1 = distribution_param_1, distribution_param_2 = distribution_param_2)
      risk_matrix <- relative_risk_mat_generator(risk_trans = risk_trans, contact_matrix)
    } 
    
    if(gT_diff == "vaccine"){
      risk_trans = c(1, output[i,"relative_inf"])
      gT_group2 = gT_proposal(mean_difference = output[i,"mean_diff"], distribution_type = output[i,"distribution_type"], distribution_param_1 = distribution_param_1, distribution_param_2 = distribution_param_2)
      risk_matrix <- relative_risk_mat_generator(risk_trans = risk_trans, contact_matrix, susceptibility_vec = susceptibility_vec)
    }
    
    if(gT_diff == "trunc"){
      risk_trans = c(1, output[i,"relative_inf"])
      gT_group2 = gT_proposal_trunc(truncation_point = output[i,"relative_inf"], distribution_type = output[i,"distribution_type"], distribution_param_1 = distribution_param_1, distribution_param_2 = distribution_param_2)
      risk_matrix <- relative_risk_mat_generator(risk_trans = risk_trans, contact_matrix)
    }
    
    
    
    #gT_weighted = rowSums(t(t(cbind(gT_group1, gT_group2))*c(1-output[i, "population_size"],output[i,"population_size"])*c(risk_trans)))
    #gT_weighted_normalised = gT_weighted / sum(gT_weighted)
    
    output[i,c(10:14)] <- R_solver_nd(output[i,"w"], gT_list=cbind(gT_group1, gT_group2), risk_matrix=risk_matrix, day_granularity)
    
    gT_weighted = rowSums(t(t(cbind(gT_group1, gT_group2))*c(output[i, "k1"],output[i,"k2"])*c(risk_trans)))
    gT_weighted_normalised = gT_weighted / sum(gT_weighted)
    
    output[i,"overall_R_weighted"] <-  R_solver_nd(output[i,"w"], gT_list=cbind(gT_weighted_normalised), risk_matrix=matrix(c(1)), day_granularity)["overall_R"]
  }
  
  output_final <- output %>% rowwise() %>%
    mutate(R_simple = R_solver_nd(w, gT_list = cbind(gT_group1), risk_matrix = matrix(c(1)), day_granularity)[1], 
           R_relative = overall_R / R_simple, 
           R_relative_weighed = overall_R/overall_R_weighted,
           k1_normalised = k1/(k1+k2),
           k2_normalised = k2/(k1+k2)) %>%
    mutate(growth_rate = as.factor(as.character(w)), population_size = as.factor(population_size))
  
  return(output_final)
}

## for trunc, relative inf should be 1
## for non-trunc, relative inf should be 0.5 ??

ptm <- proc.time()
R_by_w_gT_trunc_simple <- R_by_gT_A_w_simple(w_vec=round(seq(-0.3,0.3,0.15),3), distribution_type="gamma", distribution_param_1=10.56, distribution_param_2=1.85, 
                                             associativity_vec = c(1), population_sizes_vec = c(0.2,0.5,0.8),
                                             gT_diff = "trunc", mean_diff_vec = seq(1), relative_inf_vec = round(seq(0.01,0.99,0.02),3),
                                             day_granularity=0.01, log_associativity = T)

R_by_gT_rel_simple <- R_by_gT_A_w_simple(w_vec=round(seq(-0.3,0.3,0.15),3), distribution_type="gamma", distribution_param_1=10.56, distribution_param_2=1.85,
                                         associativity_vec = c(1), population_sizes_vec = c(0.2,0.5,0.8), gT_diff = "relative", relative_inf_vec = round(seq(0.01, 0.99, 0.02),3), 
                                         mean_diff_vec = c(0.5,1,2), day_granularity=0.01, log_associativity = T)

R_by_gT_rel_A_simple <- R_by_gT_A_w_simple(w_vec=round(seq(-0.3,0.3,0.15),3), distribution_type="gamma", distribution_param_1=10.56, distribution_param_2=1.85,
                                           associativity_vec = c(-20:20), population_sizes_vec = c(0.2,0.5,0.8),
                                           gT_diff = "vaccine", mean_diff_vec = seq(0.25, 1.0, 0.25), relative_inf_vec = (seq(0.25, 1.0, 0.25))^2,
                                           day_granularity=0.01, log_associativity = F, susceptibility_vec = c(1,0.3)) %>%
  mutate(scenario = ifelse(mean_diff^2 == relative_inf, mean_diff, NA)) %>%
  filter(is.na(scenario)==FALSE)




R_isol_no_isol <- ggplot(data = R_by_w_gT_trunc_simple %>% rename(isolating_population = population_size) %>% arrange(w) %>% mutate(growth_rate=factor(growth_rate, levels=c(-0.3, -0.15, 0, 0.15,0.3))), aes(x=relative_inf, y=R_relative)) +
  geom_line(aes(color=growth_rate, linetype=isolating_population)) +
  scale_linetype_manual(values=c("dotted", "dashed","solid")) +
  theme_bw() +
  labs(y = expression(R["adjusted"] / R["unadjusted"]), x = "Proportion of infectiousness passed at point of isolation") +
  annotate("text", x=0, y=max(R_by_w_gT_trunc_simple$R_relative), label="E", fontface="bold", size=5)

isol_no_isol_plot <- ggarrange(trunc_plots, R_isol_no_isol)

ggsave(plot=isol_no_isol_plot, filename=paste("graphs/error_R_isol_no_isol_11.png"), width=10, height=5)


R_rel <- ggplot(R_by_gT_rel_simple %>% filter(mean_diff != 1) %>% rename(asymptomatic_population = population_size) %>% arrange(w) %>% mutate(growth_rate = factor(growth_rate, levels=c(-0.3, -0.15, 0 , 0.15, 0.3))), aes(x=relative_inf, y=R_relative)) +
  geom_line(aes(color=growth_rate, linetype=asymptomatic_population)) +
  scale_linetype_manual(values=c("dotted", "dashed","solid")) +
  facet_wrap(~mean_diff, labeller = label_both) +
  scale_x_continuous(labels=mult_format()) +
  theme_bw() +
  theme(strip.text.x = element_blank()) +
  labs(y = expression(R["adjusted"] / R["unadjusted"]), x = "Relative infectiousness of asymptomatics") +
  geom_text(data = data.frame(x = c(0,0), y=rep(max(R_by_gT_rel_simple$R_relative),2), label=c("C","D"), mean_diff = c(0.5,2)), aes(x=x, y=y, label=label), fontface="bold", size=5)

symp_asymp_plot <- ggarrange(symp_asymp_gt, R_rel)

ggsave(plot=symp_asymp_plot, filename=paste("graphs/error_R_symp_asymp_11.png"), width=10, height=5)


R_vacc_no_vacc_A <- ggplot(R_by_gT_rel_A_simple %>% filter(population_size %in% c(0.2,0.5,0.8), mean_diff != 1) %>% rename(vaccinated_population = population_size) %>% arrange(w) %>% mutate(growth_rate = factor(growth_rate, levels=c(-0.3, -0.15, 0 , 0.15, 0.3))), aes(x=associativity, y=R_relative)) +
  geom_line(aes(color=growth_rate, linetype=vaccinated_population)) +
  scale_linetype_manual(values=c("dotted", "dashed", "solid")) +
  facet_wrap(~scenario, labeller = label_both, nrow=1) +
  theme_bw() +
  scale_x_continuous(breaks = seq(-20,20,20)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank()) +
  ylab(expression(R["adjusted"] / R["unadjusted"])) +
  geom_text(data = data.frame(x = rep(min(R_by_gT_rel_A_simple$associativity),3), y=rep(max(R_by_gT_rel_A_simple$R_relative),3), label=c("D","E","F"), scenario = as.factor(c(0.25,0.5,0.75))), aes(x=x, y=y, label=label), fontface="bold", size=5)


# R_vacc_no_vacc_R_diff <- ggplot(R_by_gT_rel_A_simple %>% rename(vaccinated_population = population_size, R_unvaccinated = R1, R_vaccinated = R2) %>% arrange(w) %>% filter(associativity %in% c(-19:19), mean_diff != 1) %>% mutate(growth_rate = factor(growth_rate, levels=c(-0.3, -0.15, 0 , 0.15, 0.3))) %>% tidyr::gather(metric, "value", c(10:14,19,20)) %>% filter(growth_rate == 0.15, metric %in% c("overall_R", "R_vaccinated", "R_unvaccinated"), vaccinated_population %in% c(0.2,0.5,0.8)) %>% mutate(metric = factor(metric, levels=c("overall_R", "R_vaccinated", "R_unvaccinated"))), aes(x=associativity, y=value)) +
#   geom_line(aes(color = metric, linetype = vaccinated_population)) +
#   scale_linetype_manual(values=c("dotted", "dashed", "solid")) +
#   facet_wrap(~scenario, labeller = label_both, nrow=1) +
#   theme_bw() +
#   scale_x_continuous(breaks = seq(-20,20,20)) +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         strip.text.x = element_blank()) +
#   ylab("R") +
#   geom_text(data = data.frame(x = rep(min(R_by_gT_rel_A_simple$associativity),3), y=rep(max(R_by_gT_rel_A_simple$k2_normalised),3), label=c("G","H","I"), scenario = as.factor(c(0.25,0.5,0.75))), aes(x=x, y=y, label=label), fontface="bold", size=5)
# 
#   #xlab(expression(dissassortative %->% homogeneous %->% assortative)) 

R_vacc_no_vacc_k_diff_0.3 <- ggplot(R_by_gT_rel_A_simple %>% rename(vaccinated_population = population_size, k_unvaccinated = k1_normalised, k_vaccinated = k2_normalised) %>% arrange(w) %>% filter(associativity %in% c(-19:19), mean_diff != 1) %>% mutate(growth_rate = factor(growth_rate, levels=c(-0.3, -0.15, 0 , 0.15, 0.3))) %>% tidyr::gather(metric, "value", c(10:14,19,20)) %>% filter(growth_rate == 0.3, metric %in% c("k_unvaccinated", "k_vaccinated"), vaccinated_population %in% c(0.2,0.5,0.8)) %>% mutate(metric = factor(metric, levels=c( "k_vaccinated", "k_unvaccinated"))), aes(x=associativity, y=value)) +
  geom_line(aes(color = metric, linetype = vaccinated_population)) +
  scale_linetype_manual(values=c("dotted", "dashed", "solid")) +
  scale_color_manual(values = c("blue","red")) +
  facet_wrap(~scenario, labeller = label_both, nrow=1) +
  theme_bw() +
  scale_x_continuous(breaks = seq(-20,20,20)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank()) +
  ylab("k") +
  geom_text(data = data.frame(x = rep(min(R_by_gT_rel_A_simple$associativity),3), y=rep(max(R_by_gT_rel_A_simple$k1_normalised),3), label=c("G","H","I"), scenario = as.factor(c(0.25,0.5,0.75))), aes(x=x, y=y, label=label), fontface="bold", size=5) + 
  geom_text(data = data.frame(associativity=10, value=1, label="r = 0.3"), aes(x=associativity, y=value, label=label))


R_vacc_no_vacc_k_diff_0.3_neg <- ggplot(R_by_gT_rel_A_simple %>% rename(vaccinated_population = population_size, k_unvaccinated = k1_normalised, k_vaccinated = k2_normalised) %>% arrange(w) %>% filter(associativity %in% c(-19:19), mean_diff != 1) %>% mutate(growth_rate = factor(growth_rate, levels=c(-0.3, -0.15, 0 , 0.15, 0.3))) %>% tidyr::gather(metric, "value", c(10:14,19,20)) %>% filter(growth_rate == -0.3, metric %in% c("k_unvaccinated", "k_vaccinated"), vaccinated_population %in% c(0.2,0.5,0.8)) %>% mutate(metric = factor(metric, levels=c("k_vaccinated", "k_unvaccinated"))), aes(x=associativity, y=value)) +
  geom_line(aes(color = metric, linetype = vaccinated_population)) +
  scale_color_manual(values = c("blue","red")) +
  scale_linetype_manual(values=c("dotted", "dashed", "solid")) +
  facet_wrap(~scenario, labeller = label_both, nrow=1) +
  theme_bw() +
  scale_x_continuous(breaks = seq(-20,20,20)) +
  theme(axis.text.x = element_blank(),
        strip.text.x = element_blank()) +
  xlab(expression(dissassortative %->% homogeneous %->% assortative)) +
  ylab("k") +
  geom_text(data = data.frame(x = rep(min(R_by_gT_rel_A_simple$associativity),3), y=rep(max(R_by_gT_rel_A_simple$k1_normalised),3), label=c("J","K","L"), scenario = as.factor(c(0.25,0.5,0.75))), aes(x=x, y=y, label=label), fontface="bold", size=5)+ 
  geom_text(data = data.frame(associativity=10, value=1, label="r = -0.3"), aes(x=associativity, y=value, label=label))



#vaccine_plot2 <- ggarrange(vaccine_impact, R_vacc_no_vacc_A,nrow=2)
#ggsave(plot=vaccine_plot2, filename=paste("graphs/error_R_vacc_no_vacc_4.png"), width=15, height=10)


vaccine_plot3 <- ggarrange(vaccine_impact, R_vacc_no_vacc_A, R_vacc_no_vacc_k_diff_0.3, R_vacc_no_vacc_k_diff_0.3_neg, nrow=4)
ggsave(plot=vaccine_plot3, filename=paste("graphs/error_R_vacc_no_vacc_zoom_11.png"), width=10, height=12)


R_by_gT_rel_A_test <- R_by_gT_rel_A_simple %>% rename(vaccinated_population = population_size) %>% arrange(w) %>% filter(associativity != -20) %>% mutate(growth_rate = factor(growth_rate, levels=c(-0.3, -0.15, 0 , 0.15, 0.3))) %>% filter(growth_rate == 0.15, vaccinated_population == 0.5, associativity %in% c(-19,0,19), mean_diff == 0.5)


## Application to UK situation 
shape = 10.56
rate = 1.85

shape/rate
shape/(rate^2)

normalise <- function(vector) return(vector/sum(vector))

weighted_gt_generator <- function(distribution_type, distribution_param_1, distribution_param_2, mean_difference=1, truncation_point=0.999, risk_trans, pop_vec, isol=FALSE, asymp=FALSE, plot){
  gt_no_isol <- gT_generator("gamma", distribution_param_1, distribution_param_2)$gT_dist
  gt_asymp   <-  gT_proposal(mean_difference, "gamma", distribution_param_1, distribution_param_2)
  gt_isol <- gT_proposal_trunc(truncation_point, "gamma", distribution_param_1, distribution_param_2)
  
  if(isol == FALSE && asymp == FALSE) gt_weighted = gt_no_isol
  
  if(isol == TRUE && asymp == FALSE){
    gt_weighted = normalise(c(t(pop_vec*risk_trans)) %*% matrix(c(gt_no_isol, gt_isol), nrow=2, byrow=TRUE))
  }
  
  if(isol == FALSE && asymp == TRUE){
    gt_weighted = normalise(c(t(pop_vec*risk_trans)) %*% matrix(c(gt_no_isol, gt_asymp), nrow=2, byrow=TRUE))
  }
  
  if(isol == TRUE && asymp == TRUE){
    gt_weighted = normalise(c(t(pop_vec*risk_trans)) %*% matrix(c(gt_no_isol, gt_isol, gt_asymp), nrow=3, byrow=TRUE))
  }
  
  gt_weighted <- c(rowSums(matrix(gt_weighted[-length(gt_weighted)], ncol=100, byrow=TRUE)))
  
  if(plot ==TRUE) plot(gt_weighted, type="l", xlim=c(0,20))
  
  return(gt_weighted)
}

deaths_data <- read.csv("data/UK_deaths_time_series.csv")

deaths_data_smoothed <- deaths_data %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

ggplot(deaths_data_smoothed, aes(x=date, y=deaths)) +
  geom_line()


gt_control    <- weighted_gt_generator(distribution_type="gamma", distribution_param_1=shape, distribution_param_2=rate, mean_difference=1, truncation_point=1, risk_trans=c(1), pop_vec=c(1), isol=FALSE, asymp=FALSE, plot=TRUE)
gt_isol_optim <- weighted_gt_generator(distribution_type="gamma", distribution_param_1=shape, distribution_param_2=rate, mean_difference=1, truncation_point=0.3, risk_trans=c(1,0.3), pop_vec=c(1-0.63*0.75,0.63*0.75), isol=TRUE, asymp=FALSE, plot=TRUE)
gt_isol_pesim <- weighted_gt_generator(distribution_type="gamma", distribution_param_1=shape, distribution_param_2=rate, mean_difference=1, truncation_point=0.7, risk_trans=c(1,0.7), pop_vec=c(1-0.63*0.25,0.63*0.25), isol=TRUE, asymp=FALSE, plot=TRUE)
gt_asym_optim <- weighted_gt_generator(distribution_type="gamma", distribution_param_1=shape, distribution_param_2=rate, mean_difference=0.5, truncation_point=1, risk_trans=c(1,0.5), pop_vec=c(0.63, 0.37), isol=FALSE, asymp=TRUE, plot=TRUE)
gt_asym_pesim <- weighted_gt_generator(distribution_type="gamma", distribution_param_1=shape, distribution_param_2=rate, mean_difference=  2, truncation_point=1, risk_trans=c(1,2),   pop_vec=c(0.63, 0.37), isol=FALSE, asymp=TRUE, plot=TRUE)
gt_isol_asym_optim <- weighted_gt_generator(distribution_type="gamma", distribution_param_1=shape, distribution_param_2=rate, mean_difference=0.5, truncation_point=0.3, risk_trans=c(1, 0.3, 0.5), pop_vec=c(0.63*0.25 ,0.63*0.75, 0.37), isol=TRUE, asymp=TRUE, plot=TRUE)
gt_isol_asym_pesim <- weighted_gt_generator(distribution_type="gamma", distribution_param_1=shape, distribution_param_2=rate, mean_difference=2,   truncation_point=0.7, risk_trans=c(1, 0.7,   2), pop_vec=c(0.63*0.75, 0.63*0.25, 0.37), isol=TRUE, asymp=TRUE, plot=TRUE)

gt_asymp_control_1sd_up <- weighted_gt_generator(distribution_type="gamma", distribution_param_1=shape*(1.1^2)/0.9, distribution_param_2=rate*1.1/0.9, mean_difference=1, truncation_point=0.99999, risk_trans=c(1), pop_vec=c(1), isol=FALSE, asymp=FALSE, plot=TRUE)
gt_asymp_control_1sd_down <- weighted_gt_generator(distribution_type="gamma", distribution_param_1=shape*(0.9^2)/1.1, distribution_param_2=rate*0.9/1.1, mean_difference=1, truncation_point=0.99999, risk_trans=c(1), pop_vec=c(1), isol=FALSE, asymp=FALSE, plot=TRUE) 

plot(gt_control, type="l", xlim=c(0,10))
lines(gt_isol_optim, type="l")
lines(gt_isol_pesim, type="l")
lines(gt_asym_optim, type="l", col="red")
lines(gt_asym_pesim, type="l")
lines(gt_isol_asym_optim, type="l")
lines(gt_isol_asym_pesim, type="l")

R_asymp.control <- estimate_R(incid = deaths_data$deaths,
                              method = "non_parametric_si",
                              config = make_config(list(si_distr = c(0,gt_control))))$R[,c("Quantile.0.25(R)","Median(R)", "Quantile.0.975(R)")] %>%
  set_colnames(c("asymp.control.lower", "asymp.control.median", "asymp.control.upper"))

R_isol.control <- estimate_R(incid = deaths_data$deaths,
                             method = "non_parametric_si",
                             config = make_config(list(si_distr = c(0,gt_control))))$R[,c("Quantile.0.25(R)","Median(R)", "Quantile.0.975(R)")] %>%
  set_colnames(c("isol.control.lower", "isol.control.median", "isol.control.upper"))

R_isol.asymp.control <- estimate_R(incid = deaths_data$deaths,
                                   method = "non_parametric_si",
                                   config = make_config(list(si_distr = c(0,gt_control))))$R[,c("Quantile.0.25(R)","Median(R)", "Quantile.0.975(R)")] %>%
  set_colnames(c("isol_asymp.control.lower", "isol_asymp.control.median", "isol_asymp.control.upper"))


R_isol.optim <- estimate_R(incid = deaths_data$deaths,
                           method = "non_parametric_si",
                           config = make_config(list(si_distr = c(0,gt_isol_optim))))$R[,c("Quantile.0.25(R)","Median(R)", "Quantile.0.975(R)")] %>%
  set_colnames(c("isol.optimistic.lower", "isol.optimistic.median", "isol.optimistic.upper"))


R_isol.pesim <- estimate_R(incid = deaths_data$deaths,
                           method = "non_parametric_si",
                           config = make_config(list(si_distr = c(0,gt_isol_pesim))))$R[,c("Quantile.0.25(R)","Median(R)", "Quantile.0.975(R)")]%>%
  set_colnames(c("isol.pessimistic.lower", "isol.pessimistic.median", "isol.pessimistic.upper"))


R_asymp.optim <- estimate_R(incid = deaths_data$deaths,
                            method = "non_parametric_si",
                            config = make_config(list(si_distr = c(0,gt_asym_optim))))$R[,c("Quantile.0.25(R)","Median(R)", "Quantile.0.975(R)")]%>%
  set_colnames(c("asymp.optimistic.lower", "asymp.optimistic.median", "asymp.optimistic.upper"))


R_asymp.pesim <- estimate_R(incid = deaths_data$deaths,
                            method = "non_parametric_si",
                            config = make_config(list(si_distr = c(0,gt_asym_pesim))))$R[,c("Quantile.0.25(R)","Median(R)", "Quantile.0.975(R)")]%>%
  set_colnames(c("asymp.pessimistic.lower", "asymp.pessimistic.median", "asymp.pessimistic.upper"))


R_isol.asymp.optim <- estimate_R(incid = deaths_data$deaths,
                                 method = "non_parametric_si",
                                 config = make_config(list(si_distr = c(0,gt_isol_asym_optim))))$R[,c("Quantile.0.25(R)","Median(R)", "Quantile.0.975(R)")]%>%
  set_colnames(c("isol_asymp.optimistic.lower", "isol_asymp.optimistic.median", "isol_asymp.optimistic.upper"))


R_isol.asymp.pesim <- estimate_R(incid = deaths_data$deaths,
                                 method = "non_parametric_si",
                                 config = make_config(list(si_distr = c(0,gt_isol_asym_pesim))))$R[,c("Quantile.0.25(R)","Median(R)", "Quantile.0.975(R)")]%>%
  set_colnames(c("isol_asymp.pessimistic.lower", "isol_asymp.pessimistic.median", "isol_asymp.pessimistic.upper"))

t_start <- seq(as.Date("2020-02-29"), as.Date("2020-02-29")+nrow(R_asymp.control)-1, 1)

overall_df <- cbind(t_start, R_isol.control, R_asymp.control, R_isol.asymp.control, R_isol.optim, R_isol.pesim, R_asymp.optim, R_asymp.pesim, R_isol.asymp.optim, R_isol.asymp.pesim)

overall_df_final <- overall_df %>% tidyr::gather("metric", "value", 2:ncol(overall_df)) %>% rowwise() %>%
  mutate(scenario = str_split(metric,"\\.")[[1]][1],
         sensitivity = str_split(metric,"\\.")[[1]][2], 
         interval = str_split(metric,"\\.")[[1]][3]) %>%
  select(-metric) %>%
  tidyr::spread(interval, "value") %>%
  group_by(scenario, sensitivity) #%>%

R_UK_symp_plot <- ggplot(overall_df_final %>% filter(scenario != "isol"), aes(x = t_start)) +
  geom_line(aes(y=median, color=sensitivity)) +
  scale_color_manual(values=c("blue", "green", "red")) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=sensitivity), alpha=0.3) +
  scale_fill_manual(values=c("blue", "green", "red")) +
  #geom_ribbon(aes(ymin=control.lower, ymax=control.upper), alpha=0.1) +
  facet_wrap(~scenario, nrow=2) +
  theme_bw() +
  scale_y_continuous(limits=c(0,2)) +
  theme(axis.title.x = element_blank()) +
  labs(y = "Derived R") +
  ggtitle("n-SARS-CoV2 infection, United Kingdom") +
  theme(strip.text.x = element_blank()) +
  geom_text(data = data.frame(x=rep(min(overall_df_final$t_start),2), y=rep(2,2), label=c("A","B"), scenario = overall_df_final %>% ungroup() %>% filter(scenario != "isol") %>% select(scenario) %>% unique()), aes(x=x,y=y, label=label), fontface="bold", size=5)

ggsave(plot=R_UK_symp_plot, filename=paste("graphs/R_UK_symp_plot6.png"), width=10, height=5)

R_UK_isolation_plot <- ggplot(overall_df_final %>% filter(scenario=="isol"), aes(x = t_start)) +
  geom_line(aes(y=median, color=sensitivity)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=sensitivity), alpha=0.1) +
  scale_color_manual(values=c("blue", "green", "red")) +
  scale_fill_manual(values=c("blue", "green", "red")) +
  facet_wrap(~scenario) +
  theme_bw() +
  scale_y_continuous(limits=c(0,2)) +
  theme(axis.title.x = element_blank()) +
  labs(y = "Derived R") +
  ggtitle("n-SARS-CoV2 infection, United Kingdom") +
  theme(strip.text.x = element_blank(),
        legend.title = element_blank()) +
  geom_text(data = data.frame(x=rep(min(overall_df_final$t_start)), y=rep(2), label=c("B"), scenario = "isol" %>% unique()), aes(x=x,y=y, label=label), fontface="bold", size=5)



## Application to Guinea EVD
missed_cases <- 0.64

mu <- 15.3
sd <- 9.3

alpha <- (mu/sd)^2
beta <- (mu/sd^2)

gt_control_evd <-  weighted_gt_generator(distribution_type="gamma", distribution_param_1=2.7066, distribution_param_2=0.1769, mean_difference=1, risk_trans=c(1), pop_vec=c(1), isol=FALSE, asymp=FALSE, plot=TRUE)
gt_hosp_evd <-     weighted_gt_generator(distribution_type="gamma", distribution_param_1=2.7066, distribution_param_2=0.1769, mean_difference=1, truncation_point=0.565, risk_trans=c(1,0.565), pop_vec=c(missed_cases, 1-missed_cases), isol=TRUE, asymp=FALSE, plot=TRUE)

plot(gt_control_evd, type="l")
lines(gt_hosp_evd, type="l")

evd_dates2 <- data.frame(date = seq(as.Date("2014-01-01"), as.Date("2015-09-25"), 1))

evd_data2 <- read.csv("data/guinea_ebola_cases.csv") %>%
  mutate(date = as.Date(DateOnsetInferred, format="%d/%m/%Y")) %>%
  select(c("Country", "date")) %>%
  filter(Country == "Guinea") %>%
  na.omit() %>%
  group_by(date) %>%
  summarise(count=n())

evd_data_frame <- left_join(evd_dates2, evd_data2, "date") %>%
  mutate(count = ifelse(is.na(count)==TRUE, 0, count)) %>% filter(date >= as.Date("2014-03-01"), date <= as.Date("2015-07-01"))

R_control_evd <- estimate_R(incid = evd_data_frame$count,
                            method = "non_parametric_si",
                            config = make_config(list(si_distr = c(0,gt_control_evd))))$R[,c("Quantile.0.25(R)","Median(R)", "Quantile.0.975(R)")] %>%
  set_colnames(c("control.lower", "control.median", "control.upper"))

R_isol_evd <- estimate_R(incid = evd_data_frame$count,
                         method = "non_parametric_si",
                         config = make_config(list(si_distr = c(0,gt_hosp_evd))))$R[,c("Quantile.0.25(R)","Median(R)", "Quantile.0.975(R)")] %>%
  set_colnames(c("isol.lower", "isol.median", "isol.upper"))

plot(R_control_evd$control.lower, type="l") 


t_start_evd <- evd_data_frame$date[1:(nrow(evd_data_frame)-7)]

overall_df_evd <- cbind(t_start_evd, R_control_evd, R_isol_evd) %>% 
  tidyr::gather("metric", "R", c(2:7)) %>% rowwise() %>%
  mutate(scenario = str_split(metric, "\\.")[[1]][1],
         bound = str_split(metric, "\\.")[[1]][2]) %>%
  select(-metric) %>%
  tidyr::spread(bound, R) %>%
  mutate(lower = ifelse(scenario == "control", NA, lower),
         upper = ifelse(scenario == "control", NA, upper))

R_EVD_Guinea_plot <- ggplot(overall_df_evd, aes(x=as.Date(t_start_evd))) + 
  geom_line(aes(y=median, color = scenario)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = scenario), alpha=0.2) +
  #scale_fill_manual(values=c("blue", "red")) + 
  ylim(c(0,4)) + 
  scale_fill_manual(values = c("blue", "orange")) +
  scale_color_manual(values = c("blue", "orange")) +
  theme_bw() + 
  labs(x="date", y="Derived R") +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) +
  scale_x_date(date_labels =  "%b %Y") +
  ggtitle("Ebola Virus Disease, Guinea 2014")+
  geom_text(data = data.frame(x=rep(min(overall_df_evd$t_start_evd)), y=rep(4), label=c("A"), scenario = "isol" %>% unique()), aes(x=x,y=y, label=label), fontface="bold", size=5)


applied_isol_plot <- ggarrange(R_EVD_Guinea_plot,R_UK_isolation_plot, nrow=2)  
ggsave(plot=applied_isol_plot, filename=paste("graphs/R_isolation_UK_Guinea6.png"), width=10, height=5)


## Extension application with uncertain serial interval
R_asymp.control
R_asymp.optim
R_asymp.pesim

R_asymp.control.unc_up <- estimate_R(incid = deaths_data$deaths,
                              method = "non_parametric_si",
                              config = make_config(list(si_distr = c(0,gt_asymp_control_1sd_up))))$R[,c("Quantile.0.25(R)","Median(R)", "Quantile.0.975(R)")] %>%
  set_colnames(c("asymp.control_up.lower", "asymp.control_up.median", "asymp.control_up.upper"))

R_asymp.control.unc_down <- estimate_R(incid = deaths_data$deaths,
                              method = "non_parametric_si",
                              config = make_config(list(si_distr = c(0,gt_asymp_control_1sd_down))))$R[,c("Quantile.0.25(R)","Median(R)", "Quantile.0.975(R)")] %>%
  set_colnames(c("asymp.control_down.lower", "asymp.control_down.median", "asymp.control_down.upper"))

t_start <- seq(as.Date("2020-02-29"), as.Date("2020-02-29")+nrow(R_asymp.control)-1, 1)

overall_df_unc <- cbind(t_start, R_asymp.control.unc_up, R_asymp.control.unc_down, R_asymp.optim, R_asymp.pesim)

overall_df_final_unc <- overall_df_unc %>% tidyr::gather("metric", "value", 2:ncol(overall_df_unc)) %>% rowwise() %>%
  mutate(scenario = str_split(metric,"\\.")[[1]][1],
         sensitivity = str_split(metric,"\\.")[[1]][2], 
         interval = str_split(metric,"\\.")[[1]][3]) %>%
  select(-metric) %>%
  tidyr::spread(interval, "value") %>%
  group_by(scenario, sensitivity) #%>%

R_UK_symp_plot_unc <- ggplot(overall_df_final_unc , aes(x = t_start)) +
  geom_line(aes(y=median, color=sensitivity)) +
  scale_color_manual(values=c("blue","orange", "green", "red")) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=sensitivity), alpha=0.1) +
  scale_fill_manual(values=c("blue","orange", "green", "red")) +
  #geom_ribbon(aes(ymin=control.lower, ymax=control.upper), alpha=0.1) +
  facet_wrap(~scenario, nrow=2) +
  theme_bw() +
  scale_y_continuous(limits=c(0,2)) +
  theme(axis.title.x = element_blank()) +
  labs(y = "Derived R") +
  ggtitle("n-SARS-CoV2 infection, United Kingdom") +
  theme(strip.text.x = element_blank()) 

ggsave(plot=R_UK_symp_plot_unc, filename=paste("graphs/R_UK_symp_plot_unc6.png"), width=10, height=5)
















#applications_plot <- ggarrange(R_EVD_Guinea_plot, R_UK_plot, nrow=2)
#ggsave(plot=applications_plot, filename=paste("graphs/R_Guinea_UK_plot.png"), width=15, height=10)

#####

gt_no_isol <- gT_generator("gamma", 10.56, 1.85)$gT_dist
gt_isol <- gT_proposal_trunc(0.2,"gamma", 10.56, 1.85)
gt_asymp <- gT_proposal(0.5, "gamma", 10.56, 1.85)

risk_matrix_UK <- relative_risk_mat_generator(risk_trans = c(1,0.5,0.25),contact_matrix = matrix(c(0.33,0.33,0.33,
                                                                                                   0.33,0.33,0.33,
                                                                                                   0.33,0.33,0.33), nrow=3, byrow=T))

UK_data <- read.csv(file="PhD/Data/covid_cases_data") %>% 
  select(date, cases.daily) %>% arrange(desc(row_number())) %>%
  mutate(date = as.Date(date), 
         cases.daily = as.numeric(cases.daily), 
         cases.weekly = rollsum(cases.daily,7, fill=NA)) %>%
  mutate(growth.rate.weekly = (cases.weekly - lag(cases.weekly))/lag(cases.weekly),
         growth.rate.daily = (cases.daily - lag(cases.daily))/lag(cases.daily)) %>% 
  filter(date >= as.Date("2020-03-15"), date <= as.Date("2021-01-14")) %>% rowwise() %>%
  mutate(R.simple.daily = R_solver_nd(w = growth.rate.daily, gT_list = cbind(gt_no_isol), risk_matrix = matrix(1), day_granularity = 0.01)[1],
         R.simple.weekly = R_solver_nd(w = growth.rate.weekly, gT_list = cbind(gt_no_isol), risk_matrix = matrix(1), day_granularity = 0.01)[1],
         R.simple.isol.weekly = R_solver_nd(w = growth.rate.weekly, gT_list = cbind(gt_isol), risk_matrix = matrix(1), day_granularity = 0.01)[1],
         R.true.daily = R_solver_nd(w = growth.rate.daily, gT_list = cbind(gt_no_isol, gt_isol, gt_asymp), risk_matrix = risk_matrix_UK, day_granularity = 0.01)[1],
         R.true.weekly = R_solver_nd(w = growth.rate.weekly, gT_list = cbind(gt_no_isol, gt_isol, gt_asymp), risk_matrix = risk_matrix_UK, day_granularity = 0.01)[1]) 

UK_data_clean <- UK_data %>%
  tidyr::gather("metric","value", 6:10)

UK_vaccine_projection <- read.csv("PhD/20210113-100503-f12b78ef/projected_vaccinations.csv") %>%
  mutate(date = as.Date(date),
         scenario = as.factor(scenario))

UK_projection <- read.csv("PhD/20210113-100503-f12b78ef/projected_infections.csv") %>%
  mutate(date = as.Date(date),
         growth.rate.daily = (infections - lag(infections))/lag(infections)) %>%
  filter(is.na(growth.rate.daily) == FALSE) 


UK_infections_vaccinations <- full_join(UK_projection, UK_vaccine_projection, by=c("date", "scenario"))

R_inference_vaccines <- function(data_frame, distribution_type, distribution_param_1, distribution_param_2, one_dose_eff=0.48, two_dose_eff=0.6){
  gt_no_vacc <- gT_generator(distribution_type, distribution_param_1, distribution_param_2)$gT_dist
  gt_one_vacc <- gT_proposal(mean_difference = sqrt(1-one_dose_eff), distribution_type, distribution_param_1, distribution_param_2)
  gt_two_vacc <- gT_proposal(mean_difference = sqrt(1-two_dose_eff), distribution_type, distribution_param_1, distribution_param_2)
  
  new_df <- data_frame %>% na.omit() %>% rowwise() %>%
    mutate(R.simple.daily = R_solver_nd(w = growth.rate.daily, gT_list = cbind(gt_no_isol), risk_matrix = matrix(1), day_granularity = 0.01)[1],
           R.multi.daily =  R_solver_nd(w = growth.rate.daily, gT_list = cbind(gt_no_isol, gt_one_vacc, gt_two_vacc), 
                                        risk_matrix = relative_risk_mat_generator(risk_trans = c(1, 1-one_dose_eff, 1-two_dose_eff), contact_matrix = matrix(c(1-prop_1_dose-prop_2_dose,1-prop_1_dose-prop_2_dose,1-prop_1_dose-prop_2_dose,
                                                                                                                                                               prop_1_dose, prop_1_dose, prop_1_dose,
                                                                                                                                                               prop_2_dose, prop_2_dose, prop_2_dose), nrow=3, byrow=T)), day_granularity = 0.01)[1])
  
  return(new_df)
  
}

R_inference_symp_isol <- function(data_frame, mean_difference = 0.5, truncation_point = 0.5, prop_isol = 0.5, prop_symp = 0.63, distribution_type, distribution_param_1, distribution_param_2){
  gt_symp_no_isol <- gT_generator(distribution_type, distribution_param_1, distribution_param_2)$gT_dist
  
  new_df <- data_frame %>% na.omit() %>% rowwise() %>%
    mutate(R.simple.daily = R_solver_nd(w = growth.rate.daily, gT_list = cbind(gt_no_isol), risk_matrix = matrix(1), day_granularity = 0.01)[1])
  
  gT_list_symp_asymp = cbind(gt_symp_no_isol, gt_asymp)
  contact_matrix_symp_asymp = matrix(rep(c(prop_symp, 1-prop_symp),2),nrow=2, byrow=TRUE)
  
  gT_list_isol_no_isol = cbind(gt_symp_no_isol, gt_isol)
  contact_matrix_isol_no_isol = matrix(rep(c(prop_symp*(1-prop_isol), 1-prop_symp*(1-prop_isol)), 2), nrow=2, byrow=T)
  
  gT_list_symp_isol = cbind(gt_symp_no_isol, gt_isol, gt_asymp)
  contact_matrix_symp_isol = matrix(rep(c(prop_symp*(1-prop_isol), prop_symp*prop_isol, (1-prop_symp)),3), nrow=3, byrow=T)
  
  new_df <- data_frame %>% na.omit() %>% rowwise() %>%
    mutate(R.simple.daily = R_solver_nd(w = growth.rate.daily, gT_list = cbind(gt_no_isol), risk_matrix = matrix(1), day_granularity = 0.01)[1],
           R.symp.asymp   = R_solver_nd(w = growth.rate.daily, gT_list = gT_list_symp_asymp, risk_matrix = relative_risk_mat_generator(risk_trans = c(1, mean_difference), contact_matrix = contact_matrix_symp_asymp), day_granularity = 0.01)[1],
           R.isol.no.isol = R_solver_nd(w = growth.rate.daily, gT_list = gT_list_isol_no_isol, risk_matrix = relative_risk_mat_generator(risk_trans = c(1, truncation_point), contact_matrix = contact_matrix_isol_no_isol), day_granularity = 0.01)[1],
           R.symp.isol    = R_solver_nd(w = growth.rate.daily, gT_list = gT_list_symp_isol, risk_matrix = relative_risk_mat_generator(risk_trans = c(1, truncation_point, mean_difference), contact_matrix = contact_matrix_symp_isol), day_granularity = 0.01)[1])
  
  
  return(new_df)
  
}


R_output_historic <- R_inference_symp_isol(data_frame=UK_projection %>% filter(scenario=="control"), mean_difference = 0.5, truncation_point = 0.5, prop_isol = 0.5, prop_symp = 0.63, distribution_type="gamma", distribution_param_1=10.56, distribution_param_2=1.85) %>% ungroup() %>%
  mutate(rolling.simple=rollapply(R.simple.daily,7,mean, fill=NA),
         rolling.symp.asymp=rollapply(R.symp.asymp,7,mean, fill=NA),
         rolling.isol.no.isol=rollapply(R.isol.no.isol,7,mean, fill=NA),
         rolling.symp.isol=rollapply(R.symp.isol,7,mean, fill=NA)) %>%
  tidyr::gather("multi.scenario", "rolling.multi", 11:13) #%>%
tidyr::gather("simple", )

ggplot(UK_projection, aes(x=date, y=infections))+
  geom_line()

ggplot(data=R_output_historic, aes(x=date, y=infections))+
  geom_line()

comparitor_plot <- ggplot(data=R_output_historic, aes(x=date)) +
  geom_line(aes(y=rolling.multi, color=multi.scenario)) +
  geom_line(aes(y=rolling.simple, color="simple.scenario")) +
  scale_color_manual(values = c("red", "green", "blue", "black")) +
  theme_bw() +
  facet_wrap(~multi.scenario) + 
  labs(y = "Estimated R") +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) 

ratio_plot <- ggplot(data=R_output_historic, aes(x=date, y=rolling.multi/rolling.simple)) +
  geom_line(aes(color=multi.scenario)) +
  scale_color_manual(values = c("red", "green", "blue")) +
  theme_bw() +
  facet_wrap(~multi.scenario) + 
  labs(y = "R.true/R.simple") + 
  theme(legend.title = element_blank(),
        legend.position = "none")

UK_historic_plot <- ggarrange(comparitor_plot, ratio_plot, nrow=2)

ggsave(plot=UK_historic_plot, filename=paste("PhD/Graphs/new_graphs/error_R_UK_renewal.png"), width=10, height=7)



R_output_UK <- R_inference_vaccines(data_frame=UK_infections, distribution_type="gamma", distribution_param_1=10.56, distribution_param_2=1.85, one_dose_eff=0.48, two_dose_eff=0.6)  

infections_plot <- ggplot(UK_infections_vaccinations, aes(x=date, y=infections)) + 
  geom_line(aes(color=scenario)) + 
  theme_bw()

scenarios_plot <- ggplot(UK_infections_vaccinations %>% filter(date >= "2021-01-13"), aes(x=date, y=infections)) + 
  geom_line(aes(color=scenario)) + 
  theme_bw()

R_plot <- ggplot(R_output_UK %>% filter(date >= "2021-01-15", infections >100), aes(x = date, y=R.simple.daily)) +
  geom_line(aes(color = scenario)) +
  theme_bw()


R_ratio_plot <- ggplot(R_output_UK %>% filter(date >= "2021-01-15", infections >100), aes(x = date, y=R.simple.daily/R.multi.daily)) +
  geom_line(aes(color = scenario)) +
  theme_bw()


ggarrange(scenarios_plot, R_plot, nrow=2)

#







for(i in 1:nrow(evd_data)) evd_data$max[i] = max(evd_data$total_cases[1:i])

evd_data_cleaned <- evd_data %>%
  mutate(cases_clean = pmax(total_cases, max)) 

evd_data_daily <- data_frame(report_date=seq(as.Date("2014-03-31"), as.Date("2016-03-03"), 1)) %>%
  left_join(., evd_data_cleaned) 

for(i in 2:nrow(evd_data_daily)){
  if(is.na(evd_data_daily$cases_clean[i]) == TRUE){
    lower_value <- max(evd_data_daily$cases_clean[1:(i-1)])
    upper_value <- min(evd_data_daily$cases_clean[(i+1):nrow(evd_data_daily)],na.rm = TRUE)
    
    lower_value_index <- which(evd_data_daily$cases_clean == lower_value)[length(which(evd_data_daily$cases_clean == lower_value))]
    upper_value_index <- which(evd_data_daily$cases_clean == upper_value)[1]
    
    incriment <- (upper_value/lower_value)^(1/(upper_value_index - lower_value_index))
    
    evd_data_daily$cases_clean[i] <- as.numeric(lower_value * incriment)
  }
}

evd_data_daily_inc <- evd_data_daily %>%
  mutate(cases_inc = cases_clean - lag(cases_clean),
         rolling_cases_inc = as.numeric((rollapply(cases_inc, 7, mean, fill=NA)))) %>%
  filter(is.na(rolling_cases_inc) == FALSE)


gt1 <- gT_generator("gamma", 5, 1)$gT_dist
gt2 <- gT_proposal(0.5, "gamma", 5, 1)

contact_matrix_dis <- contact_matrix_generator(population_2_size = 0.8, associativity_vec = seq(-20,20,1), i=2)
contact_matrix_homo <- contact_matrix_generator(population_2_size = 0.8, associativity_vec = seq(-20,20,1), i=21)
contact_matrix_30 <- contact_matrix_generator(population_2_size = 0.8, associativity_vec = seq(-20,20,1), i=30)
contact_matrix_35 <- contact_matrix_generator(population_2_size = 0.8, associativity_vec = seq(-20,20,1), i=35)
contact_matrix_ass <- contact_matrix_generator(population_2_size = 0.8, associativity_vec = seq(-20,20,1), i=40)

risk_matrix_dis <- relative_risk_mat_generator(risk_trans=c(1,0.0625), contact_matrix = contact_matrix_dis) 
risk_matrix_homo <- relative_risk_mat_generator(risk_trans=c(1,0.0625), contact_matrix = contact_matrix_homo) 
risk_matrix_30 <- relative_risk_mat_generator(risk_trans=c(1,0.0625), contact_matrix = contact_matrix_30) 
risk_matrix_35 <- relative_risk_mat_generator(risk_trans=c(1,0.0625), contact_matrix = contact_matrix_35) 
risk_matrix_ass <- relative_risk_mat_generator(risk_trans=c(1,0.0625), contact_matrix = contact_matrix_ass) 

R_solver_nd(w=0.15, gT_list = cbind(gt1, gt2), risk_matrix = relative_risk_mat_generator(risk_trans=c(1,0.0625), contact_matrix = contact_matrix_dis),  day_granularity = 0.01)
R_solver_nd(w=0.15, gT_list = cbind(gt1, gt2), risk_matrix = relative_risk_mat_generator(risk_trans=c(1,0.0625), contact_matrix = contact_matrix_homo),  day_granularity = 0.01)
R_solver_nd(w=0.15, gT_list = cbind(gt1, gt2), risk_matrix = relative_risk_mat_generator(risk_trans=c(1,0.0625), contact_matrix = contact_matrix_30),  day_granularity = 0.01)
R_solver_nd(w=0.15, gT_list = cbind(gt1, gt2), risk_matrix = relative_risk_mat_generator(risk_trans=c(1,0.0625), contact_matrix = contact_matrix_35),  day_granularity = 0.01)
R_solver_nd(w=0.15, gT_list = cbind(gt1, gt2), risk_matrix = relative_risk_mat_generator(risk_trans=c(1,0.0625), contact_matrix = contact_matrix_ass),  day_granularity = 0.01)

risk_matrix_homo * 1.862
risk_matrix_30 * 2.002
risk_matrix_35 * 2.009
risk_matrix_ass * 2.011



c(17:07)

R_by_w_gT_trunc <- R_by_gT_A_w(w_vec=round(seq(-0.3,0.3,0.15),3), distribution_type="gamma", distribution_param_1=10.56, distribution_param_2=1, NPI = "trunc",
                               associativity_vec = c(1), population_sizes_vec = c(0.2,0.8),
                               gT_diff = "trunc", mean_diff_vec = seq(1), relative_inf_vec = round(seq(0.01,0.99,0.02),3),
                               day_granularity=0.01, log_associativity = T)


ptm2 <- proc.time()
R_by_gT_rel <- R_by_gT_A_w(w_vec=round(seq(-0.3,0.3,0.15),3), distribution_type="gamma", distribution_param_1=10.56, distribution_param_2=1, NPI="non",
                           associativity_vec = c(1), population_sizes_vec = c(0.2,0.8), gT_diff = "relative", relative_inf_vec = round(seq(0.01, 0.99, 0.02),3), 
                           mean_diff_vec = c(0.5,1,2), day_granularity=0.01, log_associativity = T)


ptm3 <- proc.time()

R_by_gT_rel_A <- R_by_gT_A_w(w_vec=round(seq(-0.3,0.3,0.15),3), distribution_type="gamma", distribution_param_1=5, distribution_param_2=1, NPI="vaccine",
                             associativity_vec = c(-20:20), population_sizes_vec = c(0.2,0.5,0.8),
                             gT_diff = "relative", mean_diff_vec = seq(0.25, 1.0, 0.25), relative_inf_vec = (seq(0.25, 1.0, 0.25))^2,
                             day_granularity=0.01, log_associativity = F) %>%
  mutate(scenario = ifelse(mean_diff^2 == relative_inf, mean_diff, NA)) %>%
  filter(is.na(scenario)==FALSE)

w_R_by_base_R <- function(w=0.1, distribution_type="gamma", distribution_param_1 =10.56, distribution_param_2 = 1, day_granularity = 0.01, 
                          R_group1, mean_diff_vec, associativity_vec, relative_inf_vec, population_sizes_vec, diff_type){
  
  output_storage_w <- array(NA, dim=c(length(population_sizes_vec), length(associativity_vec), length(mean_diff_vec), length(relative_inf_vec)), dimnames=list(population_sizes_vec,associativity_vec,mean_diff_vec, relative_inf_vec))
  output_storage_R <- array(NA, dim=c(length(population_sizes_vec), length(associativity_vec), length(mean_diff_vec), length(relative_inf_vec)), dimnames=list(population_sizes_vec,associativity_vec,mean_diff_vec, relative_inf_vec))
  
  for(i in 1:length(population_sizes_vec)){
    for(j in 1:length(associativity_vec)){
      for(k in 1:length(mean_diff_vec)){
        for(l in 1:length(relative_inf_vec)){
          output_storage <- w_R_solver_base_R(R_group1, mean_diff=mean_diff_vec[k], rel_inf=relative_inf_vec[l], population_sizes=population_sizes_vec[i], associativity=associativity_vec[j], diff_type=diff_type, distribution_type, distribution_param_1, distribution_param_2, day_granularity = 0.01)
          output_storage_w[i,j,k,l] <- output_storage$w
          output_storage_R[i,j,k,l] <- output_storage$R
        }
      }
    }
  }
  
  output_w <- as.data.frame.table(output_storage_w) %>% set_colnames(c("population_size", "associativity", "mean_diff","relative_inf", "w")) %>% as.data.frame() %>%
    mutate(population_size = as.numeric(as.character(population_size)), 
           associativity = as.numeric(as.character(associativity)),
           mean_diff = as.numeric(as.character(mean_diff)),
           relative_inf = as.numeric(as.character(relative_inf)))
  
  output_R <- as.data.frame.table(output_storage_R) %>% set_colnames(c("population_size", "associativity", "mean_diff","relative_inf", "R")) %>% as.data.frame() %>%
    mutate(population_size = as.numeric(as.character(population_size)), 
           associativity = as.numeric(as.character(associativity)),
           mean_diff = as.numeric(as.character(mean_diff)),
           relative_inf = as.numeric(as.character(relative_inf)),
    )
  
  output <- list()
  output[["output_w"]] <- output_w
  output[["output_R"]] <- output_R
  
  return(output)
}


a <- w_R_by_base_R(w=0.1, distribution_type="gamma", distribution_param_1 =10.56, distribution_param_2 = 1, day_granularity = 0.01, 
                   R_group1=3, mean_diff_vec=seq(0.1,1,0.1), associativity_vec=seq(0,1,0.25), relative_inf_vec=seq(0.2,0.8,0.2), population_sizes_vec=seq(0.2,0.8,0.2), diff_type="rel")



ggplot(a , aes(x=mean_diff, y=relative_inf)) +
  scale_color_gradient(low="yellow", high="red") +
  geom_contour(aes(z=R, color=..level..)) +
  geom_text_contour(aes(z=R), skip=0.1) +
  #facet_wrap(~distribution) +
  theme_bw() +
  scale_y_continuous(labels=mult_format()) +
  labs(title = "test", y = "Relative infectiousness asymptomatics", x = "mean gT difference") +
  theme(legend.position = "none")

output_plot <- function(output, x_var, y_var, z_var=NULL, type, facet_vec, var_filter_vec, var_filter_value_vec,  x_lab, y_lab, title_lab){
  if(length(var_filter_vec) != 0 ) for(i in 1:length(var_filter_vec)) output <- output %>% filter(!!sym(var_filter_vec[i]) == var_filter_value_vec[i])
  
  p <- ggplot(data=output %>% filter(relative_inf != 0, relative_inf != 1), aes(x=!!sym(x_var), y=!!sym(y_var))) +
    theme_bw() + 
    labs(x = x_lab, y=y_lab, title=title_lab) 
  
  if("line" %in% type) p <- p + geom_line()
  if("contour" %in% type) p <- p +  geom_contour(aes(z=!!sym(z_var), color=..level..), breaks=brks) #+ geom_dl(aes(label=..levels..), method="bottom.pieces", stat="contour", breaks=brks)
  if(length(facet_vec)==1) p <- p + facet_wrap(facet_vec[1])
  if(length(facet_vec)==2) p  <- p + facet_grid(row=vars(!!sym(facet_vec[1])), col=vars(!!sym(facet_vec[2])))  
  
  print(plotting_df)
  print(p)
}

output_plot(output=a, x_var="associativity", y_var="relative_inf", facet_vec = c(), z_var="R", type="contour", var_filter_vec=c("population_size","mean_diff"), var_filter_value_vec=c(0.4, 0.7), x_lab="associativity", y_lab="R", title_lab="Inferred R by associativity")
output_plot(output=a, x_var="associativity", y_var="relative_inf", facet_vec = c("population_size"), z_var="R", type="contour", var_filter_vec=c("mean_diff"), var_filter_value_vec=c(0.7), x_lab="associativity", y_lab="R", title_lab="Inferred R by associativity")
output_plot(output=a, x_var="associativity", y_var="relative_inf", facet_vec = c("population_size","mean_diff"), z_var="R", type="contour", var_filter_vec=c(), var_filter_value_vec=c(), x_lab="associativity", y_lab="R", title_lab="Inferred R by associativity")



gT_gamma1 <- gT_generator(distribution_type="gamma", distribution_param_1=5, distribution_param_2=1)$gT_dist
gT_gamma2 <- gT_proposal_trunc(truncation_point = 0.15, distribution_type="gamma", distribution_param_1=5, distribution_param_2=1, normalise=TRUE)
risk_trans2 =c(1,0.15)
gT_gamma3 <- gT_proposal_trunc(truncation_point = 0.75, distribution_type="gamma", distribution_param_1=5, distribution_param_2=1, normalise=TRUE)
risk_trans3 =c(1,0.76)

contact_matrix=matrix(c(0.5,0.5,0.5,0.5), nrow=2, byrow=TRUE)

risk_matrix2 <- relative_risk_mat_generator(risk_trans2, contact_matrix)
risk_matrix3 <- relative_risk_mat_generator(risk_trans3, contact_matrix)

R_solver_nd(w=0.1, gT_list=cbind(gT_gamma1,gT_gamma2),risk_matrix = risk_matrix2, day_granularity = 0.01)$overall_R
R_solver_nd(w=0.1, gT_list=cbind(gT_gamma1,gT_gamma3),risk_matrix = risk_matrix3, day_granularity = 0.01)$overall_R

w_solver_nd(R=1.55, gT_list=cbind(gT_gamma1,gT_gamma2),risk_matrix = risk_matrix2, day_granularity = 0.01)
w_solver_nd(R=1.55, gT_list=cbind(gT_gamma1,gT_gamma3),risk_matrix = risk_matrix3, day_granularity = 0.01)






weighted_average_offspring <- function(R, risk_matrix, contact_matrix){
  total_offspring <- R*colSums(risk_matrix)
  return(total_offspring)
}

R_by_gT_rel <- function(w, distribution_type, distribution_param_1, distribution_param_2, mean_diff_min= 0.5, mean_diff_max = 2, mean_diff_by = 0.1, risk_symp_by = 0.1, day_granularity=0.01){
  gt_group1 <- gT_generator(distribution_type=distribution_type,distribution_param_1=distribution_param_1,distribution_param_2=distribution_param_2)$gT_dist
  
  risk_symp_trials <- seq(0,1,risk_symp_by)
  risk_matrix_set <- array(NA, dim=c(2,2,length(risk_symp_trials)))
  for(i in 1:dim(risk_matrix_set)[3]) risk_matrix_set[,,i] <- relative_risk_mat_generator(c(risk_symp_trials[i],1-risk_symp_trials[i]),contact_matrix)
  
  mean_difference_trials <- seq(mean_diff_min,mean_diff_max,mean_diff_by)
  gt_set <- matrix(NA, nrow=length(gt_group1), ncol=length(mean_difference_trials))
  for(j in 1:ncol(gt_set)) gt_set[,j] = gT_proposal(mean_difference=mean_difference_trials[j],distribution_type=distribution_type, distribution_param_1 = distribution_param_1,distribution_param_2 = distribution_param_2)
  
  output_storage <- matrix(NA, nrow=length(risk_symp_trials), ncol=length(mean_difference_trials))
  
  for(i in 1:length(risk_symp_trials)){
    print(i)
    for(j in 1:length(mean_difference_trials)){
      output_storage[i,j] <- R_solver_nd(w, cbind(gt_group1, gt_set[,j]), risk_matrix_set[,,i], day_granularity)$overall_R
    }
  }
  colnames(output_storage) <- mean_difference_trials
  rownames(output_storage) <- risk_symp_trials
  
  output_storage <- as.data.frame(output_storage) %>%
    mutate(distribution = distribution_type)
  
  return(output_storage)
}


R_gT_A_graph1 <- R_by_gT_A()
R_gT_A_graph2 <- R_by_gT_A(w=0.1, distribution_type="gamma", distribution_param_1=5, distribution_param_2=1, 
                           associativity_vec = c(1), population_sizes_vec = seq(0.01, 0.99, 0.05),
                           gT_diff = "trunc", relative_inf_vec = seq(0.01, 0.99, 0.05), day_granularity=0.01)
R_gT_A_graph3 <- R_by_gT_A(w=0.1, distribution_type="gamma", distribution_param_1=5, distribution_param_2=1, 
                           associativity_vec = seq(0.1,1,0.05), population_sizes_vec = c(0.5),
                           gT_diff = "relative", mean_diff_vec = seq(0.5, 2, 0.2), relative_inf_vec = seq(0, 1, 0.05),
                           day_granularity=0.01)



a <- R_by_gT_A_prepare(R_gT_A_graph1)
b <- R_by_gT_A_prepare(R_gT_A_graph2)
d <- R_by_gT_A_prepare(R_gT_A_graph3)

ggplot(a, aes(x=mean_diff, y=relative_inf)) +
  scale_color_gradient(low="yellow", high="red") +
  geom_contour(aes(z=R, color=..level..)) +
  geom_text_contour(aes(z=R), skip=0.1) +
  #facet_wrap(~distribution) +
  theme_bw() +
  scale_y_continuous(labels=mult_format()) +
  labs(title = title, y = "Relative infectiousness asymptomatics", x = "mean gT difference") +
  theme(legend.position = "none")

ggplot(b, aes(x=relative_inf, y=population_size)) +
  #geom_tile(aes(fill=R)) +
  scale_color_gradient(low="yellow", high="red") +
  geom_contour(aes(z=R, color=..level..)) +
  geom_text_contour(aes(z=R), skip=0.1) +
  #facet_wrap(~distribution) +
  theme_bw() +
  labs(title = title, y = "Population size", x = "mean gT difference") +
  theme(legend.position = "none")

ggplot(d %>% filter(mean_diff == 0.5, relative_inf == 0.25), aes(x=associativity, y=R)) +
  geom_line() +
  theme_bw() 

R_gT_A_prepared_heatmap <- function(output_prepared){
  p <- ggplot(data = output_prepared, aes(x=associativity, y=population_size)) +
    #geom_tile(aes(fill=R)) +
    facet_wrap(mean_diff~peak_height) +
    #scale_color_gradient(low="yellow", high="red") +
    geom_contour(aes(z=R)) +
    geom_text_contour(aes(z=R), skip=0.1) +
    theme_bw() +
    theme(legend.position = "none")
  
  print(p)    
}

R_gT_A_prepared_heatmap(a)

output_plot_prepare <- function(output){
  output_formatted <- as.data.frame(output) %>%
    tibble::rownames_to_column(var="symp_risk") %>%
    tidyr::gather(mean_gt_diff, "R", 2:(ncol(.)-1)) %>%
    mutate(symp_risk = as.numeric(symp_risk)/100, mean_gt_diff = as.numeric(mean_gt_diff), R = as.numeric(R)) 
  
  return(output_formatted)
}

output_plotting_function <- function(output_set, title="R by mean difference in generation time and relative symptomatic risk", x_lab="Mean gT difference", y_lab="Infectiousness symptomatics/infectiousness asymptomatics"){
  output_plot <- do.call(rbind, lapply(output_set, output_plot_prepare))
  
  p <- ggplot(output_plot, aes(x=mean_gt_diff, y=symp_risk)) +
    #geom_tile(aes(fill=R)) +
    scale_color_gradient(low="yellow", high="red") +
    geom_contour(aes(z=R, color=..level..)) +
    geom_text_contour(aes(z=R), skip=0.1) +
    facet_wrap(~distribution) +
    theme_bw() +
    scale_y_continuous(labels=mult_format()) +
    labs(title = title, y = y_lab, x = x_lab) +
    theme(legend.position = "none")
  
  ggsave(plot=p, filename=paste("PhD/Graphs/new_graphs/inferred_R_symp_asymp.png"), width=5, height=5)
  
  print(p)
  
}

R_by_gT_trunc <- function(w, distribution_type, distribution_param_1, distribution_param_2, proportion_isol_min = 0.1, proportion_isol_max = 0.9, proportion_isol_by = 0.1, trunc_min= 0.1, trunc_max = 0.9, trunc_diff_by = 0.1, normalise, day_granularity=0.01){
  gt_no_isol <- gT_generator(distribution_type=distribution_type,distribution_param_1=distribution_param_1,distribution_param_2=distribution_param_2)$gT_dist
  
  proportion_isolating <- seq(proportion_isol_min, proportion_isol_max, proportion_isol_by)
  proportion_matrix_set <- array(NA, dim=c(2,2,length(proportion_isolating)))
  
  truncation_trials <- seq(trunc_min,trunc_max,trunc_diff_by)
  gt_isol <- matrix(NA, nrow=length(gt_no_isol), ncol=length(truncation_trials))
  for(j in 1:ncol(gt_isol)) gt_isol[,j] = gT_proposal_trunc(truncation_point=truncation_trials[j],distribution_type=distribution_type, distribution_param_1 = distribution_param_1,distribution_param_2 = distribution_param_2, normalise=TRUE)
  
  output_storage <- matrix(NA, nrow=length(proportion_isolating), ncol=length(truncation_trials))
  
  for(i in 1:nrow(output_storage)){
    for(j in 1:ncol(output_storage)){
      risk_trans = c(1/(1+truncation_trials[j]),truncation_trials[j]/(1+truncation_trials[j]))
      contact_matrix = matrix(c(1-proportion_isolating[i],1-proportion_isolating[i],
                                proportion_isolating[i],proportion_isolating[i]), nrow=2, byrow=TRUE)
      
      risk_matrix <- relative_risk_mat_generator(risk_trans, contact_matrix)
      
      output_storage[i,j] <- R_solver_nd(w, cbind(gt_no_isol, gt_isol[,j]), risk_matrix=risk_matrix, day_granularity = day_granularity)
      
      cat("\n")
      print("risk_trans = ", quote=F)
      print(risk_trans)
      
      print(c("contact_matrix ="),quote=F)
      print(contact_matrix)
      
      risk_matrix <- relative_risk_mat_generator(risk_trans, contact_matrix)
      
      print(c("risk_matrix ="),quote=F)
      print(risk_matrix)
      
      print(c("R ="),quote=F)
      print(output_storage[i,j])
      
      print(c("implied_risk_matrix ="),quote=F)
      print(output_storage[i,j] * risk_matrix)
      
      cat("\n")
      
    }
  }
  
  colnames(output_storage) <- truncation_trials
  rownames(output_storage) <- proportion_isolating
  
  output_storage <- as.data.frame(output_storage) %>%
    mutate(distribution = distribution_type) %>%
    mutate(proportion_isolating = proportion_isolating)
  
  return(output_storage)
}

output_plot_prepare_trunc <- function(output, distribution_type, fit_type=NULL){
  output_formatted <- as.data.frame(output) %>%
    tidyr::gather(infection_passed_isolation, "R", 2:ncol(.)-2) %>%
    mutate(proportion_isolating = as.numeric(proportion_isolating), 
           infection_passed_isolation = as.numeric(infection_passed_isolation), 
           R = as.numeric(R)) 
  
  return(output_formatted)
}

output_plotting_function_trunc <- function(output_set, title="Inferred R when r=0.1", x_lab="Proportion infectivity passed at isolation", y_lab="Proportion isolating"){
  output_plot <- do.call(rbind, lapply(output_set, output_plot_prepare_trunc))
  
  p <- ggplot(output_plot, aes(x=infection_passed_isolation, y=proportion_isolating)) +
    #geom_tile(aes(fill=R)) +
    geom_contour(aes(z=R, color=..level..)) +
    scale_color_gradient(low="yellow", high="red") +
    geom_text_contour(aes(z=R),) +
    facet_wrap(~distribution) +
    theme_bw() +
    labs(title = title, y = y_lab, x = x_lab) +
    theme(legend.position = "none")
  
  print(p)
  
  return(p)
  
  ggsave(plot=p, filename=paste("PhD/Graphs/new_graphs/inferred_R_by_prop_isol_inf.png"), width=5, height=5)
}

w_by_gT_trunc <- function(R, distribution_type, distribution_param_1, distribution_param_2, proportion_symp_by = 0.1, trunc_min = 0.1, trunc_max = 0.9, trunc_diff_by = 0.1, normalise, day_granularity=0.01){
  gt_no_isol <- gT_generator(distribution_type=distribution_type,distribution_param_1=distribution_param_1,distribution_param_2=distribution_param_2)$gT_dist
  
  proportion_isolating <- seq(0.1, 0.9, proportion_symp_by)
  proportion_matrix_set <- array(NA, dim=c(2,2,length(proportion_isolating)))
  
  truncation_trials <- seq(trunc_min,trunc_max,trunc_diff_by)
  gt_isol <- matrix(NA, nrow=length(gt_no_isol), ncol=length(truncation_trials))
  for(j in 1:ncol(gt_isol)) gt_isol[,j] = gT_proposal_trunc(truncation_point=truncation_trials[j],distribution_type=distribution_type, distribution_param_1 = distribution_param_1,distribution_param_2 = distribution_param_2, normalise=FALSE)
  
  output_storage <- matrix(NA, nrow=length(proportion_isolating), ncol=length(truncation_trials))
  
  for(i in 1:nrow(output_storage)){
    print(c("i = ",i), quote=F)
    for(j in 1:ncol(output_storage)){
      contact_matrix = matrix(c(1-proportion_isolating[i],1-proportion_isolating[i],
                                proportion_isolating[i],proportion_isolating[i]), nrow=2, byrow=TRUE)
      
      risk_trans <- c(1, truncation_trials[j])
      
      cat("\n")
      print("risk_trans = ", quote=F)
      print(risk_trans)
      
      print(c("contact_matrix ="),quote=F)
      print(contact_matrix)
      
      risk_matrix <- relative_risk_mat_generator(risk_trans, contact_matrix)
      
      print(c("risk_matrix ="),quote=F)
      print(risk_matrix)
      
      det(risk_matrix)
      output_storage[i,j] <- w_solver_nd(R=R, gT_list = cbind(gt_no_isol, gt_isol[,j]), risk_matrix=risk_matrix, day_granularity = 0.01)
      cat("\n")
      print(c("w = ",output_storage[i,j]), quote=F)
      cat("\n")
      cat("\n")
    }
  }
  
  colnames(output_storage) <- truncation_trials
  rownames(output_storage) <- proportion_isolating
  
  output_storage <- as.data.frame(output_storage) %>%
    mutate(distribution = distribution_type) %>%
    mutate(proportion_isolating = proportion_isolating)
  
  return(output_storage)
  
}

w_output_plot_prepare_trunc <- function(output, distribution_type, fit_type=NULL){
  output_formatted <- as.data.frame(output) %>%
    tidyr::gather(infection_passed_isolation, "w", 2:ncol(.)-2) %>%
    mutate(proportion_isolating = as.numeric(proportion_isolating), infection_passed_isolation = as.numeric(infection_passed_isolation), w = as.numeric(w)) 
  
  return(output_formatted)
}

w_output_plotting_function_trunc <- function(output_set, title="Growth rate by R of next generation matrix", x_lab="Proportion infectivity passed at isolation", y_lab="Proportion isolating"){
  output_plot <- do.call(rbind, lapply(output_set, w_output_plot_prepare_trunc))
  
  p <- ggplot(output_plot, aes(x=infection_passed_isolation, y=proportion_isolating)) +
    geom_tile(aes(fill=w)) +
    scale_fill_gradient(low="yellow", high="red") +
    geom_contour(aes(z=w),color="black") +
    geom_text_contour(aes(z=w),) +
    facet_wrap(~distribution) +
    theme_bw() +
    labs(title = title, y = y_lab, x = x_lab)
  
  print(p)
}

output_1g <- R_by_gT_rel(w=0.1, distribution_type="gamma", distribution_param_1=5, distribution_param_2=1, mean_diff_by = 0.01, risk_symp_by = 0.01, day_granularity=0.01)
output_3gR <- R_by_gT_trunc(w=0.1, distribution_type="gamma", distribution_param_1=5, distribution_param_2=1, trunc_min= 0.01, trunc_max = 0.99, trunc_diff_by = 0.1, proportion_isol_min = 0.1 , proportion_isol_max = 0.99, proportion_isol_by = 0.01, normalise=FALSE, day_granularity=0.01)
output_3gw <- w_by_gT_trunc(R=3, distribution_type="gamma", distribution_param_1=5, distribution_param_2=1, trunc_min= 0.01, trunc_max = 0.99, trunc_diff_by = 0.01, normalise=FALSE)


output_plotting_function(output_set=list(output_1g), title="R by relative generation time and infection risk", x_lab="Mean gT difference", y_lab="Infectiousness symptomatics/infectiousness asymptomatics")
output_3gr_plot <- output_plotting_function_trunc(output_set=list(output_3gR), x_lab="Proportion infectivity passed at isolation", y_lab="Proportion isolating")

ggsave(plot = grid.arrange(inf_0.15_50, output_3gr_plot, inf_0.75_50, nrow=1),filename=paste("PhD/Graphs/new_graphs/inferred_R_by_prop_isol_inf_wt_gT_dist.png"), width=15, height=5)

overall_infectious_curve_plotter <- function(prop_isol, truncation_point, distribution_type, distribution_param_1, distribution_param_2){
  gT_no_isol <- gT_generator(distribution_type, distribution_param_1, distribution_param_2)$gT_dist
  gT_isol <- gT_proposal_trunc(truncation_point, distribution_type, distribution_param_1, distribution_param_2, normalise=FALSE)
  
  risk_matrix <- relative_risk_mat_generator(risk_trans = c(1,truncation_point), contact_matrix=matrix(c(1-prop_isol, 1-prop_isol, 
                                                                                                         prop_isol, prop_isol), nrow=2, byrow=T))
  gT_matrix <- cbind(gT_no_isol, gT_isol)
  
  output <- data.frame(day = seq(0,100,0.01),
                       inf_no_isol = rowSums(t(colSums(risk_matrix) * t(gT_matrix)))) %>%
    mutate(inf_no_isol = inf_no_isol/sum(inf_no_isol))
  
  p <- ggplot(data=output, aes(x=day, y = inf_no_isol*100)) +
    geom_line() + 
    xlim(c(0,25)) +
    theme_bw() +
    labs(x = "day", y = "probability density", title = paste("Weighted infectiousness distribution" , sep=""))
  
  print(p)
  
  return(p)
  ggsave(plot=p, filename=paste("PhD/Graphs/new_graphs/",metric,"weighted_gT_distribution.png"), width=10, height=5)
}

inf_0.15_50 <- overall_infectious_curve_plotter(prop_isol=0.5, truncation_point = 0.15, distribution_type="gamma", distribution_param_1=5, distribution_param_2 = 1)
inf_0.75_50 <-overall_infectious_curve_plotter(prop_isol=0.5, truncation_point = 0.75, distribution_type="gamma", distribution_param_1=5, distribution_param_2 = 1)

w_solver_isolation <- function(R_no_isol, truncation_point, proportion_isolating, distribution_type, distribution_param_1, distribution_param_2, day_granularity=0.01){
  gT_no_isol <- gT_generator(distribution_type=distribution_type,distribution_param_1=distribution_param_1,distribution_param_2=distribution_param_2)$gT_dist
  gT_isol <- gT_proposal_trunc(truncation_point, distribution_type, distribution_param_1, distribution_param_2)
  
  R_isol <- R_no_isol * truncation_point
  
  contact_matrix <- matrix(c(1-proportion_isolating, 1-proportion_isolating,
                             proportion_isolating,   proportion_isolating),nrow=2,byrow=T)
  
  
  R_matrix <- t(c(R_no_isol, R_isol) * t(contact_matrix))
  R <- eigen(R_matrix)$values[1]
  
  risk_matrix <- relative_risk_mat_generator(risk_trans = c(1, truncation_point), contact_matrix = contact_matrix)
  
  w <- w_solver_nd(R=R, gT_list=cbind(gT_no_isol,gT_isol), risk_matrix = risk_matrix, day_granularity = day_granularity)
  output <- list()
  output[["R"]] <- R
  output[["w"]] <- w
  return(output)
  
}

w_solver_isolation(R_no_isol = 1.76, truncation_point = 0.75, proportion_isolating = 0.5, distribution_type = "gamma", distribution_param_1 =10.56, distribution_param_2 = 1)
w_solver_isolation(R_no_isol = 2.69, truncation_point = 0.15, proportion_isolating = 0.5, distribution_type = "gamma", distribution_param_1 =10.56, distribution_param_2 = 1)


w_by_prop_isolating <- function(R_no_isol, truncation_point_min=0.1, truncation_point_max=0.9, truncation_point_by=0.1, proportion_isolating_min=0.1, proportion_isolating_max=0.9, proportion_isolating_by=0.1, distribution_type, distribution_param_1, distribution_param_2){
  truncation_trials <- seq(truncation_point_min, truncation_point_max, truncation_point_by)
  proportion_trials <- seq(proportion_isolating_min, proportion_isolating_max, proportion_isolating_by)
  
  output_w <- matrix(NA, nrow=length(proportion_trials), ncol=length(truncation_trials))
  output_R <- matrix(NA, nrow=length(proportion_trials), ncol=length(truncation_trials))
  
  for(i in 1:nrow(output_w)){
    print(truncation_trials[i])
    for(j in 1:ncol(output_w)){
      results <- w_solver_isolation(R_no_isol = R_no_isol, truncation_point = truncation_trials[j], proportion_isolating = proportion_trials[i], distribution_type=distribution_type, distribution_param_1 = distribution_param_1, distribution_param_2 = distribution_param_2)
      output_w[i,j] <- results$w
      output_R[i,j] <- results$R
    }
  }
  
  output <- list()
  
  output_w <- as.data.frame(output_w) %>%
    set_colnames(truncation_trials) %>%
    mutate(proportion_isolating = proportion_trials, 
           distribution = distribution_type)
  
  output_R <- as.data.frame(output_R) %>%
    set_colnames(truncation_trials) %>%
    mutate(proportion_isolating = proportion_trials, 
           distribution = distribution_type)
  
  output[["w"]] <- output_w
  output[["R"]] <- output_R
  
  return(output)
}

output_5 <- w_by_prop_isolating(R_no_isol=3, truncation_point_min=0.01, truncation_point_max=0.99, truncation_point_by=0.01, proportion_isolating_min=0.01, proportion_isolating_max=0.99, proportion_isolating_by=0.01, distribution_type="gamma", distribution_param_1=5, distribution_param_2=1)

w_by_prop_output_prepare <- function(output){
  output_w_formatted <- as.data.frame(output$w) %>%
    tidyr::gather(infection_passed_isolation, "w", 1:(ncol(.)-2)) %>%
    mutate(proportion_isolating = as.numeric(proportion_isolating), infection_passed_isolation = as.numeric(infection_passed_isolation), w = as.numeric(w)) 
  
  output_R_formatted <- as.data.frame(output$R) %>%
    tidyr::gather(infection_passed_isolation, "R", 1:(ncol(.)-2)) %>%
    mutate(proportion_isolating = as.numeric(proportion_isolating), infection_passed_isolation = as.numeric(infection_passed_isolation), R = as.numeric(R)) 
  
  output_formatted <- left_join(output_w_formatted, output_R_formatted, by=c("proportion_isolating", "distribution", "infection_passed_isolation")) %>%
    tidyr::gather(metric, "value", (ncol(.)-1):ncol(.))
  
  return(output_formatted)
}

output_5_prepared <- w_by_prop_output_prepare(output_5)

w_by_prop_trunc_plotting_function <- function(output_prepared, title="", x_lab="", y_lab=""){
  p <- ggplot(data=output_prepared %>% filter(metric=="w"), aes(x=infection_passed_isolation, y=proportion_isolating, color=value)) +
    #geom_tile(aes()) +
    #scale_fill_gradient(low="yellow", high="red") +
    geom_contour(aes(z=value, color=..level..)) +
    scale_color_gradient(low="green", high="red") +
    geom_text_contour(aes(z=value)) +
    theme_bw() +
    labs(title = "Epidemic growth rate, r", x = "Proportion infectivity passed at isolation", y = "Proportion isolating") +
    theme(axis.title = element_blank(), 
          legend.position = "none") + 
    xlim(0,1) + ylim(0,1)
  
  q <- ggplot(data=output_prepared %>% filter(metric=="R"), aes(x=infection_passed_isolation, y=proportion_isolating, color=value)) +
    #geom_tile(aes()) +
    #scale_fill_gradient(low="yellow", high="red") +
    geom_contour(aes(z=value,color=..level..)) +
    scale_color_gradient(low="green", high="red") +
    geom_text_contour(aes(z=value)) +
    theme_bw() +
    labs(title = "Reproduction number, R", x = "Proportion infectivity passed at isolation") +
    theme(axis.title=element_blank(), 
          legend.position = "none") + 
    xlim(0,1) + ylim(0,1)
  
  r <- grid.arrange(p, q, nrow=1, bottom="Proportion infectivity passed at isolation", left="Proportion isolating")
  
  ggsave(plot=r, filename=paste("PhD/Graphs/new_graphs/",metric,"with_trunc_prop.png"), width=10, height=5)
  
  return(r)
}

w_by_prop_trunc_plotting_function(output_prepared = output_5_prepared)
w_by_prop_trunc_plotting_function(output_prepared = output_5_prepared, metric_value = "R", title = "R by isolation adherence and isolation time", x_lab = "Infection passed at time of isolation", y_lab = "Proportion isolating")

w_by_prop_plotting_function <- function(output_prepared, trunc_vec, metric_value = "w", title="Epidemic growth rate by isolation point and proportion isolating", x_lab="Proportion isolating", y_lab="Epidemic growth rate"){
  output_prepared1 <- output_prepared %>%
    dplyr::filter(infection_passed_isolation %in% trunc_vec, metric == metric_value) %>%
    mutate(infection_passed_isolation = factor(infection_passed_isolation, levels=rev(trunc_vec)))
  
  p <- ggplot(data = output_prepared1, aes(x=proportion_isolating, y=value)) +
    geom_line(aes(colour = infection_passed_isolation)) +
    theme_bw()+
    #scale_color_manual(guide=guide_legend(reverse=T)) +
    labs(title = title, y = y_lab, x = x_lab, color="infectivity passed\nat isolation")
  
  return(p)
}

w1 <- w_by_prop_plotting_function(output_prepared = output_5_prepared, trunc_vec=round(seq(0.1,0.9,0.2),2), metric_value = "w", title="Epidemic growth rate by isolation point and proportion isolating", x_lab="Proportion isolating", y_lab="Epidemic growth rate")

w_by_trunc_plotting_function <- function(output_prepared, prop_vec, metric_value = "w", title="Epidemic growth rate by isolation point and proportion isolating", x_lab="Infectivity passed at isolation", y_lab="Epidemic growth rate"){
  output_prepared2 <- output_prepared %>% mutate(proportion_isolating = as.factor(proportion_isolating))
  
  p <- ggplot(data = output_prepared2 %>% filter(proportion_isolating %in% prop_vec, metric == metric_value), aes(x=infection_passed_isolation, y=value)) +
    geom_line(aes(colour = proportion_isolating)) +
    theme_bw()+
    labs(title = title, y = y_lab, x = x_lab, color="proportion\nisolating")
  
  return(p)
}

w2 <- w_by_trunc_plotting_function(output_prepared = output_5_prepared, prop_vec=seq(0.1,0.9,0.2), metric_value = "w", title="Epidemic growth rate by isolation point and proportion isolating", x_lab="Infection passed isolation", y_lab="Epidemic growth rate")


grid.arrange(w1, w2, ncol=2)








weighted_average_offspring(R = 3, risk_matrix = risk_matrix_0.01, contact_matrix = matrix(c(0.3,0.3,0.7,0.7), nrow=2, byrow=T))

gT_gamma_1 <- gT_generator(distribution_type="gamma", distribution_param_1 =10.56,distribution_param_2 = 1)$gT_dist

gT_gamma_0.15 <- gT_proposal_trunc(truncation_point=0.15, distribution_type="gamma", distribution_param_1=5, distribution_param_2=1, plot=TRUE, normalise = TRUE)
gT_gamma_0.76 <- gT_proposal_trunc(truncation_point=0.76, distribution_type="gamma", distribution_param_1=5, distribution_param_2=1, plot=TRUE, normalise = TRUE)

risk_matrix_0.15 <- relative_risk_mat_generator(risk_trans = c(1,0.15), contact_matrix = matrix(c(0.5,0.5,0.5,0.5), nrow=2, byrow=T))
risk_matrix_0.76 <- relative_risk_mat_generator(risk_trans = c(1,0.76), contact_matrix = matrix(c(0.5,0.5,0.5,0.5), nrow=2, byrow=T))

w_solver_nd(R=1.55, gT_list = cbind(gT_gamma_1, gT_gamma_0.15), risk_matrix=risk_matrix_0.15, day_granularity = 0.01)
w_solver_nd(R=1.55, gT_list = cbind(gT_gamma_1, gT_gamma_0.76), risk_matrix=risk_matrix_0.76, day_granularity = 0.01)

R_solver_nd(w=0.1, gT_list = cbind(gT_gamma_1, gT_gamma_0.76), risk_matrix=risk_matrix_0.76, day_granularity = 0.01)
R_solver_nd(w=0.1, gT_list = cbind(gT_gamma_1, gT_gamma_0.15), risk_matrix=risk_matrix_0.15, day_granularity = 0.01)

a <- numerical_outbreak_simulation(R=1.55, gT_list=cbind(gT_gamma_1, gT_gamma_0.15), risk_matrix = risk_matrix_0.15)
b <- numerical_outbreak_simulation(R=1.55, gT_list=cbind(gT_gamma_1, gT_gamma_0.76), risk_matrix = risk_matrix_0.76)

library(plotly)
library(reshape2)
df <- melt(volcano)

p <- ggplot(df, aes(Var1, Var2, colour=stat(level))) +
  geom_contour(aes(z=value)) 




pop_isol_optim_check <- weighted_gt_generator(distribution_type = "gamma", distribution_param_1 = 10.56, distribution_param_2 = 1.85, mean_difference = 0.5, truncation_point = 0.3, risk_trans = c(1, 0.3), pop_vec = c(1-0.63* 0.75, 0.63* 0.75), isol=TRUE, asymp=FALSE, plot=TRUE, aggregate_first = FALSE)
lines(gt_isol_no_isol_optim, col="red")


gt_no_isol <- gT_generator("gamma", shape, rate)$gT_dist

gt_isol_optim <- gT_proposal_trunc(0.3,"gamma", shape, rate)
gt_isol_pesim <- gT_proposal_trunc(0.7,"gamma", shape, rate)

gt_no_isol_daily <- vector(length = 20)
gt_isol_optim_daily <- vector(length = 20)

for(i in 1:length(gt_no_isol_daily)) gt_no_isol_daily[i] = sum(gt_no_isol[((i-1)*100):(i*100)])
for(i in 1:length(gt_isol_optim_daily)) gt_isol_optim_daily[i] = sum(gt_no_isol[((i-1)*100):(i*100)])



gt_asymp_optim <- gT_proposal(0.5, "gamma", shape, rate)
gt_asymp_pesim <- gT_proposal(2, "gamma", shape, rate)

risk_trans_isol_optim = c(1,0.3,1)
risk_trans_isol_pesim= c(1,0.7,1)

risk_trans_asymp_optim = c(1,1,0.5)
risk_trans_asymp_pesim= c(1,1,2)

risk_trans_isol_asymp_optim = c(1,0.3,0.5)
risk_trans_isol_asymp_pesim= c(1,0.7,2)

pop_isol_optim <- 0.63* 0.75
pop_no_isol_optim <- 1 - pop_isol_optim

pop_isol_pesim <- 0.63 * 0.25
pop_no_isol_pesim <- 1 - pop_isol_pesim

pop_symp <- 0.63
pop_asymp <- 0.37


gt_isol_no_isol_optim <- (risk_trans_isol_optim[1]*gt_no_isol*pop_no_isol_optim + risk_trans_isol_optim[2]*gt_isol_optim*pop_isol_optim)/sum(risk_trans_isol_optim[1]*gt_no_isol*pop_no_isol_optim + risk_trans_isol_optim[2]*gt_isol_optim*pop_isol_optim)
gt_isol_no_isol_pesim <- (risk_trans_isol_pesim[1]*gt_no_isol*pop_no_isol_pesim + risk_trans_isol_pesim[2]*gt_isol_pesim*pop_isol_pesim)/sum(risk_trans_isol_pesim[1]*gt_no_isol*pop_no_isol_pesim + risk_trans_isol_pesim[2]*gt_isol_pesim*pop_isol_pesim)

gt_symp_asymp_optim <- (risk_trans_asymp_optim[1]*gt_no_isol*pop_symp + risk_trans_asymp_optim[3]*gt_asymp_optim*pop_asymp)/sum(risk_trans_asymp_optim[1]*gt_no_isol*pop_symp + risk_trans_asymp_optim[3]*gt_asymp_optim*pop_asymp)
gt_symp_asymp_pesim <- (risk_trans_asymp_pesim[1]*gt_no_isol*pop_symp + risk_trans_asymp_pesim[3]*gt_asymp_pesim*pop_asymp)/sum(risk_trans_asymp_pesim[1]*gt_no_isol*pop_symp + risk_trans_asymp_pesim[3]*gt_asymp_pesim*pop_asymp)

gt_isol_asymp_optim <- (risk_trans_isol_asymp_optim[1]*gt_no_isol*(1-pop_asymp-pop_isol_optim) + risk_trans_isol_asymp_optim[2]*gt_isol_optim*pop_isol_optim + risk_trans_isol_asymp_optim[3]*gt_asymp_optim*pop_asymp)/sum(risk_trans_isol_asymp_optim[1]*gt_no_isol*(1-pop_asymp-pop_isol_optim) + risk_trans_isol_asymp_optim[2]*gt_isol_optim*pop_isol_optim + risk_trans_isol_asymp_optim[3]*gt_asymp_optim*pop_asymp)
gt_isol_asymp_pesim <- (risk_trans_isol_asymp_pesim[1]*gt_no_isol*(1-pop_asymp-pop_isol_pesim) + risk_trans_isol_asymp_pesim[2]*gt_isol_pesim*pop_isol_pesim + risk_trans_isol_asymp_pesim[3]*gt_asymp_pesim*pop_asymp)/sum(risk_trans_isol_asymp_pesim[1]*gt_no_isol*(1-pop_asymp-pop_isol_pesim) + risk_trans_isol_asymp_pesim[2]*gt_isol_pesim*pop_isol_pesim + risk_trans_isol_asymp_pesim[3]*gt_asymp_pesim*pop_asymp)

gt_no_isol_daily <- vector(length = 20)
gt_isol_no_isol_optim_daily <- vector(length = 20)
gt_isol_no_isol_pesim_daily <- vector(length = 20)
gt_symp_asymp_optim_daily <- vector(length = 20)
gt_symp_asymp_pesim_daily <- vector(length = 20)
gt_isol_asymp_optim_daily <- vector(length = 20)
gt_isol_asymp_pesim_daily <- vector(length = 20)


for(i in (1:length(gt_no_isol_daily))) gt_no_isol_daily[i] = sum(gt_no_isol[((i-1)*100+1):(i*100)])
for(i in (1:length(gt_isol_no_isol_optim_daily))) gt_isol_no_isol_optim_daily[i] = sum(gt_isol_no_isol_optim[((i-1)*100+1):(i*100)])
for(i in (1:length(gt_isol_no_isol_pesim_daily))) gt_isol_no_isol_pesim_daily[i] = sum(gt_isol_no_isol_pesim[((i-1)*100+1):(i*100)])
for(i in (1:length(gt_symp_asymp_optim_daily))) gt_symp_asymp_optim_daily[i] = sum(gt_symp_asymp_optim[((i-1)*100+1):(i*100)])
for(i in (1:length(gt_symp_asymp_pesim_daily))) gt_symp_asymp_pesim_daily[i] = sum(gt_symp_asymp_pesim[((i-1)*100+1):(i*100)])
for(i in (1:length(gt_isol_asymp_optim_daily))) gt_isol_asymp_optim_daily[i] = sum(gt_isol_asymp_optim[((i-1)*100+1):(i*100)])
for(i in (1:length(gt_isol_asymp_pesim_daily))) gt_isol_asymp_pesim_daily[i] = sum(gt_isol_asymp_pesim[((i-1)*100+1):(i*100)])

gt_isol_no_isol_optim2 <- (risk_trans_isol_optim[1]*gt_no_isol_daily*pop_no_isol_optim + risk_trans_isol_optim[2]*gt_isol_optim_daily*pop_isol_optim)/sum(risk_trans_isol_optim[1]*gt_no_isol_daily*pop_no_isol_optim + risk_trans_isol_optim[2]*gt_isol_optim_daily*pop_isol_optim)



plot(gt_no_isol_daily, type="l", ylim=c(0,0.4))
lines(gt_isol_no_isol_optim_daily, type="l", col="red")
lines(gt_isol_no_isol_optim2, type="l")

lines(gt_isol_no_isol_pesim_daily, type="l", col="green")
lines(gt_symp_asymp_optim_daily, type="l", col="red")
lines(gt_symp_asymp_pesim_daily, type="l", col="red")
lines(gt_isol_asymp_optim_daily, type="l", col="red")
lines(gt_isol_asymp_pesim_daily, type="l", col="red")
