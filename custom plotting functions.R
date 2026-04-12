
seroprevalence.fit.custom<- function(FOIfit,                   
                              individual_samples = 0,
                              age_class = 10,
                              YLIM=1,
                              colorribbon,
                              colorline,
                              ...){
  
  plots  <- NULL
  data= FOIfit$data
  
  chains <- rstan::extract(FOIfit$fit)
  se = FOIfit$model$se 
  sp = FOIfit$model$sp 
  A <- FOIfit$data$A
  latest_sampling_year <- max(FOIfit$data$sampling_year)
  years <- seq(1,A)
  
  index.plot=0
  unique.categories = data$unique.categories
  sorted.year = sort.int(unique(FOIfit$data$sampling_year),index.return = TRUE)
  Y=0
  for(sampling_year in sorted.year$x ){
    Y=Y+1
    for(cat in unique.categories){ 
      
      if(length(unique.categories)==1){
        #  title =  paste0("Sampling year: ", sampling_year)
        title =  ""
      }
      if(length(unique.categories)>1){
        # title= paste0('Category: ',cat," Sampling year: ", sampling_year)
        title= paste0('Category: ',cat)
      }
      
      index.plot=index.plot+1
      age_group = data$age_group[which(data$sampling_year ==  sampling_year)][1]
      w = which(data$sampling_year ==  sampling_year & data$category==cat, arr.ind = TRUE)[,1]
      subdat = subset(data,sub = w)
      
      # compute the proportion of seropositive
      P=chains$P[,,sorted.year$ix[Y], 1]
      d = data$categoryindex[w]
      p1=proportions.index(d)
      
      M=dim(chains$P)[1] 
      Pinf=matrix(0, nrow = M, ncol=FOIfit$data$A)
      
      # infection probability weighted on the categories      
      for(i in 1:length(p1$index)){
        Pinf =  Pinf+  p1$prop[i]*( se-(se+sp-1)*chains$P[,,sorted.year$ix[Y],p1$index[i]] ) 
      }
      
      par_out <- apply(Pinf, 2, function(x)c(mean(x), quantile(x, probs=c(0.025, 0.975))))
      par_out[par_out>YLIM]= YLIM # set to the upper limits for plotting
      
      # X axis
      years.plotted =  seq(latest_sampling_year-sampling_year+1, dim(chains$P)[2])
      years.plotted.normal= years.plotted-min(years.plotted)+1
      meanFit <- data.frame(x = years.plotted.normal, y = par_out[1,years.plotted ])
      
      # create the envelope
      xpoly <- (c(years.plotted.normal, rev(years.plotted.normal)))
      ypoly <-  c(par_out[3,years.plotted ], rev(par_out[2,years.plotted ]))
      DataEnvelope = data.frame(x = xpoly, y = ypoly)
      
      # histogram  of data
      histdata <- sero.age.groups.kc(dat = subdat,age_class = age_class,YLIM=YLIM) 
      histdata$new_labels_text <- as.character(histdata$new_labels)
      
      last_row_index <- tail(which(complete.cases(histdata)), 1)
      if(length(last_row_index)>0){
        histdata = histdata[1:last_row_index,]
      }
      
      #max.age= max(histdata$age)
      max.age= max(histdata$true_mean_age)
      
      DataEnvelope =  subset(DataEnvelope, x<=max.age)
      meanFit = subset(meanFit, x<=max.age )
      
      # plot the mean and 95% credible interval of the seroprevalence
      p <- ggplot2::ggplot() + 
        ggplot2::geom_polygon(data=DataEnvelope, ggplot2::aes(x, y), fill=colorribbon) + 
        ggplot2::geom_line(data = meanFit, ggplot2::aes(x = x, y = y), linewidth = 1, color =colorline)
      
      # if plot individual runs of the chain
      if(individual_samples>0){
        Index_samples <- sample(nrow(Pinf), individual_samples)
        for (i in Index_samples){
          ind_foi <-  data.frame(x = years.plotted.normal,y = Pinf[i, years.plotted])
          p <- p + ggplot2::geom_line(data = ind_foi, ggplot2::aes(x = x, y = y), linewidth = 0.8, colour = "#bbbbbb", alpha = 0.6)
        }
      }
      
      p <- p  +
        scale_x_continuous(breaks=histdata$new_breaks,
                           labels=histdata$new_labels_text, 
                           expand = c(0,0),
                           limits = c(0,max.age+0.2))+
        scale_y_continuous(expand = c(0,0), limits = c(0,YLIM)) +
        theme_minimal() +
        geom_point(data = histdata, aes(x=true_mean_age, y=mean))  + # Changed to true_mean_age from age
        geom_segment(data=histdata, aes(x=true_mean_age,y=lower, xend= true_mean_age, yend=upper))+ # Changed to true_mean_age from age
        ggplot2::xlab("Age (years)") + 
        ggplot2::ylab("Seropositivity") + 
        theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.85, size = 10), # 11
              axis.text.y = element_text(size=10), # 11
              text=element_text(size=11)) + # 13
        #ylim(0,YLIM)+
        ggtitle(title)
      
      
      plots[[index.plot]] <- p 
      plots[[index.plot]]$category <- cat 
      plots[[index.plot]]$year <- sampling_year 
      
    }
  }
  return(plots)
  
}

proportions.index <- function(d){
  b.x=c()
  b.y=c()
  ii=0
  for(i in unique(d)){
    ii=ii+1
    b.x[ii]  = i  
    b.y[ii]  = sum(d==i)/length(d)
  }
  return(list(index=b.x,prop = b.y ))
  
}

# get the seroprevalence (mean and 95%CI) for each age group
sero.age.groups.kc <- function(dat,age_class,YLIM){
  
  #  age_categories <- seq(from = 0, to = min(dat$A, max(dat$age)), by = age_class)
  
  if(age_class<= max(dat$age_at_sampling)){
    age_categories <- seq(from = 0, to =  max(dat$age_at_sampling), by = age_class)
  }else{
    age_categories <- seq(from = 0, to =  max(dat$age_at_sampling), by = max(dat$age_at_sampling))
  }
  
  age_bin <- sapply(dat$age, function(x) tail(which(x-age_categories >= 0), 1L)) # find the closest element
  S <- as.integer(as.logical(dat$Y)) 
  S1 <- sapply(1:length(age_categories), function(x) length(which(age_bin==x)) )
  S2 <- sapply(1:length(age_categories), function(x) sum(S[which(age_bin==x)] ))
  C <- (rbind((age_categories[1:length(age_categories)-1]), (age_categories[2:length(age_categories)]-1)))
  
  df = data.frame(x=age_categories,y=S2/S1)
  
  G=matrix(NA,nrow =  dim(df)[1], ncol=3)
  
  for(j in seq(1,length(S1))){
    if(S1[j]>1){
      B= binom::binom.confint(x=S2[j],n = S1[j],methods = "exact")
      G[j,1]=B$lower
      G[j,2]=B$upper
      G[j,3]=B$mean
    }
  }
  
  G[which(G >YLIM)] =YLIM
  mean_age =  c( (age_categories[1:length(age_categories)-1] +age_categories[2:length(age_categories)])/2, age_categories[length(age_categories)] ) 
  
  C <- (rbind((age_categories[1:length(age_categories)-1]), (age_categories[2:length(age_categories)]-1)))
  
  if(sum(C[1, ] - C[2, ]) == 0 ){ # means that the age categories are each 1 year long
    histo_label <- append(format(C[1, ]), paste("≥", tail(age_categories, n = 1), sep = ""))
  } else{
    histo_label <- append(apply(format(C), 2, paste, collapse = "-"), paste("≥", tail(age_categories, n = 1), sep = ""))
  }
  
    new_histo_label <- c(age_categories[-length(age_categories)], paste("≥", tail(age_categories, n = 1), sep = ""))
  
  
  # Add a column for average age within each age class
  dff <- data.frame(age = dat$age_at_sampling, age_cat = age_bin)
  
  dff <- dff |>
    group_by(age_cat) |>
    summarise(true_mean_age = round(mean(age), digits = 0))
  
  histdata <- data.frame(age = mean_age,
                         mean = G[,3],
                         lower = G[,1],
                         upper = G[, 2],
                         labels = factor(histo_label, levels=histo_label),
                         age_cat = 1:length(mean_age))
  
  histdata <- histdata |>
    left_join(dff)
  
  histdata$new_breaks <- age_categories
  histdata$new_labels <- factor(new_histo_label, levels=new_histo_label)
  
  return(histdata)
  
  
}


