##############################################
#
# asw.cluster.collection_2.10.1.R
#
##############################################
#
# Version:   2.10.1
# Date:      29.05.2022
# Author:    Andreas Glenz (a.glenz@psychologie.uzh.ch)
#
# Purpose:
# Collection of functions for determining faultlines for groups with n people and up to n subgroups.
# Group data has to be prepared as data-frames with numerical (or dummy-coded categrial) atributes in
# columns and individual data records in rows.
#
# Version 1.101:
# Some minor bugfixes 
#    Columns in the $long format are now convertet to "numeric" values instead of factors.
#    Van Knippenbergs measure can now be calculated for datasets with only two attributes
# 
# Version 1.102:
# Changed asw_cluster.agglomerative function to return a fl-value of zero instead of "NA" for completely homogeneous teams
#
# Version 1.103i:
# Individual-Level Faultlines
#
# Version 1.104ip:
# Parallel processing option. Specify "cores=<number of parallel calculation threads>". Applies to datasets that include multiple groups, separated by the variable 
# specified by the "group.par"-parameter. Output is omitted in parallel mode. Note that R may be "freezed" during parallel calculations.
# 
# Version 1.105ip:
# FL-level (team/individual) as new option in silhouette width function
# 
# Version 1.2:
# New silhouette.with algorithm which takes responds to subgroup heterogenity
# New function "s.outliers" to detect outlier-members
# New function ingroups.ind which detects individual member's ingroups
# 
# Version 1.21:
# Extended output of summary()$long for individual Faultlines 
#
# Version 1.23: 
# Modified asw.cluster routine for increased efficiency
#
# Version 1.3:
# Added option (usesghomo) to calculate silhouette widths that are sensitive to subgroup homogeneity
#
# Version 1.4:
# Added option (by.attr) to calculate silhouette widths separately for each attribute and apply weighting factors afterwards.
# This enables non-variance attributes to reduce faultline-strengths.
# Added "gower"-Metric, which uses Manhattan Distances instead of Euclidean Distances
#
# Version 1.5:
# Modified silhouette.width calculation for by.attr option to correctly use attribute-relevance (weights)
# Added by.attr calculation for individual mode 
#
# Version 1.51:
# Bug fixes in individual mode with by.attr=T
#
# Version 1.52:
# when determining individual faultlines, procedure does not stop at local maximum anymore, but goes through entire team
#
# Version 1.53:
# "Demo"-Mode for particular methods (ASW,...) based on version 1.52 to demonstrate calculation process
#
# Version 1.54:
# by.attr-mode implemented for individual-level faultlines
#
# Version 1.55:
# individual faultline calculation: if "by.attr" option is set to TRUE, attributes are standardized by their range
# and Euclidean distances between individuals are replaced by the weighted average of singe-attribute distances.#
#
# Version 1.56:
# Small change in ingroup.ind procedure that causes group members to be included in the focal member's ingroup during
# ingroup growing phase, even if they don't cause an immediate increase of the focal member's individual silhouette width by themselves.
#
# Version 1.57:
# outliers calculated based on percentile of distances between group members (in this version only - discontinued)
# 
# Version 1.58:
# "All-zero"-dummyvariables removed from faultline calculation
#
# Version 1.581:
# Minor bug fix in faultlines.default
#
# Version 2.0:
# Published Version
#
# Version 2.01:
# Bug Fixes in fuctions fau.dist.mahal, faultlines.calc and fau.numerator.mahal (Thanks to Marian Dragt for correcting these
#
# Version 2.02:
# Minor bug fix in code
#
# Version 2.03:
# Load package parallels upon call of faultlines() function
#
# Version 2.10:
# Attribute weight adjusted for nominal attributes, when by.attr=TRUE
#
# Version 2.10.1:
# Eliminate dependency on package "QuantPsyc"
#
# Version 2.10.2:
# Minor bug fix in code
##############################################


# Load required Libraries
require(nFactors)
require(psych)
require(flexmix)
require(nnet)
require(parallel)



s.outliers <- function(d.w,k,alpha,df)
{
	
	n <- nrow(d.w)
	
	h <- vector(mode="numeric",length=n)
	for (i in 1:n)
	{
		d.sorted <- (cbind(c(1:n),d.w[,i]))[-i,]
		d.sorted <- d.sorted[order(d.sorted[,2]),]
		
		h[i] <- mean(d.sorted[1:k,2]^2)
				
	}
	
	# compare with chi-square distribution

	return(h > (qchisq(alpha,df)))
	
}

rescale.den <- function(x) {
	n <- length(x)
	#range <- max(x,na.rm=TRUE) - min(x,na.rm=TRUE)
	#a.num <- sample(1:range,n,rep=TRUE)
	a.num <- x

	#layout(matrix(c(1:3)))
	a.den <- density(a.num,na.rm=TRUE)
	#plot(a.den)
	#plot(a.num,rep(0,n))



	s <- 0

	r <- rank(a.num,ties.method="first")
	a.num <- sort(a.num)
	a.num.1 <- a.num
	a.num.1[1] <- 0
	for (i in 2:n) {
		dif <- (a.num[i] - a.num[i - 1]) * (approx(a.den$x, a.den$y, a.num[i])$y)
		#cat(a.num[i] - a.num[i - 1], " --> ", dif, "\n")
		s <- s + dif
		a.num.1[i] <- s
	}
	#plot(a.num.1,rep(0,n))
	
	for (i in 1:n)
	{
		a.num[i] <- a.num.1[r[i]]
	}
	
	# range adjustment
	a.num <- a.num*(1/max(a.num))
	
	return(a.num)
}



ingroups.ind <- function(d,d.w,d.l,data,usesghomo=FALSE,weights=rep(1,length(d.l)),by.attr=FALSE)
{
   
   n <- nrow(d)
   if (n > 2)
   {
   gmat <- matrix(1,nrow=n,ncol=n) # matrix that holds individual in-/outgroups by row

   outliers <- s.outliers(d.w,k=1,alpha=0.9,ncol(data) - 1)   # detect outliers

   for (i in 1:n)
   {
	   # sort group-members by their distance to member i
	   d.sorted <- (cbind(c(1:n),d.w[,i]))[-i,]
	   d.sorted <- d.sorted[order(d.sorted[,2]),]
       
       g <- rep(2,n)    # initial outgroup (2's denote outgroup members)
       g[i] <- 1        # focal group-member is the first one in his ingroup
       gmat[i,] <- g    # store in-/outgroup assoc. in gmatrix
       s <- c()         # s is vector of maximized individual silhouette widths (currently not used)
       
       # by.attr=TRUE uses weighted average of single-attribute distances to calculate silhouette-widths, 
       # instead of Euclidean distances of attribute-vaues multiplied by weighting factors
       
       # silhouette width for outlier condition
       if (by.attr == FALSE) s.i.max <- silhouette.width(d,g,usesghomo=usesghomo)[i] else s.i.max <- silhouette.width.by.attr(d.l,g,usesghomo=usesghomo,weights=weights)[i]
       
       if (is.na(outliers[i]) | outliers[i] == FALSE)  
       {
          g[d.sorted[1,1]] <- 1   # place outgroup member nearest to i in i's ingroup
          gmat[i,] <- g           # store in gmat
          
          # determine initial silhouette width (2-member ingroup)
          if (by.attr == FALSE) s.i.max <- silhouette.width(d,groups=g,usesghomo=usesghomo)[i] else s.i.max <- silhouette.width.by.attr(d.l,g,usesghomo=usesghomo,weights=weights)[i]
          
          # add additional outgroup members to i's ingroup
          g.i <- g
          for (o in d.sorted[-1,1])
          {  	  
    	     g.i[o] <- 1  # g.i is i's tentative in-/outgroup distribution
    	     # silhouette-width with g.i
    	     if (by.attr == FALSE) s.i <- silhouette.width(d,groups=g.i,individual=i,usesghomo=usesghomo)[i] else s.i <- silhouette.width.by.attr(d.l,g.i,usesghomo=usesghomo,weights=weights)[i]

             if (s.i >= s.i.max) # do we have a new max?
    	     {
    	     	   # register new max
    	     	   g <- g.i  
    	     	   gmat[i,] <- g
    	     	   s.i.max <- s.i
	    
    	     } # else break (removed in V1.52)
          }
          s <- c(s,s.i.max)
       }


   }
   } else gmat <- 0
   return(gmat)
}


#################################
# function: var.ratio.ws.bs
#################################
# Version:    1.01
# Date:       18.01.2012
# Author:     A. Glenz
#
# Purpose:
# Calculate variability ratio as ratio of the within-subgroup sum-of-squares 
# over the between-subgroup sum-of-squares
#
# Called by:             command-line
# Calls:                 
#
# Parameters:
#            data:       Data set as data.frame with attributes in colums and 
#                        data-records in rows.
# Return value:          variability ratio
#
# Version 1.01:
# Separation of group association vector from data matrix
################################# 
var.ratio.ws.bs <- function(data,groups)
{
	    centroids <- by(data,groups,mean)   # subgroup centroids
	    centroids <- do.call(rbind,centroids)          # convert to data frame	    
	    
	    # sum of squares ratio
	    qs.matrix <- data
	    for (i in 1:nrow(data)) qs.matrix[i,] <- (data[i,] - centroids[groups[i],])^2
	    qs.w <- sum(qs.matrix)
	    
	    ct <- sapply(data,mean)
	    qs.b <- 0
	    
	    for (k in unique(groups)) qs.b <- qs.b +
	                                    sum((centroids[k,] - ct)^2 *
	                                    (sum(groups == k)))
	    ratio <- qs.w / qs.b
	     
	    return(ratio)
}





#################################
# function: fau_numerator
#################################
# Version:    1.01
# Date:       18.01.2012
# Author:     A. Glenz
#
# Purpose:
# Calculate between-subgroup sum-of-squares as the numerator of Thatcher et al.'s 
# formula for Fau, based on Euclidean distances between individuals.
#
# Called by:             command-line
# Calls:  
#
# Parameters:
#            data:       Data set as data.frame with attributes in colums and 
#                        data-records in rows. 
#            means:      Attribute's mean values
#
# Return value:          Between-subgroups sums-of-squares
#
# Version 1.01:
# Separation of group association vector from data matrix
################################# 
fau.numerator <- function(data,groups,means)
{
   p <- ncol(data)
   
   sum_j <- 0
   for (j in 1:p)
   {
      x.j. <- means[j]
      sum_k <- 0
      for (k in unique(groups))
      {
         x.jk <- mean(data[groups == k,j])
         n <- sum(groups == k)
         sum_k <- sum_k + (n * (x.jk - x.j.)^2)
      }      
      sum_j <- sum_j + sum_k
   }
   return(as.numeric(sum_j))
}


   
#################################
# function: fau_denominator
#################################
# Version:    1.01
# Date:       18.01.2012
# Author:     A. Glenz
#
# Purpose:
# Calculate total subgroup sum-of-squares as the denominator of Thatcher et al.'s 
# formula for Fau, based on Euclidean distances between individuals.
#
# Called by:             command-line
# Calls:  
#
# Parameters:
#            data:       Data set as data.frame with attributes in colums and 
#                        data-records in rows. 
#
# Return value:          Total subgroups sums-of-squares
#
# Version 1.01:
# Separation of group association vector from data matrix
#################################   
fau.denominator <- function(data)  
{
      
   sum.x <- 0
   means <- sapply(data,mean)
   
   for (i in 1:nrow(data))  sum.x <- sum.x + (data[i,] - means)^2 
   
   return(sum(sum.x))	   
}





#################################
# function: random.assignment
#################################
# Version:    1.01
# Date:       18.01.2012
# Author:     A. Glenz
#
# Purpose:
# Calculate random partition (subgroup assignment) of a given group
#
# Called by:             asw_cluster.random
# Calls:  
#
# Parameters:
#            data:       Data set as data.frame with attributes in colums and 
#                        data-records in rows. 
#            ng:         Number of subgroups
#
# Return value:          Vector of n integer values, containing randomized subgroup
#                        assignments of a n-sized group's members
#
# Version 1.01:
# Separation of group association vector from data matrix
#################################
random.assignment <- function(data,ng)
{
   groups <- c()
   while (length(unique(groups)) != ng) groups <- sample(1:ng,nrow(data),replace=TRUE)
   return(groups)
}




#################################
# function: reord.groups
#################################
# Version:    1.0
# Date:       22.12.2011
# Author:     A. Glenz
#
# Purpose:
# Re-arrange and sort subgroup assignment values in a way that
# there are no gaps between existing subgroup id's
#
# Called by:             asw_cluster.random
#                        asw_cluster.agglomerative
# Calls:  
#
# Parameters:
#            groups:     Vector of n unsorted integer values, representing
#                        subgroup assignments of an n-sized group 
#
# Return value:          Sorted and re-arranged vector of n integer values, subgroup
#                        assignments of a n-sized group's members
#
#################################
reord.groups <- function(groups)
{
   g <- 1
   i <- 1
   gnew <- groups
   ng <- length(unique(groups))
   while ((g <= ng) && (i <= length(gnew)))
   {
      if (gnew[i] > g)
      {
   	     x <- gnew[i]
         gnew[gnew == g] <- -1
         gnew[gnew == x] <- g
         gnew[gnew == -1] <- x
         g <- g+1
      } else if (gnew[i] == g) g <- g+1
   
      i <- i+1
   }
   return(gnew)
} 




#################################
# function: fau.dist   #
#################################
# Version:    1.01
# Date:       18.01.2012
# Author:     A. Glenz
#
# Purpose:
# Determine faultline distance as the average euclidean distance
# between centroids of a given group partition
#
# Parameters:
#            data:    Data set as data.frame with attributes in colums 
#                     and data-records in rows. 
#            groups:  Group-partition with integer values 
#                     representing the individua's assigned subgroup id
#
# Return value:       faultline.distance, as suggested by Bezrukova et al. (2009)
#
# Version 1.01:
# Separation of group association vector from data matrix
#################################
fau.dist <- function(data,groups)
{
    groups_uq <- unique(groups) # subgroups in data
    ng <- length(groups_uq)
    
    if (ng > 1)
    {
       # calculate subgroup centroids (vectors of means of each variable)
       centroids <- matrix(nrow=max(groups_uq),ncol=ncol(data))
       for (i in groups_uq)
       {
          centroids[i,] <- apply(data[groups == i,],2,mean)
       }
    
       fau.d <- 0 # initialize fau-distance
       for (i in groups_uq)
       {
          for (j in groups_uq)
          {
             if(i < j) 
             {
          	    fau.d <- fau.d + sqrt(sum((centroids[i,] - centroids[j,])^2)) # euclidean distance of centroid-pair
             }
          }
       }
       fau.d <- fau.d / (ng*(ng-1)/2) # average distance
    } else fau.d <- 0
    return(fau.d)      
}



#################################
# function: fau.dist.mahal   #
#################################
# Version:    1.0
# Date:       5.04.2012
# Author:     A. Glenz
#
# Purpose:
# Determine faultline distance as the average euclidean distance
# between centroids of a given group partition
#
# Parameters:
#            data:    Data set as data.frame with attributes in colums 
#                     and data-records in rows. 
#            groups:  Group-partition with integer values 
#                     representing the individua's assigned subgroup id
#            S:       Attribute's inverted variance-/covariance matrix
#
# Return value:       faultline.distance, as suggested by Bezrukova et al. (2009),
#                     but based on Mahalanobis distances
#
#################################
fau.dist.mahal <- function(data,groups,S)
{
    groups_uq <- unique(groups) # subgroups in data
    ng <- length(groups_uq)

    
    if (ng > 1)
    {
       # calculate subgroup centroids (vectors of means of each variable)
       centroids <- matrix(nrow=max(groups_uq),ncol=ncol(data))
       for (i in groups_uq)
       {
          centroids[i,] <- apply(data[groups == i,],2,mean)
       }
    
       fau.d <- 0 # initialize fau-distance
       for (i in groups_uq)
       {
          for (j in groups_uq)
          {
             if(i < j) 
             {
          	    fau.d <- fau.d + sqrt(dist.m(centroids[i,], centroids[j,],S)) # mahalonobis distance of centroid-pair
             }
          }
       }
       fau.d <- fau.d / (ng*(ng-1)/2) # average distance
    } else fau.d <- 0
    return(fau.d)      
}




#################################
# function: fau_cluster.mahal   #
#################################
# Version:    1.01
# Date:       18.01.2012
# Author:     A. Glenz
#
# Purpose:
# Determine subgroup partition with the highest beween-subgroup
# sum-of-squares, which is the numerator in the FAU measure
# according to Thatcher et al. (2003). 
# Algorithm uses incremental one-by-one clustering. 
# Calculation is based on Mahalanobis distances.
#
# Parameters:
#            data:    Data set as data.frame with attributes in colums 
#                     and data-records in rows. 
#            groups:  Group-partition with integer values 
#                     representing the individua's assigned subgroup id
#            gmax:    maximum number of subgroups allowed
#            S:       Attribute's inverted variance-/covariance matrix
#
# Return value:       Data vector containing the maximum value for the
#                     between-subgroup sum-of-squares in the 
#                     first place, followed by the subgroup partition 
#                     (n values) for this value.
#
# Version 1.01:
# Separation of group association vector from data matrix
#################################
fau_cluster.mahal <- function(data,groups,gmax,S)
{
   # start with minimal value 
   maxf <- -1    
   
   
   # n = group size
   n <- nrow(data)
   
   # nv = number of variables in data
   nv <- ncol(data)
   
   # initialize stop-condition
   imax <- 1   
   
   # loop as long as fau is growing
   while (imax > 0)    
   {
         # establish stop-condition
         imax <- 0      

         for (i in 1:n)  # for each person
         {
   	        for (newgroup in c(gmax:1)) # group to put people into
   	        {
   	             
   	           if (groups[i] != newgroup) # if group changes
   	           { 	              
    	              # g_new is a temporary variable for the partition
    	              g_new <- groups   
   	              
   	              # put person i into group "newgroup"
   	              g_new[i] <- newgroup  
   	              
   	                	              
   	              # consider only partitions with 2 and more subgroups 
   	              if (length(unique(g_new)) >= 2)
   	              {

   	                 # maximizing numerator is enough, because denominator is
   	                 # constant for app partitions
   	                 f <- fau.numerator.mahal(data, g_new, as.vector(sapply(data,mean)),S)
                      
   	                 
   	                 # if the new fau exceeds the previous maximum
   	                 if (f > maxf) 
   	                 {
                               # store new maximized value (highscore)
   	  	                      maxf <- f 
   	  	                      

   	  	                      # store id of person whose new assignment
   	  	                      # caused the new maximal value
   	  	                      imax <- i 
   	  	                      
   	  	                      # store better subgroup id for person i
   	  	                      newgroupmax <- newgroup  
   	  	                
   	  	            } 
                   }
   	           }
   	        } 
         }
         if (imax > 0) # did we faund a new fau maximum?
         { 
         	   # write new subgroup assignment back to data set
         	   groups[imax] <- newgroupmax
         	   
         } 
   }
   maxf <- fau.numerator.mahal(data, groups, as.vector(sapply(data,mean)),S)
   return(c(maxf,groups))                                               
}


#################################
# function: fau_cluster         #
#################################
# Version:    1.01
# Date:       18.10.2012
# Author:     A. Glenz
#
# Purpose:
# Determine subgroup partition with the highest beween-subgroup
# sum-of-squares, which is the numerator in the FAU measure
# according to Thatcher et al. (2003). 
# Algorithm uses incremental one-by-one clustering. 
# Calculation is based on Euclidean distances.
#
# Parameters:
#            data:    Data set as data.frame with attributes in colums 
#                     and data-records in rows. 
#            groups:  Group-partition with integer values 
#                     representing the individua's assigned subgroup id
#            gmax:    maximum number of subgroups allowed
#
# Return value:       Data vector containing the maximum value for the
#                     between-subgroup sum-of-squares in the 
#                     first place, followed by the subgroup partition 
#                     (n values) for this value.
#
# Version 1.01:
# Separation of group association vector from data matrix
#################################
fau_cluster <- function(data, groups, gmax)
{
   # start with minimal value 
   maxf <- -1    
   
   # n = group size
   n <- nrow(data)
   
   # nv = number of variables in data
   nv <- ncol(data)
   
   # initialize stop-condition
   imax <- 1   
   
   # loop as long as fau is growing
   while (imax > 0)    
   {
         # establish stop-condition
         imax <- 0      

         for (i in 1:n)  # for each person
         {
   	        for (newgroup in c(gmax:1)) # group to put people into
   	        {
   	             
   	           if (groups[i] != newgroup) # if group changes
   	           { 	              
    	              # g_new is a temporary variable for the partition
    	              g_new <- groups   
   	              
   	              # put person i into group "newgroup"
   	              g_new[i] <- newgroup  

   	                 	              
   	              # consider only partitions with 2 and more subgroups 
   	              if (length(unique(g_new)) >= 2)
   	              {

   	                 # maximizing numerator is enough, because denominator is
   	                 # constant for app partitions
   	                 f <- fau.numerator(data, g_new, as.vector(sapply(data,mean)))

   	                 
   	                 # if the new fau exceeds the previous maximum
   	                 if (f > maxf) 
   	                 {
                               # store new maximized value (highscore)
   	  	                      maxf <- f 


   	  	                      # store id of person whose new assignment
   	  	                      # caused the new maximal value
   	  	                      imax <- i 
   	  	                      
   	  	                      # store better subgroup id for person i
   	  	                      newgroupmax <- newgroup  
   	  	                
   	  	            } 
                   }
   	           }
   	        } 
         }
         if (imax > 0) # did we faund a new fau maximum?
         { 
         	   # write new subgroup assignment back to data set
         	   groups[imax] <- newgroupmax
         	   
         } 
   }
   maxf <- fau.numerator(data, groups, as.vector(sapply(data,mean)))
   return(c(maxf,groups))                                               
}



#################################
# function: fau_cluster.random
#################################
# Version:    1.01
# Date:       18.01.2012
# Author:     A. Glenz
#
# Purpose:
# Invoke fau-maximizing clustering procedure starting at random start.
# configurations. This function is very similar to asw_cluster.random,
# but here, the value of Thatcher's FAU (according to Thatcher, 2003)
# is maximized. ASW is used to determine the optimal number of
# subgroups.
#
# 
# configurations. 
#
# Called by:             command-line
# Calls:                 S.mahal
#                        C.cor
#                        reord.groups
#                        fau_cluster
#                        fau_cluster.mahal
#                        fau.denominator
#                        fau.denominator.mahal
#                        distmatrix.e
#                        distmatrix.m
#                        silhouette.width
#                        random.assignment
#
# Parameters:
#            data:       Data set as data.frame with attributes in colums and 
#                        data-records in rows. 
#            ng:         Maximum/exact number of subgroups to detect, depending 
#                        on value of numgrp.
#            maxnuc:     Number of unsuccessful attempts to find a new local 
#                        maximum for asw
#            verbose:    Verbosity level from 0=silent to 3=talkative
#            metric:     if metric = "mahal", partitions will be evaluated based 
#                        on Mahalanobis-distances
#                        if metric = "corr", data rows (attributes) will be 
#                        weighted with 1 - cm, where cm denotes the column mean 
#                        of the pxp correlation matrix, and p denotes the
#                        number of attributes in the data set 
#                        For all other values of metric, including NA, partitions 
#                        will be evaluated based on Euclidean distances
#            numgrp:     if set to "auto" (the default), the optimal number of 
#                        subgroups between 2 and ng will be detected. For all 
#                        other values, exactly ng subgroups will be determined.
#
# Return vaue:           Data vector containing the maximum value for asw in the 
#                        first place, followed by the subgroup partition (n values) 
#                        for this value.
#
# Version 1.01:
# Separation of group association vector from data matrix
#################################
fau_cluster.random <- function(data,ng,maxnuc=5,verbose=1,metric='euklid',numgrp='auto')
{
   # Initialize...	   
   maxasw <- -1
   # current number of unsuccessful attempts to find new local maximum
   nuc <- 0 
   maxfaugroups <- rep(1,nrow(data)) # initial subgroup associations
   res <- list() # initialize result list

   
  
   # determine group range
   ifelse (numgrp == 'auto', gmin <- 2, gmin <- ng) 
   gmax <- ng
   
      
   if (metric == 'mahal')
   {
      S <- S.mahal(data)
      D <- distmatrix.m(data,S)
      f.denom <- fau.denominator.mahal(data,as.vector(sapply(data,mean)),S)
   } else
   {
      if (metric == 'corr')
      {
         C <- C.cor(data) * diag(ncol(data))
         data <- as.matrix(data) %*% C
      }
      D <- distmatrix.e (data)
      f.denom <- fau.denominator(data)
   }
     
   ng <- gmax 

   
   while (nuc < maxnuc)
   {
      groups <- random.assignment(data,ng) # new random group assignment
      #cat(data$groupID,"   ")
      ng.out <- ng
      if (nrow(data) > 1)
      {
      	if (metric == 'mahal')
      	{
      	   fcl <- fau_cluster.mahal(data,groups,ng,S)
      	} else  fcl <- fau_cluster(data,groups,ng)
      } else fcl <- c(0,1)
      
      f <- round(fcl[1],digits=10)
      groups.x <- reord.groups(fcl[-1])
      asw <- round(mean(silhouette.width(D,groups.x)), digits=10)
      #cat(groups.x,"--> fau =",faud,"\n")
      newmax <- TRUE # initialize stop condition
      #cat(asw,maxasw,"\n")
      if (length(res) > 0)
      {
         for (i in 1:length(res))
         {
            # we've already found this local maximum before
            if (f == res[[i]][2]) newmax <- FALSE 
         }	
      }	
   
      if (newmax)
      {
         
         groups <- groups.x
         d <- 1
         if (asw >= maxasw)
         {
         	   maxasw <- asw
         	   maxaswng <- ng
         	   maxfaugroups <- groups.x
         	   maxfau <- f
            if (length(res) > 0)
            {
               for (i in 1:length(res))
               {
                  rf <- as.numeric(res[[i]][2])
                  if ((maxaswng == as.numeric(res[[i]][3])) & (rf > maxfau))  
                  {
                     maxfau <- rf
                     maxfaugroups <- res[[i]][1][[1]]
                  }
               }	
            }	         	   
         } else
         {
            if (f > maxfau)
            {
               maxfau <- f
               maxfaugroups <- groups.x
               maxdist <- d
               if (ng < gmax)
               {
            	    ng <- ng + 1
            	    if (ng > gmax) ng <- gmax
               }	 
            }
         }
         
         
         res[[length(res)+1]] <- list(groups.x,f,length(unique(groups.x)))

         fau <- f / f.denom
         
         if ((verbose > 2) && (nuc > 0)) cat("\n")
         if (verbose > 1) cat("New local maximum: Fau*Dist =",fau,
                              " Maxfau =",maxfau / f.denom,
                              " Maxasw=",maxasw,
                              " maxaswng=",maxaswng,
                              " Groups =",groups.x,
                              " Fau =",fau/d,
                              " Dist =",d,
                              " ASW:",asw,
                              " # groups:",ng.out,
                              "->",length(unique(groups.x)),
                              "\n"); flush.console()
         nuc <- 0
         newmax <- FALSE
      } else
      {
         nuc <- nuc+1
         if (verbose > 2) cat("Already known local maximum,  # groups:",
                              length(unique(groups.x)),
                              "(attempt",nuc,"of",maxnuc,") \n")
         if((nuc == maxnuc) & (ng > gmin))
         {           
            nuc <- 0 
            ng <- ng - 1          
         }
      }
   }
   groups <- maxfaugroups
   maxasw <- mean(silhouette.width(D,groups)) 
   maxfau <- maxfau / f.denom  
   if (verbose > 0) cat("Maximum Fau:",maxfau,
                        "  Groups:",groups,
                        " ASW:",maxasw,
                        " # groups:",length(unique(groups)),
                        "\n")
   return(c(maxfau,groups,maxasw))
}





#################################
# function: asw_cluster         #
#################################
# Version:    1.01
# Date:       19.01.2012
# Author:     A. Glenz
#
# Purpose:
# Determine subgroup partition with the highest asw-value by 
# incremental one-by-one clustering. 
#
# Called by:           asw_cluster.random
#                      asw_cluster.agglomerative
# Calls:               silhouette.width 
#
# Parameters:
#            data:    Data set as data.frame with attributes in colums and 
#                     data-records in rows. 
#            groups:  Group-partition with integer values 
#                     representing the individua's assigned subgroup id
#            gmax:    maximum number of subgroups allowed
#            D:       n x n Distance matrix, containing the pairwise distances 
#                     between group memebers. Distances may base on various 
#                     distance metrics and can be weighted.
#
# Return value:       Data vector containing the maximum value for asw in 
#                     the first place, followed by the subgroup partition 
#                     (n values) for this value.
#
# Version 1.01:
# Separation of group association vector from data matrix
#################################
asw_cluster <- function(data,groups,gmax,D,individual="avg",usesghomo=FALSE,weights=NA,by.attr=FALSE)
{
   maxasw <- -1    # start with minimal value for asw 
   
   
   
   n <- length(groups)
   nv <- ncol(data)
   
   # allow up to gmax subgroups
   nofgroups <- min(gmax,max(groups)+1)
   
   # initialize stop-condition
   imax <- 1                                                  
   
   # loop as long as target value is growing
   while (imax > 0)                                         
   {
         # establish stop-condition
         imax <- 0  
         
         # for each person
         for (i in 1:n)                                      
         {
   	        # for each subgroup to put people into
   	        for (newgroup in c(nofgroups:1))                
   	        {
   	             
   	           # if group assignment changes...
   	           if (groups[i] != newgroup)             
   	           { 	              
    	              # g_new is a temporary variable for the new subgroup assignment
    	              g_new <- groups                     
   	              
   	              # put person i into group "newgroup"
   	              g_new[i] <- newgroup  
                    
   	                 	                 
   	              # if the new assignment produces more than two subgroups
   	              if (length(unique(g_new)) >= 2)           
   	              {
    	                 # calculate silhouette widths...
    	                 if (by.attr == TRUE) sil <- silhouette.width.by.attr(D,g_new,individual=individual,usesghomo=usesghomo,weights=weights) 
    	                 else sil <- silhouette.width(D,g_new,individual=individual,usesghomo=usesghomo)  
    	                      
   	                 if (individual == "avg")
   	                 {
   	                 	 # ... and asw as their mean value
   	                    asw <- mean(sil)
   	                 } else
   	                 {
   	                    asw <- sil[individual]
   	                 }
   	                 
   	                 # alternative version which excludes zero-silhouettes from mean                       
   	                 #ifelse(sum(sil != 0) > 0, asw <- sum(sil) / sum(sil != 0), asw <- 0)  
                    
                      
   	                 # if the new value for asw reaches or exceeds the previous maximum
   	                 if (asw > maxasw)                      
   	                 {      
   	  	                  
   	  	                   # store new maximized asw value (highscore). 
   	  	                   # Note that this line is for clarification only and could be omitted.
   	  	                   maxasw <- asw                   

   	  	                   
   	  	                   # store person whose new group assignment caused the maximal value       
   	  	                   imax <- i   
   	  	                   
   	  	                   # store better group for person i                     
   	  	                   newgroupmax <- newgroup         	                
   	  	             } # end if
                   } # end if
   	           } # end if
   	        } # end for
         } # end for

         # did we find a new fau maximum in this round?
         if (imax > 0)                                      
         { 
         	   # write new optimal partition back to grouping vector...
         	   groups[imax] <- newgroupmax               
         	                             
         }             	                
   } # end while
   # cat(maxasw,"   ",groups,"\n")
   # maxasw <- mean(silhouette.width(D,groups))

   return(c(maxasw,groups))  # that's it...                                                                  
} # end function


#################################
# function: asw_cluster.random
#################################
# Version:    1.01
# Date:       19.01.2012
# Author:     A. Glenz
#
# Purpose:
# Invoke asw-maximizing clustering procedure starting at random start 
# configurations. 
#
# Called by:             command-line
# Calls:                 S.mahal
#                        distmatrix.e
#                        distmatrix.m
#                        C.cor
#                        reord.groups
#                        asw_cluster
#                        random.assignment
#
# Parameters:
#            data:       Data set as data.frame with attributes in colums and 
#                        data-records in rows. 
#            ng:         Maximum/exact number of subgroups to detect, depending 
#                        on value of numgrp.
#            maxnuc:     Number of unsuccessful attempts to find a new local 
#                        maximum for asw
#            verbose:    Verbosity level from 0=silent to 3=talkative
#            metric:     if metric = "mahal", partitions will be evaluated based 
#                        on Mahalanobis-distances
#                        if metric = "corr", data rows (attributes) will be 
#                        weighted with 1 - cm, where cm denotes the column mean 
#                        of the pxp correlation matrix, and p denotes the
#                        number of attributes in the data set 
#                        For all other values of metric, including NA, partitions 
#                        will be evaluated based on Euclidean distances
#            numgrp:     if set to "auto" (the default), the optimal number of 
#                        subgroups between 2 and ng will be detected. For all 
#                        other values, exactly ng subgroups will be determined.
#
# Return vaue:           Data vector containing the maximum value for asw in the 
#                        first place, followed by the subgroup partition (n values) 
#                        for this value.
#
# Version 1.01:
# Separation of group association vector from data matrix
#################################
asw_cluster.random <- function(data,ng,maxnuc=5,verbose=1,metric='euclid',numgrp='auto')
{
   # Initialize...	   
   maxasw <- -1
   
   # current number of unsuccessful attempts to find new local maximum
   nuc <- 0                                           
   
   # initial subgroup associations
   maxaswgroups <- rep(1,nrow(data))  
   
   # initialize result list                     
   res <- list()                                      
   
   
  
   # determine group range
   ifelse (numgrp == 'auto', gmin <- 2, gmin <- ng) 
   gmax <- ng
   
   # determine distance matrix according to selected metric   
   if (metric == 'mahal')
   {
      S <- S.mahal(data)
      D <- distmatrix.m(data,S)
   } else
   {
      if (metric == 'corr')
      {
         C <- C.cor(data) * diag(ncol(data))
         data <- as.matrix(data) %*% C
      }
      D <- distmatrix.e (data)
   }
     
 

   # loop over selected number of attempts
   while (nuc < maxnuc)                               
   {
      # new random group assignment
      groups <- random.assignment(data,ng)      

      # numer of subgroups (for output messages)
      ng.out <- ng                                    
      
      # if data set contains data for more than 1 individual...
      if (nrow(data) > 1)                             
      {
      	# ...call clustering procedure...
      	fcl <- asw_cluster(data,groups,ng,D)                
      } else fcl <- c(0,1) # ... or otherwise report constant result                           
      
      # sort and arrange subgroup numbers
      groups.x <- reord.groups(fcl[-1])              
      
      asw <- round(fcl[1], digits=10)                

      # initialize break condition
      newmax <- TRUE                                  
      
      # if more than 1 solution determined so far
      if (length(res) > 0)                            
      {
         # search through already found solutions
         for (i in 1:length(res))                     
         {
            # we've already found this local maximum before
            if (asw == res[[i]][2]) newmax <- FALSE   
         }	
      }	
   
      # do we have a new local maximum?
      if (newmax)                                     
      {
         # write subgroup partition to data set
         groups <- groups.x                       

         if (asw > maxasw)
         {
            # store maximum asw value
            maxasw <- asw                             
            
            # store subgroup partition for maximum asw
            maxaswgroups <- groups.x                    

            # if we have not yet reached the maximum number of allowed subgroups...
            if (ng < gmax)                            
            {
            	 # ...the number of subgroups increases by 1...
            	 ng <- ng + 1                           
            	 
            	 # ...but must not exceed gmax
            	 if (ng > gmax) ng <- gmax              
            }	 
         }

         # store result for this random start configuration
         res[[length(res)+1]] <- list(groups.x,asw)     
                           

         if ((verbose > 2) && (nuc > 0)) cat("\n")  
         if (verbose > 1) cat("New local maximum: Groups =",groups.x," ASW:", asw,
                              " # groups:",ng.out,
                              "->",length(unique(groups.x)),
                              "\n"); flush.console()
         
         # reset counter of unsuccessful attempts
         nuc <- 0                                     
         newmax <- FALSE                              
      } else # we have found a previously detected local maximum again                                         
      {
         # increase number of unsuccessful attempts
         nuc <- nuc+1                                 
         if (verbose > 2) cat("Already known local maximum,  # groups:",
                              length(unique(groups.x)),
                              "(attempt",nuc,"of",maxnuc,") \n")
         
         # we did not find a new maximum with this number of subgroups...
         if((nuc == maxnuc) & (ng > gmin))            
         {           
            # ...so we try again with a number of subgroups...
            nuc <- 0    
            
            # ...decreased by 1                              
            ng <- ng - 1                              
         }
      }
   }
   
   # write subgroup partition with maximum asw value back to data set
   groups <- maxaswgroups                       

   
   if (verbose > 0) cat("Groups:",groups,
                        "   ASW:",maxasw,
                        " # groups:",length(unique(groups)),
                        "\n")
   #return(c(maxasw,groups))
   list(asw = maxasw,
        groups = groups)
}






#################################
# function: fau.numerator.mahal
#################################
# Version:    1.01
# Date:       19.01.2012
# Author:     A. Glenz
#
# Purpose:
# Calculate between-subgroup sum-of-squares as the numerator of Thatcher et al.'s 
# formula for Fau, based on Mahalanobis-distances between individuals.
#
# Called by:             
# Calls:  
#
# Parameters:
#            data:       Data set as data.frame with attributes in colums and 
#                        data-records in rows. 
#            groups:     Group-partition with integer values 
#                        representing the individua's assigned subgroup id
#            means:      Attribute's mean values
#            S:          Attribute's inverted variance-/covariance matrix
#
# Return value:          Between-subgroups sums-of-squares
#
# Version 1.01:
# Separation of group association vector from data matrix
#################################
fau.numerator.mahal <- function(data,groups,means,S)
{
   sum_k <- 0
   for (k in unique(groups))
   {
         x.jk <- sapply(data[groups == k, ],mean)
         n <- sum(groups == k)
         sum_k <- sum_k + (dist.m(x.jk,means,S) * n)
   }      
   return(as.numeric(sum_k))
}



#################################
# function: fau.denominator.mahal
#################################
# Version:    1.01
# Date:       19.01.2012
# Author:     A. Glenz
#
# Purpose:
# Calculate total subgroup sum-of-squares as the denominator of Thatcher et al.'s 
# formula for Fau, based on Mahalanobis-distances between individuals.
#
# Called by:             
# Calls:  
#
# Parameters:
#            data:       Data set as data.frame with attributes in colums and 
#                        data-records in rows. 
#            means:      Attribute's mean values
#            S:          Attribute's inverted variance-/covariance matrix
#
# Return value:          Total subgroups sums-of-squares
#
# Version 1.01:
# Separation of group association vector from data matrix
#################################
fau.denominator.mahal <- function(data,means,S)  
{
   n <- nrow(data)
   
   sum_j <- 0
   for (j in 1:n)
   {
         x.j <- data[j,]
         sum_j <- sum_j + dist.m(x.j,means,S)
   }      
   return(sum_j)
}






#################################
# function: dist.g
#################################
# Version:    1.0
# Date:       22.12.2011
# Author:     A. Glenz
#
# Purpose:
# Calculate distance between two individuals according to gower metric (manhattan)
#
# Called by:             command-line
# Calls:  
#
# Parameters:
#            x:          Vector of length p containing the individual X's
#                        values for all p attributes
#            y:          Vector of length p containing the individual Y's
#                        values for all p attributes
#
# Return value:          Manhattan distance between X and Y
#
#################################
dist.g <- function(x,y)
{
   x <- as.numeric(x)
   y <- as.numeric(y)
   d <- sum(abs(x-y))   
   return(d)
} 


#################################
# function: dist.c
#################################
# Version:    1.0
# Date:       22.12.2011
# Author:     A. Glenz
#
# Purpose:
# Calculate Euclidean distance between two individuals, with attributes
# weighted by their intercorrelations
#
# Called by:             
# Calls:  
#
# Parameters:
#            x:          Vector of length p containing the individual X's
#                        values for all p attributes
#            y:          Vector of length p containing the individual Y's
#                        values for all p attributes
#            C:          Vector of length p, containing  the value of 1 - cm, where 
#                        cm denotes the column mean of the pxp correlation matrix, and p denotes the
#                        number of attributes in the data set 
#
# Return value:          Weighted Euclidean distance between X and Y
#
#################################
dist.c <- function(x,y,C)
{
   x <- as.numeric(x)
   y <- as.numeric(y)
   d <- ((x-y)^2 %*% C)   
   return(d)
}




#################################
# function: dist.m
#################################
# Version:    1.0
# Date:       22.12.2011
# Author:     A. Glenz
#
# Purpose:
# Calculate Mahalanobis distance between two individuals
#
# Called by:             distmatrix.m
#                        fau.denominator.mahal
# Calls:  
#
# Parameters:
#            x:          Vector of length p containing the individual X's
#                        values for all p attributes
#            y:          Vector of length p containing the individual Y's
#                        values for all p attributes
#            S:          Attribute's inverted variance-/covariance matrix
#            weights:    Vector of length p, containing weighting factors for the p attributes
#
# Return value:          Mahalanobis distance between X and Y
#
#################################
dist.m <- function(x,y,S,weights=NA)
{
    if (is.na(sum(weights)))
    {
    	   weights <- rep(1,ncol(S))
    }  
   x <- as.numeric(x)
   y <- as.numeric(y)
   d <- abs(((x-y) * weights) %*% S %*% ((x-y) * weights))
   return(d)
} 


#################################
# function: dist.e
#################################
# Version:    1.0
# Date:       22.12.2011
# Author:     A. Glenz
#
# Purpose:
# Calculate Euclidean distance between two individuals
#
# Called by:             command-line
# Calls:  
#
# Parameters:
#            x:          Vector of length p containing the individual X's
#                        values for all p attributes
#            y:          Vector of length p containing the individual Y's
#                        values for all p attributes
#
# Return value:          Euclidean distance between X and Y
#
#################################
dist.e <- function(x,y)
{
   x <- as.numeric(x)
   y <- as.numeric(y)
   d <- sum((x-y)^2)   
   return(d)
} 


#################################
# function: C.cor
#################################
# Version:    1.0
# Date:       22.12.2011
# Author:     A. Glenz
#
# Purpose:
# Calculate Vector of weights to determine correlation-weighted distance 
# between two individuals
#
# Called by:             command-line
#                        asw_cluster.agglomerative
# Calls:  
#
# Parameters:
#            data:       Data set as data.frame with attributes in colums and 
#                        data-records in rows. 
#
# Return value:          Vector of weights to determine correlation-weighted 
#                        distance between X and Y
#
#################################
C.cor <- function(data)
{
   return(1 - (apply(abs(cor(data) - diag(ncol(data))),2,sum)  / (ncol(data) - 1)))
}




#################################
# function: rel.dem
#################################
# Version:    1.0
# Date:       09.02.2012
# Author:     A. Glenz
#
# Purpose:
# Calculate Vector of Euclidean distances of group members
# to the group's centroid to represent a relational demography measure
#
# Called by:             command-line
#                        
# Calls:  
#
# Parameters:
#            data:          Data set as data.frame with attributes in colums and 
#                           data-records in rows. 
#            exclude.self:  Logical (TRUE/FALSE), indicating whether the individual's 
#                           attribute value should be excluded from the calculation
#                           of the group's centroids (default=TRUE)
#            weights:       Vector (length=p) of weighting factors for the p attributes
#
# Return value:             Vector of Euclidean distances
#
#################################
rel.dem <- function(data, exclude.self=FALSE, weights=NA)
{
   n <- nrow(data)
  
   # appy weighting factors
   if (!is.na(weights) && (is.vector(weights)) && (length(weights) == ncol(data)))
   {
       data <- t(t(data) * weights)
   }
   
   if (exclude.self)
   {   
   	  centroid <- sapply(data,mean)
       diff_2 <- (t(t(data) - centroid))^2
       distances <- sqrt(apply(diff_2,1,sum))
   } else
   {
   	  distances <- vector(length=nrow(data),mode="numeric")
   	  for (i in 1:n)
   	  {
   	      centroid <- sapply(data[-i,],mean)
   	      distances[i] <- sqrt(sum((data[i,] - centroid)^2)) 
   	  } 
   }
   
   Dscores <- data * 0   
   for (attr in 1:ncol(data))
   {
   	   for (i in 1:n)
   	   {
   	   	   Dscores[i,attr] <- sqrt(sum((data[,attr] - data[i,attr])^2) / n)
   	   }
   }   
      
   return(list(distances=distances,dscores=Dscores))
}



#################################
# function: S.mahal
#################################
# Version:    1.0
# Date:       22.12.2011
# Author:     A. Glenz
#
# Purpose:
# Calculate inverted variance-/covariance-matrix
# needed to determine Mahalanobis distance
# between two individuals
#
# Called by:             command-line
#                        asw_cluster.agglomerative
#                        asw_cluster.random
# Calls:  
#
# Parameters:
#            data:       Data set as data.frame with attributes in colums and 
#                        data-records in rows. 
#
# Return value:          Inverted p x p variance-/covariance matrix of the p
#                        attributes in the data set
#################################
S.mahal <- function(data,dummycols=NA)
{
   #require(MASS)
   
   p <- ncol(data)
   S <- matrix(nrow=p,ncol=p)
   
   S.cov <- cov(data)
   if (!is.na(dummycols)) 
   {
      S.var <- diag(S.cov)
      S.cov[dummycols,dummycols] <- 0
      diag(S.cov) <- S.var   
   }
   
   if (abs(det(S.cov)) > 1e-17)
   {
      S <- solve(S.cov)
   } else S <- matrix(0,p,p)
   
   # use moore-penrose pseudo inversion to avoid errors when S is not invertible
   # S <- ginv(S.cov)

   
   return(S)
}



#################################
# function: distmatrix.m
#################################
# Version:    1.0
# Date:       22.12.2011
# Author:     A. Glenz
#
# Purpose:
# Calculate n x n Mahalanobis distance matrix, containing the pairwise Mahalanobis
# distances between the n group members
#
# Called by:             command-line
#                        asw_cluster.random
#                        asw_cluster.agglomerative
# Calls:                 dist.m
#
# Parameters:
#            data:       Data set as data.frame with attributes in colums and 
#                        data-records in rows. 
#            S:          Attribute's inverted variance-/covariance matrix
#            weights:    Vector of length p, containing weighting factors for the p attributes
#
# Return value:          Matrix of Mahalanobis distances
#################################
distmatrix.m <- function(data,S,weights=NA)
{   
    n <- nrow(data)
    D <- matrix(ncol=n,nrow=n)
    
    for (i in 1:n)
    {
       for (j in 1:n)
       {
          D[i,j] <- sqrt(dist.m(data[i,],data[j,],S,weights))
       }
    }
    
    return(D)
}



#################################
# function: distmatrix.e
#################################
# Version:    1.0
# Date:       22.12.2011
# Author:     A. Glenz
#
# Purpose:
# Calculate n x n Euclidean distance matrix, containing the pairwise Euclidean
# distances between the n group members
#
# Called by:             command-line
#                        asw_cluster.random
#                        asw_cluster.agglomerative
# Calls:  
#
# Parameters:
#            data:       Data set as data.frame with attributes in colums and 
#                        data-records in rows. 
#            weights:    Vector of length p, containing weighting factors for the p attributes
#
# Return value:          Matrix of Euclidean distances
#################################
distmatrix.e <- function(data,weights=NA)
{

    if (is.na(sum(weights)))
    {
    	   weights <- rep(1,ncol(data))
    }
    
    n <- nrow(data)
    D <- matrix(ncol=n,nrow=n)
    
    for (i in 1:n)
    {
       for (j in 1:n)
       {
          D[i,j] <- sqrt(sum(((data[i,] * weights) - (data[j,] * weights))^2))
       }
    }
    
    return(D)
}


# calculate absolute differences of values separately for each attribute
distmatrix.e.list <- function(data)
{

    DL <- list()
    for (i in 1:ncol(data))
    {
    	DL[[i]] <- abs(outer(data[,i],data[,i],"-"))
    }
    
    return(DL)
}



#################################
# function: distmatrix.g
#################################
# Version:    1.0
# Date:       22.12.2011
# Author:     A. Glenz
#
# Purpose:
# Calculate n x n Manhattan distance matrix (gower metric), containing the pairwise Manhattan
# distances between the n group members
#
# Called by:             command-line
#                        asw_cluster.random
#                        asw_cluster.agglomerative
# Calls:  
#
# Parameters:
#            data:       Data set as data.frame with attributes in colums and 
#                        data-records in rows. 
#            weights:    Vector of length p, containing weighting factors for the p attributes
#
# Return value:          Matrix of Manhattan distances
#################################
distmatrix.g <- function(data,weights=NA)
{

    if (is.na(sum(weights)))
    {
    	   weights <- rep(1,ncol(data))
    }
    
    n <- nrow(data)
    D <- matrix(ncol=n,nrow=n)
    
    for (i in 1:n)
    {
       for (j in 1:n)
       {
          D[i,j] <- (sum(abs((data[i,] * weights) - (data[j,] * weights))))
       }
    }
    
    return(D)
}

distmatrix.g.list <- function(data)
{

    DL <- list()
    for (i in 1:ncol(data))
    {
    	DL[[i]] <- abs(outer(data[,i],data[,i],"-"))
    }
    
    return(DL)
}

#################################
# function: distmatrix.w
#################################
# Version:    1.0
# Date:       22.10.2014
# Author:     A. Glenz
#
# Purpose:
# Calculate n x n matrix with weighted averages of single-attribute distances
# between the n group members
#
# Called by:             command-line
#                        asw_cluster.random
#                        asw_cluster.agglomerative
# Calls:  
#
#
# Return value:          Matrix of distances
#################################
distmatrix.w <- function(data,weights=NA)
{

    if (is.na(sum(weights)))
    {
    	   weights <- rep(1,ncol(data))
    }
    
    n <- nrow(data)
    D <- matrix(0,ncol=n,nrow=n)
    ranges <- diff(apply(data,2,range))
    
    
    for (attr in 1:ncol(data))
    {
         # standardize differences by range and multiply with weighting factor, sum up over all attributes
         if (ranges[attr] > 0) D <- D + weights[attr] * abs(outer(data[,attr],data[,attr],"-")) / ranges[attr] 
    }    
    # and calculate weighted average (divide by sum of weighting factors)
    D <- D / sum(weights)
    
    return(D)
}

#################################
# function: generate.random
#################################
# Version:    1.01
# Date:       18.01.2012
# Author:     A. Glenz
#
# Purpose:
# Generate random data set with specified subgroup structure.
#
# Called by:             command-line
# Calls:  
#
# Parameters:
#            ng:         Number of subgroups
#            range:      Upper limit for subgroup-centroids (lower limit is 0)
#            gsize:      Group Size, number of individuals in group 
#            nv:         Number of variables (attributes) for data set 
#            nvh:        Number of variables with low standard-deviation
#                        within subgroups (sd_hg). Must be lower or equal
#                        to nv.
#            sd_hg:      mean standard-deviation of homogemeus subgroups
#            sd_ihg:     mean standard-deviation for inhomogeneous subgroups
#                        (should be higher than sd_hg)
#
# Return value:          Data frame of size gsize x (nv), containing the nv attributes
#                        in columns. Vector containing the proposed
#                        subgroup partition as integer values
#
# Version 1.01:
# Separation of group association vector from data matrix
#################################
generate.random <- function(ng,range,gsize,nv,nhv,sd_hg,sd_ihg,teams=1)
{
   team.col <- 0; if (teams > 1) team.col <- 1  # add additional column for team-id, if there are 2 ore more teams
   data.all <- c()

   groups <- sort(as.vector(rep(1:ng,round((gsize/ng)+1)))[1:gsize])

   for (team in 1:teams)
   {
      data <- as.data.frame(matrix(ncol=nv,nrow=gsize))

      for (var in 1:nv)
      {
         gsample <- sample(ng,ng)
         for (subgroup in 1:ng)
         {
            if(nhv >= var)
            { 
      	        sdx <- sd_hg  
      	        meanx <- range/ng*(gsample[subgroup])
            } else 
            {
      	        sdx <- sd_ihg
      	        meanx <- 0
            }
            data[groups == subgroup,var] <- rnorm(sum(groups == subgroup),meanx,sdx)       
         }   
      }
      if (teams > 1)
      {
      	 data <- cbind(rep(team,gsize),data)
      	 names(data)[1] <- "teamid"
      }
      data.all <- rbind(data.all,data)
   }
   #return(data)
   list(data = data.all,
        groups = groups)
}




#################################
# function: cal_hara
#################################
# Version:    1.01
# Date:       18.01.2012
# Author:     A. Glenz
#
# Purpose:
# Calculate Calinski-Harabasz (1974) criterion for determining the
# optimal number of clusters
#
# Called by:             command-line
# Calls:  
#
# Parameters:
#            data:       Data set as data.frame with attributes in colums and 
#                        data-records in rows. 
#            groups:     Group-partition with integer values 
#                        representing the individua's assigned subgroup id
#
# Return value:          Calinski-Harabasz criterion
#
# Version 1.01:
# Separation of group association vector from data matrix
#################################
cal_hara <- function(data,groups)
{
   n <- length(groups)
   k <- length(unique(groups))
      
   # initialize attribute centroids matrix
   centroids <- matrix(nrow=max(groups),ncol=ncol(data))
   
   # caclulate cluster centroids
   for (group in unique(groups)) 
   {
   	  centroids[group,] <- sapply(data[groups == group,],mean) 
   }
   
   # calculate sum of within-cluster sum of squares
   w_k <- 0	  
   for (x in 1:n)
   {
      # sum of squared differences between cluster centroid and person x's measures
      qs.within <- sum((centroids[groups[x],] - as.numeric(data[x,]))^2) 
      w_k <- w_k + qs.within
   }   
   
   # calculate sum of between-cluster sum of squares
   # grand means = variable means over the whole sample
   means <- sapply(data,mean) 
   
   # sum of suared differences between clusters centroids and grand mean
   qs.between <- apply((t(t(centroids) - means))^2,1,sum)  
   
   # determine number of individuals in each cluster
   cluster.sizes <- describe.by(groups,groups,mat=TRUE)$n 
   
   # multiply qs.between with the number of individuals in each 
   # cluster and sum over all clusters 
   b_k <- sum(qs.between * cluster.sizes) 
   
   # now calculate calinski-harabasz index according to Calinski and Harabasz (1974)
   ch <- (b_k / (k - 1)) / (w_k / (n-k))
   
   return(ch)
}


#################################
# function: silhouette.width
#################################
# Version:    1.1
# Date:       21.07.2014
# Author:     A. Glenz
#
# Purpose:
# Calculate silhouette widths of a pre-clustered group according to
# Rousseeuw (1987)
#
# Called by:             asw_cluster
#                        fau_cluster.random
#                        fau_cluster.agglomerative
# Calls:  
#
# Parameters:
#            distmatrix: n x n distance matrix, containing the pairwise
#                        distances between all members of an n-sized group 
#            groups:     Vector of n integer values, representing the partitioning 
#                        of an n-sized group
#
# Return value:          Vector of length n, conatining the silhouette widths of
#                        all n members of the given group
#
# Version 1.1:
# added option to consider subgroup homogeneity
#################################

silhouette.width <- function(distmatrix,groups,individual="avg",usesghomo=FALSE,...)
{
 
   n <- length(groups)
   u <- sort(unique(groups))
   k <- length(u)

   
   
   # use squared distances
   #distmatrix <- distmatrix^2
   
   #cat("\n")
   s <- vector(length=n)
   for (i in 1:n)
   {
      s[i] <- -99
      if (sum(groups == groups[i]) > 1) # i is not an isolate
      {
      
      
      	a_i <- 0
      	b_i <- Inf
      	h <- c()
      	h.i <- c()
      	d.i <- c()
      	for (j in u)
      	{
          g <- which(groups == j)  # all group members of group j
          g <- g[g != i]  # all group members of group j, except of member i
          g.i <- c(g,i)  # all group members of group j, and member i
          n.j <- length(g)
          n.i <- n.j + 1
          
          d.i.j <- 0
          if (n.j > 0) d.i.j <- sum(distmatrix[g,i]) / n.j # average distance from i to all elements of j
          d.i <- c(d.i,d.i.j)
          
          h.j <- 0
          if (n.j > 1) h.j <- sum(distmatrix[g,g]) / (n.j * (n.j - 1)) # average distance between all elements of j, except i --> subgroup heterogentiy
          h <- c(h,h.j)
          
          h.i.j <- 0
          if (n.i > 1) h.i.j <- sum(distmatrix[g.i,g.i]) / (n.i * (n.i - 1)) # average distance between all elements of j, incl. i --> subgroup heterogenity
          h.i <- c(h.i,h.i.j)
      	}
      
      
      	dh.i <- (h.i - h) # dh.i is change in mean distances of each cluster, when i is included

      	a <- which(u == groups[i])  # which group (entry of group list u) does i belong to (ingroup)?
    	
      	d.i.inf <- d.i
      	d.i.inf[a] <- Inf
      	
      	b <- which.min(d.i.inf)  


        a.i <- d.i[a] + usesghomo * h[b] 
        b.i <- d.i[b] + usesghomo * h[a]
      	           

        #cat(a.i,b.i,"\n")
      	m <- max(a.i,b.i) # silhouette width denominator

      	if (m != 0) s[i] <- (b.i - a.i) / m
      	#if (i == n) print(paste(i,s[i],a.i,b.i))
      }
   }
   
   s[s == -99] <- 0


   if (k == 1) s <- (rep(-1,n))
   
   return(s)
}

silhouette.width.by.attr <- function(distmatrix.list,groups,individual="avg",usesghomo=FALSE,weights=rep(NA,length(distmatrix.list)),...)
{


   n <- length(groups)
   u <- sort(unique(groups))
   k <- length(u)
   
   weights <- as.vector(weights)
   
   if (is.na(weights)) weights <- rep(1,length(distmatrix.list))
   weights[is.na(weights)] <- 1
   
   #cat("\n")
   s <- matrix(nrow=n,ncol=length(distmatrix.list))
   for (i in 1:n)
   {
      s[i,] <- rep(-99,length(distmatrix.list))
      if (sum(groups == groups[i]) > 1) # i is not an isolate
      {
        a.i <- c()
        b.i <- c()
        for (attr in 1:length(distmatrix.list))
        { 
      

      	   h <- c()
      	   h.i <- c()
      	   d.i <- c()
      	   for (j in u)
      	   {
              g <- which(groups == j)  # all group members of group j
              g <- g[g != i]  # all group members of group j, except of member i
              g.i <- c(g,i)  # all group members of group j, and member i
              n.j <- length(g)
              n.i <- n.j + 1
          
              d.i.j <- 0
              if (n.j > 0) d.i.j <- sum(distmatrix.list[[attr]][g,i]) / n.j # average distance from i to all elements of j
              d.i <- c(d.i,d.i.j)
          
              h.j <- 0
              if (n.j > 1) h.j <- sum(distmatrix.list[[attr]][g,g]) / (n.j * (n.j - 1)) # average distance between all elements of j, except i --> subgroup heterogentiy
              h <- c(h,h.j)

      	   }
      
      

      	   a <- which(u == groups[i])  # which group (entry of group list u) does i belong to (ingroup)?
    	
      	   d.i.inf <- d.i
      	   d.i.inf[a] <- Inf
      	
      	   b <- which.min(d.i.inf)  


        	
           a.i[attr] <- d.i[a] + usesghomo * h[b] 
           b.i[attr] <- d.i[b] + usesghomo * h[a]      	
      	}           

        #cat(a.i,b.i,"\n")
      	#m <- max(a.i,b.i) # silhouette width denominator
      	m <- (a.i - b.i >= 0) * a.i + (a.i - b.i < 0) * b.i  # maximum of each element pair of a.i and b.i

        s[i,m != 0] <- (b.i[m != 0] - a.i[m!= 0]) / m[m != 0]
      	#if (i == n) print(paste(i,s[i],a.i,b.i))
      }
      
   }
   
   s[s == -99] <- 0
   s.x <<- s
   s <- apply(t(t(s) * weights),1,sum)
   s <- s / sum(weights)
    
   if (k == 1) s <- (rep(-1,n))
   
   return(s)
}



#################################
# function: categorize
#################################
# Version:    1.04
# Date:       01.03.2012
# Author:     A. Glenz
#
# Purpose
# Calculate categorzized version of data set with continuous attributes
#
# Called by:             command-line
# Calls:  
#
#
# Parameters:
#            data:       Data set as data.frame with continuous attributes in colums and 
#                        data-records in rows. 
#            boundaries: List of p vectors, where p denotes the number of attributes, containing 
#                        each attribute's values of boundaries for categorization. 
#                        "NA" must be provided for attributes that are already in a categorial form.  
#
# Return value:          vcategorized data-frame
#
# Version 1.01:
# Separation of group association vector from data matrix
#
# Version 1.02:
# Fixed a bug causing wrong categorizatino, if the number of categories
# differed between attributes
#
# Version 1.03:
# New parameter "ranges" allows data range to be specified for each
# separately
#
# Version 1.04:
# Specify categorization boundaries for each attribute instead of attribute ranges
#################################
categorize <- function(data,boundaries)
{

   if (length(boundaries) != ncol(data)) stop(paste("Boundaries were specified for ",length(boundaries)," attributes (but ",ncol(data)," are needed).",sep=""))
   
   data.c <- data
   for (att in 1:ncol(data))
   {
       if (!is.na(boundaries[[att]]))
       {
       	 for (boundary in 1:(length(boundaries[[att]])))
          {
              data.c[data[,att] >= sort(boundaries[[att]])[boundary],att] <- boundary + 1
              data.c[data[,att] < min(boundaries[[att]]),att] <- 1           
          }
      } 
   }
   return(data.c)
}



#################################
# function: fau.shaw
#################################
# Version:    1.0
# Date:       22.12.2011
# Author:     A. Glenz
#
# Purpose:
# Calculate Shaw's FLS according to Shaw (2004)
#
# Called by:             command-line
# Calls:  
#
# Parameters:
#            data.c:     Data set as data.frame with categorized attributes in colums and 
#                        data-records in rows. 
#            categories: Vector of length p, containing the number of categories for each of
#                        the p categorial attributes.
#
# Return value:          Faultline value representing Shaw's FLS.
#
#################################
fau.shaw <- function(data.c,categories=NULL)
{
  nv <- ncol(data.c)

  categories.calculated <- as.vector(sapply(lapply(data.c,unique),length))
  #cat(categories.calculated,"\n")

  if (is.null(categories)) categories <- categories.calculated

  if ((sum((categories - categories.calculated) < 0)) > 0)
  {
   cat.v <- which((categories - categories.calculated) < 0)[1]
   stop(paste("Variable",cat.v,"has",
               categories.calculated[cat.v],
              "or more categories."))
  }

  if (length(categories) != nv) stop("Length of categories-vector 
                                      should match the number of 
                                      variables in the dataset")

  #if (min(categories) < 2) stop("Each variable should have at least 2 categories")


  FLS <- 0
  # Calculate internal alignment (IA)
  for (var1 in 1:nv)
  {
   IA_var1 <- 0
   for (cat1 in unique(data.c[,var1]))	 
   {
      IA_var2 <- 0
      for (var2 in c(1:nv)[-var1])
      {
          o <- vector(length=categories[var2])  # o: vector for observed values
          e <- o	# e: vector for expected values
          p <- rep(0,categories[var2])			
          noalign <- p
          p[1] <- sum(data.c[,var1] == cat1)
          uni_var2 <- unique(data.c[,var2])
          for (i2 in 1:categories[var2])
          {
       	      ifelse(i2 <= length(uni_var2),cat2 <- uni_var2[i2],cat2 <- NA)
       	      o[i2] <- sum((data.c[,var1]==cat1) & (data.c[,var2]==cat2))
       	      e[i2] <- sum(data.c[,var1] == cat1) / categories[var2]
       	      
       	 }	
       	 o[is.na(o)] <- 0
       	 e[is.na(e)] <- 0
          o_var1 <- sum(o)
      
          for (x in 1:o_var1) noalign[which.min(noalign)] <- noalign[which.min(noalign)] + 1
          
          IA_obs <- sum((o-e)^2/e)
          IA_perfect <- sum((p-e)^2/e)
          IA_noalign <- sum((noalign-e)^2/e)
          MaxDiff <- IA_perfect - IA_noalign
          ifelse(MaxDiff != 0,
                 IA_var1_var2 <- ((IA_obs - IA_noalign)/MaxDiff), 
                 IA_var1_var2 <- 0)
          IA_var1 <- IA_var1 + IA_var1_var2        
      }
       
   }
   IA_var1 <- IA_var1 / (categories[var1] * (nv - 1))
   #cat("IA_var1_total:",IA_var1,"\n")




   # Calculate Cross-Group Alignment Index (CGAI)
   CGAI_var1 <- 0
   var2 <- var1+1
   while (var2 <= nv)
   {
       #cat("var1:",var1," var2:",var2,"\n")
       
       obs <- matrix(ncol=categories[var2],nrow=categories[var1])
       
       i1 <- 0
       for (cat1 in unique(data.c[,var1]))
       {
          i2 <- 0
          i1 <- i1 + 1
          for (cat2 in unique(data.c[,var2]))
          {
              i2 <- i2 + 1
              obs[i1,i2] <- sum((data.c[,var1]==cat1) & (data.c[,var2]==cat2))
              n_var1 <- apply(obs,2,sum)
              n_var2 <- apply(obs,1,sum)
              n_var1[is.na(n_var1)] <- 0
              n_var2[is.na(n_var2)] <- 0
          }
       }
       obs[is.na(obs)] <- 0
       
       if (length(unique(data.c[,var1])) > 1) 
       {
          l <- 0
          CP <- vector(length=((categories[var1]) * (categories[var1] - 1) / 2))
          Wt <- CP

          for (i in 1:(categories[var1]-1))
          {
        	
             for (j in (i+1):categories[var1])
             {
                l <- l+1

             
                for (x in 1:categories[var2])
                {
                   CP[l] <- CP[l] + (obs[i,x] * obs[j,x])
                }
                Wt[l] <- (n_var2[i] * n_var2[j])
             
                ifelse(Wt[l] != 0, CP[l] <- CP[l] / Wt[l], CP[l] <- 0)
             
                #cat("CP:",CP[l],"\n")
             }
          }
       
          W_denom <- sum(Wt)
          ifelse((W_denom != 0),Wt <- Wt / W_denom,Wt <- Wt * 0)
       
          CGAI_var1 <- CGAI_var1 + sum(CP * Wt)
          #cat("CGAI:",CGAI_var1,"\n")
       }
       var2 <- var2 + 1	
   }
   CGAI_var1 <- CGAI_var1 / (nv-1)	
   FLS_var1 <- IA_var1 * (1 - CGAI_var1)
   FLS <- FLS + FLS_var1
  }
  FLS <- FLS / nv


  return(FLS)
}



#################################
# function: fau.dawson
#################################
# Version:    1.05
# Date:       25.10.2012
# Author:     A. Glenz
#
# Purpose:
# Calculate Faultline Strengths according to v. Knippenberg and Dawson (2011)
#
# Called by:             command-line
#
# Parameters:
#            data:       Data set as data.frame with numerical attributes in colums and 
#                        data-records in rows. 
#
# Return value:          Van Knippenberg & Dawson's Faultline Strength measure.
#
# Version 1.01:
# Separation of group association vector from data matrix
#
# Version 1.02:
# Replaced multiple R-squared by numtiple correlation
#
# Version 1.03:
# Add multinomial logistic regression for nominal attributes.
# Add square root of Nagelkerke's Pseudo-R2 as replacement for
# multiple correlation, when criterion variable is nominal scaled.
# Use dummy-variables for nominal-scaled predictors.
#
# Version 1.04:
# Fixed bug when invoking multinomial regression with constant criterion
# 
# Versino 1.05:
# Set faultline strength to zero, if dataset contains constant attributes
#################################
fau.dawson <- function(data,scale=rep("numeric",ncol(data)))
{

   
   # convert nominal attributes to numbers
   for (i in which(scale == "nominal"))
   {
   	  data[,i] <- as.numeric(as.factor(data[,i]))
   }
     
     
   # result vector (for multiple correlation)
   mc <- c(1:ncol(data))   

   result = as.data.frame(matrix(ncol=ncol(data), nrow=ncol(data) + 1))
   rownames(result)=c("Criterion",
            paste("Predictor",1:(ncol(data) - 1),sep=""),
            "Multiple Correlation")
   names(result)=c(paste("Model",1:(ncol(data)),sep=""))

   variables <- c(1:ncol(data))

   
   # sample size
   n <- nrow(data)
   
   # loop over all attributes in data set
   for (i in variables)   
   {
      av <- i # criterion	      				
      uv <- c(1:ncol(data))[-i] # predictor(s)
   
      # write criterion-name into result-matrix
      result[1,i] <- names(data)[i]  
   
      # Build new dataset wit dummycoded nominal variables (if we have them)

      # write criterion variable
      data.d <- data[av]
      d.col <- 1
      names(data.d)[d.col] <- names(data)[av]
      
      	
      # add predictors
      for (uvar in uv)
      {      	    
      	    if (scale[uvar] == "numeric")
            {
            	d.col <- d.col + 1
            	data.d <- cbind(data.d,data[,uvar])
            	names(data.d)[d.col] <- names(data)[uvar] 
            }	else # attribute is nominal-scaled
            {
                cats <- sort(unique(data[,uvar]))
                
                # dummycoding
                dummies <- c()
                for (category in cats) # we use a dummy variable for each of the variable's categories (while knowing that we could have omitted one of them)
                {
                   dummyvar <- rep(0,nrow(data))
                   dummyvar[data[,uvar] == category] <- 1 
                   dummies <- cbind(dummies,dummyvar)
                }
                d.col <- d.col + ncol(dummies)
                dummies <- as.data.frame(dummies)
                names(dummies) <- paste(names(data)[uvar],cats,sep=".")
            	data.d <- cbind(data.d,dummies)           	
            }     	    	
      }
   
   
   
   
      # Build formula for regression
      
      regr.formula <- paste("as.numeric(",names(data.d)[1],") ~ ", sep="")  
      regr.formula.null <- paste(regr.formula, "1")

      
      # build up regression formula
      for (j in 2:ncol(data.d))    # loop over all predictors  
      {   	    
   	    # no "+" for the first predictor
   	    ifelse((j == 2), pred.string <- "", pred.string <- " + ")   
   	    pred.string <- paste(pred.string,"as.numeric(",names(data.d)[j],")", sep="")
   	    regr.formula <- paste(regr.formula, pred.string)
      }
      
      # set names in result matrix
      loop <- 1
      for (j in variables[-i])
      {
      	loop <- loop + 1
      	result[loop,i] <- names(data[j])
      }
      
      # Multiple correlation from regression-Output
      # for use with lm():   mc[i] <- sqrt(summary(glm(regr.formula, data = data))$r.squared)    
      # The following is for use with glm()
      
            
      # Use multiple linear regression for numeric dependent variables
      if (length(unique(data[,av])) == 1)  # Criterion variable is constant
      {
      	 mc[i] <- 0 # According to mail-message from J. Dawson (October 24, 2012)
      } else if (all(apply(cbind(data[,uv]),2,var) == 0)) # all predictors are constant, while criterion has non-zero variance
      {
      	 mc[i] <- 0 
      }
      else if (scale[av] == "numeric") 
      {
      	 mod <- glm(regr.formula, data = data.d)
      	 mc[i] <- cor(data.d[,1], predict(mod))
      } else # Use multinomial logistic regression for nominal scaled dependent variables if (scale[av] == "nominal") 
      {

      	 attach(data.d)
      	 mod <- multinom(regr.formula)
      	 mod.null <- multinom(regr.formula.null)
      	 detach(data.d)
      	 
      	 #### calculate Nagelkerke's Pseudo R-square
      	 # First extract Log Likelihood from null- and full models
      	 ll.null <- logLik(mod.null)[1]
         ll.mod <- logLik(mod)[1]

         RN.numerator <- 1 - exp(-2/n * (ll.mod - ll.null)) # entspricht Cox & Snell's Pseudo-R2
         RN.denominator <- 1 - exp(2/n * ll.null)
         mc[i] <- sqrt(RN.numerator / RN.denominator)
         
      }
      
   }

   result[ncol(data) + 1,] <- mc

   x <- list()
   x$result <- result
   x$fau <- prod(mc)
   return(x)
}



#################################
# function: fau.trezzini
#################################
# Version:    1.0
# Date:       22.12.2011
# Author:     A. Glenz
#
# Purpose:
# Calculate Trezzinis's PMD_cat according to Trezzini (2008)
#
# Called by:             command-line
# Calls:  
#
# Parameters:
#            data.c:     Data set as data.frame with categorized attributes in colums and 
#                        data-records in rows. 
#            categories: Vector of length p, containing the number of categories for each of
#                        the p categorial attributes.
#
# Return value:          Faultline value representing Shaw's FLS.
#
#################################
fau.trezzini <- function(data.c,categories=NULL)
{
  nv <- ncol(data.c)  # number of attributes
  l <- nrow(data.c)   # length of data set (corresponds to group size)

  # determine categories from data set
  categories.calculated <- as.vector(sapply(lapply(data.c,unique),length))

  # use calculated categories if none were supplied as a parameter of this function
  if (is.null(categories)) categories <- categories.calculated

  # plausibility check
  if ((sum((categories - categories.calculated) < 0)) > 0)
  {
   cat.v <- which((categories - categories.calculated) < 0)[1]
   stop(paste("Variable",cat.v,"has",
              categories.calculated[cat.v],
              "or more categories."))
  }

  if (length(categories) != nv) stop("Length of categories-vector 
                                     should match the number of 
                                     variables in the dataset")
                                     
                            
   # convert to numeric
   for (i in 1:nv)
   {
      data.c[,i] <- as.numeric(data.c[,i])
   }
 
   # generate list of attribute occurrencies
   c_list <- list()
   for (p in 1:nv)
   {
       c_list[[p]] <- unique(data.c[,p])
   }
   
   # generate matrix of attribute combinations
   c_matrix <- as.matrix(expand.grid(c_list)) 
   
   # count number of different attribute combinations
   n <- nrow(c_matrix)
   
   # determine relative prevalences of attribute_combinations in data.c
   c_occ <- vector(length=n)
   for (c_row in 1:n)
   {
       c_occ[c_row] <- sum(apply(abs(rep(as.matrix(c_matrix[c_row,]),l) 
                           - t(data.c)),2,sum) == 0)    
   }
   c_occ <- c_occ / l
   
   # determine PMD measure
   pmd <- 0
   for (i in 1:n)
   {
     for (j in 1:n)
       {
           d_ij <- sum(c_matrix[i,] != c_matrix[j,]) / nv
           pmd <- pmd + ((c_occ[i] + c_occ[j]) * c_occ[i] * c_occ[j] * d_ij)
       
       }
   }
   return(pmd)
}


#################################
# function: fau.gibson
#################################
# Version:    1.01
# Date:       9.02.2012
# Author:     A. Glenz
#
# Purpose:
# Calculate Subgroup Strength according to Gibson's (2003) attribute Overlap index.
#
# Called by:             command-line
#
# Parameters:
#            data:       Data set as data.frame with numerical attributes in colums and 
#                        data-records in rows. 
#            scale:      Vector of length p, where p is the number of attributes in
#                        the data set, containing strings that denote the type of attributes
#                        (either "nominal" or "numeric")
#
# Return value:          Gibson's Subgroup Strength measure.
#
# Version 1.01:
# Added parameter 'scale' to distinguish ratio-scaled attributes from
# categorial attributes
#################################
# Gibson's Subgroup Strength
fau.gibson <- function(data,scale=rep("numeric",ncol(data)))
{
   n <- nrow(data)
   nv <- ncol(data)
   p <- (n*(n-1)/2)
   o <- matrix(nrow=p, ncol=nv)
   
   
   for (k in 1:nv)
   {
      index <- 0
      for (i in 1:(n-1))
      {
         for (j in (i+1):n)
         {
             index <- index + 1
             if (scale[k] == "nominal")
             {
                o[index,k] <- sum(data[i,k] == data[j,k])
             } else
             {
                o[index,k] <- min(data[i,k],data[j,k]) / max(data[i,k],data[j,k]) 
             }
         }
      }
   }
   overlap <- apply(o,1,sum)
   
   overlap.std <- sd(overlap)
   
   return(overlap.std)
   
}

#################################
# function: fau.lcca
#################################
# Version:    1.0
# Date:       26.19.2012
# Author:     A. Glenz
#
# Purpose:
# Calculate faultlines using latent class cluster analysis provided by the "flexmix" package
#
# Called by:             command-line
#
# Parameters:
#            data:       Data set as data.frame with numerical attributes in colums and 
#                        data-records in rows. 
#
# Return value:          Vector containing 
#                        - Faultline Strength (mean posterior sg-membership probability)
#                        - Individual posterior sg-membership probability
#                        - Number of Subgroups
#                        - Member-to-Subgroup association
#
#################################
# LCCA
fau.lcca <- function(data,maxgroups=6)
{
   # set number of repetitions with random start configurations
   nreps <- 10
   
   
   n <- nrow(data)
   cont <- list(minprior = 1 / n, verbose = 0)
   model <- tryCatch(stepFlexmix(as.matrix(data) ~ 1, k = 1:maxgroups, nrep = nreps, model = FLXmclust(diagonal=FALSE), control = cont),error=function(e) NULL)
   #model <- tryCatch(stepFlexmix(as.matrix(data) ~ 1, k = 1:maxgroups, nrep = nreps, model = FLXMCmvnorm(diagonal=FALSE), control = cont),error=function(e) NULL)



   best.model <- getModel(model, which="BIC")

   p <- posterior(best.model)
   groups <- reord.groups(apply(p,1,which.max))
   i.post <- apply(p,1,max)

   pe <- p * log2(p)
   pe[is.na(pe)] <- 0
   entropy <- 1 - mean(-1 * apply(pe,1,sum),na.rm=TRUE) # Shannon's Entropy Measure   
   
   
   cat("Entropy:",entropy,"   EIC:",EIC(best.model),"\n")
   
   list(fs = entropy,
         i.post = i.post,
         groups = groups,
         num_sg = length(unique(groups)))
}









#################################
# function: fau.thatcher
#################################
# Version:    1.01
# Date:       4.10.2012
# Author:     A. Glenz
#
# Purpose:
# Calculate Subgroup Strength a suggested by Thatcher et al. (2003)
#
# Called by:             command-line
#
# Parameters:
#            data:       Data set as data.frame with numerical attributes in colums and 
#                        data-records in rows. 
#
# Return value:          Thatcher's Faultline Strength measure.
#
# Version 1.01:
# Bug fixed, which aused function to crash, if there were no variance in all attributes
#################################
fau.thatcher <- function(data,weights = NA,metric="euclid",dummycols=NA)
{
   n <- nrow(data)
   nv <- ncol(data)
   
   if (is.na(sum(weights)))
   {
   	  weights <- rep(1,nv)
   }
  
   data <- as.data.frame(t(t(data) * weights))
  
   if (metric == 'mahal')
   {
       S <- S.mahal(data,dummycols)
       fau.denom <- fau.denominator.mahal(data,sapply(data,mean),S)
   } else
   {
       if (metric == 'corr')
       {
          C <- C.cor(data) * diag(ncol(data))
          data <- as.matrix(data) %*% C
       }
       fau.denom <- fau.denominator(data)
   }
     
     
   if (fau.denom != 0)
   {  
  
      # number of partitions  
      npart <- 2^(n - 1) - 1
   
      data.mean <- sapply(data,mean)
   
      max.fau <- 0
      for (partition in 1:npart)
      {
          p <- partition
          partition.bin <- c()
          x <- 0
          while (p >= 1)
          {
              # convert to binary value (base-2)
              x <- x + 1
              p <- p / 2
              ifelse(p == as.integer(p), char <- "0", char <- "1")
              p <- as.integer(p)
              partition.bin <- c(char,partition.bin)
          }
          groups <- c(rep("0",n-length(partition.bin)),partition.bin)
          groups <- as.numeric(groups) + 1
       
          if (metric == 'mahal')
          {
             fau <- fau.numerator.mahal(data,groups,data.mean,S) / fau.denom
          } else
          {
             fau <- fau.numerator(data,groups,data.mean) / fau.denom
          }
          if (fau > max.fau)
          {
              max.fau <- fau
              max.groups <- groups
          }      
      }
   } else
   {
   	 max.fau <- 0
   	 max.groups <- rep(1,n)
   }
   return(list(fau = max.fau,
               groups = max.groups))
}



#####################################
# function: asw_cluster.agglomerative
#####################################
# Version:    1.031
# Date:       25.07.2013
# Author:     A. Glenz
#
# Purpose:
# Invoke asw-maximizing clustering procedure starting at per-clustered solutions from ward's and 
# average-linkage agglomerative methods. 
#
# Called by:             command-line
# Calls:                 S.mahal
#                        distmatrix.e
#                        distmatrix.m
#                        C.cor
#                        reord.groups
#                        asw_cluster
#                        silhouette.width
#
# Parameters:
#            data:       Data set as data.frame with attributes in colums and data-records in rows. 
#            maxgroups:  Maximum number of subgroups allowed in partitions
#            metric:     if metric = "mahal", partitions will be evaluated based on Mahalanobis-distances
#                        if metric = "corr", data rows (attributes) will be weighted with 1 - cm, where 
#                        cm denotes the column mean of the pxp correlation matrix, and p denotes the
#                        number of attributes in the data set 
#                        For all other values of metric, including NA, partitions will be evaluated 
#                        based on Euclidean distances
#            weights:    Vector of length p, containing weighting factors for the p attributes
#            dummycols:
#            quiet:
#            ilevel:     Logical, specifying if individual silhouette-widths
#                        will be maximized
#
# Returned value:        Data vector containing the maximum value for asw in the first place, followed by 
#                        the subgroup partition (n values) for this value.
#
# Version 1.01:
# Separation of group association vector from data matrix
#
# Version 1.02:
# Output contains individual silhouette widths ($sil)
# 
# Version 1.03:
# New logical parameter 'i.level' for calculating individual faultlines
#
# Version 1.031:
# Return a asw/sil value of 0 (instead of NA), if all attributes have no variance for a given team
#####################################
asw_cluster.agglomerative <- function(data,i.level = FALSE,maxgroups=6,metric='euclid',weights=NA,dummycols=NA,quiet=FALSE,usesghomo=FALSE,by.attr=FALSE)
{
    n <- nrow(data)
    groups <- c(1:n)
     
    methods <- c("constant","ward","average_linkage")

    if (metric == 'mahal')
    {
       S <- S.mahal(data,dummycols)
       D <- distmatrix.m(data,S,weights)
       D.S <- distmatrix.m(scale(data),S,weights)
    } else
    {
       if (metric == 'corr')
       {
          C <- C.cor(data) * diag(ncol(data))
          data <- as.matrix(data) %*% C
       } else 
       {
       	  if (metric == 'gower')
          {
             D.W <- distmatrix.w(data,weights)
             D <- distmatrix.e(data,weights)
         	if (by.attr == TRUE) D.L <- distmatrix.g.list(data)
          } else # metric = euclid
          {    
             if (by.attr == TRUE) D.L <- distmatrix.e.list(data)
             
             # standardize by range
             D.W <- distmatrix.w(data,weights)
             D <- distmatrix.e(data,weights)


          }
       }
    }
    
    ifelse(by.attr == TRUE,dst <- D.L,dst <- D)  
    
    if (sum(D) > 0)
    {
       sil <- vector(length=n,mode='numeric')
       num_sg <- vector(length=n,mode='numeric')
       mbr <- vector(length=n,mode='numeric')
       adjm.all <- matrix(0,n,n)
       
       aswmax <- -1
       subgroups <- matrix(ncol=n,nrow=0)     
     
      if (i.level == FALSE)
      {
    
       # Combine Clusters using ward's algorithm
       cl1 <- 0
       cl2 <- 1
    
       ward.x <- as.matrix(c(1:n))
    
       while ((length(unique(groups)) > 1) & (sum(D) > 0) & (cl1 != cl2))
       {
         
            # initialize attribute centroids matrix
            nc <- max(groups)#
            centroids <- matrix(nrow=nc,ncol=ncol(data))

         
            # caclulate cluster centroids
            for (group in unique(groups)) 
            {
         	       centroids[group,] <- colMeans(data[groups == group,])
         	   }
         
            # calculate cluster sizes
            cl.sizes <- as.vector(tapply(groups,groups,length))
         
            # calculate distance matrix from centroids
            if (metric == 'mahal')
            {
               S <- S.mahal(data,dummycols)
               dist.cent <- distmatrix.m(centroids,S,weights)
            } else
            {
               if (metric == 'corr')
               {
                  C <- C.cor(centroids) * diag(ncol(centroids))
                  centroids <- as.matrix(centroids) %*% C
               } else
               {
               	  if (metric == 'gower')
               	  {
               	  	  dist.cent <- distmatrix.g(centroids)
               	  } else
               	  dist.cent <- distmatrix.e(centroids,weights)
               }
               
            }
         
            cl.mat <- matrix(ncol=length(cl.sizes),nrow=length(cl.sizes))
            for (i in 1:length(cl.sizes))
            {
         	      for (j in 1:length(cl.sizes))
         	      {
         	         cl.mat[i,j] <- cl.sizes[i] + cl.sizes[j]
         	      }
         	   }
            #dist.cent <- dist.cent * ((cl.sizes %*% t(cl.sizes))  / (cl.sizes %*% t(cl.sizes)))
            dist.cent <- dist.cent * ((cl.sizes %*% t(cl.sizes))  / cl.mat)
         
            dist.cent[is.na(dist.cent)] <- max(dist.cent) + 1 

            # replace diagonal elements with max. distance
            diag(dist.cent) <- max(dist.cent) + 1
         
            # find minimal distance
            dcmin <- which.min(dist.cent)
         
            # find centroids to combine
            cl1 <- ceiling(dcmin / nrow(dist.cent))
            cl2 <- dcmin - ((cl1 - 1) * nrow(dist.cent))
         
            # combine clusters
            groups[groups == cl2] <- cl1 
            groups <- reord.groups(groups)
            ward.x <- cbind(groups,ward.x) 
         
  
       } 
           
    
       # Combine clusters with average linkage method
       groups <- c(1:n)
       cl1 <- 0
       cl2 <- 1
       average.link <- as.matrix(groups)
       while ((length(unique(groups)) > 1) & (sum(D) > 0) & (cl1 != cl2))
       {
          min.avg.link <- max(D)
          for (i in unique(groups))
          {
             for (j in unique(groups)[-i])
             {
                d <- mean(as.vector(as.matrix(D[groups == i, groups == j])))
                if (d <= min.avg.link)
                {
                   min.avg.link <- d
                   cl1 <- i
                   cl2 <- j
                }
             }
          }
          # combine clusters
       groups[groups == cl2] <- cl1 
       groups <- reord.groups(groups)
       average.link <- cbind(groups,average.link)
 

       }
    
       # now call clustering with asw maximizing 
       n.w <- min(maxgroups,ncol(ward.x),ncol(average.link))

       a.l <- as.data.frame(as.matrix(average.link[,2:n.w])
           [,apply(as.matrix(abs(average.link[,2:n.w] - ward.x[,2:n.w])),2,sum) > 0])
       m <- cbind(rep(1,n),ward.x[,2:n.w],a.l)
       ifelse(ncol(a.l) > 0, a.l_names <- rep("average_linkage",ncol(a.l)), a.l_names <- "")
    
       methods <- c("constant",rep("ward",n.w - 1),a.l_names)
    
       asw <- vector(length=n.w,mode='numeric')
       sil <- vector(length=n,mode='numeric')
       num_sg <- vector(length=n,mode='numeric')
       mbr <- vector(length=n,mode='numeric')
       
       aswmax <- -1
       subgroups <- matrix(ncol=n,nrow=0)

       i.steps <- "avg"
       if (i.level == TRUE)
       {
          i.steps <- c(1:n)
       }
       
       adjm.all <- matrix(0,n,n)
       asw.all <- vector(length=n.w,mode='numeric')
       
       for (individual in i.steps)
       {
       aswmax <- -1
       groupmax <- c(1)
       
    
       for (i in 1:ncol(m))
       {
          groups <- m[,i]

          ifelse(n > 1,res.asw <- asw_cluster(data,groups,maxgroups,dst,individual,usesghomo=usesghomo,weights=weights,by.attr=by.attr),res.asw <- c(-1,1))
          g <- reord.groups(res.asw[-1])
          asw[i] <- round(res.asw[1],digits=8)

          groups.x <- g
          if (asw[i] >= aswmax)
          {

             if (asw[i] == aswmax)
             {
              
                 if ((length(unique(g)) < length(unique(groupmax))))
                 {
                    aswmax <- asw[i]
                    groupmax <- g
                    ng <- length(unique(groupmax))
                    max_i <- i
                    max_method <- methods[i]
                 } else
                 {
                    if ((length(unique(g)) == length(unique(groupmax))))
                    {
                       max_method <- unique(c(max_method,methods[i]))
                    }
                 }
             } else
             {
                aswmax <- asw[i]
                groupmax <- g
                ng <- length(unique(groupmax))
                max_i <- i
                max_method <- methods[i]
             }
          }
       
          method <- methods[i]
          if (quiet == FALSE) cat("Attempt",i," Method:",method," Asw: ",asw[i],"  Groups:",m[,i]," -->",g,"\n")     
       }

          
          if (quiet == FALSE) cat("Groups:",groupmax," ASW:",aswmax," #groups:",ng," method:",max_method,"\n")
          if (individual == "avg")
          {
          	 if (by.attr == TRUE) sil <- (silhouette.width.by.attr(D.L,groupmax,individual=individual,usesghomo=usesghomo,weights=weights))
          	 else sil <- (silhouette.width(D,groupmax,individual=individual,usesghomo=usesghomo))
          	 num_sg <- length(unique(groupmax))
             mbr <- individual
             subgroups <- groupmax
             adjm.all <- NA
          }  else
          {
          	 sil[individual] <- aswmax
          	 num_sg[individual] <- length(unique(groupmax))
          	 mbr[individual] <- individual
          	 adjm.individual <- (groupmax %*% t(groupmax) == groupmax^2) * 1
          	 subgroups <- rbind(subgroups,groupmax)
             adjm.all <- adjm.all + adjm.individual
          }
          
       } 
        
       adjm.all <- adjm.all / n

      } else # individual mode
      {
      	gmat <- ingroups.ind(D,D.W,D.L,data,usesghomo=usesghomo,weights=weights,by.attr=by.attr)
      	sil <- c()
      	if (n > 2)
      	{
      		for (focal in 1:n)
      	{
       		 groupmax <- gmat[focal,]
       		 if (by.attr == TRUE) sil.focal <- (silhouette.width.by.attr(D.L,groupmax,usesghomo=usesghomo,weights=weights))[focal]
       		 else sil.focal <- (silhouette.width(D,groupmax,usesghomo=usesghomo))[focal]
       		 
       		 sil <- c(sil,sil.focal)

          	 num_sg[focal] <- length(unique(groupmax))
          	 mbr[focal] <- focal
          	 adjm.individual <- (groupmax %*% t(groupmax) == groupmax^2) * 1
             adjm.all <- adjm.all + adjm.individual
      	}
      	adjm.all <- adjm.all / n
      	aswmax <- 0
      	subgroups <- gmat
      	} else
      	{
            aswmax <- 0
            sil <- rep(0,n)
            subgroups <- rep(1,n)
            num_sg <- 1
            mbr <- 1
            adjm.all <- matrix(0,n,n)
         
        }

     }
    
    } else # all persons are identical in all attributes
    {
        aswmax <- 0
        sil <- rep(0,n)
        if (i.level == TRUE) subgroups <- matrix(1,n,n) else subgroups <- rep(1,n)
        
        num_sg <- 1
        mbr <- 1
        adjm.all <- matrix(0,n,n)
    }
    #return(c(aswmax,groupmax))
    list(asw = mean(sil),
         sil = sil,
         groups = subgroups,
         num_sg = num_sg,
         mbr = mbr,
         adjm = adjm.all)
}







#################################
# function: partitions
#################################
# Version:    1.0
# Date:       22.12.2011
# Author:     A. Glenz
#
# Purpose:
# Determines all possible partitions of a group with n members
# and a given maximum number of subgroups
#
# Called by:             command-line
# Calls:                 partitions (recursive)
#
# Parameters:
#            n:          Group size
#            n.groups:   Maximum number of subgroups for partitions
#            x:          start position for recursion (do not use for direct calls)
#            s:          initial partition for recursion (do not use for direct calls)
#            p:          partition result table (do not use for direct calls)
#
# Return value:          Table of partitions
#################################
partitions <- function(n,n.groups=n,x=2,s=rep(1,n),p=c())
{   
    for (i in 1:(min(max(s[c(1:x-1)]) + 1,n.groups)))
    {
         s[x] <- i
         ifelse (x == n,p <- rbind(p,s),p <- partitions(n,n.groups,x+1,s,p)) 
    }
   return(p)
}





#################################
# function: convert.to.long
#################################
# Version:    1.1
# Date:       6.02.2013
# Author:     A. Glenz
#
# Purpose:
# Converts the output of the "faultlines"-function into a 'long'-format, suitable
# for multilevel-analysis
#
# Called by:             command-line
# Calls:                 
#
# Parameters:
#          res:          result dataframe from 'faultlines'-function
#
# Return value:          result table (data-frame) in long-format
#
# Version 1.1
# String-type group-ids are now accepted
#################################
convert.to.long <- function(res, ...)
{
   pars <- 3
   mts <- mbr <- list()
   n.teams <- length(res) - pars
   res.out <- c()
   for (i in 1:n.teams)
   {
   	   team <- res[[i]]$team
   	   mbr <- res[[i]]$fl.mbr
   	   if (is.null(res[[i]]$mbr_to_subgroups)) res[[i]]$mbr_to_subgroups <- NA
   	   for (mb in 1:length(mbr))
   	   {
   	   	   mbr_sg <- matrix(res[[i]]$mbr_to_subgroups,byrow=FALSE,ncol=res[[i]]$teamsize)[mb,]
   	   	   for (sg in mbr_sg)
   	   	   {
   	   	   	   
   	   	   	   if (res$i.level == TRUE)
   	   	   	   {
   	   	   	      res.out <- rbind(res.out,
   	   	   	                    c(res[[i]]$team,
   	   	   	                      res[[i]]$teamsize,
   	   	   	                      res[[i]]$fl.ind[mb],
   	   	   	                      res[[i]]$fl.mbr[mb],
   	   	   	                      sg,
   	   	   	                      res[[i]]$number_of_subgroups[mb],
   	   	   	                      sum(sg == mbr_sg),
   	   	   	                      sum(mbr_sg[mb] == mbr_sg)))
   	   	   	   	
   	   	   	   } else
   	   	   	   {
   	   	   	      res.out <- rbind(res.out,
   	   	   	                    c(res[[i]]$team,
   	   	   	                      res[[i]]$teamsize,
   	   	   	                      res[[i]]$fl.value[mb],
   	   	   	                      res[[i]]$fl.mbr[mb],
   	   	   	                      sg,
   	   	   	                      res[[i]]$number_of_subgroups[mb],
   	   	   	                      sum(sg == mbr_sg),
   	   	   	                      NA))
   	   	   	   }                   
   	   	   }
   	   }
   }
   
   
   res.out <- as.data.frame(res.out)
   names(res.out) <- c(names(res[[1]])[1:(ncol(res.out) - 2)], "subgroup_size","own_subgroup_size")
   
   res.out$teamsize <- as.numeric(as.character(res.out$teamsize))
   res.out$fl.value <- as.numeric(as.character(res.out$fl.value))
   res.out$mbr_to_subgroups <- as.numeric(as.character(res.out$mbr_to_subgroups))
   res.out$number_of_subgroups <- as.numeric(as.character(res.out$number_of_subgroups))
   res.out$subgroup_size <- as.numeric(as.character(res.out$subgroup_size))
   res.out$own_subgroup_size <- as.numeric(as.character(res.out$own_subgroup_size))
   return(res.out)
}

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

faultlines <- function(data,group.par="NA",attr.type=NA,attr.weight=NA,rescale=NA,method="asw",metric="euclid",maxgroups=6, i.level=FALSE, usesghomo=FALSE,by.attr=FALSE,cores=1, ...) UseMethod("faultlines")


 
faultlines.calc <- function(group, group.vect , data, method=method, maxgroups=maxgroups, metric=metric, attr.weight.team=attr.weight.team, attr.type = attr.type, dummycols.team=dummycols.team, dummylist, quiet=TRUE,i.level=i.level,usesghomo=FALSE,by.attr=FALSE,cores=1) 
{
        data.team <- data[group.vect == group,]
        data.team.X <<- data.team
        cat ("Group:",group,"  Groupsize:",nrow(data.team),"\n")
      
       if ((method %in% c("THATCHER","ASW","BEZRUKOVA")) & (toupper(metric) == "MAHAL" ))
       {
          zerovarcols <- (1:ncol(data.team))[apply(data.team,2,var) == 0]
          if (length(zerovarcols) > 0)
          {
                 dummycols.team <- dummycols.team[-(which(dummycols.team %in% zerovarcols))]
                 data.team <- data.team[,-zerovarcols]
                 attr.weight.team <- attr.weight.team[-zerovarcols]
                 for (col in sort(zerovarcols, decreasing=TRUE)) dummycols.team[dummycols.team > col] <- dummycols.team[dummycols.team > col] - 1
          }
       }
       
       data.team <- as.data.frame(data.team)
       
       if (method %in% c("THATCHER","ASW","BEZRUKOVA")) 
       {
       	   cols.none <- as.vector(which(apply(data.team,2,function(x) all(x==0))))  # colnubers of variables that have 0 in all rows
       	   dummycols.none <- cols.none[cols.none %in% dummycols.team] # dummyvariables that do not occur in this team
       	   #dummycols.none <- as.vector(which(apply(as.data.frame(data.team[,dummycols.team]),2,function(x) all(x==0)))) # dummyvariables that do not occur in this team
       	   if (length(dummycols.none) > 0)
       	   {
       	      nodummycols <- !(1:ncol(data.team) %in% dummycols.team)[-cols.none]
       	      data.team <- data.team[,-dummycols.none]
       	      for (dcol in rev(dummycols.none))
       	      {
       	        dummycols.team[which(dummycols.team==dcol)] <- 0
       	        dummycols.team[dummycols.team > dcol] <- dummycols.team[dummycols.team > dcol] - 1
       	        dummycols.team <- dummycols.team[dummycols.team > 0]
       	        
       	        attr <- which(unlist(lapply(dummylist,function(x) dcol %in% x)))
       	        dummylist[[attr]] <- dummylist[[attr]][-which(dummylist[[attr]]==dcol)]
       	        
       	      }
       	   }  
       	      
       	   if (by.attr)
       	   {
       	        for (attr in 1:length(dummylist))
       	        {
       	           attr.weight.team[dummylist[[attr]]] <- attr.weight.team[dummylist[[attr]]] / length(attr.weight.team[dummylist[[attr]]])
       	        }
       	   }
       	          
       	   if (length(dummycols.none) > 0) attr.weight.team <- attr.weight.team[-dummycols.none]
       	   # dummycols.team <- dummycols.team[!dummycols %in% dummycols.none]
       	   
       }	   
       

       if ((is.null(dummycols.team)) | (length(dummycols.team) == 0)) dummycols.team <- NA
       
       
       dummycols.X <<- dummycols.team
       data.X <<- data.team
       attr.weight.team.X <<- attr.weight.team
       dummylist.X <<- dummylist

       #cat(dummycols,"\n")
       if (ncol(data.team) > 1)
       {
          fl.mbr <- NA
          fl.adjm <- NA
          fl.ind <- NA    
          fl.groups <- NA      
          # choose function according to specified method
          if (method == "ASW")
          {
              fl.asw <- asw_cluster.agglomerative(data.team, maxgroups=maxgroups, metric=metric, weights=attr.weight.team, dummycols=dummycols.team, quiet=TRUE,i.level=i.level,usesghomo=usesghomo,by.attr=by.attr)

              fl <- fl.asw$asw
              fl.ind <- fl.asw$sil
              fl.groups <- fl.asw$groups 
              fl.num_sg <- fl.asw$num_sg
              fl.mbr <- fl.asw$mbr
              fl.adjm <- fl.asw$adjm
          }
       
          if (method == "THATCHER")
          {
              fl.thatcher <- fau.thatcher(data.team, weights=attr.weight.team, metric=metric, dummycols=dummycols.team)
              fl <- fl.thatcher$fau
              fl.groups <- fl.thatcher$groups
              fl.num_sg <- length(unique(fl.groups))
              fl.mbr <- "avg"
          }

          if (method == "LCCA")
          {
              options(show.error.messages=FALSE)
              fl.lcca <- tryCatch(fau.lcca(data.team,maxgroups=maxgroups),error=function(e) NULL)
              options(show.error.messages=TRUE)
              fl <- fl.lcca$fs
              fl.ind <- fl.lcca$i.post
              fl.groups <- fl.lcca$groups
              fl.num_sg <- fl.lcca$num_sg            
          }
      
          if (method == "SHAW")
          {
              fl <- fau.shaw(data.team)
          }
       
          if (method == "TREZZINI")
          {
              fl <- fau.trezzini(data.team)
          }
       
          if (method == "GIBSON")
          {
           
              fl <- fau.gibson(data.team,scale=attr.type)
          }
       
          if (method == "KNIPPENBERG")
          {
              fl <- fau.dawson(data.team,scale=attr.type)$fau
          }
       
          if (method == "BEZRUKOVA")
          {
              fl.thatcher <- fau.thatcher(data.team, weights=attr.weight.team, metric=metric, dummycols=dummycols.team)
              ifelse(metric == "mahal",
                 fl <- fl.thatcher$fau * fau.dist.mahal(data.team,fl.thatcher$groups,S.mahal(data.team,dummycols.team)),
                 fl <- fl.thatcher$fau * fau.dist(data.team,fl.thatcher$groups)) 
              fl.groups <- fl.thatcher$groups
              fl.num_sg <- length(unique(fl.groups))
              fl.mbr <- "avg"
          }
          
                
       
       } else # none of the data columns had a variance > 0
       {
            fl <- NA
            fl.groups <- NA
       }
       
       ifelse(is.na(fl.groups) || is.null(fl.groups),fl.nofsg <- NA, fl.nofsg <- fl.num_sg)
       
       if (metric == "mahal") 
       {
           	S <- S.mahal(data.team, dummycols.team)
           	D <- distmatrix.m(data.team, S)
       } else 
       {
           	D <- distmatrix.e(data.team)
       }
       
       
       return(list(team=group,teamsize=nrow(data.team),fl.value=fl,fl.mbr=fl.mbr,mbr_to_subgroups=fl.groups,number_of_subgroups=fl.nofsg,adjm=fl.adjm,fl.ind=fl.ind,distmat=D))
}


#################################
# function: faultlines
#################################
# Version:    1.05
# Date:       22.01.2013
# Author:     A. Glenz
#
# Purpose:
# Wrapper function for various faultline calculations
#
# Version 1.01
# Fixed syntax error
# Added error handling for datasets that have no variance in all attributes
# 
# Version 1.02
# Add functionality to calculate individual-leve faultlines
#
# Version 1.03
# Add standardisation feature to rescale numeric attributes
#
# Version 1.04
# Added method "lcca"
# Changed name of result variable "sil" to "fl.ind" (individual faultline value)
#
# Version 1.05
# Handling of situation when all nominal attributes in the data set have zero variance
# Standard deviation for rescaling numeric attributes is now calculated based on the
# whole sample, instead of seoarately for each single team 
#################################
faultlines.default <- function(data,group.par="NA",attr.type=NA,attr.weight=NA,rescale=NA,method="asw",metric="euclid",maxgroups=6, i.level=FALSE, cores=1, usesghomo=FALSE,by.attr=FALSE,...)
{
   # Load required Libraries (if not already loaded)
   library(nFactors)
   library(psych)
   library(flexmix)
   library(nnet) 
   library(parallel)
  
  
  
   result <- list()

   if (group.par != "NA")
   {

      if (!group.par %in% names(data)) 
      {
         stop(paste("Dataset doesn't contain a column named \'",group.par,"\'.",sep=""))
      }
      group.vect <- data[,which(names(data) == group.par)]
      groups <- sort(unique(group.vect))
      
      # group association has been saved in group.vect, so we can delete the group.par column
      data <- data[,-which(names(data) == group.par)] 
   } else
   {
      group.vect <- rep(1,nrow(data))
      groups <- 1     
   }
   
   if (!exists("i.level")) i.level <- FALSE
   if (missing(i.level)) i.level <- FALSE
  
     
   methods <- toupper(c("asw","thatcher","shaw","knippenberg","trezzini","gibson","bezrukova","lcca"))
   if (!toupper(method) %in%  methods)
   {
       stop(paste("Unknown method \'",method,"\'.",sep=""))
   }
   method <- toupper(method)
   
   if (!rescale %in% c(NA,"sd","sd.group"))
      stop(paste("Unknown rescale parameter:" ,rescale))

   
   if (!metric %in% c("euclid","mahal","gower"))
      stop(paste("Unknown metric:",metric))
      
   if ((metric == "mahal") & (!method %in% c("ASW","THATCHER","BEZRUKOVA")))
      stop(paste("Specified Method does not support Mahalanobis Distances."))
      
   if ((metric == "gower") & (!method %in% c("ASW","THATCHER","BEZRUKOVA")))
      stop(paste("Specified Method does not support Manhattan Distances (gower metric)."))
      
   if ((i.level == TRUE) & (method != "ASW"))
      stop(paste("Specified method does not support individual-level faultlines."))
      

   options(warn=-1)         
   if (all(is.na(attr.type)))
   {
      
          attr.type <- rep("numeric",ncol(data))
      
   }
   if (length(attr.type) != ncol(data))
      stop(paste("Length of attr.type does not match the number of columns in the dataset."))
   if (!all(attr.type %in% c("numeric","nominal")))
      stop(paste("Attribute types other than 'nominal' or 'numeric' have been specified."))
   if (sum(is.na(apply(as.matrix(data[,attr.type == "numeric"]),2,as.numeric))) > 0)
   {
      options(warn=0)
      stop(paste("Column(s) declared as 'numeric' could not be treated as numbers."))
   }  
    	
      	
   if (all(is.na(attr.weight)))
   {
      attr.weight <- rep(1,ncol(data))
   }

   if (length(attr.weight) != ncol(data))
      stop(paste("Length of attr.weight does not match the number of columns in the dataset."))
   if (!is.numeric(attr.weight))
      stop(paste("Attribute weights must be numeric.")) 

   if ((method %in% c("TREZZINI","SHAW")) & (any(attr.type != "nominal"))) 
      stop("Selected Method requires all attributes to be nominal-scaled")

   #if ((method %in% c("KNIPPENBERG")) & (any(attr.type != "numeric"))) 
    #  stop("The current implementation of the selected method requires all attributes to be numeric.")
           
      
   if ((any(attr.weight != 1)) & (method %in% c("TREZZINI","SHAW","KNIPPENBERG","GIBSON","LCCA")))
      cat("WARNING: Selected method does not support attribute weighting. Weighting parameter is ignored. \n")

   if ((method == "ASW") & (maxgroups < 2))
   {
       maxgroups <- 2
       cat("Note: The maximum number of subgroups has been set to the minimum value of 2 \n")
   }

   if ((method %in% c("THATCHER","BEZRUKOVA")) & (maxgroups != 2))
      cat("Note: The selected method limits the Number of subgroups to 2 \n")    
      
   # rescale numeric attributes
   if (!is.na(rescale))
   {
   	   if (rescale == "sd")
     
	   {
	          attr.numeric <- sort(which(attr.type == "numeric"),decreasing=FALSE)   
	       
	          for (attr in attr.numeric)
	          {
	             rescale.factor <- sd(data[,attr])
	             
	             if ((!is.na(rescale.factor)) & (rescale.factor != 0))
	             { 
	                data[,attr] <- data[,attr] / rescale.factor
	             }
	          }
	   }
   	   if (rescale == "sd.group")
     
	   {
          attr.numeric <- sort(which(attr.type == "numeric"),decreasing=FALSE)   
       
          for (attr in attr.numeric)
          {
             for (group in groups)
             {
                rescale.factor <- sd(data[group.vect == group,attr])
                
                if (!is.na(rescale.factor) & (rescale.factor != 0))
                { 
                   data[group.vect == group,attr] <- data[group.vect == group,attr] / rescale.factor
                }
             }
          }
	   }
   }

       

   # dummycode nominal attributes for thatcher, bezrukova and asw methods
   dummycols <- c()
   dummylist <- as.list(1:ncol(data))
   if (method %in% c("THATCHER","ASW","BEZRUKOVA"))
   {
      attr.nominal <- sort(which(attr.type == "nominal"),decreasing=TRUE)
   

      for (attr in attr.nominal)
      {
          cats <- sort(unique(data[,attr]))
          
          
       
          dummies <- c()
          for (category in cats)
          {
              dummyvar <- rep(0,nrow(data))
              dummyvar[data[,attr] == category] <- sqrt(2) / 2 
              dummies <- cbind(dummies,dummyvar)
          }
          dummies <- as.data.frame(dummies)
          
          dummylist[[attr]] <- dummylist[[attr]] + 1:ncol(dummies) - 1
          if (length(dummylist) > attr)
          {
            for (entry in (attr+1):length(dummylist))
            {
              dummylist[[entry]] <- dummylist[[entry]] + ncol(dummies) - 1
            }
          }
          
          names(dummies) <- paste(names(data)[attr],cats,sep=".")
          dummycols <- dummycols + ncol(dummies) - 1
          dummycols <- c(attr:(attr + ncol(dummies) - 1),dummycols)
       
          data.new <- c()
          data.names <- names(data)
          if (attr > 1) 
          {
       	   data.new <- cbind(data[,1:(attr-1)],dummies)
       	   names(data.new) <- c(names(data)[1:(attr-1)],names(dummies))
       	   # agl 5.9.2017: adjust attribute weights for nominal attributes by dividing them by the number of categories
       	   # in order to assure equality between 2-cat nominal and 2-cat numeric attributes
       	   #if (by.attr)
       	   #{
       	  #   if (length(cats) == 1) 
       	  #   {
       	  #     attr.weight.new <- c(attr.weight[1:(attr-1)],attr.weight[attr])
       	  #   } else attr.weight.new <- c(attr.weight[1:(attr-1)],rep(attr.weight[attr]/(length(cats)-1),length(cats)))
       	  # } else
       	  # {
       	     attr.weight.new <- c(attr.weight[1:(attr-1)],rep(attr.weight[attr],length(cats)))
       	  # }
          } else 
          {
       	   data.new <- dummies
       	   names(data.new) <- names(dummies)
       	   #if (by.attr)
       	   #{
       	  #   if (length(cats) == 1) 
       	  #   {
       	  #     attr.weight.new <- c(attr.weight[attr])
       	  #   } else attr.weight.new <- rep(attr.weight[attr]/(length(cats)-1),length(cats))
       	  # } else
       	  # {
       	     attr.weight.new <- rep(attr.weight[attr],length(cats))
       	  # }
          }
          # set weight of last dummy variable to zero in order to remove any impact of redundant information
          #if (by.attr & (length(cats) > 1)) attr.weight.new[length(attr.weight.new)] <- 0 
          
          if (attr < ncol(data)) 
          {
       	   data.new.names <- names(data.new)
       	   data.new <- cbind(data.new,data[,(attr+1):ncol(data)])
       	   names(data.new) <- c(data.new.names,names(data)[(attr+1):ncol(data)])
       	   attr.weight.new <- c(attr.weight.new,attr.weight[(attr+1):length(attr.weight)])
          }    
          data <- data.new
          

          attr.weight <- attr.weight.new

          
      }
   }
 


   fl <- NA
   fl.groups <- NA
   fl.adjm <- NA
   fl.num_sg <- NA
   fl.mbr <- NA
   fl.ind <- NA
   fl.nofsg <- NA
   distmat <- NA
   
 
# parallel processing


       attr.weight.team <- attr.weight
       dummycols.team <- dummycols
       result <- mclapply(groups,faultlines.calc,group.vect=group.vect,data=data,method=method,maxgroups=maxgroups, metric=metric, attr.weight.team=attr.weight.team, attr.type = attr.type, dummycols.team=dummycols.team, dummylist=dummylist, quiet=TRUE,i.level=i.level,mc.cores=cores,usesghomo=usesghomo,by.attr=by.attr,mc.preschedule=TRUE)

 
     

   #result <- unlist(result)
   #res <- as.data.frame(do.call(rbind,result))
   result$method=method
   result$i.level=i.level
   result$metric=metric
   class(result) <- "aswclust"
   return(result)
}


###########################################################################################################
###########################################################################################################
###########################################################################################################
# convert vector to text string
v2t <- function(x)
{

	if (all(!is.na(x)))
	{
		x.text <- ""
	    for (i in 1:length(x))
	    {
	    	x.text <- paste(x.text,x[i],sep=" ")
	    }
	} else
	{
		x.text <- "NA"
	}

    return(x.text)
}


print.aswclust <- function(x, ...)
{
    
    pars <- 3
    n.teams <- length(x) - pars

    res <- (matrix(ncol=5,nrow=0))

    
    for (i in 1:n.teams)
    {
    	i.tab <- NA
    	i.freq <- NA
    	if (any(is.null(x[[i]]$mbr_to_subgroups))) x[[i]]$mbr_to_subgroups <- NA
    	if (any(!is.na(x[[i]]$mbr_to_subgroups))) 
    	{
    		
    		i.tab <- as.data.frame(table(x[[i]]$mbr_to_subgroups))
    	    if (ncol(i.tab) == 2) i.freq <- i.tab[,2]
    	}
    	res <- rbind(res,c(x[[i]]$team,
    	                   x[[i]]$fl.value,
    	                   v2t(x[[i]]$mbr_to_subgroups),
    	                   v2t(x[[i]]$number_of_subgroups),
    	                   v2t(i.freq)))
    }
    res <- as.data.frame(res)
    names(res) = c("team","fl.value","mbr_to_subgroups","number_of_subgroups","subgroup_sizes")
    if (x$i.level == TRUE)
    {
    	res$subgroup.association <- "multiple"
    	res$subgroup.sizes <- "multiple"
    }
    
    print(res)
}

summary.aswclust <- function(object, ...)
{
	
	pars <- 3
	n.teams <- length(object) - pars
	method <- object$method
	metric <- object$metric
	level <- object$i.level
	
	ifelse(object$i.level == TRUE, level <- "individual", level <- "team")
    
  teams <- c()
  flt <- c()
  fli <- list()	
  nof_sg <- list()
  nof_sg_text <- c()
  mbr_sg <- list()
  adjm <- list()
  distmat <- list()
	for (i in 1:n.teams)
    {
    	teams <- c(teams,object[[i]]$team)
    	if (is.null(object[[i]]$fl.value)) object[[i]]$fl.value <- NA
    	flt <- c(flt,object[[i]]$fl.value)
    	nof_sg[[i]] <- object[[i]]$number_of_subgroups
        nof_sg_text <- c(nof_sg_text,v2t(nof_sg[[i]]))
    	fli[[i]] <- object[[i]]$fl.ind
    	distmat[[i]] <- data.frame(object[[i]]$distmat)
    	
    	if (object$i.level == TRUE)
    	{
    		mbr_sg[[i]] <- data.frame(object[[i]]$mbr_to_subgroups)
    		adjm[[i]] <- data.frame(object[[i]]$adjm)
    		
    	} else
    	{
    		mbr_sg[[i]] <- object[[i]]$mbr_to_subgroups  		
    		adjm[[i]] <- NA 
    		
    	}
    }


    fltab <- data.frame(team = teams, fl.value = flt, number_of_subgroups = nof_sg_text)
	res <- list(n.teams = n.teams, 
	            method = method, 
	            metric = metric,
	            level = level,
	            fltab = fltab,
	            fl.team = flt,
	            fl.ind = fli,
	            mbr_to_subgroups = mbr_sg,
	            number_of_subgroups = nof_sg,
	            adjm = adjm,
	            distmat = distmat,
	            long = convert.to.long(object))
	class(res) = "summary.aswclust"
	res
}

print.summary.aswclust <- function(x, ...)
{
    cat("Number of Teams: ")
    cat(x$n.teams,"\n\n")
    
    cat("Calculation features: \n")
    cat("\tMethod: \t",x$method,"\n")
    cat("\tLevel: \t",x$level,"\n")
    cat("\tMetric: \t", x$metric, "\n")
    cat("\n")
    
    for (team in 1:x$n.teams)
    {
    	tstr <- paste("Team ",team," (",x$fltab[team,1],"):",sep="")
    	tul <- rep("=",nchar(tstr))
    	cat("\n",tstr,"\n",sep="")
    	cat(as.character(tul),"\n",sep="")
    	
	    cat("Faultline Strength:\n")
	    print(x$fl.team[[team]])
	    cat("\n")
	    
	    cat("Individual Faultline Strengths (silhouette widths):\n")
	    print(x$fl.ind[[team]])
	    cat("\n")
	    
	    cat("Member to Subgroup Association:\n")
	    if (length(x$mbr_to_subgroups) > 0) print(x$mbr_to_subgroups[[team]]) else print(NA)
	    cat("\n")
	    
	    cat("Number of Subgroups:\n")
	    if (length(x$number_of_subgroups) > 0) print(x$number_of_subgroups[[team]]) else print(NA)
	    cat("\n")
	   
	    if (!is.na(x$adjm))
	    {
	    	cat("Subgroup Network:\n")
	    	print(x$adjm[[team]])
	    	cat("\n")
	    }
	    
	    cat("Distances:\n")
	    print(x$distmat[[team]])
	    cat("\n")
	    
    }
     
}




#################################
# function: get_weights
#################################
# Version:    1.0
# Date:       12.09.2017
# Author:     A. Glenz
#             a.glenz@psychologie.uzh.ch
#
# Purpose:
# Extract attribute weights for faultline calculation based on a given data set and a given group outcome vector
#################################
get_weights <- function(data,group.par="NA",attr.type=NA,n.rand=0,rescale=NA,method="asw",metric="euclid",maxgroups=6, i.level=FALSE, cores=1, usesghomo=FALSE,criterion,repetitions=1,quiet=FALSE)
{
  library(parallel) 
  # analyze dataset
  if (group.par=="NA" & !quiet)
    stop("Dataset must contain data from multiple groups (we recommend at least 20).")
  
  n.attr <- ncol(data)-1
  if ((n.attr + n.rand) < 6 & !quiet)
  {  
    warning(paste("We recommend at least 6 attributes (including random noise attributes). You specified only",n.attr,"data attributes and",n.rand,"random noise attributes."),immediate.=TRUE)
    ans <- readline("Do you want to continue (y/n)? ")
    if (ans != "y") 
      stop("Aborted by User")
  }  
  
  n.groups <- length(unique(data[[group.par]])) 
  if (n.groups < 20 & !quiet)
  {  
    warning(paste("We recommend that the data contains at least 20 groups. Your data has only",n.groups,"groups."),immediate.=TRUE)
    ans <- readline("Do you want to continue (y/n)? ")
    if (ans != "y") 
      stop("Aborted by User")
  }
  
  # generating random noise attributes
  if (n.rand > 0)
  {
    for (rcol in 1:n.rand)
    {
      rndvar <- rnorm(nrow(data))
      data <- cbind(data,rndvar)
    }
    names(data)[(length(names(data))-n.rand+1):length(names(data))] <- paste("rand",1:n.rand,sep="")
    attr.type <- c(attr.type,rep("numeric",n.rand))
  }
  message(paste("Attribute names are:"))
  attr.names <- names(data)[-which(names(data)==group.par)]
  print(paste(attr.names,sep=", "))
  if (!quiet) 
  {
    ans <- readline("Is this correct (y/n)? ")
    if (ans != "y") 
      stop("Aborted by User")
  }         
  
  if (length(criterion) != n.groups  & !quiet)
    stop(paste("The length of the criterion vector (",length(criterion),") does not match the number of groups (",n.groups,")",sep=""))
  
  
  
  
  ######### now do it...
  
  require(R.utils)
  
  n.spec <- n.attr + n.rand 
  
  combos <- intToBin(1:((2^n.spec)-1))
  combos <- strsplit(combos,"")
  
  starttime <- Sys.time()
  
  numcalc <- length(combos) * repetitions
  calc <- 0
  res <- list()
  for (turn in 1:repetitions)
  {
    
    fs.combos <- NA
    fs.combos <- list() 
    
    
    for (i in 1:length(combos))
    { 
      calc <- calc + 1
      if (calc > 1)
      { 
        curtime <- Sys.time()
        print(paste("Calculating combination ",calc," of ",numcalc,". EOT at ",(curtime - starttime) / (calc-1) * (numcalc-calc-1) + curtime,sep=""))
      } else
        print(paste("Calculating combination ",calc," of ",numcalc,sep=""))
      
      
      # re-generating random noise attributes
      if (n.rand > 0)
      {
        
        data <- data[,1:(n.attr+1)]
        for (rcol in 1:n.rand)
        {
          rndvar <- rnorm(nrow(data))
          data <- cbind(data,rndvar)
        }
        names(data)[(length(names(data))-n.rand+1):length(names(data))] <- paste("rand",1:n.rand,sep="")
      }
      
      x <- faultlines(data, group.par=group.par,attr.type=attr.type,cores=cores,by.attr=TRUE,method=method,attr.weight=as.numeric(combos[[i]])[c(1:(n.attr+n.rand))],quiet=quiet,metric=metric,maxgroups=maxgroups, i.level=i.level, cores=1, usesghomo=usesghomo)
      fs.combos[[i]] <- summary(x)$fltab$fl.value
      
      
    }
    
    #########  
    
    
    # calculate outcome correlations
    fs.cor = fs.cor.p <- c()
    for (i in 1:length(combos))
    {
      fs.cor[i] <- cor(fs.combos[[i]],criterion)
      fs.cor.p[i] <- cor.test(fs.combos[[i]],criterion)$p.value
    }
    
    
    
    combos.df <- as.data.frame(matrix(as.numeric(unlist(combos)),ncol=n.spec,byrow=TRUE))
    
    combos.sorted <- combos.df[order(fs.cor),]
    names(combos.sorted) <- attr.names
    
    
    numsig <- sum(fs.cor.p <= p.level)
    
    
    ############### calculate weights
    chi.pos = chi.neg <- c()
    chi.pos.p = chi.neg.p <- c()
    w.pos = w.neg <- c()
    occ.pos = occ.neg <- c()
    exp.pos = exp.neg <- c()
    
    for (a in 1:(n.attr + n.rand))
    {
      # h1: positive effects
      roi <-
        fs.cor[order(fs.cor)] > 0 &
        1:nrow(combos.sorted) > as.integer(nrow(combos.sorted) * .5)
      roi.size <- sum(roi)
      roi.expected <- floor(roi.size / 2)
      exp.pos[a] <- roi.expected
      occ.pos[a] <- sum(combos.sorted[roi, a])
      
      chi.pos[a] <- 0
      if (occ.pos[a] > 0)
      {
        ct <-
          chisq.test(c(roi.expected, sum(combos.sorted[roi, a])), p = c(.5, .5))
        chi.pos[a] <-
          round(ct$statistic * (2 * (sum(combos.sorted[roi, a]) >= roi.expected) - 1), 2)
      }
      
      
      
      chi.pos.p[a] <- 1
      if (sum(combos.sorted[roi, a]) >= roi.expected)
        chi.pos.p[a] <- ct$p.value
      
      # h2: negative effects
      roi <-
        fs.cor[order(fs.cor)] < 0 &
        1:nrow(combos.sorted) <= as.integer(nrow(combos.sorted) * .5)
      roi.size <- sum(roi)
      roi.expected <- floor(roi.size / 2)
      exp.neg[a] <- roi.expected
      occ.neg[a] <- sum(combos.sorted[roi, a])
      
      chi.neg[a] <- 0
      if (occ.neg[a] > 0)
      {
        ct <-
          chisq.test(c(roi.expected, sum(combos.sorted[roi, a])), p = c(.5, .5))
        chi.neg[a] <-
          round(ct$statistic * (2 * (sum(combos.sorted[roi, a]) >= roi.expected) - 1), 2)
      }
      
      
      
      chi.neg.p[a] <- 1
      if (sum(combos.sorted[roi, a]) >= roi.expected)
        chi.neg.p[a] <- ct$p.value
    }
    
    
    # calculate relative attribte weights
    w.pos <- chi.pos * 0
    if (any(chi.pos > 0))
      w.pos <- chi.pos / max(chi.pos)
    
    w.neg <- chi.neg * 0
    if (any(chi.neg > 0))
      w.neg <- chi.neg / max(chi.neg)
    
    
    chi.res <-
      data.frame(
        Attribute = names(combos.sorted)[1:(n.attr + n.rand)],
        exp.pos = format(exp.pos, nsmall = 0),
        occ.pos = format(occ.pos, nsmall = 0),
        chi2.pos = format(chi.pos, nsmall = 2),
        w.pos = w.pos,
        exp.neg = format(exp.neg, nsmall = 0),
        occ.neg = format(occ.neg, nsmall = 0),
        chi2.neg = format(chi.neg, nsmall = 2),
        w.neg = w.neg
      )
    
    
    
    colnames(chi.res) <-
      c(
        "Attribute",
        "exp_pos",
        "obs_pos",
        "chi2_pos",
        "weight_pos",
        "exp_neg",
        "obs_neg",
        "chi2_neg",
        "weight_neg"
      )
    
    
    
    
    res[[turn]] <- chi.res
  }
  
  if (repetitions > 1)
  {
    for (i in 1:length(res))
    {
      for (col in 2:ncol(res[[i]]))
      {
        res[[i]][,col] <- as.numeric(as.character(res[[i]][,col]))
      }
    }
    sum.res <- res[[1]] * 0
    for (i in 1:length(res))
    {
      sum.res <- sum.res + res[[i]]
    }
    
    
    
    # perform chi-statistics on cumulated observations
    chi.pos = chi.neg <- c()
    chi.pos.p = chi.neg.p <- c()
    w.pos = w.neg <- c()
    occ.pos = occ.neg <- c()
    exp.pos = exp.neg <- c()
    
    for (a in 1:(n.attr + n.rand))
    {
      # h1: positive effects
      exp.pos[a] <- sum.res[a,2]
      occ.pos[a] <- sum.res[a,3]
      chi.pos[a] <- 0
      
      if (occ.neg[a] > 0)
      {
        ct <- chisq.test(c(exp.pos[a], occ.pos[a]), p = c(.5, .5))
        chi.pos[a] <-
          round(ct$statistic * (occ.pos[a] >= exp.pos[a]), 2)
      }
      
      chi.pos.p[a] <- 1
      if (occ.pos[a] >= exp.pos[a])
        chi.pos.p[a] <- ct$p.value
      
      # h2: negative effects
      exp.neg[a] <- sum.res[a,6]
      occ.neg[a] <- sum.res[a,7]
      chi.neg[a] <- 0
      
      if (occ.neg[a] > 0)
      {
        ct <- chisq.test(c(exp.neg[a], occ.neg[a]), p = c(.5, .5))
        chi.neg[a] <-
          round(ct$statistic * (occ.neg[a] >= exp.neg[a]), 2)
      }
      
      chi.neg.p[a] <- 1
      if (occ.neg[a] >= exp.neg[a])
        chi.neg.p[a] <- ct$p.value
    }
    
    
    
    # calculate relative attribte weights
    w.pos <- chi.pos * 0
    if (any(chi.pos > 0))
      w.pos <- chi.pos / max(chi.pos)
    
    w.neg <- chi.neg * 0
    if (any(chi.neg > 0))
      w.neg <- chi.neg / max(chi.neg)
    
    
    chi.res <-
      data.frame(
        Attribute = names(combos.sorted)[1:(n.attr + n.rand)],
        exp.pos = format(exp.pos, nsmall = 0),
        occ.pos = format(occ.pos, nsmall = 0),
        chi2.pos = format(chi.pos, nsmall = 2), 
        w.pos = w.pos,
        exp.neg = format(exp.neg, nsmall = 0),
        occ.neg = format(occ.neg, nsmall = 0),
        chi2.neg = format(chi.neg, nsmall = 2), 
        w.neg = w.neg
      )
    
    
    
    colnames(chi.res) <-
      c(
        "Attribute",
        "exp_pos",
        "obs_pos",
        "chi2_pos",
        "weight_pos",
        "exp_neg",
        "obs_neg",
        "chi2_neg",
        "weight_neg"
      )
    
  }
  chi.res$weight_pos[chi.res$weight_pos < 0] <- 0
  chi.res$weight_neg[chi.res$weight_neg < 0] <- 0
  
  # normalize weights
  if (n.rand > 0)
  {   
    w.pos <- chi.res$weight_pos
    w.pos[(length(w.pos) - n.rand + 1):length(w.pos)] <- 0
    if (any(w.pos > 0)) w.pos <- w.pos / max(w.pos) 
    chi.res$weight_pos <- w.pos
    
    w.neg <- chi.res$weight_neg
    w.neg[(length(w.neg) - n.rand + 1):length(w.neg)] <- 0
    if (any(w.neg > 0)) w.neg <- w.neg / max(w.neg) 
    chi.res$weight_neg <- w.neg
    
  }
  return(chi.res)
}


