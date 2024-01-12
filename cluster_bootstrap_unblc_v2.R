data("PetersenCL", package = "sandwich")

# construct an unbalanced panel with different numbers of firms every year. 
datalist <- list()
for(yr in 1:7){
  datalist[[yr]] <- subset(PetersenCL, (firm %in% c(1:(2*yr)) & year == yr))
}
data <- bind_rows(datalist)

cluster_var <- "year"

cluster_names = unique(data[,cluster_var])

# get the different cluster sizeS. This is necessary to cluster bootstrapping with clusters of different sizes. 
sizes <- table(data[,cluster_var])

# don't bother the size of every cluster. Use the average cluster size only. 
avg_cl_size <- mean(sizes)
# sample size of original data
N <- nrow(data)
# Number of draws, for the bootstrap replicate data sets to be as large as original data, on average. 
n_draws <- N/avg_cl_size

#par_list <- list(n_draws = n_draws)

ran.gen_cluster <- function(original_data, arg_list){
  
  cl_boot_dat <- NULL
  
  # for later: non-unique names of clusters (repeated when there is more than one obs. in a cluster) 
  nu_cl_names <- as.character(original_data[,cluster_var]) 
  
  # sample, in the vector of names of clusters, N/avg_cl_size
  sample_cl <- sample(x = cluster_names, 
                      size = n_draws, 
                      replace = TRUE) 
  
  # because of replacement, some names are sampled more than once
  # we need to give them a new cluster identifier, otherwise a cluster sampled more than once 
  # will be "incorrectly treated as one large cluster rather than two distinct clusters" (by the fixed effects) (Cameron and Miller, 2015)    
  sample_cl_tab <- table(sample_cl)
  
  for(n in 1:max(sample_cl_tab)){ # from 1 to the max number of times a name was sampled bc of replacement
    # vector to select obs. that are within the sampled clusters. 
    names_n <- names(sample_cl_tab[sample_cl_tab == n])
    sel <- nu_cl_names %in% names_n
    
    # select data accordingly to the cluster sampling (duplicating n times observations from clusters sampled n times)
    clda <- original_data[sel,][rep(seq_len(sum(sel)), n), ]
    
    #identify row names without periods, and add ".0" 
    row.names(clda)[grep("\\.", row.names(clda), invert = TRUE)] <- paste0(grep("\\.", row.names(clda), invert = TRUE, value = TRUE),".0")
    
    # add the suffix due to the repetition after the existing cluster identifier. 
    clda[,cluster_var] <- paste0(clda[,cluster_var], sub(".*\\.","_",row.names(clda)))
    
    # stack the bootstrap samples iteratively 
    cl_boot_dat <- rbind(cl_boot_dat, clda) 
  }
  return(cl_boot_dat)
}

#the returned data ARE NOT the same dimension as input data. It only has the same dimension on average over bootstrap replicates. 
test_boot_d <- ran.gen_cluster(original_data = data)
dim(test_boot_d)
dim(data)

# test new clusters are not duplicated (correct if anyDuplicated returns 0)
base::anyDuplicated(test_boot_d[,c("firm","year")])

# custom estimation function
est_fun <- function(est_data){
  
  est <- lm(as.formula("y ~ x"), est_data)
  
  # statistics we want to evaluate the variance of:
  return(est$coefficients)
}

# Run the bootstrap (see ?boot::boot for more details on the arguments)
set.seed(1234)
boot(data = data, 
     statistic = est_fun, 
     ran.gen = ran.gen_cluster,
     mle = list(), # this argument cannot be left empty 
     sim = "parametric",
     parallel = "no",
     R = 400)

# Test that it computes the same standard error as sandwich::vcovBS and compare with the asymptotic solution of sandwich::vcovCL
set.seed(1234)
sdw_bs <- vcovBS(lm(as.formula("y ~ x"), data), cluster = ~year, R=400, type = "xy")#
sqrt(sdw_bs["x","x"])

# Compare with asymptotic solution
sdw_cl <- vcovCL(lm(as.formula("y ~ x"), data), cluster = ~year)
sqrt(sdw_cl["x","x"])


