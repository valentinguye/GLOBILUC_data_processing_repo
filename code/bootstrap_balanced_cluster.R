library(boot)
library(sandwich)

## Make some necessary objects
# unbalanced panel data
data("PetersenCL", package = "sandwich")
data <- PetersenCL

# list of parameters related to the dataset, and the clustering variable
cluster_var <- "firm"

# names and numbers of clusters of size s
par_list <- list(cluster_variable = cluster_var, 
                 cluster_names = unique(data[,cluster_var]),
                 number_clusters = length(unique(data[,cluster_var])))

ran.gen_cluster_blc <- function(original_data, arg_list){
  # to store 
  cl_boot_dat <- list()
  
  # such that we don't have to call it from the list every time
  cluster_var <- arg_list[["cluster_variable"]]
  
  # non-unique names of clusters (repeated when there is more than one obs. in a cluster) 
  nu_cl_names <- as.character(original_data[,cluster_var]) 
  
  # sample, in the vector of names of clusters as many draws as there are clusters, with replacement
  sample_cl <- sample(x = arg_list[["cluster_names"]], 
                      size = arg_list[["number_clusters"]], 
                      replace = TRUE) 
  
  # because of replacement, some names are sampled more than once
  # we need to give them a new cluster identifier, otherwise a cluster sampled more than once 
  # will be "incorrectly treated as one large cluster rather than two distinct clusters" (by the fixed effects) (Cameron and Miller, 2015)    
  sample_cl_tab <- table(sample_cl)
  
  for(n in 1:max(sample_cl_tab)){ # from 1 to the max number of times a cluster was sampled bc of replacement
    # vector to select obs. that are within the clusters sampled n times. 
    # seems slightly faster to construct the names_n object beforehand 
    names_n <- names(sample_cl_tab[sample_cl_tab == n])
    sel <- nu_cl_names %in% names_n

    # select data accordingly to the cluster sampling (duplicating n times observations from clusters sampled n times)
    clda <- original_data[sel,][rep(seq_len(sum(sel)), n), ]
    
    #identify row names without periods, and add ".0" 
    row.names(clda)[grep("\\.", row.names(clda), invert = TRUE)] <- paste0(grep("\\.", row.names(clda), invert = TRUE, value = TRUE),".0")
    
    # add the suffix due to the repetition after the existing cluster identifier. 
    clda[,cluster_var] <- paste0(clda[,cluster_var], sub(".*\\.","_",row.names(clda)))
    
    # stack the bootstrap samples iteratively 
    cl_boot_dat[[n]] <- clda
  }
  return(bind_rows(cl_boot_dat))
}

#test that the returned data are the same dimension as input
test_boot_d <- ran.gen_cluster_blc(original_data = data,
                                   arg_list = par_list)

dim(test_boot_d)
dim(data)
# test new clusters are not duplicated (correct if anyDuplicated returns 0)
base::anyDuplicated(test_boot_d[,c("firm","year")])


# test that it computes the same standard error as sandwich::vcovBS, for statistic being a regression coefficient 
est_fun <- function(est_data){
  
  est <- lm(as.formula("y ~ x"), est_data)
  
  # statistics we want to evaluate the variance of:
  return(est$coefficients)
}
# see ?boot::boot for more details on these arguments
set.seed(1234)
boot(data = data, 
     statistic = est_fun, 
     ran.gen = ran.gen_cluster_blc,
     mle = par_list,
     sim = "parametric",
     parallel = "no",
     R = 400)

set.seed(1234)
sdw_bs <- vcovBS(lm(as.formula("y ~ x"), PetersenCL), cluster = ~firm, R=400)#
sqrt(sdw_bs["x","x"])

sdw_cl <- vcovCL(lm(as.formula("y ~ x"), data), cluster = ~firm)
sqrt(sdw_cl["x","x"])

ran.gen_blc <- function(original_data, arg_list){
  rowids <- row.names(original_data)
  
  new_rowids <- sample(x = rowids, 
                       size = nrow(original_data), 
                       replace = TRUE)
  
  return(original_data[new_rowids,])
}

set.seed(1234)
boot(data = data, 
     statistic = est_fun, 
     ran.gen = ran.gen_blc,
     mle = list(),
     sim = "parametric",
     parallel = "no",
     R = 400)