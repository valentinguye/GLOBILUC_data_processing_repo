# methodology comes from:
# https://stats.stackexchange.com/questions/202916/cluster-boostrap-with-unequally-sized-clusters/202924#202924

library(boot)
library(sandwich)

## Make some necessary objects
# unbalanced panel data
data("PetersenCL", package = "sandwich")
data <- subset(PetersenCL, !(firm %in% c(1:10) & year == 10))

cluster_var <- "firm"

# get the different cluster sizeS. This is necessary to cluster bootstrapping with clusters of different sizes. 
sizes <- table(data[,cluster_var])
u_sizes <- sort(unique(sizes))

# names and numbers of clusters of every sizes
cl_names <- list()
n_clusters <- list()
for(s in u_sizes){
  cl_names[[s]] <- names(sizes[sizes == s])
  n_clusters[[s]] <- length(cl_names[[s]])
}

## Design the bootstrap sampling function. 
# It will tell boot::boot how to sample data at each replicate 
ran.gen_cluster <- function(original_data, arg_list){
  
  cl_boot_dat <- NULL

  # non-unique names of clusters (repeated when there is more than one obs. in a cluster) 
  nu_cl_names <- as.character(original_data[,cluster_var]) 
  
  for(s in arg_list[["unique_sizes"]]){
    # sample, in the vector of names of clusters of size s, as many draws as there are clusters of that size, with replacement
    sample_cl_s <- sample(arg_list[["cluster_names"]][[s]], 
                          arg_list[["number_clusters"]][[s]], 
                          replace = TRUE) 
    
    # because of replacement, some names are sampled more than once
    sample_cl_s_tab <- table(sample_cl_s)
    
    # we need to give them a new cluster identifier, otherwise a cluster sampled more than once 
    # will be "incorrectly treated as one large cluster rather than two distinct clusters" (by the fixed effects) (Cameron 2015)    
    
    for(n in 1:max(sample_cl_s_tab)){ # from 1 to the max number of times a name was sampled bc of replacement
      # vector to select obs. that are within the sampled clusters. 
      sel <- nu_cl_names %in% names(sample_cl_s_tab[sample_cl_s_tab == n])

      # replicate the selected data n times
      clda <- original_data[sel,][rep(seq_len(nrow(original_data[sel,])), n), ]
      
      #identify row names without periods, and add ".0" 
      row.names(clda)[grep("\\.", row.names(clda), invert = TRUE)] <- paste0(grep("\\.", row.names(clda), invert = TRUE, value = TRUE),".0")
      
      # add the suffix due to the repetition after the existing cluster identifier. 
      clda[,cluster_var] <- paste0(clda[,cluster_var], sub(".*\\.","_",row.names(clda)))
      
      # stack the bootstrap samples iteratively 
      cl_boot_dat <- rbind(cl_boot_dat, clda) 
    }
  }  
  return(cl_boot_dat)
}

#test that the returned data are the same dimension as input
test_boot_d <- ran.gen_cluster(original_data = data,
                               arg_list = list(unique_sizes = u_sizes,
                                               cluster_names = cl_names,
                                               number_clusters = n_clusters))
dim(test_boot_d)
dim(data)
# test new clusters are not duplicated (correct if anyDuplicated returns 0)
base::anyDuplicated(test_boot_d[,c("firm","year")])


## Custom ("black-box") estimation function
est_fun <- function(est_data){
  
  est <- lm(as.formula("y ~ x"), est_data)
  
  # statistics we want to evaluate the variance of:
  return(est$coefficients)
}

# see ?boot::boot for more details on these arguments
boot(data = data, 
      statistic = est_fun, 
      ran.gen = ran.gen_cluster,
      mle = list(unique_sizes = u_sizes, 
                 cluster_names = cl_names, 
                 number_clusters = n_clusters),
      sim = "parametric",
      parallel = "no",
      R = 200)

sdw_cl <- vcovCL(lm(as.formula("y ~ x"), data), cluster = ~firm)
sqrt(sdw_cl["x","x"])


## compare with sandwich::vcovBS, on balanced panel data
set.seed(1234)
custom_bs <- boot(data = PetersenCL, 
                  statistic = est_fun, 
                  ran.gen = ran.gen_cluster,
                  mle = list(unique_sizes = u_sizes, 
                             cluster_names = cl_names, 
                             number_clusters = n_clusters),
                  sim = "parametric",
                  parallel = "no",
                  R = 200)
sd(custom_bs$t[,2])

# by hand computation of SE changes nothing 
stat_bar <- mean(custom_bs$t)
df <- as.data.frame(x = list(stat_repl = custom_bs$t, stat_bar = stat_bar))
df <- dplyr::mutate(df, dev2 = (stat_repl - stat_bar)^2)
B <- nrow(df)
custom_bs_se <- sqrt((1/(B-1))*sum(df$dev2))
 
custom_bs_se == sd(custom_bs$t)

sdw_cl <- vcovCL(lm(as.formula("y ~ x"), PetersenCL), cluster = ~firm)
sqrt(sdw_cl["x","x"])


set.seed(1234)
sdw_bs <- vcovBS(lm(as.formula("y ~ x"), PetersenCL), cluster = ~firm, R=200)#
sqrt(sdw_bs["x","x"])


# compare standard error estimate from the two methods
all.equal(sqrt(sdw_bs["x","x"]), sd(custom_bs$t))






identical(as.double(sqrt(sdw_bs["x","x"])), as.double(sd(custom_bs$t[,2])), 
          ignore.bytecode = T, 
          ignore.environment = T,
          ignore.srcref = T)


sdw_cl <- vcovCL(est, cluster = ~firm)
sqrt(sdw_cl["x","x"])


