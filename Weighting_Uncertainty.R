# This script calculates the uncertainty of the weighted product derived in Weighting.R (or both the uncertainty and the weighted product)
# The uncertainty standard deviation reflects the true discrepancy of the hybrid product with respect to observations


# The method is based on the ensemble dependence transformation process presented in Bishop and Abramowitz (2013) and also explained in Hobeichi et al. (2018)
# I advice the users to refer to the methods section in 
# Hobeichi, S., Abramowitz, G., Evans, J., and Ukkola, A.: Derived Optimal Linear Combination Evapotranspiration (DOLCE): 
# a global gridded synthesis ET estimate, Hydrol. Earth Syst. Sci., 22, 1317-1336, https://doi.org/10.5194/hess-22-1317-2018, 2018.



# le_df is a dataframe containing observed latent heat flux from various flux tower sites as well as modelled latent heat from several models 
# every row corresponds to a unique site-time step
# var_names a vector of colnames (in le_df) corresponding to the participating models


### remove all the records (rows) where any of the participating models or observations are not available  
le_df <- le_df[rowSums(is.na(le_df[ , c(var_names,"LE_obs")])) == 0, ] 

### dataframe containing the modeled latent heat
model <- le_df[, var_names]

### vector containing observed latent heat
obs <- le_df[,"LE_obs"] 


### calculate the bias correction terms (bc1, bc2,... bcn)
le_df_error <- model - obs ### subtract the observation from the values of the models
bc_term <- colMeans(le_df_error, na.rm=TRUE) ## calculate the mean difference for each model 

### bias correct the models 

# initialization 
model_bc <- model
for(i in 1:ncol(model_bc)){
  
  # bias correct each model
  model_bc[,i]<- model_bc[,i]-bc_term[i] 
}


### Calculate the weights

error_df <- model_bc - obs

# error covariance matrix
M_matrix<-data.matrix(error_df)
M_cov<- cov(M_matrix) 

# number of participating models
model_count <- ncol(M_cov) 

# Unit column
unit_col <- matrix(data=1, nrow=model_count) 

# inverse of the covariance matrix
M_cov_inv <- solve(M_cov) 

# Unit column transposed
unit_transpose <- t(unit_col) 

### calculate optimal weights
weights <- (M_cov_inv %*% unit_col)/(unit_transpose %*% M_cov_inv %*% unit_col)[1,1]

### dataframe containing the weights and the bias correction terms
w_df <- data.frame(w=weights, bc=bc_term)
w_df$raster_names <- raster_names ## raster names of the models being merged

### Calculate the weighted vector by applying the optimal weights on the bias corrected models
weighted_vec <- rep(0,nrow(model_bc))
for(i in 1:ncol(model_bc)){
  weighted_vec<- weighted_vec + weights[i]*model_bc[,i]
}

### use the same variable names as in the method documentation to make the script more understandable
x_df <- model_bc

mu <- weighted_vec

w <- weights

######## Calculate alpha - parameter ####################################
k <- length(var_names)  ## number of participating models

# check if all the weights are positive
all_positive <- TRUE
for(i in 1: k)
  if(w[i] <0)
    all_positive <- FALSE

if(all_positive){
  alpha <- 1
}else alpha <- 1 - k*min(w)  ## In case there is a negative weight, the smallest weight is turned to zero

w_tilda <- sapply(w, function(x) (x-min(w))/alpha)
########## x_bar ######################################################## 
x_bar <- rowMeans(x_df)

########## z_df #########################################################
z_df <- x_df
for(i in 1:k)
  z_df[,i] <- x_bar + alpha*(x_df[,i]-x_bar)

########## Calculate (Se)^2 #############################################
se2 <- sum((mu - obs)^2)/(length(obs)-1)
########## Calculate Beta ###############################################

# 1. calculate the denonminator Beta
den_term <- z_df
for(i in 1:k)
  den_term[,i] <- z_df[,i] - mu

den_term <- sapply(den_term, function(x) x^2)
den_term_sum_col <- Get_weighted_vec(den_term, w_tilda)
den_term_sum_col_row <- sum(den_term_sum_col)/length(obs)

# 2. calculate Beta
beta <- sqrt(se2/den_term_sum_col_row)

##########################################################################


### If the weighted product has been already created in Weighting.R read it here:
#mu_global <- raster(path_to_hybrid_product)


### It the weighted product has not been already created, calcuate it in this step
# apply the weights and bias terms summarized in w_df to derive the hybrid product: hybrid_raster =  w1(model_raster1-bc1) + ... + wn(model_rastern-bcn)
for(i in 1:nrow(w_df)){
  
  temp <- row_df$w[i] *(raster(w_df$raster_names[i]) - w_df$bc[i])
  
  if(i ==1){
    mu_global <- temp
  }else{
    mu_global <- mu_global + temp
  }
  i <- i+1
}

writeRaster(mu_global, "hybrid_raster.tif", format="GTiff", overwrite=TRUE)

### Read the models and bias correct them in x1, x2, x3, .. xn 
for(i in 1:nrow(w_df)){
  
  assign(paste0("x",i), raster(w_df$raster_names[i]))
}


############# Calculate the mean of the bias corrected models #############
first <- TRUE
for(i in 1:nrow(w_df)){
  
  
  # bias correct the models
  temp_m_bc <- eval(parse(text=paste0("x",i," - w_df$bc[i]")))
  
  assign(paste0("x",i), temp_m_bc)
  
  # Calculate the sum of the bias corrected models
  if(first){
    sum_temp_m_bc <- temp_m_bc
    first <- FALSE
  }else{
    sum_temp_m_bc <- sum_temp_m_bc + temp_m_bc
  }
}

# mean of bias corrected models
x_bar_global <- sum_temp_m_bc/k

########## Now we want to calculate the transformed models, we will use the same alpha and Beta

# initialization
counter <- 1
first <- TRUE

for(i in 1:nrow(w_df)){
  
  ### calculate the transformed models x_tilda1, x_tilda2, ... x_tildan
  temp_x_tilda <- eval(parse(text=paste0("mu_global + beta*(x_bar_global+ alpha* (x",i,"- x_bar_global)-mu_global)")))
  assign(paste0("x_tilda",i), temp_x_tilda)
  
  
  ### Calculate the variance of the transformed models around the weighted product
  temp <- eval(parse(text=paste0("x_tilda",i,"- mu_global")))
  assign(paste0("sigma",i), temp)
  
  temp <- eval(parse(text=paste0("sigma",i,"^2")))
  assign(paste0("sigma",i), temp)
  
  temp <- eval(parse(text=paste0("sigma",i,"* w_tilda[counter]")))
  
  
  counter <- counter+1
  if(first){
    sigma <- temp
    first <- FALSE
    
  }else{
    sigma <- sigma + temp
    
  }
}

#####  Now we can calculate the uncertainty standard deviation of the weighted product
sigma <- sqrt(sigma)

writeRaster(sigma, "sigma_raster.tif", format="GTiff", overwrite=TRUE)
