
#### Optimally combine an ensemble of global products so that the RMSD of the resultant hybrid product and in-situ observation is minimized
#### Function Get_weights_and_bias returns a dataframe containing optimal weights (w1, w2, ..wn) and bias correction terms (bc1, bc2, .. bcn) of participating models
#### Users can then derive a hybrid product as: hybrid = w1(model1-bc1) + w2(model2-bc2) + ..... wn(modeln-bcn)

### The Weights ensure that the RMSD of derived hybrid product with respect to the observation is minimized
### The weights account for dependence between participating product. Error dependence is defined by error covariance


## le_df is a dataframe containing observed latent heat flux from various flux tower sites as well as modelled latent heat from several models
## every row corresponds to a unique site-time step
## var_names a vector of colnames (in le_df) corresponding to the participating models
Get_Weights_and_bias <- function(var_names, le_df){

  # remove all the records (rows) where any of the participating models or observations are not available
  le_df <- le_df[rowSums(is.na(le_df[ , c(var_names,"LE_obs")])) == 0, ]

  # dataframe containing the modeled latent heat
  model <- le_df[, var_names]

  # vector containing observed latent heat
  obs <- le_df[,"LE_obs"]


  ### calculate the bias correction terms (bc1, bc2,... bcn)
  le_df_error <- model - obs ### subtract the observation from the values of the models
  bc_term <- colMeans(le_df_error, na.rm=TRUE) ## calculate the mean difference for each model

  ### bias correct the models

  # initialization
  model_bc <- model
  for(i in 1:ncol(model_bc)){

     # subtract the mean difference
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

  # calculate optimal weights
  weights <- (M_cov_inv %*% unit_col)/(unit_transpose %*% M_cov_inv %*% unit_col)[1,1]

  ### dataframe containing the weights and the bias correction terms
  w_df <- data.frame(w=weights, bc=bc_term)


  return(w_df)

}


## calculate the weights and bias correction terms for participating models. Observed data (LE_obs) and modelled estimates(var_names) are summarized in le_df
w_df <- Get_Weights_and_bias(var_names, le_df)

## Add the Name for the models' rasters being merged
w_df$raster_names <- raster_names


## apply the weights and bias terms summarized in w_df to derive the hybrid product: hybrid_raster =  w1(model_raster1-bc1) + ... + wn(model_rastern-bcn)
for(i in 1:nrow(w_df)){

  temp <- row_df$w[i] *(raster(w_df$raster_names[i]) - w_df$bc[i])

if(i ==1){
  hybrid_raster <- temp
}else{
  hybrid_raster <- hybrid_raster + temp
}
i <- i+1
}

writeRaster(hybrid_raster, "hybrid_raster.tif", format="GTiff", overwrite=TRUE)
