

SFA.TopFeaturesPerFac <- function(fmat,
                                  top_features)
{
  indices_mat=matrix(0,dim(fmat)[2],top_features);
  for(k in 1:dim(fmat)[2]){
    indices_mat[k,] <- order(abs(fmat[,k]), decreasing=TRUE)[1:top_features];
  }
  return(indices_mat)
}

