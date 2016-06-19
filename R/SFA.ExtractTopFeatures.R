

##############  SFA  extract top features  #########################

SFA.ExtractTopFeatures <- function(fmat,
                               top_features = 10,
                               options=c("min", "max"),
                               mult.annotate=FALSE)
{
    score <- lapply(1:dim(fmat)[2], function(n) {
      out <- t(apply(fmat, 1, function(x){
    #    y=x[n] *log(x[n]/x) + x - x[n];
        y = (x[n]-x)^2;
        return(y)
      }));
      return(out)
    })
  
  
  indices_mat=matrix(0,dim(fmat)[2],top_features);
  
  if(mult.annotate==FALSE){
  
  if(dim(fmat)[2]==2){
    for(k in 1:dim(fmat)[2])
    {
      temp_mat <- score[[k]][,-k];
      if(options=="min"){
        vec <- apply(as.matrix(temp_mat), 1, function(x) min(x))}
      if(options=="max"){
        vec <- apply(as.matrix(temp_mat), 1, function(x) max(x))}
      #vec <- temp_mat;
      ordered_kl <- order(vec, decreasing = TRUE);
      counter <- 1
      flag <- counter
      while(flag <= top_features)
      {
        if(counter > dim(fmat)[1]){
          indices_mat[k,(flag:top_features)]=NA;
          break
        }
        if(which.max(fmat[ordered_kl[counter],])==k){
          indices_mat[k, flag] <- ordered_kl[counter];
          flag <- flag + 1;
          counter <- counter + 1;}
        else{
          counter <- counter + 1;
        }
      }
    }
    
  }else{
    for(k in 1:dim(fmat)[2])
    {
      temp_mat <- score[[k]][,-k];
      if(options=="min"){
        vec <- apply(temp_mat, 1, function(x) min(x))}
      if(options=="max"){
        vec <- apply(temp_mat, 1, function(x) max(x))}
      
      ordered_kl <- order(vec, decreasing = TRUE);
      counter <- 1
      flag <- counter
      while(flag <= top_features)
      {
        if(counter > dim(fmat)[1]){
          indices_mat[k,(flag:top_features)]=NA;
          break
        }
        
        if(which.max(fmat[ordered_kl[counter],])==k){
          indices_mat[k, flag] <- ordered_kl[counter];
          flag <- flag + 1;
          counter <- counter + 1;
        } else {
          counter <- counter + 1;
        }
      }
    }
  }
  
  return(indices_mat);
  }
  
  if(mult.annotate==TRUE){
    
    if(dim(fmat)[2]==2){
      for(k in 1:dim(fmat)[2])
      {
        temp_mat <- score[[k]][,-k];
        if(options=="min"){
          vec <- apply(as.matrix(temp_mat), 1, function(x) min(x))}
        if(options=="max"){
          vec <- apply(as.matrix(temp_mat), 1, function(x) max(x))}
        #vec <- temp_mat;
        ordered_kl <- order(vec, decreasing = TRUE);
        counter <- 1
        flag <- counter
        while(flag <= top_features)
        {
          if(counter > dim(fmat)[1]){
            indices_mat[k,(flag:top_features)]=NA;
            break
          }
            indices_mat[k, flag] <- ordered_kl[counter];
            flag <- flag + 1;
            counter <- counter + 1;
        }
      }
    } else {
        for(k in 1:dim(fmat)[2])
        {
          temp_mat <- score[[k]][,-k];
          if(options=="min"){
            vec <- apply(temp_mat, 1, function(x) min(x))}
          if(options=="max"){
            vec <- apply(temp_mat, 1, function(x) max(x))}
          
          ordered_kl <- order(vec, decreasing = TRUE);
          counter <- 1
          flag <- counter
          while(flag <= top_features)
          {
            if(counter > dim(fmat)[1]){
              indices_mat[k,(flag:top_features)]=NA;
              break
            }
              indices_mat[k, flag] <- ordered_kl[counter];
              flag <- flag + 1;
              counter <- counter+1;
          }
        }
      }
    return(indices_mat)
  }
}

