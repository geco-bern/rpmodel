# solves quadratic equation given by y = a*x^2 + bx + c
# Based on MAESTRA equivalent (B. Medlyn)

# - larger root
QUADP <- function(A,B,C){
  
  if (any(is.na(c(A,B,C)))){
    return(NA)
  } else {
    if((B^2 - 4*A*C) < 0){
      warning("IMAGINARY ROOTS IN QUADRATIC")
      return(0)
    }
    
    if(identical(A,0)){
      if(identical(B,0)){
        return(0)
      } else {
        return(-C/B)
      }
    } else {
      return((- B + sqrt(B^2 - 4*A*C)) / (2*A))
    }
  }
  
}

# - smaller root
QUADM <- function(A,B,C){
  
  if (any(is.na(c(A,B,C)))){
    return(NA)
  } else {
    if((B^2 - 4*A*C) < 0){
      warning("IMAGINARY ROOTS IN QUADRATIC")
      return(0)
    }
    
    if(identical(A,0)){
      if(identical(B,0)){
        return(0)
      } else {
        return(-C/B)
      }
    } else {
      return((- B - sqrt(B^2 - 4*A*C)) / (2*A))
    } 
  }
  
}
