bIErrorCheck <-
function(x,k,y,Z=NULL,nCov=0,plot0=FALSE,nIts=50000,nBurn=10000,propGamma1Pi=0.7,cx=10,pix=0.5){
  
  # check x,k, and y are vectors of the same length
  if(!(length(x)==length(k) & length(k)==length(y))){stop("Error: Vectors x, k and y must be vectors of the same length.")}
  
  # Check if Z is null
  # Check that nCov is not 0 if z is not NULL
  # Check that the number of rows in Z is equal to the length of x (and thus k and y)
  # Check that there are either 0, 1, or 2 covariates
  # Check that the number of covariates is equal to the number of columns of Z
  if(!is.null(Z)){
    if(nCov==0){stop("Error: If nCov is zero, then Z must be NULL.")}
    if(!(nCov==0 | nCov==1 | nCov==2)){stop("Error: nCov must be 0, 1 or 2.")}
    if(nCov!=ncol(Z)){stop("Error: nCov must equal the number of elements in Z.")}
    if(nrow(Z)!=length(x)){stop("Error: The number of rows in Z should equal the length of vector x.")}
  }else{
    if(nCov!=0){stop("Error: If Z is NULL, then nCov must be 0.")}
  }

  # Check propGamma1Pi is in 0 to 1
  if(!(propGamma1Pi<= 1 & propGamma1Pi>=0)){stop("Error: propGamma1pi must be between 0 and 1.")}
  
  # Check pix is in 0 to 1
  if(!(pix<= 1 & pix>=0)){stop("Error: pix must be between 0 and 1.")}
  
  # Check cx > 1
  if(!(cx>0.5)){stop("Error: cx must be greater than 0.5.")}
  
  # Check nIts and nBurn are positive integers
  if(!(nIts%%1==0 & nIts>0)){stop("Error: nIts must be a positive integer.")}
  if(!(nBurn%%1==0 & nBurn>0)){stop("Error: nBurn must be a positive integer.")}
  if(nIts<nBurn){stop("Error: nIts must be >= nBurn.")}

  ###################
}
