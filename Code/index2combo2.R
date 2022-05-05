#Gets the sigma-tuple corresponding to the Ith element of the list 0,1,...,choose(n,sigma)-1.  Supply I (index), n (number of objects), and sigma (number to draw)
args = commandArgs(trailingOnly = TRUE) 
#print(args)
for (i in 1:length(args)) { #check that inputs are of the form x = y
  if (regexpr('^[A-Za-z0-9_,.]{1,}(?=.*=)={1}?(?!.*=)',args[i],perl=TRUE)[1] == -1) { 
    stop(sprintf('Input %s must be of the form x = y',i)) 
  }
  if (regexpr('.+=[[:blank:]]*+([0-9]+)',args[i],perl=TRUE)[1] == 1) {
    eval(parse(text = gsub('(.+)=[[:blank:]]*+([0-9]+)',"\\1=as.numeric(\\2)",args[i],perl = TRUE)))
  }
}

if (!exists("allow.repeats")) { #positive integer for TRUE, 0 or negative integer for false 
  allow.repeats = FALSE
} else {
  allow.repeats = as.logical(as.numeric(allow.repeats))
}

if (allow.repeats) {
  n <- n + sigma - 1
}

# print(I)
# print(n)
# print(sigma)
# print(allow.repeats)

Ind_array <- array(I,dim=c(1,sigma+1)) #initate an array of the linear positions of the multi-index in sigma, sigma-1, ..., 1 dimensional tables
multiindex <- array(0,dim=c(1,sigma+1))

if (sigma == 1) {
  multiindex[2] <- I + 1 #multiindex corresponds to linear index when choosing a single element per draw
} else {
  for (i in seq(1,sigma,by = 1)) { #loop through elements of the multi-index i1, i2, ... , i_sigma
    ind_temp <- 0
    k <- 0
    while (ind_temp <= 0) {
      #subtract the k layers of the sigma-i+1-dimensional table not containing the multiindex
      ind_temp <- (choose(n - multiindex[i],sigma-i+1) - Ind_array[i]) - choose(n-(multiindex[i] + k),sigma-i+1)
      k <- k + 1
    }
    #print(c(i,k,ind_temp))
    Ind_array[,(i+1):(sigma+1)] <- choose(n-(multiindex[i] + k-1),sigma-(i+1)+1) - ind_temp  #value of the linear index in a table of one fewer dimensions
    multiindex[i+1] <- multiindex[i] + k - 1 #record the layer k at which the linear index occurs, in sigma-i+1 dimensions
  }
}
if (allow.repeats) {
  multiindex <- multiindex - seq_along(multiindex) + 1 + 1
}
#paste(lapply(multiindex[,2:(sigma+1)],function(X) as.character(X)),collapse=",")
#print(Ind_array)
cat(multiindex[,2:(sigma+1)],fill = TRUE) #print the tuple, starting from 0, corresponding to the multi-index 1,2,...,sigma
  