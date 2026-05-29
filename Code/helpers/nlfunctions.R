###################################################################
# This library includes the following functions:
#   f.owner: obtains the Nevo (2001) ownership matrix from a vector of firm IDs
#   f.share: obtains conditional and unconditional shares given delta, nesting parameter
#   f.share2: obtains conditional shares given unconditional shares
#   f.dqddelta: obtains derivatives with respect to mean valuations as function of delta
#   f.condDqddelta: obtains conditional derivatives with respect to mean valuations as function of delta
#   f.elas: obtain matrix of price elasticities given delta, price coefficient
# All code is written for the 2-level nested logit, which collapses to a 1-level nested
#   logit with some values of the nesting parameters. The math is worked out in the accompanying
#   "Notes on the Nested Logit Demand Model" by Mansley, Miller, Ryan, Weinberg.
# The model is not fully general: the highest level demarcation is assumed to be the
#   outside good vs. all inside good.
# Updated August 20, 2019.
###################################################################

f.owner = function(fNum){
  num_prod = length(fNum)
  if (any(fNum!=mean(fNum)) ){
    owner = model.matrix(1:num_prod~-1+as.factor(fNum))
    owner=tcrossprod(owner)
  } else {
    owner = matrix(1,length(fNum),length(fNum))
  }
  return(owner)
}

f.share <- function(delta,sigma,group){
    emv   = exp(delta/(1-sigma[1]))
    I_hg  = (1-sigma[1])*log(aggregate(x=emv,by=list(group[,1]),FUN=sum)$x)
    temp  = unique(cbind(exp(I_hg[group[,1]]/(1-sigma[2])  ),group[,2]))
    I_g   = (1-sigma[2])*log(aggregate(x=temp[,1],by=list(temp[,2]),FUN=sum  )$x)
    I     = log(1+sum(exp(I_g)))
    sj_hg = emv / exp(I_hg[group[,1]]/(1-sigma[1]))
    sh_g  = exp(I_hg / (1-sigma[2])) / exp(I_g[unique(group)[,2]]/(1-sigma[2]))
    sg    = (exp(I_g) / exp(I) )
    sj    = sj_hg * sh_g[group[,1]] * sg[group[,2]]
    s0    = 1-sum(sj)
    return(list(sj=sj,sj_hg=sj_hg,sh_g=sh_g,sg=sg,s0=s0))
}

f.share2 = function(sj,group){
  sg <- aggregate(x=sj,by=list(group[,2]),FUN=sum)$x
  temp <- aggregate(x=sj,by=list(group[,1]),FUN=sum)$x
  sh_g <- unique(temp[group[,1]] / sg[group[,2]])
  sj_hg <- sj / temp[group[,1]]
  s0    = 1-sum(sj)
  return(list(sj=sj,sj_hg=sj_hg,sh_g=sh_g,sg=sg,s0=s0))
}

f.dqddelta = function(delta,sigma,group,cal,sj){
  dlist <- f.condDqddelta(delta=delta,sigma=sigma,group=group,cal=cal,sj=sj)
  if (cal==0){
    slist <- f.share(delta=delta,sigma=sigma,group=group)
    num_prod <- length(delta)
  } else {
    slist <- f.share2(sj=sj,group=group)
    num_prod <- length(sj)
  }
  dqddelta <- matrix(0,num_prod,num_prod)
  term1 <- dlist$dsj_hg * slist$sh_g[group[,1]]   * slist$sg[group[,2]]
  term2 <- slist$sj_hg  * dlist$dsh_g[group[,1],] * slist$sg[group[,2]]
  term3 <- slist$sj_hg  * slist$sh_g[group[,1]]   * dlist$dsg[group[,2],]
  dqddelta <- term1 + term2 + term3
  return(dqddelta)
}

f.condDqddelta = function(delta,sigma,group,cal,sj){
  if (cal==0){
    slist = f.share(delta=delta,sigma=sigma,group=group)
    num_prod <- length(delta)
  } else {
    slist <- f.share2(sj=sj,group=group)
    num_prod <- length(sj)
  }
  dsj_hg = f.owner(group[,1])
  dsj_hg[dsj_hg==TRUE]  = -(1/(1-sigma[1]))*tcrossprod(slist$sj_hg,slist$sj_hg)[dsj_hg==TRUE]
  diag(dsj_hg) = (1/(1-sigma[1]))*slist$sj_hg*(1-slist$sj_hg)
  sameSubGroup = t(matrix(group[,1],num_prod,length(unique(group[,1]))))==matrix(unique(group[,1]),length(unique(group[,1])),num_prod)
  temp = matrix(unique(group[,2]),length(unique(group[,1])),num_prod)
  sameGroup    = t(matrix(group[,2],num_prod,length(unique(group[,1]))))==temp[ order(temp[,1]), ]
  combine      = sameSubGroup+sameGroup
  dsh_g = matrix(0,dim(combine)[1],dim(combine)[2])
  dsh_g[combine==2] =  (1/(1-sigma[2]))*tcrossprod(slist$sh_g*(1-slist$sh_g),slist$sj_hg)[combine==2]
  dsh_g[combine==1] = -(1/(1-sigma[2]))*tcrossprod(slist$sh_g,slist$sj_hg*slist$sh_g[group[,1]])[combine==1]
  dsg=t(matrix(group[,2],num_prod,length(unique(group[,2]))))==matrix(unique(group[,2]),length(unique(group[,2])),num_prod)
  dsg[dsg==FALSE] = -tcrossprod(slist$sg,slist$sj)[dsg==FALSE]
  dsg[dsg==TRUE]  =  tcrossprod(1-slist$sg,slist$sj)[dsg==TRUE]
  return(list(dsj_hg=dsj_hg,dsh_g=dsh_g,dsg=dsg))
}

f.elas = function(price,alpha,delta,sigma,group){
  dqdp  = alpha*f.dqddelta(delta=delta,sigma=sigma,group=group,cal=0,sj=0)
  quant = f.share(delta=delta,sigma=sigma,group=group)$sj
  elas  = dqdp * tcrossprod(1/quant,price)
  return(elas)
}
