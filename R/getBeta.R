
getBeta = function(
y, # nomalised gene expression
gplvm, # output of gplvm 
g0,  # expanded genotype vector (n cells x 1)
af, # allele frequency of g0
deltaG = 0.001
){
    g = (g0-2*af)/sqrt(2*af*(1-af))
    
    omega2 = gplvm$Param$omega2
    theta = gplvm$Param$theta
    
    tA = rbind(gplvm$Param$Alpha,1)[,1]
    tW = cbind(gplvm$Data$W,gplvm$Data$Z%*%gplvm$Param$zeta)
    yt = y - tW%*%tA
    
    # kernel for g
    K1 = updateKernelSE(gplvm$Param$Xi[[2]][,-c(3, 4)],
                        reduceDim1(gplvm$Param$Ta[[2]][,-c(3, 4)],4),
                        gplvm$Param$rho[[2]][-c(3, 4)])
    
    tZ = cbind(gplvm$Data$Z, g*K1$Knm, g)
    Dinv = dbind(gplvm$Param$Dinv, K1$K/theta/deltaG, 1/deltaG)
    Phiinv = Dinv + t(tZ/omega2)%*%tZ
    
    ZtOinvy = c(t(tZ/omega2)%*%yt)
    beta = Solve(Phiinv,ZtOinvy)
    
    betas = c((cbind(K1$Knm,1)/sqrt(2*af*(1-af)))%*%rev(rev(beta)[1:(ncol(K1$Knm)+1)]))
    betas
}




