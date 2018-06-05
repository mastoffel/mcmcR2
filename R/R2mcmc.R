#' Calculates marginal R2 with CI for MCMCglmm models.
#'
#' @param mod MCMCglmm model
#' @param type At the moment just "marginal"
#' @param family At the moment just "gaussian"
#'
#'
#' @return
#' Returns a list with two objects
#' \item{partR2} Mean, Mode and HPD interval of R2
#' \item{R2_chain} R2 MCMC chain
#'
#' @details At the moment, the actualy function arguments are defined in the function, except for formula.
#'
#' @examples
#' seals <- data(seal_data)
#' phylo <- data(seal_phylo)
#'
#' inv_phylo <- inverseA(seal_phylo, nodes="TIPS",scale=FALSE)$Ainv
#' prior<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))
#'
#' mod_gen <- MCMCglmm(TPM70_ratio ~ num_alleles_mean +  obs_het_mean + prop_low_afs_mean,
#'           random=~tip_label, nodes = "TIPS",
#'           family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
#'           data=seal_data,nitt=1100,burnin=100,thin=10)
#'
#' out <- R2mcmc(mod_gen, type = "marginal")
#' out <- R2mcmc(mod_gen, type = "conditional")
#' out$partR2
#'
#'
#' @export
#'
#'

R2mcmc <- function(mod, type = "marginal", family = "gaussian"){
   # if (type != "marginal") stop("At the moment, there is just the marginal R2")
    if (family != "gaussian") stop("At the moment just gaussian errors are supported")
    # Shinichis answer on Researchgate
    mVarF <- var(as.vector(apply(mod$Sol,2,mean) %*% t(as.matrix(mod$X))))  #t(as.matrix(mod$X))))
    # R2 <- mVarF/(mVarF+sum(apply(mod$VCV,2,mean)))

    # alternative with crebile intervals
    mcmc_chain_length <- nrow(mod$VCV)
    vmVarF<-numeric(mcmc_chain_length)

    vmVarF <- vapply(1:mcmc_chain_length, function(x) out <- var(as.vector(mod$Sol[x,] %*% t(as.matrix(mod$X)))), #t(as.matrix(mod$X))))
        FUN.VALUE = numeric(length(mcmc_chain_length)))

    if (type == "marginal"){
        R2m <- vmVarF/(vmVarF+ rowSums(mod$VCV)) # include here all random effects plus errors
    } else if (type == "conditional"){
        # Delete last column in VCV which contains residuals (observational level variance)
        R2m <- (vmVarF + rowSums(mod$VCV[, -ncol(mod$VCV), drop=FALSE]))/(vmVarF+ rowSums(mod$VCV))

    }
    outR2m <- R2m
    class(R2m) <- "mcmc"

    # R2m<-vmVarF/(vmVarF+mod$VCV[,1]+mod$VCV[,2])
    #mean(R2m)
    #posterior.mode(R2m)
    #data.frame(HPDinterval(R2m), row.names = NULL)

    out <- list("partR2" = data.frame("meanR2" = mean(R2m),"medianR2" = median(R2m), "modeR2" = posterior.mode(R2m), data.frame(HPDinterval(R2m)), row.names = NULL),
        "R2m_chain" = outR2m)
}
