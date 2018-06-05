#' Calculates structure coefficients and commonality coefficients for a gaussian phylogenetic MCMCglmm.
#'
#' @param mod MCMCglmm model
#' @param partvars Predictors for the analysis
#' @param data data frame
#' @param type Either "marginal" or "conditional"
#' @param inv_phylo inverse phylogenetic covariance matrix calcd by MCMCglmm::inverseA
#' @param prior prior for MCMCglmm
#' @param nitt Number of iterations
#' @param burnin Burnin
#' @param thin Thinning interval
#'
#'
#' @return
#' Returns a list with two objects
#' \item{partR2} Mean, Mode and HPD interval of R2
#' \item{R2_chain} R2 MCMC chain
#'
#' @details Not very sophisticated at the moment, just works with this model.
#'
#' @examples
#' seals <- data(seal_data)
#' phylo <- data(seal_phylo)
#'
#' inv_phylo <- inverseA(seal_phylo, nodes="TIPS",scale=FALSE)$Ainv
#' prior<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))
#'
#' mod_gen <- MCMCglmm(TPM70_ratio ~ num_alleles_mean +  obs_het_mean, # prop_low_afs_mean,
#'           random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
#'           family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
#'           data=seal_data,nitt=1100,burnin=100,thin=10)
#'
#' out <- partR2(mod = mod_gen, partvars = c("num_alleles_mean", "obs_het_mean"),
#'               data = seal_data, inv_phylo = inv_phylo, prior = prior, nitt = 1100,
#'               burnin = 100, thin = 10)
#'
#' @export
#'
#'


# partition the R2 into variation unique and common to the predictors
partR2 <- function(mod, partvars = NULL, data = NULL, type = "marginal", inv_phylo = NULL, prior = NULL,
                    nitt=10000,burnin=1000, thin=50){

    if (is.null(partvars)) stop("partvars has to contain the variables for the commonality analysis")

    chain_length <- nrow(mod$VCV)

    # calculate structure coefficients -----------------------------

    calc_struc_coef <- function(mcmc_iter, partvar, mod){
        resp <- MCMCglmm::predict.MCMCglmm(mod, it = mcmc_iter)
        mod_mat <- as.matrix(mod$X)
        pred <- grep(partvar, colnames(mod_mat))
        out <- cor(resp,as.matrix(mod$X)[, pred])
        #out <- sqrt(summary(lm(resp~data[[partvar]]))$r.squared) # should be correct but still to check
        out
    }
    calc_struc_coef_full <- function(partvar, mod){
        all_coef <- sapply(1:chain_length, calc_struc_coef,  partvar, mod)
        class(all_coef) <- "mcmc"
        out <- data.frame("pred" = partvar,  medianSC = stats::median(all_coef), HPDinterval(all_coef)) #modeSC = MCMCglmm::posterior.mode(all_coef),
        row.names(out) <- NULL
        out
    }

    # calculate all structure coefficients
    all_SC <- do.call(rbind, lapply(partvars, calc_struc_coef_full, mod))

    # calculate common and unique R2 --------------------------------
    model_formula <- mod$Fixed$formula

    # just unique effects
    all_unique_R2 <- c()
    R2_full <- R2mcmc(mod)

    all_comb <- lapply(1:(length(partvars)), function(x) combn(partvars, x)) # length(partvars - 1)
    all_comb2 <- lapply(all_comb, function(x) apply(x, 2, function(x) out <- list(x)))
    all_comb3 <- unlist(unlist(all_comb2, recursive = FALSE), recursive = FALSE)

    CI <- 0.95
    calc_CI <- function(x) {
        out <- stats::quantile(x, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
    }

    calc_partR2 <- function(var_to_red, R2_full, data) {
        to_del <- paste(paste("-", var_to_red, sep= ""), collapse = " ")
        new_formula <- update.formula(model_formula, paste(". ~ . ", to_del, sep=""))

        mod_red <- MCMCglmm(new_formula,
            random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
            family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
            data=data,nitt=nitt,burnin=burnin,thin=thin)

        R2_red <- R2mcmc(mod_red)
        # R2 mdedian
        R2_diff <- R2_full$partR2$medianR2 - R2_red$partR2$medianR2
        # R2 CI
        R2_diff_vec <- calc_CI(R2_full$R2_chain - R2_red$R2_chain)
        out <- c("R2median" = R2_diff, "CIlow" = R2_diff_vec[1], "CIhigh" = R2_diff_vec[2])
    }

    R2_out <- as.data.frame(do.call(rbind, lapply(all_comb3, calc_partR2, R2_full, data)))
    names(R2_out) <- c("medianR2", "lower", "upper")
    # all_vars <- lapply(all_comb3, function(x) gsub('(a|e|i|o|u)', '', x))
    all_comb_names <- unlist(lapply(all_comb3, function(x) paste(x, collapse = " & ")))
    # get full model estimates
    R2_full_mod <- data.frame("combinations" = "full model", R2_full$partR2[c(2,4,5)]) # substract mean

    out <- rbind(R2_full_mod, data.frame("combinations" = all_comb_names, R2_out))
    out_full <- list("R2" = out, "SC" = all_SC)
}
