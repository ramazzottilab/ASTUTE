#' Perform the ASTUTE framework to estimate K expression signatures associated to a set of somatic alterations provided as input.
#'
#' @examples
#' set.seed(12345)
#' data(datasetExample)
#' resExample <- ASTUTE( alterations = datasetExample$alterations, 
#'                       expression = datasetExample$expression, 
#'                       regularization = TRUE, 
#'                       nboot = NA, 
#'                       num_processes = NA, 
#'                       verbose = FALSE )
#'
#' @title ASTUTE
#' @param alterations Input binary alterations matrix.
#' @param expression Input log2(x+1)-transformed normalized expression matrix.
#' @param regularization Boolean. If TRUE, perform regularization via elastic net with LASSO penalty.
#' @param nboot Number of bootstrap sampling (minimum 2) to be performed for a robust estimation of the expression signatures. 
#' If bootstrap does not need to be executed, this parameter needs to be set to either NA or NULL.
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode,
#' this parameter needs to be set to either NA or NULL.
#' @param verbose Boolean. If TRUE, print status messages.
#' @return A list with the discovered expression signatures. It includes 7 elements:
#'              input_data: list providing the input data (i.e., alterations and expression data);
#'              bootstrap: results of the inference by bootstrap (i.e., alpha alterations matrix, beta matrix, and intercept estimates);
#'              parameters: list with the paremeters used for the inference (i.e., regularization TRUE/FALSE and nboot);
#'              goodness_fit: goodness of fit estimated as the cosine similarity comparing observations and predictions;
#'              fold_changes: log2 fold changes estimates;
#'              pvalues: p-values estimates;
#'              qvalues: p-values estimates corrected for false discovery.
#' @export ASTUTE
#' @importFrom glmnet glmnet cv.glmnet 
#' @importFrom lsa cosine
#' @importFrom parallel clusterExport clusterEvalQ clusterSetRNGStream detectCores makeCluster parLapply stopCluster
#' @importFrom stats coef p.adjust runif t.test
#'
ASTUTE <- function( alterations, expression, regularization = TRUE, nboot = 100, num_processes = Inf, verbose = TRUE ) {

    # perform the analysis
    if(is.na(nboot) || is.null(nboot)) { # do not perform bootstrap

        # estimate the mutational signatures
        res <- signaturesEstimation(alterations = alterations, expression = expression, 
                regularization = regularization, verbose = verbose)

        # get the results of the inference
        alpha <- res$inference$alpha
        beta <- res$inference$beta
        baseline <- res$inference$intercept

        # evaluate the goodness of fit
        goodness_fit <- rep(NA, nrow(alpha))
        names(goodness_fit) <- rownames(alpha)
        expression <- expression[,colnames(beta)]
        predicted_exp <- t(t((alpha %*% beta)) + baseline)
        for(i in names(goodness_fit)) {
            obs_s <- as.numeric(expression[i,])
            pred_s <- as.numeric(predicted_exp[i,])
            goodness_fit[i] <- as.numeric(cosine(obs_s,pred_s))
        }

        # compute fold changes
        fold_changes <- matrix(NA, nrow = ncol(beta), ncol = nrow(beta))
        rownames(fold_changes) <- colnames(beta)
        colnames(fold_changes) <- rownames(beta)
        pvalues <- fold_changes
        qvalues <- pvalues
        intercept <- rep(NA, ncol(beta))
        names(intercept) <- colnames(beta)
        for(i in 1:nrow(fold_changes)) {
            baseline_estimate <- as.numeric(baseline[i])
            effect_estimate <- baseline_estimate + as.numeric(beta[,i])
            baseline_estimate_abs <- ((2^baseline_estimate)-1)
            effect_estimate_abs <- ((2^effect_estimate)-1)
            min_val <- min(baseline_estimate_abs,effect_estimate_abs)
            if(min_val < 0) { # expression counts must be at minimum 0
                baseline_estimate_abs <- baseline_estimate_abs + abs(min_val)
                effect_estimate_abs <- effect_estimate_abs + abs(min_val)
            }
            intercept[i] <- baseline_estimate_abs
            fold_changes[i,] <- log(((effect_estimate_abs+1)/(baseline_estimate_abs+1)),base=2)
        }
        fold_changes <- cbind(intercept,fold_changes)
        colnames(fold_changes)[1] <- "BASELINE"

    }
    else { # perform bootstrap

        # estimate the mutational signatures
        res <- signaturesBootstrap(alterations = alterations, expression = expression, 
                regularization = regularization, nboot = nboot, num_processes = num_processes, 
                verbose = verbose)

        # get the results for each bootstrap estimate
        alpha <- res$input_data$alterations
        valid_genes <- colnames(res$bootstrap$beta[[1]])
        for(i in 2:length(res$bootstrap$beta)) {
            valid_genes <- intersect(valid_genes,colnames(res$bootstrap$beta[[i]]))
        }
        beta <- matrix(0, nrow = ncol(alpha), ncol = length(valid_genes))
        rownames(beta) <- colnames(alpha)
        colnames(beta) <- sort(unique(valid_genes))
        baseline <- rep(0, length(valid_genes))
        names(baseline) <- colnames(beta)
        for(i in 1:length(res$bootstrap$beta)) {
            beta <- beta + res$bootstrap$beta[[i]][,colnames(beta)]
            baseline <- baseline + res$bootstrap$intercept[[i]][names(baseline)]
        }
        beta <- beta / length(res$bootstrap$beta)
        baseline <- baseline / length(res$bootstrap$beta)

        # evaluate the goodness of fit
        goodness_fit <- rep(NA, nrow(alpha))
        names(goodness_fit) <- rownames(alpha)
        expression <- expression[,colnames(beta)]
        predicted_exp <- t(t((alpha %*% beta)) + baseline)
        for(i in names(goodness_fit)) {
            obs_s <- as.numeric(expression[i,])
            pred_s <- as.numeric(predicted_exp[i,])
            goodness_fit[i] <- as.numeric(cosine(obs_s,pred_s))
        }

        # compute fold changes
        if (verbose) {
            message("Estimating fold change for each gene...", "\n")
        }
        fold_changes <- matrix(0, nrow = ncol(beta), ncol = nrow(beta))
        rownames(fold_changes) <- colnames(beta)
        colnames(fold_changes) <- rownames(beta)
        pvalues <- matrix(NA, nrow = ncol(beta), ncol = nrow(beta))
        rownames(pvalues) <- colnames(beta)
        colnames(pvalues) <- rownames(beta)
        intercept <- rep(0, ncol(beta))
        names(intercept) <- colnames(beta)
        for(i in 1:nrow(fold_changes)) {
            curr_pvalues <- matrix(NA, nrow = length(res$bootstrap$beta), ncol = nrow(beta))
            rownames(curr_pvalues) <- 1:nrow(curr_pvalues)
            colnames(curr_pvalues) <- rownames(beta)
            for(j in 1:length(res$bootstrap$beta)) {
                baseline_estimate <- as.numeric(res$bootstrap$intercept[[j]][names(baseline)[i]])
                effect_estimate <- baseline_estimate + as.numeric(res$bootstrap$beta[[j]][,colnames(beta)[i]])
                baseline_estimate_abs <- ((2^baseline_estimate)-1)
                effect_estimate_abs <- ((2^effect_estimate)-1)
                min_val <- min(baseline_estimate_abs,effect_estimate_abs)
                if(min_val < 0) { # expression counts must be at minimum 0
                    baseline_estimate_abs <- baseline_estimate_abs + abs(min_val)
                    effect_estimate_abs <- effect_estimate_abs + abs(min_val)
                }
                intercept[i] <- intercept[i] + baseline_estimate_abs
                curr_fc <- log(((effect_estimate_abs+1)/(baseline_estimate_abs+1)),base=2)
                fold_changes[i,] <- fold_changes[i,] + curr_fc
                curr_pvalues[j,] <- curr_fc
            }
            for(j in colnames(curr_pvalues)) {
                curr_p <- t.test(x = as.numeric(curr_pvalues[,j]), mu = 0.00, alternative = "two.side")$p.value
                pvalues[i,j] <- curr_p
            }
            if (verbose) {
                message((i/nrow(fold_changes)), "\n")
            }
        }
        fold_changes <- cbind(intercept,fold_changes)
        colnames(fold_changes)[1] <- "BASELINE"
        fold_changes <- fold_changes / length(res$bootstrap$beta)
        qvalues <- pvalues
        for(j in colnames(pvalues)) {
            qvalues[,j] <- p.adjust(p = pvalues[,j], method = "fdr")
        }

    }
    rm(res)
    gc(verbose = FALSE)

    # save the results
    input_data <- list(alterations = alterations, expression = expression)
    inference <- list(alpha = alpha, beta = beta, intercept = baseline)
    parameters <- list(regularization = regularization, nboot = nboot)
    results <- list(input_data = input_data, inference = inference, parameters = parameters, goodness_fit = goodness_fit, 
        fold_changes = fold_changes, pvalues = pvalues, qvalues = qvalues)

    # return the results
    return(results)

}

# Perform a robust estimation by bootstrap of K expression signatures associated to a set of somatic alterations provided as input.
#
# @param alterations Input binary alterations matrix.
# @param expression Input log2(x+1)-transformed normalized expression matrix.
# @param regularization Boolean. If TRUE, perform regularization via elastic net with LASSO penalty.
# @param nboot Number of bootstrap sampling (minimum 2) to be performed for a robust estimation of the expression signatures.
# @param num_processes Number of processes to be used during parallel execution. To execute in single process mode,
# this parameter needs to be set to either NA or NULL.
# @param verbose Boolean. If TRUE, print status messages.
# @return A list with the discovered expression signatures. It includes 3 elements:
#              input_data: list providing the input data (i.e., alterations and expression data);
#              bootstrap: results of the inference by bootstrap (i.e., alpha alterations matrix, beta matrix, and intercept estimates);
#              parameters: list with the paremeters used for the inference (i.e., regularization TRUE/FALSE and nboot).
#
signaturesBootstrap <- function( alterations, expression, regularization = TRUE, nboot = 100, num_processes = Inf, verbose = TRUE ) {

    # check the input parameters
    alterations <- as.matrix(alterations)
    expression <- as.matrix(expression)
    if (nboot < 2) {
        warning("The minimum value of nboot must be 2...")
        nboot <- 2
    }

    # setting up parallel execution
    parallel <- NULL
    close_parallel <- FALSE
    if (is.na(num_processes) || is.null(num_processes)) {
        num_processes <- 1
    } else if (num_processes == Inf) {
        cores <- as.integer((detectCores() - 1))
        num_processes <- min(cores, nboot)
    } else {
        num_processes <- min(num_processes, nboot)
    }
    if(num_processes > 1) {
        parallel <- makeCluster(num_processes, outfile = "")
        close_parallel <- TRUE
        res_clusterEvalQ <- clusterEvalQ(parallel, library("glmnet", warn.conflicts = FALSE,
            quietly = TRUE, verbose = FALSE))
        clusterExport(parallel, varlist = c("signaturesEstimation"), envir = environment())
        clusterExport(parallel, varlist = c("alterations", "expression"), envir = environment())
        clusterExport(parallel, varlist = c("regularization", "nboot"), envir = environment())
        clusterSetRNGStream(parallel, iseed = round(runif(1) * 1e+05))
    }

    if (verbose) {
        message("Performing expression signatures estimation...", "\n")
        if (num_processes > 1) {
            message("Executing ", num_processes, " processes in parallel...", "\n")
        }
    }

    # perform the inference
    if (num_processes == 1) { # performing the inference sequentially
        results <- list()
        for (iteration in seq_len(nboot)) {
            if (verbose) {
                message("Performing bootstrap iteration ", iteration, " out of ",
                    nboot, "...", "\n")
            }
            curr_entries <- sample(x = 1:nrow(alterations), size = nrow(alterations), replace = TRUE)
            boot_alterations <- alterations[curr_entries,,drop=FALSE]
            boot_expression <- expression[curr_entries,,drop=FALSE]
            res <- signaturesEstimation(alterations = boot_alterations, expression = boot_expression, 
                regularization = regularization, verbose = FALSE)
            res <- res$inference
            results[[iteration]] <- res
        }
    } else { # performing the inference in parallel
        results <- parLapply(parallel, seq_len(nboot), function(iteration) {
            if (verbose) {
                message("Performing bootstrap iteration ", iteration, " out of ",
                    nboot, "...", "\n")
            }
            curr_entries <- sample(x = 1:nrow(alterations), size = nrow(alterations), replace = TRUE)
            boot_alterations <- alterations[curr_entries,,drop=FALSE]
            boot_expression <- expression[curr_entries,,drop=FALSE]
            res <- signaturesEstimation(alterations = boot_alterations, expression = boot_expression, 
                regularization = regularization, verbose = FALSE)
            res <- res$inference
            return(res)
        })
    }

    # get the results
    alpha <- list()
    beta <- list()
    intercept <- list()
    for(i in 1:length(results)) {
        alpha[[i]] <- results[[i]][["alpha"]]
        beta[[i]] <- results[[i]][["beta"]]
        intercept[[i]] <- results[[i]][["intercept"]]
    }

    # save the results
    input_data <- list(alterations = alterations, expression = expression)
    bootstrap <- list(alpha = alpha, beta = beta, intercept = intercept)
    parameters <- list(regularization = regularization, nboot = nboot)
    results <- list(input_data = input_data, bootstrap = bootstrap, parameters = parameters)

    # close parallel
    if (close_parallel) {
        stopCluster(parallel)
    }
    rm(close_parallel)
    rm(parallel)
    gc(verbose = FALSE)

    # return the results
    return(results)

}

# Perform the estimation of K expression signatures associated to a set of somatic alterations provided as input.
#
# @param alterations Input binary alterations matrix.
# @param expression Input log2(x+1)-transformed normalized expression matrix.
# @param regularization Boolean. If TRUE, perform regularization via elastic net with LASSO penalty.
# @param verbose Boolean. If TRUE, print status messages.
# @return A list with the discovered expression signatures. It includes 3 elements:
#              input_data: list providing the input data (i.e., alterations and expression data);
#              inference: results of the inference (i.e., alpha alterations matrix, beta matrix, and intercept estimates);
#              parameters: list with the paremeters used for the inference (i.e., regularization TRUE/FALSE).
#
signaturesEstimation <- function( alterations, expression, regularization = TRUE, verbose = TRUE ) {
    
    # performing the analysis
    if (verbose) {
        message("Performing expression signatures estimation...", "\n")
    }

    # initialize alpha with alterations data
    # alpha is an N x K binary matrix, where rows (N) represent the patients and columns (K) represent the K alterations
    alpha <- as.matrix(alterations)

    # initialize beta with an empty matrix
    # beta is a K x M matrix, where rows represent expression signatures for K alterations, and columns represent the 
    # M considered genes (usually a selected set of genes of interest)
    beta <- matrix(NA, nrow = ncol(alpha), ncol = ncol(expression))
    rownames(beta) <- colnames(alpha)
    colnames(beta) <- colnames(expression)

    # initialize intercept estimates
    intercept <- rep(NA,ncol(beta))
    names(intercept) <- colnames(beta)

    # perform expression signatures estimation
    if (regularization) { # in this case, fit a generalized linear model via elastic net with LASSO penalty to estimate beta
        for (i in seq_len(ncol(beta))) {
            x <- as.numeric(expression[,i])
            beta[,i] <- tryCatch({ # try to perform fit
                fit <- cv.glmnet(alpha, x, type.measure = "mse", nfolds = 10, nlambda = 100, family = "gaussian")
                fit <- as.numeric(coef(fit,s=fit$lambda.min))
                intercept[i] <- fit[1]
                res <- fit[-1]
                res
            }, error = function(err) { # if the fit fails, discard the current gene
                if (verbose) {
                    message("Fit failed for gene ",colnames(beta)[i],": discarding the gene.","\n")
                }
                return(NA)
            })
        }
    } else { # in this case, fit a generalized linear model without regularization to estimate beta
        for (i in seq_len(ncol(beta))) {
            x <- as.numeric(expression[,i])
            beta[,i] <- tryCatch({ # try to perform fit
                fit <- glmnet(alpha, x, lambda = 0, family = "gaussian")
                fit <- as.numeric(coef(fit,s=fit$lambda))
                intercept[i] <- fit[1]
                res <- fit[-1]
                res
            }, error = function(err) { # if the fit fails, discard the current gene
                if (verbose) {
                    message("Fit failed for gene ",colnames(beta)[i],": discarding the gene.","\n")
                }
                return(NA)
            })
        }
    }

    # remove any invalid gene
    invalid_genes <- which(is.na(beta),arr.ind=TRUE)
    if (nrow(invalid_genes)>0) {
        invalid_genes <- unique(invalid_genes[,"col"])
        expression <- expression[,-invalid_genes]
        beta <- beta[,-invalid_genes]
        intercept <- intercept[-invalid_genes]
    }

    # return the assigned signatures
    input_data <- list(alterations = alterations, expression = expression)
    inference <- list(alpha = alpha, beta = beta, intercept = intercept)
    parameters <- list(regularization = regularization)
    results <- list(input_data = input_data, inference = inference, parameters = parameters)
    return(results)

}
