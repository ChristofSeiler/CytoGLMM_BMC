# function for cluster experiments
experiment = function(seed = 1,                 # random seed
                      n_markers = 10,           # number of markers
                      n_true = 2,               # number of differential markers
                      bias = 0.5,               # biased induced by pipe confounding
                      n_samples = 16,            # donors per group
                      paired = TRUE,            # paired experimental design
                      n_cells = 1000,           # per sample
                      pi = 0.05,                # zero inflation probability
                      beta_treatment = log(15), # log of mean in treatment
                      beta_control = log(10),   # log of mean in control
                      rho_b = 0.1,              # cell level correlation
                      rho_u = 0.1,              # donor level correlation
                      sigma_b = 1,              # cell level standard deviation
                      sigma_u = 1,              # donor level standard deviation
                      fdr = 0.05                # FDR control
) {
  
  if((n_samples %% 2) != 0) stop("n_samples need to be even")
  
  # simulation parameters
  set.seed(seed)
  
  if(paired) {
    n_donors = n_samples/2
  } else {
    n_donors = n_samples
  }
  
  # patient information
  donor = rep(1:n_donors, each = n_cells)
  if(paired) {
    donor %<>% rep(2)
    condition = c(
      rep("treatment", length(donor)/2),
      rep("control", length(donor)/2)
    )
  } else {
    condition = c(
      rep("treatment", length(donor)/2),
      rep("control", length(donor)/2)
    )
  }
  df = tibble(donor, condition)
  
  # generate protein counts
  protein_names = paste0("m", str_pad(1:n_markers, width = 2, pad = "0"))
  rcov = function(rho, sigma) {
    corr = rho^toeplitz(0:(n_markers-1))
    sigma_vec = rep(sigma, n_markers)
    diag(sigma_vec) %*% corr %*% diag(sigma_vec)
  }
  rcov_block = function(rho, sigma) {
    corr = diag(1, nrow = n_markers)
    corr_act = rho^toeplitz(0:(n_true-1))
    corr_notact = rho^toeplitz(0:(n_markers-n_true-1))
    corr[1:n_true,1:n_true] = corr_act
    corr[(n_true+1):n_markers,(n_true+1):n_markers] = corr_notact
    sigma_vec = rep(sigma, n_markers)
    diag(sigma_vec) %*% corr %*% diag(sigma_vec)
  }
  Sigma_b = rcov(rho_b, sigma_b) # cell level variability
  Sigma_u = rcov(rho_u, sigma_u) # donor level variability
  b = mvrnorm(n = nrow(df), mu = rep(0, n_markers), Sigma_b)
  u = mvrnorm(n = n_donors, mu = rep(0, n_markers), Sigma_u)
  u = u[donor, ]
  beta = matrix(beta_control, nrow = nrow(b), ncol = n_markers)
  beta[,1:n_true] = ifelse(condition == "treatment", beta_treatment, beta_control)
  log_lambda = beta + b + u
  # pipe confounding
  if(bias > 0) {
    idx = (n_true+1):(n_true*2)
    log_lambda[,idx] = log_lambda[,idx] + bias*log_lambda[,1:n_true]
  }
  lambda = exp(log_lambda)
  y = rpois(length(lambda), lambda)
  dim(y) = dim(lambda)
  I = rbinom(n = length(y), size = 1, prob = pi)
  I = matrix(I, ncol = n_markers)
  y = y*(1-I)
  colnames(y) = protein_names
  df %<>% bind_cols(as_tibble(y))
  df$condition %<>% factor(levels = c("control",
                                      "treatment"))
  df %<>% dplyr::mutate_at(protein_names, function(x) asinh(x/5))
  
  # # debugging: marginal expressions
  # ggplot(df, aes(m01, fill = condition)) +
  #   geom_histogram(binwidth = 0.25, alpha = 0.5, position = "identity")
  # ggplot(df, aes(m02, fill = condition)) +
  #   geom_histogram(binwidth = 0.25, alpha = 0.5, position = "identity")
  # ggplot(df, aes(m01, m02, color = condition)) +
  #   geom_density_2d()
  
  # compute false discovery proportion and power
  compute_stats = function(fit, method) {
    non_null = protein_names[1:n_true]
    null = protein_names[(n_true+1):n_markers]
    tb = summary(fit, method = method)
    tb %<>% mutate(discovery = ifelse(pvalues_adj <= fdr, TRUE, FALSE))
    a = tb %>% filter(protein_name %in% null & discovery) %>% nrow
    R = tb %>% filter(discovery) %>% nrow
    fdp = a/R
    if(a == 0 & R == 0) fdp = 0
    b = tb %>% filter(protein_name %in% non_null & discovery) %>% nrow
    N1 = n_true
    power = b/N1
    list(fdp = fdp, power = power)
  }
  
  # mixed model
  glmm_fit = CytoGLMM::cytoglmm(df, 
                                protein_names = protein_names,
                                condition = "condition", 
                                group = "donor",
                                num_cores = 1)
  
  # bootstrap
  glm_fit = CytoGLMM::cytoglm(df,
                              num_boot = 1000,
                              num_cores = 1,
                              protein_names = protein_names,
                              condition = "condition", 
                              group = "donor")
  
  # Citrus
  compute_stats_citrus = function() {
    
    # fit lasso on median marker expressions
    df_median = df %>%
      group_by(donor, condition) %>%
      summarise_at(protein_names, median)
    X = as.matrix(df_median[,protein_names])
    cond = as.numeric(df_median$condition)-1
    cv_fits = cv.glmnet(x = X, y = cond, family = "binomial", nfolds = nrow(X))
    
    # # selet the model with fdr
    # # define lambda sequence
    # lambda_sequence = 10^seq(-5, 1, 0.05)
    # cv_fits = cv.glmnet(x = X, y = cond, family = "binomial", lambda = lambda_sequence, nfolds = nrow(X))
    # # workaround for bug in glmnet
    # perms = 100
    # zeros = list()
    # for(i in 1:perms) {
    #   tryCatch(
    #     error = function(cnd) {
    #       paste0("--", conditionMessage(cnd), "--")
    #     },
    #     zeros[[i]] <- cv.glmnet(x = X, y = sample(cond), family = "binomial", lambda = lambda_sequence, nfolds = nrow(X))$nzero
    #   )
    # }
    # zeros %<>% bind_rows()
    # median_zeros = apply(zeros, MARGIN = 2, function(x) median(x, na.rm = TRUE))
    # fdr_hat = median_zeros/cv_fits$nzero
    # model_id = which.max(fdr_hat < fdr)
    
    # calculate stats
    non_null = protein_names[1:n_true]
    null = protein_names[(n_true+1):n_markers]
    tb = tibble(protein_name = protein_names, beta = 0)
    # if(cv_fits$nzero[model_id] > 0) {
    #   fit = glmnet(x = X, y = cond, family = "binomial", lambda = cv_fits$lambda[model_id])
    #   tb %<>% mutate(beta = as.numeric(fit$beta))
    # }
    fit = glmnet(x = X, y = cond, family = "binomial", lambda = cv_fits$lambda.min)
    tb %<>% mutate(beta = as.numeric(fit$beta))
    tb %<>% mutate(discovery = ifelse(beta != 0, TRUE, FALSE))
    a = tb %>% filter(protein_name %in% null & discovery) %>% nrow
    R = tb %>% filter(discovery) %>% nrow
    fdp = 0
    if(R > 0) fdp = a/R
    b = tb %>% filter(protein_name %in% non_null & discovery) %>% nrow
    N1 = n_true
    power = b/N1
    list(fdp = fdp, power = power)
    
  }

  # return a table with one row per method
  tb_info = tibble(
    seed, n_markers, n_true, bias, n_samples, paired, n_cells, pi, 
    beta_treatment, beta_control, rho_b, rho_u, sigma_b, sigma_u, fdr
  )
  bind_rows(
    bind_cols(tb_info, method = "GLMM-BH", compute_stats(glmm_fit, method = "BH")),
    bind_cols(tb_info, method = "GLMM-BY", compute_stats(glmm_fit, method = "BY")),
    bind_cols(tb_info, method = "GLM-BH", compute_stats( glm_fit, method = "BH")),
    bind_cols(tb_info, method = "GLM-BY", compute_stats( glm_fit, method = "BY")),
    bind_cols(tb_info, method = "Citrus", compute_stats_citrus())
  )
  
}
