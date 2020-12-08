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
  
  # regression calibration
  n = 2*n_donors
  k = n_cells
  df_splits = df %>% group_by(donor, condition) %>% group_split()
  Ws = lapply(df_splits, function(df) dplyr::select(df, protein_names))
  Sigma_w = lapply(Ws, function(W) cov(W)) %>% simplify2array()
  Sigma_uu = apply(Sigma_w, c(1,2), sum)/n
  Wbar = lapply(Ws, function(W) colMeans(W)) %>% simplify2array() %>% t()
  Sigma_wbar = cov(Wbar)
  mu_x = mu_w = colMeans(Wbar)
  X = apply(Wbar, 1, function(wbar) {
    mu_x + (Sigma_wbar - Sigma_uu) %*% ginv(Sigma_wbar) %*% (wbar - mu_w)
  }) %>% t
  colnames(X) = protein_names
  df_info = df %>% group_by(donor, condition) %>% tally()
  df_x = bind_cols(df_info, as_tibble(X))
  
  # permutation mixel model
  n_perm = 1000
  glm_wrap = function(df) {
    df %<>% mutate(y = ifelse(condition == "control", 0, 1))
    formula_str = paste("y ~", paste(protein_names, collapse = " + "))
    glm(formula_str, family = "binomial", data = df)
  }
  cytopermute = function() {
    glm_permute = function(i) {
      df_x$condition %<>% sample()
      glm_fit = glm_wrap(df_x)
      tibble(protein_name = protein_names,
             coef = glm_fit$coefficients[protein_names],
             perm = i)
    }
    tb_coef = lapply(1:n_perm, glm_permute) %>% bind_rows()
    tb_coef %<>% pivot_wider(names_from = protein_name, values_from = coef)
    glm_fit = glm_wrap(df_x)
    out = NULL
    out$tb_coef = tb_coef
    out$glm_fit = glm_fit
    class(out) = "cytopermute"
    out
  }
  summary.cytopermute = function(fit, method = "BH") {
    coef_obsv = fit$glm_fit$coefficients[protein_names]
    pvalues_unadj = sapply(protein_names, function(x) {
      (sum(abs(fit$tb_coef[x]) > abs(coef_obsv[x])) + 1) / (nrow(fit$tb_coef) + 1)
    })
    df_pvalues = tibble(protein_name = protein_names, 
                        pvalues_unadj = pvalues_unadj)
    df_pvalues %<>% mutate(pvalues_adj = p.adjust(pvalues_unadj,
                                                  method = method))
    df_pvalues = df_pvalues[order(df_pvalues$pvalues_unadj),]
    df_pvalues
  }
  glm_fit_perm = cytopermute()
  
  # bootstrap
  glm_fit = CytoGLMM::cytoglm(df,
                              num_boot = 1000,
                              num_cores = 1,
                              protein_names = protein_names,
                              condition = "condition", 
                              group = "donor")
  
  # # GLM on median markers
  # df_median = df %>% 
  #   group_by(donor, condition) %>% 
  #   summarize_at(protein_names, median)
  # df_tally  = df %>% 
  #   group_by(donor, condition) %>%
  #   tally()
  # df_median %<>% add_column(n = df_tally$n)
  # # select some proteins
  # foci_select = foci(Y = as.numeric(df_median$condition), 
  #                    X = df_median[,protein_names], 
  #                    num_features = round(nrow(df_median)/2), 
  #                    stop = FALSE, numCores = 1)
  # protein_selected = foci_select$selectedVar$names
  # formula_str = paste("condition ~", paste(protein_selected, collapse = " + "))
  # glm_median_fit = glm(formula = formula_str, family = binomial(), data = df_median, 
  #                      weights = df_median$n)
  
  # return a table with one row per method
  tb_info = tibble(
    seed, n_markers, n_true, bias, n_samples, paired, n_cells, pi, 
    beta_treatment, beta_control, rho_b, rho_u, sigma_b, sigma_u, fdr
  )
  bind_rows(
    bind_cols(tb_info, method = "GLMM-BH", compute_stats(glmm_fit, method = "BH")),
    bind_cols(tb_info, method = "GLMM-BY", compute_stats(glmm_fit, method = "BY")),
    bind_cols(tb_info, method = "GLM-PERM-BH", compute_stats( glm_fit_perm, method = "BH")),
    bind_cols(tb_info, method = "GLM-PERM-BY", compute_stats( glm_fit_perm, method = "BY")),
    bind_cols(tb_info, method = "GLM-BH", compute_stats( glm_fit, method = "BH")),
    bind_cols(tb_info, method = "GLM-BY", compute_stats( glm_fit, method = "BY"))
  )
  
}
