# functions for plotting

plot_dag = function(dag, title_str) {
  expand_plot = function(expand_x = expansion(c(.1, .1)),
                         expand_y = expansion(c(.1, .1))) {
    list(
      ggplot2::scale_x_continuous(expand = expand_x),
      ggplot2::scale_y_continuous(expand = expand_y)
    )
  }
  ggdag(dag, layout = "circle") + 
    geom_dag_point() +
    geom_dag_text() +
    theme_dag_blank() +
    expand_plot(expand_x = expansion(c(.4, .4)), 
                expand_y = expansion(c(.4, .4))) +
    ggtitle(title_str) +
    theme(plot.title = element_text(hjust = 0.5))
}

plot_cor = function(tb_summary, .paired = TRUE, .sigma_b, .sigma_u) {
  
  tb_summary %<>% 
    filter(sigma_b == .sigma_b & sigma_u == .sigma_u) %>%
    filter(paired == .paired)
  
  if(.paired) {
    title_str = bquote("Power ( paired /"~sigma[B]==.(.sigma_b)~"/"~sigma[U]==.(.sigma_u)~")")
  } else {
    title_str = bquote("Power ( unpaired /"~sigma[B]==.(.sigma_b)~"/"~sigma[U]==.(.sigma_u)~")")
  }
  p_power = ggplot(tb_summary %>% filter(paired == .paired) %>% mutate(rho_u = -rho_u), 
         aes(rho_b, power, color = method)) +
    geom_line() + 
    geom_point() +
    facet_wrap(~rho_u, labeller = label_bquote(rho[U] == -.(rho_u)), 
               ncol = length(unique(tb_summary$rho_u))) + 
    ggtitle(title_str) +
    scale_color_few() +
    scale_x_continuous(trans = "reverse", breaks = unique(tb_summary$rho_u)) +
    theme(axis.title.y = element_blank(), panel.spacing.x = unit(1, "lines")) +
    xlab(bquote(rho[B])) +
    ylim(0, 1)
  
  if(.paired) {
    title_str = bquote("FDR ( paired /"~sigma[B]==.(.sigma_b)~"/"~sigma[U]==.(.sigma_u)~")")
  } else {
    title_str = bquote("FDR ( unpaired /"~sigma[B]==.(.sigma_b)~"/"~sigma[U]==.(.sigma_u)~")")
  }
  p_fdr = ggplot(tb_summary %>% filter(paired == .paired) %>% mutate(rho_u = -rho_u),
         aes(rho_b, FDR, color = method)) +
    geom_hline(yintercept = unique(tb_summary$fdr), linetype = "dashed", alpha = 0.5) +
    geom_line() + 
    geom_point() +
    facet_wrap(~rho_u, labeller = label_bquote(rho[U] == -.(rho_u)), 
               ncol = length(unique(tb_summary$rho_u))) + 
    ggtitle(title_str) +
    scale_color_few() +
    scale_x_continuous(trans = "reverse", breaks = unique(tb_summary$rho_u)) +
    theme(axis.title.y = element_blank(), panel.spacing.x = unit(1, "lines")) +
    xlab(bquote(rho[B])) +
    ylim(0, 0.25)
  
  plot_grid(p_power, p_fdr, nrow = 2)
  
}

plot_sig = function(tb_summary, .paired = TRUE, .rho_b, .rho_u) {
  
  tb_summary %<>% 
    filter(rho_b == .rho_b & rho_u == .rho_u) %>%
    filter(paired == .paired)
  
  if(.paired) {
    title_str = bquote("Power ( paired /"~rho[B]==.(.rho_b)~"/"~rho[U]==.(.rho_u)~")")
  } else {
    title_str = bquote("Power ( unpaired /"~rho[B]==.(.rho_b)~"/"~rho[U]==.(.rho_u)~")")
  }
  p_power = ggplot(tb_summary, aes(sigma_b, power, color = method)) +
    geom_line() + 
    geom_point() +
    facet_wrap(~sigma_u, labeller = label_bquote(sigma[U] == .(sigma_u)), 
               ncol = length(unique(tb_summary$sigma_u))) + 
    ggtitle(title_str) +
    scale_color_few() +
    scale_x_continuous(breaks = unique(tb_summary$sigma_b)) +
    theme(axis.title.y = element_blank(), panel.spacing.x = unit(1, "lines")) +
    xlab(bquote(sigma[B])) +
    ylim(0, 1)
  
  if(.paired) {
    title_str = bquote("FDR ( paired /"~rho[B]==.(.rho_b)~"/"~rho[U]==.(.rho_u)~")")
  } else {
    title_str = bquote("FDR ( unpaired /"~rho[B]==.(.rho_b)~"/"~rho[U]==.(.rho_u)~")")
  }
  p_fdr = ggplot(tb_summary, aes(sigma_b, FDR, color = method)) +
    geom_hline(yintercept = unique(tb_summary$fdr), linetype = "dashed", alpha = 0.5) +
    geom_line() + 
    geom_point() +
    facet_wrap(~sigma_u, labeller = label_bquote(sigma[U] == .(sigma_u)), 
               ncol = length(unique(tb_summary$sigma_b))) + 
    ggtitle(title_str) +
    scale_color_few() +
    scale_x_continuous(breaks = unique(tb_summary$sigma_b)) +
    theme(axis.title.y = element_blank(), panel.spacing.x = unit(1, "lines")) +
    xlab(bquote(sigma[B])) +
    ylim(0, 0.25)
  
  plot_grid(p_power, p_fdr, nrow = 2)
  
}

plot_bias = function(tb_summary, num_cells = 1000) {
  
  tb_summary_select = tb_summary %>% filter(n_markers == 25)
  ggplot(tb_summary_select %>% 
           mutate(paired = if_else(paired, "paired", "unpaired")), 
         aes(bias, value, shape = type, color = method)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = 0.8, linetype = "dashed", alpha = 0.5) +
    geom_line() + 
    geom_point() +
    facet_wrap(~paired, labeller = label_value, ncol = 2) + 
    ggtitle(paste0("Paired and Unpaired with n_cells = ", num_cells)) +
    ylab("probability") +
    scale_color_manual(values = c("#5DA5DA", "#FAA43A")) +
    ylim(0, 1)
  
}

plot_experiment = function(tb_summary, num_cells = 1000) {
  # paired
  tb_summary_select = tb_summary %>% filter(
    paired == TRUE &
      effect_size == 1.5 &
      rho_b == 0.5 & 
      rho_u == 0.5 & 
      n_true == 5 &
      n_cells == num_cells & 
      bias %in% c(0, 1)
  )
  p_paired = ggplot(tb_summary_select, aes(n_markers, value, shape = type, color = method)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = 0.8, linetype = "dashed", alpha = 0.5) +
    geom_line() + 
    geom_point() +
    facet_wrap(~bias, labeller = label_both, 
               ncol = nlevels(factor(tb_summary_select$bias))) + 
    ggtitle(paste0("Paired with n_cells = ", num_cells)) +
    ylab("probability") +
    scale_color_manual(values = c("#5DA5DA", "#FAA43A")) +
    ylim(0, 1)
  
  # unpaired
  tb_summary_select = tb_summary %>% filter(
    paired == FALSE &
      effect_size == 30 &
      rho_b == 0.5 & 
      rho_u == 0.5 & 
      n_true == 5 &
      n_cells == num_cells & 
      bias %in% c(0, 1)
  )
  p_unpaired = ggplot(tb_summary_select, aes(n_markers, value, shape = type, color = method)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = 0.8, linetype = "dashed", alpha = 0.5) +
    geom_line() + 
    geom_point() +
    facet_wrap(~bias, labeller = label_both, 
               ncol = nlevels(factor(tb_summary_select$bias))) + 
    ggtitle(paste0("Unpaired with n_cells = ", num_cells)) +
    ylab("probability") +
    scale_color_manual(values = c("#5DA5DA", "#FAA43A")) +
    ylim(0, 1)
  
  # arrange the three plots in a single row
  prow = plot_grid(
    p_paired + theme(legend.position="none"),
    p_unpaired + theme(legend.position="none"),
    align = 'vh',
    hjust = -1,
    nrow = 1
  )
  
  # extract the legend from one of the plots
  legend = get_legend(
    # create some space to the left of the legend
    p_paired + theme(legend.box.margin = margin(0, 0, 0, 12))
  )
  
  # add the legend to the row we made earlier. Give it one-third of 
  # the width of one plot (via rel_widths).
  plot_grid(prow, legend, rel_widths = c(3, 0.6))
}
