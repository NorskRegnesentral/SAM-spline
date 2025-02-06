library(ggplot2)
library(data.table)
library(patchwork)
library(here)

# ==============================================================================
# Paths and other stuff
# ==============================================================================
sam_dir = here::here()
image_dir = file.path(sam_dir, "raw_data", "figures")
if (!dir.exists(image_dir)) dir.create(image_dir)
results_dir = file.path(sam_dir, "raw_data", "results")

# ==============================================================================
# Load all the helper functions
# ==============================================================================
function_path = file.path(sam_dir, "scripts", "helper_functions.R")
source(function_path)

# ==============================================================================
# Load all the result data and create some plots
# ==============================================================================
result_files = list.files(results_dir, full.names = TRUE)

# Load all the model fits
fits = lapply(result_files, function(x) readRDS(x)$fits)
names(fits) = sub("(.+).rds", "\\1", basename(result_files))

# Plot all of the estimated Q/sigma/omega params
plots = list()
for (i in seq_along(result_files)) {
  res = readRDS(result_files[i])
  name = sub("(.+).rds", "\\1", basename(result_files[i]))
  plots[[name]] = plot_fits(res$fits) +
    labs(title = name)
}
pdf(file.path(image_dir, "params.pdf"))
for (plot in plots) print(plot)
dev.off()

# Plot estimates of SSB for each fish stock
ssb_df = lapply(
  seq_along(fits),
  FUN = function(i) {
    df = lapply(
      seq_along(fits[[i]]),
      function(j) {
        fit = fits[[i]][[j]]
        index = which(names(fit$sdrep$value) == "logssb")
        if (length(index) == 0) {
          return(data.table(name = names(fits[[i]])[j], stock = names(fits)[i]))
        }
        logssb = fit$sdrep$value[index]
        ci = logssb + fit$sdrep$sd[index] %o% c(-1.96, 1.96)
        data.table(
          value = exp(logssb),
          lower = exp(ci[, 1]),
          upper = exp(ci[, 2]),
          y = fit$data$years,
          name = names(fits[[i]])[j],
          stock = names(fits)[i]
        )
      })
    rbindlist(df, fill = TRUE)
  }
)
ssb_df = rbindlist(ssb_df)
ssb_df = ssb_df[name != "base"]
plots = list()
for (s in unique(ssb_df$stock)) {
  plots[[s]] = ggplot(ssb_df[stock == s]) +
    geom_ribbon(aes(x = y, ymin = lower, ymax = upper, fill = name), alpha = .3) +
    geom_line(aes(x = y, y = value, group = name, col = name)) +
    labs(col = "", x = "Year", fill = "", y = "SSB", title = s) +
    theme_light() +
    theme(
      strip.text = element_text(colour = "black"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      text = element_text(size = 16),
      legend.position = "top"
    )
}
pdf(file.path(image_dir, "ssb.pdf"))
for (plot in plots) print(plot)
dev.off()

# Plot Q/sigma/omega params in a pretty format, for the paper
plots = list()
for (i in seq_along(result_files)) {
  res = readRDS(result_files[i])
  name = sub("(.+).rds", "\\1", basename(result_files[i]))
  plots[[name]] = local({
    model_names = c("official", "maximal", "spline1", "spline2")
    new_names = c("Official", "Maximal", "Spline1", "Spline2")
    fits = res$fits[names(res$fits) %in% model_names]
    for (j in seq_along(model_names)) {
      if (model_names[j] %in% names(fits)) names(fits)[names(fits) == model_names[j]] = new_names[j]
    }
    areas = c("Baltic Sea", "North Sea", "North-East Arctic", "Faroe Plateau",
              "Celtic Sea", "Widely distributed")
    area_codes = c("BS", "NS", "NEA", "F", "CSWS", "WMS")
    area_code = strsplit(name, "_")[[1]][1]
    area = areas[area_codes == area_code]
    stock = strsplit(name, "_")[[1]][2]
    toupper_first = function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
    stock = toupper_first(stock)
    plot = plot_fits(fits, pretty = TRUE) +
      geom_line(aes(alpha = fit_type)) +
      labs(alpha = "", col = "", size = "", x = "$a$", y = "") +
      theme(legend.position = "top") +
      scale_y_continuous(breaks = -30:20) +
      scale_alpha_manual(values = c(0, rep(1, length(model_names))) * .5) +
      theme(text = element_text(size = 15))
    plot
  })
}
for (i in seq_along(result_files)) {
  res = readRDS(result_files[i])
  name = sub("(.+).rds", "\\1", basename(result_files[i]))
  plots[[paste0(name, "_", 2)]] = local({
    model_names = c("official")
    new_names = c("Official")
    fits = res$fits[names(res$fits) %in% model_names]
    for (j in seq_along(model_names)) {
      if (model_names[j] %in% names(fits)) names(fits)[names(fits) == model_names[j]] = new_names[j]
    }
    areas = c("Baltic Sea", "North Sea", "North-East Arctic", "Faroe Plateau",
              "Celtic Sea", "Widely distributed")
    area_codes = c("BS", "NS", "NEA", "F", "CSWS", "WMS")
    area_code = strsplit(name, "_")[[1]][1]
    area = areas[area_codes == area_code]
    stock = strsplit(name, "_")[[1]][2]
    toupper_first = function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
    stock = toupper_first(stock)
    plot = plot_fits(fits, pretty = TRUE) +
      geom_line(aes(alpha = fit_type)) +
      labs(alpha = "", col = "", size = "", x = "$a$", y = "") +
      guides(col = "none", size = "none", alpha = "none") +
      scale_y_continuous(breaks = -30:20, minor_breaks = NULL) +
      scale_alpha_manual(values = c(0, 1, 1, 1) * .5) +
      theme(text = element_text(size = 15))
    plot
  })
}
for (i in seq_along(result_files)) {
  res = readRDS(result_files[i])
  name = sub("(.+).rds", "\\1", basename(result_files[i]))
  plots[[paste0(name, "_", 3)]] = local({
    model_names = c("official", "maximal", "spline1", "spline2")
    new_names = c("Official", "Maximal", "Spline1", "Spline2")
    fits = res$fits[names(res$fits) %in% model_names]
    for (j in seq_along(model_names)) {
      if (model_names[j] %in% names(fits)) names(fits)[names(fits) == model_names[j]] = new_names[j]
    }
    areas = c("Baltic Sea", "North Sea", "North-East Arctic", "Faroe Plateau",
              "Celtic Sea", "Widely distributed")
    area_codes = c("BS", "NS", "NEA", "F", "CSWS", "WMS")
    area_code = strsplit(name, "_")[[1]][1]
    area = areas[area_codes == area_code]
    stock = strsplit(name, "_")[[1]][2]
    toupper_first = function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
    stock = toupper_first(stock)
    plot = plot_fits(fits, pretty = TRUE) +
      geom_line(aes(alpha = fit_type)) +
      labs(alpha = "", col = "", size = "", x = "$a$", y = "", fill = "") +
      theme(legend.position = "top") +
      scale_y_continuous(breaks = -30:20) +
      scale_alpha_manual(values = c(0, rep(1, length(model_names))) * .5) +
      geom_ribbon(aes(ymin = value - 1.96 * sd, ymax = value + 1.96 * sd, fill = fit_type), alpha = .2, col = NA) +
      theme(text = element_text(size = 15))
    plot
  })
}
plot_tikz(
  plot = plots[grep("(BS_plaice)|(NEA_cod)|(CSWS_plaice)", names(plots))],
  width = 9,
  height = 5,
  file = file.path(image_dir, "params_selected.pdf")
)

cod_official_param_plot = local({
  file = grep("NEA_cod", result_files, value = TRUE)
  res = readRDS(file)
  name = sub("(.+).rds", "\\1", basename(file))
  model_names = c("official")
  new_names = c("Official")
  fits = res$fits[names(res$fits) %in% model_names]
  for (j in seq_along(model_names)) {
    if (model_names[j] %in% names(fits)) names(fits)[names(fits) == model_names[j]] = new_names[j]
  }
  areas = c("Baltic Sea", "North Sea", "North-East Arctic", "Faroe Plateau",
            "Celtic Sea", "Widely distributed")
  area_codes = c("BS", "NS", "NEA", "F", "CSWS", "WMS")
  area_code = strsplit(name, "_")[[1]][1]
  area = areas[area_codes == area_code]
  stock = strsplit(name, "_")[[1]][2]
  toupper_first = function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  stock = toupper_first(stock)
  plot = plot_fits(fits, pretty = TRUE) +
    geom_line(aes(alpha = fit_type)) +
    labs(alpha = "", col = "", size = "", x = "$a$", y = "") +
    guides(col = "none", size = "none", alpha = "none") +
    scale_y_continuous(breaks = -30:20, minor_breaks = NULL) +
    scale_alpha_manual(values = c(0, 1, 1, 1) * .5) +
    theme(text = element_text(size = 15))
  plot
})
plot_tikz(
  plot = cod_official_param_plot,
  width = 9,
  height = 5,
  file = file.path(image_dir, "cod_params.pdf")
)


# Do the same once more, but with all three chosen fish stocks on the same page
myplots = plots[grep("(BS_plaice)|(NEA_cod)|(CSWS_plaice)", names(plots))]
myplot = patchwork::wrap_plots(
  myplots[c(3, 2, 1)],
  nrow = 3,
  guides = "collect",
  tag_level = "new"
) &
  theme(
    legend.position = "top",
    strip.text = element_text(size = rel(.7)),
    axis.text.x = element_text(size = rel(.8)),
    axis.text.y = element_text(size = rel(.8))
  ) &
  plot_annotation(tag_levels = "A", tag_suffix = ")")
plot_tikz(
  plot = myplot,
  width = 8,
  height = 10,
  file = file.path(image_dir, "params_selected2.pdf")
)

# Plot SSB for BS_plaice in a pretty format
ff = fits[[grep("BS_plaice", names(fits))]]
model_names = c("official", "maximal", "spline1", "spline2", "official2")
new_names = c("Official", "Maximal", "Spline1", "Spline2", "Official 2")
ssb_df = lapply(
  seq_along(ff),
  function(i) {
    index = which(names(ff[[i]]$sdrep$value) == "logssb")
    logssb = ff[[i]]$sdrep$value[index]
    ci = logssb + ff[[i]]$sdrep$sd[index] %o% c(-1.96, 1.96)
    data.table(
      value = exp(logssb),
      lower = exp(ci[, 1]),
      upper = exp(ci[, 2]),
      y = ff[[i]]$data$years,
      name = names(ff)[i]
    )
  }
)
ssb_df = rbindlist(ssb_df)
ssb_df = ssb_df[name %in% model_names]
ssb_df[, let(name = factor(name, levels = model_names, labels = new_names))]
plot = ggplot(ssb_df) +
  geom_ribbon(aes(x = y, ymin = lower, ymax = upper, fill = name), alpha = .3) +
  geom_line(aes(x = y, y = value, group = name, col = name)) +
  labs(col = "", x = "Year", fill = "", y = "SSB") +
  scale_y_continuous(
    breaks = seq(0, 1e5, by = 2e4),
    labels = paste0("$", c(0, paste(seq(2, 8, by = 2), "\\cdot 10^4"), "10^5"), "$")
  ) +
  theme_light() +
  theme(
    strip.text = element_text(colour = "black"),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
    text = element_text(size = 16),
    legend.position = "top"
  )
plot_tikz(
  plot = plot,
  width = 9,
  height = 5,
  file = file.path(image_dir, "bs-plaice_ssb.pdf")
)

# Plot all SSB estimates for all model fits from the cross/forward-evaluation
cv_ssb = lapply(result_files, function(x) {
  message(x)
  data = readRDS(x)
  df = lapply(
    X = data$crossval,
    FUN = function(y) {
      logssb = attr(y, "misc_info")$logssb
      years = unique(attr(y, "misc_info")$aux[, 1])
      data.table(
        logssb = logssb,
        year = years,
        cv_year = unique(y$year),
        fit_type = unique(y$fit_type),
        eval_type = "crossval",
        stock = sub("\\.rds", "", basename(x))
      )
    })
  df = rbindlist(df, fill = TRUE)
  df2 = lapply(
    X = data$forward,
    FUN = function(y) {
      logssb = attr(y, "misc_info")$logssb
      years = unique(attr(y, "misc_info")$aux[, 1])
      data.table(
        logssb = logssb,
        year = years,
        cv_year = unique(y$year),
        fit_type = unique(y$fit_type),
        eval_type = "forward",
        stock = sub("\\.rds", "", basename(x))
      )
    })
  df2 = rbindlist(df2, fill = TRUE)
  rbind(df, df2)
})
cv_ssb = rbindlist(cv_ssb)[fit_type != "base"]
plots = list()
for (s in unique(cv_ssb$stock)) {
  plots[[s]] = ggplot(cv_ssb[stock == s]) +
    geom_point(aes(x = year, y = exp(logssb))) +
    facet_grid(eval_type~fit_type)
}
pdf(file.path(image_dir, "cv-ssb.pdf"))
for (plot in plots) print(plot)
dev.off()

# Count the number of disagreeing SSB values for the different stocks for BS_plaice_2022
cv_ssb[eval_type == "crossval" & stock == "BS_plaice_2022" & fit_type == "spline1" & year == 2010, sort(exp(logssb))]
cv_ssb[eval_type == "crossval" & stock == "BS_plaice_2022" & fit_type == "spline2" & year == 2010, sort(exp(logssb))]
cv_ssb[eval_type == "crossval" & stock == "BS_plaice_2022" & fit_type == "maximal" & year == 2010, sort(exp(logssb))]
cv_ssb[eval_type == "crossval" & stock == "BS_plaice_2022" & fit_type == "official" & year == 2010, sort(exp(logssb))]


# ==============================================================================
# Examine RMSE values
# ==============================================================================

# Load all RMSE values from the cross/forward-validation
eval = lapply(result_files, function(x) {
  message(x)
  data = readRDS(x)
  eval = rbind(
    rbindlist(data$forward, fill = TRUE)[, let(eval_type = "forward")],
    rbindlist(data$crossval, fill = TRUE)[, let(eval_type = "crossval")],
    fill = TRUE
  )
  eval$stock = sub("(.+).rds", "\\1", basename(x))
  eval$region = strsplit(eval$stock[1], "_")[[1]][1]
  n_par = sapply(data$fits, function(x) length(x$opt$par))
  eval$n_par = n_par[eval$fit_type]
  eval
})
eval = rbindlist(eval, use.names = TRUE)

# Create a nice latex-table that describes all the 17 fish stocks
latex_table = local({
  toupper_first = function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  areas = c("Baltic Sea", "North Sea", "North-East Arctic", "Faroe Plateau",
            "Celtic Sea", "Widely distributed")
  area_codes = c("BS", "NS", "NEA", "F", "CSWS", "WMS")
  top = "\\begin{tabular}{llllll}"
  top = c(top, "\\toprule")
  top = c(top, paste(
    paste("Fish stock", "Area", "Years catch", "Years survey", "Ages", "Data source", sep = " & "),
    "\\\\"
  ))
  top = c(top, "\\midrule")
  top = paste(top, collapse = "\n")
  table = NULL
  for (i in seq_along(result_files)) {
    stock = strsplit(names(fits)[i], "_")[[1]][2]
    if (stock == "bw") stock = "Blue whiting"
    area_code = strsplit(names(fits)[i], "_")[[1]][1]
    area = areas[area_code == area_codes]
    data_source = readRDS(result_files[i])$webname
    age_range = paste(
      max(0, min(fits[[i]][[1]]$data$minAgePerFleet)),
      max(fits[[i]][[1]]$data$maxAgePerFleet),
      sep = " - "
    )
    aux = fits[[i]][[1]]$data$aux
    survey_years = paste(
      min(aux[aux[, 2] != 1, 1]),
      max(aux[aux[, 2] != 1, 1]),
      sep = " - "
    )
    catch_years = paste(
      min(aux[aux[, 2] == 1, 1]),
      max(aux[aux[, 2] == 1, 1]),
      sep = " - "
    )
    table = c(table, paste(
      toupper_first(stock), area, catch_years, survey_years, age_range, data_source, sep = " & "
    ))
  }
  table = paste(paste(table, collapse = " \\\\\n"), "\\\\")
  table = gsub("_", "\\\\_", table)
  bottom = c("\\bottomrule\n\\end{tabular}")
  full_table = paste(top, table, bottom, "\n", sep = "\n")
  full_table
})
cat(latex_table)

# Examine which models that converged for the different data subsets
convergence_stats = unique(eval[, .(
  year, fit_type, convergence, eval_type, stock, region
)])[, let(
  all_convergence = all(convergence)
), by = .(year, eval_type, stock, region)][, .(
  n_convergence = sum(convergence),
  n_all_convergence = sum(all_convergence),
  n_total = .N
), by = .(fit_type)][, let(percentage = round(n_convergence / n_total, 3) * 100)][]
convergence_stats
convergence_stats$n_all_convergence[1] / convergence_stats$n_total[1]

# Only keep results from the data subsets where all the competing models converged
eval[, let(all_convergence = all(convergence)), by = c("year", "stock", "eval_type")]
eval = eval[all_convergence == TRUE]

# Create a data.table containing all the score data
score = local({
  eval = copy(eval)[, let(fleet_type = ifelse(fleet == 1, "catch", "survey"))]
  score = eval[, .(
    rmse = sqrt(mean((exp(obs) - exp(pred))^2, na.rm = TRUE)),
    conditional_rmse = sqrt(mean((exp(obs) - conditional_pred)^2, na.rm = TRUE))
    #conditional_rmse2 = sqrt(mean((exp(obs) / sum(exp(obs), na.rm = TRUE) - conditional_pred / sum(conditional_pred, na.rm = TRUE))^2, na.rm = TRUE)),
  ), by = c("fit_type", "fleet_type", "stock", "eval_type")]
  score[, let(
    rmse = rmse / rmse[fit_type == "official"],
    conditional_rmse = conditional_rmse / conditional_rmse[fit_type == "official"]
  ), by = c("fleet_type", "stock", "eval_type")]
  score = score[fit_type != "official"]
  score
})

# Remove the stocks where it makes no sense to do conditional prediction,
# due to a large amounts of missing observations
score[stock %in% c("F_saithe_2023", "F_haddock_2023"), let(conditional_rmse = NA)]

# Print some stats for the RMSE scores. How often do the different models beat the
# official model? I.e., how often do they have a standardised RMSE < 1
message("\nMean forecast rmse:")
score[eval_type == "forward", .(
  good_percentage = round(mean(conditional_rmse < 1, na.rm = TRUE), digits = 2)
),
by = c("fit_type", "fleet_type")] |>
  dcast(fleet_type ~ fit_type, value.var = "good_percentage")
message("\nRMSE:")
score[, .(
  good_percentage = round(mean(rmse < 1, na.rm = TRUE), digits = 2)
),
by = c("fit_type", "fleet_type", "eval_type")] |>
  dcast(eval_type + fleet_type ~ fit_type, value.var = "good_percentage")


# Create a box plot with standardised RMSE values for all 17 fish stocks,
# and also highlight RMSE values from the three stocks that are given special
# focus in the paper.
plot_df = rbind(
  score[, .(fit_type, fleet_type, stock, rmse, tag = eval_type)],
  score[eval_type == "forward" & fleet_type == "catch", .(fit_type, fleet_type, stock, rmse = conditional_rmse, tag = "conditional")]
)
plot_df = plot_df[fit_type != "base"]
plot_df[, let(
  fleet_type = factor(fleet_type, levels = c("catch", "survey"), labels = c("Catch", "Survey")),
  fit_type = factor(
    fit_type,
    levels = c("maximal", "spline1", "spline2"),
    labels = c("Maximal", "Spline1", "Spline2")
  ),
  tag = factor(
    tag,
    levels = c("crossval", "forward", "conditional"),
    labels = c("Cross-validation", "Forward-validation", "Conditional forward-validation")
  )
)]
plot = ggplot(plot_df) +
  geom_boxplot(aes(x = as.numeric(fit_type), y = rmse, group = fit_type), outlier.size = 1) +
  scale_x_continuous(breaks = seq_along(levels(plot_df$fit_type)), labels = levels(plot_df$fit_type)) +
  facet_grid(fleet_type ~ tag) +
  geom_hline(
    data = unique(plot_df[, .(tag, fleet_type, y = 1)]),
    aes(yintercept = y),
    linetype = "dashed",
    col = "darkgray"
  ) +
  labs(x = "Model", y = "Standardised RMSE") +
  theme_light() +
  theme(
    strip.text = element_text(colour = "black"),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")
  ) +
  geom_point(
    data = plot_df[stock == "NEA_cod_2023"],
    aes(x = as.numeric(fit_type), y = rmse),
    col = "deepskyblue",
    size = 2
  ) +
  geom_point(
    data = plot_df[stock == "CSWS_plaice_2022"],
    aes(x = as.numeric(fit_type), y = rmse),
    col = "indianred1",
    size = 2
  ) +
  geom_point(
    data = plot_df[stock == "BS_plaice_2022"],
    aes(x = as.numeric(fit_type), y = rmse),
    col = "gold",
    size = 2
  )

plot_tikz(
  plot = plot,
  width = 9,
  height = 5,
  file = file.path(image_dir, "model_evaluation.pdf")
)

# ==============================================================================
# Plot the penalty prior
# ==============================================================================

f = function(x, k, d) exp(d * (k - x)) / (1 + exp(d * (k - x)))

k = 5
x = seq(k - 4, k + 4, length.out = 500)
d = c(
  #10^seq(-1, 3, by = 0.5)
  #.5, 1, 5, 10, 50, 100, 500, 1000
  2^(-1:8)
  #seq(.1, 1, by = .1),
  #10^seq(-1, 3.4, length.out = 100)
  #seq(2, 10, by = 1),
  #seq(20, 100, by = 10),
  #seq(200, 2000, by = 100)
)
df = expand.grid(x = x, d = d)
df$y = f(x = df$x, k = k, d = df$d)

plot = ggplot(df) +
  geom_line(aes(x = x, y = y, group = d, col = factor(d)), alpha = .9) +
  scale_color_viridis_d() +
  #scale_color_viridis_d(
  #  breaks = d,
  #  labels = paste0("$10^{", log10(d), "}$")
  #) +
  #scale_color_viridis_c(
  #  breaks = -1:3,
  #  labels = paste0("$", c(0.1, 1, 10, 100, 1000), "$")
  #) +
  theme_light() +
  theme(
    #axis.title.y = element_text(angle = 0, vjust = .5),
    text = element_text(size = 15)
  ) +
  scale_x_continuous(
    breaks = (k-5):(k + 5),
    labels = c(paste0("$K - ", 5:1, "$"), "$K$", paste0("$K + ", 1:5, "$")),
    expand = c(0, 0)
  ) +
  labs(x = "$\\rho$", y = "$\\pi(\\rho; K, \\delta)$", col = "$\\delta$")
  
plot_tikz(
  plot = plot,
  width = 9,
  height = 5,
  file = file.path(image_dir, "penalty_prior.pdf")
)

