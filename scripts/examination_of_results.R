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
# Load all the result data
# ==============================================================================
result_files = list.files(results_dir, full.names = TRUE)

# Load all the model fits
fits = lapply(result_files, function(x) readRDS(x)$fits)
names(fits) = sub("(.+).rds", "\\1", basename(result_files))

sapply(fits, function(x) sapply(x, class) == "sam") |>
  apply(1, sum)

# ==============================================================================
# Create a lot of different plots
# ==============================================================================

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

cod_final_param_plot = local({
  file = grep("NEA_cod", result_files, value = TRUE)
  res = readRDS(file)
  name = sub("(.+).rds", "\\1", basename(file))
  model_names = c("final")
  new_names = c("Final")
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
  plot = cod_final_param_plot,
  width = 9,
  height = 5,
  file = file.path(image_dir, "cod_params.pdf")
)

# Plot Q/sigma/omega params in a pretty format, for the paper
plots = list()
for (i in seq_along(result_files)) {
  res = readRDS(result_files[i])
  name = sub("(.+).rds", "\\1", basename(result_files[i]))
  plots[[name]] = local({
    model_names = c("final", "maximal", "spline")
    new_names = c("Final", "Maximal", "Spline")
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
  file = file.path(image_dir, "params_selected.pdf")
)

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

# Plot SSB for BS_plaice in a pretty format
ff = fits[[grep("BS_plaice", names(fits))]]
# ff$final2 = local({
#   conf = ff$final$conf
#   conf$keyLogFpar[2, ] = c(0:5, -1)
#   conf$keyLogFpar[3, ] = c(6:11, -1)
#   par = stockassessment::defpar(ff$final$data, conf)
#   my_sam_fit(data, conf, par)
# })
model_names = c("final", "maximal", "spline", "final2")
new_names = c("Final", "Maximal", "Spline", "Modified final")
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
    breaks = seq(0, 2e5, by = 2e4),
    labels = paste0("$", c(0, paste(seq(2, 8, by = 2), "\\cdot 10^4"), "10^5", paste(seq(1.2, 2, by = .2), "\\cdot 10^5")), "$")
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

# Plot estimates of Fbar for each fish stock
fbar_df = lapply(
  seq_along(fits),
  FUN = function(i) {
    df = lapply(
      seq_along(fits[[i]]),
      function(j) {
        fit = fits[[i]][[j]]
        index = which(names(fit$sdrep$value) == "logfbar")
        if (length(index) == 0) {
          return(data.table(name = names(fits[[i]])[j], stock = names(fits)[i]))
        }
        logfbar = fit$sdrep$value[index]
        ci = logfbar + fit$sdrep$sd[index] %o% c(-1.96, 1.96)
        data.table(
          value = exp(logfbar),
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
fbar_df = rbindlist(fbar_df)
plots = list()
for (s in unique(fbar_df$stock)) {
  plots[[s]] = ggplot(fbar_df[stock == s]) +
    geom_ribbon(aes(x = y, ymin = lower, ymax = upper, fill = name), alpha = .3) +
    geom_line(aes(x = y, y = value, group = name, col = name)) +
    labs(col = "", x = "Year", fill = "", y = "Fbar", title = s) +
    theme_light() +
    theme(
      strip.text = element_text(colour = "black"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      text = element_text(size = 16),
      legend.position = "top"
    )
}
pdf(file.path(image_dir, "fbar.pdf"))
for (plot in plots) print(plot)
dev.off()

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
  df2 = rbindlist(df2, fill = TRUE, use.names = TRUE)
  rbind(df, df2)
})
cv_ssb = rbindlist(cv_ssb, use.names = TRUE)[fit_type != "base"]
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
cv_ssb[eval_type == "crossval" & stock == "BS_plaice_2022" & fit_type == "spline" & year == 2010, sort(exp(logssb))]
cv_ssb[eval_type == "crossval" & stock == "BS_plaice_2022" & fit_type == "maximal" & year == 2010, sort(exp(logssb))]
cv_ssb[eval_type == "crossval" & stock == "BS_plaice_2022" & fit_type == "final" & year == 2010, sort(exp(logssb))]

# ==============================================================================
# Examine simulation results
# ==============================================================================

sim_data = lapply(
  X = result_files,
  FUN = function(x) {
    data = readRDS(x)
    out = data$simulations
    if (is.null(out)) return(NULL)
    out$stock = sub("(.+).rds", "\\1", basename(x))
    out$region = strsplit(out$stock[1], "_")[[1]][1]
    out
  })
sim_data = rbindlist(sim_data)

tmp = unique(sim_data[, .(
  m = mean(pred), s = mean(sd)
),
by = c("model", "sim_nr", "simulation_model", "stock", "convergence")]
)
tmp[, let(bad_fit = (is.na(convergence) | convergence != 0 | is.nan(m) | is.nan(s)))]
tmp[, let(all_good = all(!bad_fit)), by = c("stock", "sim_nr", "simulation_model")]

tmp[, sum(bad_fit), by = "model"]
unique(tmp[, .(all_good, stock, sim_nr, simulation_model)])[, .(
  n_good = sum(all_good),
  n_bad = sum(!all_good),
  .N
)]

tmp[, .(
  n_crash = sum(is.na(convergence)),
  n_no_convergence = sum(convergence != 0, na.rm = TRUE),
  n_nan = sum(is.nan(m) | is.nan(s)),
  n_bad_total = sum(convergence != 0 | is.na(convergence) | is.nan(m) | is.nan(s), na.rm = TRUE)
  ),
by = "model"]

tmp[, .(n_bad = sum(bad_fit)), by = c("model", "simulation_model")] |>
  dcast(simulation_model ~ model)

tmp[, .(n_bd = sum(bad_fit)), by = c("model", "stock")] |>
  dcast(stock ~ model)

sim_data = merge(sim_data, tmp[, .(model, sim_nr, simulation_model, stock, all_good)], all.y = TRUE)
sim_data = sim_data[all_good == TRUE]

sim_eval = sim_data[, .(
  rmse = sqrt(mean((exp(pred) - exp(truth))^2)),
  rmse2 = sqrt(mean((pred - truth)^2))
), by = c("model", "stock", "simulation_model", "tag")]
sim_eval[, let(
  rmse = rmse / rmse[model == simulation_model],
  rmse2 = rmse2 / rmse2[model == simulation_model]
),
by = c("stock", "simulation_model", "tag")]
sim_eval = sim_eval[model != simulation_model]
sim_eval = melt(sim_eval, measure.vars = c("rmse", "rmse2"))
sim_eval[variable == "rmse", let(tag = ifelse(tag == "logssb", "SSB", "$\\bar F$"))]
sim_eval[variable == "rmse2", let(tag = ifelse(tag == "logssb", "log SSB", "$\\log \\bar F$"))]
sim_eval[, let(
  model = factor(
    model,
    levels = c("final", "maximal", "spline"),
    labels = c("Final", "Maximal", "Spline")
  ),
  simulation_model = factor(
    simulation_model,
    levels = c("final", "maximal", "spline"),
    labels = c("Final", "Maximal", "Spline")
  )
)]

plot = ggplot(sim_eval) +
  geom_boxplot(aes(x = simulation_model, y = value, col = model)) +
  facet_wrap(~ tag, nrow = 1) +
  geom_hline(yintercept = 1) +
  lims(y = c(0.75, 1.65)) +
  theme_light() +
  labs(x = "True model", y = "Standardised RMSE", col = "Fitted model") +
  theme(
    strip.text = element_text(colour = "black"),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
    text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, vjust = .5),
    legend.position = "top"
  )

plot_tikz(
  plot = plot,
  file = file.path(image_dir, "simulation_rmse.pdf"),
  width = 10,
  height = 3.8
)
  
probs = seq(.05, .95, by = .05)

coverage_data = list()
for (prob in probs) {
  z = qnorm(1 - (1 - prob) / 2)
  coverage_data[[length(coverage_data) + 1]] = sim_data[, .(
    coverage = mean((pred - sd * z) <= truth & truth <= (pred + sd * z)),
    prob = prob
  ),
  by = c("model", "simulation_model", "stock", "tag")]
  #by = c("model", "simulation_model")]
}
coverage_data = rbindlist(coverage_data)
coverage_data[, let(
  tag = factor(
    tag,
    levels = c("logssb", "logfbar"),
    labels = c("SSB", "$\\bar F$")
  ),
  model = factor(
    model,
    levels = c("final", "maximal", "spline"),
    labels = paste("Fitted model:\n", c("Final", "Maximal", "Spline"))
  ),
  simulation_model = factor(
    simulation_model,
    levels = c("final", "maximal", "spline"),
    labels = paste("True model:", c("Final", "Maximal", "Spline"))
  )
)]

plot = ggplot(coverage_data[tag == "SSB"]) +
  geom_boxplot(aes(x = prob, y = coverage, group = prob)) +
  geom_abline(slope = 1, intercept = 0) +
  coord_equal() +
  theme_light() +
  theme(
    strip.text = element_text(colour = "black"),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
    strip.text.y = element_text(angle = 0)
  ) +
  labs(x = "$p$", y = "Coverage probability") +
  scale_x_continuous(
    breaks = seq(0, 1, by = .25),
    labels = c("0", "0.25", "0.50", "0.75", "1"),
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = .25),
    labels = c("0", "0.25", "0.50", "0.75", "1"),
    limits = c(0, 1)
  ) +
  facet_grid(model ~ simulation_model)
  #facet_wrap(~model)

plot_tikz(
  plot = plot,
  file = file.path(image_dir, "simulation_coverage.pdf"),
  width = 6,
  #width = 11,
  height = 5.5
)


# ==============================================================================
# Create a nice latex-table that describes all the 17 fish stocks
# ==============================================================================

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
convergence_stats$n_all_convergence[1]
convergence_stats$n_all_convergence[1] / convergence_stats$n_total[1]

eval[, .N, by = "fit_type"]
dcast(eval[, .N, by = c("fit_type", "stock")], stock ~ fit_type, value.var = "N")

unique(eval[, .(
  year, fit_type, convergence, eval_type, stock, region
)])[, let(
  all_convergence = all(convergence)
), by = .(year, eval_type, stock, region)][, .(
  n_convergence = sum(convergence),
  n_all_convergence = sum(all_convergence),
  n_total = .N
), by = .(stock, fit_type)][, let(percentage = round(n_convergence / n_total, 3) * 100)][, .(stock, percentage, fit_type)] |>
  dcast(stock ~ fit_type, value.var = "percentage")

# log(a + 1) and log(a + .1), and especially log(a + .01), is worse for spline convergence than sqrt(a)

# Only keep results from the data subsets where all the competing models converged
eval[, let(all_convergence = all(convergence)), by = c("year", "stock", "eval_type")]
eval = eval[all_convergence == TRUE]

# Create a data.table containing all the score data
score = local({
  eval = copy(eval)[, let(fleet_type = ifelse(fleet == 1, "catch", "survey"))]
  score = eval[, .(
    rmse = sqrt(mean((exp(obs) - exp(pred))^2, na.rm = TRUE)),
    rmse2 = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
    conditional_rmse = sqrt(mean((exp(obs) - exp(conditional_pred))^2, na.rm = TRUE)),
    conditional_rmse2 = sqrt(mean((obs - conditional_pred)^2, na.rm = TRUE))
  ), by = c("fit_type", "fleet_type", "stock", "eval_type")]
  score[, let(
    rmse = rmse / rmse[fit_type == "final"],
    rmse2 = rmse2 / rmse2[fit_type == "final"],
    conditional_rmse = conditional_rmse / conditional_rmse[fit_type == "final"],
    conditional_rmse2 = conditional_rmse2 / conditional_rmse2[fit_type == "final"]
  ), by = c("fleet_type", "stock", "eval_type")]
  score = score[fit_type != "final"]
  score
})

# Remove the stocks where it makes no sense to do conditional prediction,
# due to a large amounts of missing observations
score[stock %in% c("F_saithe_2023", "F_haddock_2023"), let(conditional_rmse = NA, conditional_rmse2 = NA)]

# Print some stats for the RMSE scores. How often do the different models beat the
# final model? I.e., how often do they have a standardised RMSE < 1
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
message("\nlog-RMSE:")
score[, .(
  good_percentage = round(mean(rmse2 < 1, na.rm = TRUE), digits = 2)
),
by = c("fit_type", "fleet_type", "eval_type")] |>
  dcast(eval_type + fleet_type ~ fit_type, value.var = "good_percentage")


# Create a box plot with standardised RMSE values for all 17 fish stocks,
# and also highlight RMSE values from the three stocks that are given special
# focus in the paper.
plot_df = rbind(
  score[, .(fit_type, fleet_type, stock, rmse, tag = eval_type)],
  score[eval_type == "forward" & fleet_type == "catch", .(
    fit_type, fleet_type, stock, rmse = conditional_rmse, tag = "conditional"
  )],
  score[, .(fit_type, fleet_type = paste0("log-", fleet_type), stock, rmse = rmse2, tag = eval_type)],
  score[eval_type == "forward" & fleet_type == "catch", .(
    fit_type, fleet_type = paste0("log-", fleet_type), stock, rmse = conditional_rmse2, tag = "conditional"
  )]
)
plot_df[, let(
  fleet_type = factor(
    fleet_type,
    levels = c("catch", "survey", "log-catch", "log-survey")[c(1, 3, 2, 4)],
    labels = c("$\\hat C$", "$\\hat I$", "$\\log \\hat C$", "$\\log \\hat I$")[c(1, 3, 2, 4)]
  ),
  fit_type = factor(
    fit_type,
    levels = c("maximal", "spline"),
    labels = c("Maximal", "Spline")
  ),
  tag = factor(
    tag,
    levels = c("crossval", "forward", "conditional"),
    labels = c("Cross-validation", "Forward-validation", "Conditional\nforward-validation")
  )
)]
plot = ggplot(plot_df) +
  geom_boxplot(aes(x = as.numeric(fit_type), y = rmse, group = fit_type), outlier.size = 1) +
  scale_x_continuous(breaks = seq_along(levels(plot_df$fit_type)), labels = levels(plot_df$fit_type)) +
  facet_grid(fleet_type ~ tag, scales = "free_y") +
  scale_y_continuous(breaks = seq(0, 2, by = .1)) +
  #facet_grid(tag ~ fleet_type) +
  geom_hline(
    data = unique(plot_df[, .(tag, fleet_type, y = 1)]),
    aes(yintercept = y),
    linewidth = 1,
    linetype = "dashed",
    col = "darkgray"
  ) +
  labs(x = "Model", y = "Standardised RMSE") +
  theme_light() +
  theme(
    strip.text = element_text(colour = "black"),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
    strip.text.y = element_text(angle = 0)
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
  width = 7,
  height = 5,
  file = file.path(image_dir, "model_evaluation.pdf")
)


# # ==============================================================================
# # Plot the penalty term for the penalty parameter
# # ==============================================================================
# 
# f = function(x, k, d) exp(d * (k - x)) / (1 + exp(d * (k - x)))
# 
# f2 = function(x, k, d) d * (k - x) - log(1 + exp(d * (k - x)))
# 
# f3 = function(x, k, d) {
#   out = d * (k - x)
#   i1 = which(out >= 50)
#   i2 = which(out < 50)
#   if (any(i1)) out[i1] = 0
#   if (any(i2)) out[i2] = out[i2] - log(1 + exp(out[i2]))
#   out
# }
# 
# k = 5
# x = seq(k - 1, k + 2, length.out = 500)
# d = 2^(-1:8)
# df = expand.grid(x = x, d = d)
# df$y = f3(x = df$x, k = k, d = df$d)
# 
# plot = ggplot(df) +
#   geom_line(aes(x = x, y = y, group = d, col = factor(d)), alpha = .9) +
#   scale_color_viridis_d() +
#   theme_light() +
#   theme(
#     text = element_text(size = 15)
#   ) +
#   scale_x_continuous(
#     breaks = (k-5):(k + 5),
#     labels = c(paste0("$K - ", 5:1, "$"), "$K$", paste0("$K + ", 1:5, "$")),
#     expand = c(0, 0)
#   ) +
#   labs(x = "$\\rho$", y = "$e(\\rho; K, \\delta)$", col = "$\\delta$")
#   
# plot_tikz(
#   plot = plot,
#   width = 9,
#   height = 5,
#   file = file.path(image_dir, "penalty_penalisation.pdf")
# )
# 
