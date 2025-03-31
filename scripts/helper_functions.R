library(stockassessment)
library(ggplot2)
library(data.table)
library(scales)
library(tikzDevice)

# This function is essentially a tryCatch wrapper around stockassessment::sam.fit,
# which returns NULL if sam.fit crashes
my_sam_fit = function(data,
                      conf,
                      parameters,
                      map = list(),
                      verbose = FALSE,
                      ...) {
  fit = tryCatch({
    stockassessment::sam.fit(
      data = data,
      conf = conf,
      parameters = parameters,
      map = map,
      silent = !verbose,
      trace = ifelse(verbose, 1, 0),
      ...
    )},
    error = function(e) {
      message("stockassessment::sam.fit crashed with error message:\n", e$message)
      NULL
    })
  fit
}

# This function performs cross-validation and forward-validation for a list of
# different sam objects. The cross/forward-validation is performed in parallel
do_evaluation = function(fits,
                         eval_years,
                         min_obs_years = 5,
                         n_cores = 8,
                         forward = FALSE,
                         eval_data = NULL) {
  start_time = Sys.time()
  combinations = expand.grid(year = eval_years, fit_nr = seq_along(fits))
  res = parallel::mclapply(
    X = seq_len(nrow(combinations)),
    mc.cores = n_cores,
    mc.preschedule = FALSE,
    FUN = function(i) {
      fit_nr = combinations$fit_nr[i]
      year = combinations$year[i]
      if (!is.null(eval_data)) {
        fit_type = names(fits)[fit_nr]
        all_eval_years = sapply(eval_data, function(x) unique(x$year))
        all_eval_fit_types = sapply(eval_data, function(x) unique(x$fit_type))
        if (any(fit_type == all_eval_fit_types & year == all_eval_years)) {
          index = which(fit_type == all_eval_fit_types & year == all_eval_years)
          return(eval_data[[index]])
        }
      }
      if (forward) {
        fit = forward_validate(
          year = year,
          data = fits[[fit_nr]]$data,
          conf = fits[[fit_nr]]$conf,
          min_obs_years = min_obs_years
        )

        if (!is.null(fit)) {
          all_data = fits[[fit_nr]]$data
          catch_mean_weight = all_data$catchMeanWeight[, , 1]
          catch_mean_weight = catch_mean_weight[row.names(catch_mean_weight) == as.character(year), ]
          obs = as.data.frame(cbind(all_data$aux, obs = exp(all_data$logobs)))
          obs = obs[obs$year == year & obs$fleet == 1, ]
          obs = obs[order(obs$age), ]
          catch_mean_weight = catch_mean_weight[as.integer(names(catch_mean_weight)) %in% obs$age]
          stopifnot(nrow(obs) == length(catch_mean_weight))
          obs$est_weight = obs$obs * catch_mean_weight
          tot_weight = sum(obs$est_weight)
          conditional_pred = tryCatch(
          {
            set.seed(1)
            # This is a forecast for the values at the end of the year
            fc = forecast(
              fit = fit,
              savesim = TRUE,
              catchval.exact = tot_weight
            )
            pred_mean = apply(fc[[1]]$catchatage, 1, mean)
            pred_mean
          },
          error = function(e) NA
          )
          all_ages = sort(unique(fit$data$aux[, "age"]))
          conditional_pred = data.table(
            conditional_pred = conditional_pred[all_ages %in% obs$age],
            age = obs$age,
            fleet = obs$fleet
          )
        }
      } else {
        fit = cross_validate(
          years = year,
          data = fits[[fit_nr]]$data,
          conf = fits[[fit_nr]]$conf
        )
      }
      if (is.null(fit)) {
        message("fit_nr: ", fit_nr, "/", length(fits),
                ", year: ", which(year == eval_years), "/", length(eval_years),
                ", FAILED")
        return(data.table(fit_type = names(fits)[fit_nr], year = year, convergence = FALSE, convergence_sam = FALSE))
      }
      misc_info = list(
        logN = fit$pl$logN,
        logssb = fit$sdrep$value[names(fit$sdrep$value) == "logssb"],
        aux = fit$data$aux,
        par = tryCatch(get_age_dependent_params(fit), error = function(e) NULL)
      )
      if (is.null(misc_info$par)) warning("get_age_dependent_params failed for fit_nr ", fit_nr, ", year ", year)
      year_index = which(fit$prediction$year == year)
      res = data.table(
        fit_type = names(fits)[fit_nr],
        convergence = (fit$opt$nlminb_convergence == 0),
        convergence_sam = (fit$opt$convergence == 0),
        year = year,
        obs = fit$prediction$obs,
        pred = fit$prediction$pred,
        fleet = fit$prediction$fleet[year_index],
        age = fit$prediction$age[year_index]
      )
      if (forward) {
        res = merge(res, conditional_pred, all.x = TRUE)
        res$conditional_pred
      }
      attr(res, "misc_info") = misc_info
      time_passed = difftime(Sys.time(), start_time)
      message("fit: ", fit_nr, "/", length(fits),
              ", year: ", which(year == eval_years), "/", length(eval_years),
              ", time passed: ", round(time_passed, 2), " ", attr(time_passed, "units"))
      res
    })
  res
}

# Compute the spline matrices X and S, and place them inside conf.
# The matrices are computed using the recipe in recipe_func().
# The X matrix should have one column for each spline parameter and one row for each observation.
# The S matrices should have one row and one column for each spline parameter. Furthermore,
# the actual conf$logFparSplineS/conf$varObsSplineS objects should be lists, an each element of the
# list should be a sparse matrix. S_rank should be a vector of equal length to the number of S matrices,
# describing the rank of each S matrix.
add_spline_info_to_conf = function(conf,
                                   recipe_func,
                                   data,
                                   logFpar = TRUE,
                                   varObs = TRUE,
                                   ...) {
  if (logFpar) {
    # Extract the SAM key, which tells us which age groups that are available for the
    # catch-at-age estimates and for each of the survey fleet types. Update it to reflect
    # the fact that each logFpar parameter should differ, since we are using a spline to
    # model logFpar
    key = local({
      key = t(conf[["keyLogFpar"]])
      key[key != -1] = seq_len(sum(key != -1)) - 1
      t(key)
    })
    tmp = recipe_func(conf, key, data, ...)
    conf[["logFparSplineX"]] = tmp[["X"]]
    conf[["logFparSplineS"]] = tmp[["S"]]
    conf[["logFparSplineSrank"]] = tmp[["S_rank"]]
    conf[["keyLogFpar"]] = key
    # If S_rank is zero, then there is nothing to penalise. Trying to penalise when
    # there is nothing to penalise will lead to overflow errors in c++
    if (sum(conf[["logFparSplineSrank"]]) == 0) conf[["useLogFparSplinePenalty"]] = 0
  }
  if (varObs) {
    # Extract the SAM key, which tells us which age groups that are available for the
    # catch-at-age estimates and for each of the survey fleet types. Update it to reflect
    # the fact that each logFpar parameter should differ, since we are using a spline to
    # model logFpar
    key = local({
      key = t(conf[["keyVarObs"]])
      key[key != -1] = seq_len(sum(key != -1)) - 1
      t(key)
    })
    tmp = recipe_func(conf, key, data, ...)
    conf[["varObsSplineX"]] = tmp[["X"]]
    conf[["varObsSplineS"]] = tmp[["S"]]
    conf[["varObsSplineSrank"]] = tmp[["S_rank"]]
    conf[["keyVarObs"]] = key
    # If S_rank is zero, then there is nothing to penalise. Trying to penalise when
    # there is nothing to penalise will lead to overflow errors in c++
    if (sum(conf[["varObsSplineSrank"]]) == 0) conf[["usevarObsSplinePenalty"]] = 0
  }
  # Add the function for updating the conf as an attribute to the conf. This way, if we
  # want to do e.g. forward-validation, and we remove some data points, it is easy to
  # automatically update the spline matrices so that they reflect the changes in the data
  attr(conf, "update_spline_info_in_conf_func") = function(x, data) {
    add_spline_info_to_conf(
      conf = x,
      data = data,
      recipe_func = recipe_func,
      logFpar = logFpar,
      varObs = varObs,
      ...
    )
  }
  conf
}

# This function takes in a list of sam objects, and plots parameter estimates for
# logFpar and logSdLogObs for all the sam objects
plot_fits = function(fits, pretty = FALSE) {
  if (is.null(names(fits))) names(fits) = paste0("fit_", seq_along(fits))
  df = lapply(
    seq_along(fits),
    function(i) {
      if (!is(fits[[i]], "sam")) return(NULL)
      out = get_age_dependent_params(fits[[i]])
      out$fit_type = names(fits)[i]
      out
    })
  df = do.call(rbind, df)
  all_fit_types = sort(unique(df$fit_type))
  n_fits = length(all_fit_types)
  if (n_fits == 1) {
    cols = "gray"
    sizes = 2.5
  } else {
    has_official_model = any(tolower(all_fit_types) == "official")
    if (has_official_model) {
      official_index = which(tolower(all_fit_types) == "official")
      all_fit_types = c(all_fit_types[official_index], all_fit_types[-official_index])
      cols = c("gray", scales::hue_pal()(n_fits - 1))
      sizes = c(2.5, rep(1, n_fits - 1))
    } else {
      cols = scales::hue_pal()(n_fits)
      sizes = rep(1, n_fits)
    }
  }
  df$fit_type = factor(df$fit_type, levels = all_fit_types)
  if (length(unique(df$fleet)) > 2) {
    fleets = sort(unique(df$fleet))
    df$fleet = factor(
      df$fleet,
      levels = fleets,
      labels = c("Catch", paste("Survey", tail(fleets, -1) - 1))
    )
  } else {
    fleets = sort(unique(df$fleet))
    df$fleet = factor(
      df$fleet,
      levels = fleets,
      labels = c("Catch", "Survey")
    )
  }
  if (pretty) {
    df$tag = factor(df$tag, levels = c("logFpar", "varObs"), labels = c("$\\log Q_{a, j}$", "$\\log \\sigma_a/\\log \\omega_{a, j}$"))
  }
  plot = ggplot(df, aes(x = a, y = value, col = fit_type, group = fit_type)) +
    geom_point(aes(size = fit_type)) +
    scale_size_manual(values = sizes) +
    facet_grid(tag ~ fleet, scales = "free") +
    labs(y = "Value", col = "Model", x = ifelse(pretty, "$a$", "Age"), size = "Model") +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    theme_light() +
    theme(
      strip.text = element_text(colour = "black"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")
    )
  if (diff(range(df$a)) > 8) {
    plot = plot + scale_x_continuous(breaks = seq(-100, 100, by = 2), minor_breaks = NULL)
  } else {
    plot = plot + scale_x_continuous(breaks = seq(-100, 100, by = 1), minor_breaks = NULL)
  }
  plot
}

# This function is used as a `recipe_func` in the add_spline_info_to_conf() function.
# Create a "cs" spline with a suitable penalty matrix, using mgcv
get_spline1_conf_parts = function(conf, key, data, k = Inf, ...) {
  # Preallocate lists for the spline matrices X and S, which we will compute for
  # each row of the key matrix
  X = S = vector("list", nrow(key))
  # We need to keep track of the rank of S. At the moment it is zero, and then
  # we will update S_rank as S is made
  S_rank = 0
  for (i in seq_len(nrow(key))) {
    # Find out which age groups that should be given parameters for this row of the key matrix
    age_group_index = which(key[i, ] != -1)
    # If there are no age groups, then skip to the next row of the key matrix
    if (!any(age_group_index)) next
    ages = age_group_index + conf$minAge - 1
    log_ages = log(ages + 1)
    n_ages = length(ages)
    if (n_ages < 3) {
      # If the catch-at-age estimates or survey indices of interest only contains two unique
      # age groups, or less, then it makes no sense to fit a spline to the data.
      # In that case, we give unique parameter values to each age group, and we do not
      # perform any penalisation of the parameter values.
      X[[i]] = diag(n_ages)
      S[[i]] = matrix(0, n_ages, n_ages)
    } else {
      # Create the spline matrices X and S, using mgcv
      spline_model = mgcv::gam(
        formula = y ~ s(log_age, bs = "cs", k = min(n_ages, k)),
        data = data.frame(y = 1, log_age = log_ages),
        fit = FALSE
      )
      X[[i]] = spline_model$X
      S[[i]] = spline_model$smooth[[1]]$S[[1]]
      # Reweight the elements of S that belong to the youngest age groups
      weights = pmin(1, exp(-4 + seq_len(ncol(S[[i]]))))
      S[[i]] = diag(weights) %*% S[[i]] %*% diag(weights)
      # Add zero columns/rows to S[[i]] to ensure that it has the same number of
      # rows and columns as the number of parameters. In this case, we know that
      # S contains a row/column for all spline parameters, but not for the intercept term,
      # which is the first parameter. We therefore add one column/row of zeros to S
      S[[i]] = rbind(0, cbind(0, S[[i]]))
      # Update S_rank to account for the added S matrix
      S_rank = S_rank + spline_model$smooth[[1]]$rank
    }
  }
  res = list(
    X = as.matrix(Matrix::bdiag(X[lengths(X) > 0])),
    S = list(as(Matrix::bdiag(S[lengths(S) > 0]), "TsparseMatrix")),
    S_rank = S_rank
  )
  res
}

# This function is used as a `recipe_func` in the add_spline_info_to_conf() function.
# Create a "bs" spline with a suitable penalty matrix, using mgcv
get_spline2_conf_parts = function(conf, key, data, k = Inf, ...) {
  # Preallocate lists for the spline matrices X and S, which we will compute for
  # each row of the key matrix
  X = S = vector("list", nrow(key))
  # We need to keep track of the rank of S. At the moment it is zero, and then
  # we will update S_rank as S is made
  S_rank = 0
  for (i in seq_len(nrow(key))) {
    # Find out which age groups that should be given parameters for this row of the key matrix
    age_group_index = which(key[i, ] != -1)
    # If there are no age groups, then skip to the next row of the key matrix
    if (!any(age_group_index)) next
    ages = age_group_index + conf$minAge - 1
    log_ages = log(ages + 1)
    n_ages = length(ages)
    if (n_ages < 3) {
      # If the catch-at-age estimates or survey indices of interest only contains two unique
      # age groups, or less, then it makes no sense to fit a spline to the data.
      # In that case, we give unique parameter values to each age group, and we do not
      # perform any penalisation of the parameter values.
      X[[i]] = diag(n_ages)
      S[[i]] = matrix(0, n_ages, n_ages)
    } else {
      # We cannot use a cubic spline if there are four or less distinct age groups. So if n_ages <= 4,
      # we use a quadratic spline instead. The first element of m describes the order of the spline,
      # i.e. 2 for quadratic and 3 for cubic. The second element of m describes which derivative we penalise,
      # i.e. 2 means that we penalise the second derivative of the spline.
      m = if (n_ages <= 4) c(2, 2) else c(3, 2)
      # Create the spline matrices X and S, using mgcv
      spline_model = mgcv::gam(
        formula = y ~ s(log_age, bs = "bs", k = min(n_ages, k), m = m),
        data = data.frame(y = 1, log_age = log_ages),
        fit = FALSE
      )
      X[[i]] = spline_model$X
      S[[i]] = spline_model$smooth[[1]]$S[[1]]
      # Reweight the elements of S that belong to the youngest age groups
      weights = pmin(1, exp(-4 + seq_len(ncol(S[[i]]))))
      S[[i]] = diag(weights) %*% S[[i]] %*% diag(weights)
      # Add zero columns/rows to S[[i]] to ensure that it has the same number of
      # rows and columns as the number of parameters. In this case, we know that
      # S contains a row/column for all spline parameters, but not for the intercept term,
      # which is the first parameter. We therefore add one column/row of zeros to S
      S[[i]] = rbind(0, cbind(0, S[[i]]))
      # Update S_rank to account for the added S matrix
      S_rank = S_rank + spline_model$smooth[[1]]$rank
    }
  }
  res = list(
    X = as.matrix(Matrix::bdiag(X[lengths(X) > 0])),
    S = list(as(Matrix::bdiag(S[lengths(S) > 0]), "TsparseMatrix")),
    S_rank = S_rank
  )
  res
}

# The actual function for doing cross-validation, which is called from inside
# the do_evaluation function
cross_validate = function(years, data, conf) {

  # Find out which parts of the data to set to NA
  missing_index = which(data$aux[, "year"] %in% years)
  stopifnot(any(missing_index))
  removed_obs = data$logobs[missing_index]
  data$logobs[missing_index] = NA

  # Find out which parts of the data are equal to NA
  # (some observations might have been NA before we started doing all this)
  na_index = which(is.na(data$logobs))

  # Perform inference for the remaining data
  fit = my_sam_fit(
    data = data,
    conf = conf,
    parameters = stockassessment::defpar(data, conf),
    newtonsteps = 0
  )

  if (is.null(fit)) return(NULL)

  # Add a matrix with predicted values and standard deviations, and true observed values,
  # for the removed observations
  fit$prediction = as.data.frame(cbind(
    data$aux[missing_index, , drop = FALSE],
    obs = removed_obs,
    pred = fit$pl$missing[na_index %in% missing_index]
  ))

  fit
}

# The actual function for doing forward-validation, which is called from inside
# the do_evaluation function
forward_validate = function(year,
                            data,
                            conf,
                            min_obs_years = 5,
                            n_pred_years = 1) {

  fleets_to_keep = fleets_to_drop = NULL
  for (f in sort(unique(data$aux[, "fleet"]))) {
    fleet_years = unique(data$aux[data$aux[, "fleet"] == f, "year"])
    fleet_years = fleet_years[fleet_years < year]
    if (length(fleet_years) >= min_obs_years) {
      fleets_to_keep = c(fleets_to_keep, f)
    } else {
      fleets_to_drop = c(fleets_to_drop, f)
    }
  }

  update_spline_info_in_conf_func = attr(conf, "update_spline_info_in_conf_func")

  if (!is.null(fleets_to_drop)) {
    data = stockassessment::reduce(data, fleet = fleets_to_drop, conf = conf)
    conf = attr(data, "conf")
    if (!is.null(update_spline_info_in_conf_func)) conf = update_spline_info_in_conf_func(conf, data)
  }

  years_to_drop = data$years[data$years >= year + n_pred_years]
  if (length(years_to_drop) == 0) years_to_drop = NULL

  if (!is.null(years_to_drop)) {
    data = stockassessment::reduce(data, year = years_to_drop, conf = conf)
    conf = attr(data, "conf")
    if (!is.null(update_spline_info_in_conf_func)) conf = update_spline_info_in_conf_func(conf, data)
  }

  # Find out which parts of the data to set to NA
  missing_index = which(data$aux[, "year"] >= year)
  stopifnot(any(missing_index))
  removed_obs = data$logobs[missing_index]
  data$logobs[missing_index] = NA

  # Find out which parts of the data are equal to NA
  # (some observations might have been NA before we started doing all this)
  na_index = which(is.na(data$logobs))

  # Perform inference for the remaining data
  fit = my_sam_fit(
    data = data,
    conf = conf,
    parameters = stockassessment::defpar(data, conf),
    newtonsteps = 0
  )

  if (is.null(fit)) return(NULL)

  # Add a matrix with predicted values and standard deviations, and true observed values,
  # for the removed observations
  fit$prediction = as.data.frame(cbind(
    data$aux[missing_index, , drop = FALSE],
    obs = removed_obs,
    pred = fit$pl$missing[na_index %in% missing_index]
  ))

  fit
}

# Helper function for extracting a data.frame
# with logFpar/logSdLogObs parameter estimates from a sam model fit
get_age_dependent_params = function(fit) {
  info = list(
    logFpar = c("logFpar", "keyLogFpar", "useLogFparSpline"),
    varObs = c("logSdLogObs", "keyVarObs", "useVarObsSpline"))
  out = list()
  for (i in seq_along(info)) {
    out[[i]] = list()
    par_index = which(names(fit$sdrep$value) == paste0("paraset.", info[[i]][[1]]))
    stopifnot(length(par_index) > 0)
    par = fit$sdrep$value[par_index]
    par_sd = fit$sdrep$sd[par_index]
    key = fit$conf[[info[[i]][2]]]
    is_parametric = as.logical(fit$conf[[info[[i]][[3]]]])
    if (is_parametric) {
      count = 0
      for (j in seq_len(nrow(key))) {
        for (l in seq_len(ncol(key))) {
          if (key[j, l] != -1) {
            key[j, l] = count
            count = count + 1
          }
        }
      }
    }
    for (j in seq_len(nrow(key))) {
      for (l in seq_len(ncol(key))) {
        if (key[j, l] != -1) {
          out[[i]][[length(out[[i]]) + 1]] = data.frame(
            value = par[key[j, l] + 1],
            sd = par_sd[key[j, l] + 1],
            a = fit$conf$minAge + l - 1,
            fleet = j,
            parametric = is_parametric,
            tag = names(info)[i])
        }
      }
    }
  }
  do.call(rbind, unlist(out, recursive = FALSE))
}

# Turn an R plot into a beautiful pdf made by LaTeX and TikZ,
# using the tikzDevice package
plot_tikz = function(plot = NULL,
                     expression = NULL,
                     file = "Rplots.pdf",
                     extra_packages = NULL,
                     tex_engine = c("pdflatex", "lualatex"),
                     ...) {
  expression = substitute(expression)
  if (is.null(plot) && is.null(expression)) {
    stop("Either `plot` or `expression` must be non-NULL")
  }

  # Create a temporary file for the tikz-output
  tmp = tempfile(tmpdir = getwd())
  # Clean up after yourself on early interrupt
  on.exit(unlink(tmp))

  # I am nor sure about what is going on here...
  opt = options()
  on.exit(options(opt), add = TRUE) #Reset global options on exit
  tikzDevice::setTikzDefaults(overwrite = FALSE)

  # Extract default tex usepackages and extra packages
  tex_packages = options()$tikzLatexPackages
  extra_packages = c(extra_packages, "bm", "amsmath", "amssymb")
  extra_packages = paste0("\\usepackage{", extra_packages, "}\n")
  tex_packages = unique(c(tex_packages, extra_packages))

  # Open a device for creating a tex-file
  tikzDevice::tikz(tmp, standAlone = TRUE, packages = tex_packages, ...)
  # Call dev.off() on exit in case of interruptions
  current_device = dev.cur()
  on.exit({
    if (current_device %in% dev.list()) dev.off(current_device)
  }, add = TRUE)

  # Plot something into the tex-file
  if (!is.null(plot)) {
    if (any(class(plot) %in% c("gg", "ggplot", "patchwork"))) {
      print(plot)
    } else {
      for (p in plot) print(p)
    }
  } else {
    eval(expression)
  }

  # Finish the creation of the tex-file
  dev.off()

  # Compile to pdf using lualatex
  system2(tex_engine[1], shQuote(tmp))

  # Copy pdf file to final destination
  file.copy(paste0(tmp, ".pdf"), file, overwrite = TRUE)

  # Clean up all temporary files
  tmp_filename = tail(strsplit(tmp, "/")[[1]], 1)
  files_to_clean = grep(tmp_filename, list.files(full.names = TRUE), value = TRUE)
  unlink(files_to_clean)

  # Silently return a TRUE to show that the function completed successfully
  invisible(TRUE)
}
