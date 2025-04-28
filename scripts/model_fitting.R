library(stockassessment)
library(RhpcBLASctl)
library(here)

# ==============================================================================
# Paths and other stuff
# ==============================================================================

# The minimum number of years required to not delete a survey fleet during forward-validation
min_obs_years_forward = 4
# The number of parallel cores to be used for cross/forward-validation
n_cores = 10

n_sims = 100

rel.tol = 1e-8

age_transform = function(a) sqrt(a)

# Path to the directory where we store our results
sam_dir = here::here()
results_dir = file.path(sam_dir, "raw_data", "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Should we overwrite the results if they already exist?
overwrite_fits = FALSE # Overwrite model fits for each of the 17 fish stocks?
overwrite_eval = FALSE # Overwrite cross/forward validation results?
overwrite_sim = FALSE # Overwrite results from the simulation study?

# Ensure that only one thread is used by OMP and BLAS
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)

# ==============================================================================
# Load all the helper functions
# ==============================================================================
function_path = file.path(sam_dir, "scripts", "helper_functions.R")
source(function_path)

# ==============================================================================
# Modelling
# ==============================================================================

# This is the list of all fish stocks we want to fit our models to.
# webname is the name of the data set at stockassessment.org.
# filename is the name we use when saving our results for that fish stock
stocks = list(
  list(
    filename = "NEA_cod_2023.rds",
    webname = "NEA_cod_2023_final_run"
  ),
  list(
    filename = "NEA_saithe_2022.rds",
    webname = "NEAsaithe_2022_v3"
  ),
  list(
    filename = "NS_whiting_2023.rds",
    webname = "NSwhiting_2023"
  ),
  list(
    filename = "NS_sole_2022.rds",
    webname = "Sole20_24_2022vs21"
  ),
  list(
    filename = "NS_plaice_2022.rds",
    webname = "plaice_final_10fix"
  ),
  list(
    filename = "NS_haddock_2023.rds",
    webname = "NShaddock_WGNSSK2023_Run1"
  ),
  list(
    filename = "BS_plaice_2022.rds",
    webname = "ple.27.21-23_WGBFAS_2023_ALT_v1"
  ),
  list(
    filename = "BS_herring_2022.rds",
    webname = "GoR_BP_v2.2.3qF_s"
  ),
  list(
    filename = "BS_cod_2021.rds",
    webname = "WBcod22Fsq"
  ),
  list(
    filename = "CSWS_whiting_2023.rds",
    webname = "whg.7b-ce-k_WGCSE22_RevRec_2023"
  ),
  list(
    filename = "CSWS_cod_2022.rds",
    webname = "Cod_7ek_2023"
  ),
  list(
    filename = "CSWS_haddock_2023.rds",
    webname = "HAD7bk_2023_final"
  ),
  list(
    filename = "CSWS_plaice_2022.rds",
    webname = "Ple.7fg.2022.main"
  ),
  list(
    filename = "F_saithe_2023.rds",
    webname = "fsaithe-NWWG-2023"
  ),
  list(
    filename = "F_haddock_2023.rds",
    webname = "NWWG2023_faroehaddock"
  ),
  list(
    filename = "F_ling_2023.rds",
    webname = "lin.27.5b_wgdeep2023_final"
  ),
  list(
    filename = "WMS_bw_2023.rds",
    webname = "BW-2023"
  )
)

# Loop over all the fish stock data sets and fit our four competing models to the data
start_time = Sys.time()
for (i in seq_along(stocks)) {

  # Print the progress so far
  time_passed = Sys.time() - start_time
  message(
    "\n==============================================================================\n",
    "Start on fish stock nr. ", i, " / ", length(stocks),
    "\nTime passed since start: ", round(as.numeric(time_passed), 2), " ", attr(time_passed, "units"),
    "\n==============================================================================\n"
  )

  # Check if the results already exist, and if we should use them or not
  out_path = file.path(results_dir, stocks[[i]]$filename)
  if (file.exists(out_path)) {
    res = readRDS(out_path)
  } else {
    res = list(fits = list(), forward = NULL, crossval = NULL, simulations = NULL)
  }
  if (overwrite_eval) res$forward = res$crossval = NULL
  if (overwrite_fits) res$fits = list()

  # Fit the final model to the data.
  # This is done using the fitfromweb function, to download the final model configuration
  # from stockassessment.org. Then we use the refit function to estimate all model parameters
  if (!"final" %in% names(res$fits)) {
    message("Start fitting the final model")
    res$fits$final = stockassessment::fitfromweb(stocks[[i]]$webname, character.only = TRUE)
    res$fits$final = stockassessment:::refit(res$fits$final, silent = TRUE, trace = 0)
  }

  data = res$fits$final$data
  conf = res$fits$final$conf

  # Fit the maximal model to the data.
  # This is done by simply changing the keyLogFpar and keyVarObs objects of
  # the conf so that every non-negative element has a unique value
  if (!"maximal" %in% names(res$fits)) {
    message("Start fitting the maximal model")
    res$fits$maximal = local({
      # Update conf$keyLogFpar and conf$keyVarObs
      for (name in c("keyLogFpar", "keyVarObs")) {
        count = 0
        key = conf[[name]]
        for (i in seq_len(nrow(key))) {
          for (j in seq_len(ncol(key))) {
            if (key[i, j] != -1) {
              key[i, j] = count
              count = count + 1
            }
          }
        }
        conf[[name]] = key
      }
      # Get the default initial values for the current combination of data and conf
      par = stockassessment::defpar(data, conf)
      # Estimate model parameters
      fit = my_sam_fit(data, conf, par, rel.tol = rel.tol)
      # If the model fitting failed, we just return a list containing the data and conf
      # necessary for trying to estimate the parameters. This is later used
      # when performing cross/forward-validation for different subsets of the data
      if (is.null(fit)) fit = list(data = data, conf = conf)
      fit
    })
  }

  # Fit the spline model to the data.
  # This is done by computing the spline matrices X and S, and adding them to the conf,
  # in addition to setting/changing some other conf options.
  # The computation of X and S happens inside the add_spline_info_to_conf() function
  if (!"spline" %in% names(res$fits)) {
    message("Start fitting the spline model")
    res$fits$spline = local({
      # Compute and add all necessary spline matrices to the conf
      conf = add_spline_info_to_conf(
        conf = conf,
        recipe_func = get_spline_conf_parts,
        modify_S = modify_S,
        age_transform = age_transform,
        data = data
      )

      # This tells us that we should use splines for logFpar
      conf$useLogFparSpline = 1
      # This tells us that we should penalise the splines for logFpar
      conf$useLogFparSplinePenalty = 1
      # This tells us that we should use splines for logSdLogObs
      conf$useVarObsSpline = 1
      # This tells us that we should penalise the splines for logSdLogObs
      conf$useVarObsSplinePenalty = 1
      # This tells us that the logFpar spline coefficients should be modelled as latent parameters
      conf$latentLogFparSpline = 1
      # This tells us the max value of the log of the smoothness penalty parameters before the
      # smoothness penalty prior should start to heavily penalise large smoothness penalty values
      conf$splinePenaltyMax = 7
      conf$splinePenaltySpeed = 100

      # Get the default initial values for the current combination of data and conf
      par = stockassessment::defpar(data, conf)
      # Estimate model parameters
      fit = my_sam_fit(data, conf, par, rel.tol = rel.tol)
      # If the model fitting failed, we just return a list containing the data and conf
      # necessary for trying to estimate the parameters. This is later used
      # when performing cross/forward-validation for different subsets of the data
      if (is.null(fit)) fit = list(data = data, conf = conf)

      fit
    })
  }

  message("Done! Classes: ", paste(sapply(res$fit, class), collapse = ", "), "\n")

  saveRDS(res, out_path)
}

# Perform forward- and cross-validation
start_time = Sys.time()
for (i in seq_along(stocks)) {
  
  # Print the progress so far
  time_passed = Sys.time() - start_time
  message(
    "\n==============================================================================\n",
    "Start on fish stock nr. ", i, " / ", length(stocks),
    "\nTime passed since start: ", round(as.numeric(time_passed), 2), " ", attr(time_passed, "units"),
    "\n==============================================================================\n"
  )

  out_path = file.path(results_dir, stocks[[i]]$filename)
  res = readRDS(out_path)
  data = res$fits[[1]]$data
  if (overwrite_eval) res$forward = res$crossval = NULL

  # Perform the forward-validation study
  res$forward = do_evaluation(
    res$fits,
    eval_years = tail(sort(data$years), floor(length(data$years) / 3)),
    min_obs_years = min_obs_years_forward,
    n_cores = n_cores,
    forward = TRUE,
    rel.tol = rel.tol,
    eval_data = res$forward
  )

  # Perform the cross-validation study
  res$crossval = do_evaluation(
    res$fits,
    eval_years = tail(sort(data$years), -1),
    n_cores = n_cores,
    forward = FALSE,
    rel.tol = rel.tol,
    eval_data = res$crossval
  )

  saveRDS(res, out_path)
}

start_time = Sys.time()
for (i in seq_along(stocks)) {
  
  # Print the progress so far
  time_passed = Sys.time() - start_time
  message(
    "\n==============================================================================\n",
    "Start on fish stock nr. ", i, " / ", length(stocks),
    "\nTime passed since start: ", round(as.numeric(time_passed), 2), " ", attr(time_passed, "units"),
    "\n==============================================================================\n"
  )

  out_path = file.path(results_dir, stocks[[i]]$filename)
  res = readRDS(out_path)

  if (overwrite_sim || is.null(res$simulations)) {
    res$simulations = list()
    sim_start_time = Sys.time()
    set.seed(1)
    for (j in seq_along(res$fits)) {
      if (class(res$fits[[j]]) != "sam") next
      sims = simulate(res$fits[[j]], n_sims)
      logssb = res$fits[[j]]$sdrep$value[names(res$fits[[j]]$sdrep$value) == "logssb"]
      logfbar = res$fits[[j]]$sdrep$value[names(res$fits[[j]]$sdrep$value) == "logfbar"]

      res$simulations[[j]] = parallel::mclapply(
        X = seq_len(n_sims),
        mc.cores = n_cores,
        mc.preschedule = FALSE,
        FUN = function(l) {
          my_fits = lapply(
            X = res$fits,
            FUN = function(x) {
              par = stockassessment::defpar(sims[[l]], x$conf)
              my_sam_fit(sims[[l]], x$conf, par, rel.tol = rel.tol)
            })
          out = lapply(
            X = names(res$fits),
            FUN = function(name) {
              out = data.table(
                model = name,
                year = rep(res$fits[[name]]$data$years, 2),
                sim_nr = l,
                simulation_model = names(res$fits)[j],
                truth = c(logssb, logfbar),
                tag = rep(c("logssb", "logfbar"), c(length(logssb), length(logfbar)))
              )
              if (is.null(my_fits[[name]])) return(out)
              x = my_fits[[name]]
              ssb_index = which(names(x$sdrep$value) == "logssb")
              fbar_index = which(names(x$sdrep$value) == "logfbar")
              out[, let(
                convergence = x$opt$nlminb_convergence,
                pred = x$sdrep$value[c(ssb_index, fbar_index)],
                sd = x$sdrep$sd[c(ssb_index, fbar_index)]
              )]
              out
            })
          out = rbindlist(out, fill = TRUE)
          time_passed = difftime(Sys.time(), sim_start_time)
          message(
            "stock nr: ", i, "/", length(stocks),
            ", fit_nr: ", j, "/", length(res$fits),
            ", sim: ", l, "/", n_sims,
            ", time passed: ", round(time_passed, 2), " ", attr(time_passed, "units")
          ) 
          out
        }
      )
      res$simulations[[j]] = rbindlist(res$simulations[[j]], fill = TRUE)

    }
    res$simulations = rbindlist(res$simulations, fill = TRUE)
  }

  saveRDS(res, out_path)
}
