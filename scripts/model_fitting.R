library(stockassessment)
library(RhpcBLASctl)
library(here)

# ==============================================================================
# Paths and other stuff
# ==============================================================================

# Path to the directory where we store our results
sam_dir = here::here()
results_dir = file.path(sam_dir, "raw_data", "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# The minimum number of years required to not delete a survey fleet during forward-validation
min_obs_years_forward = 4
# The number of parallel cores to be used for cross/forward-validation
n_cores = 10

# Should we overwrite the results if they already exist?
overwrite_fits = FALSE # Override model fits for each of the 17 fish stocks?
overwrite_eval = FALSE # Override cross/forward validation results?

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

# Loop over all the fish stock data sets, fit our four competing models to the data,
# and then perform cross-validation and forward-validation
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
    res = list(fits = list(), forward = NULL, crossval = NULL)
  }
  if (overwrite_eval) res$forward = res$crossval = NULL
  if (overwrite_fits) res$fits = list()


  # Fit the official model to the data.
  # This is done using the fitfromweb function, to download the official model configuration
  # from stockassessment.org. Then we use the refit function to estimate all model parameters
  if (!"official" %in% names(res$fits)) {
    message("Start fitting the official model")
    res$fits$official = stockassessment::fitfromweb(stocks[[i]]$webname, character.only = TRUE)
    res$fits$official = stockassessment:::refit(res$fits$official, silent = TRUE, trace = 0)
  }

  data = res$fits$official$data
  conf = res$fits$official$conf

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
      fit = my_sam_fit(data, conf, par)
      # If the model fitting failed, we just return a list containing the data and conf
      # necessary for trying to estimate the parameters. This is later used
      # when performing cross/forward-validation for different subsets of the data
      if (is.null(fit)) fit = list(data = data, conf = conf)
      fit
    })
  }

  # Fit the spline1 model to the data.
  # This is done by computing the spline matrices X and S, and adding them to the conf,
  # in addition to setting/changing some other conf options.
  # The computation of X and S happens inside the add_spline_info_to_conf() function
  if (!"spline1" %in% names(res$fits)) {
    message("Start fitting spline 1")
    res$fits$spline1 = local({
      # Compute and add all necessary spline matrices to the conf
      conf = add_spline_info_to_conf(
        conf = conf,
        recipe_func = get_spline1_conf_parts,
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
      fit = my_sam_fit(data, conf, par)
      # If the model fitting failed, we just return a list containing the data and conf
      # necessary for trying to estimate the parameters. This is later used
      # when performing cross/forward-validation for different subsets of the data
      if (is.null(fit)) fit = list(data = data, conf = conf)

      fit
    })
  }

  # Fit the spline2 model to the data.
  # This is done by computing the spline matrices X and S, and adding them to the conf,
  # in addition to setting/changing some other conf options.
  # The computation of X and S happens inside the add_spline_info_to_conf() function
  if (!"spline2" %in% names(res$fits)) {
    message("Start fitting spline 2")
    res$fits$spline2 = local({
      # Compute and add all necessary spline matrices to the conf
      conf = add_spline_info_to_conf(
        conf = conf,
        recipe_func = get_spline2_conf_parts,
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
      fit = my_sam_fit(data, conf, par)
      # If the model fitting failed, we just return a list containing the data and conf
      # necessary for trying to estimate the parameters. This is later used
      # when performing cross/forward-validation for different subsets of the data
      if (is.null(fit)) fit = list(data = data, conf = conf)

      fit
    })
  }

  # Perform the forward-validation study
  res$forward = do_evaluation(
    res$fits,
    eval_years = tail(sort(data$years), floor(length(data$years) / 3)),
    min_obs_years = min_obs_years_forward,
    n_cores = n_cores,
    forward = TRUE,
    eval_data = res$forward
  )

  # Perform the cross-validation study
  res$crossval = do_evaluation(
    res$fits,
    eval_years = tail(sort(data$years), -1),
    n_cores = n_cores,
    forward = FALSE,
    eval_data = res$crossval
  )

  # The fitted SAM models require a lot of memory to save, because they contain
  # an entire R environment that is full of large objects. We remove this R environment
  # before saving the SAM models
  for (j in seq_along(res$fits)) {
    res$fits[[j]]$obj = NULL
  }

  saveRDS(res, out_path)
}

