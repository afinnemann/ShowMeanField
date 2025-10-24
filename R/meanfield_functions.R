

#' Computes mean-field solutions
#'
#' This function computes the mean-field solutions for various models as specified. It allows for the customization
#' of the beta parameter space, external field parameter space, and the selection of the mean-field model to solve.
#' Valid model options include "Ising", "Percolation", "han", "zero_corrected", "independent_a", "Blume_capel",
#' "Potts", and "Continuous_01".
#'
#' @param alpha_sim Numeric vector, external field parameter space.
#' @param beta_sim Numeric vector, beta parameter space.
#' @param average_density Numeric, average density for the model.
#' #' @param MFA_to_solve Character string, specifies the mean-field model to solve.
#'   Options:
#'   \itemize{
#'     \item \code{"Ising"}: Classical Ising mean-field with binary spins,
#'       \code{alpha} as external field, \code{beta} as inverse temperature.
#'     \item \code{"Percolation"}: Site percolation mean-field,
#'       \code{alpha} shifts activation, \code{beta} controls occupation probability.
#'     \item \code{"han"}: Han’s variant mean-field with centered external field.
#'     \item \code{"zero_corrected"}: Mean-field with correction for zero-density terms.
#'     \item \code{"independent_a"}: Independent-site approximation,
#'       \code{alpha} acts as penalty term, \code{beta} scales density effect.
#'     \item \code{"Blume_capel"}: Spin-1 Blume–Capel model,
#'       \code{alpha} is crystal-field parameter, \code{beta} inverse temperature.
#'     \item \code{"Potts"}: q-state Potts model; here \code{alpha} represents the number
#'       of categories (\eqn{q}, must be integer > 2), \code{beta} is inverse temperature.
#'     \item \code{"Continuous_01"}: Continuous-spin Ising variant on [0,1],
#'       \code{alpha} shifts mean activity, \code{beta} controls steepness.
#'       #'     \item \code{"Continuous_-11"}: Continuous-spin Ising variant on [-1,1],
#'       \code{alpha} shifts mean activity, \code{beta} controls steepness.
#'   }
#' @return A list or data frame with the simulation results.
#' @examples
#' MFA_sim(alpha_sim = seq(-3,3, length.out = 100), beta_sim = seq(0,3, length.out = 100), average_density = 2, MFA_to_solve = "Ising")
#' @export

MFA_sim <- function(
    alpha_sim = NA,
    beta_sim = NA,
    average_density = 2,
    gamma = 0.1,
    MFA_to_solve = c("Ising", "Percolation",
                     #"han", "zero_corrected", "independent_a",
                     "Blume_capel", "Potts", "Continuous_01", "Continuous_-11",
                     "Blume_capel_inf","Blume_capel_0inf")
){

  # Validate choice
  allowed <- c("Ising","Percolation",
               #"han","zero_corrected","independent_a",
               "Blume_capel","Potts","Continuous_01","Continuous_-11","Blume_capel_inf", "Blume_capel_0inf")
  if (!(MFA_to_solve %in% allowed)) stop("Incorrect mean-field selected")

  if (0 %in% beta_sim) stop("Error: beta has to be larger than zero.")

  if (MFA_to_solve == "Potts") {
    if (any(alpha_sim < 2) || any(alpha_sim %% 1 != 0)) {
      stop("Alpha indicates the number of categories for the Potts model and must be integers > 2")
    }
  }

  # Root functions (mean-field consistency equations)
  if ((!is.na(alpha_sim[1]) & !is.na(beta_sim[1])) & !is.na(average_density[1])) {

    if (MFA_to_solve == "Ising") {
      root_MFA <- function (x, average_density, alpha, beta)
        ((exp(2 * beta * (average_density * x + alpha)) - 1) /
           (1 + exp(2 * beta * (average_density * x + alpha)))) - x

    } else if (MFA_to_solve == "Percolation") {
      root_MFA <- function (x, average_density, alpha, beta)
        (exp(beta * (average_density * x + alpha)) /
           (1 + exp(beta * (average_density * x + alpha)))) - x

    } else if (MFA_to_solve == "han") {
      root_MFA <- function (x, average_density, alpha, beta)
        (exp(beta * (average_density * x + (alpha - 1) / 2)) /
           (1 + exp(beta * (average_density * x + (alpha - 1) / 2)))) - x

    } else if (MFA_to_solve == "zero_corrected") {
      root_MFA <- function (x, average_density, alpha, beta)
        (exp(beta * (average_density * x + alpha - 0.5 * average_density)) /
           (1 + exp(beta * (average_density * x + alpha - 0.5 * average_density)))) - x

    } else if (MFA_to_solve == "independent_a") {
      root_MFA <- function (x, average_density, alpha, beta)
        (exp(beta * average_density * x - alpha - 0.5 * average_density) /
           (1 + exp(beta * average_density * x - alpha - 0.5 * average_density))) - x

    } else if (MFA_to_solve == "Blume_capel") {
      root_MFA <- function (x, average_density, alpha, beta)
        (exp(beta * ((average_density * x) - alpha)) -
           exp(beta * ((-average_density * x) - alpha))) /
        (exp(beta * (average_density * x - alpha)) +
           exp(beta * (-average_density * x - alpha)) + 1) - x

    } else if (MFA_to_solve == "Potts") {
      root_MFA <- function (x, average_density, alpha, beta) {
        eps   <- 1e-8
        denom <- alpha - 1
        a     <- ifelse(abs(denom) < eps, sign(denom) * eps + eps, denom)
        z     <- (beta * (1 - alpha * x) / a)
        z     <- pmax(pmin(z, 700), -700)          # prevent overflow
        p     <- 1 / (1 + a * exp(z))
        p     <- pmin(pmax(p, 1e-8), 1 - 1e-8)     # clamp to (0,1)
        p - x
      }

    } else if (MFA_to_solve == "Continuous_01") {
      root_MFA <- function (x, average_density, alpha, beta) {
        theta <- beta * (alpha + average_density * x)


            val <-   ((1/theta) * (((theta * exp(theta))/(expm1(theta)))-1))-x
            #val <- #(1/(-expm1(-theta)) - 1/theta) - x
      }

    } else if (MFA_to_solve == "Continuous_-11") {
      root_MFA <- function (x, average_density, alpha, beta) {
        theta <- beta * (average_density * x + alpha)
        val <- (cosh(theta) / sinh(theta)) - (1/theta)   # μ = coth(θ) − 1/θ
        val - x
      }

    } else if (MFA_to_solve == "Blume_capel_0inf") {
      if (any(alpha_sim <= 0)) stop("'Blume_capel_0inf' needs alpha (gamma) > 0.")
      terf <- function(z) 2 * pnorm(z * sqrt(2), lower = FALSE)
      root_MFA <- function (x, average_density, alpha, beta) {

        theta <- beta * (average_density * x + 1)
        gamma = alpha

        phi <- (1 / (2 * sqrt(pi * gamma))) * exp(-(theta^2) / (4 * gamma))

        psi <- 1 + terf(x / (2 * sqrt(gamma)))

        mu  <- phi / psi + theta / (2 * gamma)
        mu - x
      }

    }else if (MFA_to_solve == "Blume_capel_inf") {
      #if (any(alpha_sim <= 0)) stop("For 'Blume_capel_0inf', alpha (gamma) must be > 0 to avoid division by zero and keep the model well-defined.")
      root_MFA <- function (x, average_density, alpha, beta) {
        theta <- beta * (average_density * x + 0.3) #this zero is external field
        delta <- alpha #alpha becomes delta here
        val <- (theta / (2 * delta)) * exp((theta^2) / (2 * delta))
        val - x
      }
    }
    # run solver
    result <- MFA_solver(
      alpha_sim       = alpha_sim,
      beta_sim        = beta_sim,
      average_density = average_density,
      root_MFA        = root_MFA,
      MFA_to_solve    = MFA_to_solve
    )
  } else {
    stop("Missing parameters for MFA simulation.")
  }

  return(result)
}



#' Internal function that helps solving the mean field
#'
#' This internal function assists in solving the mean field by employing simulations based on specified parameters.
#' It is primarily used within the MFA_sim function.
#'
#' @param root_MFA Function, the root mean-field function to be solved.
#' @param alpha_sim Numeric, external field parameter space.
#' @param beta_sim Numeric, beta parameter space.
#' @param average_density Numeric, average density for the model.
#' @param MFA_to_solve Character, specifies the mean-field model being solved.
#' @return A data frame with the simulation results, including stability and value of each root found.
#' @export
#' @examples
#' # Example usage not provided as this function is intended for internal use by MFA_sim()


MFA_solver <- function(root_MFA,
                       alpha_sim = 0,
                       beta_sim = 1,
                       average_density = 4,
                       MFA_to_solve = MFA_to_solve
                       ){


  root_MFA_as.list = list(root_MFA)

    result =  parSim::parSim(
      ### SIMULATION CONDITIONS
      alpha = alpha_sim,
      beta = beta_sim,
      average_density = average_density,
      root_MFA_as.list = root_MFA_as.list,
      MFA_to_solve = MFA_to_solve,
      reps = 1, # repetitions per condition
      write = FALSE, # Writing to a file
      nCores = 1, # Number of cores to use
      expression = {


        solve_MFA <- function(x)  root_MFA_as.list[[1]](x, average_density, alpha, beta)

        if (MFA_to_solve == "Blume_capel_inf"){
          interval = c(-2, 2)
        #} else if (MFA_to_solve == "Blume_capel_0inf"){
        #  interval = c(0, 1)
        } else if (MFA_to_solve == "Blume_capel_0inf"){
        interval = c(0, 10)
          }else {
          interval = c(-1, 1)
        }


        all_roots <- try(rootSolve::uniroot.all(solve_MFA, interval = interval), silent = TRUE)
        if (inherits(all_roots, "try-error") || length(all_roots) == 0L) return(data.frame(stability = NA_character_, value = NA_real_))


        # Determine the stability of each root by testing if the function is positive or negative prior to the root
        res <- data.frame(stability = sapply(all_roots, function(x) ifelse(solve_MFA(x - 0.001) > 0, "stable", "unstable")),
                          value = all_roots)
        #res <- data.frame(value = all_roots)
      }
    )
    result$value = as.numeric(as.character(result$value))
    return(result)
}



#' Plots a MFA_sim() object
#'
#' This function creates a 3D plot of the results from the MFA_sim() function. It allows customization of camera angles
#' and axis labels for detailed visualization.
#'
#' @param result The output of the MFA_sim() function, expected to be a data frame.
#' @param camera_x Numeric, controls the camera angle of the x-axis.
#' @param camera_y Numeric, controls the camera angle of the y-axis.
#' @param camera_z Numeric, controls the camera angle of the z-axis.
#' @param x_label Character, label of the x-axis.
#' @param y_label Character, label of the y-axis.
#' @param z_label Character, label of the z-axis.
#' @return A 3D plotly object.
#' @export
#' @examples
#' result <- MFA_sim(alpha_sim = seq(-3,3, length.out = 100), beta_sim = seq(0,3, length.out = 100), average_density = 2, MFA_to_solve = "Ising")
#' plot_MFA(result)
#'
plot_MFA <- function(result,
                     camera_x = 2,
                     camera_y = 0,
                     camera_z = 0,
                     y_label = "Alpha",
                     x_label = "Beta",
                     z_label = "Mean field"){

  result <- result |> dplyr::filter(!is.na(value))

  #Coloring based on stability and magnitude
  reds <- grDevices::colorRampPalette(c("lightpink", "darkred"))(100) #color pallete for unstable
  purp <- grDevices::colorRampPalette(c("lavender", "darkorchid4"))(100) #color palette for stable

  # Creating a color mapping function
  color_reds <- scales::col_numeric(palette = reds, domain = c(min(result$value), max(result$value)))
  color_purp <- scales::col_numeric(palette = purp, domain = c(min(result$value), max(result$value)))

  # Apply the color function based on stability
  result$color <- ifelse(result$stability == "stable", color_purp(result$value), color_reds(result$value))

  plot <- plotly::plot_ly(z = ~result$value,
                  x = ~result$beta,
                  y = ~result$alpha,
                  mode="markers",
                  type = "scatter3d",
                  marker=list(color=result$color, #~result$stability,
                              size = 6,
                              #colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                              showscale=TRUE)
                  #line=list(width=2,color='DarkSlateGrey'))
  ) %>%
    #layout(scene = list(aspectratio = list(x = 1, y = 1, z = 0.6))) %>%
    plotly::layout(scene = list(xaxis = list(title = x_label),
                        yaxis = list(title = y_label),
                        zaxis = list(title = z_label),
                        aspectratio = list(x = 1, y = 1, z = 0.6),
                        camera = list(eye = list(x = camera_x, y = camera_y, z = camera_z)))) %>%
    plotly::hide_colorbar()

  return(plot)
}


#' Saves the output of plot_MFA() as a PDF file
#'
#' This function saves the 3D plot generated by `plot_MFA()` as a PDF file. It requires the setup of a Python environment
#' using `reticulate`. Before using this function, ensure Python is configured by running `reticulate::use_miniconda('r-reticulate')`
#' and `reticulate::py_run_string("import sys")`. These steps are necessary to enable R to interface with Python for saving
#' the plot because the plot saving functionality relies on Python libraries.
#'
#' @param MFA_plot The plotly object output from `plot_MFA()` to be saved.
#' @param width Numeric, the width of the saved plot in pixels. Defaults to 600.
#' @param height Numeric, the height of the saved plot in pixels. Defaults to 600.
#' @param saved_image_name Name of the saved file ending with ".pdf".
#' @export
#' @examples
#' # First, generate a plot with plot_MFA() function.
#' result <- MFA_sim(alpha_sim = seq(-3,3, length.out = 100), beta_sim = seq(0,3, length.out = 100), average_density = 2, MFA_to_solve = "Ising")
#' MFA_plot <- plot_MFA(result)
#' # Then, save the generated plot as a PDF file.
#' reticulate::use_miniconda('r-reticulate')
#' reticulate::py_run_string("import sys")
#' save_plot_MFA(MFA_plot)

save_plot_MFA <- function(MFA_plot,
                          saved_image_name,
                          width = 600,
                          height = 600){
  reticulate::use_miniconda('r-reticulate')
  reticulate::py_run_string("import sys")
  plotly::save_image(p = MFA_plot, file = saved_image_name, width = width, height = height)

}




