

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
#' @param MFA_to_solve Character vector, specifies the mean-field model to solve. 
#'        Options: "Ising", "Percolation", "han", "zero_corrected", "independent_a", "Blume_capel", "Potts", "Continuous_01".
#' @return A list or data frame with the simulation results.
#' @examples
#' MFA_sim(alpha_sim = seq(-3,3, length.out = 100), beta_sim = seq(0,3, length.out = 100), average_density = 2, MFA_to_solve = "Ising")
#' @export

MFA_sim <- function(alpha_sim = NA,
                    beta_sim = NA,
                    average_density = 2,
                    MFA_to_solve = c("Ising", "Percolation", "han", 
                                     "zero_corrected", "independent_a",
                                     "Blume_capel", "Potts",  "Continuous_01")){
  
  
  #Errors
  if (!(MFA_to_solve %in% c("Ising", "Percolation", "han", 
                            "zero_corrected", "independent_a",
                            "Blume_capel", "Potts",  "Continuous_01"))) stop("Incorrect mean-field selected")
  
  
  #The root function is zero when the mean field is equal to the mean field.
  if(MFA_to_solve != "custom"){ #custom MFA?
    if ((!is.na(alpha_sim[1]) & !is.na(beta_sim[1])) & !is.na(average_density[1])){ #custom parameters? 
      
      
      if (MFA_to_solve == "Ising"){
        root_MFA = function (x, average_density, alpha, beta) ((exp(2 * beta * (average_density * x + alpha)) - 1) / (1 + exp(2 * beta * (average_density * x + alpha)))) - x
      } else if (MFA_to_solve == "Percolation"){
        root_MFA = function (x, average_density, alpha, beta) (exp(beta * (average_density * x + alpha)) / (1 + exp(beta * (average_density * x + alpha)))) - x
      } else if (MFA_to_solve == "han"){
        root_MFA = function (x, average_density, alpha, beta) (exp(beta * (average_density * x + (alpha - 1) / 2)) / (1 + exp(beta * (average_density * x + (alpha - 1) / 2)))) - x
      } else if (MFA_to_solve == "zero_corrected"){
        root_MFA = function (x, average_density, alpha, beta) (exp(beta * (average_density * x + alpha - 0.5 * average_density)) / (1 + exp(beta * (average_density * x + alpha - 0.5 * average_density)))) - x
      } else if (MFA_to_solve == "independent_a"){
        root_MFA = function (x, average_density, alpha, beta) (exp(beta * average_density * x - alpha - 0.5 * average_density) / (1 + exp(beta * average_density * x - alpha - 0.5 * average_density))) - x
      } else if (MFA_to_solve == "Blume_capel"){
        #For Blume capel, beta is the coupling constant, and alpha is the costs to presence. Temp is assumed to be 1
        root_MFA = function (x, average_density, alpha, beta) (exp(beta * ((average_density * x) - alpha)) - exp(beta * (( - average_density * x) - alpha))) / (exp(beta * ((average_density * x - alpha))) + exp( beta * ((- average_density * x - alpha))) + 1) - x
      } else if (MFA_to_solve == "Potts"){
        #For Blume capel, beta is the coupling constant, and alpha is the costs to presence. Temp is assumed to be 1
        root_MFA = function (x, average_density, alpha, beta) (1 /(1 + (alpha - 1.0001) * exp((beta * (1 - alpha * x)/(alpha - 1.0001))))+0.0001) - x
      } else if (MFA_to_solve == "Continuous_01"){
        #For Blume capel, beta is the coupling constant, and alpha is the costs to presence. Temp is assumed to be 1
        root_MFA = function (x, average_density, alpha, beta){
          # Define theta
          theta = beta * (alpha + average_density * x) 
          
          # continuous Ising in form: 0 = equation - x
          result = ((1 / theta) * (1 - (theta * exp(-theta)) / (1 - exp(-theta)))) - x
          
          return(result)
        }
      }
      
      #default MFA with custom parameters
      result <- MFA_solver(alpha_sim = alpha_sim ,
                           beta_sim = beta_sim,
                           average_density = average_density,
                           root_MFA = root_MFA,
                           MFA_to_solve = MFA_to_solve)
      
    }else{
      #load one of the default data sets to avoid simulating
      
      if (MFA_to_solve == "Blume_capel"){
        result <- read.csv("Blume_capel_default.csv")
      }
    }
  }else{
    #Custom H with custom parameters  
    
    result <- MFA_sim(alpha_sim = ifelse(is.null(custom_pars[[1]]),NA,custom_pars[[1]]),
                      beta_sim = ifelse(is.null(custom_pars[[2]]),NA,custom_pars[[2]]),
                      average_den = ifelse(is.null(custom_pars[[3]]),NA,custom_pars[[3]]),
                      root_MFA = custom_MFA,
                      MFA_to_solve = MFA_to_solve)
    
  }
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
                       MFA_to_solve = MFA_to_solve){
  
  
  root_MFA_as.list = list(root_MFA)
  
  if(MFA_to_solve != "custom"){
    
    result =  parSim::parSim(
      ### SIMULATION CONDITIONS
      alpha = alpha_sim,
      beta = beta_sim, 
      average_density = average_density,
      root_MFA_as.list = root_MFA_as.list,
      reps = 1, # repetitions per condition
      write = FALSE, # Writing to a file
      nCores = 1, # Number of cores to use
      expression = {
        
        #print(alpha)
        #print(beta)
        solve_MFA <- function(x)  root_MFA_as.list[[1]](x, average_density, alpha, beta)
        #print(beta * (alpha + average_density * x)+ 0.0001)
        # Find all roots of the root function in the interval [-1,1]
        all_roots <- rootSolve::uniroot.all(solve_MFA,interval = c(-1, 1))
        
        # Determine the stability of each root by testing if the function is positive or negative prior to the root
        res <- data.frame(stability = sapply(all_roots, function(x) ifelse(solve_MFA(x - 0.001) > 0, "stable", "unstable")),
                          value = all_roots)
        #res <- data.frame(value = all_roots)
      }
    )
    result$value = as.numeric(as.character(result$value))
    return(result)
    
    
  }else if(MFA_to_solve == "custom"){
    result =  parSim(
      ### SIMULATION CONDITIONS
      alpha = alpha_sim,
      beta = beta_sim, 
      average_density = average_density,
      root_MFA_as.list = root_MFA_as.list,
      #root_MFA = root_MFA(),
      reps = 1, # repetitions per condition
      write = FALSE, # Writing to a file
      nCores = 1, # Number of cores to use
      expression = {
        
        #solve_MFA <- function(x)  root_MFA_as.list[[1]](x, param = list(alpha,beta,average_density))
        solve_MFA <- function(x)  root_MFA_as.list[[1]](x, average_density = average_density,alpha = alpha, beta = sim_df)
        
        # Find all roots of the root function in the interval [-1,1]
        all_roots <- uniroot.all(solve_MFA,interval = c(-10, 10))
        
        # Determine the stability of each root by testing if the function is positive or negative prior to the root
        res <- data.frame(stability = sapply(all_roots, function(x) ifelse(solve_MFA(x - 0.001) > 0, "stable", "unstable")),
                          value = all_roots)
      }
    )
    result$value = as.numeric(as.character(result$value))
    return(result)
  }
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




