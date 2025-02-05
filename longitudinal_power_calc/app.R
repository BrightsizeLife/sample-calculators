# ============================================================================
# Title    : Mixed Model Power Explorer
# Purpose  : Demonstrate random effects with correlation, plus Type S & Type M 
#            with error bars for both Frequentist & Bayesian approaches.
# Author   : Derek
# Date     : YYYY-MM-DD
# ----------------------------------------------------------------------------
library(shiny)
library(ggplot2)
library(MASS)      # for mvrnorm to get correlated random effects
library(lme4)
library(rstanarm)
library(dplyr)

# We'll let rstanarm run multiple cores if available
options(mc.cores = parallel::detectCores())

# ----------------------------------------------------------------------------
# 1. Minimal Polis-inspired style (inline CSS)
# ----------------------------------------------------------------------------
app_css <- "
/* RESET & GLOBAL */
* {
  margin: 0;
  padding: 0;
  box-sizing: border-box;
}
:root {
  --color-bg: #ffffff;
  --color-text: #60665f;
  --color-link: #007bff;
  --color-heading: #333333;
  --font-main: 'Space Mono', monospace;
}
body {
  background: var(--color-bg);
  color: var(--color-text);
  font-family: var(--font-main);
  font-size: 16px;
  line-height: 1.5;
}
main, .container {
  max-width: 800px;
  margin: 0 auto;
  padding: 2rem;
}
h1 {
  font-size: 2.5rem;
  font-weight: 400;
  margin: 1rem 0 0.5rem;
  color: var(--color-heading);
}
h2, h3, h4 {
  font-weight: 700;
  color: var(--color-heading);
  margin: 1rem 0 0.5rem;
}
a {
  color: var(--color-link);
  text-decoration: none;
}
a:hover {
  text-decoration: underline;
}
.btn {
  font-family: var(--font-main);
  background: var(--color-link);
  color: #fff;
  border: none;
  padding: 0.5rem 1rem;
  cursor: pointer;
  text-decoration: none;
  transition: background 0.2s;
}
.btn:hover {
  background: #0056b3;
}
footer {
  margin-top: 2rem;
  padding-top: 1rem;
  border-top: 1px solid #ccc;
  text-align: center;
  font-size: 0.9rem;
  color: #888;
}
"

# ----------------------------------------------------------------------------
# 2. Simulation Function with Correlated Random Effects
# ----------------------------------------------------------------------------
simulate_data <- function(J      = 30,
                          K      = 4,
                          mu_a   = 4.8,
                          mu_b   = -5,
                          sigma_a= 1.3,
                          sigma_b= 0.7,
                          corr_ab= 0.0,   # correlation between intercept & slope
                          gamma1 = 2.0,
                          beta2  = 0.5) {
  # time points
  time_vals  <- rep(seq(0, 1, length.out = K), J)
  person_ids <- rep(seq_len(J), each = K)
  
  # assign half the subjects to treatment
  treatment_assignment <- sample(rep(0:1, length.out = J))
  treatment_vec <- treatment_assignment[person_ids]
  
  # correlated random intercept & slope
  #  Cov matrix: [[sigma_a^2, corr * sigma_a*sigma_b],
  #               [corr*sigma_a*sigma_b, sigma_b^2]]
  # means = (mu_a, mu_b)
  Sigma <- matrix(c(sigma_a^2, corr_ab*sigma_a*sigma_b,
                    corr_ab*sigma_a*sigma_b, sigma_b^2), 2,2)
  re_draws <- MASS::mvrnorm(J, mu = c(mu_a, mu_b), Sigma = Sigma)
  a_person <- re_draws[,1]  # intercept
  b_person <- re_draws[,2]  # slope
  
  # outcome
  y_vals <- a_person[person_ids] +
    b_person[person_ids] * time_vals +
    gamma1 * treatment_vec +
    beta2  * (time_vals * treatment_vec) +
    rnorm(J*K, 0, 1.0)
  
  data.frame(
    person    = factor(person_ids),
    time      = time_vals,
    treatment = factor(treatment_vec),
    y         = y_vals
  )
}

# ----------------------------------------------------------------------------
# 3. Parameter Mapping
# ----------------------------------------------------------------------------
# time => "time"
# treatment => "treatment1"
# time×treatment => "time:treatment1"
map_parameter_name <- function(choice) {
  if (choice == "Time effect") {
    return("time")
  } else if (choice == "Treatment effect") {
    return("treatment1")
  } else {
    return("time:treatment1")
  }
}
get_true_val <- function(param, mu_b, gamma1, beta2) {
  if (param == "time") {
    return(mu_b)
  } else if (param == "treatment1") {
    return(gamma1)
  } else {
    return(beta2)
  }
}

# ----------------------------------------------------------------------------
# 4. Frequentist: Power + Type M
# ----------------------------------------------------------------------------
# Type M = mean(|estimate - trueVal|)
freq_metrics_for_j <- function(J,
                               K        = 4,
                               N        = 200,
                               param    = "time:treatment1",
                               mu_a     = 4.8,
                               mu_b     = -5,
                               sigma_a  = 1.3,
                               sigma_b  = 0.7,
                               corr_ab  = 0,
                               gamma1   = 2.0,
                               beta2    = 0.5) {
  sig_indicator <- logical(N)
  abs_errors    <- numeric(N)
  sign_error    <- numeric(N)
  
  true_val <- get_true_val(param, mu_b, gamma1, beta2)
  
  for (i in seq_len(N)) {
    df_sim <- simulate_data(
      J       = J, K = K,
      mu_a    = mu_a, mu_b = mu_b,
      sigma_a = sigma_a, sigma_b = sigma_b,
      corr_ab = corr_ab,
      gamma1  = gamma1, beta2 = beta2
    )
    
    fit   <- lmer(y ~ time * treatment + (1 | person), data = df_sim)
    coefs <- summary(fit)$coefficients
    row_idx <- which(rownames(coefs) == param)
    
    if (length(row_idx) == 0) {
      sig_indicator[i] <- FALSE
      abs_errors[i]    <- NA
      sign_error[i]    <- NA
      next
    }
    
    est <- coefs[row_idx, "Estimate"]
    se  <- coefs[row_idx, "Std. Error"]
    
    z_val <- est / se
    p_val <- 2 * (1 - pnorm(abs(z_val)))
    sig_indicator[i] <- (p_val < 0.05)
    
    # sign error if estimate * true_val < 0 => opposite sign
    sign_error[i] <- if ((est * true_val) < 0) 1 else 0
    
    # Type M = abs(est - true_val)
    abs_errors[i] <- abs(est - true_val)
  }
  
  power_val <- mean(sig_indicator, na.rm = TRUE)
  
  # Binomial approx for power
  p <- power_val
  se_p <- sqrt(p * (1 - p) / N)
  power_low  <- max(0, p - 1.96 * se_p)
  power_high <- min(1, p + 1.96 * se_p)
  
  mean_abs  <- mean(abs_errors, na.rm = TRUE)
  low_q     <- quantile(abs_errors, 0.025, na.rm = TRUE)
  high_q    <- quantile(abs_errors, 0.975, na.rm = TRUE)
  
  # Type S
  type_s_val <- mean(sign_error, na.rm = TRUE)
  n_sign  <- sum(sign_error, na.rm = TRUE)
  n_total <- sum(!is.na(sign_error))
  type_s_low  <- qbeta(0.025, n_sign + 1, n_total - n_sign + 1)
  type_s_high <- qbeta(0.975, n_sign + 1, n_total - n_sign + 1)
  
  c(
    Power      = power_val,
    Power_low  = power_low,
    Power_high = power_high,
    
    TypeM      = mean_abs,
    M_low      = low_q,
    M_high     = high_q,
    
    TypeS      = type_s_val,
    TypeS_low  = type_s_low,
    TypeS_high = type_s_high
  )
}



# ----------------------------------------------------------------------------
# 5. Bayesian: Type S + Type M
# ----------------------------------------------------------------------------
bayes_metrics_for_j <- function(J,
                                K       = 4,
                                N       = 50,
                                param   = "time:treatment1",
                                mu_a    = 4.8,
                                mu_b    = -5,
                                sigma_a = 1.3,
                                sigma_b = 0.7,
                                corr_ab = 0,
                                gamma1  = 2.0,
                                beta2   = 0.5) {
  
  type_s_vec <- numeric(N)
  type_m_vec <- numeric(N)
  
  true_val <- get_true_val(param, mu_b, gamma1, beta2)
  
  for (i in seq_len(N)) {
    df_sim <- simulate_data(J = J, K = K,
                            mu_a   = mu_a, mu_b   = mu_b,
                            sigma_a= sigma_a, sigma_b= sigma_b,
                            corr_ab= corr_ab,
                            gamma1 = gamma1, beta2  = beta2)
    
    fit_bayes <- suppressWarnings(
      stan_glmer(
        y ~ time * treatment + (1 | person),
        data    = df_sim,
        family  = gaussian(),
        chains  = 2,
        cores   = 2,       # parallel each fit
        iter    = 1000,
        refresh = 0
      )
    )
    post_draws <- as.matrix(fit_bayes)
    col_idx    <- which(colnames(post_draws) == param)
    
    if (length(col_idx) == 0) {
      type_s_vec[i] <- NA
      type_m_vec[i] <- NA
      next
    }
    draws_param <- post_draws[, col_idx]
    
    # Type S: fraction on the opposite side of 0 from true effect
    if (true_val > 0) {
      type_s_vec[i] <- mean(draws_param < 0)
    } else {
      type_s_vec[i] <- mean(draws_param > 0)
    }
    
    # Type M = average abs distance from true
    type_m_vec[i] <- mean(abs(draws_param - true_val))
  }
  
  tS_mean <- mean(type_s_vec, na.rm = TRUE)
  tS_low  <- quantile(type_s_vec, 0.025, na.rm = TRUE)
  tS_high <- quantile(type_s_vec, 0.975, na.rm = TRUE)
  
  tM_mean <- mean(type_m_vec, na.rm = TRUE)
  tM_low  <- quantile(type_m_vec, 0.025, na.rm = TRUE)
  tM_high <- quantile(type_m_vec, 0.975, na.rm = TRUE)
  
  c(TypeS_mean = tS_mean,
    TypeS_low  = tS_low,
    TypeS_high = tS_high,
    TypeM_mean = tM_mean,
    TypeM_low  = tM_low,
    TypeM_high = tM_high)
}

# ----------------------------------------------------------------------------
# 6. Shiny UI
# ----------------------------------------------------------------------------
ui <- fluidPage(
  tags$head(
    tags$style(HTML(app_css)),
    tags$link(
      rel  = "stylesheet",
      href = "https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&display=swap"
    )
  ),
  
  titlePanel("Mixed Model Power Explorer"),
  
  tabsetPanel(
    
    # ------------------------------------------------------------------------
    # A) Effect Size Explorer
    # ------------------------------------------------------------------------
    tabPanel(
      "Effect Size Explorer",
      fluidRow(
        column(
          width = 4,
          h4("Random Effects"),
          sliderInput("sigma_a", "SD of Intercepts (sigma_a):",
                      min = 0, max = 5, value = 1.3, step = 0.1),
          sliderInput("sigma_b", "SD of Slopes (sigma_b):",
                      min = 0, max = 5, value = 0.7, step = 0.1),
          sliderInput("corr_ab", "Correlation (Intercept, Slope):",
                      min = -1, max = 1, value = 0, step = 0.1),
          p("Positive correlation => subjects with higher intercepts also have steeper slopes. ",
            "Negative => higher intercept, shallower slope. Zero => independent.")
        ),
        column(
          width = 4,
          h4("Fixed Effects"),
          sliderInput("mu_a", "Baseline (mu_a):",
                      min = 0, max = 10, value = 4.8, step = 0.1),
          sliderInput("mu_b", "Time effect (Control) (mu_b):",
                      min = -10, max = 10, value = -5, step = 0.1),
          p("If mu_b < 0, control decreases from time=0 to time=1. ",
            "If > 0, it increases.")
        ),
        column(
          width = 4,
          h4("Treatment & Interaction"),
          sliderInput("gamma1", "Treatment @ Baseline (gamma1):",
                      min = -10, max = 10, value = 2, step = 0.1),
          sliderInput("beta2",  "Time × Treatment (beta2):",
                      min = -10, max = 10, value = 0.5, step = 0.1),
          p("Positive beta2 => slope in treatment is mu_b + beta2, i.e. more steep if beta2>0.")
        )
      ),
      fluidRow(
        column(
          width = 5,
          h4("Random Effects Distribution"),
          plotOutput("random_effects_plot", height = "250px")
        ),
        column(
          width = 5,
          h4("Trajectory Preview (No Noise)"),
          plotOutput("trajectory_plot", height = "250px")
        )
      ),
      fluidRow(
        column(
          width = 10,
          offset = 1,
          h4("Current Parameters"),
          verbatimTextOutput("model_params")
        )
      )
    ),
    
    # ------------------------------------------------------------------------
    # B) Frequentist Tab
    # ------------------------------------------------------------------------
    tabPanel(
      "Frequentist",
      sidebarLayout(
        sidebarPanel(
          numericInput("f_n1",    "Min sample size:",   value = 20, min = 2),
          numericInput("f_n2",    "Max sample size:",   value = 100, min = 2),
          numericInput("f_step",  "Step size:",         value = 10, min = 1),
          numericInput("f_nsims", "Num sims per size:", value = 200, min = 10),
          
          selectInput("f_param",
                      "Parameter to Explore:",
                      choices = c("Time effect", "Treatment effect", "Time × Treatment"),
                      selected = "Time × Treatment"),
          
          actionButton("run_freq", "Run Frequentist", class = "btn")
        ),
        mainPanel(
          plotOutput("freq_plot_power"),
          plotOutput("freq_plot_s"),
          plotOutput("freq_plot_m")
        )
      )
    ),
    
    # ------------------------------------------------------------------------
    # C) Bayesian Tab
    # ------------------------------------------------------------------------
    tabPanel(
      "Bayesian",
      sidebarLayout(
        sidebarPanel(
          tags$strong("Warning: Bayesian calculations can be slow."),
          numericInput("b_n1",    "Min sample size:",   value = 20, min = 2),
          numericInput("b_n2",    "Max sample size:",   value = 60, min = 2),
          numericInput("b_step",  "Step size:",         value = 10, min = 1),
          numericInput("b_nsims", "Num sims per size:", value = 50, min = 10),
          
          selectInput("b_param",
                      "Parameter to Explore:",
                      choices = c("Time effect", "Treatment effect", "Time × Treatment"),
                      selected = "Time × Treatment"),
          
          actionButton("run_bayes", "Run Bayesian", class = "btn btn-danger")
        ),
        mainPanel(
          plotOutput("bayes_plot_s"),
          plotOutput("bayes_plot_m")
        )
      )
    ),
    
    # ------------------------------------------------------------------------
    # D) About
    # ------------------------------------------------------------------------
    tabPanel(
      "About",
      fluidRow(
        column(
          width = 10, offset = 1,
          h3("Understanding This Dashboard"),
          p("This dashboard shows how sample size (J), random effects (with possible correlation), ",
            "and fixed effects (time, treatment, time×treatment) influence your ability to estimate ",
            "and detect parameters in a linear mixed model with repeated measurements over time."),
          tags$ul(
            tags$li("Frequentist Power: Probability of p<0.05 for the chosen parameter."),
            tags$li("Frequentist Type M: Mean absolute error |estimate - true|, with error bars showing 2.5% & 97.5% quantiles."),
            tags$li("Bayesian Type S: Probability that the posterior sign is wrong (Gelman & Carlin, 2014)."),
            tags$li("Bayesian Type M: Mean of |posterior draw - true|, then averaged across simulations, with 2.5%-97.5% error bars.")
          ),
          p("Adjust 'sigma_a' and 'sigma_b' to control how much random intercepts and slopes vary across subjects, ",
            "and 'corr_ab' if you expect correlation between individuals who start higher/lower and how they evolve over time."),
          p("For 'time effect' (mu_b), 'treatment effect' (gamma1), or 'interaction' (beta2), you can specify any combination ",
            "and see how detection or estimation error changes."),
          p("Type M was changed to mean absolute error from a simpler notion of confidence/credible interval width. ",
            "If the random‐effects correlation is nonzero, you'll see a different distribution in the 'Random Effects Distribution' panel."),
          h3("Sample Size Clarification"),
          p("The 'sample size' J in this app is the total number of subjects.
          By default, we randomly assign half to treatment and half to control.
          For example, if J=50, it means ~25 control and ~25 treatment. 
          If you want J to represent subjects per group, you could double it."),
          h3("Random Effects & Correlation Impact"),
          p("When 'sigma_a' (random intercept SD) or 'sigma_b' (random slope SD) 
          is large, there's more subject-to-subject variability, making it 
          harder to detect the fixed (average) effect. 
          The correlation 'corr_ab' indicates whether individuals 
          with higher intercepts also tend to have steeper or shallower slopes.
          All of this extra variation can lower power and raise Type M 
          (mean absolute error), because the model must disentangle larger random fluctuations."),
          h4("References & Credits"),
          tags$ul(
            tags$li("Gelman, A., & Carlin, J. (2014). Beyond power calculations: Assessing Type S and Type M errors. ",
                    "Perspectives on Psychological Science, 9(6), 641-651."),
            tags$li("Thanks to the R community for lme4, rstanarm, and shiny!")
          )
        )
      )
    )
  ),
  
  tags$footer(
    HTML(
      "<p>Check out the source code on <a href='https://github.com/BrightsizeLife/sample-calculators' target='_blank'>GitHub</a>. 
       Contributions and feedback are welcome!</p>"
    )
  )
)

# ----------------------------------------------------------------------------
# 7. Server
# ----------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # ---- Show current random/fixed effect parameters
  output$model_params <- renderPrint({
    cat("sigma_a =", input$sigma_a, "(SD of intercepts)\n")
    cat("sigma_b =", input$sigma_b, "(SD of slopes)\n")
    cat("corr_ab =", input$corr_ab, "(correlation of intercept & slope)\n")
    cat("mu_a    =", input$mu_a,    "(baseline intercept)\n")
    cat("mu_b    =", input$mu_b,    "(control slope)\n")
    cat("gamma1  =", input$gamma1,  "(treatment effect @ baseline)\n")
    cat("beta2   =", input$beta2,   "(time × treatment slope)\n")
  })
  
  # ---- Plot random effects distribution
  # We'll draw e.g. 300 samples from the correlated distribution, just to illustrate.
  output$random_effects_plot <- renderPlot({
    set.seed(123)  # for reproducible example
    Sigma <- matrix(c(input$sigma_a^2, input$corr_ab * input$sigma_a * input$sigma_b,
                      input$corr_ab * input$sigma_a * input$sigma_b, input$sigma_b^2), 2,2)
    re_draws <- mvrnorm(300, mu = c(input$mu_a, input$mu_b), Sigma = Sigma)
    df_re    <- data.frame(Intercept = re_draws[,1], Slope = re_draws[,2])
    
    ggplot(df_re, aes(x = Intercept, y = Slope)) +
      geom_point(alpha = 0.3, color = "#007bff") +
      coord_fixed() +
      labs(
        title = "Random Intercept & Slope Distribution",
        x     = "Random Intercept (a_person)",
        y     = "Random Slope (b_person)"
      ) +
      theme_minimal(base_size = 14)
  })
  
  
  output$trajectory_plot <- renderPlot({
    # Create a sequence of times from 0 to 1
    times <- seq(0, 1, length.out = 20)
    
    # Control’s mean: mu_a + mu_b * t
    control_y <- input$mu_a + input$mu_b * times
    
    # Treatment’s mean: (mu_a + gamma1) + (mu_b + beta2)*t
    treat_y   <- (input$mu_a + input$gamma1) + (input$mu_b + input$beta2) * times
    
    # Combine into a small data frame for plotting
    df_plot <- data.frame(
      time  = rep(times, 2),
      group = rep(c("Control", "Treatment"), each = length(times)),
      y     = c(control_y, treat_y)
    )
    
    ggplot(df_plot, aes(x = time, y = y, color = group)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c("Control"="#999999", "Treatment"="#E69F00")) +
      labs(
        title = "Mean Trajectories (Ignoring Random Effects)",
        x = "Time (0 to 1)",
        y = "Mean Outcome"
      ) +
      theme_minimal(base_size = 14)
  })
  
  # ---------------- FREQUENTIST REACTIVE ----------------
  freq_results <- eventReactive(input$run_freq, {
    param_name <- map_parameter_name(input$f_param)
    J_vals <- seq(input$f_n1, input$f_n2, by = input$f_step)
    
    withProgress(message = "Running Frequentist Sims...", value = 0, {
      do.call(rbind, lapply(J_vals, function(Ji) {
        incProgress(1 / length(J_vals), detail = paste("J =", Ji))
        
        out <- freq_metrics_for_j(
          J       = Ji,
          N       = input$f_nsims,
          param   = param_name,
          mu_a    = input$mu_a,
          mu_b    = input$mu_b,
          sigma_a = input$sigma_a,
          sigma_b = input$sigma_b,
          corr_ab = input$corr_ab,
          gamma1  = input$gamma1,
          beta2   = input$beta2
        )
        
        data.frame(
          J         = Ji,
          Power     = out["Power"],
          Power_low = out["Power_low"],
          Power_high= out["Power_high"],
          
          TypeM     = out["TypeM"],
          M_low     = out["M_low"],
          M_high    = out["M_high"],
          
          TypeS     = out["TypeS"],
          TypeS_low = out["TypeS_low"],
          TypeS_high= out["TypeS_high"]
        )
      }))
    })
  })
  
  
  output$freq_plot_power <- renderPlot({
    df <- freq_results()
    req(df)
    
    ggplot(df, aes(x = J, y = Power)) +
      geom_line(color = "#2C3E50", size = 1) +
      geom_point(color = "#2C3E50", size = 2) +
      geom_errorbar(
        aes(ymin = Power_low, ymax = Power_high),
        width = 0.5, color = "#2C3E50"
      ) +
      ylim(0,1) +
      labs(
        title = "Frequentist Power vs. Sample Size",
        x     = "Sample Size (J)",
        y     = "Power (p<0.05)"
      ) +
      theme_minimal(base_size = 14)
  })
  
  
  
  output$freq_plot_s <- renderPlot({
    df <- freq_results()
    req(df)
    
    ggplot(df, aes(x = J, y = TypeS)) +
      geom_line(color = "#9b59b6", size = 1) +
      geom_point(color = "#9b59b6", size = 2) +
      geom_errorbar(
        aes(ymin = TypeS_low, ymax = TypeS_high),
        width = 0.5, color = "#9b59b6"
      ) +
      ylim(0,1) +
      labs(
        title = "Frequentist Type S (Wrong Sign) vs. Sample Size",
        x     = "Sample Size (J)",
        y     = "Type S"
      ) +
      theme_minimal(base_size = 14)
  })
  
  
  
  output$freq_plot_m <- renderPlot({
    df <- freq_results()
    req(df)
    
    ggplot(df, aes(x = J, y = TypeM)) +
      geom_line(color = "#1abc9c", size = 1) +
      geom_point(color = "#1abc9c", size = 2) +
      geom_errorbar(
        aes(ymin = M_low, ymax = M_high),
        width = 0.5, color = "#1abc9c"
      ) +
      labs(
        title = "Frequentist Type M (Mean Abs Error) vs. Sample Size",
        x     = "Sample Size (J)",
        y     = "Type M Error"
      ) +
      theme_minimal(base_size = 14)
  })
  
  
  # ---------------- BAYESIAN REACTIVE ----------------
  bayes_results <- eventReactive(input$run_bayes, {
    param_name <- map_parameter_name(input$b_param)
    J_vals     = seq(input$b_n1, input$b_n2, by = input$b_step)
    
    withProgress(message = "Running Bayesian Sims...", value = 0, {
      do.call(rbind, lapply(J_vals, function(Ji) {
        incProgress(1 / length(J_vals), detail = paste("J =", Ji))
        
        out <- bayes_metrics_for_j(
          J       = Ji,
          N       = input$b_nsims,
          param   = param_name,
          mu_a    = input$mu_a,
          mu_b    = input$mu_b,
          sigma_a = input$sigma_a,
          sigma_b = input$sigma_b,
          corr_ab = input$corr_ab,
          gamma1  = input$gamma1,
          beta2   = input$beta2
        )
        
        data.frame(
          J           = Ji,
          TypeS_mean  = out["TypeS_mean"],
          TypeS_low   = out["TypeS_low"],
          TypeS_high  = out["TypeS_high"],
          TypeM_mean  = out["TypeM_mean"],
          TypeM_low   = out["TypeM_low"],
          TypeM_high  = out["TypeM_high"]
        )
      }))
    })
  })
  

  
  
  output$bayes_plot_s <- renderPlot({
    df <- bayes_results()
    req(df)
    ggplot(df, aes(x = J, y = TypeS_mean)) +
      geom_line(color = "firebrick", size = 1) +
      geom_point(color = "firebrick", size = 2) +
      geom_errorbar(
        aes(ymin = TypeS_low, ymax = TypeS_high),
        width = 0.5, color = "firebrick"
      ) +
      ylim(0,1) +
      labs(
        title = "Bayesian Type S (Wrong Sign) vs. Sample Size",
        x     = "Sample Size (J)",
        y     = "Type S"
      ) +
      theme_minimal(base_size = 14)
  })
  
  output$bayes_plot_m <- renderPlot({
    df <- bayes_results()
    req(df)
    ggplot(df, aes(x = J, y = TypeM_mean)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(color = "darkgreen", size = 2) +
      geom_errorbar(
        aes(ymin = TypeM_low, ymax = TypeM_high),
        width = 0.5, color = "darkgreen"
      ) +
      labs(
        title = "Bayesian Type M (Abs Distance) vs. Sample Size",
        x     = "Sample Size (J)",
        y     = "Type M Error"
      ) +
      theme_minimal(base_size = 14)
  })
}

# ----------------------------------------------------------------------------
# 8. shinyApp
# ----------------------------------------------------------------------------
shinyApp(ui, server)




