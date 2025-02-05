Mixed Model Power Explorer 🚀
An interactive Shiny dashboard for simulating and analyzing mixed‐effects models.
You can vary:

Sample size
Random‐effects standard deviations (and correlation!)
Fixed‐effects (Time, Treatment, and Interaction)
Frequentist vs. Bayesian inference
And get Power, Type S, Type M (defined as mean absolute error), and more, complete with error bars and visualizations.

Features ✨
Effect Size Explorer

Adjust intercept (mu_a), time slope (mu_b), treatment effect at baseline (gamma1), and time×treatment interaction (beta2).
Control random‐effects SDs (sigma_a, sigma_b) and the correlation between intercept & slope.
Visualize the random‐effects distribution in a scatter plot.
Frequentist Tab

Power: Probability that 
𝑝
<
0.05
p<0.05 for the chosen parameter (Time, Treatment, or Interaction).
Type M: Mean absolute estimation error (
∣
𝜃
^
−
𝜃
true
∣
∣ 
θ
^
 −θ 
true
​
 ∣) across simulations.
Error bars (2.5% & 97.5% quantiles) show variability across simulations.
Bayesian Tab

Type S (Wrong Sign): Fraction of posterior draws that cross zero in the wrong direction (Gelman & Carlin, 2014).
Type M: Average absolute distance of the posterior draws from the true parameter.
Again, error bars display the 2.5% & 97.5% quantiles across repeated simulations.
About

Briefly explains how to interpret the metrics and parameters.
Quick Start 🏎
Clone or download this repo:
bash
Copy
Edit
git clone https://github.com/BrightsizeLife/sample-calculators.git
Install the required R packages in R or RStudio:
r
Copy
Edit
install.packages(c("shiny", "ggplot2", "lme4", "rstanarm", "MASS"))
Open the app.R file in RStudio (or your preferred environment).
Run the Shiny app:
r
Copy
Edit
shiny::runApp("path/to/app.R")
Explore in your browser: adjust parameters, run Frequentist or Bayesian simulations, and see the results!
Interpreting the Parameters 🧭
Random Effects

sigma_a: Variation in random intercepts across subjects.
sigma_b: Variation in random slopes across subjects.
corr_ab: Correlation between intercept & slope. Positive means those with higher intercepts also tend to have steeper slopes (or less negative).
Fixed Effects

mu_a: Baseline (control) intercept.
mu_b: “Time effect” in the control group’s slope from time=0 to time=1. Negative => it decreases over time.
gamma1: Treatment effect at baseline. If positive, the treatment group starts higher than control at time=0.
beta2: Additional slope for treatment. If control slope = mu_b, treatment slope = mu_b + beta2.
Parameter to Analyze (Time, Treatment, Time×Treatment)

You can see how the chosen parameter’s power (Frequentist) or Type S/M (Bayesian) changes as you vary sample size, effect sizes, and random effects.
Metrics Explained 📊
Power (Frequentist): Probability that 
𝑝
<
0.05
p<0.05.
Type M (Frequentist): Mean of 
∣
𝜃
^
−
𝜃
true
∣
∣ 
θ
^
 −θ 
true
​
 ∣ across simulations.
Type S (Bayesian): Probability that the posterior draws cross zero incorrectly, i.e., the sign is wrong relative to zero (Gelman & Carlin, 2014).
Type M (Bayesian): Average 
∣
𝜃
draw
−
𝜃
true
∣
∣θ 
draw
​
 −θ 
true
​
 ∣ in the posterior, aggregated over simulations.
We draw error bars from the distribution of these metrics over repeated simulations (2.5% & 97.5% quantiles).

Helpful Hints 💡
Parallel Processing: We set cores=2 in stan_glmer. Increase or decrease as your CPU allows.
Large N or wide J ranges can be slow. Start small for quick tests.
Beta2 Interpretation:
Additive slope difference. If control slope = mu_b, treatment slope = mu_b + beta2.
If you want, e.g., “50% steeper,” that implies (mu_b + beta2)/mu_b ≈ 1.5.
References 📚
Gelman & Carlin (2014). Beyond power calculations: Assessing Type S (sign) and Type M (magnitude) errors. Perspectives on Psychological Science, 9(6), 641-651.
Cohen (1988). Statistical Power Analysis for the Behavioral Sciences. Routledge.
Stan Development Team: rstanarm: Bayesian applied regression modeling via Stan. https://mc-stan.org/rstanarm/
Contributing 🤝
We welcome issues, feature requests, and pull requests!
Please visit our GitHub repo to collaborate.

Enjoy exploring your mixed‐model simulations—and have fun running your own power and Type S/Type M analyses! 🏆
