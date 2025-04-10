function marginal_r_squared=Marginal_R2(tbl, equation)
% Calculate the total variance of FC
total_variance = var(tbl.FC);
% Fit a model with only fixed effects (ignoring random effects)
lme_fixed = fitlme(tbl, equation);
% Get the fitted values from the fixed-effects-only model
fitted_values = predict(lme_fixed);
% Calculate the variance explained by fixed effects
fixed_effects_variance = var(fitted_values);
% Calculate the residual variance (variance explained by random effects)
random_effects_variance = total_variance - fixed_effects_variance;

% Compute Marginal R² (fixed effects only)
marginal_r_squared = fixed_effects_variance / total_variance;