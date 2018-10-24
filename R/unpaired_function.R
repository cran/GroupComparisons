
#' @title Unpaired Parametric/Non-Parametric Group Comparisons
#' @description Receives two vectors, computes best function for unpaired group comparison (t-test, Mann-Whitney),
#' and reports the findings (mean/median, standard deviation, test statistic, p-value, effect size) in APA format
#' (Field, A. (2013). Discovering statistics using IBM SPSS statistics. New York, NY: SAGE.).
#' @param vec1 A vector of numbers
#' @param vec2 A vector of numbers
#' @return This function returns a sentence summarizing the findings and reporting them in APA format (effect size included)
#' @examples
#' dt <- mtcars
#' vector1 <- dt$mpg
#' vector2 <- dt$hp
#' Group_Test <- Group_Comparison_Unpaired(vector1, vector2)
#' Group_Test
#' @export


Group_Comparison_Unpaired = function(vec1, vec2) {

  # Save inputs as Vec1 and Vec2
  Vec1 = vec1
  Vec2 = vec2

  # Define means, medians, and standard deviations
  # Vector 1
  mean1 = mean(Vec1)
  sd1 = sd(Vec1)
  med1 = median(Vec1)
  # Vector 2
  mean2 = mean(Vec2)
  sd2 = sd(Vec2)
  med2 = median(Vec2)

  # Check both vectors for normality
  # Vector1
  Vec1_Shap = shapiro.test(Vec1)
  Vec1_Shap_p = Vec1_Shap$p.value
  # Vector2
  Vec2_Shap = shapiro.test(Vec2)
  Vec2_Shap_p = Vec2_Shap$p.value

  # Test vector 1 for normality
  if (Vec1_Shap_p < 0.05) {
    print("Vector1 Shapiro-Wilk test p-value is < 0.05: Data is not normally distributed.")
  } else {
    print("Vector1 Shapiro-Wilk test p-value is > or equal to 0.05: Data is normally distributed.")
  }

  # Test vector 2 for normality
  if (Vec2_Shap_p < 0.05) {
    print("Vector2 Shapiro-Wilk test p-value is < 0.05: Data is not normally distributed.")
  } else {
    print("Vector2 Shapiro-Wilk test p-value is > or equal to 0.05: Data is normally distributed.")
  }

  # Determine which test to use (parametric or non-parametric)
  if (Vec1_Shap_p >= 0.05 & Vec2_Shap_p >= 0.05) {
    print("Both vectors are normally distributed: A t-test will be used to compare these groups.")
  } else {
    print("One or more vectors are not normally distributed: Mann-Whitney U test will be used to compare these two groups.")
  }

  # Check the equality of variances
  #library(car)
  Equ_Var = leveneTest(Vec1, Vec2)
  Equ_Var_p = Equ_Var$`Pr(>F)`[1] # if p < 0.05, variances are not equal

  # Print the outcome
  if (Equ_Var_p < 0.05) {
    print("Levene's test of variances reveals unequal variances in Vector1 and Vector2.")
  } else {
    print("Levene's test of variances reveals equal variances in Vector1 and Vector2.")
  }

  # Compute the tests
  # If both groups are normal and variances are equal
  if (Vec1_Shap_p >= 0.05 & Vec2_Shap_p >= 0.05 & Equ_Var_p >= 0.05) {
    Group_Comparison = t.test(Vec1, Vec2, var.equal = TRUE)
    # If both groups are normal and variances are unequal
  } else if (Vec1_Shap_p >= 0.05 & Vec2_Shap_p >= 0.05 & Equ_Var_p < 0.05) {
    Group_Comparison = t.test(Vec1, Vec2, var.equal = FALSE)
    # If at least one group is not normal, use wilcox test
  } else {
    Group_Comparison = wilcox.test(Vec1, Vec2, exact = FALSE)
  }

  # Get group comparison p-values to make reporting easier
  if (Group_Comparison$p.value < 0.001) {
    Group_Comparison_p = '< 0.001'
  } else {
    Group_Comparison_p = Group_Comparison$p.value
  }

  # Compute effect size
  # for for t-test
  if (Vec1_Shap_p >= 0.05 & Vec2_Shap_p >= 0.05) {
    diff = abs(mean1-mean2) # mean absolute difference
    sdpool = sqrt(((sd1^2)+(sd2^2))/2) # pooled standard deviation
    ES = abs(diff/sdpool) # effect size
    # for wilcoxon test
  } else {
    Z = abs(qnorm(Group_Comparison$p.value))
    n = length(Vec1) + length(Vec2)
    r = Z/sqrt(n)
  }

  # Report the outcome of these tests
  # if Using t-test and there is a significant difference:
  if (Vec1_Shap_p >= 0.05 & Vec2_Shap_p >= 0.05 & Group_Comparison$p.value < 0.05) {
    return(paste("There was a significant difference in mean values between vector1 (M =", round(mean1, 2),
                 ", SD =", round(sd1,2),
                 ") and Vector2 (M =", round(mean2, 2),
                 ", SD = ", round(sd2, 2),
                 "), t(", round(Group_Comparison$parameter, 2), ") =", round(Group_Comparison$statistic,2),
                 ", p =", Group_Comparison_p,", d = ", round(ES, 2), "."))
    # if Using t_test and there IS NOT a significant difference:
  } else if (Vec1_Shap$p.value >= 0.05 & Vec2_Shap$p.value >= 0.05 & Group_Comparison$p.value >= 0.05) {
    return(paste("There was  not a significant difference in mean values between vector1 (M =", round(mean1, 2),
                 ", SD =", round(sd1,2),
                 " and Vector2 (M =", round(mean2, 2),
                 "), SD = ", round(sd2, 2),
                 ", t(", round(Group_Comparison$parameter, 2), ") =", round(Group_Comparison$statistic,2),
                 ", p =", Group_Comparison_p,", d = ", round(ES, 2), "."))
    # if using wilcoxon signed rank test and there is a significant difference:
  } else if (Vec1_Shap$p.value >= 0.05 | Vec2_Shap$p.value >= 0.05 & Group_Comparison$p.value < 0.05) {
    return(paste("There ws a significant difference in median values between Vector1 (Mdn =", round(med1, 2),
                 ") and Vector2 (Mdn =", round(med2, 2),
                 "), U =", round(Group_Comparison$statistic, 2),
                 ", p =", Group_Comparison_p, ", r =", round(r, 2), "."))
  }
  # End function
}


