# Repeatability for time-to-event data

It is common in animal behavior and cognition research to evaluate consistent individual differences (repeatability) in performance measured as the time to accomplish some event (latency to do a behavior or the number of trials to meet a pre-set criterion). Some individuals may never do the target behavior or meet the criterion and failing to account for the right-censored nature of these data is potentially problematic. Survival analyses explicitly account for censored data but there is currently no method for using this technique to quantify repeatability. Through a collaboration funded by the SQuID fellowship, we address this gap in our ability to statistically quantify individual differences in animal behavior.

### Contents
 - In our manuscript: "Repeatability and intra-class correlations from time-to-event data: towards a standardized approach" we describe this issue in more detail and present a solution.

 - We compiled [supplementary materials](https://kelseybmccune.github.io/Time-to-Event_Repeatability/Supplementary-materials.html) that consist of R code guiding the user through several worked examples from openly accessible data sets. In addition, we present results from a simulation study demonstrating the robustness of our proposed solution in relation to similar methods and varying data inputs.

 - We wrote a function in R that estimates the repeatability, as well as a confidence interval and p-value, from time-to-event data using survival analysis and the residual variance estimator we present in the manuscript.

