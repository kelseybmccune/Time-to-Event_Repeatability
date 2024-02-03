# Repeatability for time-to-event data

It is common in animal behavior to evaluate consistent individual differences (repeatability) in performance on time-to-event based tasks, like the latency to do a behavior or the number of trials to meet a pre-set criterion. Some individuals may never do the target behavior or meet the criterion and failing to account for the right-censored nature of these data is potentially problematic. Survival analyses explicitly account for censored data but there is currently no method for using this technique to quantify repeatability. Through a collaboration funded by the SQuID fellowship, we aim to fill this gap in our ability to statistically quantify individual differences in animal behavior.

## Goals
  - Determine methods for variance partitioning from mixed-effects survival models or similar analytical techniques that account for censored data so that repeatability can be quantified

  - Develop the capacity of the rptR package for integrating survival analysis. For this piece, we will also likely involve SQuID member Dr. Holger Schielzeth, as well as Dr. Martin Stoffel

  - Publish the methods for quantifying repeatability from time-to-event data and a comparison of the performance of this new method relative to results from the current approach of fitting Poisson or Gaussian models to data with ceiling values (rather than censored values)
