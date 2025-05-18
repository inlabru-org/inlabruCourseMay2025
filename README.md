
<!-- README.md is generated from README.Rmd. Please edit that file -->

# inlabruCourseMay2025

<!-- badges: start -->

<!-- badges: end -->

This is the course materials page and package for the [Edinburgh
`inlabru`
course](https://centreforstatistics.maths.ed.ac.uk/events/upcoming-events/inlabru),
19 and 20 May 2025. Follow the menu links to information and tutorials,
as and when they become available.

## Event Overview

Bayesian Latent Gaussian Models (LGMs) are closely related to
Generalized Additive Models (GAMs), offering Bayesian estimation and
uncertainty quantification for spatial and spatio-temporal models. The
[INLA](https://www.r-inla.org/) and
[inlabru](https://inlabru-org.github.io/inlabru/) R packages combine
these Gaussian process models with numerical optimization and
integration techniques, in a fast and flexible analysis toolkit. The
taught part of the course will provide an overview of LGM theory and the
INLA/inlabru methods and software, while the hands-on sessions will make
sure the attendees will be ready to start doing spatial LGM modelling in
R as soon as the course is over.

## Topics covered

Basics of latent Gaussian process models in the Bayesian spatial
statistics context. The principles of the INLA method for fast Baysian
inference, and inlabru extensions for non-linear models. The inlabu
package principles and interface. Building spatial and spatio-temporal
models for point-referenced, spatially aggregated, and point pattern
observations. Computing and assessing posterior predictions and
visualisation. Diagnosing modelling problems.

## Instructor

Prof [Finn Lindgren](https:://www.maths.ed.ac.uk/~flindgre) is Chair of
Statistics in the School of Mathematics, University of Edinburgh.

His research focuses on spatial and spatio-temporal stochastic models,
environmetrics, and computational methods and software. Among many
others, he co-authored the influential paper “An Explicit Link Between
Gaussian Fields and Gaussian Markov Random Fields: The Stochastic
Partial Differential Equation Approach,” published in the Journal of the
Royal Statistical Society: Series B. Professor Lindgren has contributed
to the development of several R packages, including INLA for Bayesian
latent Gaussian models and inlabru, a user-friendly interface for INLA
with additional features

## Learning outcomes

- Understand the basic theory underpinning spatial latent Gaussian
  process models and Bayesian inference
- Use inlabru to fit various spatial models to data, including
  point-referenced, aggregated, and point pattern data.
- Be able to compute posterior predictions.
- Assess and compare models

## Daily timetable (tentative)

### Monday 19/5

- 09:30 - 10:30 Lecture (Spatial modelling with random fields)
- 10:30 - 11:00 coffee break
- 11:00 - 12:30 Lecture/Hands on session ([Introduction to
  `INLA`/`inlabru`/`fmesher`](articles/random_fields_1d.html))
- 12:30 - 13:30 Lunch
- 13:30 - 15:00 Lecture/Hands on session ([Spatial models for
  point-referenced data](articles/random_fields_2d.html))
- 15:00 - 15:30 break
- 15:30 - 17:00 Lecture/Hands on session ([Aggregated counts and
  non-linear predictors](articles/aggregated_counts_1d.html))
  ([Non-separable space-time](articles/nonseparable_spacetime.html))

### Tuesday 20/5

- 09:30 - 10:30 Lecture (Point process models)
- 10:30 - 11:00 coffee break
- 11:00 - 12:30 Lecture/Hands on session ([Spatial
  covariates](articles/lgcp_2d_covars.html))
  ([Space-time](articles/lgcp_2d_spatiotemporal.html))
- 12:30 - 13:30 Lunch
- 13:30 - 15:00 Lecture/Hands on session ([Multi-likelihood models;
  hurdle
  models](https://inlabru-org.github.io/inlabru/articles/zip_zap_models.html))
- 15:00 - 15:30 break
- 15:30 - 17:00 Lecture/Hands on session ([Predictive model
  assessment](https://inlabru-org.github.io/inlabru/articles/prediction_scores.html))
  ([Spatially varying
  coefficients](https://inlabru-org.github.io/inlabru/articles/zip_zap_models.html))

## Target audience

Anyone with some statistics training who is aware of spatial data and
the basics of additive models could benefit from attending. Fields where
this may be most popular are: ecology, geosciences, ecology,
epidemiology, public health, psychology, econometrics.

## Assumed knowledge

Attendees should be comfortable with using R. They should understand
linear and additive models, though this can be intuitive and doesn’t
have to be mathematically rigorous. They do not need to have used INLA
or inlabru before.
