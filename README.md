# Lion Movement and Habitat Selection During a Flooding Event

## Overview

This project investigates how lions (*Panthera leo*) responded to a major flooding event in Kruger National Park, South Africa, in January 2026.

GPS tracking data from free-ranging lions were combined with environmental and spatial datasets to examine changes in movement distances, habitat selection, diel activity patterns and landscape use before, during and after the flood.

The analyses use linear mixed-effects models, generalised additive mixed models and integrated step-selection functions to investigate how lions responded to rapidly changing environmental conditions.

## Study area and data

The study was conducted in the central region of Kruger National Park, South Africa.

GPS tracking data were collected from six lions and analysed across three consecutive 10-day periods:

- **Before:** 31 December 2025 – 9 January 2026
- **During:** 10–19 January 2026
- **After:** 20–29 January 2026

The analysis combines GPS telemetry with weather station observations, elevation data, river and road networks, and vegetation cover derived from Dynamic World satellite products.

All spatial analyses use projected coordinates to calculate movement distances and environmental characteristics.

## Analyses

### Movement distances

Hourly and daily movement distances were calculated from successive GPS locations.

Linear mixed-effects models were used to investigate differences in daily movement distances between flood periods.

Lion identity was included as a random effect, with temporal autocorrelation and differences in residual variance accounted for in the model structure.

### Elevation and river use

Elevation and distance to the nearest river were extracted for GPS locations.

Daily mean elevation and river distances were compared between flood periods using mixed-effects models.

These analyses examined whether lions changed their use of different landscape features during and after flooding.

### Integrated step-selection analysis

Integrated step-selection functions (iSSFs) were used to investigate habitat selection while accounting for movement characteristics.

Observed movement steps were compared with randomly generated available steps.

Environmental predictors included:

- Elevation
- Distance to rivers
- Distance to roads
- Woody vegetation cover

Interactions with flood period were used to investigate whether habitat selection changed before, during and after flooding.

Step length and turning angle terms were included to account for movement characteristics.

### Diel movement patterns

Generalised additive mixed models (GAMMs) were used to investigate changes in hourly movement patterns.

Separate cyclic smooth functions were fitted for each flood period, allowing movement patterns to vary across the 24-hour cycle.

Lion identity and temporal autocorrelation were incorporated into the model.

Morning and evening movement peaks were identified from model predictions, with uncertainty estimated through simulations of model coefficients.

### Rainfall and movement

Hourly rainfall measurements were examined in relation to hourly movement distances.

Mixed-effects models were used to investigate whether the relationship between rainfall and movement differed between flood periods.

Temporal autocorrelation and diel variation were accounted for in the analysis.

### Road use

Distances between GPS locations and the nearest road were calculated using a spatial road network.

Daily road distances were compared between flood periods.

Step-selection models were also fitted to examine whether selection relative to roads changed during flooding.

### Woody vegetation selection

Woody vegetation cover was calculated within a 50 m radius using satellite-derived land-cover data.

Observed movement steps were compared with available steps to estimate selection for woody vegetation during each flood period.

Individual-lion models were fitted, and selection coefficients were summarised across animals to investigate population-level patterns.

### Spatial visualisation

Spatial datasets were processed and mapped to visualise environmental conditions across the study area.

These included vegetation probability layers, landscape features and GPS-derived spatial summaries.

## R packages

The analysis uses several R packages, including:

- `dplyr` and `tidyr` for data preparation
- `lubridate` for date and time processing
- `sf` and `terra` for spatial analysis
- `amt` for movement and step-selection analysis
- `survival` for conditional logistic regression
- `nlme`, `lme4` and `lmerTest` for mixed-effects modelling
- `mgcv` for generalised additive mixed models
- `emmeans` for model contrasts
- `ggplot2` for visualisation
- `MASS` for simulation-based uncertainty estimation

Additional packages used in individual sections are specified in the script.

## Running the analysis

The analysis is contained in:

`Lion_Flood_Master_Analysis_Clean.R`

The script requires the original GPS, weather, raster and spatial data files.

Input files include:

- Combined GPS collar and environmental data
- Weather station observations
- Elevation raster
- River shapefile
- Road network GeoPackage
- Dynamic World woody vegetation and probability rasters

Place the input files in the working directory and run the script in R.

Optional analyses can be enabled or disabled using the settings at the beginning of the script.

## Outputs

The script produces model summaries, statistical comparisons, spatial datasets and figures.

Outputs include:

- Daily movement distance comparisons
- Elevation and river-use summaries
- Habitat selection coefficients
- Road-use and road-selection analyses
- Diel movement predictions and confidence intervals
- Morning and evening movement peak estimates
- Rainfall and movement relationships
- Woody vegetation selection summaries
- Vegetation and spatial maps

Selected fitted models are saved as RDS files, with numerical results exported as CSV files and figures saved in vector formats.

## Ecological relevance

Extreme rainfall and flooding can rapidly alter habitat accessibility, landscape connectivity and resource availability.

Understanding how large carnivores respond to these events can provide insight into behavioural flexibility and the effects of short-term environmental disturbances on wildlife movement.

This project demonstrates how high-resolution telemetry and environmental data can be combined to investigate animal responses to extreme weather events.

The methods are applicable to other ecological questions involving movement behaviour, habitat selection, environmental disturbance and landscape use.

## Limitations

The analysis is based on a relatively small number of GPS-tracked individuals and a single flooding event.

Observed responses may reflect individual differences, local environmental conditions and other factors that were not explicitly included in the models.

Habitat selection estimates are conditional on the available movement steps generated by the step-selection procedure.

The analyses identify associations between flood period, environmental conditions and movement behaviour rather than establishing causal relationships.

## Data availability

The original GPS telemetry and associated spatial datasets are not included in this repository.

The analysis script is provided to demonstrate the statistical and spatial methods used in the study.

Reproducing the analysis requires access to the original datasets and supporting spatial files.
