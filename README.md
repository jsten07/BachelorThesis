# BachelorThesis
## Comparing spatial determinants of urban growth from 1990 to 2014 between three European cities

### Abstract
Driven by a fast increasing population and urbanization the world’s cities expanded rapidly over the last decades and this trend will not stop in the following decades. To better understand and control the urban growth and thereby reduce the damages on the environment, it is helpful to know the driving factors of the expansion. This thesis aims to examine the spatial determinants of the urban growth of the three cities Dresden, Seville and Krakow for the time frame from 1990 to 2014 and compares them with regard to the questions where the urban growth occurred, how the determinants differ between the case study cities and why they differ. The investigated determinants are built-up density, population density, slope, distance to roads, distance to the city center, distance to train stations, distance to the airport, distance to the major river and the land use, divided in the 6 categories artificial, crop, pasture, open land, forest and water. They were examined by using a logistic regression model and the models were evaluated with the ROC value. The drivers showed different influences depending on the city. Each of the factors, except the distance to the river, had a significant effect on the probability of change to a built-up area for at least one of the cities. Built-up density was the only determinant which influenced the probability of change for each of the cities. Considering the differences between the cities, more location-based reasons were found than reasons depending on the cities’ histories. Especially areas that are closer to already built-up areas show a higher probability of change, possibly caused by lower costs for their development. The thesis thus showed distinct differences between the determinants of urban growth in various cities. To better understand these factors and their differences, future studies should consider further cities, determinants and time frames and compare the results with each other and the ones achieved in this research.

#### Workflow / use of code
* First data observation and preparation was done with QGIS. All other steps are performed using R.
* R functions for the data preparation are implemented in *preparation_funtions.R*.
* The preparation functions were used to prepare the data of the particular cities in *process_data.R* by creating a raster stack for every city.
* Different logisitc regression models can be calculated with the functions in *regression.R*.
  * Models including all variables calculated with the *glm* and *train* functions.
  * Models that only consist of variables that increase the models performances calculated with *ffs* and *bss* functions from the *CAST* package.
* Various sampling functions are implemented in the *sampling.R* section. For the thesis a stratified random sampling is used.
* The correlation between the variables was detected within the *interaction_variables.R*.
* Visualisation of the data and results was done in the *visualisation.R* part. 
  * Initial data
  * Regression coefficients in tables
  * Boxplots 
