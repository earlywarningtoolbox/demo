library(shiny)
library(ggplot2)

  shinyUI(pageWithSidebar(

headerPanel("Early Warning Signals Toolbox - Quick Detection Analysis"),
sidebarPanel(
selectInput(inputId = 'timeseries', "Choose a time-series:",
choices = c("simulated - overharvested resource", "real-world - climate data", "User data")),

selectInput(inputId = 'analysis',
label = 'Analysis',
choices = c("Indicator trend analysis", "Trend significance analysis","Potential analysis"),
selected = "Indicator trend analysis"),

# Display this only with analysis = "Indicator trend analysis"
conditionalPanel(condition = "input.analysis == 'Indicator trend analysis'",

sliderInput(inputId = "winsize",
label = "Sliding window size (fraction of time series):",
min = 5, max = 100, value = 50, step = 0.2),

selectInput(inputId = 'detrending',
label = 'Detrending',
choices = c("no", "gaussian", "linear", "first-diff"),
selected = "no"),

checkboxInput('logtransform', 'Logarithmize'),

checkboxInput('interpolate', 'Interpolate'),
                 
submitButton("Update View")
),

# Display this only with analysis = "Trend significance analysis"
conditionalPanel(condition = "input.analysis == 'Trend significance analysis'",

sliderInput(inputId = "winsize",
label = "Sliding window size (fraction of time series):",
min = 5, max = 100, value = 50, step = 0.2),

selectInput(inputId = 'detrending',
label = 'Detrending',
choices = c("no", "gaussian", "linear", "first-diff"),
selected = "no"),

checkboxInput('logtransform', 'Logarithmize'),

checkboxInput('interpolate', 'Interpolate'),

selectInput(inputId = "boots",
label = 'Number of surrogate sets',
choices = c(50, 100, 200, 1000),
selected = 50),
                 
# selectInput(inputId = 's_level',
# label = 'Level of significance',
# choices = c(0.1, 0.05, 0.01),
# selected = 0.05),

radioButtons(inputId = 's_level',
label = "Level of significance:",
choices = c(0.1, 0.05, 0.01),
selected = 0.05),

                 
helpText("Note: depending on the length of the time series
          the estimation of the trends from the surrogate
          datasets may be slow"),
                 
submitButton("Update View")
),

# Display this only with analysis = "Potential analysis"
conditionalPanel(condition = "input.analysis == 'Potential analysis'",

sliderInput(inputId = "detection.threshold",
label = "Threshold for local minima detection",
min = 0, max = 0.5, value = 0.002, step = 0.001),

sliderInput(inputId = 'grid.size',
label = 'Grid size',
min = 10, max = 200, value = 25, step = 5),

sliderInput(inputId = 'cutoff',
label = 'Cutoff for visualizing the potential landscape',
min = 0, max = 1, value = 0.5, step = 0.01),
                 
submitButton("Update View")                 
)

),

mainPanel(
plotOutput(outputId = "plot", height = "500px")
)
)
)






