# Rocker R studio image with R 4.4 as the base
FROM rocker/rstudio:4.4

# Install the necessary packages for the analysis
RUN R -e "install.packages(c('here','tidyr','ComplexHeatmap'))"

# Create the directory for RStudio preferences
RUN mkdir -p /home/rstudio/.config/rstudio

# Set Rstudio preferences using a .json file
RUN echo '{"editor_theme": "Material", "rainbow_parentheses": true, "highlight_r_function_calls": true}' \
  > /home/rstudio/.config/rstudio/rstudio-prefs.json




