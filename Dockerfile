# docker image (tidyverse)
FROM rocker/verse:4.4.0

# Update and install standard utilities and Circos dependencies
RUN apt update && apt install -y \
    man \
    wget \
    python3 \
    python3-pip \
    libgd-dev \
    libxml2-dev \
    libcairo2-dev \
    libpango1.0-dev \
    perl \
    cpanminus \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install Conda and set up environment
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# add conda to the path
ENV PATH="/opt/conda/bin:$PATH"

# update conda and create the new conda environment (with -y for non-interactive)
RUN conda update -y conda 
RUN conda create -y -n run-env python=3.11

# Correctly initialize the run-env environment in .bashrc
RUN echo 'source /opt/conda/etc/profile.d/conda.sh && conda activate run-env' >> ~/.bashrc

# add new conda environment to the path
ENV PATH="/opt/conda/envs/run-env/bin:$PATH"

# Activate the newly created conda run environment and install packages using pip3
RUN conda run -n run-env pip install jupyter pysam

# Install Circos and its Perl dependencies
RUN cpanm Clone Config::General Font::TTF::Font GD Math::Bezier Math::Round \
    Math::VecStat Readonly Regexp::Common SVG Set::IntSpan Statistics::Basic \
    && wget --no-check-certificate http://circos.ca/distribution/circos-0.69-9.tgz \
    && tar xvfz circos-0.69-9.tgz \
    && mv circos-0.69-9 /opt/circos \
    && rm circos-0.69-9.tgz
ENV PATH="/opt/circos/bin:${PATH}"

# Install the necessary R packages
RUN R -e "install.packages(c('here','NatParksPalettes','drc','UpSetR','circlize'), repos='https://cloud.r-project.org/')"
RUN R -e "install.packages(c('ggridges','ggforce','cowplot','colorRamp2','gt','ggsignif'))"
RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_github('jokergoo/ComplexHeatmap')"

# Set password for rstudio user
RUN echo "rstudio:prame" | chpasswd

# Create the directory for RStudio preferences
RUN mkdir -p /home/rstudio/.config/rstudio

# Set Rstudio preferences using a .json file
RUN echo '{"editor_theme": "Material", "rainbow_parentheses": true, "highlight_r_function_calls": true}' \
  > /home/rstudio/.config/rstudio/rstudio-prefs.json

# Create a script to start both RStudio Server and Jupyter Notebook
RUN echo '#!/bin/bash\n'\
'source /opt/conda/etc/profile.d/conda.sh\n'\
'conda activate run-env\n'\
'# Start RStudio Server in the background\n'\
'/usr/lib/rstudio-server/bin/rserver &\n'\
'# Start Jupyter Notebook as rstudio user\n'\
'su - rstudio -c "source /opt/conda/etc/profile.d/conda.sh && conda activate run-env && jupyter-notebook --ip=0.0.0.0 --port=8888 --no-browser" &\n'\
'# Keep the container running\n'\
'tail -f /dev/null' > /usr/local/bin/start_services.sh

# Make the script executable
RUN chmod +x /usr/local/bin/start_services.sh

# Ensure permissions for rstudio user
RUN chown -R rstudio:rstudio /home/rstudio

# Set the working directory
WORKDIR /home/rstudio/work

# Set the default command to run the services
CMD ["/usr/local/bin/start_services.sh"]
