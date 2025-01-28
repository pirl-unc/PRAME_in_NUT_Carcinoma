#!/bin/bash

winpty docker run -it --rm \
        -v $(pwd):/home/rstudio/work \
	-e PASSWORD=prame \
	-p 8787:8787 \
	-t prame



