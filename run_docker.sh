#!/bin/bash

docker run -it --rm \
        -v $(pwd):/home/rstudio/work \
	-p 8787:8787 \
	-p 8888:8888 \
	-t prame



