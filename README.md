# PRAME in NUT Carcinoma

This repository contains all scripts used to complete analysis and make figures for our pre-print (currently under review): PRAME Epitopes are T-Cell Immunovulnerabilities in BRD4::NUTM1 Initiated NUT Carcinoma (https://www.biorxiv.org/content/10.1101/2025.03.07.642090v1) 

To get started, clone this repository (can also use SSH keys if you have that set up)
```
git clone https://github.com/pirl-unc/PRAME_in_NUT_Carcinoma.git 
```

To build the Docker container
```
bash build_docker.sh
```

To run the Docker container
```
bash run_docker.sh
```

The docker container runs both Rstudio and Python/Jupyter notebook.
- To access Rstudio, navigate to localhost:8787 in your browser. **username: rstudio, password: prame**
- To access Jupyter notebook go to http://<host_machine_ip>:8888/tree?token=<token> . Where host machine IP is the IP address of the machine running the Docker container, and the token can be found in the container log in the bash shell you started the container in.There will be a link you can click to in the log that will take you to the Jupyter Notebook. 

Great!! You're in!!
