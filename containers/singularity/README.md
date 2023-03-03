## Create a singularity image that includes R libraries needed for this pipeline

This directory contains a singularity definition file to create a singularity image containing the R libraries needed by scripts in this pipeline.  It is based on the [rocker/r-ver docker image](https://hub.docker.com/r/rocker/r-ver) and adds on tidyverse packages and a couple additional packages.  

#### To modify and create the singularity image

Modify the `r_tools.def` file as appropriate

```
# create the singularity image for R tools
sudo singularity build r_tools.sif r_tools.def
# now sign the image.  Uses linked github email.  Passphrase != access token.
singularity sign r_tools.sif

# create the singularity image for python tools
sudo singularity build python_tools.sif python_tools.def
# now sign the image.  Uses linked github email.  Passphrase != access token.
singularity sign python_tools.sif
```


#### To push this image to the singularity library 

I first had to create a singularity library account and login as described [here](https://sylabs.io/guides/latest/user-guide/cloud_library.html?highlight=push#overview)

Do a remote login to singularity cloud library
```
# this uses access token created at sylabs cloud library site.  These expire every few months so may need to create a new one.
singularity remote login
```

The can push images
```
# modify version number as appropriate
singularity push -D "A singularity image containing tidyverse and a few other R packages" r_tools.sif library://stenglein-lab/r_tools/r_tools:1.0.0

# push the python tools image
singularity push -D "A singularity image containing python 3.10.6 and several modules" python_tools.sif library://stenglein-lab/python_tools/python_tools:1.0.0
```

