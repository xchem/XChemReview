# XChemReview
A shiny app to collect the responses from XChemExperiments

## Setup
### Requirements
* Podman/Docker (The scripts int this are setup to use podman but you can replace podman with docker)
* Access to a database containing various information about xchem experiments.
* In this example a kubernetes cluster...

### Installation
If rebuilding the image configure the `kube.yaml`, `Dockerfile` and `rebuild_image.sh` correctly, specifying the repository you want to store you image into. Then run `rebuild_image.sh` to clean up images, rebuild and then redeploy.

If you have made changes to the `tjgorrie/nglshiny` or `app.R` run `restart_pod.sh` to restart the kubernetes pod.

## Running
While connected to the VPN you should be able to access xchemreview.diamond.ac.uk and use the website.

## Using staging
Click on buttons to make things happen, it shouldn't crash! If it does, email me..

