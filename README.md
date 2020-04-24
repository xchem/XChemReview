# staging
Shiny App to collect responses for XChem Experiments

## Setup
### Requirements
* Python3 environment with django, pandas, xchem-db and psycopg-binary installed
* Relevant setting.py and setup_django.py files containing keys to xchem-db in the root of the directory
* Docker installed and running :)

### Installation
Configure `recontain.sh` to use the python env you want and then run:
```
bash recontain.sh
``` 
Which should populate the data file with structures we want to stage, build and save the docker to the root of the directory and finally load the containered app. 

Alternatively you can do it step by step: 

Populate the Data folder with the summary of crystals that are to be staged.
```
<venv>/bin/python3 getDatafd.py
```

Use docker to contain and load the shinyApp
```
docker build -t staging . 
docker save 'staging' > staging.tar
docker load --input staging.tar
```

## Running
To Run the docker on a server simply run, configure ports accordingly in `app.R` and in `Dockerfile`
```
docker run -p 3838:3838 staging
```

## Using staging
Click on buttons to make things happen, it shouldn't crash!

## Extracting Responses from container?
To do...
