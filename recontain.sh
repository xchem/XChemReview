# Update stager and docker
#echo 'Deleting Old Docker'
#docker images -a | grep 'my-shiny-app' | awk '{print $3}' | xargs docker rmi -f
echo 'Grabbing New Data'
~/anaconda3/envs/staging/bin/python3 getDatafd.py
echo 'Recontainerising Data'
docker build -t staging . 
docker save 'staging' > staging.tar
echo 'Loading Docker'
docker load --input staging.tar

# Run With
# docker run -p 3838:3838 staging  
