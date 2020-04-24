# Update stager and docker
#echo 'Deleting Old Docker'
docker images -a | grep 'staging' | awk '{print $3}' | xargs docker rmi -f
docker images -a | grep 'none' | awk '{print $3}' | xargs docker rmi -f
echo 'Grabbing New Data'
~/anaconda3/envs/staging/bin/python3 getDatafd.py
echo 'Recontainerising Data'
docker build --rm -t staging . 
docker save 'staging' > staging.tar
echo 'Loading Docker'
docker load --input staging.tar
# Run With
# docker run --rm -p 3838:3838 staging  
