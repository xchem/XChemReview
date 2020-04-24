# Update stager and docker
#echo 'Deleting Old Docker'
#docker images -a | grep 'my-shiny-app' | awk '{print $3}' | xargs docker rmi -f
#echo 'Grabbing New Data'
#python3 getDatafd.py
echo 'Recontainerising Data'
docker build -t my-shiny-app . 
docker save 'my-shiny-app' > my-shiny-app.tar
echo 'Loading Docker'
docker load --input my-shiny-app.tar

# Run With
# docker run -p 3838:3838 my-shiny-app  
