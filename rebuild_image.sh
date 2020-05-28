#!/bin/bash
module load argus
podman image prune -a
podman build --rm -t xchemstructurereview .
x=$(podman images | grep -e 'xchem' | awk '{print $3}')
podman tag $x gcr.io/diamond-privreg/xchemapps/xchemstructurereview:latest
podman push gcr.io/diamond-privreg/xchemapps/xchemstructurereview:latest 
x=$(kubectl get pods | grep -e 'xchem' | awk '{print $1}')
kubectl delete pod $x
echo 'Finished'
