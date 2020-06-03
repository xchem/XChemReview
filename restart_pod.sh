#!/bin/bash
module load argus
x=$(kubectl get pods | grep -e 'xchem' | awk '{print $1}')
kubectl delete pod $x
kubectl get pods
