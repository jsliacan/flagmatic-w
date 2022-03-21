docker build -t flagmatic . 
docker run -it -u $(id -u):$(id -g) -v /home/cspiegel/git-repos/flagmatic/notebooks:/home/sage flagmatic