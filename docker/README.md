# How to make the Docker container

Before building the image, you must have a Linux binary of the `pear` executable in this directory so that it can be added to the image.

Next, ensure Docker daemon is running, e.g., on lytic:

```
$ sudo service docker start
```

To create a Docker image, either run:

```
$ make img
```

or the equivalent:

```
$ sudo docker build --tag=hurwitzlab/funcadelic:0.1.0 .
```

You can get an interactive terminal using to verify that all is well:

```
$ make shell
```

**TODO**: We need something here about testing the container, running samples, etc..

Once you have everything working, login to Docker Hub and "push" the image to the "hurwitzlab" account:

```
$ make push
```

or

```
$ docker push hurwitzlab/funcadelic:0.1.0
```

The image must be pushed in order to build the Singularity container from it.

## Author

Ken Youens-Clark <kyclark@arizona.edu>
