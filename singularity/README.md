# Building the Singularity container

First, you must build/test/push the Docker container to the "hurwitzlab" account. 
See the "../docker" directory.

Next run the following command to build the Singularity container from the base Docker image:

```
$ make img
```

## Author

Ken Youens-Clark <kyclark@arizona.edu>
