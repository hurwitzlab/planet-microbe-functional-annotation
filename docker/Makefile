TAG = hurwitzlab/funcadelic:0.1.0
DOCKER = sudo docker

img:
	$(DOCKER) build --tag=$(TAG) .

shell:
	$(DOCKER) run -it $(TAG) bash

run:
	$(DOCKER) run -v `cd ../data && pwd`:/data $(TAG) pipeline.py -c /data/configs/docker.txt

push:
	$(DOCKER) push $(TAG)
