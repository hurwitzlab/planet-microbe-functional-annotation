TAG = hurwitzlab/funcadelic:0.1.0
DOCKER = sudo docker

img:
	$(DOCKER) build --tag=$(TAG) .

shell:
	$(DOCKER) run -it $(TAG) bash

push:
	$(DOCKER) push $(TAG)
