.PHONY = run clean img shell

SINGULARITY = /usr/local/bin/singularity
IMG = funcadelic-0.1.0.img
DEF = image.def

run:
	sudo $(SINGULARITY) run $(IMG)

clean:
	rm -f $(IMG)

img: clean
	sudo $(SINGULARITY) build $(IMG) $(DEF)

shell:
	sudo $(SINGULARITY) shell --writable -B $(shell pwd):/tmp $(IMG)
