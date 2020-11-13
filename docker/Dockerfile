FROM python:3.9.0-buster
RUN apt-get -y update
RUN apt-get install -y git vim

WORKDIR /app
RUN git clone git@github.com:hurwitzlab/planet-microbe-functional-annotation.git
# && python3 -m pip install -r /app/biofx_python/requirements.txt

CMD ["python3", "--version"]
