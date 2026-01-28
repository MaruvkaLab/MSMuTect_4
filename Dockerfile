FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y samtools gcc g++ libopenblas-dev liblapack-dev procps


COPY . /MSMuTect4
WORKDIR /MSMuTect4

RUN pip install --upgrade pip setuptools wheel && pip install -r requirements.txt
RUN python3 -c "import scipy"
RUN bash build_cython.sh

# Set working directory
ENTRYPOINT ["/MSMuTect4/msmutect.sh"]
