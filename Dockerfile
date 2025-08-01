FROM ubuntu:20.04

# Avoid interactive prompts (e.g., for tzdata)
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y python3.8 python3-pip samtools libbz2-dev liblzma-dev
RUN ln -sf /usr/bin/python3.8 /usr/bin/python
RUN ln -sf /usr/bin/pip3 /usr/bin/pip
RUN mkdir MSMuTect4

COPY ./src MSMuTect4/src
COPY ./msmutect.sh MSMuTect4
COPY ./msmutect.sh MSMuTect4
COPY ./build.sh MSMuTect4
COPY ./rename.sh MSMuTect4
COPY ./setup.py MSMuTect4
COPY ./requirements.txt MSMuTect4
COPY ./LICENSE MSMuTect4
#RUN cd MSMuTect4 && ls
RUN cd MSMuTect4 && pip3 install -r requirements.txt && bash build.sh


# Set working directory
WORKDIR /app

#RUN ls /
#RUN ls /MSMuTect4/
#RUN ls /MSMuTect4/msmutect.sh
# Default command
ENTRYPOINT ["/MSMuTect4/msmutect.sh"]
