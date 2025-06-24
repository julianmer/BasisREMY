FROM python:3.11.12-slim

# install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        octave \
        liboctave-dev \
        build-essential \
        make \
        libglib2.0-0 \
        tk \
        libxext6 \
        libxrender1 \
        libsm6 \
        libx11-6 \
        libxft2 \
        libxinerama1 \
        libxcursor1 \
        libxrandr2 \
        libfreetype6 \
        libfontconfig1 \
        xvfb \
        x11-utils \
        xauth \
        x11-apps && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install Octave toolboxes
RUN octave --eval "pkg install -forge control signal" && \
    octave --eval "pkg load control; pkg load signal; savepath"

# set environment for oct2py to find octave
ENV OCTAVE_EXECUTABLE=octave

# install Python requirements
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# copy BasisREMY code
COPY . /app
WORKDIR /app

# default command to run your script
CMD ["python", "main.py"]