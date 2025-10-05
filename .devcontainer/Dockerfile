FROM python:3.9-slim

# Install system dependencies including headless rendering support
RUN apt-get update && apt-get install -y \
    git \
    build-essential \
    gcc \
    g++ \
    make \
    cmake \
    pkg-config \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    libegl1-mesa-dev \
    freeglut3-dev \
    xvfb \
    x11-utils \
    libx11-dev \
    libxrandr-dev \
    libxinerama-dev \
    libxcursor-dev \
    libxi-dev \
    curl \
    vim \
    sudo \
    && rm -rf /var/lib/apt/lists/*

# Set up virtual display for headless rendering
ENV DISPLAY=:99

# Create workspace directory
WORKDIR /workspace

# Set up permissions
RUN chmod 777 /workspace

# Create a script to start Xvfb automatically
RUN echo '#!/bin/bash\nXvfb :99 -screen 0 1024x768x24 > /dev/null 2>&1 &\nexec "$@"' > /usr/local/bin/start-xvfb
RUN chmod +x /usr/local/bin/start-xvfb

# Set entrypoint to start virtual display
ENTRYPOINT ["/usr/local/bin/start-xvfb"]
CMD ["bash"]