FROM ubuntu:rolling

# Aggiorna i certificati e il sistema
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && update-ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Installa i pacchetti necessari per OpenFHE
RUN apt-get update && apt-get install -y \
    build-essential cmake git wget libssl-dev nano less


CMD ["/bin/bash"]



# Con Mount:
# docker run -v "${PWD}/src:/app/target" -v "${PWD}/openfhe-development:/usr/local/include/openfhe" -it tiziana2 bash
