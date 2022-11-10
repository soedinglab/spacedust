FROM debian:stable-slim as spacedust-builder

RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    build-essential cmake xxd git zlib1g-dev libbz2-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/spacedust
ADD . .

RUN mkdir -p build_sse/bin && mkdir -p build_avx/bin

WORKDIR /opt/spacedust/build_sse
RUN cmake -DHAVE_SSE4_1=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
    make -j $(nproc --all) && make install;

WORKDIR /opt/spacedust/build_avx
RUN cmake -DHAVE_AVX2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
    make -j $(nproc --all) && make install;

FROM debian:stable-slim
MAINTAINER Milot Mirdita <milot@mirdita.de>

RUN apt-get update && apt-get upgrade -y && apt-get install -y \
     gawk bash grep wget libstdc++6 libgomp1 zlib1g libbz2-1.0 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=spacedust-builder /opt/spacedust/build_sse/bin/spacedust /usr/local/bin/spacedust_sse42
COPY --from=spacedust-builder /opt/spacedust/build_avx/bin/spacedust /usr/local/bin/spacedust_avx2
RUN echo '#!/bin/bash\n\
if $(grep -q -E "^flags.+avx2" /proc/cpuinfo); then\n\
    exec /usr/local/bin/spacedust_avx2 "$@"\n\
else\n\
    exec /usr/local/bin/spacedust_sse42 "$@"\n\
fi' > /usr/local/bin/spacedust
RUN chmod +x /usr/local/bin/spacedust

ENTRYPOINT ["/usr/local/bin/spacedust"]
