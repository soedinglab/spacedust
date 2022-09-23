FROM debian:stable-slim as clustersearch-builder

RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    build-essential cmake xxd git zlib1g-dev libbz2-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/clustersearch
ADD . .

RUN mkdir -p build_sse/bin && mkdir -p build_avx/bin

WORKDIR /opt/clustersearch/build_sse
RUN cmake -DHAVE_SSE4_1=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
    make -j $(nproc --all) && make install;

WORKDIR /opt/clustersearch/build_avx
RUN cmake -DHAVE_AVX2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
    make -j $(nproc --all) && make install;

FROM debian:stable-slim
MAINTAINER Milot Mirdita <milot@mirdita.de>

RUN apt-get update && apt-get upgrade -y && apt-get install -y \
     gawk bash grep wget libstdc++6 libgomp1 zlib1g libbz2-1.0 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=clustersearch-builder /opt/clustersearch/build_sse/bin/clustersearch /usr/local/bin/clustersearch_sse42
COPY --from=clustersearch-builder /opt/clustersearch/build_avx/bin/clustersearch /usr/local/bin/clustersearch_avx2
RUN echo '#!/bin/bash\n\
if $(grep -q -E "^flags.+avx2" /proc/cpuinfo); then\n\
    exec /usr/local/bin/clustersearch_avx2 "$@"\n\
else\n\
    exec /usr/local/bin/clustersearch_sse42 "$@"\n\
fi' > /usr/local/bin/clustersearch
RUN chmod +x /usr/local/bin/clustersearch

ENTRYPOINT ["/usr/local/bin/clustersearch"]
