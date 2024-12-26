ARG APP=spacedust
ARG downloader=foldseek_${TARGETARCH}_downloader
ARG FOLDSEEK_COMMIT=9c31b42b4abe3f1506fd34817d6aede307091b01

FROM scratch AS foldseek_amd64_downloader
ARG FOLDSEEK_COMMIT
WORKDIR /opt/build
ONBUILD ADD https://mmseqs.com/archive/${FOLDSEEK_COMMIT}/foldseek-linux-avx2.tar.gz  .
ONBUILD ADD https://mmseqs.com/archive/${FOLDSEEK_COMMIT}/foldseek-linux-sse41.tar.gz  .

FROM scratch AS foldseek_arm64_downloader
ARG FOLDSEEK_COMMIT
WORKDIR /opt/build
ONBUILD ADD https://mmseqs.com/archive/${FOLDSEEK_COMMIT}/foldseek-linux-arm64.tar.gz  .

FROM $downloader AS downloader

FROM --platform=$BUILDPLATFORM debian:stable-backports as builder
ARG TARGETARCH
ARG APP
ARG FOLDSEEK_COMMIT

RUN dpkg --add-architecture $TARGETARCH \
    && apt-get update \
    && apt-get install -y \
      build-essential curl xxd git cmake \
      zlib1g-dev libbz2-dev libatomic1 \
      crossbuild-essential-$TARGETARCH zlib1g-dev:$TARGETARCH libbz2-dev:$TARGETARCH \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/build
ADD . .

RUN if [ "$TARGETARCH" = "arm64" ]; then \
      mkdir -p build_$TARGETARCH/src; \
      cd /opt/build/build_$TARGETARCH; \
      CC=aarch64-linux-gnu-gcc CXX=aarch64-linux-gnu-g++ cmake -DHAVE_ARM8=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j $(nproc --all); \
      mv src/${APP} /opt/build/${APP}_arch; \
      touch /opt/build/${APP}_sse41 /opt/build/${APP}_avx2; \
    else \
      mkdir -p build_sse41/src && mkdir -p build_avx2/src; \
      cd /opt/build/build_sse41; \
      cmake -DHAVE_SSE4_1=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j $(nproc --all); \
      mv src/${APP} /opt/build/${APP}_sse41; \
      cd /opt/build/build_avx2; \
      cmake -DHAVE_AVX2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
      make -j $(nproc --all); \
      mv src/${APP} /opt/build/${APP}_avx2; \
      touch /opt/build/${APP}_arch; \
    fi

ADD https://raw.githubusercontent.com/steineggerlab/foldseek/${FOLDSEEK_COMMIT}/util/foldseek_wrapper.sh binaries/foldseek

COPY --from=downloader /opt/build/* .
RUN mkdir binaries; \
    if [ "$TARGETARCH" = "arm64" ]; then \
      for i in foldseek; do \
        if [ -e "${i}-linux-arm64.tar.gz" ]; then \
          cat ${i}-linux-arm64.tar.gz | tar -xzvf- ${i}/bin/${i}; \
          mv -f -- ${i}/bin/${i} binaries/${i}; \
        fi; \
      done; \
    else \
      for i in foldseek; do \
        for j in sse41 avx2; do \
          if [ -e "${i}-linux-${j}.tar.gz" ]; then \
            cat ${i}-linux-${j}.tar.gz | tar -xzvf- ${i}/bin/${i}; \
            mv -f -- ${i}/bin/${i} binaries/${i}_${j}; \
          fi; \
        done; \
      done; \
    fi; \
    chmod -R +x binaries;

FROM debian:stable-slim
ARG TARGETARCH
ARG APP

RUN apt-get update && apt-get install -y \
      gawk bash grep libstdc++6 libgomp1 libatomic1 zlib1g libbz2-1.0 wget tar aria2 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/build/${APP}_arch /opt/build/${APP}_sse41 /opt/build/${APP}_avx2 /opt/build/binaries /usr/local/bin/
ADD util/${APP}_wrapper.sh /usr/local/bin/entrypoint
RUN if [ "$TARGETARCH" = "arm64" ]; then rm -f /usr/local/bin/entrypoint; ln -s /usr/local/bin/${APP}_arch /usr/local/bin/entrypoint; fi

ENTRYPOINT ["/usr/local/bin/entrypoint"]

