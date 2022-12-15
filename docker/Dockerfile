
ARG FROM_IMAGE=vnv_base
FROM ${FROM_IMAGE}

ARG COMMIT_HASH=unknown
LABEL  vnv.version.heat=${COMMIT_HASH}

ENV VNV_CONFIG=${VNV_CONFIG}:${SOURCE_DIR}/heat/build/vnv.__registration__
ENV VNV_CONFIG=${VNV_CONFIG}:${SOURCE_DIR}/heat/docker/config/vnv.__registration__

COPY . /${SOURCE_DIR}/heat  

RUN  cd ${SOURCE_DIR}/heat \
     && mkdir build \
     && cd build \
     && cmake -DInjection_DIR=${VNV_DIR}/lib/cmake .. \
     && make 



