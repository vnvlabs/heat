ARG FROM_IMAGE=vnv_base
FROM ${FROM_IMAGE}:latest

ARG COMMIT_HASH=unknown
LABEL  vnv.version.heat=${COMMIT_HASH}


COPY . /${SOURCE_DIR}/heat  

RUN  cd ${SOURCE_DIR}/heat \
     && mkdir build \
     && cd build \
     && cmake -DInjection_DIR=${VNV_DIR}/lib/cmake .. \
     && make 


ENV VNV_CONFIG=${VNV_CONFIG}:${SOURCE_DIR}/heat/build/vnv.__registration__
