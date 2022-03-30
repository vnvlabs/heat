ARG FROM_IMAGE=vnv_base
FROM ${FROM_IMAGE}:latest

#Install the simple application example. 
COPY . ./heat  
RUN  cd heat && mkdir build && cd build && cmake -DInjection_DIR=${VNV_DIR}/lib/cmake .. && make 

ARG COMMIT_HASH=unknown
LABEL vnv.version.heat=${COMMIT_HASH}

ENV VNV_CONFIG=${VNV_CONFIG}:/heat/.vnv-config.json
