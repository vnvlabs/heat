
ARG FROM_IMAGE=vnv_base
FROM ${FROM_IMAGE}


COPY . /${SOURCE_DIR}/heat  

RUN  cd ${SOURCE_DIR}/heat \
     && mkdir build \
     && cd build \
     && cmake -DInjection_DIR=${VNV_DIR}/lib/cmake .. -DCMAKE_INSTALL_PREFIX=${SOFTWARE_DIR}/heat \
     && make \
     && make install  \
     && ls \
     && cp -r ../vnv/gui ${SOFTWARE_DIR}/heat/gui

#Remove the source
RUN rm -r ${SOURCE_DIR}/heat 

ENV VNV_CONFIG=${VNV_CONFIG}:${SOFTWARE_DIR}/heat/vnv.__registration__
ENV VNV_CONFIG=${VNV_CONFIG}:${SOFTWARE_DIR}/heat/gui/vnv.__registration__

