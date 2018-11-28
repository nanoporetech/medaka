FROM ubuntu:16.04 
RUN apt-get update && apt-get install -y python3-setuptools curl 
RUN curl https://bootstrap.pypa.io/get-pip.py -o /root/get-pip.py
RUN python3 /root/get-pip.py
RUN pip3 install --upgrade pip

COPY for_docker /root/medaka
RUN cd /root/medaka && pip3 install --prefer-binary -r requirements.txt && python3 setup.py install

CMD medaka_consensus -i basecalls.fa.gz -d draft.fa -o medaka_consensus -t 8 
