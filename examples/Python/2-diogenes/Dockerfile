FROM python:3

RUN pip3 install numpy
RUN pip3 install pandas
RUN pip3 install sklearn

RUN mkdir -p /scripts/
COPY scripts/* /scripts/
WORKDIR /scripts/

CMD [ "python", "./glucose_dynamics.py" ]
