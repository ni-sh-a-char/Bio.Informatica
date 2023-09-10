FROM python:3.7-slim-buster

WORKDIR /app

COPY . /app/

RUN apt-get update && \
    pip install --upgrade pip && \
    apt-get install -y libglib2.0-0 libsm6 libxrender1 libxext6 && \
    apt-get install -y openjdk-11-jre-headless && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    pip install -r requirements.txt

EXPOSE 80

RUN mkdir -p ~/.streamlit && \
    cp config.toml ~/.streamlit/config.toml && \
    cp credentials.toml ~/.streamlit/credentials.toml

WORKDIR /app

ENTRYPOINT ["streamlit", "run"]
CMD ["Bio.Informatica.py", "--server.enableCORS=false", "--server.enableWebsocketCompression=false", "--server.enableXsrfProtection=false"]
