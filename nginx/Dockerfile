FROM nginx:1.21-alpine

RUN rm /etc/nginx/conf.d/default.conf
RUN mkdir /blast/
RUN mkdir /blast/reciprocal_blast/
RUN mkdir /blast/reciprocal_blast/assets
COPY ./celery_blast/static /blast/reciprocal_blast/assets
COPY ./nginx/nginx.conf /etc/nginx/conf.d