FROM nginx:1.21-alpine@sha256:529db430e042ecef071f2e88267cee6da18f8ab44d66a0c44348886fdc2e60fc

RUN rm /etc/nginx/conf.d/default.conf
RUN mkdir -p /blast/reciprocal_blast/assets /blast/reciprocal_blast/media /var/cache/nginx /var/run && \
    chown -R nginx:nginx /blast /var/cache/nginx /var/run
COPY ./celery_blast/static /blast/reciprocal_blast/assets
COPY ./nginx/nginx.conf /etc/nginx/conf.d
RUN chown -R nginx:nginx /blast/reciprocal_blast/assets /etc/nginx/conf.d

USER nginx
