#!/bin/bash
psql -U protwis -h localhost -d protwis -c 'drop schema public cascade; create schema public;'
python3 manage.py makemigrations
python3 manage.py migrate
psql -U protwis -h localhost -o protwis < /home/vagrant/protwis.sql;