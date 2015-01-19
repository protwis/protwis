PROTWIS
====

The next generation of GPCRDB and related systems

INSTALLATION
---

*A vagrant configuration with all dependecies is in development!*

The system is based on **[Django][1] 1.7** and **[Python][2] 3**

On Ubuntu 14.04 you will need the following packages:

* postgresql (install with apt)
* postgresql-contrib (install with apt)
* python3-psycopg2 (install with apt)
* python3-pip (install with apt)
* django (install with pip3)
* django-tastypie (install with pip3)

Clone this repository:

    git clone https://bitbucket.org/gpcr/protwis.git
    cd protwis

Copy the example settings file and name it `settings.py`:

    mv protwis/settings_example.py protwis/settings.py

Create a postgresql database and enter the credentials to access it in `settings.py`

To create database tables, run the following from within the repository root:
    
    python3 manage.py migrate

Start the development webserver:

    python3 manage.py runserver 0.0.0.0:8000

Access the server at `localhost:8000`

[1]: http://djangoproject.com
[2]: http://python.org