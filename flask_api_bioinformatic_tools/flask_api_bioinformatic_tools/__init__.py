import os
from flask import Flask

''' Application Factory
Main place for configuration, registration and other setups.
__name__ = application name
instance_relative_config=True = configuration files relative to "instance" folder
'''
def create_app(test_config=None):
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SECRET_KEY='dev1', #should be overriden with a random value if deployed
    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile('config.py', silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    from . import fmafft
    app.register_blueprint(fmafft.fmafft)
    app.add_url_rule('/',endpoint='index')


    return app