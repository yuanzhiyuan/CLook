#encoding = utf8
__author__ = 'yzy'
import time
import CLook.config as config
from werkzeug import secure_filename
import os
from CLook import app
from flask import request,render_template,session,redirect,abort,url_for
# import CLook.model.user as db_user
# from CLook.controller.auth import requires_auth
import CLook.model.projects as db_projects
from CLook.controller.category import *
import CLook.model.files as db_files
import CLook.model.relation as db_relation

@app.route('/')
def index():
    categoriedProjects = db_projects.Projects().getCategoriesProjects()
    return render_template('test/index-template.html',category_map=category_map,categoriedProjects=categoriedProjects)



        # return str(li).replace("'",'"')










# print app.aaaa
# print app.config['aaaa']