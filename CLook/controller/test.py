#encoding = utf8
__author__ = 'yzy'
import time
import CLook.config as config
from werkzeug import secure_filename
import os
from CLook import app
from flask import request,render_template,session,redirect,abort,url_for,make_response,send_file,send_from_directory
# import CLook.model.user as db_user
# from CLook.controller.auth import requires_auth
import CLook.model.projects as db_projects
import CLook.model.files as db_files
import CLook.model.relation as db_relation

@app.route('/test')
def download_test():
    return send_file('upload/4.png',as_attachment=True)


@app.route('/download')
def download():
    csv = """"REVIEW_DATE","AUTHOR","ISBN","DISCOUNTED_PRICE"
"1985/01/21","Douglas Adams",0345391802,5.95
"1990/01/12","Douglas Hofstadter",0465026567,9.95
"1998/07/15","Timothy ""The Parser"" Campbell",0968411304,18.99
"1999/12/03","Richard Friedman",0060630353,5.95
"2004/10/04","Randel Helms",0879755725,4.50"""
    # We need to modify the response, so the first thing we
    # need to do is create a response out of the CSV string
    response = make_response(csv)
    # This is the key: Set the right header for the response
    # to be downloaded, instead of just printed on the browser
    response.headers["Content-Disposition"] = "attachment; filename=books.csv"
    return response

@app.route('/look')
def look():
    return render_template('test/human38.html')

@app.route('/file/static/<string:path>')
def get_static_file(path):

    # return send_from_directory('static',path)
    # return url_for('static',filename=path)
    return send_file('gene/'+path)