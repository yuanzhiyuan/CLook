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


@app.route('/matrix')
def draw_matrix():
    return render_template('test/draw_matrix.html')

@app.route('/matrix/data',methods=['POST'])
def get_matrix_data():
    if request.method=='POST':

        # print app.config['APP_ROOT']
    # with open('/upload')
        upload_dir = os.path.join(app.config['APP_ROOT'],'upload')
        to_send = ''
        max_val = 0
        min_val = 1000
        val_li = []
        # with open(os.path.join(upload_dir,'chr1_500KB_norm.txt')) as f:
        # with open(os.path.join(upload_dir,'chr1_2.5MB_norm.txt')) as f:
        with open(os.path.join(upload_dir,'chr1_250KB_norm.txt')) as f:
            current_line_str = f.next()
            for j in range(996):
                for i in range(0,j+1):
                    # current_line_str = f.next()
                    current_line_li = current_line_str.split('\t')
                    if j!= int(current_line_li[1]) or i!=int(current_line_li[0]):
                        to_send += '0'
                    else:
                        to_send += str(int(float(current_line_li[2][:-1])))
                        val_li.append(int(float(current_line_li[2][:-1])))
                        max_val = max(int(float(current_line_li[2][:-1])),max_val)
                        min_val = min(int(float(current_line_li[2][:-1])),min_val)
                        print current_line_str
                        current_line_str = f.next()
                    to_send += '&'
            to_send += '$'
            to_send += str(max_val)
            to_send += '&'
            to_send += str(min_val)
            to_send += '&'
            val_li.sort()
            x_percent_point = int(len(val_li)*0.80)
            x_percent_value = val_li[x_percent_point]
            to_send += str(x_percent_value)
        return to_send





