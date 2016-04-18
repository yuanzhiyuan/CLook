__author__ = 'yuan'
from CLook import app
from CLook.controller.tools import *
from flask import request,render_template,session,redirect,url_for,abort
from werkzeug import secure_filename
import os
import time

@app.route('/upload',methods=['POST'])
def upload_file():
    print 'upload file enter'
    if request.method=='POST':
        session['abc'] = 'abc'
        addLock()
        print 'upload clicked'
        file = request.files['file']
        print file
        if file and allowed_file(file.filename):
            former_filename = secure_filename(file.filename)
            filename_pre = int(time.time())
            encoded_filename = str(filename_pre) + former_filename
            location = os.path.join(app.config['UPLOAD_FOLDER'],encoded_filename)


            if not session.get('uploaded_file_location'):
                session['uploaded_file_location'] = []
            print location
            session['uploaded_file_location'].append(location)
            print 'upload session',session.get('uploaded_file_location')
            print session.has_key('uploaded_file_location')
            file.save(location)

            # print session['uploaded_file_location']
            # return redirect(url_for('upload_file',filename=encoded_filename))
        print 'success'
        releaseLock()
        return 'success'
    else:
        return abort(400)

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.',1)[1] in app.config['ALLOWED_EXTENSIONS']
# print getFullTree(2)
