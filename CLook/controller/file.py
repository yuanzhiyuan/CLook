__author__ = 'yuan'
from CLook import app
from CLook.controller.tools import *
from flask import request,session,render_template,abort,url_for,make_response,send_file
import CLook.model.relation as db_relation
import CLook.model.files as db_files
import CLook.model.projects as db_projects

@app.route('/addFile',methods=['GET','POST'])
def addFile():
    if request.method=='GET':
        if session.get('uploaded_file_location'):
            session.pop('uploaded_file_location')
        return render_template('test/add_file.html'
        )
    elif request.method=='POST':
        print session['abc']
        checkLock()
        name = request.form['name']

        location = None


        print 'session',session.get('uploaded_file_location')
        if session.get('uploaded_file_location'):
            location = '&'.join(session['uploaded_file_location'])
            session.pop('uploaded_file_location')
        project = request.form['project']
        type = request.form['type']
        file_type = request.form['file_type']
        author = request.form['author']
        under = request.form['under']
        infomation = request.form['infomation']




        # addFile must return the id of file just added,in order to update the location db.
        ok = db_files.Files().addFile(name,location,project,type,file_type,author,infomation)
        if ok==False:
            return 'fail'
        db_relation.Relation().addRelation(int(ok),int(under),int(author))

        # update the n_file of model projects
        n_files = db_files.Files().get_n_files(project)
        db_projects.Projects().update_n_files(project,n_files)


        return 'success'
    else:
        return abort(500)


@app.route('/file/<int:fileid>')
def view_file(fileid):
    fileObj = db_files.Files().getFileByFileid(fileid)
    if fileObj:
        file_location_str = fileObj.location
        file_location_list = file_location_str.split('&')
        filename_url_map = {}
        for f in file_location_list:
            encoded_filename = f.split('/')[-1]
            if len(encoded_filename) <= len('1457660751'):
                continue
            decoded_filename = encoded_filename[len('1457660751'):]
            filename_url_map[decoded_filename] = encoded_filename


        return render_template('test/view_file.html',fileObj=fileObj,filename_url_map=filename_url_map,comments = fileObj.comments.strip())

    else:
        return abort(404)


# app.add_url_rule('/upload/<path:filename>',endpoint='attachment',build_only=True)
# @app.route('/download/<path:filename>')
# def Create_url2(filename):
#     return url_for('static',filename=filename)
@app.route('/file/download/<string:encoded_filename>')
def download_file(encoded_filename):
    if encoded_filename and len(encoded_filename):
        full_filename = 'upload/'  +  encoded_filename
        return send_file(full_filename,as_attachment=True,attachment_filename=encoded_filename[10:])
    else:
        return abort(404)


