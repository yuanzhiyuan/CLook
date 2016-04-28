__author__ = 'yuan'
from CLook import app
from CLook.controller.tools import *
from flask import request,render_template,session,redirect,url_for,abort
import CLook.model.user as db_user
import CLook.model.files as db_files





@app.route('/login',methods=['GET','POST'])
def login():
    if request.method=='GET':
        return render_template('test/login.html')
    else:
        username=request.form['username']
        password=request.form['password']
        if not username:
            return 'error3'
        if not password:
            return 'error2'
        if not db_user.User().validate(username=username,password=password):
            return 'error1'

        userObj=db_user.User().getUserObj(username=username)
        if not userObj:
            return 'error4'

        session['id']=userObj.id
        session['username']=username
        session['password']=password
        lasttime=int(time.time())
        lastip=request.remote_addr
        # return lastip+lasttime
        db_user.User().updateInfo(username=username,lasttime=lasttime,lastip=lastip)
        return 'success'



@app.route('/user/add',methods=['GET','POST'])

def addUser():
    if request.method=='GET':
        return render_template('test/manage_user.html')
    else:



        username=request.form['username']
        password=request.form['password']
        repassword=request.form['repassword']
        if username=='':
            return 'please input username'
        if password=='':
            return 'please input password'
        if repassword=='':
            return 'please validate password'
        if password!=repassword:
            return 'validate fail'
        if db_user.User().getUserObj(username=username):
            return 'username exists'
        db_user.User().addUser(username=username,password=password)
        return 'success'



@app.route('/user/listfiles',methods=['GET','POST'])
def listFiles():
    if request.method=='GET':
        userid = session['id']
        files_obj = db_files.Files().getFilesByUserid(userid)
        return render_template('test/list_files.html',files_obj=files_obj)

@app.route('/user/list',methods=['GET','POST'])
def listUsers():
    if request.method=='GET':

        users_obj = db_user.User().getAllUsers()
        for i in users_obj:
            print i.id
        return render_template('test/list_users.html',users_obj=users_obj)



@app.route('/user/changepwd',methods=['GET','POST'])

def changepwd():
    if request.method == 'GET':
        return render_template('test/manage_infomation.html')
    else:



        username=session['username']
        oldpassword=request.form['oldpassword']
        if not(username!='' and oldpassword!=''):
            return 'please input username/password'
        if not db_user.User().validate(username,oldpassword):
            return 'wrong password'

        newpassword=request.form['newpassword']
        renewpassword=request.form['renewpassword']
        if not newpassword!='':
            return 'please input new password'
        if not renewpassword!='':
           return 'please validate your password'
        if not newpassword==renewpassword:
            return 'validate fail'
        if db_user.User().updateInfo(id=session['id'],password=newpassword):
            return 'success'
        else:
            return 'error!'