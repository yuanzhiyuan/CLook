from CLook import app
from flask import request,render_template,abort
import CLook.model.projects as db_projects
import CLook.model.files as db_files
import CLook.model.relation as db_relation
@app.route('/addProject',methods=['GET','POST'])
def addProject():
    if request.method=='GET':
        return render_template('test/add_project.html')
    elif request.method=='POST':
        name=request.form['name']
        category = request.form['category']
        species = request.form['species']
        cell_type = request.form['cell_type']
        manager_no = request.form['manager_no']
        infomation = request.form['infomation']
        ok = db_projects.Projects().addProject(name,category,species,cell_type,manager_no,infomation)
        if not ok:
            return 'fail'
        else:
            return 'success'
    else:
        return abort(500)


@app.route('/exp',methods=['GET'])
@app.route('/exp/<int:id>',methods=['GET'])
def exp(id=0):
    if request.method=='GET':

        projectInfomationName = db_projects.Projects().getValueByColumnName(id,'infomation','name')
        if False in projectInfomationName:
            return abort(500)
        rawDatas = db_files.Files().getRawDataByProjectid(id)
        if not rawDatas:
            rawDatas = []
# dataUsersRows = db_files.Files().getRawDataUsersByProjectid(projectid)
        rawDataUsers = db_files.Files().getRawDataUsersByProjectid(id)
        id_thisProjectFileObj_map = {}
        for row in rawDataUsers:
            id_thisProjectFileObj_map[int(row.id)] = []
            id_thisProjectFileObj_map[int(row.id)].append(row.location)
            id_thisProjectFileObj_map[int(row.id)].append(row.name)

        return render_template('test/Experiment-template.html',files=id_thisProjectFileObj_map,projectid=id,projectInfomationName=projectInfomationName,rawDatas=rawDatas)

@app.route('/exp/<int:projectid>/getTree',methods=['GET','POST'])
def getFullTree(projectid):
    if request.method=='POST':
    # if 1:
        dataUsersRows = db_files.Files().getRawDataUsersByProjectid(projectid)
        id_thisProjectFileObj_map = {}
        for row in dataUsersRows:
            id_thisProjectFileObj_map[str(row.id)] = []
            id_thisProjectFileObj_map[str(row.id)].append(row.location.encode('utf8'))
            id_thisProjectFileObj_map[str(row.id)].append(row.name.encode('utf8'))
        # maybe duplicate

        if len(dataUsersRows)==0:
            return 'no data'
        userList = map(lambda x:x.author,dataUsersRows)
        print userList
        userSet = set(userList)
        tree_structure = {}
        for author in userSet:
            tree_structure[str(author)] = {}
            rowsOfThisAuthor = filter(lambda row:row.author==author,dataUsersRows)
            # if not rowsOfThisAuthor:
                # continue
            # fileList = map(lambda row:row.id,rowsOfThisAuthor)
            # relationRows = db_relation.Relation().getRelationRowsByFileList(author,fileList)

            rawDataList = filter(lambda row:row.type==0,rowsOfThisAuthor)

            #record each level's files
            thisAuthorRawDataList = map(lambda a:a.id,rawDataList)
            thisAuthorTree = {}
            for rawdata in thisAuthorRawDataList:
                thisAuthorTree[str(rawdata)] = getTree(1,rawdata,author)
            tree_structure[str(author)] = thisAuthorTree

        rst = {}
        rst['file_map'] = id_thisProjectFileObj_map
        rst['tree_structure'] = tree_structure


        return str(rst).replace("'",'"').replace('L','')

            # 4 types of files
            # for i in range(1,4):
            #     for relationRow in relationRows:
            #         if relationRow.ancestor in thisLevelFileList and relationRow.depth==1:
            #             tree_structure[author][relationRow.ancestor][relationRow.descendant]

# level:
# 0 - author
# 1 - raw data
# 2 - qc report


def getTree(level,parentid,authorid):
    if level == 3:
        childList = db_relation.Relation().getChildList(parentid,authorid)
        rst = {}
        for child in childList:
            subtree = {}
            rst[str(child)] = subtree
        return rst
    else:
        childList = db_relation.Relation().getChildList(parentid,authorid)
        rst = {}
        if parentid==18:
            print childList
        for child in childList:
            subtree = getTree(level+1,child,authorid)
            rst[str(child)] = subtree
        return rst









