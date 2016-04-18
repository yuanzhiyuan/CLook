__author__ = 'yzy'
from common import *
from hashlib import md5
import time
class Projects(Base):
    __tablename__='projects'

    id=Column(Integer,primary_key=True)
    name=Column(String(30))
    category=Column(Integer)
    species=Column(String(30))
    cell_type=Column(String(30))
    n_files=Column(Integer,default=0)
    start_date = Column(Integer)
    update_date = Column(Integer)
    manager = Column(Integer)
    infomation = Column(Text)



    # if success,return new project's id
    def addProject(self,name,category,species,cell_type,manager,infomation):
        projectsObj = Projects(name=name,category=int(category),species=species,cell_type=cell_type,manager=int(manager),infomation=infomation,start_date=int(time.time()),update_date=int(time.time()))
        if projectsObj:
            session.add(projectsObj)
            session.commit()
            return projectsObj.id
        else:
            return False

    # get categoried project list
    # e.g.: rst[0] = [rows of RNA-seq]
    #       rst[1] = [rows of chip-seq]
    #         2...3...
    def getCategoriesProjects(self):
        categoriedProjectsList = []
        for i in range(4):
            categoryid = i
            ithCategoryProjects = session.query(Projects).filter(Projects.category==categoryid).all()

            categoriedProjectsList.append(ithCategoryProjects)
        return categoriedProjectsList

    def getValueByColumnName(self,id,*listOfColumns):
        selectedRow = session.query(Projects).filter(Projects.id==id).first()
        rst = []
        for column in listOfColumns:
            if hasattr(selectedRow,column):

                rst.append(eval('selectedRow'+'.'+column))
            else:
                rst.append(False)
        return rst

    def update_n_files(self,projectid,n_files):
        handle = session.query(Projects).filter(Projects.id==projectid)
        if handle:
            handle.update({'n_files':n_files})
            session.commit()
            return True
        else:
            return False

    # def updateInfo(self,id=None,username=None,password=None,lastip=None,lasttime=None):
    #     if id:
    #         userObj=session.query(User).filter(User.id==id)
    #         if password:
    #             changeDict={'password':md5(password).hexdigest()}
    #         elif lastip and lasttime:
    #             changeDict={'lastip':lastip,'lasttime':lasttime}
    #         userObj.update(changeDict)
    #         session.commit()
    #         return True
    #     elif username:
    #         userObj=session.query(User).filter(User.username==username)
    #         if password:
    #
    #             changeDict={'password':md5(password).hexdigest()}
    #         elif lastip and lasttime:
    #             changeDict={'lastip':lastip,'lasttime':lasttime}
    #         userObj.update(changeDict)
    #         session.commit()
    #         return True
    # def validate(self,username,password):
    #
    #     user=session.query(User).filter(User.username==username).first()
    #     if user:
    #         if user.password==md5(password).hexdigest():
    #             return True
    #         else:
    #             return False
    #     else:
    #         return False
    # def getUserObj(self,username=None,id=None):
    #     if username!=None:
    #         user=session.query(User).filter(User.username==username).first()
    #         if user:
    #             return user
    #         else:
    #             return False
    #     elif id!=None:
    #         user=session.query(User).filter(User.id==id).first()
    #         if user:
    #             return user
    #         else:
    #             return False
    #
    # def addUser(self,username,password):
    #     userObj=User(username=username,password=md5(password).hexdigest())
    #     if userObj:
    #         session.add(userObj)
    #         session.commit()
    #         return True
    #     else:
    #         return False
    #
    # def getAllUsers(self):
    #     return session.query(User)
    #
    # def deleteUser(self,userid):
    #     userObj=self.getUserObj(id=userid)
    #     if userObj:
    #         session.delete(userObj)
    #         session.commit()
    #         return True
    #     else:
    #         return False
    #
    # def getUserList(self):
    #     users=self.getAllUsers()
    #     user_list=[]
    #     for user in users:
    #         user_list.append(user.username)
    #     return user_list





Base.metadata.create_all(engine)
# p = Projects()
# print p.getValueByColumnName(2,'infomatio','infomation','id','name')
# p.addProject('a',1,'b','c',1,'d')
# p.getCategoriesProjects()
# print p.getCategoriesProjects()[1]
# print p.getCategoriesProjects()[2]
# print p.getCategoriesProjects()[3]


