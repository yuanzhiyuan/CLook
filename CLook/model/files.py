__author__ = 'yzy'
from common import *
from hashlib import md5
import time
class Files(Base):
    __tablename__='files'

    id=Column(Integer,primary_key=True)
    name=Column(String(30))
    location=Column(String(255))
    project=Column(Integer)
    type=Column(Integer)
    file_type=Column(String(30))
    author = Column(Integer)
    date = Column(Integer)
    comments = Column(Text)

    def addFile(self,name,location,project,type,file_type,author,infomation):
        filesObj = Files(name=name,location=location,project=int(project),type=int(type),file_type=file_type,author=int(author),comments=infomation,date=int(time.time()))
        if filesObj:
            session.add(filesObj)
            session.commit()
            return filesObj.id
        else:
            return False

    def getRawDataByProjectid(self,projectid):
        selectedRow = session.query(Files).filter(Files.project==projectid , Files.type==0).all()
        if not selectedRow:
            return False
        return selectedRow

        # get users of rawdata in a project
    def getRawDataUsersByProjectid(self,projectid):
        selectedRows = session.query(Files).filter(Files.project==projectid).all()
        if not selectedRows:
            return []
        return selectedRows

    def get_n_files(self,projectid):
        handle = session.query(Files).filter(Files.project==projectid).all()
        if not handle:
            return 0
        else:
            return len(handle)

    def getFileByFileid(self,fileid):
        handle = session.query(Files).filter(Files.id==fileid).first()
        return handle

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

# f = Files()
#
# print f.getRawDataByProjectid(100)

