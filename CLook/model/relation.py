__author__ = 'yzy'
from common import *
from hashlib import md5
import time
class Relation(Base):
    __tablename__='relation'

    ancestor=Column(Integer,primary_key=True)
    descendant=Column(Integer,primary_key=True)
    author=Column(Integer,primary_key=True)
    depth=Column(Integer)



    def addRelation(self,child,parent,author):

        # add direct relation
        if -1 == parent:
            relationObj = Relation(ancestor=parent,descendant=child,depth=0,author=author)
            session.add(relationObj)
            session.commit()
            return True
        relationObj = Relation(ancestor=parent,descendant=child,depth=1,author=author)
        session.add(relationObj)

        # add undirect relation
        undirectRows = session.query(Relation).filter(Relation.descendant==parent , author==author).all()
        if undirectRows:
            for row in undirectRows:
                ancestorid = row.ancestor
                depth = row.depth
                relationToAdd = Relation(ancestor=ancestorid,descendant=child,depth=depth+1,author=author)
                session.add(relationToAdd)
        session.commit()
        return True
    def getRelationRowsByFileList(self,author,fileList):
        handle = session.query(Relation).filter(or_(Relation.ancestor in fileList , Relation.descendant in fileList)).all()
        if not handle:
            return False
        return handle

    def getChildList(self,parentid,authorid):
        handle = session.query(Relation).filter(Relation.ancestor == parentid , Relation.depth == 1 , Relation.author == authorid).all()
        if handle:
            return map(lambda a:a.descendant,handle)
        else:
            return []


    # def addProject(self,category,species,cell_type,manager,infomation):
    #     projectsObj = Projects(category=category,species=species,cell_type=cell_type,manager=int(manager),infomation=infomation,start_date=int(time.time()),update_date=int(time.time()))
    #     if projectsObj:
    #         session.add(projectsObj)
    #         session.commit()
    #         return projectsObj.id
    #     else:
    #         return False

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




