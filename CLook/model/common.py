__author__ = 'yzy'
import CLook.config as config

from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String,DateTime,Text,desc,or_
from sqlalchemy.orm import sessionmaker

engine = create_engine('mysql://{}:{}@localhost:3306/{}?charset=utf8'.format(config.DB_USERNAME,config.DB_PASSWORD,config.DB_DB),pool_recycle=7200,encoding = 'utf-8', echo=False)
Session = sessionmaker(bind=engine)
session = Session()
Base = declarative_base()
