__author__ = 'yzy'
from flask import *
from CLook import config
import os

UPLOAD_FOLDER = '/home/yuan/PycharmProjects/CLook/CLook/upload'
ALLOWED_EXTENSIONS = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif','html','rmvb'])
app = Flask(__name__)
APP_ROOT = os.path.dirname(os.path.abspath(__file__))
app.config['APP_ROOT'] = APP_ROOT
# print APP_ROOT
app.secret_key = 'fuck'
app.config['aaaa'] = 'fuck u'
app.config['test_matrix'] = [
	        ('2.5MB',99 +1),
	        ('1MB',249 +1),
	        ('500KB',498+1),
	        ('250KB',996+1),
	        ('100KB',2492+1),
	        ('50KB',4984+1),
	        ('25KB',9969+1),
	        ('10KB',24922+1)
            ]
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['ALLOWED_EXTENSIONS'] = ALLOWED_EXTENSIONS
import CLook.controller.index
import CLook.controller.test
import CLook.controller.file
import CLook.controller.project
import CLook.controller.upload
import CLook.controller.tools
import CLook.controller.user

# import CLook.controller.user
# import CLook.controller.article
# import CLook.controller.category
# import CLook.controller.tools
# import CLook.controller.ueditor

# if not config.DEBUG:
    # import CLook.controller.errorHandler
    # import CLook.controller.log
