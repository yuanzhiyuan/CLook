__author__ = 'yzy'
from flask import *
from CLook import config


UPLOAD_FOLDER = '/home/yuan/PycharmProjects/CLook/CLook/upload'
ALLOWED_EXTENSIONS = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif','html','rmvb'])
app = Flask(__name__)

app.secret_key = 'fuck'
app.config['aaaa'] = 'fuck u'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['ALLOWED_EXTENSIONS'] = ALLOWED_EXTENSIONS
import CLook.controller.index
import CLook.controller.test
import CLook.controller.file
import CLook.controller.project
import CLook.controller.upload
import CLook.controller.tools

# import CLook.controller.user
# import CLook.controller.article
# import CLook.controller.category
# import CLook.controller.tools
# import CLook.controller.ueditor

# if not config.DEBUG:
    # import CLook.controller.errorHandler
    # import CLook.controller.log
