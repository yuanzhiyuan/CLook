########
# Database
########
DB_HOST='127.0.0.1'
DB_PORT='3306'
DB_DB='CLook'
DB_USERNAME='root'
DB_PASSWORD='002899'

########
#MailServer
########
# Admin=['kongkongyzt@gmail.com']
Admin=['707699544@qq.com','kongkongyzt@gmail.com','1138465446@qq.com']
MailHost="smtp.qq.com"
MailPort = 25
MailUser="364115318"
MailPass="haojingjie521"
MailAddr="364115318@qq.com"

########
# Run-Time
########
DEBUG=True


########
#Memerycache
########
Memerycache=False



logPath='./CLook/log/error.log'

########
# Upload
########

UPLOAD_FOLDER = '/home/yuan/PycharmProjects/CLook/CLook/upload'
ALLOWED_EXTENSIONS = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif','html'])