#encoding = utf8
__author__ = 'yzy'
import time
import CLook.config as config
from werkzeug import secure_filename
import os
from CLook import app
from flask import request,render_template,session,redirect,abort,url_for,make_response,send_file,send_from_directory
# import CLook.model.user as db_user
# from CLook.controller.auth import requires_auth
import CLook.model.projects as db_projects
import CLook.model.files as db_files
import CLook.model.relation as db_relation

@app.route('/test')
def download_test():
    return send_file('upload/4.png',as_attachment=True)


@app.route('/download')
def download():
    csv = """"REVIEW_DATE","AUTHOR","ISBN","DISCOUNTED_PRICE"
"1985/01/21","Douglas Adams",0345391802,5.95
"1990/01/12","Douglas Hofstadter",0465026567,9.95
"1998/07/15","Timothy ""The Parser"" Campbell",0968411304,18.99
"1999/12/03","Richard Friedman",0060630353,5.95
"2004/10/04","Randel Helms",0879755725,4.50"""
    # We need to modify the response, so the first thing we
    # need to do is create a response out of the CSV string
    response = make_response(csv)
    # This is the key: Set the right header for the response
    # to be downloaded, instead of just printed on the browser
    response.headers["Content-Disposition"] = "attachment; filename=books.csv"
    return response

@app.route('/look')
def look():
    return render_template('test/human38.html')

@app.route('/file/static/<string:path>')
def get_static_file(path):

    # return send_from_directory('static',path)
    # return url_for('static',filename=path)
    return send_file('gene/'+path)


@app.route('/matrix')
def draw_matrix():
    return render_template('test/draw_matrix.html')

@app.route('/matrix/data',methods=['POST'])
def get_matrix_data():
    if request.method=='POST':

        # print app.config['APP_ROOT']
    # with open('/upload')
        upload_dir = os.path.join(app.config['APP_ROOT'],'upload')
        to_send = ''
        max_val = 0
        min_val = 1000
        val_li = []
        bin_sz_li = app.config['test_matrix']
        # with open(os.path.join(upload_dir,'chr1_500KB_norm.txt')) as f:
        # with open(os.path.join(upload_dir,'chr1_2.5MB_norm.txt')) as f:

        nth_test = 2

        with open(os.path.join(upload_dir,'chr1_{name}_norm.txt'.format(name=bin_sz_li[nth_test][0]))) as f:
            current_line_str = f.next()
            for j in range(bin_sz_li[nth_test][1]):
                for i in range(0,j+1):
                    # current_line_str = f.next()
                    current_line_li = current_line_str.split('\t')
                    if j!= int(current_line_li[1]) or i!=int(current_line_li[0]):
                        to_send += '0'
                    else:
                        to_send += str(int(float(current_line_li[2][:-1])))
                        val_li.append(int(float(current_line_li[2][:-1])))
                        max_val = max(int(float(current_line_li[2][:-1])),max_val)
                        min_val = min(int(float(current_line_li[2][:-1])),min_val)
                        # print current_line_str
                        current_line_str = f.next()
                    to_send += '&'
            to_send += '$'
            to_send += str(max_val)
            to_send += '&'
            to_send += str(min_val)
            to_send += '&'
            val_li.sort()
            x_percent_point = int(len(val_li)*0.80)
            x_percent_value = val_li[x_percent_point]
            to_send += str(x_percent_value)
            to_send += '$'
            to_send += str(bin_sz_li[nth_test][1])
        return to_send




@app.route('/matrix/data/test1',methods=['POST'])
def get_data():
    if request.method=='POST':
        bin_sz_li = app.config['test_matrix']

        if request.form['type'] == 'init_data':
            # get matrix rows number and nine position point
            to_send = []
            upload_dir = os.path.join(app.config['APP_ROOT'],'upload')
            val_li = []
            nth_test = 2
            with open(os.path.join(upload_dir,'chr1_{name}_norm.txt'.format(name=bin_sz_li[nth_test][0]))) as f:
                for line in f:
                    line_li = line.split('\t')
                    line_val = int(float(line_li[2][:-1]))
                    val_li.append(line_val)
            val_li.sort()
            x_percent_point = int(len(val_li)*0.9)
            x_percent_value = val_li[x_percent_point]
            to_send.append(str(bin_sz_li[nth_test][1]))
            to_send.append(str(x_percent_value))
            return '&'.join(to_send)


@app.route('/matrix/data/get',methods=['POST'])
def get_block_data():
    if request.method=='POST':


        bin_sz_li = app.config['test_matrix']
        nth_test = 2


        x1 = int(float(request.form['x1']))
        y1 = int(float(request.form['y1']))
        x2 = int(float(request.form['x2']))
        y2 = int(float(request.form['y2']))

        upload_dir = os.path.join(app.config['APP_ROOT'],'upload')


        if x1>x2 or y1 > y2:
            return 'error 0'
        if (x2-x1-y2+y1)!=0:
            # ensure the block is square
            return 'error 1'


        nrows = x2-x1+1

        f = open(os.path.join(upload_dir,'chr1_{name}_norm.txt'.format(name=bin_sz_li[nth_test][0])))
        # remember to close f
        if x2 <= y1:
            rst = easiest_condition(x1,y1,x2,y2,f)
        # easiest condition: the block is above the diagonal line
        # with open(os.path.join(upload_dir,'chr1_{name}_norm.txt'.format(name=bin_sz_li[nth_test][0]))) as f:
        #     for line in f:


        elif y2 <= x1:
            # just opposite the easiest condition
            rst = easiest_condition(y1,x1,y2,x2,f)

            for i in range(nrows):
                for j in range(i):
                    k = i*nrows + j
                    k_ = j*nrows + i
                    # print k,k_
                    # print k,k_
                    rst[k_],rst[k] = rst[k],rst[k_]

        elif (x2 > y1) and (x1 < y1):
            rst = diagnal_cross_condition(x1,y1,x2,y2,f)
            # the block across the diagonal line, but mainly above the diagonal line
        elif (x1 < y2) and (x2 > y2):
            rst = diagnal_cross_condition(y1,x1,y2,x2,f)

            for i in range(nrows):
                for j in range(i):
                    k = i*nrows + j
                    k_ = j*nrows + i
                    # print k,k_
                    # print k,k_
                    rst[k_],rst[k] = rst[k],rst[k_]
            # opposite the former condition
        else:
            #x1==y1 and x2==y2
            rst = diagnal_cross_condition(x1,y1,x2,y2,f)

        f.close()

        rst = map(str,rst)
        # find the xpercent point
        to_send = '&'.join(rst)
        # rst.sort()
        # x_percent_point = int(len(rst)*0.9)
        # x_percent_value = rst[x_percent_point]
        # to_send += '$'
        # to_send += str(x_percent_value)

        return to_send




def parse_line(line):
    line_li = line.split('\t')
    line_x = int(line_li[0])
    line_y = int(line_li[1])
    line_v = int(float(line_li[2][:-1]))
    return line_x,line_y,line_v



# the index of (x,y) in square(x1,y1,x2,y2) is (y-y1)*(x2-x1+1)+(x-x1)
def get_array_index(x,y,x1,y1,x2,y2):
    return (y-y1)*(x2-x1+1)+(x-x1)

#use in condition that diagnal line pass the square
# return array's order is : order of knowing the value
def get_folded_part_of_square(x1,y1,x2,y2):
    rst = []
    for i in range(y1+1,x2+1):
        for j in range(y1,i):
            rst.append((i,j))
    return rst


def easiest_condition(x1,y1,x2,y2,f):
    rst = []

    # cur_x and cur_y point at the position to be given value
    cur_x = x1
    cur_y = y1

    #put the cur_x and cur_y,   and make sure line_y >= cur_y
    # for line in f:
    #     line_x,line_y,line_v = parse_line(line)[:]
    #     if line_y == y1:
    #         if line_x == x1:
    #             rst.append(line_v)
    #             cur_x += 1
    #             if cur_x > x2:
    #                 cur_x = x1
    #                 cur_y += 1
    #                 if cur_y > y2:
    #                     return rst
    #             break
    #         elif x1<line_x<=x2:
    #             rst.extend([0]*(line_x-x1))
    #             rst.append(line_v)
    #             cur_x = line_x + 1
    #             cur_y = line_y
    #             if cur_x > x2:
    #                 cur_x = x1
    #                 cur_y += 1
    #                 if cur_y > y2:
    #                     return rst
    #             break
    #         elif line_x > x2:
    #             rst.extend([0]*(x2-x1+1))
    #             cur_x = x1
    #             cur_y = line_y + 1
    #             if cur_y > y2:
    #                 return rst
    #             break
    #
    #
    #
    #     elif line_y > y1 and line_y < y2:
    #         pass
    #     elif line_y == y2:
    #         pass
    #     elif line_y > y2:
    #         pass

    for line in f:
        line_x,line_y,line_v = parse_line(line)[:]
        if line_x > x2 and line_y > y2:
            break

        if line_x == cur_x and line_y == cur_y:
            rst.append(line_v)
            cur_x += 1
            if cur_x > x2:
                cur_x = x1
                cur_y +=1
                if cur_y > y2:
                    return rst
        elif cur_x < line_x <= x2 and line_y == cur_y:
            rst.extend([0]*(line_x-cur_x))
            rst.append(line_v)
            cur_x = line_x + 1
            if cur_x > x2:
                cur_x = x1
                cur_y += 1
                if cur_y > y2:
                    return rst
        elif line_x > x2 and line_y == cur_y:

            rst.extend([0]*(x2-cur_x+1))
            cur_x = x1
            cur_y += 1
            if cur_y > y2:
                return rst

        elif line_y > cur_y:
            if line_y > y2:
                rst.extend([0]*(x2-cur_x+1 + (y2-cur_y)*(x2-x1+1)))
                return rst
            elif line_y <= y2:
                if line_x >= x1 and line_x <= x2:
                    rst.extend([0]*(x2-cur_x+1 + (line_y-cur_y-1)*(x2-x1+1) + line_x-x1))
                    rst.append(line_v)
                    cur_x = line_x + 1
                    cur_y = line_y
                    if cur_x > x2:
                        cur_x = x1
                        cur_y += 1
                        if cur_y > y2:
                            return rst
                elif line_x < x1:
                    rst.extend([0]*(x2-cur_x+1 + (line_y-cur_y-1)*(x2-x1+1)))
                    cur_x = x1
                    cur_y = line_y
                elif line_x > x2:
                    rst.extend([0]*(x2-cur_x+1 + (line_y-cur_y)*(x2-x1+1)))
                    cur_x = x1
                    cur_y = line_y + 1
                    if cur_y > y2:
                        return rst
    return rst



def diagnal_cross_condition(x1,y1,x2,y2,f):
    rst = []

    # cur_x and cur_y point at the position to be given value
    cur_x = x1
    cur_y = y1

    folded_part = get_folded_part_of_square(x1,y1,x2,y2)


    for line in f:
        line_x,line_y,line_v = parse_line(line)[:]
        if line_x > x2 and line_y > y2:
            break







        if line_x == cur_x and line_y == cur_y:
            rst.append(line_v)


            if (line_y,line_x) in folded_part:
                to_put_index = get_array_index(line_y,line_x,x1,y1,x2,y2)
                rst[to_put_index] = line_v


            cur_x += 1
            if cur_x > x2:
                cur_x = x1
                cur_y +=1
                if cur_y > y2:
                    return rst
        elif cur_x < line_x <= x2 and line_y == cur_y:
            rst.extend([0]*(line_x-cur_x))
            rst.append(line_v)

            if (line_y,line_x) in folded_part:
                to_put_index = get_array_index(line_y,line_x,x1,y1,x2,y2)
                rst[to_put_index] = line_v

            cur_x = line_x + 1
            if cur_x > x2:
                cur_x = x1
                cur_y += 1
                if cur_y > y2:
                    return rst
        elif line_x > x2 and line_y == cur_y:

            rst.extend([0]*(x2-cur_x+1))
            cur_x = x1
            cur_y += 1
            if cur_y > y2:
                return rst

        elif line_y > cur_y:
            if line_y > y2:
                rst.extend([0]*(x2-cur_x+1 + (y2-cur_y)*(x2-x1+1)))
                return rst
            elif line_y <= y2:
                if line_x >= x1 and line_x <= x2:
                    rst.extend([0]*(x2-cur_x+1 + (line_y-cur_y-1)*(x2-x1+1) + line_x-x1))
                    rst.append(line_v)

                    if (line_y,line_x) in folded_part:
                        to_put_index = get_array_index(line_y,line_x,x1,y1,x2,y2)
                        rst[to_put_index] = line_v


                    cur_x = line_x + 1
                    cur_y = line_y
                    if cur_x > x2:
                        cur_x = x1
                        cur_y += 1
                        if cur_y > y2:
                            return rst
                elif line_x < x1:
                    rst.extend([0]*(x2-cur_x+1 + (line_y-cur_y-1)*(x2-x1+1)))
                    cur_x = x1
                    cur_y = line_y
                elif line_x > x2:
                    rst.extend([0]*(x2-cur_x+1 + (line_y-cur_y)*(x2-x1+1)))
                    cur_x = x1
                    cur_y = line_y + 1
                    if cur_y > y2:
                        return rst



        # check if next iterator line_x>line_y
        if line_x > line_y:
            # put cur_x and cur_y to next lie
            # and take points that line_x > line_y placeholder
            rst.extend([0]*(x2-cur_y))
            cur_x = x1
            cur_y += 1
            if cur_y > y2:
                return rst
    return rst




def test_all_condition(x1,y1,x2,y2,nth_test):



    bin_sz_li = app.config['test_matrix']





    upload_dir = os.path.join(app.config['APP_ROOT'],'upload')


    if x1>x2 or y1 > y2:
        return 'error 0'
    if (x2-x1-y2+y1)!=0:
        # ensure the block is square
        return 'error 1'

    rst = []

    f = open(os.path.join(upload_dir,'chr1_{name}_norm.txt'.format(name=bin_sz_li[nth_test][0])))
    # remember to close f
    if x2 <= y1:
        rst = easiest_condition(x1,y1,x2,y2,f)
        # easiest condition: the block is above the diagonal line
        # with open(os.path.join(upload_dir,'chr1_{name}_norm.txt'.format(name=bin_sz_li[nth_test][0]))) as f:
        #     for line in f:


    elif y2 <= x1:
        # just opposite the easiest condition
        rst = easiest_condition(y1,x1,y2,x2,f)
    elif (x2 > y1) and (x1 < y1):
        rst = diagnal_cross_condition(x1,y1,x2,y2,f)
        # the block across the diagonal line, but mainly above the diagonal line
    elif (x1 < y2) and (x2 > y2):
        rst = diagnal_cross_condition(y1,x1,y2,x2,f)
        # opposite the former condition
    else:
        #x1==y1 and x2==y2
        rst = diagnal_cross_condition(x1,y1,x2,y2,f)

    f.close()
    return rst





# nth_test = 2
# x1 = 200
# y1 = 250
# x2 = 300
# y2 = 350
#
# a = test_all_condition(x1,y1,x2,y2,nth_test)
# print len(a)










        # print app.config['APP_ROOT']
    # with open('/upload')
    #     upload_dir = os.path.join(app.config['APP_ROOT'],'upload')
    #     to_send = ''
    #     max_val = 0
    #     min_val = 1000
    #     val_li = []
    #     bin_sz_li = [
	 #        ('2.5MB',99),
	 #        ('1MB',249),
	 #        ('500KB',498),
	 #        ('250KB',996),
	 #        ('100KB',2492),
	 #        ('50KB',4984),
	 #        ('25KB',9969),
	 #        ('10KB',24922)
    #         ]
    #     # with open(os.path.join(upload_dir,'chr1_500KB_norm.txt')) as f:
    #     # with open(os.path.join(upload_dir,'chr1_2.5MB_norm.txt')) as f:
    #
    #     nth_test = 2
    #
    #     with open(os.path.join(upload_dir,'chr1_{name}_norm.txt'.format(name=bin_sz_li[nth_test][0]))) as f:
    #         current_line_str = f.next()
    #         for j in range(bin_sz_li[nth_test][1]):
    #             for i in range(0,j+1):
    #                 # current_line_str = f.next()
    #                 current_line_li = current_line_str.split('\t')
    #                 if j!= int(current_line_li[1]) or i!=int(current_line_li[0]):
    #                     to_send += '0'
    #                 else:
    #                     to_send += str(int(float(current_line_li[2][:-1])))
    #                     val_li.append(int(float(current_line_li[2][:-1])))
    #                     max_val = max(int(float(current_line_li[2][:-1])),max_val)
    #                     min_val = min(int(float(current_line_li[2][:-1])),min_val)
    #                     print current_line_str
    #                     current_line_str = f.next()
    #                 to_send += '&'
    #         to_send += '$'
    #         to_send += str(max_val)
    #         to_send += '&'
    #         to_send += str(min_val)
    #         to_send += '&'
    #         val_li.sort()
    #         x_percent_point = int(len(val_li)*0.80)
    #         x_percent_value = val_li[x_percent_point]
    #         to_send += str(x_percent_value)
    #         to_send += '$'
    #         to_send += str(bin_sz_li[nth_test][1])
    #     return to_send
