/**
 * Created by yuan on 3/1/16.
 */


function add_project(infomation){
    var name = $('#name').val();
    var category = $('#category').val();
    var species = $('#species').val();
    var cell_type = $('#cell_type').val();
    var manager_no = $('#manager').val();

    $.post('/addProject',
        {
            name:name,
            category:category,
            species:species,
            cell_type:cell_type,
            manager_no:manager_no,
            infomation:infomation
        },
        function(data,status){
            if(data=='success'){
                alert('succeed');
                location.reload();
            }else{
                alert('fail');
                location.reload();
            }
        }
    )
}

function add_file(infomation){
    console.log('bbb')
    var name = $('#name').val();
//    var location = $('#location').val();
    var project = $('#project').val();
    var type = $('#type').val();
    var file_type = $('#file_type').val();
    var author = $('#author').val();
//    add to tree location
    var under = $('#under').val();

    $.post('/addFile',
        {
            name:name,
//            location:location,
            project:project,
            type:type,
            file_type:file_type,
            author:author,
            under:under,
            infomation:infomation

        },
        function(data,status){
            if(data=='success'){
                alert('succeed');
                location.reload();
            }else{
                alert('fail');
                location.reload();
            }
        }
    )
}


// 封装console.log
    function printLog(title, info) {
        window.console && console.log(title, info);
    }

    // ------- 配置上传的初始化事件 -------
    function uploadInit () {
        console.log('afas')
        // this 即 editor 对象
        var editor = this;
        // 编辑器中，触发选择图片的按钮的id
        var btnId = editor.customUploadBtnId;
        // 编辑器中，触发选择图片的按钮的父元素的id
        var containerId = editor.customUploadContainerId;

        //实例化一个上传对象
        var uploader = new plupload.Uploader({
            browse_button: btnId,  // 选择文件的按钮的id
            url: '/upload',  // 服务器端的上传地址
            flash_swf_url: '../../plupload/Moxie.swf',
            sliverlight_xap_url: '../../plupload/Moxie.xap',
            filters: {
                mime_types: [
                    //只允许上传图片文件 （注意，extensions中，逗号后面不要加空格）
                    { title: "图片文件", extensions: "jpg,gif,png,bmp" },
                    { title: "html files", extensions: "html" }
                ]
            }
        });

        //存储所有图片的url地址
        var urls = [];

        //初始化
        uploader.init();

        //绑定文件添加到队列的事件
        uploader.bind('FilesAdded', function (uploader, files) {
            //显示添加进来的文件名
            $.each(files, function(key, value){
                printLog('添加文件' + value.name);
            });

            // 文件添加之后，开始执行上传
            uploader.start();
        });

        //单个文件上传之后
        uploader.bind('FileUploaded', function (uploader, file, responseObject) {
            //注意，要从服务器返回图片的url地址，否则上传的图片无法显示在编辑器中
            var url = responseObject.response;
            //先将url地址存储来，待所有图片都上传完了，再统一处理
            urls.push(url);

            printLog('一个图片上传完成，返回的url是' + url);
        });

        //全部文件上传时候
        uploader.bind('UploadComplete', function (uploader, files) {
            printLog('所有图片上传完成');

            // 用 try catch 兼容IE低版本的异常情况
            try {
                //打印出所有图片的url地址
                $.each(urls, function (key, value) {
                    printLog('即将插入图片' + value);

                    // 插入到编辑器中
                    editor.command(null, 'insertHtml', '<img src="' + value + '" style="max-width:100%;"/>');
                });
            } catch (ex) {
                // 此处可不写代码
            } finally {
                //清空url数组
                urls = [];

                // 隐藏进度条
                editor.hideUploadProgress();
            }
        });

        // 上传进度条
        uploader.bind('UploadProgress', function (uploader, file) {
            // 显示进度条
            editor.showUploadProgress(file.percent);
        });
    }

    //這麼做是爲了在調用add_file之前就已經執行完成了/upload的動作
    function sumit_location(infomation){
        $.post('/upload',
            {},
            function(data,status){
                console.log(data);
                add_file(infomation);
            }

        )



    }
$(document).ready(function(){





   var editor = new wangEditor('infomation');
//    editor.config.customUpload = true;  // 配置自定义上传的开关
//    console.log(editor)
//    editor.config.customUploadInit = uploadInit;  // 配置上传事件，uploadInit方法已经在上面定义了
//    console.log('zzz')
    editor.create();
    var infomation = editor.$txt.html();
    $('#btn1').click(function(){
        add_project(infomation)
    });
//    $('#uploadfiles').click(function(){
//        sumit_location(infomation);
//    })
    $('#add_file_submit').click(function(){
        if($('#uploadfiles').length){

//            console.log('aaaaa');
//            $('#uploadfiles').click();
        }
        add_file(infomation);

    });

});




function login(){
//    alert('???')
//    记得写val()
    var username=$('#username').val()
    var password=$('#password').val()
    if(username==""){
        alert('请填写用户名')
        return false
    }
    if(password==""){
        alert('请填写密码')
        return false
    }
    $.post('/login',
        {
            username:username,
            password:password
        },
        function(data,status){

            if(data=='success'){
                alert('登录成功！')
                location.href='/'
            }
            else if(data=='error1'){
                alert('用户名/密码错误')
                location.reload()
            }
            else if(data=='error2'){
                alert('请输入密码')
                location.reload()
            }
            else if(data=='error4'){
                alert('系统错误')
                location.reload()
            }


        })
}


function addUser(){
    var username=$('#username').val()
    var password=$('#password').val()
    var repassword=$('#repassword').val()
    if(username==''){
        alert('请输入用户名')
        return false
    }
    if(password==''){
        alert('请输入新密码')
        return false
    }
    if(repassword==''){
        alert('请确认密码')
        return false
    }
    if(password!=repassword){
        alert('两次密码不一致')
        return false
    }

    $.post('/user/add',
        {
            username:username,
            password:password,
            repassword:repassword
        },
        function(data,status){
            if(data=='success'){
                alert('添加成功！')
                location.reload()
            }
            else{
                alert(data)
                location.reload()
            }
        })

}


function changePassword(){
    var oldpassword=$('#oldpassword').val()
    var newpassword=$('#newpassword').val()
    var renewpassword=$('#renewpassword').val()

    if(oldpassword==''){
        alert('请输入旧密码')
        return false
    }
    if(newpassword==''){
        alert('请输入新密码')
        return false
    }
    if(renewpassword==''){
        alert('请确认密码')
        return false
    }
    $.post('/user/changepwd',
        {
            oldpassword:oldpassword,
            newpassword:newpassword,
            renewpassword:renewpassword
        },
        function (data,status) {
            if(data=='success'){
                alert('修改成功')
                location.href='/login'
            }else{
                alert(data)
                location.reload()
            }

    })
}
