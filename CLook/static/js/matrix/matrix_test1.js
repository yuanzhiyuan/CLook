// Your beautiful D3 code will go here
            // var w = 500;
            // var h = 50;
            // var barPadding = 1;
            // var dataset = [5,10,15,20,25]



// var canvas = d3.select('body').append('canvas');
// canvas.attr('width')

div_sz = 200;

d3.select('#drag_container')
    .style('width',div_sz+'px')
    .style('height',div_sz+'px')
    .style('overflow','hidden')
//    .style('border','1px solid');


var Drawer = {
                w: 500,//default width of canvas
                h: 500,//default height of canvas
                s: 499,//default row/col number of canvas
                b: 1, //default canvas border width
                canvasId: 'myCanvas',
                borderColor: 'black',
                createNew: function(width,height,nrows){
                    var drawer = {};


                    drawer.canvas_bitmap = new Array((nrows)*(nrows)).fill(0);

                    Drawer.w = width;
                    Drawer.h = height;
                    Drawer.s = nrows;
//                    drawer.w = width;
//                    drawer.h = height;
//                    drawer.s = nrows;

                    drawer.makeCanvas = function(){
                    	var canvas = d3.select('#drag_container').append('canvas');
                    	canvas.attr('width',Drawer.w+'px')
                    			.attr('height',Drawer.h+'px')
                    			.attr('id',Drawer.canvasId)

                    	drawer.canvas = canvas;









                    };
                    drawer.setXPercentValue = function(x){

                      drawer.x_percent_value = parseInt(x);

                    };
                    drawer.beforeDrawElement = function(){
                    	drawer.each_cell_size = Drawer.w/Drawer.s;
                    	// var x_pos = x * each_cell_size;
                    	// var y_pos = y * each_cell_size;
                    	// var ctx = document.getElementById('myCanvas').getContext('2d');
                    	drawer.ctx = drawer.canvas.node().getContext('2d');

                    }
                    drawer.drawElement = function(y,x,value_r,value_g,value_b){
                        //这里的水平位x，纵位y。与matrix的x，y刚好相反。我们用这个函数的时候，就是用在matrix上，

                    	var x_pos = x * drawer.each_cell_size;
                    	var y_pos = y * drawer.each_cell_size;
                    	// var ctx = document.getElementById('myCanvas').getContext('2d');


                    	drawer.ctx.fillStyle = "rgb("+value_r+","+value_g+","+value_b+")";
                    	drawer.ctx.fillRect(x_pos,y_pos,drawer.each_cell_size,drawer.each_cell_size);
                    	// drawer.ctx.fillRect(x_pos+each_cell_size,y_pos,each_cell_size,each_cell_size);
                    }
                    return drawer;
                }




            }




var init_canvas = function(handleDrawer){

    $.post('/matrix/data/test1',
        {
            type:'init_data'   //获取matrix的row数，九分位点等。
        },
    function(data,status){



        var canvas_nrows = parseInt(data.split('&')[0]);
        var x_percent_value = parseInt(data.split('&')[1]);
        var each_cell_sz = 1;
        console.log(each_cell_sz);
        var canvas_width = canvas_nrows * each_cell_sz;
        var canvas_height = canvas_nrows * each_cell_sz;
//        console.log(canvas_nrows);
//        console.log(x_percent_value);
        var d = Drawer.createNew(canvas_width,canvas_height,canvas_nrows);
        d.makeCanvas();
        d.beforeDrawElement();
        d.setXPercentValue(x_percent_value);
//        console.log(canvas_nrows, d.x_percent_value);

//        var k=0;
//        var data_array = data.split('$')[0].split('&');
//        var data_max_min_array = data.split('$')[1].split('&');

//        var max_data = parseInt(data_max_min_array[0]);
//        var min_data = parseInt(data_max_min_array[1]);
//        var x_percent_value = parseInt(data_max_min_array[2]);
//        var red_scale = d3.scale.linear()
////                                .domain([min_data,max_data])
//                            .domain([0,x_percent_value])
//                            .range([0,255]);
//
//
//        for(var j=0;j<canvas_nrows;j++){
//            for(var i=0;i<=j;i++){
//
//                        d.drawElement(i,j,255,255-parseInt(red_scale(data_array[k])),255-parseInt(red_scale(data_array[k])));
//                        d.drawElement(j,i,255,255-parseInt(red_scale(data_array[k])),255-parseInt(red_scale(data_array[k])));
//                    k++;
//                // console.log(red_scale(i*j))
//            }
//        }

        $('#myCanvas').draggable({


            containment: [div_sz-canvas_width+8,div_sz-canvas_height+8,8,8],
            // 为什么要加8？ 这里的containment相对应的是整个页面的位置，而不是parent的位置。而canvas的parent（div）默认在页面的（8,8）位置。

            stop:function(event,ui){
//                console.log(d.canvas_bitmap[d.canvas_bitmap.length-1]);

//                获取当前窗口区域在整个canvas的坐标
                var cur_win_x1 = ui.offset.top * (-1)/(Drawer.w/Drawer.s);
                var cur_win_y1 = ui.offset.left * (-1)/(Drawer.w/Drawer.s);
                var cur_win_x2 = cur_win_x1 + div_sz/(Drawer.w/Drawer.s);
                var cur_win_y2 = cur_win_y1 + div_sz/(Drawer.w/Drawer.s);
//                console.log(cur_win_x1,cur_win_y1);
//                console.log(cur_win_x2,cur_win_y2);
                var scale_ = d3.scale.linear()
////                                .domain([min_data,max_data])
                            .domain([0,d.x_percent_value])
                            .range([0,255]);
                var red_scale = function(val,threshold,scaler){
                    if(val >= threshold)    return 255;
                    else    return scale_(val);
                }

                get_block_data(cur_win_x1,cur_win_y1,cur_win_x2,cur_win_y2,function(data){
                    var data_array = data.split('&');
                    var k = 0;
                    for(var j = cur_win_y1;j<=cur_win_y2;j++){
                        for(var i=cur_win_x1;i<=cur_win_x2;i++){
//                            console.log(d.canvas_bitmap[get_array_index(i,j,Drawer.s+1)]);
//                            if(d.canvas_bitmap[get_array_index(i,j, Drawer.s+1)]==0){
//                            console.log(red_scale(parseInt(data_array[k]), d.x_percent_value,scale_))
                            console.log(i,j);

                                d.drawElement(i,j,255,255-parseInt(red_scale(parseInt(data_array[k]), d.x_percent_value,scale_)),255-parseInt(red_scale(parseInt(data_array[k]), d.x_percent_value,scale_)));
                                d.canvas_bitmap[get_array_index(i,j, Drawer.s+1)]=1;
//                                console.log(i,j);
//                            }

                            k++;
                        }
                    }
                });
            }
        });

        handleDrawer(d);
    })
}


var get_block_data = function(x1,y1,x2,y2,handle_data){
    //获取某正方形中需要的数据
    $.post('/matrix/data/get',
        {
            x1:x1,
            y1:y1,
            x2:x2,
            y2:y2
        },
        function(data,status){
           //希望是排好序的data，减轻浏览器负担
            handle_data(data);
        });
}


init_canvas(function(drawer_obj){
//    alert('aaa');

    var cur_win_x1 = 0;
    var cur_win_y1 = 0;
    var cur_win_x2 = div_sz/(Drawer.w/Drawer.s);
    var cur_win_y2 = div_sz/(Drawer.w/Drawer.s);
    var scale_ = d3.scale.linear()
////                                .domain([min_data,max_data])
                            .domain([0,drawer_obj.x_percent_value])
                            .range([0,255]);
    var red_scale = function(val,threshold,scaler){
        if(val >= threshold)    return 255;
        else    return scale_(val);
    }
    get_block_data(cur_win_x1,cur_win_y1,cur_win_x2,cur_win_y2,function(data){
        var data_array = data.split('&');
        var k = 0;
        for(var j = cur_win_y1;j<=cur_win_y2;j++){
            for(var i=cur_win_x1;i<=cur_win_x2;i++){
                drawer_obj.canvas_bitmap[get_array_index(i,j,Drawer.s+1)] = 1;
                drawer_obj.drawElement(i,j,255,255-parseInt(red_scale(parseInt(data_array[k]),drawer_obj.x_percent_value,scale_)),255-parseInt(red_scale(parseInt(data_array[k]),drawer_obj.x_percent_value,scale_)));
                k++;
            }
        }
    });


//console.log(drawer_obj.w,drawer_obj.h,drawer_obj.s);


});


var get_array_index = function(x,y,sz){
    return x+y*sz;
}


//$.post('/matrix/data/get',
//    {
//        x1:100,
//        y1:150,
//        x2:200,
//        y2:250
//    },
//function(data,status){
//    matrix_data = data.split('$')[0].split('&');
//    x_percent_value = data.split('$')[1]
//    console.log(matrix_data)
//})




