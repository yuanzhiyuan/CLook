// Your beautiful D3 code will go here
            // var w = 500;
            // var h = 50;
            // var barPadding = 1;
            // var dataset = [5,10,15,20,25]



// var canvas = d3.select('body').append('canvas');
// canvas.attr('width')


d3.select('#drag_container')
    .style('width','200px')
    .style('height','200px')
    .style('overflow','hidden')
    .style('border','1px solid');


var Drawer = {
                w: 500,//default width of canvas
                h: 500,//default height of canvas
                s: 499,//default row/col number of canvas
                b: 1, //default canvas border width
                canvasId: 'myCanvas',
                borderColor: 'black',
                createNew: function(width,height,nrows){
                    Drawer.w = width;
                    Drawer.h = height;
                    Drawer.s = nrows;
                    var drawer = {};
                    drawer.makeCanvas = function(){
                    	var canvas = d3.select('#drag_container').append('canvas');
                    	canvas.attr('width',Drawer.w+'px')
                    			.attr('height',Drawer.h+'px')
                    			.attr('id',Drawer.canvasId)

                    	drawer.canvas = canvas;









                    };
                    drawer.setXPercentValue = function(x){

                      drawer.x_percent_value = x;

                    };
                    drawer.beforeDrawElement = function(){
                    	drawer.each_cell_size = Drawer.w/Drawer.s;
                    	// var x_pos = x * each_cell_size;
                    	// var y_pos = y * each_cell_size;
                    	// var ctx = document.getElementById('myCanvas').getContext('2d');
                    	drawer.ctx = drawer.canvas.node().getContext('2d');

                    }
                    drawer.drawElement = function(x,y,value_r,value_g,value_b){

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

    $.post('/matrix/data',
        {
//            type:'init_data'   //获取matrix的row数，九分位点等。
        },
    function(data,status){



        var canvas_nrows = parseInt(data.split('$')[2]);
//        var x_percent_value = parseInt(data.split('&')[1]);
        var each_cell_sz = 1;
        var canvas_width = canvas_nrows * each_cell_sz;
        var canvas_height = canvas_nrows * each_cell_sz;

        d = Drawer.createNew(canvas_width,canvas_height,canvas_nrows);
        d.makeCanvas();
        d.beforeDrawElement();
        d.setXPercentValue(x_percent_value);
//        console.log(canvas_nrows, d.x_percent_value);

        var k=0;
        var data_array = data.split('$')[0].split('&');
        var data_max_min_array = data.split('$')[1].split('&');

        var max_data = parseInt(data_max_min_array[0]);
        var min_data = parseInt(data_max_min_array[1]);
        var x_percent_value = parseInt(data_max_min_array[2]);
        var red_scale = d3.scale.linear()
//                                .domain([min_data,max_data])
                            .domain([0,x_percent_value])
                            .range([0,255]);


        for(var j=0;j<canvas_nrows;j++){
            for(var i=0;i<=j;i++){

                        d.drawElement(i,j,255,255-parseInt(red_scale(data_array[k])),255-parseInt(red_scale(data_array[k])));
                        d.drawElement(j,i,255,255-parseInt(red_scale(data_array[k])),255-parseInt(red_scale(data_array[k])));
                    k++;
                // console.log(red_scale(i*j))
            }
        }

        $('#myCanvas').draggable({


            containment: [200-canvas_width,200-canvas_height,0,0]
        });

        handleDrawer(d);
    })
}


var get_block_data = function(x1,y1,x2,y2){
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
        });
}


init_canvas(function(drawer_obj){
//    alert('aaa');



});





