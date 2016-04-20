// Your beautiful D3 code will go here
            // var w = 500;
            // var h = 50;
            // var barPadding = 1;
            // var dataset = [5,10,15,20,25]



// var canvas = d3.select('body').append('canvas');
// canvas.attr('width')

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
                    	var canvas = d3.select('body').append('canvas');
                    	canvas.attr('width',Drawer.w)
                    			.attr('height',Drawer.h)
                    			.attr('id',Drawer.canvasId)
                    			.style('border','1px solid #000000');
                    	drawer.canvas = canvas;









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


    $.post('/matrix/data',
        {},
    function(data,status){


            var canvas_width = 1000;
            var canvas_height = 1000;
            var canvas_nrows = 997;

            var d = Drawer.createNew(canvas_width,canvas_height,canvas_nrows);
            d.makeCanvas();
            d.beforeDrawElement();

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






    })