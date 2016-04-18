$(function(){
	
  $('#table_parent').children().each(function(){
  	$(this).hide();
  })
});



var view_detail = function(obj){
	
	var category_id = $(obj).attr('categoryid')

	var table_id = 'table' + category_id;

	var table_container = $('#'+table_id).parent();

	table_container.children().each(function(){
		$(this).hide();
	});
	$('#'+table_id).show();

}