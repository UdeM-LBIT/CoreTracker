$(document).ready(function() {
		$('#showadvanced').click(function() {
			$('#advanced-drop').slideToggle("fast");
		});

		var opts = {
			lines: 10,
			length: 5,
			width: 5 // The line thickness
				,
			radius: 10,
			corners: 1 // Corner roundness (0..1)
				,
			color: '#666' // #rgb or #rrggbb or array of colors
				,
			opacity: 0.25 // Opacity of the lines
				,
			trail: 60 // Afterglow percentage
				,
			className: 'spinner' // The CSS class to assign to the spinner
				,
			top: '50%' // Top position relative to parent
				,
			left: '50%' // Left position relative to parent
				,
			shadow: false // Whether to render a shadow
		}
		var spinner = new Spinner(opts);
		var spinspan = document.getElementById('spin')


		var filelist = {};
		$('input[type=file]').on('change', prepareUpload);
		var error_on_file = {}
		function prepareUpload(event) {
			var file = event.target.files[0];
			if(file.size/1000000 > 10){
				error_on_file[event.target.id] = false
				$('#show-error').text("File "+file.name+" should be less than 10MB").slideToggle().delay(8000).slideToggle();
				$(event.target).parent().toggleClass("ui-state-error", true);

			}
			else{
				filelist[event.target.id] = file;
				error_on_file[event.target.id] = file.name
				$(event.target).parent().toggleClass("ui-state-error", false);
  				
			}
			
		}

		$('form').on('submit', function(event){
			event.preventDefault();
			var msg = []
			var send = true;
			size = 0 
			for (key in error_on_file){
				if (error_on_file.hasOwnProperty(key)) size++;
				if (error_on_file[key] == false){
					send = false;
					msg[0] = "file larger than 10MB";
				}
			}

			if(size<3) {
				send = false;
				msg[0] = "required parameters not found"
			}

			valist = ['aamajthresh', 'freqthresh', 'idfilter', 'icfilter', 'gapfilter']
			for(len=valist.length, i=0; i<len; i++){
				val = valist[i]
				var val_id = "#"+val;
				$(val_id).toggleClass("ui-state-error", false);
				var valE = parseInt($(val_id).val());
				if (!(valE>0 && valE/100.0 <1)){
					send = false;

					msg[1]= "Invalid value in parameters"
					$(val_id).toggleClass("ui-state-error", true);
					setTimeout(function () {
						$(val_id).toggleClass("ui-state-error", false);
					}, 10000);
				}
			}

			
			if (send == true) {
				uploadFiles(event);	
			}
			else{
				$('#show-error').text("Could not submit, error in your form. Reasons : " + msg.join(', ')).slideToggle().delay(10000).slideToggle();
			}
			
		});

			

		function uploadFiles(event) {
			event.stopPropagation();
			spinner.spin(spinspan);
			// Create a formdata object and add the files
			var data = new FormData();

			$.each(filelist, function(key, value) {
				data.append(key, value);
			});

			$.ajax({
					url: ete_webplugin_URL + '/uploadfile',
					type: 'POST',
					data: data,
					cache: false,
					processData: false, // Don't process the files
					contentType: false, // Set content type to false as jQuery will tell the server its a query string request

					success: function(data, textStatus, jqXHR) {
						if (typeof data.error === 'undefined') {
							
							submitForm(event, data);
						
						} else {
							$('#show-error').text(data.error).slideToggle().delay(10000).slideToggle();
							spinner.stop();

						}
					},
					error: function(jqXHR, textStatus, errorThrown) {
						$('#show-error').text(textStatus).slideToggle().delay(10000).slideToggle();
						spinner.stop();
					}
		});
}

function submitForm(event, data) {
	// Create a jQuery object from the form
	$form = $(event.target);

	// Serialize the form data
	var formData = $form.serialize();


	// You should sterilise the file names
	$.each(data, function(key, value) {
		formData = formData + '&'+key+'=' + value;
	});

	$.ajax({
		url: ete_webplugin_URL + '/exectracker',
		type: 'GET',
		data: formData,
		dataType: 'json',
		success: function(data, textStatus, jqXHR) {
			if (typeof data.error === 'undefined') {
				// Success so call function to process the form
				console.log('SUCCESS: ' + data.success);
				// If this is a success,
				// add a link to where we should see the rest of the data and 
				// instruction to how to refresh every x second
			} else {
				// Handle errors here
				console.log('ERRORS: ' + data.error);
			}
		},
		error: function(jqXHR, textStatus, errorThrown) {
			// Handle errors here
			$('#show-error').text(textStatus).slideToggle().delay(10000).slideToggle();
		},
		complete: function() {
			spinner.stop();
		}
	});
}

});