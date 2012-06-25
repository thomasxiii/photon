var $body,
	$demos,
	$crane,
	$map,
	$mapPanel1,
	$mapPanel2,
	$mapPanel3,
	$mapCover,
	$coverflow,
	$toggleBtn,
	$toggleOn,
	$toggleOff,
	crane,
	craneFaces,
	cubeFaces,
	map,
	diamondFaces,
	coverflowFaces,
	shadeAmount,
	tintAmount,
	light,
	currentCover,
	renderTimer,
	isLit,
	domTransformProperty,
	cssTransformProperty,
	domTransitionProperty,
	cssTransitionProperty,
	transitionEndEvent,
	transitionEndEvents = {
		'WebkitTransition' : 'webkitTransitionEnd',
		'MozTransition'    : 'transitionend',
		'OTransition'      : 'oTransitionEnd',
		'msTransition'     : 'MSTransitionEnd',
		'transition'       : 'transitionend'
	};





$(document).ready(function() {
	$body = $('body');
	$demos = $('.demo');
	light = new Photon.Light();
	shadeAmount = .5;
	tintAmount = 0;
	coverflowFaces = [];
	cubeFaces = [];
	diamondFaces = [];
	currentCover = 0;
	renderCurrent = renderCrane;
	$toggleBtn = $('.toggle-btn');
	$toggleOn = $('.toggle .label-on');
	$toggleOff = $('.toggle .label-off');
	isLit = true;
	domTransformProperty = Modernizr.prefixed('transform');
	cssTransformProperty = domToCss(domTransformProperty);
	domTransitionProperty = Modernizr.prefixed('transition');
	cssTransitionProperty = domToCss(domTransitionProperty);
	transitionEndEvent = transitionEndEvents[domTransitionProperty];

	setupLightControls();
	setupCoverflow();
	setupCrane();
	setupMap();

	// demo menu
	$('.example-menu a').bind('click', onDemoNav);

	if(cssTransformProperty == '-webkit-transform') {
		showCrane();
	} else {
		$('.map-thumb').click();
		showMap();
		$('.crane').hide();
		$('.crane-thumb').hide();
	}
});







/*---------------------------------

	Light Controls

---------------------------------*/

function setupLightControls() {
	$('.toggle a').bind('click', toggleLight);
}

function toggleLight(e) {
	e.preventDefault();

	switch($(e.target).attr('id')) {
		case 'label-on':
			isLit = true;
			$toggleBtn.addClass('on');
			$toggleOn.addClass('current');
			$toggleOff.removeClass('current');
			$('.photon-shader').show();
			break;
		case 'label-off':
			isLit = false;
			$toggleBtn.removeClass('on');
			$toggleOn.removeClass('current');
			$toggleOff.addClass('current');
			$('.photon-shader').hide();
			break;
		case 'toggle-btn':
			isLit = !isLit;
			$toggleBtn.toggleClass('on');
			$toggleOn.toggleClass('current');
			$toggleOff.toggleClass('current');
			$('.photon-shader').toggle();
			break;
	}
}








/*---------------------------------

	Menus

---------------------------------*/

function onDemoNav(e) {
	e.preventDefault();

	var demo = $(e.target).attr('data-demo');

	$('.example-menu .current').removeClass('current');
	$(this).addClass('current');

	switch(demo) {
		case 'coverflow':
			hideCrane();
			hideMap();
			showCoverflow();
			renderCurrent = renderCoverflow;
			break;
		case 'crane':
			hideCoverflow();
			hideMap();
			showCrane();
			renderCurrent = renderCrane;
			break;
		case 'map':
			hideCoverflow();
			hideCrane();
			showMap();
			renderCurrent = renderMap;
			break;
	}

	renderCurrent();
	if(!isLit) {
		$('.photon-shader').hide();
	}
}









/*---------------------------------

	Crane

---------------------------------*/

function setupCrane() {
	$crane = $('.crane');
	crane = new Photon.FaceGroup($('.crane')[0], $('.crane .face'), .6, .1, true);
	renderCrane();
}

function renderCrane() {
	crane.render(light, true);
}

function showCrane() {
	$body.bind('mousemove', rotateCrane);
	$crane.show();
}

function hideCrane() {
	$body.unbind('mousemove', rotateCrane);
	$crane.hide();
}

function rotateCrane(e) {
	var xPer = e.pageX / $body.width();

	$(crane.element).css(cssTransformProperty, 'rotateX(-15deg) rotateY(' + (-180 + (xPer * 360)) + 'deg)');
	renderCrane();
}









/*---------------------------------

	Map

---------------------------------*/

function setupMap() {
	$map = $('.map');
	$mapPanel1 = $('.panel-1');
	$mapPanel2 = $('.panel-2');
	$mapPanel3 = $('.panel-3');
	$mapCover = $('.map-cover');
	$map.bind('click', toggleMap);

	map = new Photon.FaceGroup($('.map')[0], $('.map .face'), 1.5, .2, true);
	renderMap();
}

function toggleMap() {
	$map.toggleClass('is-open');

	$map.unbind();
	$map.bind(transitionEndEvent, stopRenderTimer);

	if(!renderTimer) {
		renderTimer = setInterval(renderMap, 34);
	}
}

function renderMap() {
	map.render(light, true, true);
}

function showMap() {
	$body.bind('mousemove', rotateMap);
	$map.show();
}

function hideMap() {
	$body.unbind('mousemove', rotateMap);
	$map.hide();
}

function rotateMap(e) {
	var xPer = e.pageX / $body.width();
	var yPer = e.pageY / $body.height();

	$mapPanel1.css(cssTransformProperty, 'rotateY(' + (178 - (138 * xPer)) + 'deg)');
	$mapCover.css(cssTransformProperty, 'rotateY(' + (178 - (138 * xPer)) + 'deg) translateZ(-2px) rotateY(180deg) translateX(240px)');
	$mapPanel3.css(cssTransformProperty, 'rotateY(' + (178 - (138 * xPer)) + 'deg)');
	$map.css(cssTransformProperty, 'rotateX(' + (40 - (yPer * 70)) + 'deg) rotateY(' + (20 - (xPer * 60)) + 'deg) rotateZ(0)');

	renderMap();
}









/*---------------------------------

	Coverflow

---------------------------------*/

function setupCoverflow() {
	$coverflow = $('.coverflow');
	var $coverflowItems = $coverflow.find('li');

	$coverflowItems.each(function(i) {
		coverflowFaces[i] = new Photon.Face($(this)[0], shadeAmount);
	});

	console.log(transitionEndEvent);
	$coverflowItems.eq(1).bind(transitionEndEvent, stopRenderTimer);
	
	setCoverTransforms();
}

function changeCover() {
	currentCover = currentCover < coverflowFaces.length - 1 ? currentCover + 1 : 0;
	setCoverTransforms(true);
}

function setCoverTransforms(animate) {
	if(!renderTimer && animate) {
		renderTimer = setInterval(renderCoverflow, 34);
	}
	for(var i = 0; i < coverflowFaces.length; i++) {
		var element = coverflowFaces[i].element;
		var offset = Math.abs(currentCover - i);
		var x = i == currentCover ? 0 : (150 + (100 * offset)) * (i < currentCover ? -1 : 1);
		var z = i == currentCover ? 0 : -200;

		var rotationY = i == currentCover ? 0 : (80 + (offset * -5)) * (i < currentCover ? 1 : -1);

		$(element).css(cssTransformProperty, 'translateX(' + x +'px) translateZ(' + z + 'px) rotateY(' + rotationY + 'deg)');
	}
}

function rotateCoverflow(e) {
	var xPer = e.pageX / $body.width();

	var newIndex = (coverflowFaces.length -1) - Math.round((coverflowFaces.length -1) * xPer);

	if(!renderTimer && newIndex != currentCover) {
		renderTimer = setInterval(renderCoverflow, 34);
		currentCover = newIndex;
	}
	for(var i = 0; i < coverflowFaces.length; i++) {
		var element = coverflowFaces[i].element;
		var offset = Math.abs(currentCover - i);
		var x = i == currentCover ? 0 : (150 + (100 * offset)) * (i < currentCover ? -1 : 1);
		var z = i == currentCover ? 0 : -200;

		var rotationY = i == currentCover ? 0 : (80 + (offset * -5)) * (i < currentCover ? 1 : -1);

		$(element).css(cssTransformProperty, 'translateX(' + x +'px) translateZ(' + z + 'px) rotateY(' + rotationY + 'deg)');
	}
}

function stopRenderTimer() {
	if(renderTimer) {
		clearInterval(renderTimer);
		renderTimer = null;
	}
}

function renderCoverflow() {	
	for(var i = 0; i < coverflowFaces.length; i++) {
		coverflowFaces[i].render(light, true);
	}
}

function hideCoverflow() {
	$coverflow.hide();
	$body.unbind();
}

function showCoverflow() {
	$coverflow.show();
	$body.bind('mousemove', rotateCoverflow);
}









/*---------------------------------

	Utilities

---------------------------------*/

function domToCss(property) {
	var css = property.replace(/([A-Z])/g, function (str, m1) {
		return '-' + m1.toLowerCase();
	}).replace(/^ms-/,'-ms-');

	return css;
}

function clamp(val, min, max) {
    if(val > max) return max;
    if(val < min) return min;
    return val;
}
