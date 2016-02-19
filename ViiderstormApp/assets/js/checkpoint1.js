function init(){
    
}

var scene = new THREE.Scene();
var renderer = new THREE.WebGLRenderer({ antialias: true });
var img = document.getElementById("render_image");
renderer.setSize( img.clientWidth, img.clientHeight );
document.getElementById("WebGLCanvas").appendChild(renderer.domElement);

var material1 = new THREE.MeshPhongMaterial({color:0x2194ce});
var material2 = new THREE.MeshPhongMaterial({color:0x4FF5ff});
var material3 = new THREE.MeshPhongMaterial({color:0xff0055});

//camera
var camera = new THREE.PerspectiveCamera( 75, 1, 0.1, 1000);

//light
var light = new THREE.PointLight(0xffffff, 1, 100);
light.position.set(1,1,4);

//plane
var planeGeometry = new THREE.PlaneGeometry(4,8);
var plane = new THREE.Mesh(planeGeometry, material1);
plane.rotation.x += 4.8;

//sphere 1
var sphere1Geometry = new THREE.SphereGeometry(0.5,50,50);
var sphere1 = new THREE.Mesh(sphere1Geometry, material2);

//sphere 2
var sphere2Geometry = new THREE.SphereGeometry(0.45,50,50);
var sphere2 = new THREE.Mesh(sphere2Geometry, material3);

scene.add(plane);
scene.add(sphere1);
scene.add(sphere2);
scene.add(light);

sphere1.position.set(0,0.1,2);
sphere2.position.set(-0.75,-0.3,1);
plane.position.set(-0.5,-1,0);
camera.position.x = 0;
camera.position.z = 4;
camera.position.y = 0;

$('#left-rotate-btn').click(function(){
    if(scene != undefined){
        scene.rotation.y += 0.1;    
    }
})

$('#right-rotate-btn').click(function(){
    if(scene != undefined){
        scene.rotation.y += -0.1;    
    }
})

$('#reset-rotate-btn').click(function(){
    if(scene != undefined){
        scene.rotation.y = 0;    
    }
})

var render = function(){
    
    requestAnimationFrame(render);
    
    /* DEBUG
    
    sphere1.position.set(0.5,0.1,2);
    sphere2.position.set(-0.25,-0.3,1);
    plane.position.set(0,-1,0);

    camera.position.x = 0.5;
    camera.position.z = 4;
    camera.position.y = 0;
    */
    
    renderer.render(scene, camera);
}

render();

window.addEventListener('resize', function(){
    renderer.setSize(img.clientHeight, img.clientWidth);
    camera.aspect = img.clientWidth/img.clientHeight;
    camera.updateProjectionMatrix();
    camera.position = new THREE.Vector3(0,0,4);
    render();
}, false)