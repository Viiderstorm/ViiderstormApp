function ImageFrame(iwidth, iheight, location){
    var self = this;
    self.imageWidth = iwidth;
    self.imageHeight = iheight;
    
    self.planeWidth = 2;
    self.planeHeight = 2;
    
    self.pixels = (function(){
        var tmp = [];
        for(var i = 0; i < self.imageHeight; i++){
            tmp.push([])
            for(var k = 0; k < self.imageWidth; k++){
                tmp[i].push(new THREE.Color( 0x000000))
            }
        }
        return tmp
    })()
    
    self.pixelWidth  = self.planeWidth / self.imageWidth;
    self.pixelHeight = self.planeHeight / self.imageHeight;
    
    self.generateImageData = function(){
        var tmp = [];
        for(var i = 0; i < self.imageHeight; i++){
            for(var k = 0; k < self.imageWidth; k++){
                
                var style = self.pixels[i][k].getStyle().replace("rgb(", "").replace(")", "")
                var values = style.split(",");
                
                tmp.push(parseInt(values[0]));
                tmp.push(parseInt(values[1]));
                tmp.push(parseInt(values[2]));
                tmp.push(255);
            }
        }
        
        return new Uint8ClampedArray(tmp);
    }
    
    self.displayImage = function(){
        
        var c = document.getElementById("renderCanvas");
        var ctx = c.getContext("2d");
        var imgData = ctx.createImageData(self.imageHeight,self.imageWidth);
        imgData.data.set(self.generateImageData());
        ctx.putImageData(imgData, 0, 0);
        
    }
    
    self.location = location;
}


function RayTracer(camera, frame, scene){
    var self = this;
    
    self.frame = frame;
    self.geometries = (function(){
        var tmp = []
        for(var i = 0; i < scene.children.length; i++){
            if(scene.children[i] instanceof THREE.Mesh){
                tmp.push(scene.children[i]);
            }
        }
        return tmp;
    })()
    
    self.lights = (function(){
        var tmp = []
        for(var i = 0; i < scene.children.length; i++){
            if(scene.children[i] instanceof THREE.PointLight){
                tmp.push(scene.children[i]);
            }
        }
        return tmp;
    })()
    
    self.rays = generateRays(camera, frame);
    
    self.traceRays = function(){
        for(var i = 0; i < self.rays.length; i++){
            for(var k = 0; k < self.rays[i].length; k++){
                
                var intersectedObject = findIntersection(self.rays[i][k], self.geometries);
                if(intersectedObject != null){
                    self.frame.pixels[i][k] = intersectedObject.material.color;
                }
                
            }
        }
        
        return self.frame;
        
    }
    
    self.render = function(){
        var imgPlane = self.traceRays();
        imgPlane.displayImage();
    }
     
}

function Ray(origin, direction){
    var self = this;
    
    self.origin = origin;
    self.direction = direction;
}

function generateRays(camera, frame){
    
    //assume camera is centered behind frame, looking through the center of frame
    var cameraLocation = camera.position.clone();
    var originalFramePixel = frame.location.clone();
    originalFramePixel.x = originalFramePixel.x - (frame.planeWidth/2) + (0.5 * frame.pixelWidth);
    originalFramePixel.y = originalFramePixel.y + (frame.planeHeight/2) - (0.5 * frame.pixelHeight);
    var currentFramePixel = originalFramePixel.clone();
   
    var rays = [];
    for(var i = 0; i < frame.imageHeight; i++){
        rays.push([]);
        for(var k = 0; k < frame.imageWidth; k++){
            
            var rayDirection = currentFramePixel.clone().sub(cameraLocation).normalize();
            
            rays[i].push(new Ray(currentFramePixel.clone(), rayDirection))
            
            currentFramePixel.x += frame.pixelWidth;
            
        }
        
        currentFramePixel = originalFramePixel.clone();
        currentFramePixel.y -= (frame.pixelHeight * (i+1));
    }
    
    return rays;
    
}


function findIntersection(ray, geometries){
    
    var intersections = [];
    
    for(var i = 0; i < geometries.length; i++){
       if(geometries[i].geometry instanceof THREE.SphereGeometry){
           var intersectPoint = getSphereIntersection(ray, geometries[i]);
           if(intersectPoint != null){
               intersections.push({
                   geometry: geometries[i],
                   intersection: intersectPoint
               });
           }
       }
       if(geometries[i].geometry instanceof THREE.PlaneGeometry){
           var intersectPoint = getPlaneIntersection(ray, geometries[i]);
           if(intersectPoint != null){
               intersections.push({
                   geometry: geometries[i],
                   intersection: intersectPoint
               });
           }
       }
    }
    
    var closestIntersectionDistance = 99999;
    var closestGeometry = null;
    for(var i = 0; i < intersections.length; i++){
        
        var distance = intersections[i].intersection.clone().sub(ray.origin).length();
        if(distance < closestIntersectionDistance){
            closestIntersectionDistance = distance;
            closestGeometry = intersections[i].geometry;
        }
        
    }
    
    return closestGeometry;
    
}

function getSphereIntersection(ray, sphere){
    
    var dx = ray.direction.x;
    var dy = ray.direction.y;
    var dz = ray.direction.z;
    
    var xc = sphere.position.x;
    var yc = sphere.position.y;
    var zc = sphere.position.z;
    
    var xo = ray.origin.x;
    var yo = ray.origin.y;
    var zo = ray.origin.z;
    
    var r = sphere.geometry.parameters.radius;
    
    var A = Math.pow(dx, 2) + Math.pow(dy, 2) + Math.pow(dz, 2);
    var B = 2 * ( (dx * (xo - xc)) + (dy * (yo - yc)) + (dz * (zo - zc)));
    var C = Math.pow(xo-xc, 2) + Math.pow(yo-yc, 2) + Math.pow(zo-zc, 2) - Math.pow(r, 2);
    
    var W1 = (-B - Math.sqrt( Math.pow(B, 2) - (4 * A * C)))/ (2 * A);
    var W2 = (-B + Math.sqrt( Math.pow(B, 2) - (4 * A * C)))/ (2 * A);
    
    var W = W1 > 0 ? W1 : W2;
    
    if(W >= 0){
        return ray.origin.clone().add(ray.direction.clone().multiplyScalar(W));
    }
    
    return null;
    
}

function getPlaneIntersection(ray, plane){
    
    var pointOnPlace = plane.geometry.vertices[0];
    var normal = plane.geometry.faces[0].normal;
    var f = normal.dot(pointOnPlace.clone().sub(new THREE.Vector3(0,0,0))) / normal.length();
    var denom = normal.dot(ray.direction)
    var W = ((normal.dot(pointOnPlace)+f) / denom)
    
    if(W > 0){
        return ray.origin.clone().add(ray.direction.clone().multiplyScalar(W));
    }
    
    return null;
    
}

function getPolygonIntersection(ray, plane){
    
    
}


function getNormalOfPlane(plane){
    var v1 = plane.geometry.vertices[plane.geometry.faces[0].a].applyMatrix4(new THREE.Matrix4().makeRotationX(plane.rotation.x));
    var v2 = plane.geometry.vertices[plane.geometry.faces[0].b].applyMatrix4(new THREE.Matrix4().makeRotationX(plane.rotation.x));
    var v3 = plane.geometry.vertices[plane.geometry.faces[0].c].applyMatrix4(new THREE.Matrix4().makeRotationX(plane.rotation.x));
    
    var normal = (v2.clone().sub(v1).cross(v3.clone().sub(v1))).normalize();
    return normal;
    
}

/* Exists already in THREE.js for Vector3 */
function getAngleBetweenVectors(vect1, vect2){
    
    // ||a . b|| = ||a||*||b||*cos(theta)
    // solve for theta and you get the angle between the vectors in radians
    return Math.acos(vect1.dot(vect2)/(vect1.length() * vect2.length()));
}

function radianToDegree(rad){
    return rad * 57.2958; //1 rad = ~57.2958 degrees
}


var scene = new THREE.Scene();
var renderer = new THREE.WebGLRenderer({ antialias: true });
//var img = document.getElementById("render_image");
renderer.setSize( 400, 400 );
document.getElementById("WebGLCanvas").appendChild(renderer.domElement);

var material1 = new THREE.MeshBasicMaterial({color:0x2194ce});
var material2 = new THREE.MeshBasicMaterial({color:0x4FF5ff});
var material3 = new THREE.MeshBasicMaterial({color:0xff0055});

//camera
var camera = new THREE.PerspectiveCamera( 75, 1, 0.1, 1000);

//light
var light = new THREE.PointLight(0xffffff, 1, 100);
light.position.set(1,1,4);

//plane
var planeGeometry = new THREE.PlaneGeometry(4,8);
var plane = new THREE.Mesh(planeGeometry, material1);
plane.geometry.applyMatrix(new THREE.Matrix4().makeRotationX(4.8));

//sphere 1
var sphere1Geometry = new THREE.SphereGeometry(0.5,50,50);
var sphere1 = new THREE.Mesh(sphere1Geometry, material2);

//sphere 2
var sphere2Geometry = new THREE.SphereGeometry(0.45,50,50);
var sphere2 = new THREE.Mesh(sphere2Geometry, material3);


var iframe = new ImageFrame(400,400,new THREE.Vector3(0,0,4));

scene.add(plane);
scene.add(sphere1);
scene.add(sphere2);
scene.add(light);

sphere1.position.set(0,0.1,2);
sphere2.position.set(-0.75,-0.3,1);
plane.geometry.applyMatrix(new THREE.Matrix4().makeTranslation(-0.5,-1,0));
camera.position.x = 0;
camera.position.z = 6;
camera.position.y = 0;


var test = new RayTracer(camera, iframe, scene);

var render = function(){
    
    
    test.render();
    camera.position.z = 4
    renderer.render(scene, camera);
}

render();

window.addEventListener('resize', function(){

}, false)