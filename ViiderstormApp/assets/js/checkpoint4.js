/* Class RayTracer
* Description: Object that is used for managing and tracing rays as well as converting data to images
* View view - The View used to view the scene with
* Scene scene - The scene to view
* Integer supersample - The number of samples per pixel
*/
function RayTracer(view, scene, supersample){
    var self = this;
    self.view        = view;
    self.scene       = scene;
    self.rays        = null;
    
    /*
    * Initializes the Ray Tracer with a set number of Random samples then renders the scene
    * String canvasId - The DOMElement ID for the rendering canvas
    * Integer supersample - The number of samples per pixel
    */
    self.renderRandom = function(canvasId, supersample){
        self.rays = self.view.getRaysRandom(supersample);
        self.render(canvasId);
    }
    
    /*
    *Initializes the Ray Tracer with a uniform sampling (1 SPP) and then renders the scene
    * String canvasId - The DOMElement ID for the rendering canvas
    */
    self.renderUniform = function(canvasId){
        self.rays = self.view.getRays();
        self.render(canvasId);
    }
    
    /*
    *Renders the scene by tracing each pixel's ray
    * String canvasID - The DOMElement ID for the rendering canvas
    */ 
    self.render      = function(canvasId){
        for(var i = 0; i < self.rays.length; i++){
            for(var k = 0; k < self.rays[i].length; k++){
                
                var colorSamples = []
                
                //supersampling ray loop
                for(var z = 0; z < self.rays[i][k].length; z++){
                    var intersection = self.scene.getIntersection(self.rays[i][k][z]);
                    if(intersection != null){
                        
                        var ambientLight =  self.scene.getAmbientLight(intersection, 0.5);
                        
                        var shadowIntersections = self.scene.getShadowIntersectionsFromPoint(intersection.point);
                        if(!arrayIsNull(shadowIntersections)){
                            
                            //ray intersection is in atleast 1 shadow, maybe more.
                            colorSamples.push(ambientLight);
                            
                        } else {
                            
                            //ray intersection not in any shadow
                            var diffuseLight  = self.scene.getDiffuseLight(intersection, 0.4);
                            var specularLight = self.scene.getSpecularLight(intersection, self.rays[i][k][z], 0.2); 
                            var tmpColor      = new THREE.Color(0x000000);
                            tmpColor.add(ambientLight);
                            tmpColor.add(diffuseLight);
                            tmpColor.add(specularLight);
                            colorSamples.push(tmpColor);
                            
                        }
                    } else {
                        
                        //ray hits background
                        colorSamples.push(self.scene.backgroundColor);
                    }
                }
                
                var finalColor = self.averageColor(colorSamples)
                if(finalColor != null){
                    self.view.imagePlane.pixels[i][k].color = finalColor;
                }
                
            }
        }
        self.displayImage(canvasId)
    } 
    
    /*
    * Averages a list of color samples
    * Color[] colorSamples - A List of Color Samples
    */
    self.averageColor = function(colorSamples){
        var colors = new THREE.Color(0x000000);
        for(var i = 0; i < colorSamples.length; i++){
            colors.add(colorSamples[i].clone().multiplyScalar(1/colorSamples.length));
        }
        return colors;
    }
    
    /*
    * Puts the render data into the canvas
    * String canvasId - The DOMElement ID for the rendering canvas
    */
    self.displayImage = function(canvasId){
        var renderImgData = self.view.imagePlane.generateImage();
        var c = document.getElementById(canvasId);
        var ctx = c.getContext("2d");
        var imgData = ctx.createImageData(self.view.imagePlane.displayImageHeight, self.view.imagePlane.displayImageWidth);
        imgData.data.set(renderImgData);
        ctx.putImageData(imgData, 0, 0);
    }
}

/* Class Ray
* THREE.Vector3 origin - origin of the ray
* THREE.Vector3 direction - direction of the ray
*/
function Ray(origin, direction){
    var self = this;
    self.origin    = origin.clone();
    self.direction = direction.clone();
}

/* Class Camera
* THREE.Vector3 origin - origin of the camera
* THREE.Vector3 lookAt - point to look at
* THREE.Vector3 up     - Up direction in relation to the scene
*/
function Camera(origin, lookAt, up){
    var self = this;
    self.origin = origin.clone();
    self.lookAt = lookAt.clone();
    self.up     = up.clone();
}

/* Class ImagePlane
* THREE.Vector3 origin - origin of the Image Plane
* THREE.Vector3 lookAt - point to look at
* THREE.Vector3 up     - Up direction in relation to the scene
* Integer displayImageHeight - height of the image plane
* Integer displayImageWidth  - width of the image plane
*/ 
function ImagePlane(origin, lookAt, up, displayImageHeight, displayImageWidth){
    var self = this;
    self.origin             = origin.clone();
    self.lookAt             = lookAt.clone().normalize();
    self.up                 = Math.abs(lookAt.clone().dot(up)) <= 0.01 ? up.clone().normalize() : (function(){throw new Error('Up direction not within Image Plane');})();
    self.down               = self.up.clone().negate();
    self.right              = lookAt.clone().cross(self.up).normalize();
    self.left               = self.right.clone().negate();
    self.displayImageHeight = displayImageHeight;
    self.displayImageWidth  = displayImageWidth;
    self.worldPixelHeight   = 2 / displayImageHeight;
    self.worldPixelWidth    = 2 / displayImageWidth;
    self.topLeftPixel       = self.origin.clone().add(self.up.clone().multiplyScalar(self.worldPixelHeight * (displayImageHeight/2))
                              .add(self.left.clone().multiplyScalar(self.worldPixelWidth * (displayImageWidth/2))));
    self.pixels             = (function(x,y,z){
        var pix             = []  
        var leftMostPoint   = new THREE.Vector3(x,y,z);
        for(var col = 0; col < self.displayImageHeight; col++){
            pix.push([]);
            for(var row = 0; row < self.displayImageWidth; row++){
                pix[col].push({
                    origin: leftMostPoint.clone().add(self.right.clone().multiplyScalar(self.worldPixelWidth * row)),
                    color: new THREE.Color(0x000000)
                })
            }
            leftMostPoint.add(self.down.clone().multiplyScalar(self.worldPixelHeight));
        }
        return pix;
    })( self.topLeftPixel.x, self.topLeftPixel.y, self.topLeftPixel.z );
    
    /*
    * Take the rendering data in the Image Plane and turn it into a format that the canvas can read
    */
    self.generateImage      = function(){
        var tmp             = [];
        for (var i = 0; i < self.displayImageHeight; i++) {
            for (var k = 0; k < self.displayImageWidth; k++) {

                var style = self.pixels[i][k].color.getStyle().replace("rgb(", "").replace(")", "")
                var values = style.split(",");

                tmp.push(parseInt(values[0]));
                tmp.push(parseInt(values[1]));
                tmp.push(parseInt(values[2]));
                tmp.push(255);
            }
        }
        return new Uint8ClampedArray(tmp);
    }
    
}

/* Class View
* A combination of the camera and image plane that makes it's manueverability much easier.
* THREE.Vector3 origin - origin of the image plane
* THREE.Vector3 lookAt - the point to look at
* THREE.Vector3 up     - The up direction in relation to the camera's view
* Integer displayImageHeight - height of the image plane
* Integer displayImageWidth  - width of the image plane
*/
function View(origin, lookAt, up, cameraDistance, displayImageHeight, displayImageWidth){
    var self = this;
    self.imagePlane    = new ImagePlane(origin, lookAt, up, displayImageHeight, displayImageWidth);
    self.camera        = new Camera(self.imagePlane.origin.clone().add(self.imagePlane.lookAt.clone().negate().multiplyScalar(cameraDistance)),self.imagePlane.origin, up);
    
    /*
    * Creates the set of ray for the set of pixels
    * Function sampleFun - function to modify the generation of the rays. If undefined, it defaults to uniform sampling.
    */
    self.getRays       = function(sampleFunc){
        var rays       = [];
        for(var i = 0; i < self.imagePlane.pixels.length; i++){
            rays.push([])
            for(var k = 0; k < self.imagePlane.pixels[i].length; k++){
                if(sampleFunc != undefined){
                    rays[i].push(sampleFunc(self.imagePlane.pixels[i][k].origin,
                                            self.imagePlane.right.clone().multiplyScalar(self.imagePlane.worldPixelWidth),
                                            self.imagePlane.down.clone().multiplyScalar(self.imagePlane.worldPixelHeight),
                                            self.camera.origin.clone()))
                } else {
                    var rayOrigin = self.imagePlane.pixels[i][k].origin.clone()
                                    .add(self.imagePlane.right.clone().multiplyScalar(self.imagePlane.worldPixelWidth/2))
                                    .add(self.imagePlane.down.clone().multiplyScalar(self.imagePlane.worldPixelHeight/2))
                    rays[i].push([new Ray(rayOrigin, rayOrigin.clone().sub(self.camera.origin).normalize())]);
                }
            }
        }
        return rays;
    }
    
    /*
    * Modifier function for generating random ray direction given the View information. Picks a random direction in the image plane for each pixel, potentially multiple times depending on the sampleNumber.
    * Integer sampleNumber - number of samples to take per pixel
    */
    self.getRaysRandom = function(sampleNumber){
        return self.getRays(function(pixelOrigin, rightVector, downVector, cameraOrigin){
            var rays          = []
            var distanceRight = rightVector.length();
            var distanceDown  = downVector.length();
            var normRight     = rightVector.clone().normalize();
            var normDown      = downVector.clone().normalize();
            while(rays.length < sampleNumber){
                var rayOrigin = pixelOrigin.clone()
                                .add(normRight.clone().multiplyScalar(distanceRight * Math.random()))
                                .add(normDown.clone().multiplyScalar(distanceDown * Math.random()))
                rays.push(new Ray(rayOrigin, rayOrigin.clone().sub(cameraOrigin).normalize()))
            }
            
            return rays;
        });
    }
}


///BEGIN SCENE OBJECTS

/* Class Sphere
* Description: Sphere object data structure. Built on top of the existing THREE.js gemoetry objects.
* THREE.Vector3 origin - origin of the sphere
* Integer radius - radius of the sphere
* Integer widthSegments - number of segments going across (left-right) when creating the vertices
* Integer heightSegments - number of segments going across (top-down) when creating the vertices
* THREE.Material material - THREE.js material object that defines the properties of the sphere
*/
function Sphere(origin, radius, widthSegments, heightSegments, material){
    var self = this;
    self.origin    = origin;
    self.radius    = radius;
    self.mesh      = (function(){
        var mesh = new THREE.Mesh(new THREE.SphereGeometry(radius, widthSegments, heightSegments), new THREE.MeshNormalMaterial())
        mesh.geometry.applyMatrix(new THREE.Matrix4().makeTranslation(origin.x, origin.y, origin.z));
        return mesh;
    })();
    self.vertices  = self.mesh.geometry.vertices;
    self.faces     = self.mesh.geometry.faces;
    self.material  = material;
    self.translate = function(magnitudes){
        self.mesh.geometry.applyMatrix(new THREE.Matrix4().makeTranslation(magnitudes.x, magnitudes.y, magnitudes.z));
        self.origin.add(magnitudes);
        self.vertices = self.mesh.geometry.vertices;
        self.faces = self.mesh.geometry.faces;
    }
    
    /*
    * Calculates intersection with the sphere given a ray. Returns null if no intersection is found.
    * Ray ray - ray to check intersection for
    */ 
    self.getIntersect   = function(ray){
        var dx = ray.direction.x;
        var dy = ray.direction.y;
        var dz = ray.direction.z;

        var xc = self.origin.x;
        var yc = self.origin.y;
        var zc = self.origin.z;

        var xo = ray.origin.x;
        var yo = ray.origin.y;
        var zo = ray.origin.z;

        var r = self.radius;

        var A = Math.pow(dx, 2) + Math.pow(dy, 2) + Math.pow(dz, 2);
        var B = 2 * ((dx * (xo - xc)) + (dy * (yo - yc)) + (dz * (zo - zc)));
        var C = Math.pow(xo - xc, 2) + Math.pow(yo - yc, 2) + Math.pow(zo - zc, 2) - Math.pow(r, 2);

        var W1 = (-B - Math.sqrt(Math.pow(B, 2) - (4 * A * C))) / (2 * A);
        var W2 = (-B + Math.sqrt(Math.pow(B, 2) - (4 * A * C))) / (2 * A);

        var W = W1 > 0 ? W1 : W2;

        if (W > (0)) { //Javascript math is bad, I gotta take into account error.
            
            var point = ray.origin.clone().add(ray.direction.clone().multiplyScalar(W)); 
            var normal = point.clone().sub(self.origin).normalize();
            return {
                point: point,
                normal: normal 
            }
                
        }

        return null;
    }
    
}


/* Class Plane
* Description: Data strucutre for the Plane
* THREE.Vector3 origin - origin of the plane
* THREE.Vector3 normal - normal of the plane
* Integer width - width of the plane
* Integer height - height of the plane
* THREE.Material material - THREE.js material object that defines the properties of the plane
*/
function Plane(origin, normal, width, height, material){
    var self = this;
    self.origin       = origin;
    self.normal       = normal;
    self.width        = width;
    self.height       = height;
    self.mesh         = (function(){
        var mesh = new THREE.Mesh(new THREE.PlaneGeometry(width, height), new THREE.MeshNormalMaterial())
        var curNorm = mesh.geometry.faces[0].normal;
        var angle = curNorm.angleTo(self.normal);
        mesh.geometry.applyMatrix(new THREE.Matrix4().makeRotationAxis( curNorm.clone().cross(self.normal), angle));
        mesh.geometry.applyMatrix(new THREE.Matrix4().makeTranslation(origin.x, origin.y, origin.z));
        return mesh;
    })();
    self.vertices     = self.mesh.geometry.vertices;
    self.faces        = self.mesh.geometry.faces;
    self.material     = material;
    self.translate    = function(magnitudes){
        self.mesh.geometry.applyMatrix(new THREE.Matrix4().makeTranslation(magnitudes.x, magnitudes.y, magnitudes.z));
        self.origin.add(magnitudes);
        self.vertices = self.mesh.geometry.vertices;
        self.faces    = self.mesh.geometry.faces;
    }
    self.rotate       = function(rMatrix){
        self.mesh.geometry.applyMatrix(rMatrix);
        self.vertices = self.mesh.geometry.vertices;
        self.faces    = self.mesh.geometry.faces;
    }
    
    /*
    * Intersection calculation for plane
    * Ray ray - ray to check intersection with
    */
    self.getIntersect = function(ray){
        var intersect = self.getIntersectPlane(ray);
        if(intersect != null && self.isPointInPlaneGeometry(intersect)){
            return {
                point: intersect,
                normal: self.normal.clone()
            }
                
        }
        return null;
    }
    
    /*
    * Helper function for ray intersection
    * Ray ray - ray to check intersection with
    */
    self.getIntersectPlane = function(ray){
        var point = self.vertices[0]
        var rayToPoint  = ray.origin.clone().sub(point);
        var denom = self.normal.clone().dot(ray.direction);
        var num = self.normal.clone().dot(rayToPoint);
        var t = -(num/denom);
        
        if(t >= 0.00000001){ //Javascript math is bad, I gotta take into account error.
            return ray.origin.clone().add(ray.direction.clone().multiplyScalar(t));
        }
        return null;
    }
    
    /*
    * Checks if the intersection point is within the plane geometry
    * THREE.Vector3 point - intersection point
    */ 
    self.isPointInPlaneGeometry = function(point){
        for(var i = 0; i < self.faces.length; ++i){
            var verts = [ self.vertices[self.faces[i].a],
                          self.vertices[self.faces[i].b],
                          self.vertices[self.faces[i].c] ];
            if(checkPointIsInTriangle(point, verts)){
                return true;
            }
        }
        return false;
    }
}

/* Class Light
* Description: Simple light data structure
* THREE.Vector3 origin - origin of the light
* THREE.Color intensity - intensity of the light
* THREE.Color color - color of the light
*/
function Light(origin, intensity, color){
    var self = this;
    self.intensity             = intensity;
    self.color                 = color;
    self.origin                = origin;
}

///END SCENE OBJECTS

///BEGIN MATERIAL OBJECTS

/* Class BasicMaterial
* Description: Simple material data structure that replaces the one that THREE.js has
* THREE.Color color - color of the material
*/
function BasicMaterial(color){
    var self = this;
    self.color    = new THREE.Color(color);
    
    /* Accessor method
    * Returns the color of the material
    * var intersection - The intersection object
    */
    self.getColor = function(intersection){
        return self.color;
    }
}

/* Class CheckerMaterial
* Description: Material data structure that creates a checkered pattern on the surface of an object
* THREE.Color color1 - first color
* THREE.Color color2 - second color
* Integer widthInChecks - width of the checkerboard in number of squares
* Integer heightInChecks - height of the checkboard in number of squares
*/
function CheckerMaterial(color1, color2, widthInChecks, heightInChecks){
    var self = this;
    self.color1            = new THREE.Color(color1);
    self.color2            = new THREE.Color(color2);
    self.wChecks           = widthInChecks;
    self.hChecks           = heightInChecks;
    
    /*
    * Accessor method
    * Returns the color of the material given the intersection object
    * var intersection - The intersection object
    */
    self.getColor          = function (intersection){
        switch(intersection.geometry.constructor.name){
            case "Sphere":
                return self.getSphereColor(intersection);
            case "Plane":
                return self.getPlaneColor(intersection);
            default:
                throw new Error("Material not defined for Scene Object: " + intersection.geometry.constructor.name); 
        }
    }
    
    /*
    * Calculates what color should be at a given intersection on a plane.
    * var intersection = The intersection object
    */
    self.getPlaneColor = function(intersection){
        var backCorner      = intersection.geometry.vertices[0].clone();
        var widthVector     = intersection.geometry.vertices[1].clone().sub(backCorner).normalize();
        var intersectVector = intersection.point.clone().sub(backCorner);
        var h               = intersectVector.length();
        var theta           = intersectVector.angleTo(widthVector);
        var width           = h * Math.cos(theta);
        var height          = h * Math.sin(theta);
        var x               = width / intersection.geometry.width;
        var y               = height / intersection.geometry.height;
        var isEvenCol       = ((Math.round(x / (1/self.wChecks))) % self.wChecks) % 2 == 0
        var isEvenRow       = ((Math.round(y / (1/self.hChecks))) % self.hChecks) % 2 == 0
        if(isEvenCol == isEvenRow){
            return self.color1;
        }
        return self.color2;
    }
    
    /*
    * Calculates what color should be at the given intersection on a sphere
    * var intersection = The intersection object
    */
    self.getSphereColor = function(intersection){
        return self.color1; 
    }  
}

//END MATERIAL OBJECTS

///BEGIN SCENE
/* Class Scene
* Description: Object that maintains the scene and tracing of the rays
*/
function Scene(){
    var self = this;
    self.geometries      = []
    self.lights          = []
    self.backgroundColor = new THREE.Color(0x000000)
    
     /* 
    * Adds an object to the scene
    * var object - object to add to the scene. Must be of type Light, Sphere, or Plane
    */
    self.add             = function(object){
        if(object instanceof Light){
            self.lights.push(object);
        }
        if(object instanceof Sphere || object instanceof Plane){
            self.geometries.push(object);
        }
    }
    
    /*
    * Finds all intersections for a ray in the scene then picks the closest one
    * Ray ray - ray to find intersection points for
    */
    self.getIntersections = function(ray){
        var intersections = [];
        for(var i = 0 ; i < self.geometries.length; ++i){
            var intersection = self.geometries[i].getIntersect(ray) 
            if(intersection != null){
                intersections.push({
                    point: intersection.point,
                    normal: intersection.normal,
                    geometry: self.geometries[i]
                });
            }
        }
        intersections.sort(function(a, b){
            var distanceA = a.point.clone().sub(ray.origin).length();
            var distanceB = b.point.clone().sub(ray.origin).length();
            return distanceA < distanceB ? -1 : (distanceA > distanceB ? 1: 0); 
        })
        return intersections;
    }
    
    /*
    * Finds a single intersection in the scene for a given ray
    * Ray ray - ray to find intersection point for
    */
    self.getIntersection = function(ray){
        var intersection = self.getIntersections(ray)[0];
        return intersection == undefined ? null : intersection;
    }
    
    /*
    * Finds the shadow ray intersection points for a shadow ray originating from a point
    * THREE.Vector3 point - intersection point or shadow ray origin point
    */
    self.getShadowIntersectionsFromPoint = function(point){
        var intersections = [];
        for(var i = 0; i < self.lights.length; i++){
            var direction = self.lights[i].origin.clone().sub(point).normalize();
            
            //Move point a tiny bit closer to the light source to deal with extremely miniscule Floating point errors
            //in the intersection calculation.
            //This keeps errors from occuring during the shadow intersection calculations.
            var scale     = point.clone().add(direction.clone().multiplyScalar(0.0000000000001));
            var sRay = new Ray(scale, self.lights[i].origin.clone().sub(point).normalize());
            intersections.push(self.getIntersection(sRay));
        }
        return intersections;
    }
    
    /*
    * Calculates the ambient light term of an intersection point within the scene
    * THREE.Vector3 rayIntersection - intersection object
    * Float constant - Ambient term constant
    */
    self.getAmbientLight = function(rayIntersection, constant){
        var comp         = 1/self.lights.length;
        var ambientColor = new THREE.Color(0x000000); 
        for(var i = 0 ; i < self.lights.length; ++i){
            ambientColor.add( self.lights[i].color.clone().multiplyScalar(comp)  );   
        }
        
        return ambientColor.multiply(rayIntersection.geometry.material.getColor(rayIntersection)).multiplyScalar(constant)
    }
    
    /*
    * Calculates the diffuse light term of an intersection point within the scene
    * THREE.Vector3 rayIntersection - intersection object
    * Float constant - diffuse term constant
    */
    self.getDiffuseLight  = function(rayIntersection, constant){
        var diffuseLight = new THREE.Color(0x000000);
        for(var i = 0; i < self.lights.length; ++i){
            var pointToLight = self.lights[i].origin.clone().sub(rayIntersection.point).normalize();
            var scalarVals   = pointToLight.clone().dot(rayIntersection.normal);
            var tmpColor     = self.lights[i].intensity.clone();
            var objectColor  = rayIntersection.geometry.material.getColor(rayIntersection);
            var finalColor   = tmpColor.multiply(objectColor).multiplyScalar(scalarVals);
            diffuseLight.add(finalColor);
        }
        return diffuseLight.multiplyScalar(constant);
    }
    
    /*
    * Calculates the specular light term of an intersection point within the scene
    * THREE.Vector3 rayIntersection - intersection object
    * Float constant - specular term constant
    */
    self.getSpecularLight = function(rayIntersection, ray, constant){
        var specularLight = new THREE.Color(0x000000);
        for(var i = 0; i < self.lights.length; i++){
            var reflection    = getReflection(rayIntersection, self.lights[i]).normalize();
            var viewDirection = ray.direction.clone().negate();
            var scalarVals    = Math.pow(reflection.dot(viewDirection), 30);
            var tmpColor      = self.lights[i].intensity.clone();
            var finalColor    = tmpColor.multiplyScalar(scalarVals);
            specularLight.add(finalColor);
        }
        return specularLight.multiplyScalar(constant);
    }
}

///END SCENE









/// BEGIN UTIL

/*
* Description: Checks that the intersection point lies within a face
* THREE.Vector3 point - Intersection point
* THREE.Vector3[] vertices- vertices of the plane's faces
*/
function checkPointIsInTriangle(point, vertices){
    var linestoVertices = []
    for(var i = 0; i < vertices.length; i++){
        linestoVertices.push(vertices[i].clone().sub(point));
    }
    
    var a,b,c;
    a = linestoVertices[0].angleTo(linestoVertices[1]) * (180/Math.PI);
    b = linestoVertices[1].angleTo(linestoVertices[2]) * (180/Math.PI);
    c = linestoVertices[2].angleTo(linestoVertices[0]) * (180/Math.PI);
    
    var sum = a + b + c;
    var epsilon = 0.1;
    if(Math.abs(Math.round(sum) - 360) < epsilon){
        return true;
    }
    return false;
}

/*
* Description: Checks whether an array is full of null values
*/
function arrayIsNull(arr){
    var count = 0;
    for(var i = 0 ; i < arr.length; i++){
        if(arr[i] == null){
            count += 1;
        }
    }
    return count == arr.length ? true : false;
}

/*
* Description: Calculates the relfection direction of a point and direction
* THREE.Vector3 rayIntersection - intersection object
* THREE.Vector3 light - light object
*/
function getReflection(rayIntersection, light){
    var pointToLight = light.origin.clone().sub(rayIntersection.point).normalize();
    var dotProd  = pointToLight.clone().dot(rayIntersection.normal.normalize()) * 2;
    var scaleNormal  = rayIntersection.normal.clone().multiplyScalar(dotProd);
    return scaleNormal.sub(pointToLight); 
}

/// END UTIL


/// INIT

var sc;

var checkWidth = 10;
var checkHeight = 25;

var sphereMat1 = new BasicMaterial(0x4FF5ff);
var sphereMat2 = new BasicMaterial(0xff0055);
var planeMat   = new CheckerMaterial(0xff0000, 0xffff00, checkWidth, checkHeight);

var sphere1 = new Sphere(new THREE.Vector3(0,0.1,2), 0.5, 50,50, sphereMat1);
var sphere2 = new Sphere(new THREE.Vector3(-0.75,-0.3,1), 0.45, 50,50, sphereMat2);
var plane1  = new Plane(new THREE.Vector3(-0.5,-1,0), new THREE.Vector3(0,1,0), 4, 20, planeMat);
var light1  = new Light(new THREE.Vector3(0.5,4,10), new THREE.Color(0xffffff), new THREE.Color(0xffffff));

sc = new Scene();
sc.add(sphere1);
sc.add(sphere2);
sc.add(plane1);
sc.add(light1);
$(document).ready(function(){
    
    
    $("#renderButtonUniform").click(function(){
        
        //$("#renderButtonUniform").remove();
        
        checkWidth  = parseInt($("#widthInput").val()) || 10;
        checkHeight = parseInt($("#heightInput").val()) || 25;
                
        planeMat.wChecks = checkWidth
        planeMat.hChecks = checkHeight
        
        var canvas  = $("#renderCanvasUniform")
        var view    = new View(new THREE.Vector3(0,0,3), new THREE.Vector3(0,0,-3).normalize(), new THREE.Vector3(0,1,0), 2, canvas.height(), canvas.width());
        
        var rt      = new RayTracer(view, sc);
        rt.renderUniform("renderCanvasUniform");
        
        
    })
    
})