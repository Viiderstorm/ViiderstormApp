// Ray Tracer
// --Ray
// --Sphere
// --Plane
// --ImagePlane - always a 2x2 grid in world coords, size of pixels per square can change, numbers of ray samples per image pixel can change 
// --Camera
// --View (ImagePlane and Camera)
// Light
// --Scene

function RayTracer(view, scene, supersample){
    var self = this;
    self.view        = view;
    self.scene       = scene;
    self.rays        = null;
    
    self.renderRandom = function(canvasId, supersample){
        self.rays = self.view.getRaysRandom(supersample);
        self.render(canvasId);
    }
    
    self.renderUniform = function(canvasId){
        self.rays = self.view.getRays();
        self.render(canvasId);
    }
    
    self.render      = function(canvasId){
        for(var i = 0; i < self.rays.length; i++){
            for(var k = 0; k < self.rays[i].length; k++){
                
                var colorSamples = []
                
                //supersampling ray loop
                for(var z = 0; z < self.rays[i][k].length; z++){
                    var intersection = self.scene.getIntersection(self.rays[i][k][z]);
                    if(intersection != null){
                        
                        var ambientLight =  self.scene.getAmbientLight(intersection, 0.333);
                        
                        var shadowIntersections = self.scene.getShadowIntersectionsFromPoint(intersection.point);
                        if(!arrayIsNull(shadowIntersections)){
                            
                            //ray intersection is in atleast 1 shadow, maybe more.
                            colorSamples.push(ambientLight);
                            
                        } else {
                            
                            //ray intersection not in any shadow
                            var diffuseLight  = self.scene.getDiffuseLight(intersection, 0.333);
                            var specularLight = self.scene.getSpecularLight(intersection, self.rays[i][k][z], 0.333); 
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
    
    self.averageColor = function(colorSamples){
        var colors = new THREE.Color(0x000000);
        for(var i = 0; i < colorSamples.length; i++){
            colors.add(colorSamples[i].clone().multiplyScalar(1/colorSamples.length));
        }
        return colors;
    }
    
    self.displayImage = function(canvasId){
        var renderImgData = self.view.imagePlane.generateImage();
        var c = document.getElementById(canvasId);
        var ctx = c.getContext("2d");
        var imgData = ctx.createImageData(self.view.imagePlane.displayImageHeight, self.view.imagePlane.displayImageWidth);
        imgData.data.set(renderImgData);
        ctx.putImageData(imgData, 0, 0);
    }
}

function Ray(origin, direction){
    var self = this;
    self.origin    = origin.clone();
    self.direction = direction.clone();
}

function Camera(origin, lookAt, up){
    var self = this;
    self.origin = origin.clone();
    self.lookAt = lookAt.clone();
    self.up     = up.clone();
}

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

function View(origin, lookAt, up, cameraDistance, displayImageHeight, displayImageWidth){
    var self = this;
    self.imagePlane    = new ImagePlane(origin, lookAt, up, displayImageHeight, displayImageWidth);
    self.camera        = new Camera(self.imagePlane.origin.clone().add(self.imagePlane.lookAt.clone().negate().multiplyScalar(cameraDistance)),self.imagePlane.origin, up);
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

function Sphere(origin, radius, widthSegments, heightSegments, material){
    var self = this;
    self.origin    = origin;
    self.radius    = radius;
    self.mesh      = (function(){
        var mesh = new THREE.Mesh(new THREE.SphereGeometry(radius, widthSegments, heightSegments), material)
        mesh.geometry.applyMatrix(new THREE.Matrix4().makeTranslation(origin.x, origin.y, origin.z));
        return mesh;
    })();
    self.vertices  = self.mesh.geometry.vertices;
    self.faces     = self.mesh.geometry.faces;
    self.translate = function(magnitudes){
        self.mesh.geometry.applyMatrix(new THREE.Matrix4().makeTranslation(magnitudes.x, magnitudes.y, magnitudes.z));
        self.origin.add(magnitudes);
        self.vertices = self.mesh.geometry.vertices;
        self.faces = self.mesh.geometry.faces;
    }
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

function Plane(origin, normal, width, height, material){
    var self = this;
    self.origin       = origin;
    self.normal       = normal;
    self.width        = width;
    self.height       = height;
    self.mesh         = (function(){
        var mesh = new THREE.Mesh(new THREE.PlaneGeometry(width, height), material)
        var curNorm = mesh.geometry.faces[0].normal;
        var angle = curNorm.angleTo(self.normal);
        mesh.geometry.applyMatrix(new THREE.Matrix4().makeRotationAxis( curNorm.clone().cross(self.normal), angle));
        mesh.geometry.applyMatrix(new THREE.Matrix4().makeTranslation(origin.x, origin.y, origin.z));
        return mesh;
    })();
    self.vertices     = self.mesh.geometry.vertices;
    self.faces        = self.mesh.geometry.faces;
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

function Light(origin, intensity, color){
    var self = this;
    self.intensity             = intensity;
    self.color                 = color;
    self.origin                = origin;
}



///END SCENE OBJECTS

///BEGIN SCENE

function Scene(){
    var self = this;
    self.geometries      = []
    self.lights          = []
    self.backgroundColor = new THREE.Color(0x000000)
    self.add             = function(object){
        if(object instanceof Light){
            self.lights.push(object);
        }
        if(object instanceof Sphere || object instanceof Plane){
            self.geometries.push(object);
        }
    }
    self.rotate           = function(rMatrix){
        
    }
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
    
    self.getIntersection = function(ray){
        var intersection = self.getIntersections(ray)[0];
        return intersection == undefined ? null : intersection;
    }
    
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
    
    self.getAmbientLight = function(rayIntersection, constant){
        var comp         = 1/self.lights.length;
        var ambientColor = new THREE.Color(0x000000); 
        for(var i = 0 ; i < self.lights.length; ++i){
            ambientColor.add( self.lights[i].color.clone().multiplyScalar(comp)  );   
        }
        
        return ambientColor.multiply(rayIntersection.geometry.mesh.material.color).multiplyScalar(constant)
    }
    self.getDiffuseLight  = function(rayIntersection, constant){
        var diffuseLight = new THREE.Color(0x000000);
        for(var i = 0; i < self.lights.length; ++i){
            var pointToLight = self.lights[i].origin.clone().sub(rayIntersection.point).normalize();
            var scalarVals   = pointToLight.clone().dot(rayIntersection.normal);
            var tmpColor     = self.lights[i].intensity.clone();
            var objectColor  = rayIntersection.geometry.mesh.material.color;
            var finalColor   = tmpColor.multiply(objectColor).multiplyScalar(scalarVals);
            diffuseLight.add(finalColor);
        }
        return diffuseLight.multiplyScalar(constant);
    }
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

function arrayIsNull(arr){
    var count = 0;
    for(var i = 0 ; i < arr.length; i++){
        if(arr[i] == null){
            count += 1;
        }
    }
    return count == arr.length ? true : false;
}

function getReflection(rayIntersection, light){
    var pointToLight = light.origin.clone().sub(rayIntersection.point).normalize();
    var dotProd  = pointToLight.clone().dot(rayIntersection.normal.normalize()) * 2;
    var scaleNormal  = rayIntersection.normal.clone().multiplyScalar(dotProd);
    return scaleNormal.sub(pointToLight); 
}

/// END UTIL


/// INIT

var sc;

$(document).ready(function(){
    
    
    $("#renderButton").click(function(){
        
        var material1 = new THREE.MeshBasicMaterial({ color: 0x2194ce });
        var material2 = new THREE.MeshBasicMaterial({ color: 0x4FF5ff });
        var material3 = new THREE.MeshBasicMaterial({ color: 0xff0055 });


        var sphere1 = new Sphere(new THREE.Vector3(0,0.1,2), 0.5, 50,50, material2);
        var sphere2 = new Sphere(new THREE.Vector3(-0.75,-0.3,1), 0.45, 50,50, material3);
        var plane1  = new Plane(new THREE.Vector3(-0.5,-1,0), new THREE.Vector3(0,1,0), 4, 20, material1);
        var light1  = new Light(new THREE.Vector3(0.5,4,10), new THREE.Color(0xffffff), new THREE.Color(0xffffff));
        //var light2  = new Light(new THREE.Vector3(-0.5,4,10), new THREE.Color(0xffffff), new THREE.Color(0xffffff));

        sc = new Scene();
        sc.add(sphere1);
        sc.add(sphere2);
        sc.add(plane1);
        sc.add(light1);
        //sc.add(light2);

        var view    = new View(new THREE.Vector3(0,0,3), new THREE.Vector3(0,0,-3).normalize(), new THREE.Vector3(0,1,0), 2, 400, 400);
        // var view2    = new View(new THREE.Vector3(0,0,-4), new THREE.Vector3(0,0,4).normalize(), new THREE.Vector3(0,1,0), 2, 400, 400);
        // var view3    = new View(new THREE.Vector3(-2,0,0), new THREE.Vector3(2,0,0).normalize(), new THREE.Vector3(0,1,0), 2, 400, 400);
        // var view4    = new View(new THREE.Vector3(2,0,0), new THREE.Vector3(-2,0,0).normalize(), new THREE.Vector3(0,1,0), 2, 400, 400);

        var rt      = new RayTracer(view, sc);
        // var rt2      = new RayTracer(view2, sc);
        // var rt3      = new RayTracer(view3, sc);
        // var rt4      = new RayTracer(view4, sc);
        //rt.renderUniform("renderCanvas")
        rt.renderRandom("renderCanvas", 30)
        // rt2.renderUniform("renderCanvas2")
        // rt3.renderUniform("renderCanvas3")
        // rt4.renderUniform("renderCanvas4")
        
    })
    
})