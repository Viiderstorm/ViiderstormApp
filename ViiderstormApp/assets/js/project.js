// Ray Tracer
// --Ray
// --Sphere
// --Plane
// --ImagePlane - always a 2x2 grid in world coords, size of pixels per square can change, numbers of ray samples per image pixel can change 
// --Camera
// --View (ImagePlane and Camera)
// --Light
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
    self.getShaderCode  = function(objectIdNumber){
        var qualifier = function(s){
            var i = "sphere" + objectIdNumber + "_"; //sphereN_
            return i + s;   
        }
        var origin    =  "uniform vec3 " + qualifier("origin") + ";\n";// uniform vec3 sphereN_origin;
        var radius    =  "uniform float " + qualifier("radius") +";\n";// uniform float sphereN_radius;
        var id        =  "uniform int " + qualifier("id") + ";\n";// uniform float sphereN_id;
        var intersect_uniforms  = {}
        intersect_uniforms[qualifier("origin")] = {type: 'v3', value: self.origin};
        intersect_uniforms[qualifier("radius")] = {type: 'f', value: self.radius};
        intersect_uniforms[qualifier("id")] = {type: 'i', value: objectIdNumber};
        
        var color_uniforms = {}
        color_uniforms[qualifier("id")] = {type: 'f', value: objectIdNumber};
        color_uniforms[qualifier("color")] = {type: 'c', value: self.mesh.material.color};
        
        var func_header = "vec4 " + qualifier("getIntersection") + "(vec3 rayOrigin, vec3 rayDirection){\n";
        var func_contentA = "\tfloat A = pow(rayDirection.x, 2.0) + pow(rayDirection.y, 2.0) + pow(rayDirection.z, 2.0);\n"
        var func_contentB = "\tfloat B = 2.0 * ((rayDirection.x * (rayOrigin.x - " + qualifier("origin") + ".x)) +" +
                                           "(rayDirection.y * (rayOrigin.y - " + qualifier("origin") + ".y)) +" +
                                           "(rayDirection.z * (rayOrigin.z - " + qualifier("origin") + ".z)));\n";
        var func_contentC = "\tfloat C = pow(rayOrigin.x - " + qualifier("origin") + ".x, 2.0) + " +
                                       "pow(rayOrigin.y - " + qualifier("origin") + ".y, 2.0) + " +
                                       "pow(rayOrigin.z - " + qualifier("origin") + ".z, 2.0) - " +
                                       "pow(" + qualifier("radius") + ", 2.0);\n";
        var func_contentW1 = "\tfloat W1 = (-B - sqrt(pow(B, 2.0) - (4.0 * A * C))) / (2.0 * A);\n";
        var func_contentW2 = "\tfloat W2 = (-B + sqrt(pow(B, 2.0) - (4.0 * A * C))) / (2.0 * A);\n";
        var func_contentW  = "\tfloat W  = W1 > 0.0 ? W1 : W2;\n";
        var func_contentif        = "\tif( W > 0.0 ){\n \t\treturn vec4(rayOrigin + (rayDirection * W), "+qualifier("id") +");\n\t} else {\n \t\treturn vec4(0,0,0,-1);\n\t}\n}\n";                               
        
        var boilerCode   = origin + radius + id;
        var functionCode =  func_header + func_contentA + 
                            func_contentB + func_contentC + func_contentW1 + 
                            func_contentW2 + func_contentW + func_contentif;
        
        return { qualifier, intersect_uniforms, boilerCode, functionCode, color_uniforms }
        
        
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
    
    self.getShaderCode = function(objectIdNumber){
        var qualifier = function(s){
            var i = "plane" + objectIdNumber + "_";
            return i + s;   
        }
        
        var faceToVec2Str = function(s){
            var p1 = "vec3(" + self.vertices[s.a].x + "," + self.vertices[s.a].y + "," + self.vertices[s.a].z + ")";
            var p2 =  "vec3(" + self.vertices[s.b].x + "," + self.vertices[s.b].y + "," + self.vertices[s.b].z + ")";
            var p3 =  "vec3(" + self.vertices[s.c].x + "," + self.vertices[s.c].y + "," + self.vertices[s.c].z + ")";
            return [p1,p2,p3].join(',');
        }
        
        var origin    =  "uniform vec3 " + qualifier("origin") + ";\n";// uniform vec3 planeN_origin;
        var normal = "uniform vec3 " + qualifier("normal") + ";\n";//uniform vec3 planeN_normal;
        var pointOnPlane = "uniform vec3 " + qualifier("pointOnPlane") + ";\n";//uniform vec2 planeN_pointOnPlane;
        var id        =  "uniform int " + qualifier("id") + ";\n";// uniform float planeN_id;
        var faces     = "struct " + qualifier("face") + "{\n" +
                        "\t vec3 point1;\n" +
                        "\t vec3 point2;\n" + 
                        "\t vec3 point3;\n" +
                        "}";
        var faceInstances = []                        
        for(var i = 0; i < self.faces.length; i++){
            faceInstances.push( qualifier("face_tri") + i);
        }
        faces += faceInstances.join(",") + ";\n";
        
        var faceDefine = "\t" +qualifier("face") + " " + qualifier("faces") + "[" + self.faces.length + "];\n"; 
        for(var i = 0; i < self.faces.length; i++){
            faceDefine += "\t" + qualifier("face_tri") + i + " = " + qualifier("face") + "(" + faceToVec2Str(self.faces[i]) + ");\n";
            faceDefine += "\t" + qualifier("faces") + "[" + i + "] = " + qualifier("face_tri") + i  + ";\n";
        }
        
        var func1 = "vec4 " + qualifier("getIntersection") + 
                    "(vec3 rayOrigin, vec3 rayDirection," + qualifier("face") + "[" + self.faces.length + "] faces){\n";
        var func1_content1 = "\tvec3 rayToPoint = rayOrigin - " + qualifier("pointOnPlane") + ";\n";
        var func1_content2 = "\tfloat denom = dot(" + qualifier("normal") + ", " + "rayDirection);\n";
        var func1_content3 = "\tfloat num   = dot(rayToPoint, " + qualifier("normal") + ");\n";
        var func1_content4 = "\tfloat t     = -(num/denom);\n";
        var func1_content5 = "\tif(t >= 0.00000001){\n" +
                             "\t\tvec4 intersectionPoint = vec4(rayOrigin + (rayDirection * t)," + objectIdNumber + ");\n" +
                             "\t\tbool triangleCheck1;\n" +
                             "\t\tbool triangleCheck2;\n" +
                             "\t\ttriangleCheck1 = " + qualifier("isPointInPlane") + "(intersectionPoint.xyz, faces[0]);\n" +
                             "\t\ttriangleCheck2 = " + qualifier("isPointInPlane") + "(intersectionPoint.xyz, faces[1]);\n" +
                             "\t\tif(triangleCheck1 || triangleCheck2){\n" +
                             "\t\t\t return intersectionPoint;\n" +
                             "\t\t}\n" +
                             "\t}\n" +
                             "\treturn vec4(0,0,0,-1);\n";
        var func1_end = "}\n"
        
        var func2 = "bool " + qualifier("isPointInPlane") + "(vec3 intersectionPoint," + qualifier("face") + " tri){\n";
        for(var i = 0; i < 3; i++){
            func2 += "\tvec3 linestoVertices" + i + " = tri.point" + (i+1) + " - intersectionPoint;\n"
        }
        for(var i = 0; i < 3; i++){
            func2 += "\tfloat sum" + i + " = degrees(acos( dot(normalize(linestoVertices" + i + "), normalize(linestoVertices" + ((i+1)%3) + "))));\n";
        }
        func2 += "\tfloat totalSum = ";
        var sums = [];
        for(var i = 0; i < 3; i++){
            sums.push("sum" + i);
        }
        func2 += sums.join("+") + ";\n"
        func2 += "\tif(abs(totalSum - 360.0) < 0.1){\n" +
                 "\t\t return true;\n" +
                 "\t} else {\n" +
                 "\t\t return false;\n" +
                 "\t}\n";
        func2_end = "}\n";
                         
        var intersect_uniforms = {};
        intersect_uniforms[qualifier("origin")] = {type: 'v3', value: self.origin};
        intersect_uniforms[qualifier("normal")] = {type: 'v3', value: self.normal};
        intersect_uniforms[qualifier("pointOnPlane")] = {type: 'v3', value: self.vertices[0]};
        intersect_uniforms[qualifier("id")] = {type: 'i', value: objectIdNumber};
        
        var color_uniforms = {};
        color_uniforms[qualifier("id")] = {type: 'f', value: objectIdNumber};
        color_uniforms[qualifier("color")] = {type: 'c', value: self.mesh.material.color};
        
        var boilerCode = origin + normal + pointOnPlane + id + faces;
        var functionCode = func2 + func2_end;
        functionCode += func1 + func1_content1 + func1_content2 + func1_content3 + func1_content4 + func1_content5 + func1_end;
        
        return {qualifier, intersect_uniforms, boilerCode, functionCode, faceDefine, color_uniforms};
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

function flattenArray(arr){
    var concat = function(a,b){
        return a.concat(b);
    }
    var tmpArr = arr.reduce(concat).reduce(concat).reduce(concat);
}

/// END UTIL


function SceneRenderer(id,width, height){
    var self = this;
    self.width = width;
    self.height = height;
    self.scRenderer = new THREE.WebGLRenderer();
    self.scCamera = (function(){
        var i = new THREE.OrthographicCamera( 1 / -2, 1 / 2, 1 / 2, 1 / -2, 1, 1000 )
        i.position.z = 1;
        i.position.y = 0;
        i.position.x = 0;
        i.up = new THREE.Vector3(0,1,0);
        return i;
    })()
    self.scScene  = (function(){
        var sc = new THREE.Scene();
        var planeGeo = new THREE.PlaneGeometry(1,1);
        var shaderMat = new THREE.ShaderMaterial({
            uniforms: {
                texture: {type: "t", value: self.frame}
            },
            vertexShader: document.getElementById("passThruShader").textContent,
            fragmentShader: document.getElementById("imageShader").textContent
        })
        var plane = new THREE.Mesh(planeGeo, shaderMat);
        plane.position.z = -1;
        plane.rotation.z = -1.5708; //Cancels the rotations from using a DataTexture initially. DataTextures and RenderTargets are not presented in GLSL the same way.
        sc.add(plane);
        return sc;
    })()
    self.frame = null;
    self.setFrame = function(nextFrame){
        self.frame = nextFrame;
        self.scScene.children[0].material.uniforms.texture.value = self.frame;
    }
    self.init  = function(){
        self.scRenderer.setSize(self.width, self.height);
        document.getElementById(id).appendChild(self.scRenderer.domElement);
    }
    self.render = function(){
        self.scRenderer.render(self.scScene, self.scCamera);
    }
    
}



/// PARALLEL

function getRayTextures(rays, samplesize){
    
    var  rayOrigins    = new Float32Array( rays.length * rays[0].length * samplesize * 4 );
    var  rayDirections = new Float32Array( rays.length * rays[0].length * samplesize * 4 );
    
    var count = 0;
    for(var row = 0; row < rays.length; row++){
        for(var col = 0; col < rays[row].length; col++){
            for(var sample = 0; sample < samplesize; sample++){
                rayOrigins[count] = rays[row][col][sample].origin.x;
                rayOrigins[count+1] = rays[row][col][sample].origin.y;
                rayOrigins[count+2] = rays[row][col][sample].origin.z;
                rayOrigins[count+3] = 1;
                rayDirections[count] = rays[row][col][sample].direction.x;
                rayDirections[count+1] = rays[row][col][sample].direction.y;
                rayDirections[count+2] = rays[row][col][sample].direction.z;
                rayDirections[count+3] = 1;
                count += 4;
            }
        }
    }
    
    var originTexture            = new THREE.DataTexture(rayOrigins, 400, 400, THREE.RGBAFormat, THREE.FloatType);
    var directionTexture         = new THREE.DataTexture(rayDirections, 400, 400, THREE.RGBAFormat, THREE.FloatType);
    originTexture.needsUpdate    = true;
    directionTexture.needsUpdate = true;
    
    return {origins: {type: 't', value: originTexture}, directions: {type: 't', value: directionTexture}}
}

/// END PARALLEL

/// INIT

var sc;
var material1 = new THREE.MeshBasicMaterial({ color: 0x2194ce });
var material2 = new THREE.MeshBasicMaterial({ color: 0x4FF5ff });
var material3 = new THREE.MeshBasicMaterial({ color: 0xff0055 });


var sphere1 = new Sphere(new THREE.Vector3(0,0.1,2), 0.5, 50,50, material2);
var sphere2 = new Sphere(new THREE.Vector3(-0.75,-0.3,1), 0.45, 50,50, material3);
var plane1  = new Plane(new THREE.Vector3(-0.5,-1,0), new THREE.Vector3(0,1,0), 4, 20, material1);
var light1  = new Light(new THREE.Vector3(0.5,4,10), new THREE.Color(0xffffff), new THREE.Color(0xffffff));

sc = new Scene();
sc.add(sphere1);
sc.add(sphere2);
sc.add(plane1);
sc.add(light1);


var view    = new View(new THREE.Vector3(0,0,3), new THREE.Vector3(0,0,-3).normalize(), new THREE.Vector3(0,1,0), 2, 400, 400);


function init(){
    
    var sc;
    var material1 = new THREE.MeshBasicMaterial({ color: 0x2194ce });
    var material2 = new THREE.MeshBasicMaterial({ color: 0x4FF5ff });
    var material3 = new THREE.MeshBasicMaterial({ color: 0xff0055 });


    var sphere1 = new Sphere(new THREE.Vector3(0,0.1,2), 0.5, 50,50, material2);
    var sphere2 = new Sphere(new THREE.Vector3(-0.75,-0.3,1), 0.45, 50,50, material3);
    var sphere3 = new Sphere(new THREE.Vector3(-0.75,0.5,1), 0.45, 50,50, material3);
    var plane1  = new Plane(new THREE.Vector3(-0.5,-1,0), new THREE.Vector3(0,1,0), 4, 20, material1);
    var light1  = new Light(new THREE.Vector3(0.5,4,10), new THREE.Color(0xffffff), new THREE.Color(0xffffff));

    sc = new Scene();
    sc.add(sphere1);
    sc.add(sphere2);
    sc.add(sphere3);
    sc.add(plane1);
    sc.add(light1);
    
    var scr = new SceneRenderer("WebGLCanvas", 400, 400);
    scr.init();
    
    var view    = new View(new THREE.Vector3(0,0,3), new THREE.Vector3(0,0,-3).normalize(), new THREE.Vector3(0,1,0), 2, 400, 400);
    var rt      = new RayTracer(view, sc);
    var gp = new GPGPU(scr, 400, 400, rt);
    
    var intersections = gp.getIntersections(getRayTextures(rt.view.getRays(), 1));
    var colors        = gp.getColors({intersections: {type: 't', value: intersections},
                                      bgColor: {type: 'c', value: new THREE.Color(0x000000)}});
    //rt.renderUniform("renderCanvasUniform");
    var arr = new Float32Array(400 * 400 * 4);
    scr.scRenderer.readRenderTargetPixels(intersections,0,0,400,400,arr );
    
    
    
    //var imgD = rt.view.imagePlane.generateImage();
    //var convD = convertData(imgD, 400, 400);
    scr.setFrame(colors);
    scr.render();
    var buff = new Float32Array(400 * 400 * 4)
    scr.scRenderer.readRenderTargetPixels(colors, 0, 0, 400, 400, buff);
    
    console.log("Hello");
        
}

function convertData(imageData, width, height){
    var arr = new Float32Array(width * height * 4);
    var count = 0;
    for(var i = 0; i < height; ++i){
        for(var k = 0; k < width; ++k){
            arr[count + 0] = imageData[count + 0]/255;
            arr[count + 1] = imageData[count + 1]/255;
            arr[count + 2] = imageData[count + 2]/255;
            arr[count + 3] = imageData[count + 3]/255;
            count += 4;
        }
    }
    var arrTexture = new THREE.DataTexture(arr, width, height, THREE.RGBAFormat, THREE.FloatType);
    arrTexture.needsUpdate = true;
    return arrTexture;
}

function genTexture(){
    
    var arr = new Float32Array(400 * 400 * 4);
    
    var count = 0;
    for(var i = 0; i < 400; ++i){
        for(var k = 0; k < 400; ++k){
            
            arr[count + 0] = (k > 200) ? 1.0 : 0.5;
            arr[count + 1] = 1.0;
            arr[count + 2] = (i > 200) ? 1.0 : 0.5;
            arr[count + 3] = 1.0;
            
            count += 4;
        }
    }
    
    
    var arrTexture = new THREE.DataTexture(arr, 400, 400, THREE.RGBAFormat, THREE.FloatType);
    arrTexture.needsUpdate = true;
    return arrTexture;
    
}