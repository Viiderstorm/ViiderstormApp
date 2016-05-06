//TODO:
// -add change listener to number of scene elements and recompile shaders if it changes
// -add change listener to scene object origins so objects can be moved around without recompiation


function GPGPU(SCR, width, height, rt){
    var self = this;
    
    //fields
    self.rayTracer          = rt;
    self.intersectionShader = null;
    self.shadowShader       = null;
    self.colorShader        = null;
    self.gpRenderer         = SCR.scRenderer;
    self.gpCamera           = (function(){
        var i      = new THREE.OrthographicCamera( 1 / -2, 1 / 2, 1 / 2, 1 / -2, 1, 1000 )
        i.position = new THREE.Vector3(0,0,1);
        return i;
    })()
    self.gpScene             = (function(){
        var sc       = new THREE.Scene();
        var planeGeo = new THREE.PlaneGeometry(1,1);
        var plane    = new THREE.Mesh(planeGeo, new THREE.MeshBasicMaterial());
        plane.position.z = -1;
        //plane.rotation.z = -1.5708;
        sc.add(plane);
        return sc;
    })()
    
    //init
    self.init                = function(){
        self.intersectionShader = self.getIntersectionShader(self.rayTracer.scene);
        self.shadowShader       = self.getShadowShader(self.rayTracer.scene);
        self.colorShader        = self.getColorShader(self.rayTracer.scene);
    }
    
    //class functions
    self.getIntersections = function(rayData){
        self.updateUniforms(self.intersectionShader, rayData );
        var intersections = self.compute(self.intersectionShader); // vec4(x,y,z, objectIndex)
        return intersections;
    }
    
    self.getShadows = function(intersectionData){
        self.updateUniforms(self.shadowShader, intersectionData );
        var shadows = self.compute(self.shadowShader);
        return shadows;
    }
    
    self.getColors        = function(data){
        self.updateUniforms(self.colorShader, data);
        var colors        = self.compute(self.colorShader);
        return colors;
    }
    
    self.compute          = function(shaderMat){
        self.gpScene.children[0].material = shaderMat;
        var output = new THREE.WebGLRenderTarget(width, height, {
            minFilter: THREE.NearestFilter,
			magFilter: THREE.NearestFilter,
			format: THREE.RGBAFormat,
			type: THREE.FloatType
        });
        self.gpRenderer.render(self.gpScene, self.gpCamera, output);
        return output;
    }
    
    self.updateUniforms        = function(shader, uniforms){
        shader.uniforms = $.extend(shader.uniforms, uniforms);
    }
    
    self.getIntersectionShader = function(scene){
        var sceneElems = scene.geometries;
        var shader = "precision mediump float;\nuniform sampler2D origins;\nuniform sampler2D directions;\nvarying vec2 vUv;\n";
        
        shader += getShaderObjectUniforms(sceneElems);
        shader += getShaderObjectStructs(sceneElems);
        
        //Closest Intersection
        shader += getClosestIntersectionFunction(sceneElems.length);
        shader += getShaderObjectIntersectionFunctions(sceneElems);
        
        //MAIN
        shader += "void main(){\n";
        shader += getShaderObjectStructInstances(sceneElems);
        shader += getShaderIntersectionArray(sceneElems);
        shader += "\tvec4 closestIntersection = findClosestIntersection(intersections);\n" +
                  "\tif(closestIntersection.w >= 0.0){\n" +
                  "\t\tgl_FragColor = closestIntersection;\n"+
                  "\t} else {\n" +
                  "\t\tgl_FragColor = vec4(0,0,0,9999999);\n"+ //use a exceptionally large number instead of negativies. Negatives don't return correctly.
                  "\t}\n" +
                  "}\n";
        
        var uniforms = {};
        for(var i = 0 ; i < sceneElems.length; i++){
            uniforms = $.extend(uniforms, getObjectIntersectionUniforms(sceneElems[i], i));
        }
        
        return new THREE.ShaderMaterial({
            uniforms: uniforms,
            vertexShader: document.getElementById("passThruShader").textContent,
            fragmentShader: shader
        });
    }
    
    self.getShadowShader = function(scene){
        var sceneElems = scene.geometries;
        var shader = "precision mediump float;\n" +
                     "uniform sampler2D intersections;\n"+
                     "varying vec2 vUv;\n";
        
        shader += getShaderObjectUniforms(sceneElems);
        shader += getShadowShaderUniforms(scene.lights);
        shader += getShaderObjectStructs(sceneElems);
        shader += getShaderObjectIntersectionFunctions(sceneElems);
        
        shader += getShadowShaderFunctions(sceneElems);
        
        shader += "void main(){\n" +
                  getShaderObjectStructInstances(sceneElems) +
                  "\tvec3 norm = getNormalForIntersection(texture2D(intersections, vUv.xy).xyzw);\n" + 
                  getShadowShaderCount(scene.lights) +
                  "\tgl_FragColor = vec4(norm,shadowCount);\n" +
                  "}\n"
        
        var uniforms = {};
        for(var i = 0 ; i < sceneElems.length; i++){
            uniforms = $.extend(uniforms, getObjectIntersectionUniforms(sceneElems[i], i));
        }
        for(var i = 0 ; i < scene.lights.length; i++){
            uniforms = $.extend(uniforms, getObjectShadowUniforms(scene.lights[i], i));
        }
        
        return new THREE.ShaderMaterial({
            uniforms: uniforms,
            vertexShader: document.getElementById("passThruShader").textContent,
            fragmentShader: shader
        });
        
        
    }
    
    self.getColorShader   = function(scene){
        var sceneElems = scene.geometries;
        var shader = "precision mediump float;\n" +
                     "uniform sampler2D intersections;\n" +
                     "uniform sampler2D normals;\n" +
                     "uniform vec3 bgColor;\n" +
                     "varying vec2 vUv;\n";
        for(var i = 0 ; i < sceneElems.length; i++){
            shader += getColorShaderUniforms(sceneElems[i], i);
        }
        shader += "void main(){\n";
        shader += "\tvec4 intersection = texture2D(intersections, vUv.xy);\n"
        shader += getColorShaderIfStatement(sceneElems);
        shader += "}\n";
        var uniforms = {}
        for(var i = 0 ; i < sceneElems.length; i++){
            uniforms = $.extend(uniforms, getObjectColorUniforms(sceneElems[i], i));
        }
        
        return new THREE.ShaderMaterial({
            uniforms: uniforms,
            vertexShader: document.getElementById("passThruShader").textContent,
            fragmentShader: shader
        });
        
    }
    
    self.init();
}

///////////////////////////
//      INTERSECTION
///////////////////////////
function getShaderObjectUniforms(objects){
    var str = "";
    for(var id = 0; id < objects.length; id++){
        if(objects[id] instanceof Sphere){
            str += "uniform vec3 sphere" + id + "_origin;\n" + 
                "uniform float sphere" + id + "_radius;\n" + 
                "uniform float sphere" + id + "_id;\n"; 
        } else if(objects[id] instanceof Plane){
            str += "uniform vec3 plane" + id + "_origin;\n" + 
                "uniform vec3 plane" + id + "_normal;\n" +
                "uniform vec3 plane" + id + "_pointOnPlane;\n" + 
                "uniform float plane" + id + "_id;\n";
        }
    }
    return str;
}

/*
    struct face{
        vec3 point1;
        vec3 point2;
        vec3 point3;    
    } planeN_triangleN, ... ;
*/
function getShaderObjectStructs(objects){
    //create struct
    var str = "struct face{\n" +
            "\t vec3 point1;\n" +
            "\t vec3 point2;\n" + 
            "\t vec3 point3;\n" +
            "} ";
    for(var i = 0; i < objects.length; i++){
        if(objects[i] instanceof Plane){         
            //define structs
            for(var k = 0; k < objects[i].faces.length; k++){
                str += "plane" + i + "_triangle" + k + 
                        ((k == objects[i].faces.length - 1) ? ";\n\n" : ",");
                
            }
        }
    }
    return str;
}

/*
    face planeN_faces[N];
    planeN_triangle0 = face(pos1,pos2,pos3);
    planeN_faces[0] = planeN_triangle0;
    .
    .
    .
    planeN_triangleN = face(pos1,pos2,pos3);
    planeN_faces[N] = planeN_triangleN;
*/
function getShaderObjectStructInstances(objects){
    var str = "";
    
    for(var id = 0 ; id < objects.length; id++){
        var geo = objects[id];
        if(geo instanceof Plane){
            
            //define array
            str = "\tface plane" + id + "_faces[" + geo.faces.length + "];\n";
                    
            //helper function                
            var faceToStr = function(s){
                var p1 = "vec3(" + geo.vertices[s.a].x + "," + geo.vertices[s.a].y + "," + geo.vertices[s.a].z + ")";
                var p2 =  "vec3(" + geo.vertices[s.b].x + "," + geo.vertices[s.b].y + "," + geo.vertices[s.b].z + ")";
                var p3 =  "vec3(" + geo.vertices[s.c].x + "," + geo.vertices[s.c].y + "," + geo.vertices[s.c].z + ")";
                return [p1,p2,p3].join(',');
            }
            
            //instantiate and place into array
            for(var k = 0; k < geo.faces.length; k++){
                //"planeN_triangleN = face(vec3,vec3,vec3);"
                str += "\tplane" + id + "_triangle" + k + " = " + 
                        "face" + "(" + faceToStr(geo.faces[k]) + ");\n";
                //"planeN_faces[N] = planeN_triangleN;"
                str += "\tplane" + id + "_faces[" + k + "] = " +
                        "plane"  + id + "_triangle" + k  + ";\n";
            }
                    
        }
    }
    return str;
}

function getShaderObjectIntersectionFunctions(objects){
    var sphere = "";
    var plane  = "";
    
    for(var i = 0 ; i < objects.length; i++){
        if(objects[i] instanceof Sphere && sphere == ""){
            sphere = "vec4 getSphereIntersection(vec3 rayOrigin, vec3 rayDirection, vec3 origin, float radius, float id){\n" +
                     "\tfloat A = pow(rayDirection.x, 2.0) + pow(rayDirection.y, 2.0) + pow(rayDirection.z, 2.0);\n" +
                     "\tfloat B = 2.0 * ((rayDirection.x * (rayOrigin.x - origin.x)) +" +
                                         "(rayDirection.y * (rayOrigin.y - origin.y)) +" +
                                         "(rayDirection.z * (rayOrigin.z - origin.z)));\n" +
                     "\tfloat C = pow(rayOrigin.x - origin.x, 2.0) + " +
                                 "pow(rayOrigin.y - origin.y, 2.0) + " +
                                 "pow(rayOrigin.z - origin.z, 2.0) - " +
                                 "pow(radius, 2.0);\n" +
                     "\tfloat W1 = (-B - sqrt(pow(B, 2.0) - (4.0 * A * C))) / (2.0 * A);\n" +
                     "\tfloat W2 = (-B + sqrt(pow(B, 2.0) - (4.0 * A * C))) / (2.0 * A);\n" +
                     "\tfloat W  = W1 > 0.0 ? W1 : W2;\n" +
                     "\tif( W > 0.0 ){\n \t\treturn vec4(rayOrigin + (rayDirection * W), id);\n" + 
                     "\t} else {\n \t\treturn vec4(0,0,0,-1);\n" + 
                     "\t}\n}\n\n";
        } else if (objects[i] instanceof Plane && plane == ""){
            plane = "bool isPointInPlane(vec3 intersectionPoint, face f){\n" +
                    "\tvec3 linestoVertices0 = f.point1 - intersectionPoint;\n" +
                    "\tvec3 linestoVertices1 = f.point2 - intersectionPoint;\n" +
                    "\tvec3 linestoVertices2 = f.point3 - intersectionPoint;\n" +
                    "\tfloat sum0 = degrees(acos( dot(normalize(linestoVertices0), normalize(linestoVertices1))));\n"+
                    "\tfloat sum1 = degrees(acos( dot(normalize(linestoVertices1), normalize(linestoVertices2))));\n"+
                    "\tfloat sum2 = degrees(acos( dot(normalize(linestoVertices2), normalize(linestoVertices0))));\n"+
                    "\tfloat totalSum = sum0 + sum1 + sum2;\n" +
                    "\tif(abs(totalSum - 360.0) < 0.1){\n" +
                    "\t\t return true;\n" +
                    "\t} else {\n" +
                    "\t\t return false;\n" +
                    "\t}\n" +
                    "}\n\n" +
                    "vec4 getPlaneIntersection(vec3 rayOrigin, vec3 rayDirection, vec3 pointOnPlane, face[2] faces, vec3 normal, float id){\n"+
                    "\tvec3 rayToPoint = rayOrigin - pointOnPlane;\n" + 
                    "\tfloat denom = dot(normal, rayDirection);\n" +
                    "\tfloat num   = dot(rayToPoint, normal);\n" +
                    "\tfloat t     = -(num/denom);\n" + 
                    "\tif(t >= 0.00000001){\n" +
                    "\t\tvec4 intersectionPoint = vec4(rayOrigin + (rayDirection * t), id);\n" +
                    "\t\tbool triangleCheck1;\n" +
                    "\t\tbool triangleCheck2;\n" +
                    "\t\ttriangleCheck1 = isPointInPlane(intersectionPoint.xyz, faces[0]);\n" +
                    "\t\ttriangleCheck2 = isPointInPlane(intersectionPoint.xyz, faces[1]);\n" +
                    "\t\tif(triangleCheck1 || triangleCheck2){\n" +
                    "\t\t\t return intersectionPoint;\n" +
                    "\t\t}\n" +
                    "\t}\n" +
                    "\treturn vec4(0,0,0,-1);\n" + 
                    "}\n\n";
        }
    }
    return sphere + plane;
}

function getClosestIntersectionFunction(numObjects){
    return "vec4 findClosestIntersection(vec4[" + numObjects +"] intersections){\n" +
            "\tvec4 max = vec4(9999999,9999999,9999999, -1.0);\n" +
            "\tfor(int i = 0; i < " + numObjects + "; i++){\n"+
            "\t\tif(intersections[i].w < 0.0){\n"+
            "\t\t\tcontinue;\n"+
            "\t\t}\n"+
            "\t\tfloat distance1 = distance(intersections[i].xyz, texture2D(origins, vUv.xy).xyz);\n" +
            "\t\tfloat distance2 = distance(max.xyz, texture2D(origins, vUv.xy).xyz);\n" +
            "\t\tif(distance1 < distance2){\n"+
            "\t\t\tmax = intersections[i].xyzw;\n" +
            "\t\t}\n" +
            "\t}\n" +
            "\treturn max;\n" +
            "}\n"
}

function getShaderIntersectionArray(objects){
    var str = "";
    str = "\tvec4 intersections[" + objects.length + "];\n"
    for(var i = 0; i < objects.length; i++){
        if(objects[i] instanceof Plane){
            str += "\tintersections[" + i + "] = getPlaneIntersection"+ 
                   "(texture2D(origins,vUv.xy).xyz,texture2D(directions,vUv.xy).xyz," +
                   "plane" + i + "_pointOnPlane.xyz," +
                   "plane" + i + "_faces," +
                   "plane" + i + "_normal," +
                   "plane" + i + "_id" +
                   ");\n";
        } else if (objects[i] instanceof Sphere){
            str += "\tintersections[" + i + "] = getSphereIntersection(texture2D(origins,vUv.xy).xyz,texture2D(directions,vUv.xy).xyz," +
                   "sphere" + i + "_origin.xyz," +
                   "sphere" + i + "_radius," +
                   "sphere" + i + "_id" +
                   ");\n";
        }
    }
    return str;
}

///////////////////////////
//        SHADOWS
///////////////////////////
function getShadowShaderUniforms(lights){
    var str = "";
    for(var i = 0 ; i < lights.length; i++){
        str += "uniform vec3 light" + i + "_origin;\n";
    }
    return str;
}

function getShadowShaderFunctions(objects){
    var str = "vec3 getNormalForIntersection(vec4 intersection){\n" 
    str += "\tvec3 norm;\n";
    var clauses = [];
    for(var i = 0 ; i < objects.length; i++){
        var qualifier = "";
        var eq = "";
        if(objects[i] instanceof Sphere){
            qualifier = "sphere" + i;
            eq = "norm = intersection.xyz - " + qualifier + "_origin.xyz;\n";
        } else if(objects[i] instanceof Plane){
            qualifier = "plane" + i;
            eq = "norm = " + qualifier + "_normal.xyz;\n";
        } else {
            continue;
        }
        var clause = "(intersection.w == " + qualifier + "_id){\n" +
                     "\t\t" + eq +
                     "\t}";
        clauses.push(clause);
    }
    str += "\tif" + clauses.join(" else if ");
    str += "\telse {\n" +
           "\t\tnorm = vec3(0.0,0.0,0.0);\n" +
           "\t}\n";
    str += "\treturn norm;\n}\n";
    
    str += "bool doesShadowRayIntersect(vec4 rayOrigin, vec3 light){\n"+
           "\tvec4 intersects[" + objects.length + "];\n" +
           "\tvec3 dir = normalize(light.xyz - rayOrigin.xyz);\n"+
           "\tvec3 newRayOrigin = rayOrigin.xyz + (dir * 0.00001);\n";
    for(var i = 0; i < objects.length; i++){
        if(objects[i] instanceof Plane){
            str += "\tface plane" + i + "_faces[2];\n" +
                   "\tplane" + i + "_faces[0] = plane" + i + "_triangle0;\n" + 
                   "\tplane" + i + "_faces[1] = plane" + i + "_triangle1;\n";
            str += "\tintersects[" + i + "] = getPlaneIntersection("+ 
                   "newRayOrigin.xyz," +
                   "dir," +
                   "plane" + i + "_pointOnPlane.xyz," +
                   "plane" + i + "_faces," +
                   "plane" + i + "_normal," +
                   "plane" + i + "_id" +
                   ");\n";
        } else if (objects[i] instanceof Sphere){
            str += "\tintersects[" + i + "] = getSphereIntersection(newRayOrigin.xyz,"+
                   "dir," +
                   "sphere" + i + "_origin.xyz," +
                   "sphere" + i + "_radius," +
                   "sphere" + i + "_id" +
                   ");\n";
        }
    }
    str += "\tfor(int i = 0; i < " + objects.length + "; i++){\n" +
           "\t\tif(intersects[i].w >= 0.0){\n" + 
           "\t\t\treturn true;\n" +
           "\t\t}\n" +
           "\t}\n" +
            "\treturn false;\n" +
            "}\n\n";
           
    
    return str;
    
    
}

function getShadowShaderCount(lights){
    var str = "\tbool shadows[" + lights.length + "];\n\tfloat shadowCount = 0.0;\n"
    for(var i = 0; i < lights.length; i++){
        str += "\tshadows[" + i + "] = doesShadowRayIntersect(texture2D(intersections, vUv.xy), light"+ i +"_origin);\n";
    }
    str += "\tfor(int i = 0; i < " + lights.length + "; i++){\n" +
           "\t\tif(shadows[i] == true){\n" + 
           "\t\t\tshadowCount += 1.0;\n" +
           "\t\t}\n" +
           "\t}\n";
    return str;
}

///////////////////////////
//     COLORS
///////////////////////////
function getColorShaderUniforms(geo, objectid){
    var str = "";
    if(geo instanceof Sphere){
        str = "uniform vec3 sphere" + objectid + "_color;\n" +
              "uniform float sphere" + objectid + "_id;\n";
    } else if (geo instanceof Plane){
        str = "uniform vec3 plane" + objectid + "_color;\n" +
              "uniform float plane" + objectid + "_id;\n";
    }
    return str;
}

function getColorShaderIfStatement(objects){
    var str = "";
    var clauses  = [];
    for(var i = 0 ; i < objects.length; i++){
        var qualifier = "";
        if(objects[i] instanceof Sphere){
            qualifier = "sphere" + i;
        } else if(objects[i] instanceof Plane){
            qualifier = "plane" + i;
        } else {
            continue;
        }
        var clause = "(intersection.w == " + qualifier + "_id && texture2D(normals, vUv.xy).w == 0.0){\n" +
                     "\t\tgl_FragColor = vec4(" + qualifier + "_color, 1.0);\n" +
                     "\t}";
        clauses.push(clause);
    }
    str = "\tif" + clauses.join(" else if ");
    str += "\telse {\n" +
           "\t\tgl_FragColor = vec4(bgColor, 1.0);\n" +
           "\t}\n";
           
    return str;
}


///////////////////////////
//          UNIFORMS
///////////////////////////
function getObjectIntersectionUniforms(geo, objectid){
    if(geo instanceof Sphere){
        var tmp = {};
        tmp["sphere"+objectid+"_origin"] = {type: 'v3', value: geo.origin};
        tmp["sphere"+objectid+"_radius"] = {type: 'f', value: geo.radius};
        tmp["sphere"+objectid+"_id"] = {type: 'f', value: objectid};
        return tmp;
    } else if (geo instanceof Plane){
        var tmp = {};
        tmp["plane"+objectid+"_origin"] = {type: 'v3', value: geo.origin};
        tmp["plane"+objectid+"_normal"] = {type: 'v3', value: geo.normal};
        tmp["plane"+objectid+"_pointOnPlane"] = {type: 'v3', value: geo.vertices[0]};
        tmp["plane"+objectid+"_id"] = {type: 'f', value: objectid};
        return tmp;
    }
    return null;
}

function getObjectShadowUniforms(geo, objectid){
    if (geo instanceof Light){
        var tmp = {};
        tmp["light"+objectid+"_origin"] = {type: 'v3', value: geo.origin};
        return tmp;
    }
    return null;
}

function getObjectColorUniforms(geo, objectid){
    if(geo instanceof Sphere){
        var tmp = {};
        tmp["sphere"+objectid+"_id"] = {type: 'f', value: objectid};
        tmp["sphere"+objectid+"_color"] = {type: 'c', value: geo.material.color1};
        return tmp;
    } else if (geo instanceof Plane){
        var tmp = {};
        tmp["plane"+objectid+"_id"] = {type: 'f', value: objectid};
        tmp["plane"+objectid+"_color"] = {type: 'c', value: geo.material.color1};
        return tmp;
    }
    return null;
}


function getRayTextures(rays, samplesize, gpgpu, width, height){
    
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
    
    var originTexture            = new THREE.DataTexture(rayOrigins, width, height, THREE.RGBAFormat, THREE.FloatType);
    var directionTexture         = new THREE.DataTexture(rayDirections, width, height, THREE.RGBAFormat, THREE.FloatType);
    originTexture.needsUpdate    = true;
    directionTexture.needsUpdate = true;
    
    var passthruFragShaderOrigins = "precision mediump float;\nuniform sampler2D origins;\n" +
                                    "varying vec2 vUv;\n" +
                                    "void main(){\n" +
                                    "\tgl_FragColor = texture2D(origins, vUv.yx).xyzw;\n" +
                                    "}\n";
    var passthruFragShaderDirections = "precision mediump float;\nuniform sampler2D directions;\n" +
                                    "varying vec2 vUv;\n" +
                                    "void main(){\n" +
                                    "\tgl_FragColor = texture2D(directions, vUv.yx).xyzw;\n" +
                                    "}\n";
    var originShaderMat = new THREE.ShaderMaterial({
        uniforms: { origins: {type: 't', value: originTexture}},
        vertexShader: document.getElementById("passThruShader").textContent,
        fragmentShader: passthruFragShaderOrigins
    });
    
    var directionShaderMat = new THREE.ShaderMaterial({
        uniforms: { directions: {type: 't', value: directionTexture}},
        vertexShader: document.getElementById("passThruShader").textContent,
        fragmentShader: passthruFragShaderDirections
    });

    var originsRT = gpgpu.compute(originShaderMat);
    var directionsRT = gpgpu.compute(directionShaderMat);                           
    
    return {origins: {type: 't', value: originsRT}, directions: {type: 't', value: directionsRT}}
}