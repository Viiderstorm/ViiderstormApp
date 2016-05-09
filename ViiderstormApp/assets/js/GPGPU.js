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
        self.intersectionShader = self.getIntersectionShader(self.rayTracer.scene, self.rayTracer);
        //self.shadowShader       = self.getShadowShader(self.rayTracer.scene);
        //self.colorShader        = self.getColorShader(self.rayTracer.scene);
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
    
    self.getIntersectionShader = function(scene, rt){
        var sceneElems = scene.geometries;
        var shader = "precision highp float;\nuniform sampler2D origins;\nuniform sampler2D directions;\nvarying vec2 vUv;\n";
        
        var glslUniforms = getShaderObjectUniforms(scene); 
        shader += glslUniforms.glsl;
        shader += getShaderObjectStructs(scene);
        shader += getShaderObjectIntersectionFunctions(scene);
        shader += getClosestIntersectionFunction(scene);
        shader += getShaderFindIntersectionFunction(scene);
        shader += getShadowShaderFunctions(scene);
        shader += getShaderRandom();
        shader += getColorShaderFunctions(scene, rt);
        
        //MAIN
        shader += "void main(){\n" +
                  "gl_FragColor = vec4(getColor(),1.0);"+                
                  "}\n";
        
        return new THREE.ShaderMaterial({
            uniforms: glslUniforms.uniforms,
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

/*
* Returns all uniform stuff required to render the scene
*/
function getShaderObjectUniforms(scene){
    var str = "";
    var uni = {};
    for(var id = 0; id < scene.geometries.length; id++){
        if(scene.geometries[id] instanceof Sphere){
            str += "uniform vec3 sphere" + id + "_origin;\n" + 
                   "uniform float sphere" + id + "_radius;\n" +
                   "uniform vec3 sphere" + id+ "_color;\n" + 
                   "uniform float sphere" + id + "_id;\n";
            var u = { }
            u["sphere"+id+"_origin"] = {type: 'v3', value: scene.geometries[id].origin};
            u["sphere"+id+"_radius"] = {type: 'f', value: scene.geometries[id].radius};
            u["sphere"+id+"_color"] = {type: 'c', value: scene.geometries[id].material.color1};
            u["sphere"+id+"_id"] = {type: 'f', value: id};
            uni = $.extend(uni, u);
        } else if(scene.geometries[id] instanceof Plane){
            str += "uniform vec3 plane" + id + "_origin;\n" + 
                "uniform vec3 plane" + id + "_normal;\n" +
                "uniform vec3 plane" + id + "_pointOnPlane;\n" + 
                "uniform vec3 plane" + id+ "_color;\n" +
                "uniform float plane" + id + "_id;\n";
            var u = { }
            u["plane"+id+"_origin"] = {type: 'v3', value: scene.geometries[id].origin};
            u["plane"+id+"_normal"] = {type: 'v3', value: scene.geometries[id].normal};
            u["plane"+id+"_pointOnPlane"] = {type: 'v3', value: scene.geometries[id].vertices[0]};
            u["plane"+id+"_color"] = {type: 'c', value: scene.geometries[id].material.color1};
            u["plane"+id+"_id"] = {type: 'f', value: id};
            uni = $.extend(uni, u);
        }
    }
    
    for(var id = 0; id < scene.lights.length; id++){
        if(scene.lights[id] instanceof Light){
            str += "uniform vec3 light" + id + "_origin;\n";
            str += "uniform vec3 light" + id + "_color;\n";
            str += "uniform vec3 light" + id + "_intensity;\n";
            var u = { }
            u["light"+id+"_origin"] = {type: 'v3', value: scene.lights[id].origin};
            u["light"+id+"_color"] = {type: 'c', value: scene.lights[id].color};
            u["light"+id+"_intensity"] = {type: 'c', value: scene.lights[id].intensity};
            uni = $.extend(uni, u);
        }
    }
    
    str += "uniform vec3 backgroundColor;\n";
    str += "uniform float mySeed;\n";
    var u = {};
    u['backgroundColor'] = {type: 'c', value: scene.backgroundColor};
    u['mySeed'] = {type: 'f', value: Math.random()*10};
    uni = $.extend(uni, u);
    
    var tmp = {
        glsl: str,
        uniforms: uni
    }
    
    return tmp;
}

/*
    Gets everything regarding the structs for the scene
*/
function getShaderObjectStructs(scene){
    
    //helper function                
    var faceToStr = function(geo, s){
        var p1 = "vec3(" + geo.vertices[s.a].x + "," + geo.vertices[s.a].y + "," + geo.vertices[s.a].z + ")";
        var p2 =  "vec3(" + geo.vertices[s.b].x + "," + geo.vertices[s.b].y + "," + geo.vertices[s.b].z + ")";
        var p3 =  "vec3(" + geo.vertices[s.c].x + "," + geo.vertices[s.c].y + "," + geo.vertices[s.c].z + ")";
        return [p1,p2,p3].join(',');
    }
    
    //create structs
    var str = "struct face{\n" +
              "\t vec3 point1;\n" +
              "\t vec3 point2;\n" + 
              "\t vec3 point3;\n" +
              "};\n\n" + 
              "struct material{\n" +
              "\tfloat ka;\n" +
              "\tfloat kd;\n" + 
              "\tfloat ks;\n" + 
              "\tfloat ke;\n" + 
              "\tfloat kr;\n" + 
              "\tfloat kt;\n" + 
              "\tfloat n;\n" +
              "};\n\n";
              
                
    for(var i = 0; i < scene.geometries.length; i++){
        if(scene.geometries[i] instanceof Plane){     
            str += "face plane" + i + "_faces[2];\n";    
            for(var k = 0; k < scene.geometries[i].faces.length; k++){
                str += "face plane" + i + "_triangle" + k + 
                       "= face(" + faceToStr(scene.geometries[i], scene.geometries[i].faces[k]) + ");\n";
            }
        }
        var mat = scene.geometries[i].material;
        var t = [
            mat.ka.toFixed(2),
            mat.kd.toFixed(2),
            mat.ks.toFixed(2),
            mat.ke.toFixed(2),
            mat.kr.toFixed(2),
            mat.kt.toFixed(2),
            mat.n.toFixed(2)
        ]
        str += "material mater" + i + " = material(" + t.join(',') + ");\n";
    }
    return str;
}

function getShaderObjectIntersectionFunctions(scene){
    var sphere = "";
    var plane  = "";
    
    for(var i = 0 ; i < scene.geometries.length; i++){
        if(scene.geometries[i] instanceof Sphere && sphere == ""){
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
        } else if (scene.geometries[i] instanceof Plane && plane == ""){
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

function getClosestIntersectionFunction(scene){
    return "vec4 findClosestIntersection(vec3 rayOrigin, vec4[" + scene.geometries.length +"] intersections){\n" +
            "\tvec4 max = vec4(9999999,9999999,9999999, -1.0);\n" +
            "\tfor(int i = 0; i < " + scene.geometries.length + "; i++){\n"+
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

function getShaderFindIntersectionFunction(scene){
    var str = "";
    str = "vec4 findIntersection(vec3 rayOrigin, vec3 rayDirection){\n" +
          "\tvec4 intersections[" + scene.geometries.length + "];\n"
    for(var i = 0; i < scene.geometries.length; i++){
        if(scene.geometries[i] instanceof Plane){
            str += "\tplane" + i + "_faces[0] = plane" + i + "_triangle0;\n" + 
                   "\tplane" + i + "_faces[1] = plane" + i + "_triangle1;\n" +
                   "\tintersections[" + i + "] = getPlaneIntersection"+ 
                   "(rayOrigin,rayDirection," +
                   "plane" + i + "_pointOnPlane.xyz," +
                   "plane" + i + "_faces," +
                   "plane" + i + "_normal," +
                   "plane" + i + "_id" +
                   ");\n";
        } else if (scene.geometries[i] instanceof Sphere){
            str += "\tintersections[" + i + "] = getSphereIntersection(rayOrigin,rayDirection," +
                   "sphere" + i + "_origin.xyz," +
                   "sphere" + i + "_radius," +
                   "sphere" + i + "_id" +
                   ");\n";
        }
    }
    
    str += "\treturn findClosestIntersection(rayOrigin.xyz, intersections);\n}\n\n";
    
    return str;
}

function getShadowShaderFunctions(scene){
    var objects = scene.geometries;
    var str = "vec3 getNormalForIntersection(vec4 intersection){\n" 
    str += "\tvec3 norm;\n";
    var clauses = [];
    for(var i = 0 ; i < objects.length; i++){
        var qualifier = "";
        var eq = "";
        if(objects[i] instanceof Sphere){
            qualifier = "sphere" + i;
            eq = "norm = normalize(intersection.xyz - " + qualifier + "_origin.xyz);\n";
        } else if(objects[i] instanceof Plane){
            qualifier = "plane" + i;
            eq = "norm = normalize(" + qualifier + "_normal.xyz);\n";
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
    str += "\treturn norm;\n}\n\n";
    
    str += "float findShadowRayIntersections(vec4 rayOrigin){\n"+
           "\tfloat numberOfShadows = 0.0;\n";
           
    for(var i = 0 ; i < scene.lights.length; i++){
        str += "\tvec3 dir" + i + "= normalize(light" + i + "_origin.xyz - rayOrigin.xyz);\n"+
               "\tvec3 newRayOrigin" + i + " = rayOrigin.xyz + (dir" + i + " * 0.00001);\n" +
               "\tif(findIntersection(newRayOrigin" + i + ", dir" + i + ").w >= 0.0){ numberOfShadows += 1.0;}\n\n";
    }
    
    str += "\treturn numberOfShadows;\n}\n\n";
           
    
    return str;  
}

function getColorShaderFunctions(scene, rt){
    var str = "";
    str += "vec3 getAmbientLight(vec4 intersection){\n"+
           "\tfloat comp = 1.0 / " + scene.lights.length.toFixed(1) + ";\n" +
           "\tvec3 sum;\n";
    var tmp = []
    for(var i = 0 ; i < scene.lights.length; i++){
        str += "\tvec3 color" + i + " = comp * light" + i + "_color;\n";
        tmp.push("color" + i); 
    }
    str += "\tsum = " + tmp.join("+") + ";\n";
    var terms = [];
    for(var i = 0 ; i < scene.geometries.length; i++){
        if(scene.geometries[i] instanceof Sphere){
            var tmp2 = "( intersection.w == sphere" + i + "_id){\n"+
                       "\t\tsum = mater" + i + ".ka * (sum * sphere" + i + "_color);\n"+
                       "\t}";
            terms.push(tmp2); 
        } else if (scene.geometries[i] instanceof Plane){
            var tmp2 = "( intersection.w == plane" + i + "_id){\n"+
                       "\t\tsum = mater" + i + ".ka * (sum * plane" + i + "_color);\n"+
                       "\t}";
            terms.push(tmp2);
        }
    }
    str += "\tif" + terms.join(" else if ") + " else {\n\t\tsum = vec3(1,1,1);\n\t}\n";
    str += "\treturn sum;\n"
    str += "}\n\n"
    
    str += "vec3 getDiffuseLight(vec4 intersection){\n"+
           "\tvec3 diffuseLight = vec3(0,0,0);\n" +
           "\tfloat constantKd = 0.0;\n" +
           "\tvec3 lightOrigins[" + scene.lights.length + "];\n" +
           "\tvec3 lightIntensities[" + scene.lights.length + "];\n";
    for(var i = 0; i < scene.lights.length; i++){
        str += "\tlightOrigins["  + i + "] = light" + i + "_origin;\n"
        str += "\tlightIntensities["  + i + "] = light" + i + "_intensity;\n"
    }
    
    var colorizer = function(type, id){
        return "( intersection.w == " + type + id + "_id){\n"+
            "\t\t\t diffuseLight += (lightIntensities[i] * " + type + id + "_color) " +
            " * dot(normalize(lightOrigins[i] - intersection.xyz).xyz, " + 
                    "getNormalForIntersection(intersection).xyz);\n" +
            "\t\t\tconstantKd = mater" + id + ".kd;\n" +
            "\t\t}";
    }
     
    str += "\tfor(int i = 0; i < " + scene.lights.length + "; i++){\n";
    var terms = [];
    for(var i = 0 ; i < scene.geometries.length; i++){
        if(scene.geometries[i] instanceof Sphere){
            var tmp2 = colorizer("sphere", i);
            terms.push(tmp2); 
        } else if (scene.geometries[i] instanceof Plane){
            var tmp2 = colorizer("plane", i);
            terms.push(tmp2);
        }
    }
    str += "\t\tif" + terms.join(" else if ") + " else {\n\t\t\tdiffuseLight = vec3(1,1,1);\n\t\t}\n";
    str += "\t}\n";
    str += "\treturn diffuseLight * constantKd;\n"
    str += "}\n\n"
    
    str += "vec3 getSpecularLight(vec4 intersection, vec3 direction){\n"+
           "\tvec3 specularLight = vec3(0,0,0);\n" +
           "\tfloat constantKs = 0.0;\n" +
           "\tvec3 lightOrigins[" + scene.lights.length + "];\n" +
           "\tvec3 lightIntensities[" + scene.lights.length + "];\n";
    for(var i = 0; i < scene.lights.length; i++){
        str += "\tlightOrigins["  + i + "] = light" + i + "_origin;\n"
        str += "\tlightIntensities["  + i + "] = light" + i + "_intensity;\n"
    }
    
    var colorizer1 = function(type, id){
        return "( intersection.w == " + type + id + "_id){\n"+
                "\t\t\tvec3 pointToLight = normalize(lightOrigins[i].xyz-intersection.xyz);\n" +
                "\t\t\tspecularLight += (lightIntensities[i] * "+
                        "pow(abs(dot(normalize(reflect(pointToLight, getNormalForIntersection(intersection))), (-direction.xyz))),"+
                        " mater" + id + ".ke));\n"+
                "\t\t\tconstantKs = mater" + id + ".ks;\n" +
                "\t\t}";
    }
     
    str += "\tfor(int i = 0; i < " + scene.lights.length + "; i++){\n";
    var terms = [];
    for(var i = 0 ; i < scene.geometries.length; i++){
        if(scene.geometries[i] instanceof Sphere){
            var tmp2 = colorizer1("sphere", i);
            terms.push(tmp2); 
        } else if (scene.geometries[i] instanceof Plane){
            var tmp2 = colorizer1("plane", i);
            terms.push(tmp2);
        }
    }
    str += "\t\tif" + terms.join(" else if ") + " else {\n\t\t\specularLight = vec3(1,1,1);\n\t\t}\n";
    str += "\t}\n";
    str += "\treturn specularLight * constantKs;\n"
    str += "}\n\n"
    
    str += "vec3 accumulateColor(vec3["+(rt.bounces+1)+"] colors, float["+(rt.bounces+1)+"] coeffs){\n"+
           "\tvec3 totalColor = colors["+(rt.bounces)+"];\n"+
           "\tfor(int i = "+(rt.bounces)+"; i > 0; i--){\n" +
           "\t\ttotalColor = colors[i-1]  + (totalColor * coeffs[i-1]);\n"+
           "\t}\n"+
           "\treturn totalColor;\n" +
            "}\n\n";
            
    str += "vec3 getTransmissionDirection(vec4 intersection, vec3 direction, vec3 normal, float refractIndex ){\n" +
           "\tfloat ni, nt, nit, DvN;\n" +
           "\tvec3 D;\n" +
           "\tvec3 N;\n" +
           "\tni = 1.0;\n" +
           "\tnt = refractIndex;\n" +
           "\tD = direction.xyz;\n" +
           "\tN = normal.xyz;\n" +
           "\tDvN = dot(-D, N);\n" +
           "\tif(DvN < 0.0){\n" +
           "\t\tni = refractIndex;\n" +
           "\t\tnt = 1.0;\n" +
           "\t\tN = -N;\n" +
           "\t\tDvN = dot(-D,N);\n"+
           "\t}" +
           "\tnit = ni/nt;\n"+
           "\tfloat discrim = sqrt(1.0 + (pow(nit,2.0) * (pow(DvN,2.0)-1.0) )  );\n" +
           "\tvec3 reflectRayDirection = (D * nit) + (N * ((DvN * nit) - discrim));\n" +
           "\tif(discrim < 0.0){\n" +
           "\t\treflectRayDirection = reflect(D, normalize(N));\n"+
           "\t}\n" +
           "\treturn reflectRayDirection.xyz;\n" +
           "}\n\n";
    
    str += "vec3 getColor(){\n" +
           "\tvec4 ray = texture2D(origins, vUv.xy);\n" +
           "\tvec3 direction = texture2D(directions, vUv.xy).xyz;\n" +
           "\tvec3 colors["+(rt.bounces+1)+"];\n" +
           "\tfloat coeffs["+(rt.bounces+1)+"];\n";
           
   for(var i= 0; i < (rt.bounces+1); i++){
       str += "\tcolors["+i+"] = vec3(0);\n"
       str += "\tcoeffs["+i+"] = 0.0;\n"
   }        
           
    str += "\tfor(int b = 0; b <= " + rt.bounces + "; b++){\n" +
           "\t\tvec4 intersection = findIntersection(ray.xyz, direction);\n"+
           "\t\tif(intersection.w < 0.0){colors[b] = backgroundColor; coeffs[b] = 1.0; break;}\n" +
           "\t\tfloat shadows = findShadowRayIntersections(intersection);\n"+
           "\t\tvec3 normal = getNormalForIntersection(intersection);\n"+
           "\t\tvec3 ambient = getAmbientLight(intersection);\n" +
           "\t\tvec3 diffuse = getDiffuseLight(intersection);\n" + 
           "\t\tvec3 specular = getSpecularLight(intersection, direction);\n" +  
           "\t\tvec3 curColor = ambient + diffuse + specular;\n" +
           "\t\tif(shadows > 0.0){curColor = ambient;}\n" +
           "\t\t";
           
    var colorizer2 = function(type, id){
        return "(intersection.w == " + type + id + "_id){\n" +
               "\t\t\tif(mater"+i+".kr > 0.0){\n" +
               "\t\t\t\tcolors[b] = curColor;\n" +
               "\t\t\t\tcoeffs[b] = mater"+i+".kr;\n" +
               "\t\t\t\tray = vec4(intersection.xyz + (normalize(-direction) * 0.0001), intersection.w);\n"+
               "\t\t\t\tdirection = reflect(direction, normalize(normal));\n"+
               "\t\t\t} else if(mater"+i+".kt > 0.0) {\n" +
               "\t\t\t\tcolors[b] = curColor;\n" +
               "\t\t\t\tcoeffs[b] = mater"+i+".kt;\n" +
               "\t\t\t\tray = vec4(intersection.xyz + (normalize(direction) * 0.000001), intersection.w);\n"+
               "\t\t\t\tdirection = getTransmissionDirection(ray, direction, normal, mater"+i+".n);\n"+
               "\t\t\t} else {\n" +
               "\t\t\t\tcolors[b] = curColor;\n" +
               "\t\t\t\tcoeffs[b] = 1.0;\n" +
               "\t\t\t\tray = vec4(intersection.xyz + (normalize(-direction) * 0.0001), intersection.w);\n"+
               "\t\t\t\tdirection = cosineWeightedDirection(dot(gl_FragCoord.xy, normal.xy)/mySeed, normal );\n"+
               "\t\t\t\tbreak;\n" +
               "\t\t\t}\n" +
               "\t\t} ";
    }
    
    var temp = [];
    for(var i = 0 ; i < scene.geometries.length; i++){
        if(scene.geometries[i] instanceof Sphere){
             var t = colorizer2("sphere", i);
             temp.push(t);
                    
        }  else if (scene.geometries[i] instanceof Plane){
             var t = colorizer2("plane", i);
             temp.push(t);
        }
    }
    str += "if" + temp.join(" else if ") + "\n";        
    str += "\t}\n" +
           "\t\treturn accumulateColor(colors, coeffs);\n" +
           "}\n\n";
    
    return str;       
}

function getShaderRandom(){
    var str =
    "float random(vec3 scale, float seed) {\n" +
    "\treturn fract(sin(dot(gl_FragCoord.xyz + seed, scale)) * 43758.5453 + seed);\n" +
    "}\n" + 
    "vec3 cosineWeightedDirection(float seed, vec3 normal) {\n" +
    "\tfloat u = random(vec3(12.9898, 78.233, 151.7182), seed);\n" +
    "\tfloat v = random(vec3(63.7264, 10.873, 623.6736), seed);\n" +
    "\tfloat r = sqrt(u);\n" +
    "\tfloat angle = 6.283185307179586 * v;\n" +
    "\tvec3 sdir, tdir;\n" +
    "\tif (abs(normal.x) < .5) {\n" +
    "\t    sdir = cross(normal, vec3(1, 0, 0));\n" +
    "\t} else {\n" +
    "\t    sdir = cross(normal, vec3(0, 1, 0));\n" +
    "\t}\n" +
    "\ttdir = cross(normal, sdir);\n" +
    "\treturn r * cos(angle) * sdir + r * sin(angle) * tdir + sqrt(1. - u) * normal;\n" +
    "}\n" +
    "vec3 uniformlyRandomDirection(float seed) {\n" +
    "\tfloat u = random(vec3(12.9898, 78.233, 151.7182), seed);\n" +
    "\tfloat v = random(vec3(63.7264, 10.873, 623.6736), seed);\n" +
    "\tfloat z = 1.0 - 2.0 * u;\n" +
    "\tfloat r = sqrt(1.0 - z * z);\n" +
    "\tfloat angle = 6.283185307179586 * v;\n" +
    "\treturn vec3(r * cos(angle), r * sin(angle), z);\n" +
    "}\n" +
    "vec3 uniformlyRandomVector(float seed) {\n" +
    "\treturn uniformlyRandomDirection(seed) * sqrt(random(vec3(36.7539, 50.3658, 306.2759), seed));\n" +
    "}\n";
    return str;

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